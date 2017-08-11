// estab exports elasticsearch fields as tab separated values
package main
import (
  "bufio"
  "bytes"
  "flag"
  "fmt"
  "io"
  "log"
  "os"
  "strconv"
  "strings"
  "sync"
  "regexp"
  "github.com/akotlar/bystro-utils/parse"
  "sort"
  "path/filepath"
  "compress/gzip"
)

const concurrency int = 16
const chromIdx int = 0
const posIdx int = 1
const idIdx int = 2
const refIdx int = 3
const altIdx int = 4
const qualIdx int = 5
const filterIdx int = 6
const infoIdx int = 7
const formatIdx int = 8

type Config struct {
  bedPath string
  vcfGlob string
}

func setup(args []string) *Config {
  config := &Config{}
  flag.StringVar(&config.bedPath, "bedPath", "", "The input file path (optional: default is stdin)")
  flag.StringVar(&config.vcfGlob, "vcfGlob", "", "The output path for the JSON output (optional)")
  
  flag.CommandLine.Parse(os.Args[1:])

  return config
}

func init() {
  log.SetFlags(0)
}

// NOTE: For now this only supports \n end of line characters
// If we want to support CLRF or whatever, use either csv package, or set a different delimiter
func main() {
  config := setup(nil)

  ranges := readBed(config.bedPath)

  seen := make(map[string]int)
  readVcf(config.vcfGlob, ranges, func(data []string) {
    if seen[data[0]] == 0 {
      fmt.Print(data[1])
      seen[data[0]]++
    }
  })

  for id, count := range seen {
    log.Printf("Saw %s %d times \n", id, count)
  }
}

func readBed (bedPath string) map[string][][]int {
  if bedPath == "" {
    log.Fatal("bedPath required");
  }

  inFh, err := os.Open(bedPath)

  if err != nil {
    log.Fatal("Couldn't open", bedPath)
  }

  ranges := make(map[string][][]int)
  // starts := make(map[string][]int)
  // ends := make(map[string][]int)

  reader := bufio.NewReader(inFh)

  var chr string

  for {
    // http://stackoverflow.com/questions/8757389/reading-file-line-by-line-in-go
    // http://www.jeffduckett.com/blog/551119d6c6b86364cef12da7/golang---read-a-file-line-by-line.html
    // Scanner doesn't work well, has buffer restrictions that we need to manually get around
    // and we don't expect any newline characters in a Seqant output body
    row, err := reader.ReadString('\n') // 0x0A separator = newline

    if err == io.EOF {
      break
    } else if err != nil {
      log.Fatal(err)
    } else if row == "" {
      // This shouldn't occur, however, in case
      continue
    }

    // remove the trailing \n or \r
    // equivalent of chomp https://groups.google.com/forum/#!topic/golang-nuts/smFU8TytFr4
    record := strings.Split(row[:len(row) - 1], "\t")

    if row[0:3] == "chr" {
      chr = record[0]
    } else {
      var buffer bytes.Buffer

      buffer.WriteString("chr")
      buffer.WriteString(record[0])

      chr = buffer.String()
    }

    start, err := strconv.Atoi(record[1])

    if err != nil {
      log.Fatal("couldn't convert", record)
    }

    end, err := strconv.Atoi(record[2])

    if err != nil {
      log.Fatal("couldn't convert", record)
    }

    if ranges[chr] == nil {
      ranges[chr] = [][]int{{start, end}}
    } else {
      ranges[chr] = append(ranges[chr], []int{start,end})
    }
  }

  for _, val := range ranges {
    sort.Slice(val, func(i, j int) bool {
      return val[i][0] < val[j][0]
    })
  }

  return ranges
}

func readVcf (vcfGlob string, ranges map[string][][]int, resultFunc func(data []string)) {
  foundHeader := false

  files, err := filepath.Glob(vcfGlob)

  if err != nil {
    log.Fatal(err)
  }

  // Read buffer
  // workQueue := make(chan string, 100)
  complete := make(chan bool)
  // Write buffer
  results := make(chan []string, 100)
  var wg sync.WaitGroup

  inFh, err := os.Open(files[0])

  if err != nil {
    log.Fatal(err)
  }

  gzRead, err := gzip.NewReader(inFh)

  if err != nil {
    log.Fatal(err)
  }

  reader := bufio.NewReader(gzRead)

  endOfLineByte, numChars, versionLine, err := parse.FindEndOfLine(reader, "")

  if err != nil {
    log.Fatal(err)
  }

  vcfMatch, err := regexp.MatchString("##fileformat=VCFv4", versionLine)

  if err != nil {
    log.Fatal(err)
  }

  if !vcfMatch {
    log.Fatal("Not a VCF file")
  }

  fmt.Println(versionLine)

  // Get the header
  for {
    // http://stackoverflow.com/questions/8757389/reading-file-line-by-line-in-go
    // http://www.jeffduckett.com/blog/551119d6c6b86364cef12da7/golang---read-a-file-line-by-line.html
    // Scanner doesn't work well, has buffer restrictions that we need to manually get around
    // and we don't expect any newline characters in a Seqant output body
    row, err := reader.ReadString(endOfLineByte) // 0x0A separator = newline

    if err == io.EOF {
      break
    } else if err != nil {
      log.Fatal(err)
    } else if row == "" {
      // This shouldn't occur, however, in case
      continue
    }

    // remove the trailing \n or \r
    // equivalent of chomp https://groups.google.com/forum/#!topic/golang-nuts/smFU8TytFr4
    record := strings.Split(row[:len(row) - numChars], "\t")

    if foundHeader == false {
      if row[0] == '#' {
        fmt.Print(row)
      }

      if record[chromIdx] == "#CHROM" {
        foundHeader = true
        break
      }
    }
  }

  if !foundHeader {
    log.Fatal("No header found")
  }

  // Now read them all off, concurrently.
  for i := 0; i < len(files); i++ {
    go processLines(files[i], ranges, results, complete)
  }

  wg.Add(1)
  go func() {
    defer wg.Done()
    for line := range results {
      resultFunc(line)
    }
  }()

  

  // Wait for everyone to finish.
  for i := 0; i < concurrency; i++ {
    <-complete
  }

  close(results)

  wg.Wait()
}

func processLines(file string, ranges map[string][][]int, results chan []string, complete chan bool) {
  inFh, err := os.Open(file)

  if err != nil {
    log.Fatal(err)
  }

  gzFile, err := gzip.NewReader(inFh)

  if err != nil {
    log.Fatal(err)
  }

  reader := bufio.NewReader(gzFile)

  var pos int
  var arr [][]int

  endOfLineByte, numChars, versionLine, err := parse.FindEndOfLine(reader, "")

  if err != nil {
    log.Fatal(err)
  }

  vcfMatch, err := regexp.MatchString("##fileformat=VCFv4", versionLine)

  if err != nil || vcfMatch == false {
    log.Fatal("Not a VCF file", file)
  }

  for {
    row, err := reader.ReadString(endOfLineByte) // 0x0A separator = newline

    if err == io.EOF {
      break
    } else if err != nil {
      log.Fatal(err)
    } else if row == "" || row[0] == '#' {
      // We may have not closed the pipe, but not have any more information to send
      // Wait for EOF
      continue
    }

    record := strings.Split(row[:len(row) - numChars], "\t")

    if record[6] != "PASS" {
      continue
    }

    if row[0:3] == "chr" {
      arr = ranges[record[0]]
    } else {
      var buffer bytes.Buffer

      buffer.WriteString("chr")
      buffer.WriteString(record[0])

      record[0] = buffer.String()

      arr = ranges[record[0]]
    }

    if arr == nil {
      continue
    }

    pos, err = strconv.Atoi(record[1])

    if err != nil {
      log.Fatal(err)
    }

    if pos < arr[0][0] || pos > arr[len(arr) - 1][1] {
      continue
    }

    for i := 0; i < len(arr); i++ {
      if pos >= arr[i][0] && pos <= arr[i][1] {
        var buffer bytes.Buffer
        buffer.WriteString(record[0])
        buffer.WriteRune('_')
        buffer.WriteString(record[1])
        buffer.WriteRune('_')
        buffer.WriteString(record[3])
        buffer.WriteRune('_')
        buffer.WriteString(record[4])

        results <- []string{buffer.String(), row}
        break
      } else if pos > arr[i][1] {
        break
      }
    }
  }

  // log.Println("Worker hit, missed this many times: ", hitCount, missCount)
  // Let the main process know we're done.
  complete <- true
}
