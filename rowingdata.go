package gorow

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
)

type Rowingdata struct {
	Name               string
	timestamp          []float64
	distance           []float64
	spm                []float64
	hr                 []float64
	pace               []float64
	power              []float64
	drivelength        []float64
	strokedistance     []float64
	drivetime          []float64
	dragfactor         []int
	strokerecoverytime []float64
	averageforce       []float64
	peakforce          []float64
	velo               []float64
	lapnr              []int
	intervaltime       []float64
	calories           []float64
	workoutstate       []int64
}

func ReadCSV(f string) {
	csvFile, _ := os.Open(f)
	defer csvFile.Close()
	reader := csv.NewReader(csvFile)
	// should be for record,err = reader.Read
	// dict[header[i]] = record[i]
	// https://gist.github.com/drernie/5684f9def5bee832ebc50cabb46c377a
	// rows = append(rows.dict)
	records, err := reader.ReadAll()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Print(records)
}
