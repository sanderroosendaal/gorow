package gorow

import (
	"compress/gzip"
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"log"
	"math"
	"net/http"
	"net/url"
	"os"
	"path/filepath"
	"reflect"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/fatih/structs"
	"github.com/xitongsys/parquet-go-source/local"
	"github.com/xitongsys/parquet-go/parquet"
	"github.com/xitongsys/parquet-go/reader"
	"github.com/xitongsys/parquet-go/writer"
)

// LbstoN convert lbs of force to Newton
const LbstoN = 4.44822

func reverseMap(m map[string]string) map[string]string {
	n := make(map[string]string)
	for k, v := range m {
		n[v] = k
	}
	return n
}

// StrokeRecord sort of dataframe
type StrokeRecord struct {
	Timestamp          float64 `rowingdata:"TimeStamp (sec)" parquet:"name=time, type=DOUBLE, encoding=PLAIN_DICTIONARY"`
	Distance           float64 `rowingdata:" Horizontal (meters)" parquet:"name=distance, type=DOUBLE"`
	Spm                float64 `rowingdata:" Cadence (stokes/min)" parquet:"name=spm, type=DOUBLE"`
	Hr                 float64 `rowingdata:" HRCur (bpm)" parquet:"name=hr, type=DOUBLE"`
	Pace               float64 `rowingdata:" Stroke500mPace (sec/500m)" parquet:"name=pace, type=DOUBLE"`
	Power              float64 `rowingdata:" Power (watts)" parquet:"name=power, type=DOUBLE"`
	Drivelength        float64 `rowingdata:" DriveLength (meters)" parquet:"name=drivelength, type=DOUBLE"`
	Strokedistance     float64 `rowingdata:" StrokeDistance (meters)" parquet:"name=distanceperstroke, type=DOUBLE"`
	Drivetime          float64 `rowingdata:" drivetime"`
	Dragfactor         int64   `rowingdata:" DragFactor"`
	Strokerecoverytime float64 `rowingdata:" StrokeRecoveryTime (ms)"`
	Workperstroke      float64 `rowingdata:" WorkPerStroke (joules)" parquet:"name=driveenergy, type=DOUBLE"`
	Averageforce       float64 `rowingdata:" AverageDriveForce (lbs)" parquet:"name=averageforce, type=DOUBLE"`
	Peakforce          float64 `rowingdata:" PeakDriveForce (lbs)" parquet:"name=peakforce, type=DOUBLE"`
	Velo               float64 `rowingdata:" Speed (m/sec)" parquet:"name=velo, type=DOUBLE"`
	Lapnr              int64   `rowingdata:" lapIdx"`
	Intervaltime       float64 `rowingdata:" ElapsedTime (sec)"`
	Calories           float64 `rowingdata:" Calories (kCal)"`
	Workoutstate       float64 `rowingdata:" WorkoutState" parquet:"name=workoutstate, type=DOUBLE"`
	Latitude           float64 `rowingdata:" latitude"`
	Longitude          float64 `rowingdata:" longitude"`
	Bearing            float64 `rowingdata:" bearing"`
	Nowindpace         float64 `rowingdata:"nowindpace" parquet:"name=nowindpace, type=DOUBLE"`
	Equivergpower      float64 `rowingdata:"Equiv erg Power" parquet:"name=equivergpower, type=DOUBLE"`
	Modelpower         float64 `rowingdata:"power (model)"`
	Modelfavg          float64 `rowingdata:"averageforce (model)"`
	Modeldrivelength   float64 `rowingdata:"drivelength (model)"`
	Vwind              float64 `rowingdata:"vwind"`
	Winddirection      float64 `rowingdata:"winddirection"`
	Vstream            float64 `rowingdata:"vstream"`
}

// GetField gets field value as float from StrokeRecord
func (s *StrokeRecord) GetField(field string) (float64, error) {
	r := reflect.ValueOf(s)
	f := reflect.Indirect(r).FieldByName(field)
	tip := f.Type().Name()
	switch tip {
	case "float64":
		return float64(f.Float()), nil
	case "int":
		return float64(f.Int()), nil
	}
	return 0, errors.New("GetField returned an invalid type")
}

// StrokeFieldMapping maps StrokeRecord field names to CSV header names
func StrokeFieldMapping(field string) (string, error) {
	s := &StrokeRecord{}
	r := reflect.TypeOf(s)
	f, ok := r.Elem().FieldByName(field)

	if !ok {
		return "", errors.New("FieldByName didn't work")
	}

	return f.Tag.Get("rowingdata"), nil
}

func getfloatrecord(s string) (float64, error) {
	return strconv.ParseFloat(strings.TrimSpace(s), 64)
}

func getintrecord(s string) (int64, error) {
	res, err := strconv.ParseInt(strings.TrimSpace(s), 10, 64)
	return res, err
}

func exists(name string) bool {
	if _, err := os.Stat(name); err != nil {
		if os.IsNotExist(err) {
			return false
		}
	}
	return true
}

func csvwriter(records [][]string, f string) (ok bool, err error) {
	// file does not exist or overwrite is set to true
	csvFile, err := os.OpenFile(f, os.O_WRONLY|os.O_CREATE, 0777)
	if err != nil {
		return false, err
	}
	defer csvFile.Close()

	// create writer
	w := csv.NewWriter(csvFile)
	// write header
	w.WriteAll(records)

	if err := w.Error(); err != nil {
		return false, err
	}
	return true, nil
}

func gzipper(f string) (ok bool, err error) {
	// file does not exist or overwrite is set to true
	// file does not exist or overwrite is set to true

	// zip it
	reader, err := os.Open(f)
	if err != nil {
		return false, err
	}
	defer reader.Close()

	filename := filepath.Base(f)
	target := fmt.Sprintf("%s.gz", f)
	writer, err := os.Create(target)
	if err != nil {
		return false, err
	}
	defer writer.Close()

	archiver := gzip.NewWriter(writer)
	archiver.Name = filename
	defer archiver.Close()

	_, err = io.Copy(archiver, reader)
	if err != nil {
		return false, err
	}

	err = os.Remove(f)
	if err != nil {
		return false, err
	}

	return true, nil
}

// WriteParquet writes data to Parquet
func WriteParquet(strokes []StrokeRecord, f string, overwrite bool, gz bool) (ok bool, err error) {
	if exists(f) && !overwrite {
		err := errors.New("File exists and overwrite was set to false")
		return false, err
	}
	fw, err := local.NewLocalFileWriter(f)
	defer fw.Close()
	if err != nil {
		return false, err
	}
	//parameters: writer, type of struct, size
	pw, err := writer.NewParquetWriter(fw, new(StrokeRecord), 1)
	if err != nil {
		return false, err
	}
	//compression type
	pw.CompressionType = parquet.CompressionCodec_UNCOMPRESSED
	if gz {
		pw.CompressionType = parquet.CompressionCodec_GZIP
	}
	// pw.RowGroupSize = 128 * 1024 * 1024 //128M

	for _, d := range strokes {
		if err = pw.Write(d); err != nil {
			return false, err
		}
	}
	if err = pw.WriteStop(); err != nil {
		return false, err
	}
	return true, nil
}

// WriteCSV writes data to file
func WriteCSV(strokes []StrokeRecord, f string, overwrite bool, gz bool) (ok bool, err error) {
	if exists(f) && !overwrite {
		err := errors.New("File exists and overwrite was set to false")
		return false, err
	}
	// create records
	var records [][]string
	names := structs.Names(&StrokeRecord{})
	var header []string
	for _, name := range names {
		rowingdataname, err := StrokeFieldMapping(name)
		if err != nil {
			return false, err
		}
		header = append(header, rowingdataname)
	}
	records = append(records, header)
	for _, stroke := range strokes {
		var record []string
		for _, name := range names {
			value, err := stroke.GetField(name)
			if err != nil {
				value = 0
			}
			record = append(record, fmt.Sprintf("%f", value))
		}
		records = append(records, record)
	}

	// gzip
	if gz {
		f = f[:len(f)-3]
		_, err := csvwriter(records, f)
		if err != nil {
			return false, err
		}
		return gzipper(f)
	}

	return csvwriter(records, f)

}

// ReadParquet reads rowing data into data frame
// from Parquet file
func ReadParquet(f string) ([]StrokeRecord, error) {
	fr, err := local.NewLocalFileReader(f)
	defer fr.Close()
	if err != nil {
		return nil, err
	}
	pr, err := reader.NewParquetReader(fr, new(StrokeRecord), 4)
	if err != nil {
		return nil, err
	}
	num := int(pr.GetNumRows())
	outp := make([]*StrokeRecord, num)

	if err := pr.Read(&outp); err != nil {
		return nil, err
	}

	pr.ReadStop()

	out := make([]StrokeRecord, num)
	for i := range out {
		out[i] = *outp[i]
	}
	return out, nil
}

// ReadCSV reads rowing data into data frame
// from CSV file or gzipped CSV file with extension .csv.gz
func ReadCSV(f string) ([]StrokeRecord, error) {
	// get extension
	ext := filepath.Ext(f)

	csvFile, err := os.Open(f)
	if err != nil {
		return []StrokeRecord{}, errors.New("ReadCSV: Unable to open file")
	}
	defer csvFile.Close()
	var rg io.Reader
	if ext == ".gz" {
		rg, err = gzip.NewReader(csvFile)
		if err != nil {
			return []StrokeRecord{}, errors.New("ReadCSV: Unable to open gzip")
		}
		// defer rg.Close()
	} else {
		rg = csvFile
	}
	reader := csv.NewReader(rg)
	// should be for record,err = reader.Read
	// dict[header[i]] = record[i]
	// https://gist.github.com/drernie/5684f9def5bee832ebc50cabb46c377a
	// rows = append(rows.dict)
	var header []string
	var rows []StrokeRecord

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return rows, errors.New("Hit CSV reader error, returning partial result")
		}
		if header == nil {
			header = record
		} else {
			var row StrokeRecord
			for i := range header {
				switch header[i] {
				case "TimeStamp (sec)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Timestamp = f
					}
				case " lapIdx":
					if f, err := getintrecord(record[i]); err == nil {
						row.Lapnr = f
					}
				case " ElapsedTime (sec)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Intervaltime = f
					}
				case " Horizontal (meters)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Distance = f
					}
				case " Stroke500mPace (sec/500m)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Pace = f
					}
				case " Cadence (stokes/min)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Spm = f
					}
				case " Cadence (strokes/min)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Spm = f
					}
				case " HRCur (bpm)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Hr = f
					}
				case " Power (watts)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Power = f
					}
				case " Calories (kCal)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Calories = f
					}
				case " Speed (m/sec)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Velo = f
					}
				case " StrokeDistance (meters)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Strokedistance = f
					}
				case " DriveLength (meters)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Drivelength = f
					}
				case " StrokeRecoveryTime (ms)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Strokerecoverytime = f
					}
				case " WorkPerStroke (joules)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Workperstroke = f
					}
				case " AverageDriveForce (lbs)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Averageforce = f
					}
				case " AverageDriveForce (N)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Averageforce = f / LbstoN
					}
				case " PeakDriveForce (lbs)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Peakforce = f
					}
				case " PeakDriveForce (N)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Peakforce = f / LbstoN
					}
				case " DragFactor":
					if f, err := getintrecord(record[i]); err == nil {
						row.Dragfactor = f
					}
				case " latitude":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Latitude = f
					}
				case " longitude":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Longitude = f
					}
				case "Nowindpace":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Nowindpace = f
					}
				case "Equiv erg Power":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Equivergpower = f
					}
				case "power (model)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Modelpower = f
					}
				case "averageforce (model)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Modelfavg = f
					}
				case "drivelength (model)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Modeldrivelength = f
					}
				case "vwind":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Vwind = f
					}
				case "winddirection":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Winddirection = f
					}
				case "vstream":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.Vstream = f
					}
				}
			}
			if row.Velo == 0 && row.Pace != 0 {
				row.Velo = 500. / row.Pace
			}
			if row.Workperstroke == 0 && row.Power != 0 && row.Spm != 0 {
				row.Workperstroke = 60. * row.Power / row.Spm
			}
			rows = append(rows, row)
		}
	}
	return rows, nil
}

// reportprogress
func postprogress(secret, progressurl string, progress int) (statuscode int) {
	postData := url.Values{}
	postData.Set("secret", secret)
	postData.Set("value", fmt.Sprintf("%d", progress))

	req, err := http.NewRequest("POST", progressurl, strings.NewReader(postData.Encode()))
	if err != nil {
		return 408
	}
	req.Header.Set("content-type", "application/x-www-form-urlencoded")

	client := &http.Client{Timeout: 5 * time.Second}
	resp, err := client.Do(req)
	if err != nil {
		return 408
	}
	defer resp.Body.Close()

	return resp.StatusCode
}

// OTWSetPower adds power for OTW rows
func OTWSetPower(
	strokes []StrokeRecord,
	c *Crew,
	rg *Rig,
	secret string,
	progressurl string,
	powermeasured bool,
	verbose bool) error {

	// a blocking channel to keep concurrency under control
	semaphoreChan := make(chan struct{}, 4)
	defer close(semaphoreChan)

	// a wait group enables the main process a wait for goroutines to finish
	wg := sync.WaitGroup{}

	// counter channel
	var counter = make(chan struct{})
	var done = make(chan struct{})

	// start a counter/timer
	if len(progressurl) > 0 {
		go func(aantal int) {
			var cntr int
			ticker2 := time.NewTicker(1 * time.Second)
			for {
				select {
				case <-counter:
					cntr++
				case <-done:
					return
				case <-ticker2.C:
					perc := 100 * cntr / (aantal - 1)
					if verbose {
						log.Printf("Percentage done: %d\n", perc)
					}
					postprogress(secret, progressurl, perc)
				}
			}
		}(len(strokes))
	} else {
		close(counter)
		close(done)
	}

	for i, stroke := range strokes {
		// increment the wait group internal counter
		wg.Add(1)

		go func(i int, stroke StrokeRecord, c *Crew, rg *Rig) {
			// block until the semaphore channel has room
			// this could also be moved out of the goroutine
			// which would make sense if the list is huge
			semaphoreChan <- struct{}{}

			c.Tempo = stroke.Spm
			res, err := PhysGetPower(
				stroke.Velo, c, rg, stroke.Bearing,
				stroke.Vwind, stroke.Winddirection, stroke.Vstream,
			)
			//res, err := PhysGetPower(stroke.velo, c, rg, 0, 0, 0, 0)
			if err == nil {

				pwr := res[0]
				frc := res[2]
				nowindp := res[3]
				if !powermeasured {
					strokes[i].Power = pwr
					strokes[i].Averageforce = frc / LbstoN
				}
				strokes[i].Nowindpace = nowindp
				strokes[i].Modelpower = pwr
				strokes[i].Modelfavg = frc / LbstoN
			}

			// tell the wait group that we be done
			wg.Done()
			// send signal to counter
			if len(progressurl) > 0 {
				counter <- struct{}{}
			}
			// clear a spot in the semaphore channel
			<-semaphoreChan
		}(i, stroke, c, rg)
	}
	// wait for all the goroutines to be done
	wg.Wait()
	if len(progressurl) > 0 {
		done <- struct{}{}
	}
	err := smoothnowindpace(strokes, 3)
	if err != nil {
		return err
	}
	return nil
}

func averageNowindpace(strokes []StrokeRecord) float64 {
	p := 0.0
	var counter int
	for _, stroke := range strokes {
		if !math.IsNaN(stroke.Nowindpace) {
			p += stroke.Nowindpace
		} else {
			counter++
		}
	}
	return p / float64(len(strokes)-counter)
}

// AveragePower calculates average power
func AveragePower(
	strokes []StrokeRecord,
) float64 {
	power := 0.0
	var counter int
	for _, stroke := range strokes {
		if !math.IsNaN(stroke.Power) {
			power += stroke.Power
		} else {
			counter++
		}
	}
	return power / float64(len(strokes)-counter)
}

func cumulativedistance(strokes []StrokeRecord) ([]float64, error) {
	var cumdist []float64
	cumdist = append(cumdist, 0)
	curdist := 0.0
	for i := 1; i < len(strokes); i++ {
		delta := strokes[i].Distance - strokes[i-1].Distance
		if delta > 0 {
			curdist += delta
			cumdist = append(cumdist, curdist)
		} else {
			cumdist = append(cumdist, curdist)
		}
	}
	return cumdist, nil
}

// UpdateLapNumbers returns workout data with lap number updated
func UpdateLapNumbers(strokes []StrokeRecord) ([]StrokeRecord, error) {
	lapnumber := 0
	currentstate := "rest"

	for _, stroke := range strokes {
		if currentstate == "rest" {
			switch stroke.Workoutstate {
			case
				5:
				{
					currentstate = "work"
					lapnumber += 1
				}
			}
		}
	}

	return strokes, nil
}

func divmodfloat(numerator float64, denominator int64) (quotient int64, remainder float64) {
	quotient = int64(numerator) / denominator // integer division, decimals are truncated
	remainder = numerator - float64(quotient)*float64(denominator)
	return
}

func formatTime(duration float64) string {
	min, sec := divmodfloat(duration, 60)
	if min > 60 {
		hour := min / 60
		min = min % 60
		str1 := fmt.Sprintf("%02d:%02d:%4.1f", hour, min, sec)
		return str1
	}
	str1 := fmt.Sprintf("00:%02d:%4.1f", min, sec)
	return str1

}

func formatPace(seconds float64) string {
	min, sec := divmodfloat(seconds, 60)
	return fmt.Sprintf("%02d:%04.1f", min, sec)
}

func workstring(
	totaldist, totaltime, avgpace, avgspm, avghr, maxhr, avgdps, avgpower float64,
	separator, symbol string,
) (string, error) {
	stri1 := symbol
	stri1 += fmt.Sprintf("%s%s%s%05.0f%s", separator, formatTime(totaltime), separator, totaldist, separator)
	stri1 += fmt.Sprintf("%s%s", formatPace(avgpace), separator)
	stri1 += fmt.Sprintf("%05.1f%s%4.1f%s%05.1f%s", avgpower, separator, avgspm, separator,
		avghr, separator)

	stri1 += fmt.Sprintf("%05.1f%s%04.1f", maxhr, separator, avgdps)

	return stri1, nil
}

// SummaryString writes summary
func SummaryString(strokes []StrokeRecord, title string, separator string) (string, error) {
	stri1 := fmt.Sprintf("Workout Summary - %s\n", title)
	stri1 += fmt.Sprintf("--%[1]sTotal%[1]s-Total----%[1]s--Avg--%[1]s-Avg-%[1]sAvg-%[1]s-Avg-%[1]s-Max-%[1]s-Avg\n", separator)
	stri1 += fmt.Sprintf("--%[1]sDist-%[1]s-Time-----%[1]s--Pace-%[1]s-Pwr-%[1]sSPM-%[1]s-HR--%[1]s-HR--%[1]s-DPS\n", separator)

	var avgpace, avgspm, avghr, maxhr, avgdps, avgpower float64
	var avgv float64

	cumdist, _ := cumulativedistance(strokes)
	totaldist := cumdist[len(strokes)-1]
	totaltime := strokes[len(strokes)-1].Timestamp - strokes[0].Timestamp

	for _, stroke := range strokes {
		avgv += 500. / stroke.Pace
		avgspm += stroke.Spm
		avghr += stroke.Hr
		avgdps += stroke.Strokedistance
		avgpower += stroke.Power
		if stroke.Hr > maxhr {
			maxhr = stroke.Hr
		}

	}

	avgv /= float64(len(strokes))
	avgpace = 500. / avgv
	avgspm /= float64(len(strokes))
	avghr /= float64(len(strokes))
	avgdps /= float64(len(strokes))
	avgpower /= float64(len(strokes))

	stri1 += fmt.Sprintf("--%s%05.0f%[1]s", separator, totaldist)
	stri1 += fmt.Sprintf("%s%s%s", formatTime(totaltime), separator, formatPace(avgpace))
	stri1 += fmt.Sprintf("%s%05.1f", separator, avgpower)
	stri1 += fmt.Sprintf("%[1]s%2.1f%[1]s%05.1f", separator, avgspm, avghr)
	stri1 += fmt.Sprintf("%s%05.1f%s%04.1f", separator, maxhr, separator, avgdps)

	return stri1, nil
}

// UpdateIntervalMetric updates intervals per metric
func UpdateIntervalMetric(
	strokes []StrokeRecord,
	metric string,
	setvalue float64,
	mode string,
	smoothwindow float64,
	activewindow [2]float64,
) ([]StrokeRecord, error) {
	// Set Active Window
	if activewindow[1] == 0 {
		activewindow[1] = strokes[len(strokes)-1].Timestamp

	}
	activewindow[0] -= strokes[0].Timestamp
	activewindow[1] -= strokes[0].Timestamp

	var metricvalues []float64
	var dtavg float64

	/// First loop to get some stats
	for i, stroke := range strokes {
		value, err := stroke.GetField(metric)
		if err != nil {
			return strokes, fmt.Errorf("Could not get value for %s in record number %v", metric, i)
		}
		metricvalues = append(metricvalues, value)
		if i >= 1 {
			dtavg += strokes[i].Timestamp - strokes[i-1].Timestamp
		}
	}
	dtavg /= float64(len(strokes))

	nrrecords := uint(smoothwindow / dtavg)
	metricvalues, _ = ewmovingaverage(metricvalues, nrrecords)

	// Second loop to set the values

	largerthantype := 5.
	smallerthantype := 3.
	if mode == "smaller" {
		largerthantype = 3.
		smallerthantype = 5.
	}

	timezero := strokes[0].Timestamp

	for i, stroke := range strokes {
		stroketime := stroke.Timestamp - timezero
		if metricvalues[i] >= setvalue && stroketime >= activewindow[0] && stroketime <= activewindow[1] {
			strokes[i].Workoutstate = largerthantype
		} else {
			strokes[i].Workoutstate = smallerthantype
		}
	}

	return strokes, nil
}

// AverageHR calculates average heart rate
func AverageHR(strokes []StrokeRecord) float64 {
	hr := 0.0
	var counter int
	for _, stroke := range strokes {
		if !math.IsNaN(stroke.Hr) {
			hr += stroke.Hr
		} else {
			counter++
		}
	}
	return hr / float64(len(strokes)-counter)
}

func geodistance(
	lat1 float64,
	lon1 float64,
	lat2 float64,
	lon2 float64) (float64, float64) {

	/* Approximate distance and bearing between two points
	   defined by lat1,lon1 and lat2,lon2
	   This is a slight underestimate but is close enough for our purposes,
	   We're never moving more than 10 meters between trackpoints

	   Bearing calculation fails if one of the points is a pole.

	*/

	// radius of earth in km
	const R = 6373.0

	// pi
	const pi = math.Pi

	lat1 = lat1 * pi / 180
	lat2 = lat2 * pi / 180
	lon1 = lon1 * pi / 180
	lon2 = lon2 * pi / 180

	dlon := lon2 - lon1
	dlat := lat2 - lat1

	a := sine(dlat/2)*sine(dlat/2) + math.Cos(lat1)*math.Cos(lat2)*sine(dlon/2)*sine(dlon/2)
	c := 2 * math.Atan2(math.Sqrt(a), math.Sqrt(1-a))

	distance := R * c

	x := sine(lon2-lon1) * math.Cos(lat2)
	y := math.Cos(lat1)*sine(lat2) - sine(lat1)*math.Cos(lat2)*math.Cos(lon2-lon1)

	tc1 := math.Atan2(x, y)

	if tc1 < 0 {
		tc1 += 2 * pi
	}

	tc1 = math.Mod(tc1, 2*pi)

	bearing := tc1 * 180 / pi

	return distance, bearing
}

// AverageSPM calculates average SPM
func AverageSPM(strokes []StrokeRecord) float64 {
	spm := 0.0
	var counter int
	for _, stroke := range strokes {
		if !math.IsNaN(stroke.Spm) {
			spm += stroke.Spm
		} else {
			counter++
		}
	}
	return spm / float64(len(strokes)-counter)
}

// AddBearing returns a stroke set with bearing
func AddBearing(strokes []StrokeRecord) {

	unfilteredbearing := make([]float64, len(strokes))

	for i := 0; i < len(strokes)-1; i++ {
		long1 := strokes[i].Longitude
		lat1 := strokes[i].Latitude
		long2 := strokes[i+1].Longitude
		lat2 := strokes[i+1].Latitude

		_, bearing := geodistance(lat1, long1, lat2, long2)

		unfilteredbearing[i] = bearing

	}

	filteredbearing, _ := ewmovingaverageboth(unfilteredbearing, 20)

	for i := range strokes {
		strokes[i].Bearing = filteredbearing[i]
	}

}

// AddStream adds river stream
func AddStream(strokes []StrokeRecord, vstream float64, unit string) {
	// foot / second
	if unit == "f" {
		vstream *= 0.3048
	}
	// knots
	if unit == "k" {
		vstream /= 1.994
	}
	// pace difference equivalent (approx)
	if unit == "p" {
		vstream *= 8 / 500
	}
	// now vstream is in m/s
	for _, stroke := range strokes {
		stroke.Vstream = vstream
	}
}

// AddWind adds wind speed and direction
func AddWind(strokes []StrokeRecord, vwind float64, winddirection float64, unit string) {
	// beaufort
	if unit == "b" {
		vwind = 0.837 * math.Sqrt(vwind*vwind*vwind)
	}
	// knots
	if unit == "k" {
		vwind *= 1.994
	}
	// km/h
	if unit == "kmh" {
		vwind /= 3.6
	}
	// mph
	if unit == "mph" {
		vwind *= 0.44704
	}
	for _, stroke := range strokes {
		stroke.Vwind = vwind
		stroke.Winddirection = winddirection
	}
}
