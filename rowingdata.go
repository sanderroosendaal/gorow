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

	"github.com/VividCortex/ewma"
	"github.com/fatih/structs"
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
	timestamp          float64 `rowingdata:"TimeStamp (sec)"`
	distance           float64 `rowingdata:" Horizontal (meters)"`
	spm                float64 `rowingdata:" Cadence (stokes/min)"`
	hr                 float64 `rowingdata:" HRCur (bpm)"`
	pace               float64 `rowingdata:" Stroke500mPace (sec/500m)"`
	power              float64 `rowingdata:" Power (watts)"`
	drivelength        float64 `rowingdata:" DriveLength (meters)"`
	strokedistance     float64 `rowingdata:" StrokeDistance (meters)"`
	drivetime          float64 `rowingdata:" drivetime"`
	dragfactor         int     `rowingdata:" DragFactor"`
	strokerecoverytime float64 `rowingdata:" StrokeRecoveryTime (ms)"`
	workperstroke      float64 `rowingdata:" WorkPerStroke (joules)"`
	averageforce       float64 `rowingdata:" AverageDriveForce (lbs)"`
	peakforce          float64 `rowingdata:" PeakDriveForce (lbs)"`
	velo               float64 `rowingdata:" Speed (m/sec)"`
	lapnr              int     `rowingdata:" lapIdx"`
	intervaltime       float64 `rowingdata:" ElapsedTime (sec)"`
	calories           float64 `rowingdata:" Calories (kCal)"`
	workoutstate       int     `rowingdata:" WorkoutState"`
	latitude           float64 `rowingdata:" latitude"`
	longitude          float64 `rowingdata:" longitude"`
	bearing            float64 `rowingdata:" bearing"`
	nowindpace         float64 `rowingdata:"nowindpace"`
	equivergpower      float64 `rowingdata:"Equiv erg Power"`
	modelpower         float64 `rowingdata:"power (model)"`
	modelfavg          float64 `rowingdata:"averageforce (model)"`
	modeldrivelength   float64 `rowingdata:"drivelength (model)"`
	vwind              float64 `rowingdata:"vwind"`
	winddirection      float64 `rowingdata:"winddirection"`
	vstream            float64 `rowingdata:"vstream"`
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

func getintrecord(s string) (int, error) {
	return strconv.Atoi(strings.TrimSpace(s))
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
						row.timestamp = f
					}
				case " lapIdx":
					if f, err := getintrecord(record[i]); err == nil {
						row.lapnr = f
					}
				case " ElapsedTime (sec)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.intervaltime = f
					}
				case " Horizontal (meters)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.distance = f
					}
				case " Stroke500mPace (sec/500m)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.pace = f
					}
				case " Cadence (stokes/min)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.spm = f
					}
				case " Cadence (strokes/min)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.spm = f
					}
				case " HRCur (bpm)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.hr = f
					}
				case " Power (watts)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.power = f
					}
				case " Calories (kCal)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.calories = f
					}
				case " Speed (m/sec)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.velo = f
					}
				case " StrokeDistance (meters)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.strokedistance = f
					}
				case " DriveLength (meters)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.drivelength = f
					}
				case " StrokeRecoveryTime (ms)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.strokerecoverytime = f
					}
				case " WorkPerStroke (joules)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.workperstroke = f
					}
				case " AverageDriveForce (lbs)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.averageforce = f
					}
				case " AverageDriveForce (N)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.averageforce = f / LbstoN
					}
				case " PeakDriveForce (lbs)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.peakforce = f
					}
				case " PeakDriveForce (N)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.peakforce = f / LbstoN
					}
				case " DragFactor":
					if f, err := getintrecord(record[i]); err == nil {
						row.dragfactor = f
					}
				case " latitude":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.latitude = f
					}
				case " longitude":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.longitude = f
					}
				case "nowindpace":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.nowindpace = f
					}
				case "Equiv erg Power":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.equivergpower = f
					}
				case "power (model)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.modelpower = f
					}
				case "averageforce (model)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.modelfavg = f
					}
				case "drivelength (model)":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.modeldrivelength = f
					}
				case "vwind":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.vwind = f
					}
				case "winddirection":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.winddirection = f
					}
				case "vstream":
					if f, err := getfloatrecord(record[i]); err == nil {
						row.vstream = f
					}
				}
			}
			if row.velo == 0 && row.pace != 0 {
				row.velo = 500. / row.pace
			}
			if row.workperstroke == 0 && row.power != 0 && row.spm != 0 {
				row.workperstroke = 60. * row.power / row.spm
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

			c.Tempo = stroke.spm
			res, err := PhysGetPower(
				stroke.velo, c, rg, stroke.bearing,
				stroke.vwind, stroke.winddirection, stroke.vstream,
			)
			//res, err := PhysGetPower(stroke.velo, c, rg, 0, 0, 0, 0)
			if err == nil {

				pwr := res[0]
				frc := res[2]
				nowindp := res[3]
				if !powermeasured {
					strokes[i].power = pwr
					strokes[i].averageforce = frc / LbstoN
				}
				strokes[i].nowindpace = nowindp
				strokes[i].modelpower = pwr
				strokes[i].modelfavg = frc / LbstoN
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

func averagenowindpace(strokes []StrokeRecord) float64 {
	p := 0.0
	var counter int
	for _, stroke := range strokes {
		if !math.IsNaN(stroke.nowindpace) {
			p += stroke.nowindpace
		} else {
			counter++
		}
	}
	return p / float64(len(strokes)-counter)
}

// AveragePower calculates average power
func AveragePower(strokes []StrokeRecord) float64 {
	power := 0.0
	var counter int
	for _, stroke := range strokes {
		if !math.IsNaN(stroke.power) {
			power += stroke.power
		} else {
			counter++
		}
	}
	return power / float64(len(strokes)-counter)
}

// AverageHR calculates average heart rate
func AverageHR(strokes []StrokeRecord) float64 {
	hr := 0.0
	var counter int
	for _, stroke := range strokes {
		if !math.IsNaN(stroke.hr) {
			hr += stroke.hr
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

	a := math.Sin(dlat/2)*math.Sin(dlat/2) + math.Cos(lat1)*math.Cos(lat2)*math.Sin(dlon/2)*math.Sin(dlon/2)
	c := 2 * math.Atan2(math.Sqrt(a), math.Sqrt(1-a))

	distance := R * c

	x := math.Sin(lon2-lon1) * math.Cos(lat2)
	y := math.Cos(lat1)*math.Sin(lat2) - math.Sin(lat1)*math.Cos(lat2)*math.Cos(lon2-lon1)

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
		if !math.IsNaN(stroke.spm) {
			spm += stroke.spm
		} else {
			counter++
		}
	}
	return spm / float64(len(strokes)-counter)
}

// AddBearing returns a stroke set with bearing
func AddBearing(strokes []StrokeRecord) {

	e := ewma.NewMovingAverage(20)
	e1 := ewma.NewMovingAverage()

	for i := 0; i < len(strokes)-1; i++ {
		long1 := strokes[i].longitude
		lat1 := strokes[i].latitude
		long2 := strokes[i+1].longitude
		lat2 := strokes[i+1].latitude

		_, bearing := geodistance(lat1, long1, lat2, long2)
		e.Add(bearing)
		e1.Add(bearing)
		if i < 20 {
			strokes[i].bearing = e1.Value()
		} else {
			strokes[i].bearing = e.Value()
		}
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
		stroke.vstream = vstream
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
		stroke.vwind = vwind
		stroke.winddirection = winddirection
	}
}
