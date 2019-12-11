package gorow

import (
	"encoding/csv"
	"io"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
)

// LbstoN convert lbs of force to Newton
const LbstoN = 4.44822

// StrokeRecord sort of dataframe
type StrokeRecord struct {
	timestamp          float64
	distance           float64
	spm                float64
	hr                 int
	pace               float64
	power              float64
	drivelength        float64
	strokedistance     float64
	drivetime          float64
	dragfactor         int
	strokerecoverytime float64
	workperstroke      float64
	averageforce       float64
	peakforce          float64
	velo               float64
	lapnr              int
	intervaltime       float64
	calories           float64
	workoutstate       int
	latitude           float64
	longitude          float64
	bearing            float64
}

func getfloatrecord(s string) (float64, error) {
	return strconv.ParseFloat(strings.TrimSpace(s), 64)
}

func getintrecord(s string) (int, error) {
	return strconv.Atoi(strings.TrimSpace(s))
}

// ReadCSV reads rowing data into data frame
func ReadCSV(f string) []StrokeRecord {
	csvFile, _ := os.Open(f)
	defer csvFile.Close()
	reader := csv.NewReader(csvFile)
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
			log.Fatal(err)
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
					if f, err := getintrecord(record[i]); err == nil {
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

					//	   DragFactor
				}
			}
			if row.velo == 0 && row.pace != 0 {
				row.velo = 500. / row.pace
			}
			rows = append(rows, row)
		}
	}
	return rows
}

// OTWSetPower adds power for OTW rows
func OTWSetPower(strokes []StrokeRecord) {
	// temporary default values
	c := NewCrew(
		75., 1.4, 30.0, 0.5,
		SinusRecovery{},
		Trapezium{x1: 0.15, x2: 0.5, h2: 0.9, h1: 1.0}, 1000., 1000.)
	rg := NewRig(0.9, 14, 2.885, 1.60, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)

	// newstrokes := strokes

	for i, stroke := range strokes {
		c.tempo = stroke.spm
		res, err := PhysGetPower(stroke.velo, c, rg, 0, 0, 0, 0)
		if err != nil {
			break // ignore value
		}
		pwr := res[0]
		frc := res[2]
		strokes[i].power = pwr
		strokes[i].averageforce = frc / LbstoN
	}
}

// AveragePower calculates average power
func AveragePower(strokes []StrokeRecord) float64 {
	power := 0.0
	for _, stroke := range strokes {
		power += stroke.power
	}
	return power / float64(len(strokes))
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
	c := 2 * math.Atan(math.Sqrt(1-a)/math.Sqrt(a))

	distance := R * c

	x := math.Sin(lon2-lon1) * math.Cos(lat2)
	y := math.Cos(lat1)*math.Sin(lat2) - math.Sin(lat1)*math.Cos(lat2)*math.Cos(lon2-lon1)

	tc1 := math.Atan(y / x)

	tc1 = math.Mod(tc1, 2*pi)

	bearing := tc1 * 180 / pi

	return distance, bearing
}

// AverageSPM calculates average SPM
func AverageSPM(strokes []StrokeRecord) float64 {
	spm := 0.0
	for _, stroke := range strokes {
		spm += stroke.spm
	}
	return spm / float64(len(strokes))
}

// AddBearing returns a stroke set with bearing
func AddBearing(strokes []StrokeRecord) {
	for i := 0; i < len(strokes)-1; i++ {
		long1 := strokes[i].longitude
		lat1 := strokes[i].latitude
		long2 := strokes[i+1].longitude
		lat2 := strokes[i+1].latitude

		_, bearing := geodistance(lat1, long1, lat2, long2)
		strokes[i].bearing = bearing
	}
}
