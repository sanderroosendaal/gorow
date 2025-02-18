package gorow

import (
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"testing"
)

const tolerance = 0.000001
const relativetolerance = 0.05

func GetTolerance(got float64, want float64, relativetolerance float64) bool {
	if want == 0.0 {
		if got == want {
			return true
		}
		return false

	}

	if math.IsNaN(want) && math.IsNaN(got) {
		return true
	}

	var diff = math.Abs((got - want) / (want))
	return diff < relativetolerance
}

func ToleranceTest(t *testing.T, got []float64, want []float64, name string) {
	if len(got) != len(want) {
		t.Errorf("Function %s did not return the expected slice length. Got %d, wanted %d",
			name, len(got), len(want))
		return
	}

	for i := range want {
		if !GetTolerance(got[i], want[i], relativetolerance) {
			t.Errorf("Function %s, element %d, expected %f, got %f",
				name, i, want[i], got[i])
		}
	}
	return
}

func TestSine(t *testing.T) {
	x := math.Pi * 0.1
	got := sine(x)
	want := math.Sin(x)
	if math.Abs(got-want) > tolerance {
		t.Errorf("sine approximation is wrong")
	}
}

func TestRebase(t *testing.T) {
	x := []float64{1,3,5}
	y := []float64{1,9,25}
	x2 := []float64{1,2,3,4,5}
	expected := []float64{1,5,9,17,25}
	result, err := rebase(x,y,x2)
	if err != nil {
		t.Errorf("TestRebase, got error %s", err)
	}
	for i:= 0; i<len(expected)-1; i++ {
		if result[i] != expected[i] {
			t.Errorf("TestRebase, %d, got %f, want %f\n",
			i, result[i], expected[i])
		}
	}
}

func TestCosine(t *testing.T) {
	x := math.Pi * 0.1
	got := cosine(x)
	want := math.Cos(x)
	if math.Abs(got-want) > tolerance {
		t.Errorf("Cosine approximation is wrong")
	}
}

func TestGetField(t *testing.T) {
	stroke := StrokeRecord{Spm: 22}
	got, _ := stroke.GetField("Spm")
	want := 22.0
	if math.Abs(got-want) > tolerance {
		t.Errorf("GetField gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}
}

func TestCSVReader(t *testing.T) {
	strokes, err := ReadCSV("testdata/testdata.csv")
	if err != nil {
		t.Errorf("CSVReader returned an error: %v", err.Error())
	}
	want := 191
	got := len(strokes)
	if want != got {
		t.Errorf("CSVReader got incorrect result. Got %d, wanted %d\n", got, want)
	}
}

func TestCSVReaderGZip(t *testing.T) {
	strokes, err := ReadCSV("testdata/otw.csv.gz")
	if err != nil {
		t.Errorf("CSVReader returned an error: %v", err.Error())
	}
	want := 990
	got := len(strokes)
	if want != got {
		t.Errorf("CSVReader got incorrect result. Got %d, wanted %d\n", got, want)
	}
}

func TestParquetReaderWriter(t *testing.T) {

	strokes, err := ReadCSV("testdata/testdata.csv")
	_, err = WriteParquet(strokes, "testdata/testdata.parquet", true, true)
	if err != nil {
		t.Errorf("WriteParquet error: %v", err.Error())
	}

	strokes2, err := ReadParquet("testdata/testdata.parquet")
	if err != nil {
		t.Errorf("ReadParquet error: %v", err.Error())
	}

	wantf := AveragePower(strokes)
	gotf := AveragePower(strokes2)

	if math.Abs(gotf-wantf) > tolerance {
		t.Errorf("WriteCSV equation gave incorrect result writing Power. Got %f, wanted %f\n",
			gotf, wantf)
	}

	strokes3, err := ReadParquet("testdata/strokedata_2.parquet.gz/part.0.parquet")
	if err != nil {
		t.Errorf("ReadParquet error: %v", err.Error())
	}

	wantf = 0.0
	gotf = AveragePower(strokes3)

	if math.Abs(gotf-wantf) > tolerance {
		t.Errorf("WriteCSV equation gave incorrect result writing Power. Got %f, wanted %f\n",
			gotf, wantf)
	}

}

func TestCSVReaderWriter(t *testing.T) {
	runtime.GOMAXPROCS(1)
	strokes, err := ReadCSV("testdata/testdata.csv")
	if err != nil {
		t.Errorf("CSVReader returned an error: %v", err.Error())
	}
	want := 191
	got := len(strokes)
	if want != got {
		t.Errorf("CSVReader got incorrect result. Got %d, wanted %d\n", got, want)
	}
	ok, err := WriteCSV(strokes, "testdata/out2.csv", true, false)
	if !ok {
		t.Errorf("CSVWriter: %v", err)
	}

	strokes2, err := ReadCSV("testdata/out2.csv")
	wantf := AveragePower(strokes)
	gotf := AveragePower(strokes2)

	if math.Abs(gotf-wantf) > tolerance {
		t.Errorf("WriteCSV equation gave incorrect result writing Power. Got %f, wanted %f\n",
			gotf, wantf)
	}

	err = os.Remove("testdata/out2.csv.gz")
	ok, err = WriteCSV(strokes, "testdata/out2.csv.gz", true, true)
	if !ok {
		t.Errorf("CSVWriter (gzip): %v", err)
	}

	strokes, err = ReadCSV("testdata/out2.csv.gz")
	wantf = AveragePower(strokes)
	gotf = AveragePower(strokes2)

	if math.Abs(gotf-wantf) > tolerance {
		t.Errorf("WriteCSV with gzip equation gave incorrect result writing Power. Got %f, wanted %f\n",
			gotf, wantf)
	}
}

func TestReverse(t *testing.T) {
	in := []float64{1, 2, 3, 4}
	want := []float64{4, 3, 2, 1}
	got := reverseslice(in)
	for i := range got {
		if !GetTolerance(got[i], want[i], 0.0001) {
			t.Errorf("Reverse %d, got %v, wanted %v", i, got[i], want[i])
		}
	}
}

func TestEWMA(t *testing.T) {
	want := []float64{
		1.000000,
		1.666667,
		2.428571,
		3.266667,
		4.266667,
	}

	in := []float64{
		1, 2, 3, 4, 5,
	}

	got, _ := ewmovingaverage(in, 3)

	for i, value := range got {
		if !GetTolerance(value, want[i], 0.001) {
			t.Errorf("EWMA item %d, got %f, wanted %f", i, value, want[i])
		}
	}

}

func TestEWMARight(t *testing.T) {
	want := []float64{
		1.733333,
		2.571429,
		3.333333,
		4.0,
	}

	in := []float64{
		1, 2, 3, 4,
	}

	got, _ := ewmovingaverageright(in, 3)

	for i, value := range got {
		if !GetTolerance(value, want[i], 0.001) {
			t.Errorf("EWMA item %d, got %f, wanted %f", i, value, want[i])
		}
	}

}

func TestEWMABoth(t *testing.T) {
	want := []float64{
		1.366667,
		2.2,
		3,
		3.8,
		4.633333,
	}

	in := []float64{
		1, 2, 3, 4, 5,
	}

	got, _ := ewmovingaverageboth(in, 3)

	for i, value := range got {
		if !GetTolerance(value, want[i], 0.001) {
			t.Errorf("EWMA item %d, got %f, wanted %f", i, value, want[i])
		}
	}

}

func TestMetricsShort(t *testing.T) {
	tss, normp, trimp, hrtss, normv, normw, spmtss, err := WorkoutMetrics(
		"testdata/test_short.csv",
		200.0,
		"male",
		167, 185, 54, 500,
	)

	got := []float64{tss, normp, trimp, hrtss, normv, normw, spmtss}
	want := []float64{36.7, 171.3, 2.6, 1.49, 2.56, 171.3, 0}

	for i, value := range got {
		if !GetTolerance(value, want[i], relativetolerance) {
			t.Errorf("Function WorkoutMetrics short file, %d, got %f, wanted %f",
				i, got[i], want[i])
		}
	}

	if err != nil {
		t.Error("Function WorkoutMetrics raised an error for short file")
	}

}

func TestMetricsShort2(t *testing.T) {
	tss, normp, trimp, hrtss, normv, normw, spmtss, err := WorkoutMetrics(
		"testdata/test_short2.csv",
		200.0,
		"male",
		167, 185, 54, 500,
	)

	got := []float64{tss, normp, trimp, hrtss, normv, normw, spmtss}
	want := []float64{0.8, 118.3, 1.71, 1.0, 3.40, 360, 1.56}

	for i, value := range got {
		if !GetTolerance(value, want[i], relativetolerance) {
			t.Errorf("Function WorkoutMetrics short file, %d, got %f, wanted %f",
				i, got[i], want[i])
		}
	}

	if err != nil {
		t.Error("Function WorkoutMetrics raised an error for short file")
	}

}


func TestMetricsShortNK(t *testing.T) {
	tss, normp, trimp, hrtss, normv, normw, spmtss, err := WorkoutMetrics(
		"testdata/nk_testdata.csv",
		200.0,
		"male",
		167, 185, 54, 500,
	)

	got := []float64{tss, normp, trimp, hrtss, normv, normw, spmtss}
	want := []float64{52.8, 152, -6.5, -3.75, 3.40, 465, 64.3}

	for i, value := range got {
		if !GetTolerance(value, want[i], relativetolerance) {
			t.Errorf("Function WorkoutMetrics NK file, %d, got %f, wanted %f",
				i, got[i], want[i])
		}
	}

	if err != nil {
		t.Error("Function WorkoutMetrics raised an error for NK file")
	}

}

func TestMetrics(t *testing.T) {
	tss, normp, trimp, hrtss, normv, normw, spmtss, err := WorkoutMetrics(
		"testdata/testdata.csv",
		200.0,
		"male",
		167, 185, 54, 500,
	)

	if err != nil {
		t.Error("Function WorkoutMetrics gave an error")
	}

	got := []float64{tss, normp, trimp, hrtss, normv, normw, spmtss}
	want := []float64{8.120, 147.529, 16.782, 9.670, 3.722, 414.346, 11.5}

	for i, value := range got {
		if !GetTolerance(value, want[i], relativetolerance) {
			t.Errorf("Function WorkoutMetrics, %d, got %f, wanted %f", i, got[i], want[i])
		}
	}
}


func TestTailWind(t *testing.T) {
	want := 0.10260604299770072
	got := tailwind(340, 0.3, 230, 0)
	if !GetTolerance(want, got, 0.001) {
		t.Errorf("Function Tailwind, got %f, wanted %f", got, want)
	}
}

func TestGeoDistance(t *testing.T) {
	lat1 := 49.2387198
	lon1 := 16.5140534
	lat2 := 49.2387182
	lon2 := 16.514050199999996

	_, bearing := geodistance(lat1, lon1, lat2, lon2)

	wantbearing := 232.55

	if !GetTolerance(bearing, wantbearing, 0.001) {
		t.Errorf("Function geodistance gave bearing %f, wanted %f", bearing, wantbearing)
	}

}

func TestOTWSetPower(t *testing.T) {

	// 1x
	var rg = NewRig(0.9, 14, 2.655, 1.6, 0.88, "scull", -0.93, 0.0822, 0.46, 1, 0.98)
	var c = NewCrew(80, 1.4, 30, 0.5, SinusRecovery{}, Trapezium{X1: 0.15, X2: 0.5, H1: 1.0, H2: 0.9}, 1000., 1000.)

	strokes, err := ReadCSV("testdata/otw.csv")
	if err != nil {
		t.Errorf("CSVReader returned an error: %v", err.Error())
	}
	strokes = strokes[100:120]
	AddBearing(strokes)
	/*
		wantslice := []float64{
			160.677160,
			135.033001,
			94.236924,
			79.676226,
			95.165019,
			87.488895,
		}

		for i, want := range wantslice {
			if !GetTolerance(strokes[i].bearing, want, 0.0001) {
				t.Errorf("Bearing %d, got %f, wanted %f", i, strokes[i].bearing, want)
			}
		}
	*/
	AddStream(strokes, 0, "f")
	AddWind(strokes, 0, 0, "m")
	fmt.Printf("Before: %.2f, %.2f, %.2f \n", AveragePower(strokes), AverageSPM(strokes), AverageHR(strokes))
	OTWSetPower(
		strokes, c, rg, "maherio",
		"http://localhost:8000/rowers/record-progress/testprogress/",
		false, true,
	)
	if err != nil {
		fmt.Printf("smoothnowindpace gave error: %v", err.Error())
	}
	fmt.Printf("After: %.2f, %.2f, %.2f \n", AveragePower(strokes), AverageSPM(strokes), AverageHR(strokes))
}

func TestInterPol3(t *testing.T) {
	x := []float64{1, 2, 4, 5, 6}
	y := []float64{1, 2, 4, 5, 6}

	want := 3.0
	got, _ := srinterpol3(x, y, 3)

	if math.Abs(got-want) > tolerance {
		t.Errorf("Interpol3 equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

}

func TestSlices(t *testing.T) {
	s1 := []float64{1, 2}
	s2 := []float64{3, 4}

	got := slicesadd(s1, s2)
	want := []float64{4, 6}

	for i := range want {
		if got[i] != want[i] {
			t.Errorf("Slices addition failed. Got %f, wanted %f\n", got[i], want[i])
		}
	}

	got = slicesadd(s1, s2, s2)
	want = []float64{7, 10}

	for i := range want {
		if got[i] != want[i] {
			t.Errorf("Slices addition failed. Got %f, wanted %f\n", got[i], want[i])
		}
	}

	v := []float64{2, 1, 3, 4, 5, 123, 1}
	smallest, biggest := sliceminmax(v)
	wantbig := 123.0
	wantsmall := 1.0

	if smallest != wantsmall {
		t.Errorf("Slices min failed. Got %f, wanted %f\n", smallest, wantsmall)
	}

	if biggest != wantbig {
		t.Errorf("Slices min failed. Got %f, wanted %f\n", biggest, wantbig)
	}

}

func TestDragEq(t *testing.T) {
	var got = DragEq(100, 4.5, 3.5, 0, 0)
	var want = 74.445191
	if math.Abs(got-want) > tolerance {
		t.Errorf("Drag equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = DragEq(100, 4.5, 3.5, 0, 1)
	want = 70.875
	if math.Abs(got-want) > tolerance {
		t.Errorf("Drag equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = DragEq(100, 4.5, 0, 0, 0)
	want = 74.445191
	if math.Abs(got-want) > tolerance {
		t.Errorf("Drag equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = DragEq(110, 4.2, 3.5, 0, 0)
	want = 69.51463720445798
	if math.Abs(got-want) > tolerance {
		t.Errorf("Drag equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

}

func TestDRecovery(t *testing.T) {
	var got = dRecovery(0.01, 4.5, 1.0, 0.1, 80., 14., 3.5, 100.)
	var want = 0.005681
	if math.Abs(got-want) > tolerance {
		t.Errorf("dRecovery equation gave incorrect result. Got %f, wanted %f.\n",
			got, want)
	}
}

func TestCrewExportImport(t *testing.T) {
	var c = NewCrew(80, 1.4, 30, 0.5, SinusRecovery{}, Trapezium{X1: 0.15, X2: 0.5, H1: 1.0, H2: 0.9}, 1000., 1000.)
	s, err := c.ToJSON()
	if err != nil {
		t.Errorf("rigging ToJSON yielded an error, %v", err.Error())
	}

	want := 156
	got := len(s)

	if want != got {
		t.Errorf("Crew ToJSON error. String length got %d, want %d", got, want)
	}

	var c2 Crew

	c2.FromJSON(s)

	wantf := c.Strokelength
	gotf := c2.Strokelength

	if math.Abs(gotf-wantf) > tolerance {
		t.Errorf("Rigging FromJSON equation gave incorrect result. Got %f, wanted %f\n",
			gotf, wantf)
	}

}

func TestCrew(t *testing.T) {
	var c = NewCrew(80, 1.4, 30, 0.5, SinusRecovery{}, Trapezium{X1: 0.15, X2: 0.5, H1: 1.0, H2: 0.9}, 1000., 1000.)
	var got = c.Strokelength
	var want = 1.4
	if math.Abs(got-want) > tolerance {
		t.Errorf("Crew stroke length incorrect. Got %f, wanted %f\n", got, want)
	}

	want = 1.2261495151704724
	got = c.vha(1, 0.3)

	if math.Abs(got-want) > tolerance {
		t.Errorf("Crew vha. Got %f, wanted %f\n", got, want)
	}

	want = 0.8155612244897958
	got = c.vcm(1, 0.3)
	if math.Abs(got-want) > tolerance {
		t.Errorf("Crew vcm. Got %f, wanted %f\n", got, want)
	}

	want = 0.8346938775510204
	got = c.vcm(1, 0.2)
	if math.Abs(got-want) > tolerance {
		t.Errorf("Crew vcm. Got %f, wanted %f\n", got, want)
	}

	want = 0.5400510204081632
	got = c.vcm(1, 0.9)
	if math.Abs(got-want) > tolerance {
		t.Errorf("Crew vcm. Got %f, wanted %f\n", got, want)
	}

	want = 0.2700255102040816
	got = c.vcm(0.5, 0.9)
	if math.Abs(got-want) > tolerance {
		t.Errorf("Crew vcm. Got %f, wanted %f\n", got, want)
	}

	var wantv = []float64{0.83469388, 0.91816327, 1.00163265}
	var vhandle = []float64{1, 1.1, 1.2}
	var xhandle = []float64{0.2, 0.2, 0.2}
	var gotv = c.vcma(vhandle, xhandle)

	for i := range wantv {
		if math.Abs(gotv[i]-wantv[i]) > tolerance {
			t.Errorf("Crew vcma, got %f, wanted %f", gotv[i], wantv[i])
		}
	}

	want = -0.7853981633974483
	got = c.vhandle(0.5, 2.0, 1.0)
	if math.Abs(got-want) > tolerance {
		t.Errorf("Crew sinus recovery error, got %f, wanted %f", got, want)
	}

	want = -0.5
	got = c.dxhandle(0.5, 2.0, 1.0)
	if math.Abs(got-want) > tolerance {
		t.Errorf("Crew flat recovery error, got %f, wanted %f", got, want)
	}
}

func TestRigZero(t *testing.T) {
	var rg = NewRig(0, 0, 0, 0, 0, Scull, 0, 0, 0, 0, 0)

	var got = rg.buitenhand()
	var want = 0.523899

	if math.Abs(got-want) > tolerance {
		t.Errorf("buitenhand equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}
}

func TestRigExportImport(t *testing.T) {
	var rg = NewRig(0.89, 14., 2.89, 1.61, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)
	s, err := rg.ToJSON()
	if err != nil {
		t.Errorf("rigging ToJSON yielded an error, %v", err.Error())
	}

	want := 130
	got := len(s)

	if want != got {
		t.Errorf("rigging ToJSON error. String length got %d, want %d", got, want)
	}

	var rg2 Rig

	rg2.FromJSON(s)

	wantf := rg.buitenhand()
	gotf := rg2.buitenhand()

	if math.Abs(gotf-wantf) > tolerance {
		t.Errorf("Rigging FromJSON equation gave incorrect result. Got %f, wanted %f\n",
			gotf, wantf)
	}
}

func TestRigExportImportFiles(t *testing.T) {
	var rg = NewRig(0.89, 14., 2.89, 1.61, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)
	s, err := rg.ToJSON()
	if err != nil {
		t.Errorf("rigging ToJSON yielded an error, %v", err.Error())
	}

	temppath := filepath.FromSlash("/tmp/dat1")
	fmt.Printf("Temporary file location: %v\n", temppath)

	err = ioutil.WriteFile(temppath, []byte(s), 0644)

	var rg2 Rig

	content, err := ioutil.ReadFile(temppath)
	if err != nil {
		t.Errorf("Rigging reading from file yielded an error")
	}

	err = rg2.FromJSON(string(content))
	if err != nil {
		t.Errorf("Rigging fromJSON an error")
	}

	wantf := rg.buitenhand()
	gotf := rg2.buitenhand()

	if math.Abs(gotf-wantf) > tolerance {
		t.Errorf("Rigging FromJSON equation gave incorrect result. Got %f, wanted %f\n",
			gotf, wantf)
	}
}

func TestRig(t *testing.T) {
	var rg = NewRig(0.89, 14., 2.89, 1.61, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)

	var got = rg.buitenhand()
	var want = 0.545856

	if math.Abs(got-want) > tolerance {
		t.Errorf("buitenhand equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = rg.bladearea()
	want = 0.1644

	if math.Abs(got-want) > tolerance {
		t.Errorf("bladearea equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = rg.dcatch()
	want = -0.713442

	if math.Abs(got-want) > tolerance {
		t.Errorf("dcatch equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = rg.oarangle(0.2)
	want = -0.614929

	if math.Abs(got-want) > tolerance {
		t.Errorf("oarangle equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = rg.spread()
	want = 0.805

	if math.Abs(got-want) > tolerance {
		t.Errorf("oarangle equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = rg.overlap()
	want = 0.17

	if math.Abs(got-want) > tolerance {
		t.Errorf("oarangle equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

}

func TestDStroke(t *testing.T) {
	var got = dStroke(0.01, 3.5, 3.4, 0.1, 80, 14, 3.5, 100)
	var want = 0.010501

	if math.Abs(got-want) > tolerance {
		t.Errorf("dStroke equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}
}

func TestBasics(t *testing.T) {
	var rg = NewRig(0.89, 14., 2.89, 1.61, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)
	var mc = 80.0
	var mb = 14.0

	var v = 3.9
	var lout = rg.Lscull - rg.Lin
	var got = vhandle(v, rg.Lin, lout, mc, mb)
	var want = 1.373323

	if math.Abs(got-want) > tolerance {
		t.Errorf("dStroke equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}
}

func TestDFootboard(t *testing.T) {
	var mc = 80.
	var mb = 14.0
	var vs1 = 3.9
	var vs2 = 4.0

	var got = deFootboard(mc, mb, vs1, vs2)
	var want = 4.706383

	if math.Abs(got-want) > tolerance {
		t.Errorf("dStroke equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}
}

func TestStroke2Seconds(t *testing.T) {
	t0 := []float64{1,3,5}
	y0 := []float64{1,9,25}
	want := []float64{1,5,9,17,25}
	twant := []float64{1,2,3,4,5}
	tgot, ygot, err := stroke2second(t0, y0)
	if err != nil {
		t.Errorf("TestRebase, got error %s", err)
	}
	for i:= 0; i<len(want)-1; i++ {
		if ygot[i] != want[i] {
			t.Errorf("TestStroke2Seconds (y), %d, got %f, want %f\n",
			i, ygot[i], want[i])
		}
		if tgot[i] != twant[i] {
			t.Errorf("TestStroke2Seconds (x), %d, got %f, want %f\n",
			i, tgot[i], twant[i])
		}
	}

}

func TestLinSpace(t *testing.T) {
	want := []float64{1., 1.02040816, 1.04081633, 1.06122449, 1.08163265,
		1.10204082, 1.12244898, 1.14285714, 1.16326531, 1.18367347,
		1.20408163, 1.2244898, 1.24489796, 1.26530612, 1.28571429,
		1.30612245, 1.32653061, 1.34693878, 1.36734694, 1.3877551,
		1.40816327, 1.42857143, 1.44897959, 1.46938776, 1.48979592,
		1.51020408, 1.53061224, 1.55102041, 1.57142857, 1.59183673,
		1.6122449, 1.63265306, 1.65306122, 1.67346939, 1.69387755,
		1.71428571, 1.73469388, 1.75510204, 1.7755102, 1.79591837,
		1.81632653, 1.83673469, 1.85714286, 1.87755102, 1.89795918,
		1.91836735, 1.93877551, 1.95918367, 1.97959184, 2.}

	got, _ := LinSpace(1, 2, 50)

	if len(got) != len(want) {
		t.Errorf("Function LinSpace did not return the expected slice. Got %d, wanted %d",
			len(got), len(want))
	}

	for i := range got {
		if math.Abs(got[i]-want[i]) > tolerance {
			t.Errorf("Function Linspace %d, got %f, wanted %f", i, got[i], want[i])
		}
	}

	// check error
	got, err := LinSpace(1, 2, 0)
	if err == nil {
		t.Error("Linspace should give error when called with length of zero")
	}
	if got != nil {
		t.Error("Linspace should give error when called with length of zero")
	}

}

func TestBladeForce(t *testing.T) {
	var want = []float64{
		1.8192760229325606,
		99.98194431474082,
		82.51865949087183,
		-98.81844828718056,
		15.208664210566168,
		-0.30068789227801074,
		0.046277403309847434,
		-0.15270692145314926,
	}

	var rg = NewRig(0.9, 14, 2.885, 1.60, 0.88, Scull,
		-0.93, 822.e-4, 0.46,
		1, 1.0)

	got, _ := BladeForce(-0.6, rg, 3.5, 100)

	if len(got) != len(want) {
		t.Errorf("Function BladeForce did not return the expected slice. Got %d, wanted %d",
			len(got), len(want))
	}

	for i := 0; i < len(want); i++ {
		if !GetTolerance(got[i], want[i], relativetolerance) {
			t.Errorf("Function BladeForce, element %d, expected %f, got %f",
				i, want[i], got[i])
		}
	}

	want = []float64{2.0803788704765536,
		99.99842385483569,
		76.4830132261339,
		-99.73887777592651,
		7.200071767014857,
		-0.14362994318676583,
		0.010368533533736905,
		-0.07206421091268539}

	got, _ = BladeForce(-0.7, rg, 4.5, 100)

	ToleranceTest(t, got, want, "BladeForce")

}

func TestEnergyBalance(t *testing.T) {

	want := []float64{
		0.070523, // 0.08950832793934493,
		3.3695083279393447,
		3.700052225266561,
		0.5223880597014925,
		561.8516594477129,
		280.92582972385645,
		0.718216020128656,
		5.030428996890504,
		2.4749540307829343,
		0.0, // 3.8866101681872163,
		0.0, // 2.499200316746598,
		2.5554749661075697,
		0.0, // 33.70434561030705,
		0.0, // 0.5970149253731345,
		12.429634017382732,
		0.8613483156138445,
	}

	c := NewCrew(
		80., 1.4, 30.0, 0.5,
		SinusRecovery{},
		Trapezium{X1: 0.15, X2: 0.5, H2: 0.9, H1: 1.0}, 1000., 1000.)
	rg := NewRig(0.9, 14, 2.885, 1.60, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)
	got, _ := EnergyBalance(350, c, rg, 3.28, 0.03, 5.0, 0.0, true)

	ToleranceTest(t, got, want, "EnergyBalance")

}

func TestStroke(t *testing.T) {
	want := []float64{
		-0.001490, // -0.0008034566368703589,
		3.405032807198155,
		3.805137532355478,
		0.49328358208955214,
		553.541979346718,
		276.770989673359,
		0.7211424362408791,
		5.012794324778716,
		2.5904571157190004,
		0.0, // 5.507898476518529,
		0.0, // 2.4779480499384134,
		2.422337209059717,
		0.0, // 92.55933370813963,
		0.0, // 0.5402985074626863,
		11.370089516607615,
		0.8768551568021985,
	}

	c := NewCrew(
		80., 1.4, 30.0, 0.5,
		SinusRecovery{},
		Trapezium{X1: 0.15, X2: 0.5, H2: 0.9, H1: 1.0}, 1000., 1000.)
	rg := NewRig(0.9, 14, 2.885, 1.60, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)
	got, _ := Stroke(350, c, rg, 3.42, 0.03, 20, 5.0, true, 0.0)

	ToleranceTest(t, got, want, "Stroke")

}

func TestConstantVelo(t *testing.T) {
	want := []float64{
		249.34555283522275,
		3.4481771832103107,
		0.5462686567164179,
		199.53058042827024,
		0.7239294989576714,
	}
	c := NewCrew(
		80., 1.4, 30.0, 0.5,
		SinusRecovery{},
		Trapezium{X1: 0.15, X2: 0.5, H2: 0.9, H1: 1.0}, 1000., 1000.)
	rg := NewRig(0.9, 14, 2.885, 1.60, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)

	got, _ := ConstantVeloFast(3.42, c, rg, 0.03, 5, 5, 100, 400, 5, 0.0, true)

	ToleranceTest(t, got, want, "ConstantVeloFast")
}

func TestConstantWatt(t *testing.T) {
	want := []float64{
		250,
		3.44,
		0.54,
		199,
		0.723,
	}

	c := NewCrew(
		80., 1.4, 30.0, 0.5,
		SinusRecovery{},
		Trapezium{X1: 0.15, X2: 0.5, H2: 0.9, H1: 1.0}, 1000., 1000.)
	rg := NewRig(0.9, 14, 2.885, 1.60, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)

	got, _ := ConstantWattFast(200, c, rg, 0.03, 5, 5, 50, 1000, 5, 0, true, 15)
	ToleranceTest(t, got, want, "ConstantWattFast")

}

func TestMaas(t *testing.T) {
	// nkfactor := 1.15
	// pwr := []float64{114, 129, 151, 175, 222}
	pwr := []float64{113, 116, 131, 127, 137, 119, 133, 133, 151, 163, 177, 195, 224}
	// spm := []float64{18.8, 22.6, 24.9, 28.1, 31.1}
	spm := []float64{17.9, 19.8, 22.0, 22.7, 24.4, 24.4, 20.2, 22.2, 24.4, 26.3, 27.9, 29.5, 32.1}
	// want := []float64{2.6, 3.0, 3.3, 3.5, 3.6}
	want := []float64{2.89, 2.95, 3.05, 3.01, 3.05, 2.69, 3.1, 3.06, 3.28, 3.38, 3.52, 3.51, 3.73}

	rg := RigCoastalMaas
	rgolympic := RigSingle
	sum := 0.0
	for i := range want {
		c := NewCrew(
			93.9, 1.4, spm[i], 0.5,
			SinusRecovery{},
			Trapezium{X1: 0.15, X2: 0.5, H2: 0.9, H1: 1.0}, 1000., 1000.)
		got, _ := ConstantWattFast(pwr[i], c, rg, 0.03, 5, 5, 50, 1000, 5., 0, true, 15)
		got2, _ := ConstantWattFast(pwr[i], c, rgolympic, 0.03, 5, 5, 50, 1000, 5., 0, true, 15)
		gotvelo := got[1]
		gotveloolympic := got2[1]
		ratio := math.Pow(gotveloolympic/gotvelo, 3.0)
		sum += ratio
		if !GetTolerance(ratio, 1.24, 0.2) {
			t.Errorf("Run %v : Maas Aero/Olympic shell velo ratio got %v, wanted %v at power %v", i, ratio, 1.24, pwr[i])
		}

	}

	avg := sum / float64(len(pwr))
	println(avg)
}

func TestPhysGetPower(t *testing.T) {
	want := []float64{
		237.,
		0.522,
		323.,
		139.,
		math.NaN(),
	}

	c := NewCrew(
		80., 1.4, 30.0, 0.5,
		SinusRecovery{},
		Trapezium{X1: 0.15, X2: 0.5, H2: 0.9, H1: 1.0}, 1000., 1000.)
	rg := NewRig(0.9, 14, 2.885, 1.60, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)

	got, _ := PhysGetPower(3.5, c, rg, 90., 2.1, 160., 0)
	ToleranceTest(t, got, want, "ConstantWattFast")
}
