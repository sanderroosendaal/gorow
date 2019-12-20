package gorow

import (
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"testing"
)

const tolerance = 0.000001
const relativetolerance = 0.05

func GetTolerance(got float64, want float64) bool {
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
		if !GetTolerance(got[i], want[i]) {
			t.Errorf("Function %s, element %d, expected %f, got %f",
				name, i, want[i], got[i])
		}
	}
	return
}

func TestGetField(t *testing.T) {
	stroke := StrokeRecord{spm: 22}
	fmt.Println(stroke.spm)
	got, _ := stroke.GetField("spm")
	want := 22.0
	if math.Abs(got-want) > tolerance {
		t.Errorf("GetField gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}
}

func TestCSVReader(t *testing.T) {
	strokes, err := ReadCSV("testdata.csv")
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
	strokes, err := ReadCSV("otw.csv.gz")
	if err != nil {
		t.Errorf("CSVReader returned an error: %v", err.Error())
	}
	want := 990
	got := len(strokes)
	if want != got {
		t.Errorf("CSVReader got incorrect result. Got %d, wanted %d\n", got, want)
	}
}

func TestCSVReaderWriter(t *testing.T) {
	strokes, err := ReadCSV("testdata.csv")
	if err != nil {
		t.Errorf("CSVReader returned an error: %v", err.Error())
	}
	want := 191
	got := len(strokes)
	if want != got {
		t.Errorf("CSVReader got incorrect result. Got %d, wanted %d\n", got, want)
	}
	ok, err := WriteCSV(strokes, "out2.csv", true, false)
	if !ok {
		t.Errorf("CSVWriter: %v", err)
	}

	strokes2, err := ReadCSV("out2.csv")
	wantf := AveragePower(strokes)
	gotf := AveragePower(strokes2)

	if math.Abs(gotf-wantf) > tolerance {
		t.Errorf("WriteCSV equation gave incorrect result writing Power. Got %f, wanted %f\n",
			gotf, wantf)
	}

	err = os.Remove("out2.csv.gz")
	ok, err = WriteCSV(strokes, "out2.csv.gz", true, true)
	if !ok {
		t.Errorf("CSVWriter (gzip): %v", err)
	}

	strokes, err = ReadCSV("out2.csv.gz")
	wantf = AveragePower(strokes)
	gotf = AveragePower(strokes2)

	if math.Abs(gotf-wantf) > tolerance {
		t.Errorf("WriteCSV with gzip equation gave incorrect result writing Power. Got %f, wanted %f\n",
			gotf, wantf)
	}
}

func TestOTWSetPower(t *testing.T) {

	// 1x
	var rg = NewRig(0.9, 14, 2.655, 1.6, 0.88, "scull", -0.93, 0.0822, 0.46, 1, 0.98)
	var c = NewCrew(80, 1.4, 30, 0.5, SinusRecovery{}, Trapezium{X1: 0.15, X2: 0.5, H1: 1.0, H2: 0.9}, 1000., 1000.)

	strokes, err := ReadCSV("otw.csv")
	if err != nil {
		t.Errorf("CSVReader returned an error: %v", err.Error())
	}
	strokes = strokes[100:120]
	AddBearing(strokes)
	AddStream(strokes, 0, "f")
	AddWind(strokes, 0, 0, "m")
	fmt.Printf("Before: %.2f, %.2f, %.2f \n", AveragePower(strokes), AverageSPM(strokes), AverageHR(strokes))
	OTWSetPower(
		strokes, c, rg, "maherio",
		"http://localhost:8000/rowers/record-progress/testprogress/",
		false, true,
	)
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
		if !GetTolerance(got[i], want[i]) {
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
		0.08950832793934493,
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
		-0.0008034566368703589,
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
