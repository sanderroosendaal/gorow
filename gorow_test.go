package gorow

import (
	"math"
	"testing"
)

const tolerance = 0.000001
const relativetolerance = 0.05

func GetTolerance(got float64, want float64) bool {
	var diff = math.Abs((got - want) / (want))
	return diff < relativetolerance
}

func TestDragEq(t *testing.T) {
	var got = DragEq(100, 4.5, 3.5, 0, 0)
	var want = 74.445191
	if math.Abs(got-want) > tolerance {
		t.Errorf("Drag equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = DragEq(100, 4.5, 3.5, 1, 1)
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

}

func TestDRecovery(t *testing.T) {
	var got = dRecovery(0.01, 4.5, 1.0, 0.1, 80., 14., 3.5, 100.)
	var want = 0.005681
	if math.Abs(got-want) > tolerance {
		t.Errorf("dRecovery equation gave incorrect result. Got %f, wanted %f.\n",
			got, want)
	}
}

func TestCrew(t *testing.T) {
	var c = NewCrew(80, 1.4, 30, 0.5, 1000., 1000.)
	var got = c.strokelength
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

	var wantv = []float64{0.83469388, 0.91816327, 1.00163265}
	var vhandle = []float64{1, 1.1, 1.2}
	var xhandle = []float64{0.2, 0.2, 0.2}
	var gotv = c.vcma(vhandle, xhandle)

	for i := range wantv {
		if math.Abs(gotv[i]-wantv[i]) > tolerance {
			t.Errorf("Crew vcma, got %f, wanted %f", gotv[i], wantv[i])
		}
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
	var lout = rg.lscull - rg.lin
	var got = vhandle(v, rg.lin, lout, mc, mb)
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

	var got = BladeForce(-0.6, rg, 3.5, 100)

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
}
