package gorow

import (
	"math"
	"testing"
)

const tolerance = 0.000001

func TestDragEq(t *testing.T) {
	var got = dragEq(100, 4.5, 3.5, 0, 0)
	var want = 74.445191
	if math.Abs(got-want) > tolerance {
		t.Errorf("Drag equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = dragEq(100, 4.5, 3.5, 1, 1)
	want = 70.875
	if math.Abs(got-want) > tolerance {
		t.Errorf("Drag equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

	got = dragEq(100, 4.5, 0, 0, 0)
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
		t.Errorf("dRecovery equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}
}

func TestRig(t *testing.T) {
	var rg = newRig(0.89, 14., 2.89, 1.61, 0.88, Scull, -0.93, 822.e-4, 0.46, 1, 1.0)
	var got = rg.buitenhand()
	var want = 0.545856

	if math.Abs(got-want) > tolerance {
		t.Errorf("buitenhand equation gave incorrect result. Got %f, wanted %f\n",
			got, want)
	}

}
