package gorow

import (
	"encoding/json"
	"math"
)

// Rig holds boat rigging parameters and has methods to manipulate them
type Rig struct {
	Lin         float64
	Mb          float64
	BladeLength float64
	Lscull      float64
	Nrowers     int32
	RowOrScull  string
	Span        float64
	CatchAngle  float64
	DragForm    float64
	_Spread     float64
	_Bladearea  float64
}

func (rg *Rig) spread() float64 {
	if rg.RowOrScull == Scull {
		return (rg.Span / 2.)
	}
	return rg._Spread
}

func (rg *Rig) overlap() float64 {
	if rg.RowOrScull == Scull {
		return 2.*rg.Lin - rg.Span
	}
	return rg.Lin - rg._Spread
}

func (rg *Rig) buitenhand() float64 {
	if rg.RowOrScull == Scull {
		return rg.Span - 2.*rg.Lin*math.Cos(rg.CatchAngle)
	}
	return rg._Spread - rg.Lin*math.Cos(rg.CatchAngle)
}

func (rg *Rig) bladearea() float64 {
	if rg.RowOrScull == Scull {
		return 2. * rg._Bladearea
	}
	return rg._Bladearea
}

func (rg *Rig) dcatch() float64 {
	return (rg.Lin * sine(rg.CatchAngle))
}

func (rg *Rig) oarangle(x float64) float64 {
	var dist = rg.dcatch() + x
	var angle = math.Asin(dist / rg.Lin)
	return (angle)
}

// ToJSON exports rig to JSON
func (rg *Rig) ToJSON() (string, error) {
	b, err := json.Marshal(rg)
	if err != nil {
		return "", err
	}
	return string(b), nil
}

// FromJSON sets rig from JSON
func (rg *Rig) FromJSON(s string) error {
	err := json.Unmarshal([]byte(s), rg)
	return err
}

// NewRig initiates a new boat rigging
func NewRig(Lin float64, Mb float64, Lscull float64,
	Span float64, spread float64, RowOrScull string,
	CatchAngle float64, bladearea float64, BladeLength float64,
	Nrowers int32, DragForm float64) *Rig {
	if Lin == 0 {
		Lin = 0.9
	}
	if Mb == 0 {
		Mb = 14
	}
	if Lscull == 0 {
		Lscull = 2.885
	}
	if Span == 0 {
		Span = 1.60
	}

	if CatchAngle == 0 {
		CatchAngle = -0.93
	}
	if bladearea == 0 {
		bladearea = 822.e-4
	}
	if BladeLength == 0 {
		BladeLength = 0.46
	}
	if Nrowers == 0 {
		Nrowers = 1
	}
	if DragForm == 0 {
		DragForm = 1.0
	}
	return &Rig{
		Lin:         Lin,
		Mb:          Mb,
		BladeLength: BladeLength,
		Lscull:      Lscull - 0.5*BladeLength,
		Nrowers:     Nrowers,
		RowOrScull:  RowOrScull,
		Span:        Span,
		CatchAngle:  CatchAngle,
		DragForm:    DragForm,
		_Bladearea:  bladearea,
	}
}

// RigEight 8+
var RigEight = NewRig(1.14, 151, 3.205, 1.6, 0.88, "row", -0.93, 0.1174, 0.545, 8, 1)

// RigQuad 4x
var RigQuad = NewRig(0.9, 52, 2.655, 1.6, 0.88, "scull", -0.93, 0.0822, 0.46, 4, 1.11)

// RigFour 4-
var RigFour = NewRig(1.14, 50, 3.205, 1.6, 0.88, "row", -0.93, 0.1174, 0.545, 4, 1)

// RigPair 2-
var RigPair = NewRig(1.14, 27, 3.205, 1.6, 0.88, "row", -0.93, 0.1174, 0.46, 2, 1.05)

// RigDouble 2x
var RigDouble = NewRig(0.9, 27, 2.655, 1.6, 0.88, "scull", -0.93, 0.0822, 0.46, 2, 1.05)

// RigSingle 1x
var RigSingle = NewRig(0.9, 14, 2.655, 1.6, 0.88, "scull", -0.93, 0.0822, 0.46, 1, 0.98)

// RigCoastalMaas coastal1x
var RigCoastalMaas = NewRig(0.9, 18, 2.655, 1.6, 0.88, "scull", -0.93, 0.0822, 0.46, 1, 1.121)

// RigCoastalLiteboat coastal1x
var RigCoastalLiteboat = NewRig(0.9, 35, 2.655, 1.6, 0.88, "scull", -0.93, 0.0822, 0.46, 1, 1.2)
