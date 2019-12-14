package gorow

import "math"

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
	return (rg.Lin * math.Sin(rg.CatchAngle))
}

func (rg *Rig) oarangle(x float64) float64 {
	var dist = rg.dcatch() + x
	var angle = math.Asin(dist / rg.Lin)
	return (angle)
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
