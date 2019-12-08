package gorow

import "math"

// Rig holds boat rigging parameters and has methods to manipulate them
type Rig struct {
	lin         float64
	mb          float64
	bladelength float64
	lscull      float64
	Nrowers     int32
	roworscull  string
	span        float64
	catchangle  float64
	dragform    float64
	_Spread     float64
	_Bladearea  float64
}

func (rg *Rig) spread() float64 {
	if rg.roworscull == Scull {
		return (rg.span / 2.)
	}
	return rg._Spread
}

func (rg *Rig) overlap() float64 {
	if rg.roworscull == Scull {
		return 2.*rg.lin - rg.span
	}
	return rg.lin - rg._Spread
}

func (rg *Rig) buitenhand() float64 {
	if rg.roworscull == Scull {
		return rg.span - 2.*rg.lin*math.Cos(rg.catchangle)
	}
	return rg._Spread - rg.lin*math.Cos(rg.catchangle)
}

func (rg *Rig) bladearea() float64 {
	if rg.roworscull == Scull {
		return 2. * rg._Bladearea
	}
	return rg._Bladearea
}

func (rg *Rig) dcatch() float64 {
	return (rg.lin * math.Sin(rg.catchangle))
}

func (rg *Rig) oarangle(x float64) float64 {
	var dist = rg.dcatch() + x
	var angle = math.Asin(dist / rg.lin)
	return (angle)
}

// NewRig initiates a new boat rigging
func NewRig(lin float64, mb float64, lscull float64,
	span float64, spread float64, roworscull string,
	catchangle float64, bladearea float64, bladelength float64,
	Nrowers int32, dragform float64) *Rig {
	if lin == 0 {
		lin = 0.9
	}
	if mb == 0 {
		mb = 14
	}
	if lscull == 0 {
		lscull = 2.885
	}
	if span == 0 {
		span = 1.60
	}

	if catchangle == 0 {
		catchangle = -0.93
	}
	if bladearea == 0 {
		bladearea = 822.e-4
	}
	if bladelength == 0 {
		bladelength = 0.46
	}
	if Nrowers == 0 {
		Nrowers = 1
	}
	if dragform == 0 {
		dragform = 1.0
	}
	return &Rig{
		lin:         lin,
		mb:          mb,
		bladelength: bladelength,
		lscull:      lscull - 0.5*bladelength,
		Nrowers:     Nrowers,
		roworscull:  roworscull,
		span:        span,
		catchangle:  catchangle,
		dragform:    dragform,
		_Bladearea:  bladearea,
	}
}
