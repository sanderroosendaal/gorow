package gorow

import (
	"math"
)

// technique functions
func vcm(vhandle, strokelength, xhandle float64) float64 {
	var vc = vhandle
	if xhandle > 0 {
		var xr = xhandle / strokelength
		vc = 0.85*vhandle - 0.75*vhandle*math.Pow(xr, 2)
	}
	return vc
}

func vcma(vhandle []float64, strokelength float64, xhandle []float64) []float64 {
	var vc = make([]float64, len(vhandle))

	for i := range vhandle {
		if vhandle[i] > 0 {
			vc[i] = vcm(vhandle[i], strokelength, xhandle[i])
		}
	}

	return vc

}

func vha(vcm, strokelength, xhandle float64) float64 {
	var xr = xhandle / strokelength
	var vh = vcm / (0.85 - 0.75*math.Pow(xr, 2))

	return vh
}

// Crew class with rower quantities
type Crew struct {
	mc           float64
	strokelength float64
	tempo        float64
	frac         float64
	// recprofile = sinusrecovery()
	// strokeprofile = trapezium(x1=0.15,x2=0.5,h2=0.9)
	// technique = technique_meas()
	maxpower float64
	maxforce float64
}

func (c *Crew) vcm(vhandle, xhandle float64) float64 {
	return vcm(vhandle, c.strokelength, xhandle)
}

func (c *Crew) vcma(vhandle, xhandle []float64) []float64 {
	return vcma(vhandle, c.strokelength, xhandle)
}

func (c *Crew) vha(vcm, xhandle float64) float64 {
	return vha(vcm, c.strokelength, xhandle)
}

// NewCrew inits Crew instance
func NewCrew(mc float64, strokelength float64, tempo float64, frac float64,
	maxpower float64, maxforce float64) *Crew {
	return &Crew{
		mc:           mc,
		strokelength: strokelength,
		tempo:        tempo,
		frac:         frac,
		maxpower:     maxpower,
		maxforce:     maxforce,
	}
}
