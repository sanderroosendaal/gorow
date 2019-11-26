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

// ForceProfile interface forceprofile
type ForceProfile interface {
	forceprofile(favg, x float64) float64
}

// StrongMiddle stroke profile with strong middle
type StrongMiddle struct {
	frac float64
}

func (s StrongMiddle) forceprofile(favg, x float64) float64 {
	var f = (s.frac * favg * math.Pi * math.Sin(math.Pi*x) / 2.) + (1.-s.frac)*favg
	return f
}

// StrongMiddle2 stroke profile, alternative with strong middle
type StrongMiddle2 struct {
	frac float64
}

func (s StrongMiddle2) forceprofile(favg, x float64) float64 {
	var f = (s.frac * favg * math.Pi * math.Sin(math.Pi*x) / 2.) + 2*(1.-s.frac)*favg*(1-x)
	return f
}

// Flat stroke profile
type Flat struct {
}

func (s Flat) forceprofile(favg, x float64) float64 {
	return favg
}

// Trapezium stroke profile
type Trapezium struct {
	h1, h2, x1, x2 float64
}

func (s Trapezium) forceprofile(favg, x float64) float64 {
	ratio := s.h1*0.5*s.x2 + s.h2*(0.5-0.5*s.x1)
	f := 0.0
	if x < s.x1 {
		f = favg * s.h1 * x / s.x1
	} else if x > s.x2 {
		f = favg * s.h2 * (1. - x) / (1. - s.x2)
	} else {
		f = (s.h1 + (s.h2-s.h1)*(x-s.x1)/(s.x2-s.x1)) * favg
	}

	f = f / ratio
	return f
}

// Trapezium2 stroke profile
type Trapezium2 struct {
	h0, h1, h2, x1, x2 float64
}

func (s Trapezium2) forceprofile(favg, x float64) float64 {
	ratio := s.h1*0.5*s.x2 + s.h2*(0.5-0.5*s.x1) + s.h0
	ratio2 := s.h1*0.5*s.x2 + s.h2*(0.5-0.5*s.x1)
	frac := ratio2 / ratio
	f := 0.0
	if x < s.x1 {
		f = favg * s.h1 * x / s.x1
	} else if x > s.x2 {
		f = favg * s.h2 * (1. - x) / (1. - s.x2)
	} else {
		f = (s.h1 + (s.h2-s.h1)*(x-s.x1)/(s.x2-s.x1)) * favg
	}

	f = frac*f + (1-frac)*favg*s.h0
	f = f / ratio

	return f
}

// FromFile stroke profile class not (yet) implemented

// StrongBegin stroke profile
type StrongBegin struct {
	frac float64
}

func (s StrongBegin) forceprofile(favg, x float64) float64 {
	f := (2*s.frac*(1.0-x) + (1. - s.frac)) * favg
	return f
}

// StrongEnd stroke profile
type StrongEnd struct {
	frac float64
}

func (s StrongEnd) forceprofile(favg, x float64) float64 {
	f := (2*s.frac*x + (1. - s.frac)) * favg
	return f
}

// Crew class with rower quantities
type Crew struct {
	mc           float64
	strokelength float64
	tempo        float64
	frac         float64
	// recprofile = sinusrecovery()
	// strokeprofile = trapezium(x1=0.15,x2=0.5,h2=0.9)
	strokeprofile ForceProfile
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
	strokeprofile ForceProfile,
	maxpower float64, maxforce float64) *Crew {
	return &Crew{
		mc:            mc,
		strokelength:  strokelength,
		strokeprofile: strokeprofile,
		tempo:         tempo,
		frac:          frac,
		maxpower:      maxpower,
		maxforce:      maxforce,
	}
}
