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

// RecoveryProfile interface
type RecoveryProfile interface {
	vhandle(vavg, trecovery, time float64) float64
	dxhandle(vavg, trecovery, time float64) float64
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

// FlatRecovery recovery profile
type FlatRecovery struct {
}

func (r FlatRecovery) vhandle(vavg, trecovery, time float64) float64 {
	return -vavg
}

func (r FlatRecovery) dxhandle(vavg, trecovery, time float64) float64 {
	return -time / trecovery
}

// SinusRecovery recovery profile
type SinusRecovery struct {
}

func (r SinusRecovery) vhandle(vavg, trecovery, time float64) float64 {
	vhandmax := -math.Pi * vavg / 2.
	vhand := vhandmax * math.Sin(math.Pi*time/trecovery)
	return vhand
}

func (r SinusRecovery) dxhandle(vavg, trecovery, time float64) float64 {
	dx := 0.5 * (math.Cos(math.Pi*time/trecovery) - 1)
	return dx
}

// SinusRecovery2 recovery profile
type SinusRecovery2 struct {
	p1           float64
	strokelength float64
}

func (r SinusRecovery2) vhandle(vavg, trecovery, time float64) float64 {
	w1 := math.Pi / r.p1
	w := w1 / trecovery
	vhandmax := w * r.strokelength / (1 - math.Cos(w*trecovery))
	vhand := -vhandmax * math.Sin(w*time)
	return vhand
}

func (r SinusRecovery2) dxhandle(vavg, trecovery, time float64) float64 {
	w1 := math.Pi / r.p1
	w := w1 / trecovery
	vhandmax := w * r.strokelength / (1 - math.Cos(w*trecovery))
	dx := vhandmax * (math.Cos(w*time) - 1) / (r.strokelength)
	return dx
}

// CosinusRecovery recovery profile
type CosinusRecovery struct {
	p1           float64
	strokelength float64
}

func (r CosinusRecovery) vhandle(vavg, trecovery, time float64) float64 {
	w1 := math.Pi / (2 * r.p1)
	w := w1 / trecovery
	vhandmax := w * r.strokelength / (math.Sin(w * trecovery))
	vhand := -vhandmax * math.Cos(w*time)
	return vhand
}

func (r CosinusRecovery) dxhandle(vavg, trecovery, time float64) float64 {
	w1 := math.Pi / (2 * r.p1)
	w := w1 / trecovery
	vhandmax := w * r.strokelength / (math.Sin(w * trecovery))
	dx := vhandmax * math.Sin(w*time) / (r.strokelength)
	return dx
}

// GenericRecovery not implemented

// CombiRecovery not implemented

// TriangleRecovery recovery profile
type TriangleRecovery struct {
	x1 float64
}

func (r TriangleRecovery) vhandle(vavg, trecovery, time float64) float64 {
	trel := time / trecovery
	if trel < r.x1 {
		vhand := -2 * vavg * trel / r.x1
		return vhand
	}
	vhand := -2 * vavg * (1. - trel) / (1 - r.x1)
	return vhand
}

func (r TriangleRecovery) dxhandle(vavg, trecovery, time float64) float64 {
	trel := time / trecovery
	if trel < r.x1 {
		dx := math.Pow(trel, 2) / r.x1
		return dx
	}
	dx := r.x1
	dx -= (math.Pow(1-trel, 2)) / (1 - r.x1)
	dx += 1 - r.x1
	return -dx
}

// RealisticRecovery not implemented

// Crew class with rower quantities
type Crew struct {
	mc           float64
	strokelength float64
	tempo        float64
	frac         float64
	// recprofile = sinusrecovery()
	recoveryprofile RecoveryProfile
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

func (c *Crew) vhandle(vavg, trecovery, time float64) float64 {
	return c.recoveryprofile.vhandle(vavg, trecovery, time)
}

func (c *Crew) dxhandle(vavg, trecovery, time float64) float64 {
	return c.recoveryprofile.dxhandle(vavg, trecovery, time)
}

func (c *Crew) forceprofile(F, x float64) float64 {
	return c.strokeprofile.forceprofile(F,x)
}

// NewCrew inits Crew instance
func NewCrew(mc float64, strokelength float64, tempo float64, frac float64,
	recoveryprofile RecoveryProfile, strokeprofile ForceProfile,
	maxpower float64, maxforce float64) *Crew {
	return &Crew{
		mc:              mc,
		strokelength:    strokelength,
		recoveryprofile: recoveryprofile,
		strokeprofile:   strokeprofile,
		tempo:           tempo,
		frac:            frac,
		maxpower:        maxpower,
		maxforce:        maxforce,
	}
}
