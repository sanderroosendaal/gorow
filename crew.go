package gorow

import (
	"encoding/json"
	"math"
)

// technique functions
func vcm(vhandle, Strokelength, xhandle float64) float64 {
	vc := vhandle
	if xhandle >= 0 {
		xr := xhandle / Strokelength
		vc = 0.85*vhandle - 0.75*vhandle*math.Pow(xr, 2)
	}
	return vc
}

func vcma(vhandle []float64, Strokelength float64, xhandle []float64) []float64 {
	var vc = make([]float64, len(vhandle))

	for i := range vhandle {
		if vhandle[i] > 0 {
			vc[i] = vcm(vhandle[i], Strokelength, xhandle[i])
		}
	}

	return vc

}

func vha(vcm, Strokelength, xhandle float64) float64 {
	var xr = xhandle / Strokelength
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
	Frac float64
}

func (s StrongMiddle) forceprofile(favg, x float64) float64 {
	var f = (s.Frac * favg * math.Pi * sine(math.Pi*x) / 2.) + (1.-s.Frac)*favg
	return f
}

// StrongMiddle2 stroke profile, alternative with strong middle
type StrongMiddle2 struct {
	Frac float64
}

func (s StrongMiddle2) forceprofile(favg, x float64) float64 {
	var f = (s.Frac * favg * math.Pi * sine(math.Pi*x) / 2.) + 2*(1.-s.Frac)*favg*(1-x)
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
	H1, H2, X1, X2 float64
}

func (s Trapezium) forceprofile(favg, x float64) float64 {
	ratio := s.H1*0.5*s.X2 + s.H2*(0.5-0.5*s.X1)
	f := 0.0
	if x < s.X1 {
		f = favg * s.H1 * x / s.X1
	} else if x > s.X2 {
		f = favg * s.H2 * (1. - x) / (1. - s.X2)
	} else {
		f = favg * (s.H1 + (s.H2-s.H1)*(x-s.X1)/(s.X2-s.X1))
	}

	f = f / ratio
	return f
}

// Trapezium2 stroke profile
type Trapezium2 struct {
	H0, H1, H2, X1, X2 float64
}

func (s Trapezium2) forceprofile(favg, x float64) float64 {
	ratio := s.H1*0.5*s.X2 + s.H2*(0.5-0.5*s.X1) + s.H0
	ratio2 := s.H1*0.5*s.X2 + s.H2*(0.5-0.5*s.X1)
	Frac := ratio2 / ratio
	f := 0.0
	if x < s.X1 {
		f = favg * s.H1 * x / s.X1
	} else if x > s.X2 {
		f = favg * s.H2 * (1. - x) / (1. - s.X2)
	} else {
		f = (s.H1 + (s.H2-s.H1)*(x-s.X1)/(s.X2-s.X1)) * favg
	}

	f = Frac*f + (1-Frac)*favg*s.H0
	f = f / ratio

	return f
}

// FromFile stroke profile class not (yet) implemented

// StrongBegin stroke profile
type StrongBegin struct {
	Frac float64
}

func (s StrongBegin) forceprofile(favg, x float64) float64 {
	f := (2*s.Frac*(1.0-x) + (1. - s.Frac)) * favg
	return f
}

// StrongEnd stroke profile
type StrongEnd struct {
	Frac float64
}

func (s StrongEnd) forceprofile(favg, x float64) float64 {
	f := (2*s.Frac*x + (1. - s.Frac)) * favg
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
	vhand := vhandmax * sine(math.Pi*time/trecovery)
	return vhand
}

func (r SinusRecovery) dxhandle(vavg, trecovery, time float64) float64 {
	dx := 0.5 * (math.Cos(math.Pi*time/trecovery) - 1)
	return dx
}

// SinusRecovery2 recovery profile
type SinusRecovery2 struct {
	P1           float64
	Strokelength float64
}

func (r SinusRecovery2) vhandle(vavg, trecovery, time float64) float64 {
	w1 := math.Pi / r.P1
	w := w1 / trecovery
	vhandmax := w * r.Strokelength / (1 - math.Cos(w*trecovery))
	vhand := -vhandmax * sine(w*time)
	return vhand
}

func (r SinusRecovery2) dxhandle(vavg, trecovery, time float64) float64 {
	w1 := math.Pi / r.P1
	w := w1 / trecovery
	vhandmax := w * r.Strokelength / (1 - math.Cos(w*trecovery))
	dx := vhandmax * (math.Cos(w*time) - 1) / (r.Strokelength)
	return dx
}

// CosinusRecovery recovery profile
type CosinusRecovery struct {
	P1           float64
	Strokelength float64
}

func (r CosinusRecovery) vhandle(vavg, trecovery, time float64) float64 {
	w1 := math.Pi / (2 * r.P1)
	w := w1 / trecovery
	vhandmax := w * r.Strokelength / (sine(w * trecovery))
	vhand := -vhandmax * math.Cos(w*time)
	return vhand
}

func (r CosinusRecovery) dxhandle(vavg, trecovery, time float64) float64 {
	w1 := math.Pi / (2 * r.P1)
	w := w1 / trecovery
	vhandmax := w * r.Strokelength / (sine(w * trecovery))
	dx := vhandmax * sine(w*time) / (r.Strokelength)
	return dx
}

// GenericRecovery not implemented

// CombiRecovery not implemented

// TriangleRecovery recovery profile
type TriangleRecovery struct {
	X1 float64
}

func (r TriangleRecovery) vhandle(vavg, trecovery, time float64) float64 {
	trel := time / trecovery
	if trel < r.X1 {
		vhand := -2 * vavg * trel / r.X1
		return vhand
	}
	vhand := -2 * vavg * (1. - trel) / (1 - r.X1)
	return vhand
}

func (r TriangleRecovery) dxhandle(vavg, trecovery, time float64) float64 {
	trel := time / trecovery
	if trel < r.X1 {
		dx := math.Pow(trel, 2) / r.X1
		return dx
	}
	dx := r.X1
	dx -= (math.Pow(1-trel, 2)) / (1 - r.X1)
	dx += 1 - r.X1
	return -dx
}

// RealisticRecovery not implemented

// Crew class with rower quantities
type Crew struct {
	Mc           float64
	Strokelength float64
	Tempo        float64
	Frac         float64
	// recprofile = sinusrecovery()
	Recoveryprofile RecoveryProfile
	// Strokeprofile = trapezium(X1=0.15,X2=0.5,H2=0.9)
	Strokeprofile ForceProfile
	// technique = technique_meas()
	Maxpower float64
	Maxforce float64
}

func (c *Crew) vcm(vhandle, xhandle float64) float64 {
	return vcm(vhandle, c.Strokelength, xhandle)
}

func (c *Crew) vcma(vhandle, xhandle []float64) []float64 {
	return vcma(vhandle, c.Strokelength, xhandle)
}

func (c *Crew) vha(vcm, xhandle float64) float64 {
	return vha(vcm, c.Strokelength, xhandle)
}

func (c *Crew) vhandle(vavg, trecovery, time float64) float64 {
	return c.Recoveryprofile.vhandle(vavg, trecovery, time)
}

func (c *Crew) dxhandle(vavg, trecovery, time float64) float64 {
	return c.Recoveryprofile.dxhandle(vavg, trecovery, time)
}

func (c *Crew) forceprofile(F, x float64) float64 {
	return c.Strokeprofile.forceprofile(F, x/c.Strokelength)
}

// NewCrew inits Crew instance
func NewCrew(Mc float64, Strokelength float64, Tempo float64, Frac float64,
	Recoveryprofile RecoveryProfile, Strokeprofile ForceProfile,
	Maxpower float64, Maxforce float64) *Crew {
	return &Crew{
		Mc:              Mc,
		Strokelength:    Strokelength,
		Recoveryprofile: Recoveryprofile,
		Strokeprofile:   Strokeprofile,
		Tempo:           Tempo,
		Frac:            Frac,
		Maxpower:        Maxpower,
		Maxforce:        Maxforce,
	}
}

// ToJSON exports crew to JSON
func (c *Crew) ToJSON() (string, error) {
	b, err := json.Marshal(c)
	if err != nil {
		return "", err
	}
	return string(b), nil
}

// FromJSON sets crew from JSON
func (c *Crew) FromJSON(s string) error {
	err := json.Unmarshal([]byte(s), c)
	return err
}
