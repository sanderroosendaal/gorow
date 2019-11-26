/*
Package gorow package to simulate the physics rowing (sport)
*/
package gorow

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

// Global constants

// water density (kg/m^3)

const rho = 999.97

// CLmax maximum lift coefficient for blade
const CLmax = 1.0

// Scull string
const Scull = "scull"

// N Number for bladeforce
const N = 50

// VecLinSpace to create a linear range as a gonum VecDense
func VecLinSpace(start float64, stop float64, N int) *mat.VecDense {
	return mat.NewVecDense(N, LinSpace(start, stop, N))
}

// LinSpace to create a linear range
func LinSpace(start float64, stop float64, N int) []float64 {
	rnge := make([]float64, N)
	var step = (stop - start) / float64(N)
	for x := range rnge {
		rnge[x] = start + step*float64(x)
	}

	return rnge
}

// Constvec creates a VecDense with a single value
func Constvec(value float64, N int) *mat.VecDense {
	rnge := make([]float64, N)
	for i := range rnge {
		rnge[i] = value
	}
	var r = mat.NewVecDense(N, rnge)
	return r
}

// DragEq calculates drag
func DragEq(displacement float64, velo float64,
	alfaref float64, doprint int, constantdrag int) float64 {
	if alfaref == 0 {
		alfaref = 3.5
	}
	var corr = alfaref / 3.5

	// usual coefficient - boat is a spheroid
	const c1 float64 = 1. / 3.
	const c2 float64 = 2. / 3.
	const d1 float64 = 0.06
	const d2 float64 = 28.0
	const d3 float64 = 0.10891
	var a1 float64
	const a3 float64 = 0

	var beam float64 = a1 + d1*math.Pow(displacement, c1)
	var boatlength = d2 * beam
	var wettedArea float64 = a3 + d3*math.Pow(displacement, c2)

	const kinvis float64 = 1.19e-6

	// var D float64 = displacement/rho
	var Re = boatlength * velo / kinvis
	var Cf float64 = 0.075 / (math.Pow((math.Log10(Re) - 2.0), 2))
	var alpha = 0.5 * rho * wettedArea * Cf
	alpha = alpha * corr
	a1 = alpha / 0.8

	if constantdrag == 1 {
		a1 = alfaref
	}

	if doprint == 1 {
		fmt.Print("------ Drag resistance data ----------------\n")
		fmt.Printf("Corr : %f\n", corr)
		fmt.Printf("Beam : %f\n", beam)
		fmt.Printf("Boat length : %f\n", boatlength)
		fmt.Printf("Wetted Area : %f\n", wettedArea)
		fmt.Printf("alpha skin : %f\n", alpha)
		fmt.Printf("alpha total : %f\n", a1)
		fmt.Print("------  Drag resistance data ---------------\n")
		fmt.Print("\n")
	}

	var W2 = a1 * math.Pow(velo, 2)

	return (W2)
}

func vboat(mc float64, mb float64, vc float64) float64 {
	return (mc * vc / (mc + mb))
}

func vhandle(v float64, lin float64, lout float64, mc float64, mb float64) float64 {
	var gamma = mc / (mc + mb)
	var vc = lin * v / (lout + gamma*lin)
	return (vc)
}

func dRecovery(dt, v, vc, dvc, mc, mb, alef, F float64) float64 {
	var dv float64
	var vb = vboat(mc, mb, vc)
	// var dvb = vboat(mc,mb,dvc)
	var Ftot = F - alef*math.Pow((v-vb), 2)
	dv = dt * Ftot / (mb + mc)
	return (dv)
}

func dStroke(dt, v, vc, dvc, mc, mb, alef, F float64) float64 {
	var dv float64
	var vb = vboat(mc, mb, vc)
	var Ftot = F - alef*math.Pow((v-vb), 2)

	dv = dt * Ftot / (mb + mc)

	return (dv)
}

func deFootboard(mc, mb, vs1, vs2 float64) float64 {
	var de float64
	var vt float64

	var vb1 = vboat(mc, mb, vs1)
	var vb2 = vboat(mc, mb, vs2)

	var vmb1 = vt - vb1
	var vmb2 = vt - vb2

	var vmc1 = vt + vs1 - vb1
	var vmc2 = vt + vs2 - vb2

	var e1 = 0.5*(mb*math.Pow(vmb1, 2)) + 0.5*mc*math.Pow(vmc1, 2)
	var e2 = 0.5*(mb*math.Pow(vmb2, 2)) + 0.5*mc*math.Pow(vmc2, 2)
	var eT = 0.5 * (mc + mb) * math.Pow(vt, 2)

	// var de2 = e2 - eT
	// var de1 = e1 - eT
	var samesign = (vs1 * vs2) > 0

	if samesign {
		de = e2 - e1
		if e2-e1 < 0 {
			de = 0
		}
	} else {
		de = e2 - eT
		if de < 0 {
			de = 0
		}
	}

	return (de)
}

// BladeForce calculates the blade slip given a handle force
func BladeForce(oarangle float64, rigging *Rig, vb, fblade float64) []float64 {
	var lin = rigging.lin
	var lscull = rigging.lscull
	var lout = lscull - lin
	var area = rigging.bladearea()

	const N = 50

	// temporary variables
	var phidot1 float64 // = 1.7498034669231506
	var FR float64      // = 99.99941828917106

	var Fprop float64
	var FL float64
	var FD float64
	var CL float64
	var CD float64
	var a float64

	var phidot0 = vb * math.Cos(oarangle) / lout
	var phidot = LinSpace(phidot0, 2*math.Abs(phidot0), N)
	var phidotv = mat.NewVecDense(N, phidot)
	var vblade = Constvec(lout, N)

	vblade.MulElemVec(vblade, phidotv)

	var u1v = make([]float64, N)
	var up = vb * math.Sin(oarangle)
	var av = make([]float64, N)
	var FRv = make([]float64, N)

	// var wg sync.WaitGroup
	// wg.Add(N)

	for i := 0; i < N; i++ {
		// go func(i int) {
		// defer wg.Done()
		var u1 = vblade.AtVec(i) - vb*math.Cos(oarangle)
		u1v[i] = u1
		var u = math.Sqrt(math.Pow(u1, 2) + math.Pow(up, 2)) // fluid velocity
		a = math.Atan(u1 / up)                               // angle of attack
		CD = 2 * CLmax * math.Pow(math.Sin(a), 2)
		CL = CLmax * math.Sin(2*a)

		FL = 0.5 * CL * rho * area * math.Pow(u, 2)
		FD = 0.5 * CD * rho * area * math.Pow(u, 2)
		FRv[i] = math.Sqrt(math.Pow(FL, 2) + math.Pow(FD, 2))
		// Fprop = FRv[i] * math.Cos(oarangle)
		av[i] = a
		// } (i)
	}

	// wg.Wait()

	phidot1 = srinterpol1(phidot, FRv, fblade)

	var vblade1 = phidot1 * lout
	var u1 = vblade1 - vb*math.Cos(oarangle)
	up = vb * math.Sin(oarangle)

	var u = math.Sqrt(math.Pow(u1, 2) + math.Pow(up, 2)) // fluid velocity
	a = math.Atan(u1 / up)                               // angle of attack

	CD = 2 * CLmax * math.Pow(math.Sin(a), 2)
	CL = CLmax * math.Sin(2.*a)

	FL = 0.5 * CL * rho * area * math.Pow(u, 2)
	FD = 0.5 * CD * rho * area * math.Pow(u, 2)
	FR = math.Sqrt(math.Pow(FL, 2) + math.Pow(FD, 2))
	Fprop = FR * math.Cos(oarangle)

	return []float64{
		phidot1, FR, Fprop, FL, FD, CL, CD, a,
	}
}

// EnergyBalance calculates one stroke with average handle force as input
func EnergyBalance(
	F float64,
	crew *Crew,
	rigging *Rig,
	v0 float64,
	dt float64,
	catchacceler float64,
) []float64 {
	var dv, vavg, vend, ratio, power float64

	vb0 := v0

	if catchacceler > 50 {
		catchacceler = 50
	}

	lin := rigging.lin
	lscull := rigging.lscull
	lout := lscull - lin
	tempo := crew.tempo
	mc := crew.mc
	mb := rigging.mb
	recprofile := crew.recoveryprofile
	d := crew.strokelength
	Nrowers := rigging.Nrowers
	dragform := rigging.dragform

	if catchacceler < 2 {
		catchacceler = 2
	}

	aantal := 1 + int(math.Round(60/(tempo*dt)))

	time := LinSpace(0, 60./tempo, aantal)

	vs := LinSpace(0, 0, aantal)
	vb := LinSpace(0, 0, aantal)
	vc := LinSpace(0, 0, aantal)

	oarangle := LinSpace(0, 0, aantal)
	xblade := LinSpace(0, 0, aantal)
	Fhandle := LinSpace(0, 0, aantal)
	Fblade := LinSpace(0, 0, aantal)
	Fprop := LinSpace(0, 0, aantal)

	Pbladeslip := LinSpace(0, 0, aantal)

	xdotdot := LinSpace(0, 0, aantal)
	zdotdot := LinSpace(0, 0, aantal)
	ydotdot := LinSpace(0, 0, aantal)

	xdot := LinSpace(v0, v0, aantal)
	ydot := LinSpace(v0, v0, aantal)
	zdot := LinSpace(v0, v0, aantal)

	Pf := LinSpace(0, 0, aantal)
	Foarlock := LinSpace(0, 0, aantal)
	Flift := LinSpace(0, 0, aantal)
	Fbldrag := LinSpace(0, 0, aantal)
	attackangle := LinSpace(0, 0, aantal)
	Clift := LinSpace(0, 0, aantal)
	Cdrag := LinSpace(0, 0, aantal)

	handlepos := 0

	// initial handle and boat velocities
	vs[0] = v0
	vb[0] = vb0
	vc[0] = ((float64(Nrowers)*mc+mb)*vs[0] - mb*vb[0]) / (float64(Nrowers) * mc)
	oarangle[0] = rigging.oarangle(0)
	xblade[0] = -lout * math.Sin(oarangle[0])

	i := 1

	vcstroke := 0
	vcstroke2 := 0

	vblade := xdot[aantal-1]

	for vcstroke < vcstroke2 {
		// blade entry loop
	}

	return []float64{1, 2}
}
