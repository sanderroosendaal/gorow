/*
Package gorow package to simulate the physics rowing (sport)
*/
package gorow

import (
	"errors"
	"fmt"
	"math"

	"gonum.org/v1/gonum/stat"
)

// Global constants

// water density (kg/m^3)

const rho = 999.97

// RhoAir air density
const RhoAir = 1.226 // kg.m3
// Cdw drag constant
const Cdw = 1.1 // for all boats, big approximation
const crewarea = 2.0
const scalepower = 0.67

// CLmax maximum lift coefficient for blade
const CLmax = 1.0

// Scull string
const Scull = "scull"

// N Number for bladeforce
const N = 50

// Rowing global parameters
const alfa = 3.06 // best fit to Kleshnev data for single
const alfaatkinson = 3.4

// slice helper functions
func sliceminmax(vs []float64) (min float64, max float64) {
	biggest, smallest := vs[0], vs[0]

	for _, v := range vs {
		if v > biggest {
			biggest = v
		}
		if v < smallest {
			smallest = v
		}
	}
	return smallest, biggest
}

func slicesadd(vs ...[]float64) []float64 {
	out := make([]float64, len(vs[0]))
	for i := 0; i < len(vs[0]); i++ {
		for _, v := range vs {
			out[i] += v[i]
		}
	}
	return out
}

// LinSpace to create a Linear range
func LinSpace(start float64, stop float64, N int) ([]float64, error) {
	if N <= 0 {
		return nil, errors.New("LinSpace: Length of slice must by > 0")
	}
	rnge := make([]float64, N)
	var step = (stop - start) / float64(N-1)
	for x := range rnge {
		rnge[x] = start + step*float64(x)
	}

	return rnge, nil
}

// ConstSpace to create a array of value
func ConstSpace(value float64, N int) []float64 {
	rnge := make([]float64, N)
	for x := range rnge {
		rnge[x] = value
	}
	return rnge
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
	var Cf float64 = 0.075 / ((math.Log10(Re) - 2.0) * (math.Log10(Re) - 2.0))
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

	var W2 = a1 * velo * velo

	return (W2)
}

func vboat(Mc float64, mb float64, vc float64) float64 {
	return (Mc * vc / (Mc + mb))
}

func vhandle(v float64, Lin float64, lout float64, Mc float64, mb float64) float64 {
	var gamma = Mc / (Mc + mb)
	var vc = Lin * v / (lout + gamma*Lin)
	return (vc)
}

func dRecovery(dt, v, vc, dvc, Mc, mb, alef, F float64) float64 {
	var dv float64
	var vb = vboat(Mc, mb, vc)
	// var dvb = vboat(Mc,mb,dvc)
	var Ftot = F - alef*(v-vb)*(v-vb)
	dv = dt * Ftot / (mb + Mc)
	return (dv)
}

func dStroke(dt, v, vc, dvc, Mc, mb, alef, F float64) float64 {
	var dv float64
	var vb = vboat(Mc, mb, vc)
	var Ftot = F - alef*(v-vb)*(v-vb)

	dv = dt * Ftot / (mb + Mc)

	return (dv)
}

func deFootboard(Mc, mb, vs1, vs2 float64) float64 {
	var de float64
	var vt float64

	var vb1 = vboat(Mc, mb, vs1)
	var vb2 = vboat(Mc, mb, vs2)

	var vmb1 = vt - vb1
	var vmb2 = vt - vb2

	var vMc1 = vt + vs1 - vb1
	var vMc2 = vt + vs2 - vb2

	var e1 = 0.5*(mb*vmb1*vmb1) + 0.5*Mc*vMc1*vMc1
	var e2 = 0.5*(mb*vmb2*vmb2) + 0.5*Mc*vMc2*vMc2
	var eT = 0.5 * (Mc + mb) * vt * vt

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
func BladeForce(oarangle float64, rigging *Rig, vb, fblade float64) ([]float64, error) {
	var Lin = rigging.Lin
	var Lscull = rigging.Lscull
	var lout = Lscull - Lin
	var area = rigging.bladearea()

	const N = 50

	// Temporary variables
	var phidot1 float64 // = 1.7498034669231506
	var FR float64      // = 99.99941828917106

	var Fprop float64
	var FL float64
	var FD float64
	var CL float64
	var CD float64
	var a float64

	var phidot0 = vb * math.Cos(oarangle) / lout
	phidot, err := LinSpace(phidot0, 2*math.Abs(phidot0), N)
	if err != nil {
		return nil, err
	}

	var vblade = make([]float64, N)

	for i, v := range phidot {
		vblade[i] = v * lout
	}

	var u1v = make([]float64, N)
	var up = vb * math.Sin(oarangle)
	var av = make([]float64, N)
	var FRv = make([]float64, N)

	// var wg sync.WaitGroup
	// wg.Add(N)

	for i := 0; i < N; i++ {
		// go func(i int) {
		// defer wg.Done()
		var u1 = vblade[i] - vb*math.Cos(oarangle)
		u1v[i] = u1
		var u = math.Sqrt(u1*u1 + up*up) // fluid velocity
		a = math.Atan(u1 / up)           // angle of attack
		CD = 2 * CLmax * math.Sin(a) * math.Sin(a)
		CL = CLmax * math.Sin(2*a)

		FL = 0.5 * CL * rho * area * u * u
		FD = 0.5 * CD * rho * area * u * u
		FRv[i] = math.Sqrt(FL*FL + FD*FD)
		// Fprop = FRv[i] * math.Cos(oarangle)
		av[i] = a

		// } (i)
	}

	// wg.Wait()

	phidot1, _ = srinterpol3(phidot, FRv, fblade)

	var vblade1 = phidot1 * lout
	var u1 = vblade1 - vb*math.Cos(oarangle)
	up = vb * math.Sin(oarangle)

	var u = math.Sqrt(u1*u1 + up*up) // fluid velocity
	a = math.Atan(u1 / up)           // angle of attack

	CD = 2 * CLmax * math.Sin(a) * math.Sin(a)
	CL = CLmax * math.Sin(2.*a)

	FL = 0.5 * CL * rho * area * u * u
	FD = 0.5 * CD * rho * area * u * u
	FR = math.Sqrt(FL*FL + FD*FD)
	Fprop = FR * math.Cos(oarangle)

	return []float64{
		phidot1, // 0
		FR,      // 1
		Fprop,   // 2
		FL,      // 3
		FD,      // 4
		CL,      // 5
		CD,      // 6
		a,       // 7
	}, nil
}

// EnergyBalance calculates one stroke with average handle force as input
func EnergyBalance(
	F float64,
	crew *Crew,
	rigging *Rig,
	v0 float64,
	dt float64,
	catchacceler float64,
	windv float64,
	dowind bool,
) ([]float64, error) {

	var dv, vavg, vend, ratio, power float64

	vb0 := v0

	if catchacceler > 50 {
		catchacceler = 50
	}

	Lin := rigging.Lin
	Lscull := rigging.Lscull
	lout := Lscull - Lin
	Tempo := crew.Tempo
	Mc := crew.Mc
	mb := rigging.Mb

	d := crew.Strokelength
	Nrowers := rigging.Nrowers
	DragForm := rigging.DragForm

	if catchacceler < 2 {
		catchacceler = 2
	}

	aantal := 1 + int(math.Round(60/(Tempo*dt)))

	time, err := LinSpace(0, 60./Tempo, aantal)
	if err != nil {
		return nil, err
	}

	vs := make([]float64, aantal)
	vb := make([]float64, aantal)
	vc := make([]float64, aantal)

	oarangle := make([]float64, aantal)
	xblade := make([]float64, aantal)
	Fhandle := make([]float64, aantal)
	Fblade := make([]float64, aantal)
	Fprop := make([]float64, aantal)

	Pbladeslip := make([]float64, aantal)

	xdotdot := make([]float64, aantal)
	zdotdot := make([]float64, aantal)
	ydotdot := make([]float64, aantal)

	xdot := ConstSpace(v0, aantal)
	ydot := ConstSpace(v0, aantal)
	zdot := ConstSpace(v0, aantal)

	Pf := make([]float64, aantal)

	Flift := make([]float64, aantal)
	Fbldrag := make([]float64, aantal)
	attackangle := make([]float64, aantal)
	Clift := make([]float64, aantal)
	Cdrag := make([]float64, aantal)

	mtotal := float64(Nrowers)*Mc + mb
	Mcrew := float64(Nrowers) * Mc

	var handlepos float64

	// initial handle and boat velocities
	vs[0] = v0
	vb[0] = vb0
	vc[0] = ((float64(Nrowers)*Mc+mb)*vs[0] - mb*vb[0]) / (float64(Nrowers) * Mc)
	oarangle[0] = rigging.oarangle(0)
	xblade[0] = -lout * math.Sin(oarangle[0])

	i := 1

	var vcstroke float64
	var vcstroke2 float64 = 1

	for vcstroke < vcstroke2 && i < aantal {
		// blade entry loop
		vhand := catchacceler * (time[i] - time[0])
		vcstroke = crew.vcm(vhand, handlepos)
		// phidot := vb[i-1] * math.Cos(oarangle[i-1])
		// vhand := phidot * Lin * math.Cos(oarangle[i-1])
		ydot[i] = vcstroke

		alfaref := alfa * DragForm
		Fdrag := DragEq(mtotal, xdot[i-1], alfaref, 0, 0)
		zdotdot[i] = -Fdrag / mtotal

		vw := windv - vcstroke - zdot[i-1]
		Fwind := 0.0
		if dowind {
			Fwind = 0.5 * crewarea * Cdw * RhoAir * math.Pow(float64(Nrowers), scalepower) * vw * math.Abs(vw)
		}
		zdotdot[i] = zdotdot[i] + Fwind/mtotal
		zdot[i] = zdot[i-1] + dt*zdotdot[i]
		xdot[i] = zdot[i] - (Mcrew/mtotal)*ydot[i]

		Fi := crew.forceprofile(F, handlepos)
		Fbladei := Fi * Lin / lout
		// fmt.Println(i, Fbladei, Fi)
		res, err := BladeForce(oarangle[i-1], rigging, vb[i-1], Fbladei)
		if err != nil {
			return nil, err
		}

		phidot2 := res[0]
		vhand2 := phidot2 * Lin * math.Cos(oarangle[i-1])
		// fmt.Println(i, phidot2, Lin, oarangle[i-1], math.Cos(oarangle[i-1]), vhand2)

		vcstroke2 = crew.vcm(vhand2, handlepos)

		vs[i] = zdot[i]
		vc[i] = xdot[i] + ydot[i]
		vb[i] = xdot[i]

		ydotdot[i] = (ydot[i] - ydot[i-1]) / dt
		xdotdot[i] = zdotdot[i] - ((Mcrew)/(mtotal))*ydotdot[i]

		handlepos += ydot[i] * dt
		Fhandle[i] = 0

		oarangle[i] = rigging.oarangle(handlepos)

		i++
	}

	// stroke
	for handlepos < d && i < len(time) {
		Fi := crew.forceprofile(F, handlepos)
		Fhandle[i-1] = Fi
		Fblade[i-1] = Fi * Lin / lout

		res, err := BladeForce(oarangle[i-1], rigging, vb[i-1], Fblade[i-1])
		if err != nil {
			return nil, err
		}
		phidot := res[0]

		Fprop[i-1] = res[2] * float64(Nrowers)
		Flift[i-1] = res[3] * float64(Nrowers)
		Fbldrag[i-1] = res[4] * float64(Nrowers)
		Clift[i-1] = res[5]
		Cdrag[i-1] = res[6]
		attackangle[i-1] = res[7]

		vhand := phidot * Lin * math.Cos(oarangle[i-1])

		vcstroke := crew.vcm(vhand, handlepos)
		Pbladeslip[i-1] = float64(Nrowers) * res[1] * (phidot*lout - vb[i-1]*math.Cos(oarangle[i-1]))

		alfaref := alfa * DragForm
		Fdrag := DragEq(mtotal, xdot[i-1], alfaref, 0, 0)
		zdotdot[i] = (Fprop[i-1] - Fdrag) / mtotal
		vw := windv - vcstroke - zdot[i-1]
		Fwind := 0.0
		if dowind {
			Fwind = 0.5 * crewarea * Cdw * RhoAir * math.Pow(float64(Nrowers), scalepower) * vw * math.Abs(vw)
		}
		zdotdot[i] = zdotdot[i] + Fwind/mtotal

		zdot[i] = zdot[i-1] + dt*zdotdot[i]

		ydot[i] = vcstroke
		xdot[i] = zdot[i] - ((Mcrew)/(mtotal))*ydot[i]

		handlepos = handlepos + vhand*dt
		vs[i] = zdot[i]
		vc[i] = xdot[i] + ydot[i]
		vb[i] = xdot[i]

		ydotdot[i] = (ydot[i] - ydot[i-1]) / dt
		xdotdot[i] = zdotdot[i] - (Mcrew/mtotal)*ydotdot[i]

		Pf[i-1] = float64(Nrowers) * Fblade[i-1] * xdot[i] * math.Cos(oarangle[i-1])

		oarangle[i] = rigging.oarangle(handlepos)

		i++
	}

	i = i - 1

	// recovery
	maxtime := 60. / Tempo
	trecovery := maxtime - time[i]
	ratio = time[i] / maxtime

	vavgrec := d / trecovery
	vcrecovery := make([]float64, aantal)

	for k := i + 1; k < aantal; k++ {
		vhand := crew.vhandle(vavgrec, trecovery, time[k]-time[i])
		vcrecovery[k] = crew.vcm(vhand, handlepos)

		alfaref := alfa * DragForm
		Fdrag := DragEq(mtotal, xdot[k-1], alfaref, 0, 0)
		zdotdot[k] = -Fdrag / mtotal

		vw := windv - vcrecovery[k] - zdot[k-1]
		Fwind := 0.0
		if dowind {
			Fwind = 0.5 * crewarea * Cdw * RhoAir * math.Pow(float64(Nrowers), scalepower) * vw * math.Abs(vw)
		}

		zdotdot[k] = zdotdot[k] + Fwind/mtotal

		zdot[k] = zdot[k-1] + dt*zdotdot[k]
		ydot[k] = vcrecovery[k]
		xdot[k] = zdot[k] - ydot[k]*Mcrew/mtotal

		vs[k] = zdot[k]
		vc[k] = xdot[k] + ydot[k]
		vb[k] = xdot[k]

		ydotdot[k] = (ydot[k] - ydot[k-1]) / dt
		xdotdot[k] = zdotdot[k] - ydotdot[k]*Mcrew/mtotal

		handlepos = d + d*crew.dxhandle(vavgrec, trecovery, time[k]-time[i])
		oarangle[k] = rigging.oarangle(handlepos)
	}

	Pq := make([]float64, aantal)
	Pqrower := make([]float64, aantal)
	Pdiss := make([]float64, aantal)
	//	Ekinb := make([]float64, aantal)
	// Ekinc := make([]float64, aantal)
	// Ef := make([]float64, aantal)
	// Eq := make([]float64, aantal)
	// Eblade := make([]float64, aantal)
	// Eqrower := make([]float64, aantal)
	// Ediss := make([]float64, aantal)
	// Eleg := make([]float64, aantal)
	// Ehandle := make([]float64, aantal)
	// Ew := make([]float64, aantal)
	var Ewmin float64
	var Pwmin float64
	Pw := make([]float64, aantal)
	Pmb := make([]float64, aantal)
	PMc := make([]float64, aantal)
	Phandle := make([]float64, aantal)
	Pleg := make([]float64, aantal)

	alfaref := alfa * DragForm

	// blade positions
	for i := 0; i < aantal; i++ {
		xdot[i] = vb[i]
		zdot[i] = vs[i]
		ydot[i] = vc[i] - vb[i]
	}

	xdotdot[1] = (xdot[1] - xdot[0]) / dt
	ydotdot[1] = (ydot[1] - ydot[0]) / dt

	xdotmean := stat.Mean(xdot, nil)

	Pwmin = DragEq(mtotal, xdotmean, alfaref, 0, 0) * xdotmean

	for i := 0; i < aantal; i++ {
		Pq[i] = Mcrew * (xdotdot[i] + ydotdot[i]) * ydot[i]
		Pw[i] = DragEq(mtotal, xdot[i], alfaref, 0, 0) * xdot[i]
		Pmb[i] = mb * xdot[i] * xdotdot[i]
		PMc[i] = Mcrew * (xdot[i] + ydot[i]) * (xdotdot[i] + ydotdot[i])
		Phandle[i] = float64(Nrowers) * Fhandle[i] * xdot[i] * math.Cos(oarangle[i])
		Pleg[i] = float64(Nrowers) * Mc * (xdotdot[i] + ydotdot[i]) * ydot[i]
		Pqrower[i] = math.Abs(Pq[i])
		Pdiss[i] = Pqrower[i] - Pq[i]

	}

	// Ekinb = cumsum(Pmb, dt)
	// Ekinc = cumsum(PMc, dt)
	// Ef = cumsum(Pf, dt)
	// Eq = cumsum(Pq, dt)
	Eblade := cumsum(Pbladeslip, dt)
	// Eqrower = cumsum(Pqrower, dt)
	Ediss := cumsum(Pdiss, dt)
	Ew := cumsum(Pw, dt)

	// Eleg = cumsum(Pleg, dt)
	// Ehandle = cumsum(Phandle, dt)

	Ewmin = Pwmin * 60. / Tempo

	Ekin0 := 0.5 * mtotal * zdot[0] * zdot[0]
	Ekinend := 0.5 * mtotal * zdot[aantal-1] * zdot[aantal-1]
	Eloss := Ekin0 - Ekinend

	// calculate vavg, vmin, vmax, energy, efficiency, power
	dv = zdot[aantal-1] - zdot[0]
	vavg = stat.Mean(xdot, nil)
	vend = zdot[aantal-1]
	_, energy := sliceminmax(slicesadd(Ew, Ediss, Eblade))

	energy -= Eloss
	_, efficiency := sliceminmax(Ew)
	efficiency -= Eloss
	efficiency = efficiency / energy
	energy = energy / float64(Nrowers)
	power = energy * Tempo / 60.
	vmin, vmax := sliceminmax(xdot)

	// calculate check
	/*
			decel = -(abs(xdotdot[index_offset:])-xdotdot[index_offset:])/2.
		    indices = decel.nonzero()
		    decelmean = mean(decel[indices]) */

	CNCheck := 0.0 // np.std(decel[indices])**2

	// RIM parameters
	RIMcheck := vmax - vmin
	RIME := 0.0 // max(cumsum(xdot-vmin)*dt)
	_, maxEw := sliceminmax(Ew)
	dragEff := Ewmin / maxEw
	/*
	   try:
	       t4 = time[index_offset+min(where(decel==0)[0])]
	       t3 = time[index_offset+max(where(decel==0)[0])]
	   except ValueError:
	       t4 = 1.0
	       t3 = t4
	*/
	// amin, _ := Sliceminmax(xdotdot[2:])
	RIMcatchE := 0.0 // -(amin/t4)
	RIMcatchD := 0.0 // t4+max(time)-t3

	catchacceler = ydotdot[aantal-1] - xdotdot[aantal-1]
	if catchacceler < 5.0 {
		catchacceler = 5.0
	}

	return []float64{
		dv,           // 0
		vend,         // 1
		vavg,         // 2
		ratio,        // 3
		energy,       // 4
		power,        // 5
		efficiency,   // 6
		vmax,         // 7
		vmin,         // 8
		CNCheck,      // 9
		RIME,         // 10
		RIMcheck,     // 11
		RIMcatchE,    // 12
		RIMcatchD,    // 13
		catchacceler, // 14
		dragEff,      // 15
	}, nil

}

// Stroke calculates a few (aantal) strokes and returns parameters averaged over those strokes
func Stroke(
	F float64,
	crew *Crew,
	rigging *Rig,
	v0 float64,
	dt float64,
	aantal int,
	catchacceler float64,
	dowind bool,
	windv float64,
) ([]float64, error) {

	var dv, vavg, vend, vmin, vmax, ratio, energy, power, eff float64

	var CNCheck, RIME, RIMcheck, RIMcatchE, RIMcatchD, DragEff float64
	tcatchacceler := catchacceler

	for i := 0; i < aantal; i++ {
		res, err := EnergyBalance(F, crew, rigging, v0, dt, tcatchacceler, windv, dowind)
		if err != nil {
			return nil, err
		}

		dv = dv + res[0]
		vend = vend + res[1]
		vavg = vavg + res[2]
		ratio = ratio + res[3]
		energy = energy + res[4]
		power = power + res[5]
		vmin = vmin + res[8]
		vmax = vmax + res[7]
		eff = eff + res[6]

		v0 = res[1]

		CNCheck = CNCheck + res[9]
		RIME = RIME + res[10]
		RIMcheck = RIMcheck + res[11]
		RIMcatchE = RIMcatchE + res[12]
		RIMcatchD = RIMcatchD + res[13]
		catchacceler = catchacceler + res[14]
		DragEff = DragEff + res[15]
		tcatchacceler = res[14]
	}

	faantal := float64(aantal)

	dv = dv / faantal
	vend = vend / faantal
	vavg = vavg / faantal
	vmin = vmin / faantal
	vmax = vmax / faantal
	ratio = ratio / faantal
	energy = energy / faantal
	power = power / faantal
	eff = eff / faantal
	CNCheck = CNCheck / faantal
	RIME = RIME / faantal
	RIMcheck = RIMcheck / faantal
	RIMcatchE = RIMcatchE / faantal
	RIMcatchD = RIMcatchD / faantal
	DragEff = DragEff / faantal
	catchacceler = catchacceler / faantal

	return []float64{
		dv,           // 0
		vend,         // 1
		vavg,         // 2
		ratio,        // 3
		energy,       // 4
		power,        // 5
		eff,          // 6
		vmax,         // 7
		vmin,         // 8
		CNCheck,      // 9
		RIME,         // 10
		RIMcheck,     // 11
		RIMcatchE,    // 12
		RIMcatchD,    // 13
		catchacceler, // 14
		DragEff,      // 15
	}, nil
}

// ConstantVeloFast calculates force and power to achieve certain boat speed
func ConstantVeloFast(
	velo float64,
	crew *Crew,
	rigging *Rig,
	timestep float64,
	aantal int,
	aantal2 int,
	Fmin float64,
	Fmax float64,
	catchacceler float64,
	windv float64,
	dowind bool,
) ([]float64, error) {

	F, err := LinSpace(Fmin, Fmax, aantal)
	if err != nil {
		return nil, err
	}
	velocity := make([]float64, aantal)

	ca := catchacceler
	dv := 1.
	vend := 4.0

	for i, ff := range F {
		for dv/vend > 0.001 {
			res, err := EnergyBalance(ff, crew, rigging, vend, timestep, ca, windv, dowind)
			if err != nil {
				return nil, err
			}
			dv = res[0]
			vend = res[1]
			ca = res[14]
		}
		res, err := Stroke(ff, crew, rigging, vend, timestep, 10, ca, dowind, windv)
		if err != nil {
			return nil, err
		}
		velocity[i] = res[2]
	}

	fres, _ := srinterpol3(F, velocity, velo)
	// fmt.Println(fres, velocity, F)
	dv = 1.0

	for dv/vend > 0.001 {
		res, err := EnergyBalance(fres, crew, rigging, vend, timestep, ca, windv, dowind)
		if err != nil {
			return nil, err
		}
		dv = res[0]
		vend = res[1]
		ca = res[14]
	}

	res, err := Stroke(fres, crew, rigging, vend, timestep, 10, ca, dowind, windv)
	if err != nil {
		return nil, err
	}
	vavg := res[2]
	ratio := res[3]
	pw := res[5]
	eff := res[6]

	return []float64{
		fres,  // 0
		vavg,  // 1
		ratio, // 2
		pw,    // 3
		eff,   // 4
	}, nil

}

// ConstantWattFast returns force, average speed given input power in Watt
func ConstantWattFast(
	watt float64,
	crew *Crew,
	rigging *Rig,
	timestep float64,
	aantal int,
	aantal2 int,
	Fmin float64,
	Fmax float64,
	catchacceler float64,
	windv float64,
	dowind bool,
	maxIterationsAllowed int,
) ([]float64, error) {

	F, err := LinSpace(Fmin, Fmax, aantal)
	if err != nil {
		return nil, err
	}
	power := make([]float64, aantal)

	ca := catchacceler
	dv := 1.
	vend := 4.0

	for i, ff := range F {
		for dv/vend > 0.001 {
			res, err := EnergyBalance(ff, crew, rigging, vend, timestep, ca, windv, dowind)
			if err != nil {
				return res, err
			}
			dv = res[0]
			vend = res[1]
			ca = res[14]
		}
		res, err := Stroke(ff, crew, rigging, vend, timestep, 10, ca, dowind, windv)
		if err != nil {
			return nil, err
		}
		power[i] = res[5]
	}

	fres, _ := srinterpol3(F, power, watt)

	dv = 1.0

	for count := 0; dv/vend > 0.001 && count <= maxIterationsAllowed; count++ {
		res, err := EnergyBalance(fres, crew, rigging, vend, timestep, ca, windv, dowind)
		if err != nil {
			dv = 0.01 * vend
			vend += dv
			continue
		}
		dv = res[0]
		vend = res[1]
		ca = res[14]
	}

	res, err := Stroke(fres, crew, rigging, vend, timestep, 10, ca, dowind, windv)
	if err != nil {
		return nil, err
	}
	vavg := res[2]
	ratio := res[3]
	pw := res[5]
	eff := res[6]

	return []float64{
		fres,  // 0
		vavg,  // 1
		ratio, // 2
		pw,    // 3
		eff,   // 4
	}, nil
}

func tailwind(
	bearing float64,
	vwind float64,
	winddirection float64,
	vstream float64,
) float64 {
	b := bearing * math.Pi / 180.
	w := winddirection * math.Pi / 180.

	vtail := -vwind*math.Cos(w-b) - vstream

	return vtail
}

// PhysGetPower Gets power and no wind pace
func PhysGetPower(
	velo float64,
	rower *Crew,
	rigging *Rig,
	bearing float64,
	vwind float64,
	winddirection float64,
	vstream float64,
) ([]float64, error) {

	tw := tailwind(bearing, vwind, winddirection, vstream)
	velowater := velo - vstream
	res, err := ConstantVeloFast(velowater, rower, rigging, 0.03, 5, 5, 50, 600, 5, tw, true)
	if err != nil {
		return nil, err
	}

	force := res[0]
	power := res[3]
	ratio := res[2]

	res2, err := ConstantWattFast(power, rower, rigging, 0.03, 5, 5, 50, 600, 5, 0, true, 15)
	if err != nil {
		return nil, err
	}
	pnowind := 500. / res2[1]

	return []float64{
		power,      // 0
		ratio,      // 1
		force,      // 2
		pnowind,    // 3
		math.NaN(), // 4
	}, nil
}
