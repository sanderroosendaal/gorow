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

// Linspace helper function to create a linear range, like np.linspace
func Linspace(start float64, stop float64, N int) *mat.VecDense {
	rnge := make([]float64, N)
	var step = (stop - start) / float64(N)
	for x := range rnge {
		rnge[x] = start + step*float64(x)
	}
	var r = mat.NewVecDense(N, rnge)
	return r
}

func dragEq(displacement float64, velo float64,
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

type rig struct {
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

func (rg *rig) spread() float64 {
	if rg.roworscull == Scull {
		return (rg.span / 2.)
	}
	return rg._Spread
}

func (rg *rig) overlap() float64 {
	if rg.roworscull == Scull {
		return 2.*rg.lin - rg.span
	}
	return rg.lin - rg._Spread
}

func (rg *rig) buitenhand() float64 {
	if rg.roworscull == Scull {
		return rg.span - 2.*rg.lin*math.Cos(rg.catchangle)
	}
	return rg._Spread - rg.lin*math.Cos(rg.catchangle)
}

func (rg *rig) bladearea() float64 {
	if rg.roworscull == Scull {
		return 2. * rg._Bladearea
	}
	return rg._Bladearea
}

func (rg *rig) dcatch() float64 {
	return (rg.lin * math.Sin(rg.catchangle))
}

func (rg *rig) oarangle(x float64) float64 {
	var dist = rg.dcatch() + x
	var angle = math.Asin(dist / rg.lin)
	return (angle)
}

func newRig(lin float64, mb float64, lscull float64,
	span float64, spread float64, roworscull string,
	catchangle float64, bladearea float64, bladelength float64,
	Nrowers int32, dragform float64) *rig {
	return &rig{
		lin:         lin,
		mb:          mb,
		bladelength: bladelength,
		lscull:      lscull - 0.5*bladelength,
		Nrowers:     Nrowers,
		roworscull:  roworscull,
		span:        span,
		catchangle:  catchangle,
		dragform:    dragform,
	}
}

func bladeForce(oarangle float64, rigging rig, vb, fblade float64) {
	var lin = rigging.lin
	var lscull = rigging.lscull
	var lout = lscull - lin

	var phidot0 = vb * math.Cos(oarangle) / lout
	var phidot = Linspace(phidot0, 2*math.Abs(phidot0), 50)
	var vblade = Linspace(lout, lout, 50)

	vblade.MulElemVec(vblade, phidot)

}
