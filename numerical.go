package gorow

import (
	"github.com/cnkei/gospline"
	//	"github.com/pa-m/sklearn/interpolate"
	"github.com/sgreben/piecewiselinear"
)

func cumsum(x []float64, m float64) []float64 {
	var y = make([]float64, len(x))
	for i := 1; i < len(x); i++ {
		y[i] = y[i-1] + m*x[i]
	}
	return y
}

func srinterpol1(x []float64, y []float64, target float64) float64 {
	// f := interpolate.Interp1d(x, y)
	f := piecewiselinear.Function{X: x, Y: y}

	var newx = LinSpace(x[0], x[len(x)-1], 100*len(x))
	var newy = LinSpace(x[0], x[len(x)-1], 100*len(x))

	var minysq = target
	var minindex = 0

	for i := range newx {
		// newy[i] = f(newx[i])
		newy[i] = f.At(newx[i])
		var ysq = (target - newy[i]) * (target - newy[i])
		if ysq < minysq {
			minysq = ysq
			minindex = i
		}
	}
	return newx[minindex]
}

func srinterpol2(x []float64, y []float64, target float64) float64 {
	var dx = x[1] - x[0]

	var newx = LinSpace(x[0], x[len(x)-1], 10*len(x))

	s := gospline.NewCubicSpline(x, y)

	var newy = s.Range(x[0], x[len(x)-1], dx/10.)

	var minysq = target
	var minindex = 0

	for i, val := range newy {
		var ysq = (target - val) * (target - val)
		if ysq < minysq {
			minysq = ysq
			minindex = i
		}
	}

	return newx[minindex]
}
