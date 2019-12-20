package gorow

import (
	"github.com/cnkei/gospline"
	//	"github.com/pa-m/sklearn/interpolate"
	"errors"

	"github.com/sgreben/piecewiselinear"
)

func cumsum(x []float64, m float64) []float64 {
	var y = make([]float64, len(x))
	for i := 1; i < len(x); i++ {
		y[i] = y[i-1] + m*x[i]
	}
	return y
}

func srinterpol3(x []float64, y []float64, target float64) (float64, error) {
	// assumes x are sorted in increasing order
	aantal := len(x)

	for i := 0; i < aantal-1; i++ {
		if y[i+1] >= target && y[i] < target {
			ratio := (target - y[i]) / (y[i+1] - y[i])
			Xf := x[i] + ratio*(x[i+1]-x[i])
			return Xf, nil
		} else if y[i+1] <= target && y[i] > target {
			ratio := (target - y[i]) / (y[i+1] - y[i])
			Xf := x[i] + ratio*(x[i+1]-x[i])
			return Xf, nil
		}
	}
	return x[aantal-1], errors.New("target outside input range")
}

func srinterpol1(x []float64, y []float64, target float64) float64 {
	// f := interpolate.Interp1d(x, y)
	f := piecewiselinear.Function{X: x, Y: y}

	var newx, _ = LinSpace(x[0], x[len(x)-1], 100*len(x))
	var newy, _ = LinSpace(x[0], x[len(x)-1], 100*len(x))

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

	var newx, _ = LinSpace(x[0], x[len(x)-1], 10*len(x))

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
