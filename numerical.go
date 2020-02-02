package gorow

import (
	"math"

	"github.com/cnkei/gospline"
	//	"github.com/pa-m/sklearn/interpolate"
	"errors"

	"github.com/sgreben/piecewiselinear"
)

func reverseslice(in []float64) []float64 {
	out := make([]float64, len(in))

	for i, value := range in {
		out[len(in)-i-1] = value
	}

	return out
}

func ewmovingaverage(in []float64, span uint) ([]float64, error) {
	var out1 = make([]float64, len(in))

	alfa := 2 / (1 + float64(span))

	for i := range in {
		if i < int(span) {
			beta := 1.0
			sumweights := 0.0
			total := 0.0
			for j := 0; j <= i; j++ {
				total += beta * in[i-j]
				sumweights += beta
				beta *= (1 - alfa)
			}
			out1[i] = total / sumweights
			continue
		}
		beta := 1.0
		sumweights := 0.0
		total := 0.0
		for j := 0; j <= int(span); j++ {
			total += beta * in[i-j]
			sumweights += beta
			beta *= (1 - alfa)
		}
		out1[i] = total / sumweights
	}
	return out1, nil
}

func ewmovingaverageright(in []float64, span uint) ([]float64, error) {
	in2 := reverseslice(in)
	out, err := ewmovingaverage(in2, span)
	return reverseslice(out), err
}

func ewmovingaverageboth(in []float64, span uint) ([]float64, error) {
	out1, err := ewmovingaverage(in, span)
	if err != nil {
		return out1, err
	}
	out2, err := ewmovingaverageright(in, span)
	if err != nil {
		return out2, err
	}
	for i := range out1 {
		out1[i] = (out1[i] + out2[i]) / 2.
	}
	return out1, nil
}

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

// X must be monotonously increasing and X2 higher frequency than X
func linearize(X []float64, Y []float64, X2 []float64) ([]float64, error) {
	newy := make([]float64, len(X2))
	for i, xpos := range X2 {
		for j := 0; j < len(X)-1; j++ {
			if X[j] <= xpos && X[j+1] > xpos {
				frac := (xpos - X[j]) / (X[j+1] - X[j])
				newy[i] = Y[j] + frac*(Y[j+1]-Y[j])
			}
		}
	}
	return newy, nil
}

func rolling(data []float64, windowsize int) ([]float64, error) {
	cma := 0.0
	out := make([]float64, len(data))
	for i, value := range data {
		if i < windowsize {
			cma = (float64(i)*cma + value) / (float64(i) + 1)
			out[i] = cma
		} else {
			cma = (float64(windowsize)*cma - data[i-windowsize] + value) / float64(windowsize)
			out[i] = cma
		}
	}
	return out, nil
}

func Sine(x float64) float64 {
	a0 := -0.10132118         // x
	a1 := 0.0066208798        // x^3
	a2 := -0.00017350505      // x^5
	a3 := 0.0000025222919     // x^7
	a4 := -0.000000023317787  // x^9
	a5 := 0.00000000013291342 // x^11

	x2 := x * x

	p11 := a5
	p9 := p11*x2 + a4
	p7 := p9*x2 + a3
	p5 := p7*x2 + a2
	p3 := p5*x2 + a1
	p1 := p3*x2 + a0

	r := (x - math.Pi) * (x + math.Pi) * p1 * x

	return r
}
