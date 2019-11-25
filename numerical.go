package gorow

import (
	"math"

	"github.com/cnkei/gospline"
)

func srinterpol1(x []float64, y []float64, target float64) float64 {
	var dx = x[1] - x[0]

	var newx = LinSpace(x[0], x[len(x)-1], 10*len(x))

	s := gospline.NewCubicSpline(x, y)

	var newy = s.Range(x[0], x[len(x)-1], dx/10.)

	var minysq = target
	var minindex = 0

	for i, val := range newy {
		var ysq = math.Pow(target-val, 2)
		if ysq < minysq {
			minysq = ysq
			minindex = i
		}
	}

	return newx[minindex]
}
