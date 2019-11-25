package gorow

import (
	"github.com/cnkei/gospline"
)

func srinterpol1(x []float64, y []float64, target float64) float64 {
	s := gospline.NewCubicSpline(y, x)
	return s.At(target)
}
