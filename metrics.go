package gorow

import (
	"errors"
	"math"
)

// WorkoutMetrics calculates workout level metrics (TSS, TRIMP and more)
func WorkoutMetrics(
	filename string,
	ftp float64,
	sex string,
	hrftp, hrmax, hrmin float64,
) (tss, normp, trimp, hrtss, normv, normw float64, err error) {
	strokes, err := ReadCSV(filename)
	if err != nil {
		err := errors.New("WorkoutMetrics - ReadCSV" + err.Error())
		return 0, 0, 0, 0, 0, 0, err
	}
	nrstrokes := len(strokes)

	duration := strokes[nrstrokes-1].timestamp - strokes[0].timestamp
	time, _ := LinSpace(0, duration, int(duration))

	origpower := make([]float64, len(strokes))
	origtime := make([]float64, len(strokes))
	orighr := make([]float64, len(strokes))
	origwps := make([]float64, len(strokes))
	origvelo := make([]float64, len(strokes))
	for i := range strokes {
		origpower[i] = strokes[i].power
		origtime[i] = strokes[i].timestamp - strokes[0].timestamp
		orighr[i] = strokes[i].hr
		origvelo[i] = strokes[i].velo
		origwps[i] = strokes[i].workperstroke
	}

	power, err := linearize(origtime, origpower, time)
	if err != nil {
		err := errors.New("WorkoutMetrics - couldn't linearize power" + err.Error())
		return 0, 0, 0, 0, 0, 0, err
	}

	hr, err := linearize(origtime, orighr, time)
	if err != nil {
		err := errors.New("WorkoutMetrics - couldn't linearize heart rate" + err.Error())
		return 0, 0, 0, 0, 0, 0, err
	}

	wps, err := linearize(origtime, origwps, time)
	if err != nil {
		err := errors.New("WorkoutMetrics - couldn't linearize heart WpS" + err.Error())
		return 0, 0, 0, 0, 0, 0, err
	}

	velo, err := linearize(origtime, origvelo, time)
	if err != nil {
		err := errors.New("WorkoutMetrics - couldn't linearize heart Boat Speed" + err.Error())
		return 0, 0, 0, 0, 0, 0, err
	}

	power, _ = rolling(power, 30)
	hr, _ = rolling(hr, 30)
	wps, _ = rolling(wps, 30)
	velo, _ = rolling(velo, 30)

	power4 := make([]float64, len(power))
	for i := range power {
		power4[i] = power[i] * power[i] * power[i] * power[i]
	}

	pwr4mean := mean(power4)

	if pwr4mean > 0 {
		normp = math.Pow(pwr4mean, 0.25)
	} else {
		normp = mean(power)
	}

	intensityfactor := normp / ftp

	tss = 100 * ((duration * normp * intensityfactor) / (3600 * ftp))

	hrr := make([]float64, len(hr))
	for i := range hr {
		hrr[i] = (hr[i] - hrmin) / (hrmax - hrmin)
	}

	f := 1.67
	if sex == "male" {
		f = 1.92
	}

	hrrftp := (hrftp - hrmin) / (hrmax - hrmin)
	trimp1hr := 60. * hrrftp * 0.64 * math.Exp(f*hrrftp)

	for _, value := range hrr {
		trimp += value * 0.64 * math.Exp(f*value) / 60.
	}

	hrtss = 100. * trimp / trimp1hr

	pp := 8.0

	w4 := make([]float64, len(wps))
	v4 := make([]float64, len(wps))
	for i := range wps {
		w4[i] = math.Pow(wps[i], pp)
		v4[i] = math.Pow(velo[i], pp)
	}
	w4mean := mean(w4)
	v4mean := mean(v4)
	normw = math.Pow(w4mean, 1./pp)
	normv = math.Pow(v4mean, 1./pp)

	return tss, normp, trimp, hrtss, normv, normw, nil
}
