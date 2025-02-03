package gorow

import (
	"errors"
	"math"
	//	"fmt"
	"gonum.org/v1/gonum/stat"
)

// WorkoutMetrics calculates workout level metrics (TSS, TRIMP and more)
func WorkoutMetrics(
	filename string,
	ftp float64,
	sex string,
	hrftp, hrmax, hrmin, wpsavg float64,
) (tss, normp, trimp, hrtss, normv, normw, spmtss float64, err error) {
	strokes, err := ReadCSV(filename)

	if err != nil {
		err := errors.New("WorkoutMetrics - ReadCSV" + err.Error())
		return 0, 0, 0, 0, 0, 0, 0, err
	}
	nrstrokes := len(strokes)

	duration := strokes[nrstrokes-1].Timestamp - strokes[0].Timestamp
	time, _ := LinSpace(0, float64(int(duration)), int(duration)+1)

	origpower := make([]float64, len(strokes))
	origtime := make([]float64, len(strokes))
	orighr := make([]float64, len(strokes))
	origwps := make([]float64, len(strokes))
	origvelo := make([]float64, len(strokes))
	origspmpower := make([]float64, len(strokes))
	for i := range strokes {
		origpower[i] = strokes[i].Power
		origtime[i] = strokes[i].Timestamp - strokes[0].Timestamp
		orighr[i] = strokes[i].Hr
		origvelo[i] = strokes[i].Velo
		origwps[i] = strokes[i].Workperstroke
		origspmpower[i] = wpsavg*strokes[i].Spm/60.0 
	}

	if len(strokes) < 30 {
		normp := stat.Mean(origpower, nil)
		
		hrmean := stat.Mean(orighr, nil)
		hrrmean := (hrmean - hrmin) / (hrmax - hrmin )
		f := 1.67
		if sex == "male" {
			f = 1.92
		}
		hrrftp := (hrftp - hrmin) / (hrmax - hrmin)
		trimp1hr := 60. * hrrftp * 0.64 * math.Exp(f*hrrftp)
		trimp := hrrmean * 0.64 * math.Exp(f*hrrmean)
		hrtss = 100. * trimp/trimp1hr
		velomean := stat.Mean(origvelo, nil)
		normv := velomean
		normw := normp

		intensityfactor := normp / ftp

		tss = 100 * ((duration * normp * intensityfactor) / (3600 * ftp))
		
		
		return tss, normp, trimp, hrtss, normv, normw, spmtss, nil
		
		//err := errors.New("WorkoutMetrics - Less than 30 ddata points")
		//return 0, 0, 0, 0, 0, 0, err
	}

	power, err := linearize(origtime, origpower, time)
	if err != nil {
		err := errors.New("WorkoutMetrics - couldn't linearize power" + err.Error())
		return 0, 0, 0, 0, 0, 0, 0, err
	}

	spmpower, err := linearize(origtime, origspmpower, time)
	if err != nil {
		err := errors.New("WorkoutMetrics - couldn't linearize spm power" + err.Error())
		return 0, 0, 0, 0, 0, 0, 0, err
	}

	hr, err := linearize(origtime, orighr, time)
	if err != nil {
		err := errors.New("WorkoutMetrics - couldn't linearize heart rate" + err.Error())
		return 0, 0, 0, 0, 0, 0, 0, err
	}

	wps, err := linearize(origtime, origwps, time)
	if err != nil {
		err := errors.New("WorkoutMetrics - couldn't linearize heart WpS" + err.Error())
		return 0, 0, 0, 0, 0, 0, 0, err
	}

	velo, err := linearize(origtime, origvelo, time)
	if err != nil {
		err := errors.New("WorkoutMetrics - couldn't linearize heart Boat Speed" + err.Error())
		return 0, 0, 0, 0, 0, 0, 0, err
	}

	power, _ = rolling(power, 30)
	hr, _ = rolling(hr, 30)
	wps, _ = rolling(wps, 30)
	velo, _ = rolling(velo, 30)
	spmpower, _ = rolling(spmpower, 30)

	power4 := make([]float64, len(power))
	for i := range power {
		power4[i] = power[i] * power[i] * power[i] * power[i]
	}

	spm4power := make([]float64, len(spmpower))
	for i := range spmpower {
		spm4power[i] = spmpower[i] * spmpower[i] * spmpower[i] * spmpower[i]
	}

	pwr4mean := stat.Mean(power4, nil)
	spm4mean := stat.Mean(spm4power, nil)

	if pwr4mean > 0 {
		normp = math.Pow(pwr4mean, 0.25)
	} else {
		normp = stat.Mean(power, nil)
	}

	normspmpower := 0.0
	if spm4mean > 0 {
		normspmpower = math.Pow(spm4mean, 0.25)
	} else {
		normspmpower = stat.Mean(spmpower, nil)
	}

	intensityfactor := normp / ftp
	intensityfactorspm := normspmpower / ftp

	tss = 100 * ((duration * normp * intensityfactor) / (3600 * ftp))
	spmtss = 100 * ((duration * normspmpower * intensityfactorspm) / (3600 * ftp))

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
	w4mean := stat.Mean(w4, nil)
	v4mean := stat.Mean(v4, nil)
	normw = math.Pow(w4mean, 1./pp)
	normv = math.Pow(v4mean, 1./pp)

	return tss, normp, trimp, hrtss, normv, normw, spmtss, nil
}
