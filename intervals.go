package gorow

import "fmt"

func cumulativedistance(strokes []StrokeRecord) ([]float64, []float64, error) {
	var cumdist []float64
	var dps []float64
	cumdist = append(cumdist, 0)
	dps = append(dps, 0)
	curdist := 0.0
	for i := 1; i < len(strokes); i++ {
		delta := strokes[i].Distance - strokes[i-1].Distance
		distanceperstroke := 0.0
		deltat := strokes[i].Timestamp - strokes[i-1].Timestamp
		if deltat > 0 && strokes[i].Spm > 0 {
			distanceperstroke = 60. * delta / (deltat * strokes[i].Spm)
		}
		if delta > 0 {
			curdist += delta
			cumdist = append(cumdist, curdist)
			dps = append(dps, distanceperstroke)
		} else {
			cumdist = append(cumdist, curdist)
			dps = append(dps, 0)
		}
	}
	return cumdist, dps, nil
}

// SetLapNumbers returns workout data with lap number updated
func SetLapNumbers(strokes *[]StrokeRecord) error {
	lapnumber := 0
	const reststate = "rest"
	const workstate = "work"

	currentstate := reststate

	for i := range *strokes {
		if currentstate == reststate {
			strokestate := (*strokes)[i].Workoutstate
			switch strokestate {
			case 1, 4, 5, 8, 9, 6, 7:
				{
					currentstate = workstate
					lapnumber++
				}
			}
		}
		if currentstate == workstate {
			strokestate := (*strokes)[i].Workoutstate
			switch strokestate {
			case 3:
				{
					currentstate = reststate
				}
			}
		}
		(*strokes)[i].Lapnr = int64(lapnumber)
	}

	return nil
}

func workrest(strokes []StrokeRecord, separator string) (string, error) {
	const workstate = "work"
	const reststate = "rest"
	cumdist, dps, err := cumulativedistance(strokes)
	if err != nil {
		return "", err
	}
	currentstate := workstate

	restmetersprevious := 0.0
	restmeters := 0.0
	resttimeprevious := 0.0
	resttime := 0.0

	workmetersprevious := 0.0
	workmeters := 0.0
	worktimeprevious := strokes[0].Timestamp
	worktime := 0.0

	workpower := 0.0
	restpower := 0.0

	restspm := 0.0
	workspm := 0.0

	restavghr := 0.0
	workavghr := 0.0
	restmaxhr := 0.0
	workmaxhr := 0.0

	restdps := 0.0
	workdps := 0.0

	for i, stroke := range strokes {
		strokestate := stroke.Workoutstate
		switch strokestate {
		case 1, 4, 5, 8, 9, 6, 7:
			{
				if currentstate == reststate {
					// switching to work
					restmeters += cumdist[i] - restmetersprevious
					resttime += stroke.Timestamp - resttimeprevious
					restmetersprevious = cumdist[i]
					workmetersprevious = cumdist[i]
					resttimeprevious = stroke.Timestamp
					worktimeprevious = stroke.Timestamp
					currentstate = workstate
				}
				// averages
				workpower += stroke.Power
				workspm += stroke.Spm
				workavghr += stroke.Hr
				if stroke.Hr > workmaxhr {
					workmaxhr = stroke.Hr
				}

				workdps += dps[i]

			}
		case 3:
			{
				if currentstate == workstate {
					workmeters += cumdist[i] - workmetersprevious
					worktime += stroke.Timestamp - worktimeprevious
					workmetersprevious = cumdist[i]
					restmetersprevious = cumdist[i]
					resttimeprevious = stroke.Timestamp
					worktimeprevious = stroke.Timestamp
					currentstate = reststate
				}
				// averages
				restpower += stroke.Power
				restspm += stroke.Spm
				restavghr += stroke.Hr
				if stroke.Hr > restmaxhr {
					restmaxhr = stroke.Hr
				}
				restdps += dps[i]
			}
		}
	}

	avgpacework := 500. * worktime / workmeters
	avgpacerest := 500. * resttime / restmeters

	workpower /= float64(len(strokes))
	restpower /= float64(len(strokes))
	workspm /= float64(len(strokes))
	restspm /= float64(len(strokes))
	workavghr /= float64(len(strokes))
	restavghr /= float64(len(strokes))

	workdps /= float64(len(strokes))
	restdps /= float64(len(strokes))

	if worktime > 0 {
		workdps = (workmeters / worktime) * 60. / workspm
	}
	if resttime > 0 {
		restdps = (restmeters / resttime) * 60. / restspm
	}

	striw, err := workstring(workmeters, worktime, avgpacework, workspm, workavghr, workmaxhr, workdps, workpower, separator, "W")
	if err != nil {
		return striw, err
	}

	strir, err := workstring(restmeters, resttime, avgpacerest, restspm, restavghr, restmaxhr, restdps, restpower, separator, "R")
	if err != nil {
		return striw, err
	}

	return striw + strir, nil
}

func intervalstats(strokes []StrokeRecord, separator string) (string, error) {
	stri := "Workout Details\n"
	stri += "#-|SDist|-Split-|-SPace-|-Pwr-|SPM-|AvgHR|MaxHR|DPS-\n"
	return stri, nil
}

// AllStats returns all stats
func AllStats(strokes []StrokeRecord, title string, separator string) (string, error) {
	summary, err := SummaryString(strokes, title, separator)
	if err != nil {
		return summary, err
	}

	workrest, err := workrest(strokes, separator)
	if err != nil {
		return summary + workrest, err
	}

	intervals, err := intervalstats(strokes, separator)
	if err != nil {
		return summary + workrest, err
	}

	allstats := summary + workrest + intervals

	return allstats, nil
}

// GetLapNumbers returns a list of unique lap numbers
func GetLapNumbers(strokes []StrokeRecord) ([]int64, error) {
	keys := make(map[int64]bool)
	list := []int64{}
	for _, stroke := range strokes {
		entry := stroke.Lapnr
		if _, value := keys[entry]; !value {
			keys[entry] = true
			list = append(list, entry)
		}
	}
	return list, nil
}

func divmodfloat(numerator float64, denominator int64) (quotient int64, remainder float64) {
	quotient = int64(numerator) / denominator // integer division, decimals are truncated
	remainder = numerator - float64(quotient)*float64(denominator)
	return
}

func formatTime(duration float64) string {
	min, sec := divmodfloat(duration, 60)
	if min > 60 {
		hour := min / 60
		min = min % 60
		str1 := fmt.Sprintf("%02d:%02d:%04.1f", hour, min, sec)
		return str1
	}
	str1 := fmt.Sprintf("00:%02d:%04.1f", min, sec)
	return str1

}

func formatPace(seconds float64) string {
	min, sec := divmodfloat(seconds, 60)
	return fmt.Sprintf("%02d:%04.1f", min, sec)
}

func workstring(
	totaldist, totaltime, avgpace, avgspm, avghr, maxhr, avgdps, avgpower float64,
	separator, symbol string,
) (string, error) {
	stri1 := fmt.Sprintf("%2s%s%05.0f%s%s%s", symbol, separator, totaldist, separator, formatTime(totaltime), separator)
	stri1 += fmt.Sprintf("%s%s", formatPace(avgpace), separator)
	stri1 += fmt.Sprintf("%05.1f%s%04.1f%s%05.1f%s", avgpower, separator, avgspm, separator,
		avghr, separator)

	stri1 += fmt.Sprintf("%05.1f%s%04.1f\n", maxhr, separator, avgdps)

	return stri1, nil
}

func intervalstring(
	nr, totaldist, totaltime, avgpace, avgspm, avghr, maxhr, avgdps, avgpower float64,
	separator string,
) (string, error) {
	stri1 := ""
	stri1 += fmt.Sprintf("%02.0f%s%05.0f%s%s%s", nr, separator, totaldist, separator, formatPace(totaltime), separator)
	stri1 += fmt.Sprintf("%s%s%05.1f%s%04.1f%s%05.1f", formatPace(avgpace), separator, avgpower, separator, avgspm, separator, avghr)
	stri1 += fmt.Sprintf("%s%3.1f%s%04.1f\n", separator, maxhr, separator, avgdps)

	return stri1, nil
}

// SummaryString writes summary
func SummaryString(strokes []StrokeRecord, title string, separator string) (string, error) {
	stri1 := fmt.Sprintf("Workout Summary - %s\n", title)
	stri1 += fmt.Sprintf("--%[1]sTotal%[1]s-Total----%[1]s--Avg--%[1]s-Avg-%[1]sAvg-%[1]s-Avg-%[1]s-Max-%[1]s-Avg\n", separator)
	stri1 += fmt.Sprintf("--%[1]sDist-%[1]s-Time-----%[1]s--Pace-%[1]s-Pwr-%[1]sSPM-%[1]s-HR--%[1]s-HR--%[1]s-DPS\n", separator)

	var avgpace, avgspm, avghr, maxhr, avgdps, avgpower float64
	var avgv float64

	cumdist, _, _ := cumulativedistance(strokes)
	totaldist := cumdist[len(strokes)-1]
	totaltime := strokes[len(strokes)-1].Timestamp - strokes[0].Timestamp

	for _, stroke := range strokes {
		avgv += 500. / stroke.Pace
		avgspm += stroke.Spm
		avghr += stroke.Hr
		avgdps += stroke.Strokedistance
		avgpower += stroke.Power
		if stroke.Hr > maxhr {
			maxhr = stroke.Hr
		}

	}

	avgv /= float64(len(strokes))
	avgpace = 500. / avgv
	avgspm /= float64(len(strokes))
	avghr /= float64(len(strokes))
	avgpower /= float64(len(strokes))

	avgdps = totaldist / (totaltime * avgspm / 60.)

	stri1 += fmt.Sprintf("--%s%05.0f%s", separator, totaldist, separator)
	stri1 += fmt.Sprintf("%s%s%s", formatTime(totaltime), separator, formatPace(avgpace))
	stri1 += fmt.Sprintf("%s%05.1f", separator, avgpower)
	stri1 += fmt.Sprintf("%s%2.1f%s%05.1f", separator, avgspm, separator, avghr)
	stri1 += fmt.Sprintf("%s%05.1f%s%04.1f\n", separator, maxhr, separator, avgdps)

	return stri1, nil
}

// UpdateIntervalMetric updates intervals per metric
func UpdateIntervalMetric(
	strokes []StrokeRecord,
	metric string,
	setvalue float64,
	mode string,
	smoothwindow float64,
	activewindow [2]float64,
) ([]StrokeRecord, error) {
	// Set Active Window
	if activewindow[1] == 0 {
		activewindow[1] = strokes[len(strokes)-1].Timestamp

	}
	activewindow[0] -= strokes[0].Timestamp
	activewindow[1] -= strokes[0].Timestamp

	var metricvalues []float64
	var dtavg float64

	/// First loop to get some stats
	for i, stroke := range strokes {
		value, err := stroke.GetField(metric)
		if err != nil {
			return strokes, fmt.Errorf("Could not get value for %s in record number %v", metric, i)
		}
		metricvalues = append(metricvalues, value)
		if i >= 1 {
			dtavg += strokes[i].Timestamp - strokes[i-1].Timestamp
		}
	}
	dtavg /= float64(len(strokes))

	nrrecords := uint(smoothwindow / dtavg)
	metricvalues, _ = ewmovingaverage(metricvalues, nrrecords)

	// Second loop to set the values

	largerthantype := 5.
	smallerthantype := 3.
	if mode == "smaller" {
		largerthantype = 3.
		smallerthantype = 5.
	}

	timezero := strokes[0].Timestamp

	for i, stroke := range strokes {
		stroketime := stroke.Timestamp - timezero
		if metricvalues[i] >= setvalue && stroketime >= activewindow[0] && stroketime <= activewindow[1] {
			strokes[i].Workoutstate = largerthantype
		} else {
			strokes[i].Workoutstate = smallerthantype
		}
	}

	return strokes, nil
}
