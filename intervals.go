package gorow

import "fmt"

func cumulativedistance(strokes []StrokeRecord) ([]float64, []float64, []float64, error) {
	var cumdist []float64
	var dps []float64
	var deltas []float64
	cumdist = append(cumdist, strokes[0].Distance)
	dps = append(dps, 0)
	deltas = append(deltas, strokes[1].Timestamp-strokes[0].Timestamp)
	curdist := strokes[0].Distance

	for i := 1; i < len(strokes); i++ {
		delta := strokes[i].Distance - strokes[i-1].Distance
		distanceperstroke := 0.0
		deltat := strokes[i].Timestamp - strokes[i-1].Timestamp
		newdelta := delta
		if delta < 0 {
			newdelta = 0.0
			// if i >= 2 {
			// 	t3 := strokes[i].Timestamp
			// 	t2 := strokes[i-1].Timestamp
			// 	t1 := strokes[i-2].Timestamp
			//
			// 	d2 := strokes[i-1].Distance
			// 	d1 := strokes[i-2].Distance
			// 	if t2 > t1 {
			// 		newdelta = (t3 - t2) * (d2 - d1) / (t2 - t1)
			// 	}
			//}
		}
		if deltat > 0 && strokes[i].Spm > 0 {
			distanceperstroke = 60. * newdelta / (deltat * strokes[i].Spm)
		}
		curdist += newdelta
		cumdist = append(cumdist, curdist+newdelta)
		dps = append(dps, distanceperstroke)
		deltas = append(deltas, deltat)
	}
	return cumdist, dps, deltas, nil
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
	cumdist, dps, deltas, err := cumulativedistance(strokes)
	if err != nil {
		return "", err
	}
	currentstate := workstate

	deltatotwork := 0.0
	deltatotrest := 0.0
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
				workpower += deltas[i] * stroke.Power
				workspm += deltas[i] * stroke.Spm
				workavghr += deltas[i] * stroke.Hr
				if stroke.Hr > workmaxhr {
					workmaxhr = stroke.Hr
				}

				workdps += deltas[i] * dps[i]
				deltatotwork += deltas[i]

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
				restpower += deltas[i] * stroke.Power
				restspm += deltas[i] * stroke.Spm
				restavghr += deltas[i] * stroke.Hr
				if stroke.Hr > restmaxhr {
					restmaxhr = stroke.Hr
				}
				restdps += deltas[i] * dps[i]
				deltatotrest += deltas[i]
			}
		}
	}

	avgpacework := 500. * worktime / workmeters
	avgpacerest := 500. * resttime / restmeters

	workpower /= deltatotwork
	restpower /= deltatotrest
	workspm /= deltatotwork
	restspm /= deltatotrest
	workavghr /= deltatotwork
	restavghr /= deltatotrest

	workdps /= deltatotwork
	restdps /= deltatotrest

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

func lapstats(strokes []StrokeRecord, lap int, separator string) (string, error) {
	const workstate = "work"
	const reststate = "rest"
	cumdist, dps, deltas, err := cumulativedistance(strokes)
	if err != nil {
		return "", err
	}
	currentstate := workstate

	if strokes[0].Workoutstate == 3 {
		currentstate = reststate
	}

	deltatot := 0.0
	restmetersprevious := 0.0
	restmeters := 0.0
	resttimeprevious := 0.0
	resttime := 0.0

	workmetersprevious := 0.0
	workmeters := 0.0
	if currentstate == workstate {
		workmeters = strokes[0].Distance
	}
	worktimeprevious := strokes[0].Timestamp
	worktime := 0.0

	workpower := 0.0

	workspm := 0.0

	workavghr := 0.0

	workmaxhr := 0.0

	workdps := 0.0

	countintervalworkstrokes := 0

	for i, stroke := range strokes {
		strokestate := stroke.Workoutstate
		lapid := stroke.Lapnr
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
				if int(lapid) == lap {
					countintervalworkstrokes++
					workpower += deltas[i] * stroke.Power
					workspm += deltas[i] * stroke.Spm
					workavghr += deltas[i] * stroke.Hr
					if stroke.Hr > workmaxhr {
						workmaxhr = stroke.Hr
					}
					workdps += deltas[i] * dps[i]
					deltatot += deltas[i]
				}

			}
		case 3:
			{
				if currentstate == workstate {
					if int(lapid) == lap {
						workmeters += cumdist[i] - workmetersprevious
						worktime += stroke.Timestamp - worktimeprevious
					}

					workmetersprevious = cumdist[i]
					restmetersprevious = cumdist[i]
					resttimeprevious = stroke.Timestamp
					worktimeprevious = stroke.Timestamp
					currentstate = reststate
				}

			}
		}
	}

	workpower /= deltatot
	workspm /= deltatot
	workavghr /= deltatot
	workdps /= deltatot

	if worktime > 0 {
		workdps = (workmeters / worktime) * 60. / workspm
	}

	workmeters -= float64(int(strokes[0].Distance))
	workmeters = float64(int(workmeters))
	avgpacework := 500. * worktime / workmeters

	stri := fmt.Sprintf("%02d%s%05.0f%s", lap, separator, workmeters, separator)
	stri += fmt.Sprintf("%s%s%s%s", formatPace(worktime), separator, formatPace(avgpacework), separator)
	stri += fmt.Sprintf("%05.1f%s%04.1f%s", workpower, separator, workspm, separator)
	stri += fmt.Sprintf("%05.1f%s%05.1f%s%04.1f", workavghr, separator, workmaxhr, separator, workdps)
	stri += "\n"
	return stri, nil
}

func intervalstats(strokes []StrokeRecord, separator string) (string, error) {
	stri := "Workout Details\n"
	stri += "#-|SDist|-Split-|-SPace-|-Pwr-|SPM-|AvgHR|MaxHR|DPS-\n"
	laps, err := GetLapNumbers(strokes)
	if err != nil {
		return stri, err
	}

	for lap := range laps {
		lapstri, err := lapstats(strokes, lap, separator)
		if err != nil {
			return stri, err
		}
		stri += lapstri
	}

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

	var avgpace, avgspm, avghr, maxhr, avgdps, avgpower, deltatot float64
	// var avgv float64

	cumdist, _, deltas, _ := cumulativedistance(strokes)
	totaldist := cumdist[len(strokes)-1]
	totaltime := strokes[len(strokes)-1].Timestamp - strokes[0].Timestamp

	for i, stroke := range strokes {
		// avgv += deltas[i] * 500. / stroke.Pace
		avgspm += deltas[i] * stroke.Spm
		avghr += deltas[i] * stroke.Hr
		avgdps += deltas[i] * stroke.Strokedistance
		avgpower += deltas[i] * stroke.Power
		if stroke.Hr > maxhr {
			maxhr = stroke.Hr
		}
		deltatot += deltas[i]

	}

	avgspm /= deltatot
	avghr /= deltatot
	avgpower /= deltatot

	avgv := totaldist / totaltime
	avgpace = 500. / avgv
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
