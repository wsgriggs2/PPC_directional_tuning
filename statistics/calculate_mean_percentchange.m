function effectsize = calculate_mean_percentchange(analysis_period, baseline)

percentchange = analysis_period - baseline;

effectsize = mean(percentchange);
