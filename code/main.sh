# Subjects to analyze
subjects=(41 42)

# Compute time-series
python compute_time_series.py "${subjects[@]}"

# Compute instabilities
python compute_instabilities.py "${subjects[@]}"

# Statistical analysis
ipython statistical_analysis.py "${subjects[@]}"
