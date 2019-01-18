sudo perf timechart record ../pdc_mini_aevol -n 50 &&
sudo perf timechart &&
gwenview output.svg & 
