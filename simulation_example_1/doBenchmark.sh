#!/bin/bash
rouge='\e[0;31m'
neutre='\e[0;m'
vert='\e[1;32m'

/bin/bash -c "cd .. ; make" &&
sudo perf timechart record ../pdc_mini_aevol -n ${1} &&
sudo perf timechart -p 'pdc_mini_aevol' --highlight 'pdc_mini_aevol' &&
gwenview output.svg 

#On vérifie que les résultats sont corrects
if [ "0" -eq $(diff stats/stats_simd_mean.csv stats/stats_simd_mean2.csv | wc -w) ]; then 
	echo -e "${vert}RESULTAT CORRECT${neutre}"; 
else 
	echo -e "${rouge}RESULTAT INCORRECT${neutre}";
fi
