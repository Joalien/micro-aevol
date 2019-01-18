#!/bin/bash
rouge='\e[0;31m'
neutre='\e[0;m'

/bin/bash -c "cd .. ; make" &&
sudo perf timechart record ../pdc_mini_aevol -n ${1} &&
sudo perf timechart -P &&
gwenview output.svg &

#On vérifie que les résultats sont corrects
if ["" -eq $(diff stats/stats_simd_mean.csv stats/stats_simd_mean2.csv)]; then 
	echo "Le résultat est correct"; 
else 
	echo -e "${rouge}RESULTAT INCORRECT${neutre}";
fi
