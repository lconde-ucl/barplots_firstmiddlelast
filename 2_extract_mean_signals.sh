#!/usr/bin/bash



source ~/.bash_profile 

for sample in human/bigwigs/*bw human/bigwigs_subtract/*bw human/bigwigs_nopseudo/*bw ; do
	for region in human/hs_bedfiles/*bed; do
		output="${sample%.bw}_$(basename $region .bed).txt"
		echo "$sample $region $output"
		bigWigAverageOverBed "$sample" "$region" "$output"
	done
done


for sample in mouse/bigwigs/*bw mouse/bigwigs_subtract/*bw mouse/bigwigs_nopseudo/*bw; do
	for region in mouse/mm_bedfiles/*bed; do
		output="${sample%.bw}_$(basename $region .bed).txt"
		echo "$sample $region $output"
		bigWigAverageOverBed "$sample" "$region" "$output"
	done
done


