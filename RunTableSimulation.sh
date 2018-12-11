# $1 is gamma, tau is scanned from 0.5001 to 3.0

output_file=Table1_$1_$2.out
if [ -e "${output_file}" ]; then
	rm ${output_file};
fi;
echo "Scanning parameters for the $1 scheme with gamma $2."
./bin/poiseuille_$1_stability --gamma $2 --tau 0.5001 --Ny 41 >> ${output_file}
parallel -j8 './bin/poiseuille_{1}_stability --gamma {2} --tau $(expr {3} / 1000).$(expr {3} % 1000 | xargs printf "%03d") --Ny 41' >> ${output_file} ::: $1 ::: $2 ::: {501..3000}
