jobs=4
if [ $# -gt 0 ]; then
	jobs=$1
fi;

if [ ! -e "Figure5.out" ]; then
	echo "\"Lattice Size\", \"NY\", \"Beta\", \"Tau\", \"Error\"" > Figure6_nonconvex.out
	parallel -j${jobs} --colsep '\t' ./bin/Couette_circle_convex --beta {1} --tau 1.0 --Ny 162 --dump-trace "Figure5_{1}.out" :::: Figure5.pars
fi;

