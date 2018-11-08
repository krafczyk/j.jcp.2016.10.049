jobs=4
if [ $# -gt 0 ]; then
	jobs=$1
fi;

echo "Left Figures"
if [ ! -e "Figure6_nonconvex.out" ]; then
	echo "\"Lattice Size\", \"NY\", \"Tau\", \"Error\"" > Figure6_nonconvex.out
	parallel -j${jobs} --colsep '\t' ./bin/taylor_vortex_nonconvex --tau {1} --Ny {2} >> Figure6_nonconvex.out :::: Figure6_nonconvex.pars
fi;

echo "Right Figures"
if [ ! -e "Figure6_convex.out" ]; then
	echo "\"Lattice Size\", \"NY\", \"Tau\", \"Error\"" > Figure6_convex.out
	parallel -j${jobs} --colsep '\t' ./bin/taylor_vortex_convex --tau {1} --Ny {2} >> Figure6_convex.out :::: Figure6_convex.pars
fi;
