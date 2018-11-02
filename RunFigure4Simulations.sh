echo "Left Figures"
if [ ! -e "Figure4_nonconvex.out" ]; then
	echo "\"Lattice Size\", \"NY\", \"Tau\", \"Error\"" > Figure4_nonconvex.out
	parallel -j4 --colsep '\t' ./bin/taylor_vortex_nonconvex --tau {1} --Ny {2} >> Figure4_nonconvex.out :::: Figure4_nonconvex.pars
fi;

echo "Right Figures"
if [ ! -e "Figure4_convex.out" ]; then
	echo "\"Lattice Size\", \"NY\", \"Tau\", \"Error\"" > Figure4_convex.out
	parallel -j4 --colsep '\t' ./bin/taylor_vortex_convex --tau {1} --Ny {2} >> Figure4_convex.out :::: Figure4_convex.pars
fi;
