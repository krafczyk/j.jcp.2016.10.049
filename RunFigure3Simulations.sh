echo "Left Figures"
if [ ! -e "Figure3_nonconvex.out" ]; then
	echo "\"Lattice Size\", \"NY\", \"Gamma\", \"Tau\", \"Error\"" > Figure3_nonconvex.out
	parallel -j4 --colsep '\t' ./bin/poiseuille_nonconvex --gamma {1} --tau {2} --Ny {3} >> Figure3_nonconvex.out :::: Figure3_nonconvex.pars
fi;

echo "Right Figures"
if [ ! -e "Figure3_convex.out" ]; then
	echo "\"Lattice Size\", \"NY\", \"Gamma\", \"Tau\", \"Error\"" > Figure3_convex.out
	parallel -j4 --colsep '\t' ./bin/poiseuille_convex --gamma {1} --tau {2} --Ny {3} >> Figure3_convex.out :::: Figure3_convex.pars
fi;
