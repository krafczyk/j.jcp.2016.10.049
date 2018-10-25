#echo "Left Figures"
#if [ ! -e "Figure4_nonconvex.out" ]; then
#	echo "\"Lattice Size\", \"NY\", \"Gamma\", \"Tau\", \"Error\"" > Figure4_nonconvex.out
#	parallel -j4 --colsep '\t' ./Poiseuille/poiseuille_nonconvex --gamma {1} --tau {2} --Ny {3} >> Figure4_nonconvex.out :::: Figure4_nonconvex.pars
#fi;

echo "Right Figures"
if [ ! -e "Figure4_convex.out" ]; then
	echo "\"Lattice Size\", \"NY\", \"Gamma\", \"Tau\", \"Error\"" > Figure4_convex.out
	parallel -j4 --colsep '\t' ./Poiseuille/poiseuille_convex --gamma {1} --tau {2} --Ny {3} >> Figure4_convex.out :::: Figure4_convex.pars
fi;
