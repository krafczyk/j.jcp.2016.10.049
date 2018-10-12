echo "Left Figures"
echo "\"Lattice Size\", \"NY\", \"Gamma\", \"Tau\", \"Error\"" > Figure3_nonconvex.dat
parallel -j4 --colsep '\t' ./Poiseuille/poiseuille_nonconvex --gamma {1} --tau {2} --Ny {3} >> Figure3_nonconvex.dat :::: Figure3_nonconvex_pars.dat
