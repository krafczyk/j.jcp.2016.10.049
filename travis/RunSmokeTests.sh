echo "" > smoke_tests.out
echo "poiseuille nonconvex test"
parallel -j1 --colsep '\t' ./bin/poiseuille_nonconvex --gamma {1} --tau {2} --Ny {3} >> smoke_tests.out :::: travis/poiseuille_nonconvex_smoke.pars
echo "poiseuille convex test"
parallel -j1 --colsep '\t' ./bin/poiseuille_convex --gamma {1} --tau {2} --Ny {3} >> smoke_tests.out :::: travis/poiseuille_convex_smoke.pars
echo "taylor vortex nonconvex test"
parallel -j1 --colsep '\t' ./bin/taylor_vortex_nonconvex --tau {1} --Ny {2} >> smoke_tests.out :::: travis/taylor_vortex_nonconvex_smoke.pars
echo "taylor vortex convex test"
parallel -j1 --colsep '\t' ./bin/taylor_vortex_convex --tau {1} --Ny {2} >> smoke_tests.out :::: travis/taylor_vortex_convex_smoke.pars
echo "couette circle convex beta test"
parallel -j1 --colsep '\t' ./bin/Couette_circle_convex --beta {1} --tau 1.0 --Ny 162 --dump-trace "smoke_tests_1.out" >> smoke_tests.out :::: travis/couette_circle_convex_1_smoke.pars
echo "couette circle nonconvex test"
parallel -j1 --colsep '\t' ./bin/Couette_circle_nonconvex --beta 0.5 --tau {1} --Ny {2} >> smoke_tests.out :::: travis/couette_circle_nonconvex_smoke.pars
echo "couette circle convex test"
parallel -j1 --colsep '\t' ./bin/Couette_circle_convex --beta 0.5 --tau {1} --Ny {2} >> smoke_tests.out :::: travis/couette_circle_convex_smoke.pars
