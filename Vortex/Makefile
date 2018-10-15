all: taylor_vortex_convex taylor_vortex_nonconvex

taylor_vortex_convex: taylor_vortex_convex.cpp
	g++ taylor_vortex_convex.cpp -o taylor_vortex_convex -I ../include

taylor_vortex_nonconvex: taylor_vortex_nonconvex.cpp
	g++ taylor_vortex_nonconvex.cpp -o taylor_vortex_nonconvex -I ../include

clean: 
	rm taylor_vortex_convex taylor_vortex_nonconvex
