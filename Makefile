#CXXFLAGS := -g
CXXFLAGS := -O3

all: bin/Couette_circle_convex bin/Couette_circle_nonconvex bin/taylor_vortex_convex bin/taylor_vortex_nonconvex bin/poiseuille_convex bin/poiseuille_nonconvex

bindir:
	@mkdir -p bin

bin/Couette_circle_convex: bindir Couette/Couette_circle_convex.cpp
	g++ $(CXXFLAGS) Couette/Couette_circle_convex.cpp -o bin/Couette_circle_convex -I include

bin/Couette_circle_nonconvex: bindir Couette/Couette_circle_nonconvex.cpp
	g++ $(CXXFLAGS) Couette/Couette_circle_nonconvex.cpp -o bin/Couette_circle_nonconvex -I include

bin/taylor_vortex_convex: bindir Vortex/taylor_vortex_convex.cpp
	g++ $(CXXFLAGS) Vortex/taylor_vortex_convex.cpp -o bin/taylor_vortex_convex -I include

bin/taylor_vortex_nonconvex: bindir Vortex/taylor_vortex_nonconvex.cpp
	g++ $(CXXFLAGS) Vortex/taylor_vortex_nonconvex.cpp -o bin/taylor_vortex_nonconvex -I include

bin/poiseuille_convex: bindir Poiseuille/poiseuille_convex.cpp
	g++ $(CXXFLAGS) Poiseuille/poiseuille_convex.cpp -o bin/poiseuille_convex -I include

bin/poiseuille_nonconvex: bindir Poiseuille/poiseuille_nonconvex.cpp
	g++ $(CXXFLAGS) Poiseuille/poiseuille_nonconvex.cpp -o bin/poiseuille_nonconvex -I include

clean:
	@rm -r bin
