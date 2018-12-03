#CXXFLAGS := -g
CXXFLAGS := -O3 -std=c++11

all: bin bin/Couette_circle_convex bin/Couette_circle_nonconvex bin/taylor_vortex_convex bin/taylor_vortex_nonconvex bin/poiseuille_convex bin/poiseuille_convex_stability bin/poiseuille_nonconvex bin/poiseuille_nonconvex_stability

bin:
	@mkdir -p bin

bin/Couette_circle_convex: Couette/Couette_circle_convex.cpp
	g++ $(CXXFLAGS) Couette/Couette_circle_convex.cpp -o bin/Couette_circle_convex -I include

bin/Couette_circle_nonconvex: Couette/Couette_circle_nonconvex.cpp
	g++ $(CXXFLAGS) Couette/Couette_circle_nonconvex.cpp -o bin/Couette_circle_nonconvex -I include

bin/taylor_vortex_convex: Vortex/taylor_vortex_convex.cpp
	g++ $(CXXFLAGS) Vortex/taylor_vortex_convex.cpp -o bin/taylor_vortex_convex -I include

bin/taylor_vortex_nonconvex: Vortex/taylor_vortex_nonconvex.cpp
	g++ $(CXXFLAGS) Vortex/taylor_vortex_nonconvex.cpp -o bin/taylor_vortex_nonconvex -I include

bin/poiseuille_convex: Poiseuille/poiseuille_convex.cpp
	g++ $(CXXFLAGS) Poiseuille/poiseuille_convex.cpp -o bin/poiseuille_convex -I include

bin/poiseuille_convex_stability: Poiseuille/poiseuille_convex_stability.cpp
	g++ $(CXXFLAGS) Poiseuille/poiseuille_convex_stability.cpp -o bin/poiseuille_convex_stability -I include

bin/poiseuille_nonconvex: Poiseuille/poiseuille_nonconvex.cpp
	g++ $(CXXFLAGS) Poiseuille/poiseuille_nonconvex.cpp -o bin/poiseuille_nonconvex -I include

bin/poiseuille_nonconvex_stability: Poiseuille/poiseuille_nonconvex_stability.cpp
	g++ $(CXXFLAGS) Poiseuille/poiseuille_nonconvex_stability.cpp -o bin/poiseuille_nonconvex_stability -I include

clean:
	@rm -r bin
