all: Couette_circle_convex Couette_circle_nonconvex

Couette_circle_convex: Couette_circle_convex.cpp
	g++ Couette_circle_convex.cpp -o Couette_circle_convex -I ../include

Couette_circle_nonconvex: Couette_circle_nonconvex.cpp
	g++ Couette_circle_nonconvex.cpp -o Couette_circle_nonconvex -I ../include

clean:
	rm Couette_circle_convex Couette_circle_nonconvex
