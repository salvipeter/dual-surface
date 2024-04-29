all: dual

LIBGEOM=../libgeom
C0COONS=../c0coons
FLAGS=-std=c++20 -Wall -pedantic -O3
INCLUDES=-I$(LIBGEOM) -I$(C0COONS)
LIBS=$(C0COONS)/c0coons.o $(C0COONS)/curves.o -L$(LIBGEOM)/release -lgeom

dual: dual.cc
	$(CXX) -o $@ $< $(FLAGS) $(INCLUDES) $(LIBS)
