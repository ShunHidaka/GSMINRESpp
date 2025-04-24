CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Wpedantic -O3
#LAFLAGS  = -lm -lgfortran -lblas -llapack
LAFLAGS  = -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas
ifdef debug
	CXXFLAGS = -std=c++17 -Wall -Wextra -Wconversion -fsanitize=address -fsanitize=undefined -Wpedantic -O0 -g3
endif

TARGET = sample1.out sample2.out
all: $(TARGET)
sample1.out: sample1.cpp gsminres_solver.cpp gsminres_util.cpp
	$(CXX) $(CXXFLAGS) $^ $(LAFLAGS) -o $@
sample2.out: sample2.cpp gsminres_solver.cpp gsminres_util.cpp
	$(CXX) $(CXXFLAGS) $^ $(LAFLAGS) -o $@
clean:
	rm -f *.out
