# Compiler and Flags
CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3
LAFLAGS  = -lm -lgfortran -lblas -llapack
#LAFLAGS  = -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas
ifdef debug
	CXXFLAGS += -g3 -fsanitize=address -fsanitize=undefined -Wpedantic -O0
endif

# Directories
INCDIR    = include
SRCDIR    = src
SAMPLEDIR = sample
BUILDDIR  = build

# Files
SRCFILES = $(SRCDIR)/gsminres_solver.cpp $(SRCDIR)/gsminres_util.cpp
HEADERS  = $(INCDIR)/gsminres_solver.hpp $(INCDIR)/gsminres_util.hpp $(INCDIR)/gsminres_blas.hpp
TARGETS  = $(BUILDDIR)/sample1.out $(BUILDDIR)/sample2.out

# Targets
all: $(TARGETS)
$(BUILDDIR)/sample1.out: $(SAMPLEDIR)/sample1.cpp $(SRCFILES) $(HEADERS)
	mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $^ $(LAFLAGS) -o $@
$(BUILDDIR)/sample2.out: $(SAMPLEDIR)/sample2.cpp $(SRCFILES) $(HEADERS)
	mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $^ $(LAFLAGS) -o $@

clean:
	rm -rf $(BUILDDIR)
