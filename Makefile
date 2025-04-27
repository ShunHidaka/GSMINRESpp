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
BUILDDIR  = bin

# Files
SRCFILES = $(SRCDIR)/gsminres_solver.cpp $(SRCDIR)/gsminres_util.cpp
HEADERS  = $(INCDIR)/gsminres_solver.hpp $(INCDIR)/gsminres_util.hpp $(INCDIR)/gsminres_blas.hpp
TARGETS  = $(BUILDDIR)/sample1.out $(BUILDDIR)/sample2.out

.PHONY: all install

# Default Targets
all: $(TARGETS)
$(BUILDDIR)/sample1.out: $(SAMPLEDIR)/sample1.cpp $(SRCFILES) $(HEADERS)
	mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $^ $(LAFLAGS) -o $@
$(BUILDDIR)/sample2.out: $(SAMPLEDIR)/sample2.cpp $(SRCFILES) $(HEADERS)
	mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $^ $(LAFLAGS) -o $@
# Clean
clean:
	rm -rf $(BUILDDIR)

# Install Target
INSTALL_DIR ?= $(HOME)/gsminres_install
install: libgsminres.a
	mkdir -p $(INSTALL_DIR)/include
	cp $(INCDIR)/*.hpp $(INSTALL_DIR)/include/
	mkdir -p $(INSTALL_DIR)/lib
	mv libgsminres.a $(INSTALL_DIR)/lib/
libgsminres.a: $(SRCFILES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $(SRCFILES) -I$(INCDIR)
	ar rcs libgsminres.a *.o
	rm -f *.o
