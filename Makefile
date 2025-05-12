# =====================================
# Derectory configurations
# =====================================
SRCDIR = src
INCDIR = include
OBJDIR = bin/obj
BINDIR = bin

$(shell mkdir -p $(BINDIR) $(OBJDIR))

# =====================================
# Compiler and Flags
# =====================================
CXX      = g++
CC       = gcc
FC       = gfortran
CXXFLAGS = -std=c++17 -O3 -Wextra -fPIC
CFLAGS   = -std=c99 -O3 -Wall -fPIC
FFLAGS   = -O3 -Wall -fPIC -J$(OBJDIR)
LDFLAGS  = 
LALIBS   = -lblas -llapack -lm # Linear Algeblic Libraries

USE_OPENMP               = 1
ENABLE_C_API             = 1
ENABLE_FORTRAN_INTERFACE = 1

ifeq ($(USE_OPENMP), 1)
	CXXFLAGS += -fopenmp
	CFLAGS   += -fopenmp
	FFLAGS   += -fopenmp
	LALIBS   += -fopenmp
endif

# =====================================
# Source files and Object files
# =====================================
SRC_CPP = src/gsminres_solver.cpp src/gsminres_util.cpp
SRC_C   = src/gsminres_c_api.cpp
SRC_F   = src/gsminres_fortran_interface.f90

OBJ_CPP = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC_CPP))
OBJ_C   = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC_C))
OBJ_F   = $(patsubst $(SRCDIR)/%.f90, $(OBJDIR)/%.o, $(SRC_F))

# =====================================
# Library and executable
# =====================================
LIB_SHARED = $(BINDIR)/libgsminres.so
LIB_STATIC = $(BINDIR)/libgsminres.a

SAMPLES = $(BINDIR)/sample1 $(BINDIR)/sample2
ifeq ($(ENABLE_C_API), 1)
	SAMPLES += $(BINDIR)/sample2_c
endif
ifeq ($(ENABLE_FORTRAN_INTERFACE), 1)
	SAMPLES += $(BINDIR)/sample1_f
endif

# =====================================
# Build Targets
# =====================================
all: $(LIB_SHARED) $(LIB_STATIC) $(SAMPLES)

$(LIB_SHARED): $(OBJ_CPP) $(if $(filter 1,$(ENABLE_C_API)),$(OBJ_C)) $(if $(filter 1,$(ENABLE_FORTRAN_INTERFACE)),$(OBJ_F))
	$(CXX) -shared -o $@ $^ $(LALIBS)
$(LIB_STATIC): $(OBJ_CPP) $(if $(filter 1,$(ENABLE_C_API)),$(OBJ_C)) $(if $(filter 1,$(ENABLE_FORTRAN_INTERFACE)),$(OBJ_F))
	ar rcs $@ $^

$(BINDIR)/sample1: sample/sample1.cpp $(LIB_SHARED)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -o $@ $< -L$(BINDIR) -lgsminres $(LALIBS)
$(BINDIR)/sample2: sample/sample2.cpp $(LIB_SHARED)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -o $@ $< -L$(BINDIR) -lgsminres $(LALIBS)
$(BINDIR)/sample1_f: sample/sample1_f.f90 $(LIB_SHARED)
ifeq ($(ENABLE_FORTRAN_INTERFACE), 1)
	$(FC) $(FFLAGS) -I$(INCDIR) -o $@ $< -L$(BINDIR) -lgsminres $(LALIBS)
else
	@echo "Fortran interface is disabled. Enable it by setting ENABLE_FORTRAN_INTERFACE = 1"
endif
$(BINDIR)/sample2_c: sample/sample2_c.c $(LIB_SHARED)
ifeq ($(ENABLE_C_API), 1)
	$(CC) $(CFLAGS) -I$(INCDIR) -o $@ $< -L$(BINDIR) -lgsminres $(LALIBS)
else
	@echo "C API is disabled. Enable it by setting ENABLE_C_API = 1"
endif

# =====================================
# Pattern rule
# =====================================
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -I$(INCDIR) -c $< -o $@

# =====================================
# Install setting
# =====================================
PREFIX     ?= $(HOME)/gsminres_install
INCDIR_DST  = $(PREFIX)/include
LIBDIR_DST  = $(PREFIX)/lib

install: $(LIB_SHARED) $(LIB_STATIC)
	@echo "Installing to $(PREFIX)"
	@mkdir -p $(INCDIR_DST) $(LIBDIR_DST)
	@cp -r $(INCDIR)/* $(INCDIR_DST)
	@cp $(LIB_SHARED) $(LIB_STATIC) $(LIBDIR_DST)
	@echo "Installed libgsminres.{so,a} and headers to $(PREFIX)"

# =====================================
# Clean up
# =====================================
clean:
	rm -rf bin

.PHONY: all clean install
