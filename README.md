# GSMINRESpp

A high-performance solver library for generalized shifted linear systems using MINRES method.

---

## Overview

This library solves multiple shifted linear systems of the form:

$$ (A + \sigma^{(m)}B)\textbf{x}^{(m)} = \textbf{b}, \quad (m=1, 2, \dots, M) $$

where $A$ is **real symmetric** or **Hermitian** matrix and $B$ is **positive definite real symmetric** or **positive definite Hermitian** matrices.

---

## Requirement

- C++17 compiler
- BLAS library
- LAPACK library (required for sample programs, but **complication will fail without it in some BLAS/LAPACK**)
- C compiler (optional, for C API)
- Fortran compiler (optional, for Fortran interface)
- `make` or `CMake`
- Python3 with `numpy` and `scipy` (used for matrix data conversion)
<!-- zrot_ が原因
/usr/bin/ld: bin/libgsminres.so: undefined reference to zrot_'
collect2: error: ld returned 1 exit status 
make: *** [Makefile:83: bin/sample2_c] エラー 1
-->

---

## Documents

- [Manual](https://github.com/ShunHidaka/GSMINRESpp/#readme)

---

## Directory Structure

```
.  
├── CMakeLists.txt                         # CMake build script
├── Makefile                               # Make build script
├── README.md                              # This file
├── bin/                                   # Output directory (created by Make)
├── build/                                 # Build directory (created by CMake)
├── cmake/                                 #
│   ├── gsminresConfig.cmake.in            # CMake configuration file
├── data/  
│   ├── check_PD.py                        # Check if matrix is positive definite
│   ├── converter.py                       # Convert Matrix Market format to CSR
├── include/  
│   ├── gsminres_blas.hpp                  #
│   ├── gsminres_c_api.h                   #
│   ├── gsminres_c_api_util.hpp            #
│   ├── gsminres_lapack.hpp                #
│   ├── gsminres_solver.hpp                #
│   ├── gsminres_util.hpp                  #
├── sample/  
│   ├── sample1.cpp                        # C++ example (Matrix Market packed format)
│   ├── sample1_f.f90                      # Fortran example
│   ├── sample2.cpp                        # C++ example (CSR format)
│   ├── sample2_c.c                        # C example
├── src/  
│   ├── gsminres_c_api.cpp                 #
│   ├── gsminres_fortran_interface.f90     #
│   ├── gsminres_solver.cpp                #
│   ├── gsminres_util.cpp                  #
```

---

## Installation

### Using CMake (recommended)
``` bash
mkdir build
cd build
cmake ..
make           # Build sample programs and libraries
make install   # Install to $HOME/gsminres_install by default
```
### Using Makefile
``` bash
# Edit the Makefile options correctly
make           # Build sample programs and libraries in bin/
make install   # Install to $HOME/gsminres_install by default
```
See the manual for detailed instructions.

---

## How to use sample programs

This library provides several example programs to demonstrate how to use the GSMINRES++ solver in different formats and languages.

### 1. `sample1.cpp`: C++ + Matrix Market Format
C++ program using packed Hermitian matrices in Matrix Market format. LAPACK routines (`zhpmv`, `zpptrf`, `zpptrs`) are used for matrix-vector multiplication and solving linear systems in inner iterations.
``` bash
./sample1 ../data/A.mtx ../data/B.mtx
```
### 2. `sample2.cpp`: C++ + Sparse CSR Format
C++ program using CSR format defined in `gsminres_util.cpp`. Built-in functions (`gsminres::util::SpMV`, `gsminres::util::cg`) are used for matrix-vector multiplication and inner solves.
``` bash
./sample2 ../data/A.csr ../data/B.csr
```
### 3. `sample1_f.f90`: Fortan + Matrix Market Format
Fortran program using packed Hermitian matrices in Matrix Market format. LAPACK routines (`zhpmv`, `zpptrf`, `zpptrs`) are used for matrix-vector multiplication and solving linear systems in inner iterations.
``` bash
./sample1_f ../data/A.mtx ../data/B.mtx
```
### 4. `sample2_c.c`: C + Sparse CSR Format
C program using CSR format matrices defined in `sample2_c.c`. Matrix-vector multiplication and CG solves are also defined in the same source file.
``` bash
./sample2_c ../data/A.csr ../data/B.csr
```

Do not forget to convert the matrices using the Python scripts in `data/`:
```bash
python3 data/converter.py A.mtx A.csr
python3 data/converter.py B.mtx B.csr
```

---

## How to link this library
``` bash
```
詳細はマニュアルを参照

---

## API Summary
See the manual.

---

## Known Issues
- OpenBLAS versions prior to 0.3.27 has bug in the `zrotg`.
  - See: https://github.com/OpenMathLib/OpenBLAS/issues/4909
  - **Workarounds**
    - Update OpenBLAS version 0.3.27 or later.
    - Use an alternative BLAS implementation (e.g., Netlib BLAS or Interl MKL).
    - Optionally, modify the source ~~~

---

## Citation
If you use this code, please cite:
``` bibtex
@article{
  author  = {},
  title   = {},
  doi     = {},
  journal = {},
  volume  = {},
  pages   = {},
  year    = {}
}
```

---

## Licnse
MIT License
