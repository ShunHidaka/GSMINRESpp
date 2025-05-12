# GSMINRESpp

A high-performance solver library for generalized shifted linear systems using MINRES method.


## Overview

This library solves multiple shifted linear systems of the form:

$$ (A + \sigma^{(m)}B)x^{(m)} = b, \quad (m=1, 2, \dots, M) $$

where $A$ is Hermitian (or real symmetric) matrix and $B$ is positive definite Hermitian (or real symmetric) matrices.


## Requirement
- C++ compiler
- BLAS library
- LAPACK library (Used in the sample programs, but 無ければコンパイルが通らない)
- C compiler (Optional, for C API)
- Fortran compiler (Optional, for Fortran interface)
- `make` or `CMake`
- Python3 with `numpy` and `scipy` (used for matrix data conversion)
<!-- zrot_ が原因
/usr/bin/ld: bin/libgsminres.so: undefined reference to zrot_'
collect2: error: ld returned 1 exit status 
make: *** [Makefile:83: bin/sample2_c] エラー 1
-->


## Documents
- Manual for the library

## Directory Structure

    .  
    ├── CMakeLists.txt                         #
    ├── Makefile                               #
    ├── README.md                              #
    ├── bin                                    # Directory Makefile での
    ├── build                                  # Directory CMake での ユーザーが創る
    ├── cmake  
    │   ├── gsminresConfig.cmake.in  
    ├── data  
    │   ├── check_PD.py                        # MTX形式の行列が正定値か調べる
    │   ├── converter.py                       # MTX形式の行列をCSR形式に変換する
    ├── include  
    │   ├── gsminres_blas.hpp                  #
    │   ├── gsminres_c_api.h                   #
    │   ├── gsminres_c_api_util.hpp            #
    │   ├── gsminres_lapack.hpp                #
    │   ├── gsminres_solver.hpp                #
    │   ├── gsminres_util.hpp                  #
    ├── sample  
    │   ├── sample1.cpp                        #
    │   ├── sample1_f.f90                      #
    │   ├── sample2.cpp                        #
    │   ├── sample2_c.c                        #
    ├── src  
    │   ├── gsminres_c_api.cpp                 #
    │   ├── gsminres_fortran_interface.f90     #
    │   ├── gsminres_solver.cpp                #
    │   ├── gsminres_util.cpp                  #


## Install
もっとも簡単な方法は
``` bash
mkdir build
cd build
cmake ..
make         # サンプルプログラムがコンパイルされる
make install # デフォルトではホームディレクトリにインストール
```
or
``` bash
# Makefile の設定を適切に変更する
make
make install
```

## Sample programs


## Usage


## API Summary
See the manual.


## Known Issues
- OpenBLAS versions prior to 0.3.27 has bug in the `zrotg`.
  - See: https://github.com/OpenMathLib/OpenBLAS/issues/4909
  - **Workarounds**
    - Update OpenBLAS version 0.3.27 or later.
    - Use an alternative BLAS implementation (e.g., Netlib BLAS or Interl MKL).
    - Optionally, modify the source ~~~

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

## Licnse
MIT License
