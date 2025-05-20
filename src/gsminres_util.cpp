/**
 * \file gsminres_util.cpp
 * \brief Implementation of Utility functions and data structures for GSMINRES++.
 * \author Shuntaro Hidaka
 */

#include "gsminres_util.hpp"
#include "gsminres_blas.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <vector>
#include <cstdlib>

namespace gsminres {
  namespace util {

    std::vector<std::complex<double>> load_matrix_from_mm(const std::string& filename, std::size_t& size) {
      // Referenced: https://math.nist.gov/MatrixMarket/mmio-c.html
      // Open file
      std::ifstream inputFile(filename);
      if (!inputFile) {
        std::cerr << "load_matrix_from_mm: [ERROR] Unable to open file " << filename << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // Analyze header (assumes matrix coordinate)
      std::string line;
      bool isReal      = false;
      bool isComplex   = false;
      bool isSymmetric = false;
      bool isHermitian = false;
      if (std::getline(inputFile, line)) {
        if (line.find("%%MatrixMarket matrix coordinate") != std::string::npos) {
          if (line.find("real")      != std::string::npos) isReal      = true;
          if (line.find("complex")   != std::string::npos) isComplex   = true;
          if (line.find("symmetric") != std::string::npos) isSymmetric = true;
          if (line.find("hermitian") != std::string::npos) isHermitian = true;           
        } else {
          std::cerr << "load_matrix_from_mm: [ERROR] Inappropriate format " << filename << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      // Skip comments
      while (std::getline(inputFile, line)) {
        if (line[0] == '%') continue;
        else                break;
      }
      // Read matrix size
      std::istringstream iss(line);
      std::size_t numRows, numCols, numVals;
      if (!(iss >> numRows >> numCols >> numVals)) {
        std::cerr << "load_matrix_from_mm: [ERROR] Failed to read matrix size from " << filename << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (numRows != numCols) {
        std::cerr << "load_matrix_from_mm: [ERROR] Matrix is not square in " << filename << std::endl;
        std::exit(EXIT_FAILURE);
      }
      size = numRows;
      // Read matrix elements
      std::vector<std::complex<double>> mat(size*(size+1)/2, std::complex<double>(0.0, 0.0));
      std::size_t row, col;
      double real, imag;
      for (std::size_t i=0; i < numVals; ++i) {
        if (isReal && isSymmetric) {
          if (!(inputFile >> row >> col >> real)) {
            std::cerr << "load_matrix_from_mm: [ERROR] Invalid matrix elements in " << filename << std::endl;
            std::exit(EXIT_FAILURE);
          }
          row -= 1; col -= 1;
          if (row <= col) {
            mat[row + col*(col+1)/2] = {real, 0.0};
          } else {
            mat[col + row*(row+1)/2] = {real, 0.0};
          }
        } else if (isComplex && isHermitian) {
          if (!(inputFile >> row >> col >> real >> imag)) {
            std::cerr << "load_matrix_from_mm: [ERROR] Invalid matrix elements in " << filename << std::endl;
            std::exit(EXIT_FAILURE);
          }
          row -= 1; col -= 1;
          if (row <= col) {
            mat[row + col*(col+1)/2] = {real, imag};
          } else {
            mat[col + row*(row+1)/2] = {real, -imag};
          }
        } else {
          std::cerr << "load_matrix_from_mm: [ERROR] Invalid matrix format in " << filename << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      return mat;
    }

    std::vector<std::complex<double>> load_vector(const std::string& filename) {
      // Open file
      std::ifstream inputFile(filename);
      if (!inputFile) {
        std::cerr << "load_vector: [ERROR] Unable to open file " << filename << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // Read the Number of Elements
      std::size_t numElements;
      if (!(inputFile >> numElements)) {
        std::cerr << "load_vector: [ERROR] Failed to read the number of elements from " << filename << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::vector<std::complex<double>> vec(numElements, std::complex<double>(0.0, 0.0));
      // Read the elements
      double real=0.0, imag=0.0;
      for (std::size_t count=0; count < numElements; ++count) {
        if (!(inputFile >> real >> imag)) {
          std::cerr << "load_vector: [ERROR] Failed to read element at " << count << " from " << filename << std::endl;
          std::exit(EXIT_FAILURE);
        }
        vec[count] = std::complex<double>(real, imag);
      }
      return vec;
    }

    std::vector<std::complex<double>> generate_ones(const std::size_t size) {
      return std::vector<std::complex<double>>(size, std::complex<double>(1.0, 0.0));
    }

    std::vector<std::complex<double>> generate_identity(const std::size_t size) {
      std::vector<std::complex<double>> vec(size*(size+1)/2, std::complex<double>(0.0, 0.0));
      for (std::size_t i=0; i < size; ++i) {
        vec[i + i*(i+1)/2] = {1.0, 0.0};
      }
      return vec;
    }
    
    // CSR functions
    CSRMat load_csr_from_csr(const std::string& filename) {
      std::ifstream inputFile(filename);
      if (!inputFile) {
        std::cerr << "load_vector: [ERROR] Unable to open file " << filename << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::string line;
      std::size_t ROWPSIZE, DATASIZE, tmp;
      while (std::getline(inputFile, line)) {
        if (line[0] == '#') { continue;}
        else                { break;}
      }
      std::istringstream iss(line);
      if (!(iss >> ROWPSIZE >> DATASIZE >> tmp)) {
        std::cerr << "load_matrix_from_mm: [ERROR] Failed to read matrix size from " << filename << std::endl;
        std::exit(EXIT_FAILURE);
      }
      CSRMat mat(ROWPSIZE, DATASIZE);
      std::size_t row, col;
      double real, imag;
      for (std::size_t i=0; i<DATASIZE; ++i) {
        if (!(inputFile >> row >> col >> real >> imag)) {
          std::cerr << "load_csr_from_csr: [ERROR] Invalid matrix elements in " << filename << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (i < ROWPSIZE) { mat.row_pointer[i] = row;}
        mat.col_indices[i] = col;
        mat.values[i] = {real, imag};
      }
      return mat;
    }

    void spmv(const CSRMat& A, const std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& y) {
      #pragma omp parallel for
      for (std::size_t i=0; i < A.matrix_size; ++i) {
        y[i] = {0.0, 0.0};
        for (std::size_t j=A.row_pointer[i]; j < A.row_pointer[i+1]; ++j) {
          y[i] += A.values[j] * x[A.col_indices[j]];
        }
      }
    }

    bool cg(const CSRMat& A, std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& b, const double tol=1e-12, const std::size_t max_iter=10000) {
      bool status = false;
      std::size_t N = A.matrix_size;
      double r0nrm = blas::dznrm2(N, b);
      std::vector<std::complex<double>> r(N), p(N), Ap(N);
      std::complex<double> alpha, beta, rr, rr_old;
      blas::zdscal(N, 0.0, x);
      blas::zcopy(N, b, 0, r, 0);
      blas::zcopy(N, r, 0, p, 0);
      rr = blas::zdotc(N, r, 0, r, 0);
      for (std::size_t i=0; i < max_iter; ++i) {
        spmv(A, p, Ap);
        alpha = rr / blas::zdotc(N, p, 0, Ap, 0);
        blas::zaxpy(N, alpha,   p, 0, x, 0);
        blas::zaxpy(N, -alpha, Ap, 0, r, 0);
        if (blas::dznrm2(N, r)/r0nrm < tol) {
          status = true;
          break;
        }
        rr_old = rr;
        rr = blas::zdotc(N, r, 0, r, 0);
        beta = rr / rr_old;
        blas::zscal(N, beta, p);
        blas::zaxpy(N, {1.0, 0.0}, r, 0, p, 0);
      }
      return status;
    }

  }  // namespace util
}  // namespace gsminres
