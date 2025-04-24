#ifndef GSMINRES_UTIL_HPP
#define GSMINRES_UTIL_HPP

#include <iostream>
#include <string>
#include <complex>
#include <vector>

namespace gsminres {
  namespace util {

    // Read MM format to "packed" "U" BLAS/LAPACK format
    std::vector<std::complex<double>> load_matrix_from_mm(const std::string& filename, std::size_t& size);
    // Read vector
    std::vector<std::complex<double>> load_vector(const std::string& filename);
    // Generate all-ones vector
    std::vector<std::complex<double>> generate_ones(const std::size_t size);
    // Generate identity matirx ("packed" "U" BLAS/LAPACK format)
    std::vector<std::complex<double>> generate_identity(const std::size_t size);

    // Structure representing a matrix in Compressed Sparse Row format
    struct CSRMat {
      std::size_t matrix_size;
      std::vector<std::size_t> row_pointer;
      std::vector<std::size_t> col_indices;
      std::vector<std::complex<double>> values;
      CSRMat(std::size_t ROWPSIZE, std::size_t DATASIZE) : matrix_size(ROWPSIZE-1),
                                                           row_pointer(ROWPSIZE, 0),
                                                           col_indices(DATASIZE, 0),
                                                           values(DATASIZE, std::complex<double>{0.0, 0.0}) {}
    };
    // Read MM format to CSR format
    //CSRMat load_csr_from_mm(const std::string& filename);
    // Read CSR format to CSR format
    CSRMat load_csr_from_csr(const std::string& filename);
    // CSR format Matrix-Vector multiplication
    void spmv(const CSRMat& A, const std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& y);
    // Conjugate Gradient method
    bool cg(const CSRMat&A, std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& b, const double tol, const std::size_t max_iter);

  }  // namespace util
}  //namespace gsminres

#endif // GSMINRES_UTIL_HPP
