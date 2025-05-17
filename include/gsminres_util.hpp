/**
 * \file gsminres_util.hpp
 * \brief Utility functions and data structures for GSMINRES++
 * \details This header provides matrix/vector loaders, generators and
 *          sparse matrix data structure (CSR format), as well as utility kernels
 *          such as sparse matrix-vector multiplication and a simple CG solver.
 */

#ifndef GSMINRES_UTIL_HPP
#define GSMINRES_UTIL_HPP

#include <iostream>
#include <string>
#include <complex>
#include <vector>

namespace gsminres {
  namespace util {

    /**
     * \brief Load a packed 'U' Hermitian or symmetric matrix from a Matrix Market file.
     * \param[in]  filename File path of the Matrix Market file.
     * \param[out] size     Size of the square matrix.
     * \return Packed upper-triangular matrix in std::vector<std::complex<double>>.
     * \note Exits the program on failure.
     */
    std::vector<std::complex<double>> load_matrix_from_mm(const std::string& filename,
                                                          std::size_t&       size);

    /**
     * \brief Load a complex vector from a plain text file.
     * \param[in] filename File path of the vector file.
     * \return std::vector<std::complex<double>>.
     * \note Exits the program on failure.
     */
    std::vector<std::complex<double>> load_vector(const std::string& filename);

    /**
     * \brief Generate a vector filled with ones.
     * \param[in] size Number of elements.
     * \return st::vector<std::complex<double>> with all entries set to 1.0+0.0i.
     */
    std::vector<std::complex<double>> generate_ones(const std::size_t size);

    /**
     * \brief Generate a packed 'U' identity matrix.
     * \param[in] size Matrix size
     * \return Packed upper-triangular identity matrix.
     */
    std::vector<std::complex<double>> generate_identity(const std::size_t size);

    /**
     * \struct CSRMat
     * \brief Struct representing a sparse matrix in CSR format
     */
    struct CSRMat {
      std::size_t matrix_size; ///< Number of rows/columns (N times N square matrix)
      std::vector<std::size_t> row_pointer; ///< Row pointer array (size = N+1)
      std::vector<std::size_t> col_indices; ///< Column index array (size = nnz)
      std::vector<std::complex<double>> values; /// Non-zero values (size = nnz)
      /**
       * \brief Construct an empty CSR matrix with given size and number of non-zero elements.
       * \param[in] ROWPSIZE Size of the row pointer array (usually N+1)
       * \param[in] DATASIZE Number of non-zero elements (nnz)
       */
      CSRMat(std::size_t ROWPSIZE, std::size_t DATASIZE) :
        matrix_size(ROWPSIZE-1),
        row_pointer(ROWPSIZE, 0),
        col_indices(DATASIZE, 0),
        values(DATASIZE, std::complex<double>{0.0, 0.0}) {}
    };

    // Read MM format to CSR format
    //CSRMat load_csr_from_mm(const std::string& filename);
    /**
     * \brief Load a sparse matrix form a custom CSR-format file.
     * \param[in] filename File path of the CSR matrix
     * \return CSR matrix.
     * \note Exits the program on failure.
     */
    CSRMat load_csr_from_csr(const std::string& filename);

    /**
     * \brief Sparse matrix-vector multiplication: y = A * x
     * \param[in]  A Matrix in CSR format.
     * \param[in]  x Input vector.
     * \param[out] y Output vector (A * x).
     */
    void spmv(const CSRMat& A, const std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& y);

    /**
     * \brief Conjugate Gradient method for solving Ax = b.
     * \param[in]     A        Coefficient matrix (CSR format).
     * \param[in,out] x        Solution vector (initial guess input, solution output).
     * \param[in]     b        Right-hand side vector.
     * \param[in]     tol      Tolerance for relative residual (default = 1e-12).
     * \param[in]     max_iter Maximum number of iterations (default = 10000).
     * \return true if converged, false otherwise.
     */
    bool cg(const CSRMat&A, std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& b, const double tol, const std::size_t max_iter);

  }  // namespace util
}  //namespace gsminres

#endif // GSMINRES_UTIL_HPP
