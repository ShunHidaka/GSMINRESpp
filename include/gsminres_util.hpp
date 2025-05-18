/**
 * \file gsminres_util.hpp
 * \brief Utility functions and data structures for GSMINRES++
 * \author Shuntaro Hidaka
 *
 * \details This header provides common utility components related to
 *          matrix and vector operations, including functions to load and generate test data,
 *          a custom sparse matrix data structure in CSR format,
 *          and basic computational kernels such as sparse matrix-vector multiplication
 *          and a simple Conjugate Gradient (CG) solver.
 *
 *          These functions are not used internally by the GSMINRES++ solver itself,
 *          but are provided for convenience in sample programs, testing, and exploratory use.
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
     * \param[in]  filename Path to the Matrix Market file.
     * \param[out] size     Size of the resulting square matrix.
     * \return `std::vector<std::complex<double>>` containing packed matrix (upper-triangular).
     * \note Exits the program on failure.
     */
    std::vector<std::complex<double>> load_matrix_from_mm(const std::string& filename,
                                                          std::size_t&       size);

    /**
     * \brief Load a complex-value vector from a plain text file.
     * \param[in] filename Path to the vector file.
     * \return `std::vector<std::complex<double>>`.
     * \note Exits the program on failure.
     */
    std::vector<std::complex<double>> load_vector(const std::string& filename);

    /**
     * \brief Generate a vector filled with ones.
     * \param[in] size Number of elements.
     * \return `std::vector<std::complex<double>>` with all entries set to `1.0+0.0i`.
     */
    std::vector<std::complex<double>> generate_ones(const std::size_t size);

    /**
     * \brief Generate a packed identity matrix in upper-triangular format.
     * \param[in] size Dimension of the square matrix
     * \return `std::vector<std::complex<double>>` storing identity matrix in packed format (upper triangle).
     */
    std::vector<std::complex<double>> generate_identity(const std::size_t size);

    /**
     * \struct CSRMat
     * \brief Struct representing a sparse matrix in Compressed Sparse Row (CSR) format.
     */
    struct CSRMat {
      std::size_t                       matrix_size; ///< Dimension of the square matrix (N).
      std::vector<std::size_t>          row_pointer; ///< Row pointer array (size = N+1).
      std::vector<std::size_t>          col_indices; ///< Column index array (size = nnz).
      std::vector<std::complex<double>> values;      ///< Non-zero values (size = nnz).
      /**
       * \brief Construct for CSRMat
       * \param[in] ROWPSIZE Size of the row pointer array (N+1).
       * \param[in] DATASIZE Number of non-zero elements (nnz).
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
     * \brief Load a sparse matrix from a file in custom CSR format.
     * \param[in] filename Path to the file containing CSR-formatted matrix.
     * \return CSR matrix object.
     * \note Exits the program on failure.
     */
    CSRMat load_csr_from_csr(const std::string& filename);

    /**
     * \brief Perform sparse matrix-vector multiplication: \f$ y = A x \f$.
     * \param[in]  A Matrix in CSR format.
     * \param[in]  x Input vector.
     * \param[out] y Output vector where result is stored.
     */
    void spmv(const CSRMat&                            A,
              const std::vector<std::complex<double>>& x,
              std::vector<std::complex<double>>&       y);

    /**
     * \brief Solve \f$ Ax=b \f$ using the Conjugate Gradient method.
     * \param[in]  A        Coefficient matrix (CSR format).
     * \param[out] x        Solution vector.
     * \param[in]  b        Right-hand side vector.
     * \param[in]  tol      Relative residual tolerance (default = 1e-12).
     * \param[in]  max_iter Maximum number of iterations (default = 10000).
     * \return true if converged, false otherwise.
     */
    bool cg(const CSRMat&                            A,
            std::vector<std::complex<double>>&       x,
            const std::vector<std::complex<double>>& b,
            const double tol, const std::size_t max_iter);

  }  // namespace util
}  //namespace gsminres

#endif // GSMINRES_UTIL_HPP
