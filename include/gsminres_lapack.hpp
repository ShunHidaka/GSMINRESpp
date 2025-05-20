/**
 * \file gsminres_lapack.hpp
 * \brief LAPACK wrapper functions for the GSMINRES++ solver.
 * \author Shuntaro Hidaka
 *
 * \details This header provides minimal C++ wrappers for selected LAPACK routines
 *          used in sample programs of the GSMINRES++ solver.
 *          In particular, it supports operations on packed Hermitian matrices
 *          stored in upper-triangular ('U') format.
 *
 *          The routines include Cholesky factorization (zpptrf), linear solve (zpptrs),
 *          and robust Givens rotation computation (zlartg) as a workaround for known
 *          OpenBLAS issues in zrot.
 */

#ifndef GSMINRES_LAPACK_HPP
#define GSMINRES_LAPACK_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cstddef>
#include <cstdlib>

extern "C" {
  void zpptrf_(char *uplo, int *n, std::complex<double> *ap, int *info);
  void zpptrs_(char *uplo, int *n, int *nrhs, std::complex<double> *ap, std::complex<double> *b, int *ldb, int *info);
  void zlartg_(std::complex<double> *f, std::complex<double> *g, double *c, std::complex<double> *s, std::complex<double> *r);
}

/**
 * \namespace gsminres::lapack
 * \brief This namespace provides C++ wrappers for LAPACK routines used in GSMINRES++ sample programs.
 * \details This namespace provides simple LAPACK routines necessary for
 *          Cholesky factorization and solving linear equations, used in sample programs.
 *
 *          In addition, it provides a wrapper for `zlartg` as a numerically robust
 *          alternative to the BLAS routine `zrotg`, which is known to produce incorrect
 *          results in certain environments such as OpenBLAS versions prior to 0.3.27.
 */
namespace gsminres {
  namespace lapack {

    /**
     * \brief Perform Cholesky factorization of a Hermitian positive-definite packed matrix.
     * \param[in]     n  Dimension of the matrix (n x n).
     * \param[in,out] ap Packed Hermitian matrix stored in upper triangle (column-major).
     *                   On output, contains the Cholesky factor.
     * \note Exits the program on failure.
     */
    inline void zpptrf(int n, std::vector<std::complex<double>>& ap) {
      char uplo = 'U';
      int info = 0;
      zpptrf_(&uplo, &n, ap.data(), &info);
      if (info != 0) {
        std::cerr << "zpptrf: [ERROR] failed INFO = " << std::to_string(info) << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    /**
     * \brief Solve Ax = b using the Cholesky factorization of a Hermitian positive-definite packed matrix.
     * \param[in]  n  Dimension of the matrix (n x n).
     * \param[in]  ap Cholesky factor of A (as computed by zpptrf).
     * \param[out] x  Solution vector on output.
     * \param[in]  b  Input right-hand side vector (copied internally to x).
     * \note Exits the program on failure.
     */
    inline void zpptrs(int n, const std::vector<std::complex<double>>& ap, std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& b){
      char uplo = 'U';
      int nrhs = 1, ldb = n, info = 0;
      blas::zcopy(n, b, 0, x, 0);
      zpptrs_(&uplo, &n, &nrhs, const_cast<std::complex<double>*>(ap.data()), x.data(), &ldb, &info);
      if (info != 0) {
        std::cerr << "zpptrf: [ERROR] failed INFO = " << std::to_string(info) << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    /**
     * \brief Compute Givens rotation parameters for complex values using LAPACK's `zlartg` routine.
     * \details This function computes the parameters of a Givens rotation matrix
     *          that eliminates the second entry of a 2-vector.
     *          It is a LAPACK-based alternative to zrotg, used in GSMINRES++
     *          to work around known bugs in older versions of OpenBLAS
     *          where zrot does not behave correctly.
     *
     * \param[in,out] f On input, the first component. On output, replaced with the resulting r.
     * \param[in]     g On input, the second component.
     * \param[out]    c The cosine of the rotation.
     * \param[out]    s The sine of the rotation (complex).
 */
    inline void zlartg(std::complex<double>& f, std::complex<double>& g,
                       double& c, std::complex<double>& s){
      std::complex<double> r;
      zlartg_(&f, &g, &c, &s, &r);
      f = r;
    }
      
  }  // namespace lapack
}  // namespace gsminres

#endif // GSMINRES_LAPACK_HPP
