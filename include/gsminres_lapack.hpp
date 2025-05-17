/**
 * \file gsminres_lapack.hpp
 * \brief LAPACK wrapper functions for the GSMINRES++ solver.
 * \details Provede minimal C++ wrapperf for selected LAPACK routines used sample programs.
 *          The wrappers support packed Hermitian 'U' matrices.
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

namespace gsminres {
  namespace lapack {

    /**
     * \brief Perform Cholesky factorization of a Hermitian positive-definite packed matrix.
     * \param[in]     n  Dimension of the matrix (n times n).
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
     * \param[in]     n  Dimension of the matrix (n x n).
     * \param[in]     ap Cholesky factor of A (as computed by zpptrf).
     * \param[in,out] x  Right-hand side vector on input, solution vector on output.
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
     * \brief Compute Givens rotation parameters for complex values using LAPACK's zlartg routine.
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
