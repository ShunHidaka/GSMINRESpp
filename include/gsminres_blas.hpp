/**
 * \file gsminres_blas.hpp
 * \brief BLAS wrappers for the GSMINRES++ solver.
 * \author Shuntaro Hidaka
 *
 * \details This header defines lightweight C++ wrapper functions for
 *          selected BLAS Level-1 and Level-2 routines such as
 *          `zaxpy`, `dznrm2`, `zdotc`, and `zhpmv`, which are used internally in GSMINRES++.
 *          The interfaces are designed for safety and usability using `std::vector`,
 *          and provide explicit control over starting offsets and memory strides
 *          for advanced vector operations.
 */

#ifndef GSMINRES_BLAS_HPP
#define GSMINRES_BLAS_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cstddef>
#include <cstdlib>

extern "C" {
  void dscal_(const int *n, const double *a, double *x, const int *incx);
  void dcopy_(const int *n, const double *x, const int *incx, double *y, const int *incy);
  void zdscal_(const int *n, const double *a, std::complex<double> *x, const int *incx);
  void zscal_(const int *n, const std::complex<double> *a, std::complex<double> *x, const int *incx);
  void zcopy_(const int *n,const std::complex<double> *x, const int *incx, std::complex<double> *y, const int *incy);
  void zaxpy_(const int *n, const std::complex<double> *alpha, const std::complex<double> *x, const int *incx, std::complex<double> *y, const int *incy);
  std::complex<double> zdotc_(const int *n, const std::complex<double> *x, const int *incx, const std::complex<double> *y, const int *incy);
  double dznrm2_(const int *n, const std::complex<double> *x, const int *incx);
  void zhpmv_(char *uplo, int *n, std::complex<double> *alpha, std::complex<double> *A, std::complex<double> *x, int *incx, std::complex<double> *beta, std::complex<double> *y, int *incy);
  void zrotg_(std::complex<double> *a, std::complex<double> *b, double *c, std::complex<double> *s);
  void zrot_(const int *n, std::complex<double> *x, const int *incx, std::complex<double> *y, const int *incy, const double *c, const std::complex<double> *s);
}

/**
 * \namespace gsminres::blas
 * \brief This namespace provides C++ wrappers for BLAS routines used in GSMINRES++.
 * \details This namespace provides simple C++ wrappers over the BLAS routines.
 *          These routines are used in GSMINRES++.
 */
namespace gsminres {
  namespace blas {

    /**
     * \brief Scale a real vector \f$ x \f$ by a real scalar \f$ a \f$.
     * \param[in]     n        Number of elements to scale.
     * \param[in]     a        Real scalar multiplier.
     * \param[in,out] x        Real vector to scale.
     * \param[in]     x_offset Starting index within the x vector (zero-base offset).
     * \param[in]     incx     Step size beteen elements in the x vector (stride).
     */
    inline void dscal(std::size_t n, double a,
                      std::vector<double>& x, std::size_t x_offset=0, std::size_t incx=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx);
      dscal_(&nn, &a, x.data()+x_offset, &ix);
    }

    /**
     * \brief Copy real vector \f$ x \f$ into real vector \f$ y \f$.
     * \param[in]  n        Number of elements to copy.
     * \param[in]  x        Real source vector.
     * \param[in]  x_offset Starting index within the x vector.
     * \param[out] y        Real destination vector.
     * \param[in]  y_offset Starting index within the y vector.
     * \param[in]  incx     Step size beteen elements in the x vector.
     * \param[in]  incy     Step size beteen elements in the y vector.
     */
    inline void dcopy(std::size_t n,
                      const std::vector<double>& x, std::size_t x_offset,
                      std::vector<double>&       y, std::size_t y_offset,
                      std::size_t incx=1, std::size_t incy=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx), iy = static_cast<int>(incy);
      dcopy_(&nn, x.data()+x_offset, &ix, y.data()+y_offset, &iy);
    }

    /**
     * \brief Scale a complex vector \f$ x \f$ by a real scalar \f$ a \f$.
     * \param[in]     n        Number of elements to scale.
     * \param[in]     a        Real scalar multiplier.
     * \param[in,out] x        Complex vector to scale.
     * \param[in]     x_offset Starting index within the x vector.
     * \param[in]     incx     Step size beteen elements in the x vector.
     */
    inline void zdscal(std::size_t n, double a,
                       std::vector<std::complex<double>>& x, std::size_t x_offset=0, std::size_t incx=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx);
      zdscal_(&nn, &a, x.data()+x_offset, &ix);
    }

    /**
     * \brief Scale a complex vector \f$ x \f$ by a complex scalar \f$ a \f$.
     * \param[in]     n        Number of elements to scale.
     * \param[in]     a        Complex scalar multiplier.
     * \param[in,out] x        Complex vector to scale.
     * \param[in]     x_offset Starting index within the x vector.
     * \param[in]     incx     Step size beteen elements in the x vector.
     */
    inline void zscal(std::size_t n, std::complex<double> a,
                      std::vector<std::complex<double>>& x, std::size_t x_offset=0, std::size_t incx=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx);
      zscal_(&nn, &a, x.data()+x_offset, &ix);
    }

    /**
     * \brief Copy complex vector \f$ x \f$ into complex vector \f$ y \f$.
     * \param[in]  n        Number of elements to copy.
     * \param[in]  x        Complex source vector.
     * \param[in]  x_offset Starting index within the x vector.
     * \param[out] y        Complex destination vector.
     * \param[in]  y_offset Starting index within the y vector.
     * \param[in]  incx     Step size beteen elements in the x vector.
     * \param[in]  incy     Step size beteen elements in the y vector.
     */
    inline void zcopy(std::size_t n,
                      const std::vector<std::complex<double>>& x, std::size_t x_offset,
                      std::vector<std::complex<double>>&       y, std::size_t y_offset,
                      std::size_t incx=1, std::size_t incy=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx), iy = static_cast<int>(incy);
      zcopy_(&nn, x.data()+x_offset, &ix, y.data()+y_offset, &iy);
    }

    /**
     * \brief Perform \f$ y = \alpha x + y \f$ for complex vector
     * \param[in]  n        Number of elements to perform.
     * \param[in]  alpha    Scalar multiplier.
     * \param[in]  x        Input vector.
     * \param[in]  x_offset Starting index within the x vector.
     * \param[out] y        Output vector (accumulated).
     * \param[in]  y_offset Starting index within the y vector.
     * \param[in]  incx     Step size beteen elements in the x vector.
     * \param[in]  incy     Step size beteen elements in the y vector.
     */
    inline void zaxpy(std::size_t n, std::complex<double> alpha,
                      const std::vector<std::complex<double>>& x, std::size_t x_offset,
                      std::vector<std::complex<double>>&       y, std::size_t y_offset,
                      std::size_t incx=1, std::size_t incy=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx), iy = static_cast<int>(incy);
      zaxpy_(&nn, &alpha, x.data()+x_offset, &ix, y.data()+y_offset, &iy);
    }

    /**
     * \brief Compute dot product of complex vectors: \f$ \sum_i \overline{x[i]} * y[i]) \f$.
     * \param[in] n        Number of elements to compute.
     * \param[in] x        First input vector.
     * \param[in] x_offset Starting index within the x vector.
     * \param[in] y        Second input vector.
     * \param[in] y_offset Starting index within the y vector.
     * \param[in]  incx    Step size beteen elements in the x vector.
     * \param[in]  incy    Step size beteen elements in the y vector.
     * \return Complex scalar result.
     */
    inline std::complex<double> zdotc(std::size_t n,
                                      const std::vector<std::complex<double>>& x, std::size_t x_offset,
                                      const std::vector<std::complex<double>>& y, std::size_t y_offset,
                                      std::size_t incx=1, std::size_t incy=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx), iy = static_cast<int>(incy);
      return zdotc_(&nn, x.data()+x_offset, &ix, y.data()+y_offset, &iy);
    }

    /**
     * \brief Compute the Euclidean norm (2-norm) of a complex vector: \f$ \|x\| \f$.
     * \param[in] n        Number of elements to compute.
     * \param[in] x        Input vector.
     * \param[in] x_offset Starting index within the x vector.
     * \param[in] incx     Step size beteen elements in the x vector.
     * \return 2-norm value (double).
     */
    inline double dznrm2(std::size_t n,
                         const std::vector<std::complex<double>>& x, std::size_t x_offset=0,
                         std::size_t incx=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx);
      return dznrm2_(&nn, x.data()+x_offset, &ix);
    }

    /**
     * \brief Hermitian packed 'U' matrix-vector multiplication: \f$ y = \alpha A x + \beta y \f$.
     * \param[in]     alpha    Scalar multiplier for A*x.
     * \param[in]     A        Packed Hermitian matrix (upper triagnle stored).
     * \param[in]     x        Input vector.
     * \param[in]     x_offset TMP
     * \param[in]     incx     TMP
     * \param[in]     beta     Scalar multiplier for y.
     * \param[in,out] y        Output vector.
     * \param[in]     y_offset TMP
     * \param[in]     incy     TMP
     */
    inline void zhpmv(std::complex<double> alpha,
                      const std::vector<std::complex<double>>& A,
                      const std::vector<std::complex<double>>& x,
                      std::complex<double> beta,
                      std::vector<std::complex<double>>&       y) {
      char uplo = 'U';
      int n = static_cast<int>(x.size()), incx = 1, incy = 1;
      zhpmv_(&uplo, &n, &alpha, const_cast<std::complex<double>*>(A.data()), const_cast<std::complex<double>*>(x.data()), &incx, &beta, y.data(), &incy);
    }

    /**
     * \brief Compute Givens rotation parameters.
     * \param[in,out] a First component, overwritten.
     * \param[in,out] b Second component, overwritten.
     * \param[out]    c Cosine value (real).
     * \param[out]    s Sine value (complex).
     */
    inline void zrotg(std::complex<double>& a, std::complex<double>& b,
                      double& c, std::complex<double>& s) {
      zrotg_(&a, &b, &c, &s);
    }

    /**
     * \brief Apply Givens rotation to vector pair (x, y).
     * \param[in]     n        Number of elements to apply.
     * \param[in,out] x        First vector, overwritten c*x+s*y.
     * \param[in]     x_offset Starting index within the y vector.
     * \param[in,out] y        Second vector, overwritten -conj(s)*x+c*y.
     * \param[in]     y_offset Starting index within the y vector.
     * \param[in]     c        Cosine value.
     * \param[in]     s        Sine value (complex).
     * \param[in]     incx     Step size beteen elements in the x vector.
     * \param[in]     incy     Step size beteen elements in the y vector.
     */
    inline void zrot(std::size_t n,
                     std::vector<std::complex<double>>& x, std::size_t x_offset,
                     std::vector<std::complex<double>>& y, std::size_t y_offset,
                     double c, std::complex<double> s,
                     std::size_t incx=1, std::size_t incy=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx), iy = static_cast<int>(incy);
      zrot_(&nn, x.data()+x_offset, &ix, y.data()+y_offset, &iy, &c, &s);
    }
  }  //namespace blas
}  // namespace gsminres

#endif // GSMINRES_BLAS_HPP
