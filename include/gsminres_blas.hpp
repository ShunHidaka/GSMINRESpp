#ifndef GSMINRES_BLAS_HPP
#define GSMINRES_BLAS_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cstddef>
#include <cstdlib>

extern "C" {
  // BLAS
  void dscal_(const int *n, const double *a, double *x, const int *incx);
  void dcopy_(const int *n, const double *x, const int *incx, double *y, const int *incy);
  void zdscal_(const int *n, const double *a, std::complex<double> *x, const int *incx);
  void zscal_(const int *n, const std::complex<double> *a, std::complex<double> *x, const int *incx);
  void zcopy_(const int *n,
              const std::complex<double> *x, const int *incx,
              std::complex<double> *y, const int *incy);
  void zaxpy_(const int *n, const std::complex<double> *alpha,
              const std::complex<double> *x, const int *incx,
              std::complex<double> *y, const int *incy);
  std::complex<double> zdotc_(const int *n,
                              const std::complex<double> *x, const int *incx,
                              const std::complex<double> *y, const int *incy);
  double dznrm2_(const int *n, const std::complex<double> *x, const int *incx);
  void zhpmv_(char *uplo, int *n, std::complex<double> *alpha, std::complex<double> *A, std::complex<double> *x, int *incx, std::complex<double> *beta, std::complex<double> *y, int *incy);
  void zrotg_(std::complex<double> *a, std::complex<double> *b, double *c, std::complex<double> *s);
  void zrot_(const int *n,
             std::complex<double> *x, const int *incx,
             std::complex<double> *y, const int *incy,
             const double *c, const std::complex<double> *s);
}

namespace gsminres {
  namespace blas {
    //
    // Real double precision BLAS
    //
    inline void dscal(std::size_t n, double a,
                      std::vector<double>& x, std::size_t x_offset=0, std::size_t incx=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx);
      dscal_(&nn, &a, x.data()+x_offset, &ix);
    }

    inline void dcopy(std::size_t n,
                      const std::vector<double>& x, std::size_t x_offset,
                      std::vector<double>&       y, std::size_t y_offset,
                      std::size_t incx=1, std::size_t incy=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx), iy = static_cast<int>(incy);
      dcopy_(&nn, x.data()+x_offset, &ix, y.data()+y_offset, &iy);
    }

    //
    // Complex double precision BLAS
    //
    inline void zdscal(std::size_t n, double a,
                       std::vector<std::complex<double>>& x, std::size_t x_offset=0, std::size_t incx=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx);
      zdscal_(&nn, &a, x.data()+x_offset, &ix);
    }

    inline void zscal(std::size_t n, std::complex<double> a,
                      std::vector<std::complex<double>>& x, std::size_t x_offset=0, std::size_t incx=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx);
      zscal_(&nn, &a, x.data()+x_offset, &ix);
    }

    inline void zcopy(std::size_t n,
                      const std::vector<std::complex<double>>& x, std::size_t x_offset,
                      std::vector<std::complex<double>>&       y, std::size_t y_offset,
                      std::size_t incx=1, std::size_t incy=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx), iy = static_cast<int>(incy);
      zcopy_(&nn, x.data()+x_offset, &ix, y.data()+y_offset, &iy);
    }

    inline void zaxpy(std::size_t n, std::complex<double> alpha,
                      const std::vector<std::complex<double>>& x, std::size_t x_offset,
                      std::vector<std::complex<double>>&       y, std::size_t y_offset,
                      std::size_t incx=1, std::size_t incy=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx), iy = static_cast<int>(incy);
      zaxpy_(&nn, &alpha, x.data()+x_offset, &ix, y.data()+y_offset, &iy);
    }

    inline std::complex<double> zdotc(std::size_t n,
                                      const std::vector<std::complex<double>>& x, std::size_t x_offset,
                                      const std::vector<std::complex<double>>& y, std::size_t y_offset,
                                      std::size_t incx=1, std::size_t incy=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx), iy = static_cast<int>(incy);
      return zdotc_(&nn, x.data()+x_offset, &ix, y.data()+y_offset, &iy);
    }

    inline double dznrm2(std::size_t n,
                         const std::vector<std::complex<double>>& x, std::size_t x_offset=0,
                         std::size_t incx=1) {
      int nn = static_cast<int>(n);
      int ix = static_cast<int>(incx);
      return dznrm2_(&nn, x.data()+x_offset, &ix);
    }

    inline void zhpmv(std::complex<double> alpha,
                      const std::vector<std::complex<double>>& A,
                      const std::vector<std::complex<double>>& x,
                      std::complex<double> beta,
                      std::vector<std::complex<double>>&       y) {
      char uplo = 'U';
      int n = static_cast<int>(x.size()), incx = 1, incy = 1;
      zhpmv_(&uplo, &n, &alpha, const_cast<std::complex<double>*>(A.data()), const_cast<std::complex<double>*>(x.data()), &incx, &beta, y.data(), &incy);
    }

    //
    // BLAS Givens rotation
    //
    inline void zrotg(std::complex<double>& a, std::complex<double>& b,
                      double& c, std::complex<double>& s) {
      zrotg_(&a, &b, &c, &s);
    }

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
