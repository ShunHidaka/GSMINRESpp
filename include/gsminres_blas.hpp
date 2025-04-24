#ifndef GSMINRES_BLAS_HPP
#define GSMINRES_BLAS_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cstdlib>

extern "C" {
  // using solver
  void dscal_(int *n, double *a, double *x, int *incx);
  void zdscal_(int *n, double *a, std::complex<double> *x, int *incx);
  void zscal_(int *n, std::complex<double> *a, std::complex<double> *x, int *incx);
  void zaxpy_(int *n, std::complex<double> *alpha, std::complex<double> *x, int *incx, std::complex<double> *y, int *incy);
  void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
  void zcopy_(int *n, std::complex<double> *x, int *incx, std::complex<double> *y, int *incy);
  std::complex<double> zdotc_(int *n, std::complex<double> *x, int *incx, std::complex<double> *y, int *incy);
  double dznrm2_(int *n, std::complex<double> *x, int *incx);
  void zrotg_(std::complex<double> *a, std::complex<double> *b, double *c, std::complex<double> *s);
  void zrot_(int *n, std::complex<double> *x, int *incx, std::complex<double> *y, int *incy, double *c, std::complex<double> *s);
  // using util
  void zhpmv_(char *uplo, int *n, std::complex<double> *alpha, std::complex<double> *A, std::complex<double> *x, int *incx, std::complex<double> *beta, std::complex<double> *y, int *incy);
  void zpptrf_(char *uplo, int *n, std::complex<double> *ap, int *info);
  void zpptrs_(char *uplo, int *n, int *nrhs, std::complex<double> *ap, std::complex<double> *b, int *ldb, int *info);
}

namespace gsminres {
  namespace blas {
    inline void dscal(double a, std::vector<double>& x) {
      int n = static_cast<int>(x.size()), incx = 1;
      dscal_(&n, &a, x.data(), &incx);
    }
    inline void zdscal(double a, std::vector<std::complex<double>>& x) {
      int n = static_cast<int>(x.size()), incx = 1;
      zdscal_(&n, &a, x.data(), &incx);
    }
    inline void zscal(std::complex<double> a, std::vector<std::complex<double>>& x) {
      int n = static_cast<int>(x.size()), incx = 1;
      zscal_(&n, &a, x.data(), &incx);
    }
    inline void zaxpy(std::complex<double> alpha, const std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& y) {
      int n = static_cast<int>(x.size()), incx = 1, incy = 1;
      zaxpy_(&n, &alpha, const_cast<std::complex<double>*>(x.data()), &incx, y.data(), &incy);
    }
    inline void dcopy(const std::vector<double>& x, std::vector<double>& y) {
      int n = static_cast<int>(x.size()), incx = 1, incy = 1;
      dcopy_(&n, const_cast<double*>(x.data()), &incx, y.data(), &incy);
    }
    inline void zcopy(const std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& y) {
      int n = static_cast<int>(x.size()), incx = 1, incy = 1;
      zcopy_(&n, const_cast<std::complex<double>*>(x.data()), &incx, y.data(), &incy);
    }
    inline std::complex<double> zdotc(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& y) {
      int n = static_cast<int>(x.size()), incx = 1, incy = 1;
      return zdotc_(&n, const_cast<std::complex<double>*>(x.data()), &incx, const_cast<std::complex<double>*>(y.data()), &incy);
    }
    inline double dznrm2(const std::vector<std::complex<double>>& x) {
      int n = static_cast<int>(x.size()), incx = 1;
      return dznrm2_(&n, const_cast<std::complex<double>*>(x.data()), &incx);
    }
    inline void zrotg(std::complex<double>& a, std::complex<double>& b, double& c, std::complex<double>& s) {
      zrotg_(&a, &b, &c, &s);
    }
    inline void zrot(std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& y, double c, std::complex<double> s) {
      int n = static_cast<int>(x.size()), incx = 1, incy = 1;
      zrot_(&n, x.data(), &incx, y.data(), &incy, &c, &s);
    }
    inline void zhpmv(std::complex<double> alpha, const std::vector<std::complex<double>>& A, const std::vector<std::complex<double>>& x, std::complex<double> beta, std::vector<std::complex<double>>& y) {
      char uplo = 'U';
      int n = static_cast<int>(x.size()), incx = 1, incy = 1;
      zhpmv_(&uplo, &n, &alpha, const_cast<std::complex<double>*>(A.data()), const_cast<std::complex<double>*>(x.data()), &incx, &beta, y.data(), &incy);
    }
  }  //namespace blas

  namespace lapack {
    inline void zpptrf(int n, std::vector<std::complex<double>>& ap) {
      char uplo = 'U';
      int info = 0;
      zpptrf_(&uplo, &n, ap.data(), &info);
      if (info != 0) {
        std::cerr << "zpptrf: [ERROR] failed INFO = " << std::to_string(info) << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    inline void zpptrs(int n, const std::vector<std::complex<double>>& ap, std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& b){
      char uplo = 'U';
      int nrhs = 1, ldb = n, info = 0;
      blas::zcopy(b, x);
      zpptrs_(&uplo, &n, &nrhs, const_cast<std::complex<double>*>(ap.data()), x.data(), &ldb, &info);
      if (info != 0) {
        std::cerr << "zpptrf: [ERROR] failed INFO = " << std::to_string(info) << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
      
  }  // namespace lapack
}  // namespace gsminres

#endif // GSMINRES_BLAS_HPP
