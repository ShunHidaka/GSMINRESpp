#ifndef GSMINRES_LAPACK_HPP
#define GSMINRES_LAPACK_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cstddef>
#include <cstdlib>

extern "C" {
  // LAPACK
  void zpptrf_(char *uplo, int *n, std::complex<double> *ap, int *info);
  void zpptrs_(char *uplo, int *n, int *nrhs, std::complex<double> *ap, std::complex<double> *b, int *ldb, int *info);
}

namespace gsminres {
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
      blas::zcopy(n, b, 0, x, 0);
      zpptrs_(&uplo, &n, &nrhs, const_cast<std::complex<double>*>(ap.data()), x.data(), &ldb, &info);
      if (info != 0) {
        std::cerr << "zpptrf: [ERROR] failed INFO = " << std::to_string(info) << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
      
  }  // namespace lapack
}  // namespace gsminres

#endif // GSMINRES_LAPACK_HPP
