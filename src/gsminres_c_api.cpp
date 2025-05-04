#include "gsminres_solver.hpp"
#include "gsminres_c_api.h"
#include "gsminres_c_api_util.hpp"
#include <complex>
#include <vector>
#include <cstring>
#include <iostream>

inline gsminres::Solver* as_solver(gsminres_handle handle) {
  return static_cast<gsminres::Solver*>(handle);
}
inline double _Complex* as_cmplx(void* ptr) {
  return reinterpret_cast<double _Complex*>(ptr);
}

using namespace gsminres_c_api_util;

extern "C" {

  gsminres_handle gsminres_create(size_t n, size_t m) {
    return new gsminres::Solver(n, m);
  }

  void gsminres_destroy(gsminres_handle handle) {
    delete as_solver(handle);
  }

  void gsminres_initialize(gsminres_handle handle,
                           void*           x,
                           const void*     b,
                           void*           w,
                           const void*     sigma,
                           const double    threshold,
                           const size_t    n,
                           const size_t    m) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<std::complex<double>> xvec = to_cpp_vector(as_cmplx(x), n*m);
    std::vector<std::complex<double>> bvec = to_cpp_vector(as_cmplx(const_cast<void *>(b)), n);
    std::vector<std::complex<double>> wvec = to_cpp_vector(as_cmplx(w), n);
    std::vector<std::complex<double>> svec = to_cpp_vector(as_cmplx(const_cast<void *>(sigma)), m);
    solver->initialize(xvec, bvec, wvec, svec, threshold);
    from_cpp_vector(xvec, as_cmplx(x));
    from_cpp_vector(wvec, as_cmplx(w));
  }

  void gsminres_glanczos_pre(gsminres_handle handle,
                             void*           u,
                             const size_t    n) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<std::complex<double>> uvec = to_cpp_vector(as_cmplx(u), n);
    //std::cout << "Before: " << uvec[0] << ", ";
    solver->glanczos_pre(uvec);
    //std::cout << "After: " << uvec[0] << std::endl;
    from_cpp_vector(uvec, as_cmplx(u));
  }

  void gsminres_glanczos_pst(gsminres_handle handle,
                             void*           w,
                             void*           u,
                             const size_t    n) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<std::complex<double>> wvec = to_cpp_vector(as_cmplx(w), n);
    std::vector<std::complex<double>> uvec = to_cpp_vector(as_cmplx(u), n);
    solver->glanczos_pst(wvec, uvec);
    from_cpp_vector(wvec, as_cmplx(w));
    from_cpp_vector(uvec, as_cmplx(u));
  }

  int gsminres_update(gsminres_handle handle,
                      void*           x,
                      const size_t    n,
                      const size_t    m) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<std::complex<double>> xvec = to_cpp_vector(as_cmplx(x), n*m);
    bool converged = solver->update(xvec);
    from_cpp_vector(xvec, as_cmplx(x));
    return converged ? 1 : 0;
  }

  void gsminres_finalize(gsminres_handle handle,
                         void*           conv_itr,
                         void*           conv_res,
                         const size_t m) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<std::size_t> itr(m);
    std::vector<double> res(m);
    solver->finalize(itr, res);
    for (size_t i = 0; i < m; ++i) {
      reinterpret_cast<int*>(conv_itr)[i]    = static_cast<int>(itr[i]);
      reinterpret_cast<double*>(conv_res)[i] = res[i];
    }
  }

  void gsminres_get_residual(gsminres_handle handle,
                             void*           res,
                             const size_t    m) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<double> rvec(m);
    solver->get_residual(rvec);
    for (size_t i=0; i < m; ++i) {
      reinterpret_cast<double*>(res)[i] = rvec[i];
    }
  }

}
