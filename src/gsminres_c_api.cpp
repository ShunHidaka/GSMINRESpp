#include "gsminres_solver.hpp"
#include "gsminres_c_api.h"
#include "gsminres_c_api_util.hpp"
#include <complex>
#include <vector>
#include <cstring>

inline gsminres::Solver* as_solver(gsminres_handle handle) {
  return static_cast<gsminres::Solver*>(handle);
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
                           double _Complex       *x,
                           const double _Complex *b,
                           double _Complex       *w,
                           const double _Complex *sigma,
                           const double          threshold,
                           const size_t          n,
                           const size_t          m) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<std::complex<double>> xvec = to_cpp_vector(x, n*m);
    std::vector<std::complex<double>> bvec = to_cpp_vector(b, n);
    std::vector<std::complex<double>> wvec = to_cpp_vector(w, n);
    std::vector<std::complex<double>> svec = to_cpp_vector(sigma, m);
    solver->initialize(xvec, bvec, wvec, svec, threshold);
    from_cpp_vector(xvec, x);
    from_cpp_vector(wvec, w);
  }

  void gsminres_glanczos_pre(gsminres_handle handle,
                             double _Complex *u,
                             const size_t    n) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<std::complex<double>> uvec = to_cpp_vector(u, n);
    solver->glanczos_pre(uvec);
    from_cpp_vector(uvec, u);
  }

  void gsminres_glanczos_pst(gsminres_handle handle,
                             double _Complex *w,
                             double _Complex *u,
                             const size_t    n) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<std::complex<double>> wvec = to_cpp_vector(w, n);
    std::vector<std::complex<double>> uvec = to_cpp_vector(u, n);
    solver->glanczos_pst(wvec, uvec);
    from_cpp_vector(wvec, w);
    from_cpp_vector(uvec, u);
  }

  int gsminres_update(gsminres_handle handle,
                      double _Complex *x,
                      const size_t    n,
                      const size_t    m) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<std::complex<double>> xvec = to_cpp_vector(x, n*m);
    bool converged = solver->update(xvec);
    from_cpp_vector(xvec, x);
    return converged ? 1 : 0;
  }

  void gsminres_finalize(gsminres_handle handle,
                         int          *conv_itr,
                         double       *conv_res,
                         const size_t m) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<std::size_t> itr(m);
    std::vector<double> res(m);
    solver->finalize(itr, res);
    for (size_t i = 0; i < m; ++i) {
      conv_itr[i] = static_cast<int>(itr[i]);
      conv_res[i] = res[i];
    }
  }

  void gsminres_get_residual(gsminres_handle handle,
                             double       *res,
                             const size_t m) {
    gsminres::Solver* solver = as_solver(handle);
    std::vector<double> rvec(m);
    solver->get_residual(rvec);
    for (size_t i=0; i < m; ++i) {
      res[i] = rvec[i];
    }
  }

}
