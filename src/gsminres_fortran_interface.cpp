// gsminre_fortran_interface.cpp
#include "gsminres_solver.hpp"
#include <complex>
#include <vector>
#include <cstdlib>
#include <cstring>

extern "C" {
  using SolverHandle = void*;

  //
  SolverHandle gsminres_create(int n, int m) {
    return new gsminres::Solver(static_cast<std::size_t>(n), static_cast<std::size_t>(m));
  }
  // 
  void gsminres_destroy(SolverHandle handle) {
    delete static_cast<gsminres::Solver*>(handle);
  }

  // Fortran array to C++ std::vector
  static std::vector<std::complex<double>> to_stdvec(const std::complex<double>* ptr, std::size_t size) {
    return std::vector<std::complex<double>>(ptr, ptr+size);
  }
  // C++ std::vector to Fortran array
  static void from_stdvec(const std::vector<std::complex<double>>& vec, std::complex<double>* ptr) {
    std::memcpy(ptr, vec.data(), vec.size() * sizeof(std::complex<double>));
  }

  void gsminre_initialize(SolverHandle handle,
                          std::complex<double>* x,
                          const std::complex<double>* b,
                          std::complex<double>* w,
                          const std::complex<double>* sigma,
                          int n, int m, double threshold) {
    gsminres::Solver* solver = static_cast<gsminres::Solver*>(handle);
    std::vector<std::complex<double>> vb = to_stdvec(b, n);
    std::vector<std::complex<double>> vw = to_stdvec(w, n);
    std::vector<std::complex<double>> vsigma = to_stdvec(sigma, m);
    solver->initialize(vb, vw, vsigma, threshold);
    from_stdvec(vw, w);
  }

  void gsminres_update(SolverHandle handle,
                       std::complex<double>* x,  // [m][n] in Fortran
                       int n, int m,
                       int* converged) {
    auto solver = static_cast<gsminres::Solver*>(handle);
    std::vector<std::vector<std::complex<double>>> xvec(m, std::vector<std::complex<double>>(n));
    for (int j = 0; j < m; ++j)
      std::memcpy(xvec[j].data(), x + j*n, sizeof(std::complex<double>) * n);
    bool done = solver->update(xvec);
    for (int j = 0; j < m; ++j)
      std::memcpy(x + j*n, xvec[j].data(), sizeof(std::complex<double>) * n);
    *converged = done ? 1 : 0;
  }

  void gsminres_get_residual(SolverHandle handle,
                             double* res, int m) {
    auto solver = static_cast<gsminres::Solver*>(handle);
    std::vector<double> r(m);
    solver->get_residual(r);
    std::memcpy(res, r.data(), sizeof(double)*m);
  }
}
