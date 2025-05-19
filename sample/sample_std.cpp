/**
 * \file sample_std.cpp
 * \brief Example of using GSMINRES++ with identity matrix as B (standard shifted MINRES).
 * \example sample_std.cpp
 * \author Shuntaro Hidaka
 *
 * \details This example solves a set of (standard) shifted linear systems of the form:
 *          \f[
 *            (A + \sigma^{(m)} I)x^{(m)} = b, \quad (m=1,\dots,M)
 *          \f]
 *          where \f$ I \f$ is identity matrix.
 *          Unlike the generalized shifted linear systems, this version assumes that
 *          \f$ B = I \f$, so the solver simplifies to the (standard) shifted MINRES method.
 *
 *          Matrices A is given in Matrix Market format (`.mtx`) and
 *          are loaded as packed `U` Hermitian matrices.
 *          The solver uses BLAS routines internally for matrix-vector multiplication (`zhpmv`).
 *
 *          A key feature of GSMINRES++ is that the user is free to implement
 *          matrix-vector multiplications.
 *          This example shows one such approach using BLAS, but any other representation
 *          or computation method can be used as long as the required computational steps
 *          are performed in accordance with the expected algorithmic flow.
 *
 * \par Usage:
 * \code
 *  $ ./sample_std ../data/A.mtx
 * \endcode
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <vector>
#include "gsminres_solver.hpp"
#include "gsminres_util.hpp"
#include "gsminres_blas.hpp"
#include "gsminres_lapack.hpp"


int main(int argc, char* argv[]) {
  std::size_t N, M;
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <MTX_file(A)>" << std::endl;
    return 1;
  }
  std::string Aname = argv[1];
  const std::vector<std::complex<double>> A = gsminres::util::load_matrix_from_mm(Aname, N);
  const std::vector<std::complex<double>> b = gsminres::util::generate_ones(N);
  std::vector<std::complex<double>> sigma(10);
  for (std::size_t m=0; m<10; ++m) {
    std::complex<double> I(0.0, 1.0);
    sigma[m] = 0.01 * std::exp(2 * M_PI * I * (m+0.5) / 10.0);
  }
  M = sigma.size();

  std::vector<std::complex<double>> x(M*N, {0.0, 0.0});
  std::vector<std::complex<double>> v(N, {0.0, 0.0}), Av(N, {0.0, 0.0});
  std::vector<std::size_t> itr(M);
  std::vector<double> res(M);

  gsminres::Solver solver(N, M);
  gsminres::blas::zcopy(N, b, 0, v, 0);
  solver.initialize(x, b, v, sigma, 1e-13);
  for (std::size_t j=0; j<10000; ++j) {
    gsminres::blas::zhpmv({1.0, 0.0}, A, v, {0.0, 0.0}, Av);
    solver.glanczos_pre(Av);
  gsminres::blas::zcopy(N, Av, 0, v, 0);
    solver.glanczos_pst(v, Av);
    if (solver.update(x)) {
      break;
    }
    solver.get_residual(res);
  }
  solver.finalize(itr, res);

  for (std::size_t j=0; j<M; ++j){
    std::vector<std::complex<double>> ans(x.begin()+j*N, x.begin()+(j+1)*N);
    std::vector<std::complex<double>> tmp(N, {0.0, 0.0});
    double tmp_nrm = 0.0;
    gsminres::blas::zhpmv({1.0, 0.0},  A, ans, {0.0, 0.0}, tmp);
    gsminres::blas::zaxpy(N,    sigma[j], ans,          0, tmp, 0);
    gsminres::blas::zaxpy(N, {-1.0, 0.0},   b,          0, tmp, 0);
    tmp_nrm = gsminres::blas::dznrm2(N, tmp);
    std::cout << std::right
              << std::setw(2) << j << " "
              << std::fixed << std::setw(10) << std::setprecision(6) << sigma[j].real() << " "
              << std::fixed << std::setw(10) << std::setprecision(6) << sigma[j].imag() << " "
              << std::setw(5) << itr[j] << " "
              << std::scientific << std::setw(12) << std::setprecision(5) << res[j] << " "
              << std::scientific << std::setw(12) << std::setprecision(5) << tmp_nrm
              << std::endl;
  }
}
