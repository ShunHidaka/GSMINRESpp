/**
 * \file sample2.cpp
 * \brief C++ example of using GSMINRES++ with CSR format input and built-in SpMV+CG.
 * \example sample2.cpp
 * \author Shuntaro Hidaka
 *
 * \details This example solves a set of generalized shifted linear systems of the form:
 *          \f[
 *            (A + \sigma^{(m)} B)x^{(m)} = b, \quad (m=1,\dots,M)
 *          \f]
 *          using the GSMINRES++ solver.
 *
 *          Matrices A and B are provided in a custom CSR format (`.csr`) and
 *          are read using the utilities in \ref gsminres_util.hpp "gsminres_util.cpp".
 *          Sparse matrix-vector multiplication and inner linear solves
 *          are performed using built-in routines (`SpMV` and `CG`).
 *
 *          The CSR files are generated from Matrix Market (.mtx) input files
 *          using \ref converter.py "Python script", which converts
 *          the Matrix Market format matrix into a custom CSR format.
 *
 *          A key feature of GSMINRES++ is that the user is free to implement
 *          matrix-vector multiplications and linear solves externally.
 *          This example shows built-in approach, but any other representation
 *          or computation method can be used as long as the required computational steps
 *          are performed in accordance with the expected algorithmic flow.
 *
 * \par Usage:
 * \code
 *  $ ./sample2 ../data/A.mtx ../data/B.mtx
 * \endcode
 */

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include "gsminres_solver.hpp"
#include "gsminres_util.hpp"
#include "gsminres_blas.hpp"


int main(int argc, char* argv[]) {
  std::size_t N, M;
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <CSR_file(A)> <CSR_file(B)>" << std::endl;
    return 1;
  }
  std::string Aname = argv[1], Bname = argv[2];
  const gsminres::util::CSRMat A = gsminres::util::load_csr_from_csr(Aname);
  const gsminres::util::CSRMat B = gsminres::util::load_csr_from_csr(Bname);
  N = A.matrix_size;
  const std::vector<std::complex<double>>     b = gsminres::util::generate_ones(N);
  std::vector<std::complex<double>> sigma(10);
  for(std::size_t i=0; i<10; i++) {
    std::complex<double> I(0.0, 1.0);
    std::complex<double> tmp = 2 * M_PI * I * (i+0.5) / 10.0;
    sigma[i] = 0.1 * std::exp(tmp);
  }
  M = sigma.size();

  std::vector<std::complex<double>> x(M*N, {0.0, 0.0});
  std::vector<std::complex<double>> w(N, {0.0, 0.0}), u(N, {0.0, 0.0});
  std::vector<std::size_t> itr(M);
  std::vector<double> res(M);

  gsminres::Solver solver(N, M);
  if (!gsminres::util::cg(B, w, b, 1e-13, 10000)) {
    std::cerr << "Failed" << std::endl;
    std::exit(1);
  }
  solver.initialize(x, b, w, sigma, 1e-13);
  for(std::size_t j=1; j<10000; ++j) {
    gsminres::util::spmv(A, w, u);
    solver.glanczos_pre(u);
    if (!gsminres::util::cg(B, w, u, 1e-13, 10000)) {
      std::cerr << "Failed" << std::endl;
      std::exit(1);
    }
    solver.glanczos_pst(w, u);
    if(solver.update(x)) {
      std::cout << "converged in " << j << std::endl;
      break;
    }
    solver.get_residual(res);
    /*
    if(j % 10 == 1) {
      std::cout << j;
      for(std::size_t j=0; j<M/2; ++j){
	std::vector<std::complex<double>> tmp1(N, {0.0, 0.0}), tmp2(N, {0.0, 0.0});
	double tmp_nrm = 0.0;
	gsminres::util::spmv(A, x[j], tmp1);
	gsminres::util::spmv(B, x[j], tmp2);
	gsminres::blas::zaxpy(sigma[j], tmp2, tmp1);
	gsminres::blas::zaxpy({-1.0, 0.0}, b, tmp1);
	tmp_nrm = gsminres::blas::dznrm2(tmp1);
	std::cout << " (" << res[j] << ", " << tmp_nrm << ")";
      }
      std::cout << std::endl;
    }
    */
  }
  solver.finalize(itr, res);

  for(std::size_t j=0; j<M; ++j){
    std::vector<std::complex<double>> ans(x.begin()+j*N, x.begin()+(j+1)*N);
    std::vector<std::complex<double>> tmp1(N, {0.0, 0.0}), tmp2(N, {0.0, 0.0});
    double tmp_nrm = 0.0;
    gsminres::util::spmv(A, ans, tmp1);
    gsminres::util::spmv(B, ans, tmp2);
    gsminres::blas::zaxpy(N, sigma[j], tmp2, 0, tmp1, 0);
    gsminres::blas::zaxpy(N, {-1.0, 0.0}, b, 0, tmp1, 0);
    tmp_nrm = gsminres::blas::dznrm2(N, tmp1);
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
