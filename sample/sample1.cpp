#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include "gsminres_solver.hpp"
#include "gsminres_util.hpp"
#include "gsminres_blas.hpp"


int main() {
  std::size_t N, M;
  //std::string Aname = "ELSES_MATRIX_BNZ30_A.mtx", Bname = "ELSES_MATRIX_BNZ30_B.mtx"; // 非正定値に近い
  //std::string Aname = "ELSES_MATRIX_DIAB18h_A.mtx", Bname = "ELSES_MATRIX_DIAB18h_B.mtx";
  //std::string Aname = "ELSES_MATRIX_VCNT900_A.mtx", Bname = "ELSES_MATRIX_VCNT900_B.mtx";
  std::string Aname = "ELSES_MATRIX_PPE3594_20160426_A.mtx", Bname = "ELSES_MATRIX_PPE3594_20160426_B.mtx";
  const std::vector<std::complex<double>>     A = gsminres::util::load_matrix_from_mm(Aname, N);
  //const std::vector<std::complex<double>>     I = gsminres::util::generate_identity(N);
  const std::vector<std::complex<double>>     B = gsminres::util::load_matrix_from_mm(Bname, N);
  const std::vector<std::complex<double>>     b = gsminres::util::generate_ones(N);
  //const std::vector<std::complex<double>> sigma = gsminres::util::load_vector("shift.txt");
  std::vector<std::complex<double>> sigma(10);
  for(std::size_t i=0; i<10; i++) {
    std::complex<double> I(0.0, 1.0);
    std::complex<double> tmp = 2 * M_PI * I * (i+0.5) / 10.0;
    sigma[i] = 0.01 * std::exp(tmp);
    std::cout << sigma[i] << std::endl;
  }
  M = sigma.size();

  std::vector<std::vector<std::complex<double>>> x(M, std::vector<std::complex<double>>(N, {0.0, 0.0}));
  std::vector<std::complex<double>> w(N, {0.0, 0.0}), u(N, {0.0, 0.0});
  std::vector<double> res(M);

  //for(std::size_t i=0; i<I.size(); i++) B[i] = B[i] + I[i];
  std::vector<std::complex<double>> Bcholesky = B;
  gsminres::lapack::zpptrf(N, Bcholesky);

  gsminres::Solver solver(N, M);
  gsminres::lapack::zpptrs(N, Bcholesky, w, b);
  solver.initialize(b, w, sigma, 1e-13);
  for(std::size_t j=1; j<10000; j++) {
    gsminres::blas::zhpmv({1.0, 0.0}, A, w, {0.0, 0.0}, u);
    solver.glanczos_pre(u);
    gsminres::lapack::zpptrs(N, Bcholesky, w, u);
    solver.glanczos_pst(w, u);
    if(solver.update(x)) {
      std::cout << "converged in " << j << std::endl;
      break;
    }
    solver.get_residual(res);

    if(j % 10 == 1) {
      std::cout << j;
      for(std::size_t j=0; j<M/2; j++){
	std::vector<std::complex<double>> tmp(N, {0.0, 0.0});
	double tmp_nrm = 0.0;
	gsminres::blas::zhpmv({1.0, 0.0},  A, x[j], {0.0, 0.0}, tmp);
	gsminres::blas::zhpmv(sigma[j],    B, x[j], {1.0, 0.0}, tmp);
	gsminres::blas::zaxpy({-1.0, 0.0}, b, tmp);
	tmp_nrm = gsminres::blas::dznrm2(tmp);
	std::cout << " (" << res[j] << ", " << tmp_nrm << ")";
      }
      std::cout << std::endl;
    }

  }
  solver.finalize();
  solver.get_residual(res);

  for(std::size_t j=0; j<M; j++){
    std::vector<std::complex<double>> tmp(N, {0.0, 0.0});
    double tmp_nrm = 0.0;
    gsminres::blas::zhpmv({1.0, 0.0},  A, x[j], {0.0, 0.0}, tmp);
    gsminres::blas::zhpmv(sigma[j],    B, x[j], {1.0, 0.0}, tmp);
    gsminres::blas::zaxpy({-1.0, 0.0}, b, tmp);
    tmp_nrm = gsminres::blas::dznrm2(tmp);
    std::cout << j << " " << res[j] << " " << tmp_nrm << std::endl;
  }

}

/*
  gsminres::Solver solver();
  線形方程式の求解;
  solver.initialize();
  for(int j=0; j<MAX_ITR; j++) {
    行列ベクトル積;
    solver.glanczos_pre();
    線形方程式の求解;
    solver.glanczos_pst();
    if(solver.update() == true)
      break;
  }
  solver.get_residual();
  solver.finalize();
*/
