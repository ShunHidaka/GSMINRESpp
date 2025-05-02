#include <iostream>
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

  std::vector<std::vector<std::complex<double>>> x(M, std::vector<std::complex<double>>(N, {0.0, 0.0}));
  std::vector<std::complex<double>> w(N, {0.0, 0.0}), u(N, {0.0, 0.0});
  std::vector<std::size_t> itr(M);
  std::vector<double> res(M);

  gsminres::Solver solver(N, M);
  if (!gsminres::util::cg(B, w, b, 1e-13, 10000)) {
    std::cerr << "Failed" << std::endl;
    std::exit(1);
  }
  solver.initialize(x, b, w, sigma, 1e-13);
  for(std::size_t j=1; j<10000; j++) {
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
      for(std::size_t j=0; j<M/2; j++){
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

  for(std::size_t j=0; j<M; j++){
    std::vector<std::complex<double>> tmp1(N, {0.0, 0.0}), tmp2(N, {0.0, 0.0});
    double tmp_nrm = 0.0;
    gsminres::util::spmv(A, x[j], tmp1);
    gsminres::util::spmv(B, x[j], tmp2);
    gsminres::blas::zaxpy(sigma[j], tmp2, tmp1);
    gsminres::blas::zaxpy({-1.0, 0.0}, b, tmp1);
    tmp_nrm = gsminres::blas::dznrm2(tmp1);
    std::cout << j << " " << itr[j] << " " << res[j] << " " << tmp_nrm << std::endl;
  }
}
