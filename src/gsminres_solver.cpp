#include "gsminres_solver.hpp"
#include "gsminres_blas.hpp"
#include <iostream>
#include <cmath>

namespace gsminres {

  Solver::Solver(std::size_t matrix_size, std::size_t shift_size)
    : iter_(1),
      matrix_size_(matrix_size),
      shift_size_(shift_size),
      r0_norm_(0.0),
      sigma_(shift_size, {0.0, 0.0}),
      alpha_(0.0),
      beta_prev_(0.0),
      beta_curr_(0.0),
      w_prev_(matrix_size, {0.0, 0.0}),
      w_curr_(matrix_size, {0.0, 0.0}),
      w_next_(matrix_size, {0.0, 0.0}),
      u_prev_(matrix_size, {0.0, 0.0}),
      u_curr_(matrix_size, {0.0, 0.0}),
      u_next_(matrix_size, {0.0, 0.0}),
      T_prev2_(1, {0.0, 0.0}),
      T_prev_( 1, {0.0, 0.0}),
      T_curr_( 1, {0.0, 0.0}),
      T_next_( 1, {0.0, 0.0}),
      Gc_(shift_size, std::array<double, 3>{0.0, 0.0, 0.0}),
      Gs_(shift_size, std::array<std::complex<double>, 3>{{{0.0,0.0}, {0.0,0.0}, {0.0,0.0}}}),
      p_prev2_(shift_size, std::vector<std::complex<double>>(matrix_size, {0.0, 0.0})),
      p_prev_( shift_size, std::vector<std::complex<double>>(matrix_size, {0.0, 0.0})),
      p_curr_( shift_size, std::vector<std::complex<double>>(matrix_size, {0.0, 0.0})),
      f_(shift_size, {1.0, 0.0}),
      h_(shift_size, 1.0),
      conv_num_(0),
      is_conv_(shift_size, 0),
      threshold_(1e-12) {
  }

  void Solver::initialize(std::vector<std::vector<std::complex<double>>>& x,
                          const std::vector<std::complex<double>>& b,
                          std::vector<std::complex<double>>& w,
                          const std::vector<std::complex<double>>& sigma,
                          const double threshold) {
    for (std::size_t i=0; i<shift_size_; i++) {
      blas::zdscal(0.0, x[i]);
    }
    r0_norm_ = std::sqrt((blas::zdotc(b, w)).real());
    blas::zcopy(w, w_curr_);
    blas::zcopy(b, u_curr_);
    blas::zdscal(1.0/r0_norm_, w_curr_);
    blas::zdscal(1.0/r0_norm_, u_curr_);
    blas::zcopy(w_curr_, w);
    blas::dscal(r0_norm_, h_);
    blas::zcopy(sigma, sigma_);
    threshold_ = threshold;
  }

  void Solver::glanczos_pre(std::vector<std::complex<double>>& u) {
    alpha_ = (blas::zdotc(w_curr_, u)).real();
    blas::zaxpy(-alpha_,     u_curr_, u);
    blas::zaxpy(-beta_prev_, u_prev_, u);
  }

  void Solver::glanczos_pst(std::vector<std::complex<double>>& w,
                            std::vector<std::complex<double>>& u) {
    beta_curr_ = std::sqrt((blas::zdotc(u, w)).real());
    blas::zdscal(1.0/beta_curr_, w);
    blas::zdscal(1.0/beta_curr_, u);
    blas::zcopy(w, w_next_);
    blas::zcopy(u, u_next_);
  }

  bool Solver::update(std::vector<std::vector<std::complex<double>>>& x) {
    for (std::size_t m=0; m<shift_size_; m++) {
      if (is_conv_[m] != 0) {
        continue;
      }
      T_prev2_[0] = 0.0;
      T_prev_[0]  = beta_prev_;
      T_curr_[0]  = alpha_ + sigma_[m];
      T_next_[0]  = beta_curr_;
      if (iter_ >= 3) {
        blas::zrot(T_prev2_, T_prev_, Gc_[m][0], Gs_[m][0]);
      }
      if (iter_ >= 2) {
        blas::zrot(T_prev_,  T_curr_, Gc_[m][1], Gs_[m][1]);
      }
      blas::zrotg(T_curr_[0], T_next_[0], Gc_[m][2], Gs_[m][2]);
      blas::zcopy(p_prev_[m], p_prev2_[m]);
      blas::zcopy(p_curr_[m], p_prev_[m]);
      blas::zcopy(w_curr_,    p_curr_[m]);
      blas::zaxpy(-T_prev2_[0], p_prev2_[m], p_curr_[m]);
      blas::zaxpy(-T_prev_[0],  p_prev_[m],  p_curr_[m]);
      blas::zscal(1.0/T_curr_[0], p_curr_[m]);
      blas::zaxpy(r0_norm_*Gc_[m][2]*f_[m], p_curr_[m], x[m]);
      f_[m] = -std::conj(Gs_[m][2]) * f_[m];
      h_[m] = std::abs(-std::conj(Gs_[m][2])) * h_[m];
      if (h_[m]/r0_norm_ < threshold_) {
        conv_num_++;
        is_conv_[m] = iter_;
        continue;
      }
      Gc_[m][0] = Gc_[m][1]; Gc_[m][1] = Gc_[m][2];
      Gs_[m][0] = Gs_[m][1]; Gs_[m][1] = Gs_[m][2];
    }
    beta_prev_ = beta_curr_;
    blas::zcopy(w_curr_, w_prev_); blas::zcopy(w_next_, w_curr_);
    blas::zcopy(u_curr_, u_prev_); blas::zcopy(u_next_, u_curr_);
    iter_++;
    if (conv_num_ >= shift_size_) {
      return true;
    }
    return false;
  }

  void Solver::finalize(std::vector<std::size_t>& conv_itr,
                        std::vector<double>&      conv_res) {
    // 当初はメモリの解放などを行う予定だったが
    // (動的な確保をおこなっていないため)不要なので収束までの反復回数と残差のノルムを返す関数とする
    conv_itr = is_conv_;
    conv_res = h_;
  }

  void Solver::get_residual(std::vector<double>& res) const {
    blas::dcopy(h_, res);
  }
}
