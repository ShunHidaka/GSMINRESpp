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
      p_prev2_(shift_size*matrix_size, {0.0, 0.0}),
      p_prev_( shift_size*matrix_size, {0.0, 0.0}),
      p_curr_( shift_size*matrix_size, {0.0, 0.0}),
      f_(shift_size, {1.0, 0.0}),
      h_(shift_size, 1.0),
      conv_num_(0),
      is_conv_(shift_size, 0),
      threshold_(1e-12) {
  }

  void Solver::initialize(std::vector<std::complex<double>>& x,
                          const std::vector<std::complex<double>>& b,
                          std::vector<std::complex<double>>& w,
                          const std::vector<std::complex<double>>& sigma,
                          const double threshold) {
    blas::zdscal(shift_size_*matrix_size_, 0.0, x);
    r0_norm_ = std::sqrt((blas::zdotc(matrix_size_, b, 0, w, 0)).real());
    blas::zcopy(matrix_size_, w, 0, w_curr_, 0);
    blas::zcopy(matrix_size_, b, 0, u_curr_, 0);
    blas::zdscal(matrix_size_, 1.0/r0_norm_, w_curr_);
    blas::zdscal(matrix_size_, 1.0/r0_norm_, u_curr_);
    blas::zcopy(matrix_size_, w_curr_, 0, w, 0);
    blas::dscal(shift_size_, r0_norm_, h_);
    blas::zcopy(shift_size_, sigma, 0, sigma_, 0);
    threshold_ = threshold;
  }

  void Solver::glanczos_pre(std::vector<std::complex<double>>& u) {
    alpha_ = (blas::zdotc(matrix_size_, w_curr_, 0, u, 0)).real();
    blas::zaxpy(matrix_size_, -alpha_,     u_curr_, 0, u, 0);
    blas::zaxpy(matrix_size_, -beta_prev_, u_prev_, 0, u, 0);
  }

  void Solver::glanczos_pst(std::vector<std::complex<double>>& w,
                            std::vector<std::complex<double>>& u) {
    beta_curr_ = std::sqrt((blas::zdotc(matrix_size_, u, 0, w, 0)).real());
    blas::zdscal(matrix_size_, 1.0/beta_curr_, w);
    blas::zdscal(matrix_size_, 1.0/beta_curr_, u);
    blas::zcopy(matrix_size_, w, 0, w_next_, 0);
    blas::zcopy(matrix_size_, u, 0, u_next_, 0);
  }

  bool Solver::update(std::vector<std::complex<double>>& x) {
    for (std::size_t m=0; m<shift_size_; m++) {
      if (is_conv_[m] != 0) {
        continue;
      }
      T_prev2_[0] = 0.0;
      T_prev_[0]  = beta_prev_;
      T_curr_[0]  = alpha_ + sigma_[m];
      T_next_[0]  = beta_curr_;
      if (iter_ >= 3) {
        blas::zrot(1, T_prev2_, 0, T_prev_, 0, Gc_[m][0], Gs_[m][0]);
      }
      if (iter_ >= 2) {
        blas::zrot(1, T_prev_,  0, T_curr_, 0, Gc_[m][1], Gs_[m][1]);
      }
      blas::zrotg(T_curr_[0], T_next_[0], Gc_[m][2], Gs_[m][2]);
      std::size_t offset = m*matrix_size_;
      blas::zcopy(matrix_size_, p_prev_, offset, p_prev2_, offset);
      blas::zcopy(matrix_size_, p_curr_, offset, p_prev_,  offset);
      blas::zcopy(matrix_size_, w_curr_,    0,   p_curr_,  offset);
      blas::zaxpy(matrix_size_, -T_prev2_[0], p_prev2_, offset, p_curr_, offset);
      blas::zaxpy(matrix_size_, -T_prev_[0],  p_prev_,  offset, p_curr_, offset);
      blas::zscal(matrix_size_, 1.0/T_curr_[0], p_curr_, offset);
      blas::zaxpy(matrix_size_, r0_norm_*Gc_[m][2]*f_[m], p_curr_, offset, x, offset);
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
    blas::zcopy(matrix_size_, w_curr_, 0, w_prev_, 0);
    blas::zcopy(matrix_size_, w_next_, 0, w_curr_, 0);
    blas::zcopy(matrix_size_, u_curr_, 0, u_prev_, 0);
    blas::zcopy(matrix_size_, u_next_, 0, u_curr_, 0);
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
    blas::dcopy(shift_size_, h_, 0, res, 0);
  }
}
