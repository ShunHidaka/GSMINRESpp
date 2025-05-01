#ifndef GSMINRES_SOLVER_HPP
#define GSMINRES_SOLVER_HPP

#include <complex>
#include <vector>
#include <array>

namespace gsminres {
  class Solver {
  public:
    // Constructor
    Solver(std::size_t matrix_size, std::size_t shift_size);
    // Deconstructor
    ~Solver() = default;

    // Initialize Generalized shifted MINRES solver
    void initialize(std::vector<std::vector<std::complex<double>>>& x, // solution vectors
                    const std::vector<std::complex<double>>& b,        // Right-hand side vector
                    std::vector<std::complex<double>>& w,              // B^{-1}b
                    const std::vector<std::complex<double>>& sigma,    // Shift values
                    const double threshold);                           // threshold
    // Preprocess of Generalized Laczos process
    void glanczos_pre(std::vector<std::complex<double>>& u);
    // Postprocess of Generalized Lanczos process
    void glanczos_pst(std::vector<std::complex<double>>& w,
                      std::vector<std::complex<double>>& u);
    // Update approximate solutions
    bool update(std::vector<std::vector<std::complex<double>>>& x);
    // Finalize Generalized shifted MINRES solver
    void finalize(std::vector<std::size_t>& conv_itr, std::vector<double>& conv_res);
    // Get residual norms
    void get_residual(std::vector<double>& res) const;

  private:
    std::size_t iter_;                        // Number of iterations
    std::size_t matrix_size_;                 // Matrix size
    std::size_t shift_size_;                  // Number of shift
    double r0_norm_;                          // Norm of the initial residual norm
    std::vector<std::complex<double>> sigma_; // Shift values
    // Generalized Lanczos process variables
    double alpha_;
    double beta_prev_, beta_curr_;
    std::vector<std::complex<double>> w_prev_, w_curr_, w_next_; // Basis vectors
    std::vector<std::complex<double>> u_prev_, u_curr_, u_next_; // Auxiliary vectors
    // Variables for update solutions
    std::vector<std::complex<double>> T_prev2_, T_prev_, T_curr_, T_next_; // BLASを使用するためにstd::vectorに
    std::vector<std::array<double, 3>> Gc_;                      // Givens rotation matrixs element "c"
    std::vector<std::array<std::complex<double>, 3>> Gs_;        // Givens rotation matrixs element "s"
    std::vector<std::vector<std::complex<double>>> p_prev2_, p_prev_, p_curr_; // Auxiliary vectors
    std::vector<std::complex<double>> f_;                         //Auxiliary variables
    std::vector<double> h_;                                       // Residual norms in Algorithm
    // Convergence-related variables
    unsigned int conv_num_;            // Number of converged systems
    std::vector<std::size_t> is_conv_; // Converged flag for each system
    double threshold_;                 // Convergence threshold
  };

}  // namespace gsminres

#endif // GSMINRES_SOLVER_HPP
