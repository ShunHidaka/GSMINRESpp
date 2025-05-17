/**
 * \file gsminres_solver.hpp
 * \brief Header file for the GSMINRES++ solver class
 * \details This file declares the Solver class, which implements
 *          the Generalized shifted MINRES method for solving
 *          multiple generalized shifted linear systems of the form
 *          (A + sigma^{(m)}B)x^{(m)} = b.
 */

#ifndef GSMINRES_SOLVER_HPP
#define GSMINRES_SOLVER_HPP

#include <complex>
#include <vector>
#include <array>

namespace gsminres {

  /**
   * \class Solver
   * \brief Generalized shifted MINRES solver class.
   * \details This class solves a set of shifted linear systems using
   *          the MINRES method and the generalized Lanczos process.
   */
  class Solver {
  public:
    /**
     * @brief Constructor.
     * @param[in] matrix_size Matrix size
     * @param[in] shift_size  Number of shifts
     */
    Solver(std::size_t matrix_size, std::size_t shift_size);

    /**
     * @brief Deconstructor.
     * @details Maybe do not need to use
     */
    ~Solver() = default;

    /**
     * \brief Set up the solver with input data and prepare for iteration.
     * \param[out]    x         Approximate solutions (size = matrix_size * shift_size).
     * \param[in]     b         Right-hand side vector (size = matrix_size).
     * \param[in,out] w         Pre-processed right-hand side B^{-1}b (size = matrix_size).
     * \param[in]     sigma     Array of shift parameters (size = shift_size).
     * \param[in]     threshold Convergence threshold for relative residuals.
     */
    void initialize(std::vector<std::complex<double>>& x,
                    const std::vector<std::complex<double>>& b,
                    std::vector<std::complex<double>>& w,
                    const std::vector<std::complex<double>>& sigma,
                    const double threshold);

    /**
     * \brief Perform the pre-processing step of the generalized Lanczos process.
     * \param[in,out] u Vector to which matrix-vector multiplication is applied u=A*w
     */
    void glanczos_pre(std::vector<std::complex<double>>& u);

    /**
     * \brief Perform the post-processing step of the generalized Lanczos process.
     * \param[in,out] w Pre-processed vector w=B^{-1}u
     * \param[in,out] u Vector which used in glanczos_pre
     */
    void glanczos_pst(std::vector<std::complex<double>>& w,
                      std::vector<std::complex<double>>& u);

    /**
     * \bried Update the approximate solutions and check convergence.
     * \param[in,out] x Solution vectors to be updated (size = matrix_size * shift_size)
     * \return true if all systems have converged, false otherwise.
     */
    bool update(std::vector<std::complex<double>>& x);

    /**
     * \brief Finalize the solver and retrieve iteration info.
     * \param[out] conv_itr Number of iterations for each shift (size = shift_size)
     * \param[out] conv_res Final residual norms in Algorithm for each shift (size = shift_size)
     */
    void finalize(std::vector<std::size_t>& conv_itr, std::vector<double>& conv_res);

    /**
     * \brief Retrieve current residual norms in Algorithm.
     * \param[out] res Residual norms in Algorithm for each shift (shift = shift_size).
     */
    void get_residual(std::vector<double>& res) const;

  private:
    // Basic algorithm parameters
    std::size_t iter_;                        ///< Number of iterations
    std::size_t matrix_size_;                 ///< Matrix size
    std::size_t shift_size_;                  ///< Number of shift

    double r0_norm_;                          ///< Norm of the initial residual norm

    std::vector<std::complex<double>> sigma_; ///< Shift values

    // Variables for generalized Lanczos process variables
    double alpha_;                 ///< alpha coeffcient
    double beta_prev_, beta_curr_; ///< beta coefficients (previous and current)
    std::vector<std::complex<double>> w_prev_, w_curr_, w_next_; ///< Lanczos basis vectors
    std::vector<std::complex<double>> u_prev_, u_curr_, u_next_; ///< Auxiliary vectors

    // Variables for updating the sollutions
    /**
     * @brief Elements of the tridiagonal matrix by Lanczos process
     * These vectors store the column-wise elements of a tridiagonal matrix (T),
     * used during the Lanczos process. Although each vector contains only one
     * element in practice, they are wrapped in 'std::vector' to be compatible
     * with my BLAS routine wropper that require vector inputs.
     */
    std::vector<std::complex<double>> T_prev2_, T_prev_, T_curr_, T_next_;
    std::vector<std::array<double, 3>>               Gc_; ///< Givens rotation matrixs element "c"
    std::vector<std::array<std::complex<double>, 3>> Gs_; ///< Givens rotation matrixs element "s"
    std::vector<std::complex<double>> p_prev2_, p_prev_, p_curr_; ///< Auxiliary vectors (shift*matrix)
    std::vector<std::complex<double>> f_; ///< Auxiliary variables
    std::vector<double> h_;               ///< Residual norms in Algorithm
    // Convergence-related variables
    unsigned int conv_num_;            ///< Number of converged systems
    std::vector<std::size_t> is_conv_; ///< Converged flag for each system
    double threshold_;                 ///< Convergence threshold
  };

}  // namespace gsminres

#endif // GSMINRES_SOLVER_HPP
