/**
 * \file gsminres_c_api.h
 * \brief C API for the GSMINRES++ solver.
 * \author Shuntaro Hidaka
 *
 * \details This file provides a C interface for the GSMINRES++ library,
 *          which solves multiple shifted linear systems of the form
 *          \f[
 *            (A + \sigma^{(m)} B)x^{(m)} = b 
 *          \f]
 *          using a generalized MINRES method.
 *
 *          The API is designed for interoperability with C and Fortran,
 *          and uses row-major layout for storing multiple solution vectors \f$ x^{(m)} \f$.
 */

#ifndef GSMINRES_C_API_H
#define GSMINRES_C_API_H

#include <stddef.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * \brief Opaque solver handle to internal C++ \t gsminres::Solver class.
   * \details This handle abstracts the internal solver object used in C++ and
   *          allows it to be manipulated via C-compatible interface.
   *          The actual structure and memory layout are hidden from the C side.
   */
  typedef void* gsminres_handle;

  /**
   * \brief Create a new GSMINRES Solver.
   * \param[in] n Matrix size.
   * \param[in] m Number of shifts.
   * \return Solver handle.
   */
  gsminres_handle gsminres_create(size_t n, size_t m);

  /**
   * \brief Destroy the solver and free memory.
   * \param handle Solver handle.
   */
  void gsminres_destroy(gsminres_handle handle);

  /**
   * \brief Initialize the solver with input data and prepare for iteration.
   * \param[in]     handle    Solver handle.
   * \param[out]    x         Approximate solutions (size = n*m, row-major).
   * \param[in]     b         Right-hand side vector (size = n).
   * \param[in,out] w         Pre-processed right-hand side \f$ B^{-1}b \f$ (size = n).
   * \param[in]     sigma     Array of shift parameters (size = m).
   * \param[in]     threshold Convergence threshold for relative residuals.
   * \param[in]     n         Matrix size.
   * \param[in]     m         Number of shifts.
   */
  void gsminres_initialize(gsminres_handle handle,
                           void            *x,
                           const void      *b,
                           void            *w,
                           const void      *sigma,
                           const double    threshold,
                           const size_t    n,
                           const size_t    m);

  /**
   * \brief Perform the pre-processing step of the generalized Lanczos process.
   * \param[in]     handle Solver handle.
   * \param[in,out] u      Vector initially containing the \f$ Aw\f$ and is update (size = n).
   * \param[in]     n      Matrix size.
   */
  void gsminres_glanczos_pre(gsminres_handle handle,
                             void            *u,
                             const size_t    n);

  /**
   * \brief Perform the post-processing step of the generalized Lanczos process
   * \param[in]     handle Solver handle.
   * \param[in,out] w      Pre-processed vector \f$ w = B^{-1} u \f$ (size = n).
   * \param[in,out] u      Vector which used in `glanczos_pre()` (size = n).
   * \param[in]     n      Matrix size.
   */
  void gsminres_glanczos_pst(gsminres_handle handle,
                             void            *w,
                             void            *u,
                             const size_t    n);

  /**
   * \brief Update the approximate solutions and check convergence.
   * \param[in]     handle Solver handle.
   * \param[in,out] x      Approximate solutions to be updated (size = n*m, row-major).
   * \param[in]     n      Matrix size.
   * \param[in]     m      Number of shifts.
   * \return 1 if all systems converged, 0 otherwise
   */
  int gsminres_update(gsminres_handle handle,
                      void            *x,
                      const size_t    n,
                      const size_t    m);

  /**
   * \brief Retrieve converged iteration and converged residual norm.
   * \details This function does not finalize or delete the solver instance.
   *          It only retrieves the number of iterations and the residual norm in Algorithm
   *          at convergence for each shifted systems.
   *          To properly clean up the solver, call `gsminres_destory` explicitly.
   *
   * \param[in]  handle   Solver handle.
   * \param[out] conv_itr Array storage the number of converged iterations (size = m).
   * \param[out] conv_res Array storage converged residual norms (size = m).
   * \param[in]  m        Number of shifts.
   */
  void gsminres_finalize(gsminres_handle handle,
                         void            *conv_itr,
                         void            *conv_res,
                         const size_t    m);

  /**
   * @brief Retrieve current residual norms in Algorithm.
   * @param[in]  handle Solver handle.
   * @param[out] res    Array storage residual norm for each system (size = m).
   * @param[in]  m      Number of shifts.
   */
  void gsminres_get_residual(gsminres_handle handle,
                             void            *res,
                             const size_t    m);

#ifdef __cplusplus
}
#endif

#endif // GSMINRES_C_API_H
