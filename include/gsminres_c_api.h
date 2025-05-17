/**
 * \file gsminres_c_api.h
 * \brief C API for the GSMINRES++ solver.
 * \detail
 */

#ifndef GSMINRES_C_API_H
#define GSMINRES_C_API_H

#include <stddef.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

  // Opaque solver handle to internal C++ gsminres::Solver class
  typedef void* gsminres_handle;

  /**
   * @brief Create a new GSMINRES Solver
   * @param[in] n Matrix size
   * @param[in] m Number of shifts
   * @return Solver handle
   */
  gsminres_handle gsminres_create(size_t n, size_t m);

  /**
   * @brief Destroy the solver and free memory
   * @param handle Solver handle
   */
  void gsminres_destroy(gsminres_handle handle);

  /**
   * @brief Initalize the solver
   * @param[in]     handle    Solver handle
   * @param[out]    x         Approximate solutions (length=n*m, row-major)
   * @param[in]     b         Right-hand side vector (length=n)
   * @param[in,out] w         B^{-1}b (length=n)
   * @param[in]     sigma     Shift parameters (length=m)
   * @param[in]     threshold Relative residual threshold for convergence
   * @param[in]     n         Matrix size
   * @param[in]     m         Number of shifts
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
   * @brief Apply pre-processing step of the Generalized Lanczos process
   * @param[in]     handle Solver handle
   * @param[in,out] u      A*w (length=n)
   * @param[in]     n      Matrix size
   */
  void gsminres_glanczos_pre(gsminres_handle handle,
                             void            *u,
                             const size_t    n);

  /**
   * @brief Apply post-processing step of the Generalized Lanczos process
   * @param[in]     handle Solver handle
   * @param[in,out] w      (length=n)
   * @param[in,out] u      B^{-1}Aw (length=n)
   * @param[in]     n      Matrix size
   */
  void gsminres_glanczos_pst(gsminres_handle handle,
                             void            *w,
                             void            *u,
                             const size_t    n);

  /**
   * @brief Update the approximate solutions
   * @param[in]     handle Solver handle
   * @param[in,out] x      Approximate solutions (length=n*m, row-major)
   * @param[in]     n      Matrix size
   * @param[in]     m      Number of shifts
   * @return 1 if all systems converged, 0 otherwise
   */
  int gsminres_update(gsminres_handle handle,
                      void            *x,
                      const size_t    n,
                      const size_t    m);

  /**
   * @brief Get converged iteration and converged residual norm
   * @param[in]  handle   Solver handle
   * @param[out] conv_itr Array storage the number of iterations until convergence (length=m)
   * @param[out] conv_res Array storage converged residual norm for each system (lenght=m)
   * @param[in]  m        Number of shifts
   */
  void gsminres_finalize(gsminres_handle handle,
                         void            *conv_itr,
                         void            *conv_res,
                         const size_t    m);

  /**
   * @brief Get current residual norms
   * @param[in]  handle Solver handle
   * @param[out] res    Array storage residual norm for each system (length=m)
   * @param[in]  m      Number of shifts
   */
  void gsminres_get_residual(gsminres_handle handle,
                             void            *res,
                             const size_t    m);

#ifdef __cplusplus
}
#endif

#endif // GSMINRES_C_API_H
