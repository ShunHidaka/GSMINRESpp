/**
 * \file gsminres_c_api_util.hpp
 * \brief Internal utility functions for GSMINRES++ C API complex data conversion.
 * \author Shuntaro Hidaka
 *
 * \details This header provides conversion routines
 *          between C-style `double _Complex` arrays and C++ `std::complex<double>` vectors.
 *
 *          It handles both ABI-compatible and ABI-incompatible cases for
 *          `std::complex<double>` depending on the definition of `GSMINRES_ASSUME_ABI_COMPATIBLE`.
 *          These functions are used internally in the C API implementation
 *          and are not intended for external use.
 */

#ifndef GSMINRES_C_API_UTIL_HPP
#define GSMINRES_C_API_UTIL_HPP

#include <complex>
#include <vector>
#include <cstring>
#include <type_traits>

namespace gsminres_c_api_util {

#ifdef GSMINRES_ASSUME_ABI_COMPATIBLE
  /**
   * \brief Convert a C-style `double _Complex` array to a C++ `std::vector`.
   * \details If `GSMINRES_ASSUME_ABI_COMPATIBLE` is defined, this performs a `reinterpret_cast`;
   *          otherwise, it copies the data byte-wise using `std::memcpy`.
   *
   * \param[in] src Pointer to input C-style complex array (size = n).
   * \param[in] n   Number of elements.
   * \return `std::vector<std::complex<double>>` of size `n`.
   */
  inline std::vector<std::complex<double>> to_cpp_vector(const double _Complex* src, std::size_t n) {
    return std::vector<std::complex<double>>(reinterpret_cast<const std::complex<double>*>(src),
                                             reinterpret_cast<const std::complex<double>*>(src) + n);
  }
#else // If not ABI Compatibility between "std::complex<double>" and "double _Complex"
  /**
   * \brief Convert a C-style `double _Complex` array to a C++ `std::vector`.
   * \details If `GSMINRES_ASSUME_ABI_COMPATIBLE` is defined, this performs a `reinterpret_cast`;
   *          otherwise, it copies the data byte-wise using `std::memcpy`.
   *
   * \param[in] src Pointer to input C-style complex array (size = n).
   * \param[in] n   Number of elements.
   * \return `std::vector<std::complex<double>>` of size `n`.  
   */
  inline std::vector<std::complex<double>> to_cpp_vector(const double _Complex* src, std::size_t n) {
    std::vector<std::complex<double>> dst(n);
    // Check std::complex<double> is "trivially copyable"
    static_assert(std::is_trivially_copyable<std::complex<double>>::value,
                  "Cannot safely memcpy std::complex<double> — it must be trivially copyable");
    std::memcpy(dst.data(), src, n*sizeof(std::complex<double>));
    return dst;
  }
#endif

  /**
   * \brief Copy a `std::vector<std::complex<double>>` to a C-style array of `double _Complex`.
   * \param[in]  src Source C++ vector.
   * \param[out] dst Destination pointer to C-style complex array.
   */
  inline void from_cpp_vector(const std::vector<std::complex<double>>& src, double _Complex* dst) {
    // Check std::complex<double> is "trivially copyable"
    static_assert(std::is_trivially_copyable<std::complex<double>>::value,
                  "Cannot safely memcpy std::complex<double> — it must be trivially copyable.");
    std::memcpy(dst, src.data(), src.size()*sizeof(std::complex<double>));
  }

} // namespace gsminres_c_api_util

#endif // GSMINRES_C_API_UTIL_HPP

/*
 if std::complex<double> is not trivially copyable:
 CHANGE
 std::memcpy(dst.data(), src, n*sizeof(std::complex<double>));
 TO
 for (std::size_t i=0; i<n; ++i){
   dst[i] = src[i];
 }
 */
