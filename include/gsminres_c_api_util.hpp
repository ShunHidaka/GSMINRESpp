#ifndef GSMINRES_C_API_UTIL_HPP
#define GSMINRES_C_API_UTIL_HPP

#include <complex>
#include <vector>
#include <cstring>
#include <type_traits>

namespace gsminres_c_api_util {

#ifdef GSMINRES_ASSUME_ABI_COMPATIBLE
  inline std::vector<std::complex<double>> to_cpp_vector(const double _Complex* src, std::size_t n) {
    return std::vector<std::complex<double>>(reinterpret_cast<const std::complex<double>*>(src),
                                             reinterpret_cast<const std::complex<double>*>(src) + n);
  }
#else // If not ABI Compatibility between "std::complex<double>" and "double _Complex"
  inline std::vector<std::complex<double>> to_cpp_vector(const double _Complex* src, std::size_t n) {
    std::vector<std::complex<double>> dst(n);
    // Check std::complex<double> is "trivially copyable"
    static_assert(std::is_trivially_copyable<std::complex<double>>::value,
                  "Cannot safely memcpy std::complex<double> — it must be trivially copyable");
    std::memcpy(dst.data(), src, n*sizeof(std::complex<double>));
    return dst;
  }
#endif

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
