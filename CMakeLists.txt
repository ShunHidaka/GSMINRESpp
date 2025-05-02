cmake_minimum_required(VERSION 3.14)
project(GSMINRESpp VERSION 1.0.0 LANGUAGES CXX C)

# =======================
# Build Configuration
# =======================

# Default build is Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
# Set default installation Directory to home directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/gsminres_install" CACHE PATH "Default path" FORCE)
endif()
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

# Language standard settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG   "-g3 -O0 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Wextra -Wpedantic")
set(CMAKE_C_FLAGS_DEBUG   "-g3 -O0 -Wall -Wextra -Wpedantic")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Wextra -Wpedantic")

# =======================
# Options
# =======================

# C API (Default ON)
option(GSMINRES_ENABLE_C_API "Building with C API support" ON)
# Complex number ABI compatibility (Default OFF)
option(GSMINRES_ASSUME_ABI_COMPATIBLE
  "Assume std::complex<double> and double _Complex are ABI-compatible" OFF)
# OpenMP (Default ON)
option(USE_OPENMP "Enable OpenMP parallelization" ON)

# ABI compatibility
if(GSMINRES_ASSUME_ABI_COMPATIBLE)
  message(STATUS
    "GSMINRES: Assuming ABI compatibility between std::complex and double _Complex (using reinterpret_cast)")
  add_compile_definitions(GSMINRES_ASSUME_ABI_COMPATIBLE)
else()
  message(STATUS "GSMINRES: Using memcpy-based conversion for ABI safety")
endif()

# OpenMP setup
if(USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found, enabling parallelization.")
    set(OPENMP_CXX_OPTION OpenMP::OpenMP_CXX)
    set(OPENMP_C_OPTION   OpenMP::OpenMP_C)
  else()
    message(WARNING "OpenMP not found, disabling parallelization.")
    set(OPENMP_CXX_OPTION "")
    set(OPENMP_C_OPTION   "")
  endif()
else()
  message(STATUS "OpenMP disabled.")
  set(OPENMP_CXX_OPTION "")
  set(OPENMP_C_OPTION   "")
endif()

# =======================
# Dependencies
# =======================

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")

# Check OpenBLAS version
find_path(OPENBLAS_INCLUDE_DIR openblas_config.h)
if(OPENBLAS_INCLUDE_DIR)
  file(READ "${OPENBLAS_INCLUDE_DIR}/openblas_config.h" OPENBLAS_CONFIG_CONTENTS)
  string(REGEX MATCH "#define OPENBLAS_VERSION \"[^\"]+\"" OPENBLAS_VERSION_LINE "${OPENBLAS_CONFIG_CONTENTS}")
  string(REGEX REPLACE "#define OPENBLAS_VERSION \"([^\"]+)\"" "\\1" OPENBLAS_VERSION "${OPENBLAS_VERSION_LINE}")
  message(STATUS "Detected OpenBLAS version: ${OPENBLAS_VERSION}")
  string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" OPENBLAS_VERSION_NUMBER "${OPENBLAS_VERSION}")
  set(OPENBLAS_REQUIRED_VERSION "0.3.27") # Require Version
  if(OPENBLAS_VERSION_NUMBER VERSION_LESS OPENBLAS_REQUIRED_VERSION)
    message(FATAL_ERROR "Detected OpenBLAS version: ${OPENBLAS_VERSION_NUMBER}
                         Required minimum version: ${OPENBLAS_REQUIRED_VERSION}
                         Error: OpenBLAS is too old.
                         Please upgrade to at least ${OPENBLAS_REQUIRED_VERSION}.")
  else()
    message(STATUS "Detected OpenBLAS version: ${OPENBLAS_VERSION_NUMBER} (OK)")
  endif()
endif()

# =======================
# Source Files
# =======================

set(SRC src/gsminres_solver.cpp
        src/gsminres_util.cpp)

# Optionally add C API
if(GSMINRES_ENABLE_C_API)
  list(APPEND SRC src/gsminres_c_api.cpp)
endif()


# =======================
# Build Targets
# =======================

# Include path
include_directories(include)

# Shared Library (libgsminres.so)
add_library(gsminres_shared SHARED ${SRC})
target_link_libraries(gsminres_shared PRIVATE ${BLAS_LIBRARIES}
                                              ${LAPACK_LIBRARIES}
                                              ${OPENMP_CXX_OPTION})
set_target_properties(gsminres_shared PROPERTIES OUTPUT_NAME gsminres
                      OUTPUT_NAME gsminres
                      VERSION     ${PROJECT_VERSION}
                      SOVERSION   1)

# Static Library (libgsminres.a)
add_library(gsminres_static STATIC ${SRC})
target_link_libraries(gsminres_static PRIVATE ${BLAS_LIBRARIES}
                                              ${LAPACK_LIBRARIES}
                                              ${OPENMP_CXX_OPTION})
set_target_properties(gsminres_static PROPERTIES OUTPUT_NAME gsminres)

# =======================
# Sample Executables
# =======================

add_executable(sample1 sample/sample1.cpp)
add_executable(sample2 sample/sample2.cpp)

# Using shared library
target_link_libraries(sample1 PRIVATE gsminres_shared ${BLAS_LIBRARIES}
                                                      ${LAPACK_LIBRARIES}
                                                      ${OPENMP_CXX_OPTION})
target_link_libraries(sample2 PRIVATE gsminres_shared ${BLAS_LIBRARIES}
                                                      ${LAPACK_LIBRARIES}
                                                      ${OPENMP_CXX_OPTION})
# Using static library
#target_link_libraries(sample1 PRIVATE gsminres_static ${BLAS_LIBRARIES}
#                                                      ${LAPACK_LIBRARIES}
#                                                      ${OPENMP_CXX_OPTION})
#target_link_libraries(sample2 PRIVATE gsminres_static ${BLAS_LIBRARIES}
#                                                      ${LAPACK_LIBRARIES}
#                                                      ${OPENMP_CXX_OPTION})

# C API sample program
if(GSMINRES_ENABLE_C_API)
  add_executable(sample2_c sample/sample2_c.c)
  set_target_properties(sample2_c PROPERTIES LINKER_LANGUAGE C)
  target_include_directories(sample2_c PRIVATE include)
  target_link_libraries(sample2_c PRIVATE gsminres_shared ${BLAS_LIBRARIES}
                                                          ${LAPACK_LIBRARIES}
                                                          ${OPENMP_C_OPTION})
endif()

# =======================
# Installation
# =======================

install(TARGETS gsminres_shared gsminres_static
        EXPORT  gsminresTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.hpp")

# Generate CMake Package
include(CMakePackageConfigHelpers)
install(EXPORT gsminresTargets
        FILE gsminresTargets.cmake
        NAMESPACE gsminres::
        DESTINATION lib/cmake/gsminres)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/gsminresConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/gsminresConfig.cmake"
  INSTALL_DESTINATION lib/cmake/gsminres)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/gsminresConfig.cmake" DESTINATION lib/cmake/gsminres)
