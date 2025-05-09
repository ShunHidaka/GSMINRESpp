cmake_minimum_required(VERSION 3.14)

# =====================================
# Language configureation
# =====================================
include(CheckLanguage)
check_language(C)
check_language(Fortran)

# C API / Fortran Interface option (default ON)
option(GSMINRES_ENABLE_C_API "Enable C API" ON)
option(GSMINRES_ENABLE_FORTRAN_INTERFACE "Enable Fortran interface (requires C API)" ON)

# Disable C API if C compiler not found
if(NOT CMAKE_C_COMPILER AND GSMINRES_ENABLE_C_API)
  message(WARNING "C compiler not found. Disabling C API.")
  set(GSMINRES_ENABLE_C_API OFF)
endif()

# Disable Fortran interface if Fortran compiler is missing OR C API is disabled
if((NOT CMAKE_Fortran_COMPILER OR NOT GSMINRES_ENABLE_C_API) AND GSMINRES_ENABLE_FORTRAN_INTERFACE)
  message(WARNING "Fortran interface requires both Fortran and C compilers. Disabling Fortran interface.")
  set(GSMINRES_ENABLE_FORTRAN_INTERFACE OFF)
endif()

# Define language list
set(LANGUAGE_LIST CXX)
if(GSMINRES_ENABLE_C_API)
  list(APPEND LANGUAGE_LIST C)
endif()
if(GSMINRES_ENABLE_FORTRAN_INTERFACE)
  list(APPEND LANGUAGE_LIST Fortran)
endif()

# =====================================
# Project declaration
# =====================================
project(GSMINRESpp VERSION 1.0.0 LANGUAGES ${LANGUAGE_LIST})

# =====================================
# Compiler setting
# =====================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG         "-g3 -O0 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE       "-O3 -DNDEBUG")

if(GSMINRES_ENABLE_C_API)
  set(CMAKE_C_STANDARD 99)
  set(CMAKE_C_STANDARD_REQUIRED ON)
  set(CMAKE_C_FLAGS_DEBUG         "-g3 -O0 -Wall -Wextra -Wpedantic")
  set(CMAKE_C_FLAGS_RELEASE       "-O3 -DNDEBUG")
endif()

if(GSMINRES_ENABLE_FORTRAN_INTERFACE)
  set(CMAKE_Fortran_STANDARD 2017)
  set(CMAKE_Fortran_STANDARD_REQUIRED ON)
  set(CMAKE_Fortran_FLAGS_DEBUG   "-g3 -O0 -Wall")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()

# =====================================
# Install prefix configuration
# =====================================
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/gsminres_install" CACHE PATH "Install path" FORCE)
  message(STATUS "Install prefix not specified. Using default: ${CMAKE_INSTALL_PREFIX}")
else()
  message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
endif()

# =====================================
# OpenMP setup
# =====================================
option(USE_OPENMP "Enable OpenMP parallelization" ON)
set(OPENMP_CXX_OPTION     "")
set(OPENMP_C_OPTION       "")
set(OPENMP_FORTRAN OPTION "")

if(USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found, enabling parallelization.")
    set(OPENMP_CXX_OPTION OpenMP::OpenMP_CXX)
    if(OpenMP_C_FOUND)
      set(OPENMP_C_OPTION OpenMP::OpenMP_C)
    endif()
    if(OpenMP_Fortran_FOUND)
      set(OPENMP_Fortran_OPTION OpenMP::OpenMP_Fortran)
    endif()
  else()
    message(WARNING "OpenMP not found, disabling parallelization.")
  endif()
else()
  message(STATUS "OpenMP disabled.")
endif()

# =====================================
# ABI compatibility option
# =====================================
option(GSMINRES_ASSUME_ABI_COMPATIBLE
  "Assume ABI compatibility between std::complex<double> and double _Complex"
  OFF)
if(GSMINRES_ASSUME_ABI_COMPATIBLE)
  add_compile_definitions(GSMINRES_ASSUME_ABI_COMPATIBLE)
endif()

# =====================================
# BLAS / LAPACK configuration
# =====================================
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

# =====================================
# Source files
# =====================================
set(SRC src/gsminres_solver.cpp
        src/gsminres_util.cpp)

# Optionally add C API / Fortran Interface
if(GSMINRES_ENABLE_C_API)
  list(APPEND SRC src/gsminres_c_api.cpp)
endif()
if(GSMINRES_ENABLE_FORTRAN_INTERFACE)
  list(APPEND SRC src/gsminres_fortran_interface.f90)
endif()

# =====================================
# Build Targets
# =====================================
# Shared library (libgsminres.so)
include_directories(include)
add_library(gsminres_shared SHARED ${SRC})
target_link_libraries(gsminres_shared PRIVATE ${BLAS_LIBRARIES}
                                              ${LAPACK_LIBRARIES}
                                              ${OPENMP_CXX_OPTION})
set_target_properties(gsminres_shared PROPERTIES
                      OUTPUT_NAME gsminres
                      VERSION     ${PROJECT_VERSION}
                      SOVERSION   1)
# Static Library (libgsminres.a)
add_library(gsminres_static STATIC ${SRC})
target_link_libraries(gsminres_static PRIVATE ${BLAS_LIBRARIES}
                                              ${LAPACK_LIBRARIES}
					      ${OPENMP_CXX_OPTION})
set_target_properties(gsminres_shared PROPERTIES OUTPUT_NAME gsminres)

# =====================================
# Sample programs
# =====================================
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
  target_link_libraries(sample2_c PRIVATE gsminres_shared m
                                                          ${BLAS_LIBRARIES}
                                                          ${LAPACK_LIBRARIES}
                                                          ${OPENMP_C_OPTION})
endif()

# Fortran Interface sample program
if(GSMINRES_ENABLE_FORTRAN_INTERFACE)
  add_executable(sample1_f sample/sample1_f.f90)
  set_target_properties(sample1_f PROPERTIES LINKER_LANGUAGE Fortran)
  target_include_directories(sample1_f PRIVATE include)
  target_link_libraries(sample1_f PRIVATE gsminres_shared ${BLAS_LIBRARIES}
                                                          ${LAPACK_LIBRARIES}
                                                          ${OPENMP_FORTRAN_OPTION})
endif()

# =====================================
# Installation
# =====================================
install(TARGETS gsminres_shared gsminres_static
        EXPORT  gsminresTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

# Generate CMake Package
include(CMakePackageConfigHelpers)
install(EXPORT gsminresTargets
        FILE gsminresTargets.cmake
        NAMESPACE gsminres::
        DESTINATION lib/cmake/gsminres)
configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/gsminresConfig.cmake.in"
                              "${CMAKE_CURRENT_BINARY_DIR}/gsminresConfig.cmake"
                              INSTALL_DESTINATION lib/cmake/gsminres)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/gsminresConfig.cmake" DESTINATION lib/cmake/gsminres)
