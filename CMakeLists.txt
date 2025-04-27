cmake_minimum_required(VERSION 3.14)
project(GSMINRESpp CXX)

# Set Install Directory (Default: HOME_directory)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/gsminres_install" CACHE PATH "default path" FORCE)
endif()

# C++17
set(CMAKE_XX_STANDARD 17)
set(CMAKE_STANDARD_REQUIRED ON)

# Include path
include_directories(include)

# Find BLAS/LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# OpenBLAS Version Check Start
find_path(OPENBLAS_INCLUDE_DIR openblas_config.h)
if(OPENBLAS_INCLUDE_DIR)
  file(READ "${OPENBLAS_INCLUDE_DIR}/openblas_config.h" OPENBLAS_CONFIG_CONTENTS)
  string(REGEX MATCH "#define OPENBLAS_VERSION \"[^\"]+\"" OPENBLAS_VERSION_LINE "${OPENBLAS_CONFIG_CONTENTS}")
  string(REGEX REPLACE "#define OPENBLAS_VERSION \"([^\"]+)\"" "\\1" OPENBLAS_VERSION "${OPENBLAS_VERSION_LINE}")
  message(STATUS "Detected OpenBLAS version: ${OPENBLAS_VERSION}")
  string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" OPENBLAS_VERSION_NUMBER "${OPENBLAS_VERSION}")
  # Require Version
  set(OPENBLAS_REQUIRED_VERSION "0.3.27")
  if(OPENBLAS_VERSION_NUMBER VERSION_LESS OPENBLAS_REQUIRED_VERSION)
    message(FATAL_ERROR "OpenBLAS version ${OPENBLAS_VERSION_NUMBER} is too old. Please upgrade to at least ${OPENBLAS_REQUIRED_VERSION}.")
  endif()
else()
  message(STATUS "OpenBLAS not detected, other BLAS/LAPACK are used.")
endif()
# OpenBLAS Version Check End

# Source files
set(SRC
   src/gsminres_solver.cpp
   src/gsminres_util.cpp
)

# Create Shared Library
add_library(gsminres_shared SHARED ${SRC})
target_link_libraries(gsminres_shared PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set_target_properties(gsminres_shared PROPERTIES OUTPUT_NAME gsminres
                      OUTPUT_NAME gsminres
                      VERSION     1.0.0
                      SOVERSION   1)
install(TARGETS gsminres_shared
        EXPORT  gsminresTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

# Create Static Library
add_library(gsminres_static STATIC ${SRC})
target_link_libraries(gsminres_static PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set_target_properties(gsminres_static PROPERTIES OUTPUT_NAME gsminres)
install(TARGETS gsminres_static
        EXPORT  gsminresTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

# Sample executable (using shared library)
add_executable(sample1 sample/sample1.cpp)
add_executable(sample2 sample/sample2.cpp)
target_link_libraries(sample1 PRIVATE gsminres_shared ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(sample2 PRIVATE gsminres_shared ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#target_link_libraries(sample1 PRIVATE gsminres_static ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#target_link_libraries(sample2 PRIVATE gsminres_static ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Install Setting
install(TARGETS gsminres_shared gsminres_static
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
  INSTALL_DESTINATION lib/cmake/gsminres
  )
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/gsminresConfig.cmake" DESTINATION lib/cmake/gsminres)

