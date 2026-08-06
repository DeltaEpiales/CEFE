cmake_minimum_required(VERSION 3.14)
project(cefe VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Fetch Pybind11
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.11.1
)
FetchContent_MakeAvailable(pybind11)

# Fetch Eigen3
set(BUILD_TESTING OFF CACHE BOOL "Disable tests" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "Disable Eigen tests" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "Disable Eigen docs" FORCE)
set(EIGEN_BUILD_BLAS OFF CACHE BOOL "Disable Eigen BLAS" FORCE)
set(EIGEN_BUILD_LAPACK OFF CACHE BOOL "Disable Eigen LAPACK" FORCE)
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)

# Fetch Spectra
FetchContent_Declare(
  spectra
  GIT_REPOSITORY https://github.com/yixuan/spectra.git
  GIT_TAG v1.0.1
)
FetchContent_MakeAvailable(spectra)

# Includes
include_directories(include)
include_directories(${eigen_SOURCE_DIR})
include_directories(${spectra_SOURCE_DIR}/include)

# Python module
pybind11_add_module(cefe_core
    src/python/pybind_module.cpp
    src/geometry/CausalDiamond.cpp
    src/core/EntropicFieldEngine.cpp
    dynamic_tensors/DynamicTensorEngine.cpp
)

target_include_directories(cefe_core PUBLIC include)
target_link_libraries(cefe_core PRIVATE pybind11::module)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cefe_core PRIVATE OpenMP::OpenMP_CXX)
endif()
