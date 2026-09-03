cmake_minimum_required(VERSION 3.14...4.2)

if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
endif()

set(CMAKE_POLICY_VERSION_MINIMUM "3.5")

project(PROJECT_REIMS
    VERSION 2.0
    DESCRIPTION "Riemann Explicit Implicit Magnet Simulator"
    LANGUAGES Fortran C)
	
# set(CMAKE_Fortran_COMPILER "K:/Tools/PortableIntel/compiler/latest/bin/ifort.exe")
# set(CMAKE_C_COMPILER "K:/Tools/VisualStudioCompiler_14.43.34808_22621/VC/Tools/MSVC/14.43.34808/bin/Hostx64/x64/cl.exe")

# Main configuration
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
include(cmake/CPM.cmake)
include(cmake/flags.cmake)

# Libraries and cmake configurations
include(cmake/mkl.cmake)
set(HDF5_USE_STATIC_LIBRARIES ON)
set(HDF5_REQUIRED_LIBS HDF5::HDF5 libhdf5 libhdf5_hl libhdf5_f90cstub libhdf5_fortran libhdf5_hl_f90cstub libhdf5_hl_fortran)
find_package(HDF5 REQUIRED COMPONENTS C HL Fortran Fortran_HL)
find_package(OpenMP)
CPMAddPackage("gh:Nicholaswogan/fortran-yaml-c#main")

# Source code
add_subdirectory(src)
# add_subdirectory(src_tests)
