
add_library(MKL_SEQ INTERFACE) # Single thread
add_library(MKL_OMP INTERFACE) # Multithread OpenMP
add_library(MKL_TBB INTERFACE) # Multithread Thread Building Blocks
message(STATUS "ONEAPI_ROOT is set to: $ENV{ONEAPI_ROOT}")
set(_mkl_lib "$ENV{ONEAPI_ROOT}/mkl/latest/lib/")

target_link_libraries(MKL_SEQ INTERFACE
   ${_mkl_lib}mkl_blas95_lp64.lib
   ${_mkl_lib}mkl_lapack95_lp64.lib
   ${_mkl_lib}mkl_intel_lp64.lib
   ${_mkl_lib}mkl_core.lib
   ${_mkl_lib}mkl_sequential.lib
)
target_include_directories(MKL_SEQ INTERFACE
   "$ENV{ONEAPI_ROOT}/mkl/latest/include/intel64/lp64"
   "$ENV{ONEAPI_ROOT}/mkl/latest/include"
)
# ------------------------------------------ THIS IS FULL DYNAMIC LIBS
target_link_libraries(MKL_OMP INTERFACE
   ${_mkl_lib}mkl_blas95_lp64.lib
   ${_mkl_lib}mkl_lapack95_lp64.lib
   ${_mkl_lib}mkl_intel_lp64.lib
   ${_mkl_lib}mkl_core.lib
   ${_mkl_lib}mkl_intel_thread.lib
   "$ENV{ONEAPI_ROOT}/compiler/latest/lib"
   "$ENV{ONEAPI_ROOT}/compiler/latest/bin"
)
target_include_directories(MKL_OMP INTERFACE
   #"$ENV{ONEAPI_ROOT}/mkl/include/intel64/lp64"
   "$ENV{ONEAPI_ROOT}/mkl/latest/include"
)
target_link_libraries(MKL_TBB INTERFACE
    ${_mkl_lib}mkl_blas95_lp64.lib
    ${_mkl_lib}mkl_lapack95_lp64.lib
    ${_mkl_lib}mkl_intel_lp64.lib
    ${_mkl_lib}mkl_core.lib
    ${_mkl_lib}mkl_tbb_thread.lib
    "$ENV{ONEAPI_ROOT}/tbb/lib/intel64_win/vc_mt/tbb.lib"
)
target_include_directories(MKL_TBB INTERFACE
   "$ENV{ONEAPI_ROOT}/mkl/include/intel64/lp64"
   "$ENV{ONEAPI_ROOT}/mkl/include"
)

