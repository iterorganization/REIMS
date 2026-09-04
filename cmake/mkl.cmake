
add_library(MKL_SEQ INTERFACE) # Single thread
add_library(MKL_OMP INTERFACE) # Multithread OpenMP
add_library(MKL_TBB INTERFACE) # Multithread Thread Building Blocks
message(STATUS "ONEAPI_ROOT is set to: $ENV{ONEAPI_ROOT}")
set(_mkl_search_paths
   "$ENV{ONEAPI_ROOT}/mkl/latest/lib"
   "$ENV{ONEAPI_ROOT}/mkl/latest/lib/intel64"
   "$ENV{ONEAPI_ROOT}/mkl/latest/lib/intel64_win"
)

function(_find_mkl_library out_var library_name)
   find_file(${out_var}
      NAMES ${library_name}.lib ${library_name}
      PATHS ${_mkl_search_paths}
      NO_DEFAULT_PATH
   )
   if(NOT ${out_var})
      message(FATAL_ERROR "Could not find Intel oneAPI MKL library '${library_name}' in: ${_mkl_search_paths}")
   endif()
endfunction()

_find_mkl_library(_mkl_intel_lp64 mkl_intel_lp64)
_find_mkl_library(_mkl_core mkl_core)
_find_mkl_library(_mkl_sequential mkl_sequential)
_find_mkl_library(_mkl_intel_thread mkl_intel_thread)
_find_mkl_library(_mkl_tbb_thread mkl_tbb_thread)

target_link_libraries(MKL_SEQ INTERFACE
   ${_mkl_intel_lp64}
   ${_mkl_core}
   ${_mkl_sequential}
)
target_include_directories(MKL_SEQ INTERFACE
   "$ENV{ONEAPI_ROOT}/mkl/latest/include/intel64/lp64"
   "$ENV{ONEAPI_ROOT}/mkl/latest/include"
)
# ------------------------------------------ THIS IS FULL DYNAMIC LIBS
target_link_libraries(MKL_OMP INTERFACE
   ${_mkl_intel_lp64}
   ${_mkl_core}
   ${_mkl_intel_thread}
)
target_include_directories(MKL_OMP INTERFACE
   #"$ENV{ONEAPI_ROOT}/mkl/include/intel64/lp64"
   "$ENV{ONEAPI_ROOT}/mkl/latest/include"
)
target_link_libraries(MKL_TBB INTERFACE
    ${_mkl_intel_lp64}
    ${_mkl_core}
    ${_mkl_tbb_thread}
    "$ENV{ONEAPI_ROOT}/tbb/lib/intel64_win/vc_mt/tbb.lib"
)
target_include_directories(MKL_TBB INTERFACE
   "$ENV{ONEAPI_ROOT}/mkl/include/intel64/lp64"
   "$ENV{ONEAPI_ROOT}/mkl/include"
)
