
add_library(MKL_SEQ INTERFACE) # Single thread
add_library(MKL_OMP INTERFACE) # Multithread OpenMP
add_library(MKL_TBB INTERFACE) # Multithread Thread Building Blocks
message(STATUS "ONEAPI_ROOT is set to: $ENV{ONEAPI_ROOT}")
set(_mkl_search_paths
   "$ENV{ONEAPI_ROOT}/mkl/latest/lib"
   "$ENV{ONEAPI_ROOT}/mkl/latest/lib/intel64"
   "$ENV{ONEAPI_ROOT}/mkl/latest/lib/intel64_win"
)
set(_compiler_search_paths
   "$ENV{ONEAPI_ROOT}/compiler/latest/lib"
   "$ENV{ONEAPI_ROOT}/compiler/latest/lib/intel64"
   "$ENV{ONEAPI_ROOT}/compiler/latest/lib/intel64_win"
)
set(_tbb_search_paths
   "$ENV{ONEAPI_ROOT}/tbb/latest/lib/intel64/vc14"
   "$ENV{ONEAPI_ROOT}/tbb/latest/lib/intel64_win/vc14"
   "$ENV{ONEAPI_ROOT}/tbb/latest/lib/intel64_win/vc_mt"
   "$ENV{ONEAPI_ROOT}/tbb/lib/intel64/vc14"
   "$ENV{ONEAPI_ROOT}/tbb/lib/intel64_win/vc14"
   "$ENV{ONEAPI_ROOT}/tbb/lib/intel64_win/vc_mt"
)

function(_find_oneapi_library out_var library_name search_paths_var)
   find_file(${out_var}
      NAMES ${library_name}.lib ${library_name}
      PATHS ${${search_paths_var}}
      NO_DEFAULT_PATH
   )
   if(NOT DEFINED ${out_var} OR NOT EXISTS "${${out_var}}")
      message(FATAL_ERROR "Could not find Intel oneAPI library '${library_name}' in: ${${search_paths_var}}")
   endif()
endfunction()

_find_oneapi_library(_mkl_intel_lp64 mkl_intel_lp64 _mkl_search_paths)
_find_oneapi_library(_mkl_core mkl_core _mkl_search_paths)
_find_oneapi_library(_mkl_sequential mkl_sequential _mkl_search_paths)
_find_oneapi_library(_mkl_intel_thread mkl_intel_thread _mkl_search_paths)
_find_oneapi_library(_mkl_tbb_thread mkl_tbb_thread _mkl_search_paths)
_find_oneapi_library(_libiomp5md libiomp5md _compiler_search_paths)
_find_oneapi_library(_tbb_lib tbb _tbb_search_paths)

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
   ${_libiomp5md}
)
target_include_directories(MKL_OMP INTERFACE
   #"$ENV{ONEAPI_ROOT}/mkl/include/intel64/lp64"
   "$ENV{ONEAPI_ROOT}/mkl/latest/include"
)
target_link_libraries(MKL_TBB INTERFACE
    ${_mkl_intel_lp64}
    ${_mkl_core}
    ${_mkl_tbb_thread}
    ${_tbb_lib}
)
target_include_directories(MKL_TBB INTERFACE
   "$ENV{ONEAPI_ROOT}/mkl/latest/include/intel64/lp64"
   "$ENV{ONEAPI_ROOT}/mkl/latest/include"
)
