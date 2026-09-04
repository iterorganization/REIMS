add_library(MKL_SEQ INTERFACE) # Single thread
add_library(MKL_OMP INTERFACE) # Multithread OpenMP
add_library(MKL_TBB INTERFACE) # Multithread Thread Building Blocks

set(_mkl_root_candidates
   "$ENV{MKLROOT}"
   "$ENV{ONEAPI_ROOT}/mkl/latest"
)

set(_mkl_root "")
foreach(_candidate IN LISTS _mkl_root_candidates)
   if(_candidate AND EXISTS "${_candidate}")
     set(_mkl_root "${_candidate}")
     break()
   endif()
endforeach()

if(NOT _mkl_root)
   message(FATAL_ERROR "Could not locate an Intel oneMKL installation. Checked MKLROOT and ONEAPI_ROOT.")
endif()

message(STATUS "Using oneMKL from: ${_mkl_root}")

set(_mkl_lib_dirs
   "${_mkl_root}/lib/intel64"
   "${_mkl_root}/lib"
)

set(_mkl_include_dirs
   "${_mkl_root}/include/intel64/lp64"
   "${_mkl_root}/include"
)

function(_reims_find_mkl_library out_var lib_name)
   unset(_resolved_lib CACHE)
   unset(_resolved_lib)
   find_library(_resolved_lib
     NAMES ${lib_name} ${lib_name}.lib
     HINTS ${_mkl_lib_dirs}
     NO_DEFAULT_PATH
   )

   if(NOT _resolved_lib)
     message(FATAL_ERROR "Could not find required oneMKL library '${lib_name}' under ${_mkl_root}.")
   endif()

   set(${out_var} "${_resolved_lib}" PARENT_SCOPE)
endfunction()

function(_reims_link_mkl target)
   set(_resolved_libs "")
   foreach(_lib_name IN LISTS ARGN)
     _reims_find_mkl_library(_resolved_lib "${_lib_name}")
     list(APPEND _resolved_libs "${_resolved_lib}")
   endforeach()
   target_link_libraries(${target} INTERFACE ${_resolved_libs})
endfunction()

_reims_link_mkl(MKL_SEQ
   mkl_intel_lp64
   mkl_core
   mkl_sequential
)
target_include_directories(MKL_SEQ INTERFACE ${_mkl_include_dirs})

# ------------------------------------------ THIS IS FULL DYNAMIC LIBS
_reims_link_mkl(MKL_OMP
   mkl_intel_lp64
   mkl_core
   mkl_intel_thread
)
target_include_directories(MKL_OMP INTERFACE ${_mkl_include_dirs})

_reims_link_mkl(MKL_TBB
   mkl_intel_lp64
   mkl_core
   mkl_tbb_thread
)
target_include_directories(MKL_TBB INTERFACE ${_mkl_include_dirs})
