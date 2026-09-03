if (MSVC)
    # Set stack size for MSVC
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:640000000")
    # doesn't work but it supposed to get rid of intel runtime
    # even hdf should be compiled with this flags so, more work....
    #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  /libs:static /threads /Qopenmp:static")
    #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /Qopenmp /MD /link /NODEFAULTLIB:libcmt /NODEFAULTLIB:libifcoremd /NODEFAULTLIB:libmmd /NODEFAULTLIB:libiomp5md /LIBPATH:\"${IFORT_COMPILER_LIB}\" libifcoremt.lib libmmt.lib libiomp5mt.lib") 
elseif (CMAKE_COMPILER_IS_GNUCXX)
    # Set stack size for GCC/MinGW
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,640000000")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size:64")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -NODEFAULTLIB:MSVCRT").0
# ------------------------------------------------------------------------------- THIS IS FULL DYNAMIC LIBS
#set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -warn:all -check:all /Qipo")# -threads -Qmkl:parallel -libs:dll") #-Qmkl:sequential
#set(CMAKE_Fortran_FLAGS_RELEASE "/Qprec-div- /O3")#"/Qprec-div- /Qipo /O3 /fp:fast=2 /arch:CORE-AVX2")#"${CMAKE_Fortran_FLAGS_RELEASE} -O3")
set(CMAKE_Fortran_FLAGS_RELEASE "/Qprec-div- /O3")
#set(CMAKE_Fortran_FLAGS_RELEASE "/fp:precise /Ofast /arch:CORE-AVX2")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -warn:all -check:all") 


# NOTES for the futue are below:

# /O3 (/Ofast) + /Qipo --> crash (already seen in 2021 for CS simulation)* 
# /fp:fast=2 (2 or 1) decreases the performance
# /arch:CORE-AVX2 TO BE CHECKED
# /fp:precise TO BE CHECKED 

# 
# aditional info:
#    compiler command:
#      /4I8 /module:"%MKLROOT%"\include\intel64/ilp64 -I"%MKLROOT%"\include
#    sufix: _lp64 (32bit integer) _ilp64 + /4I8 (32bit integer)
#

# To do something with it
# set_target_properties(app PROPERTIES
#                      RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR}
#                      RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR})

#/add_custom_command(TARGET app POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:app> ${CMAKE_SOURCE_DIR})

# TIP:
#    Add to "cmake-tools-kits.json" in Visual Studio Kits:
#       "compilers": {
#          "Fortran": "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/bin/intel64/ifort.exe"
#       },
#

# /warn:interfaces
# /real_size:64
# /traceback
# /check:bounds
# /check:stack
# /Qmkl:sequential

# /fpp
# /I"C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2018.5.274\windows\\mkl\include\intel64\lp64"
# /I"C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2018.5.274\windows\\mkl\include"
# /DCMAKE_INTDIR=\"Debug\"
# /W1
# /Qlocation,link,"C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\Common7\IDE\VC\\bin\amd64"
# /Qm64 "C:\Users\kosekj\source\repos\channel\src\Implicit.f90"
