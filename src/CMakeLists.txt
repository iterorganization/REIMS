# Files required by reims
add_executable(reims
	reims.f90
	krn_interface.f90
	krn_linear_system.f90
	krn_global_tools.f90
	krn_simulation.f90
	lib_ext_math.f90
	lib_input.f90
	lib_He_thermo.f90
	lib_friction_correlations.f90
	lib_nusselt_correlations.f90
	lib_material.f90
	lib_material_nb3sn.f90
	lib_material_nbti.f90
	lib_material_metal.f90
	lib_material_insulation.f90
	lib_hdf5_write.f90
	cmp_channel_init.f90
	cmp_channel_calc.f90
	cmp_channel_flux.f90
	cmp_channel_source_terms.f90
	cmp_strand_init.f90
	cmp_strand_flux.f90
	cmp_strand_source_terms.f90
	cmp_strand_calc.f90
	cmp_solid_init.f90
	cmp_solid_flux.f90	
	cmp_solid_calc.f90	
	cmp_mesh2D_init.f90
	cmp_mesh2D_flux.f90
	cmp_mesh2D_calc.f90
	cmp_mesh2D_hdf5_write.f90
	cmp_circulator_init.f90
	cmp_circulator_calc.f90
	cmp_junction_init.f90
	cmp_junction_calc.f90
	cmp_boundary_init.f90
	cmp_boundary_calc.f90
	cmp_FSlink_init.f90
	cmp_FSlink_calc.f90
	cmp_FFsrcLink_init.f90
	cmp_FFsrcLink_calc.f90
	cmp_SSsrcLink_init.f90
	cmp_SSsrcLink_calc.f90
	cmp_SSfluxLink.f90
)

# Libraries reqired by reims
target_link_libraries(reims
	MKL_OMP 
	OpenMP::OpenMP_Fortran
	fortran-yaml-c
)
#    HDF5::HDF5
#    hdf5::hdf5_fortran
#    hdf5::hdf5_hl_fortran
#    hdf5::hdf5_hl
#    hdf5::hdf5_hl_f90cstub
#    hdf5::hdf5_f90cstub
#    ${HDF5_ROOT}/lib/libhdf5_hl_f90cstub.lib
#    ${HDF5_ROOT}/lib/libhdf5_f90cstub.lib
    #${HDF5_ROOT}/lib/libhdf5_hl_fortran.lib
    #${HDF5_ROOT}/lib/libhdf5.lib
    #${HDF5_ROOT}/lib/libhdf5_fortran.lib
    #${HDF5_ROOT}/lib/libhdf5_hl.lib
#)
target_link_directories(reims PRIVATE "${HDF5_BUILD_DIR}/lib")
foreach(LIB_NAME ${HDF5_REQUIRED_LIBS})
    target_link_libraries(reims ${LIB_NAME})
endforeach()

get_cmake_property(_variableNames VARIABLES)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
