! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_mesh2D_hdf5_write_m       ! TODO : to be removed
    use cmp_mesh2D_init_m
    implicit none

contains

subroutine writing_HDF5_2D_casing(me,time,step)

  type(mesh2D_t),     intent(inout) :: me
  real(dp), intent(in) :: time
  integer(i_64), intent(in) :: step
           
  integer(hsize_t), dimension(1) :: maxdims_nbcel,maxdims_nbpoint,maxdims_nbconIds,maxdims_Offsets
  integer(hsize_t), dimension(1) :: maxdims_Types,maxdims_Connectivity,maxdims_Temperature,maxdimsVersion,maxdims_Time
  integer(hsize_t), dimension(1:2) :: maxdims_Points
  integer :: error
  integer(i_64) :: Nb2DCells,NbPoints,NbConIds
  integer(i_64), dimension(2) :: VersionData
  integer(i_64), dimension(:), allocatable :: Offsets,Connectivity
  integer(i_8), dimension(:), allocatable :: Types
  real(sp), dimension(:,:), allocatable :: Points
  real(dp), dimension(:), allocatable :: Temperature,TimeArr
  integer(hsize_t), dimension(1) :: dims_nbcel,dims_nbpoint,dims_nbconIds,dims_Offsets,dims_Types
  integer(hsize_t), dimension(1) :: dims_Connectivity,dims_Temperature,dims_Time
  integer(hsize_t), dimension(1:2) :: dims_Points
  integer(hid_t) :: nbcel_id,nbcel_space_id,nbpoint_id,nbpoint_space_id,nbconIds_id,Connectivity_id,Connectivity_space_id
  integer(hid_t) :: nbconIds_space_id,Offsets_id,Offsets_space_id,Types_id,Types_space_id,Points_id,Points_space_id
  integer(hid_t) :: Temperature_id,Temperature_space_id,Time_id,Time_space_id
  integer(hid_t) :: gr_id, grA_id, group_id
  integer :: nb_e, nb_n, idx_offset
  integer(hid_t)  ::  dtype,dataspace_id,dataspaceV_id,attrA_id,attrV_id
  integer(size_t) ::  size
  integer(hsize_t), dimension(1) :: dimsType,dimsVersion
  character(len=3) :: stepChar

    if(step==0) then
      call h5fcreate_f(me%filename_extless//".hdf", H5F_ACC_TRUNC_F, me%f_id_writing, error)  ! H5F_ACC_TRUNC_F deletes the file if already existing
    endif

    write(stepChar,'(i3.3)') step
    call h5gcreate_f(me%f_id_writing, me%filename_extless//"_"//stepChar//"", group_id, error)

    call h5gcreate_f(group_id, "VTKHDF", gr_id, error)

    ! Type (attribute)
    dimsType = (/1/)
    call H5Tcopy_f(H5T_C_S1,dtype,error)
    size = 16
    call H5Tset_size_f (dtype, size, error)
    call h5screate_f(H5S_SCALAR_F, dataspace_id, error)
    call h5acreate_f(gr_id, 'Type', dtype, dataspace_id, attrA_id, error)
    call h5awrite_f(attrA_id, dtype, "UnstructuredGrid", dimsType, error)
    call h5aclose_f(attrA_id, error)
    call h5sclose_f(dataspace_id, error)

    ! Version (attribute)
    VersionData(1) = 1
    VersionData(2) = 0
    dimsVersion = (/2/)
    maxdimsVersion=dimsVersion
    call h5screate_simple_f(1, dimsVersion, dataspaceV_id, error, maxdimsVersion)
    call h5acreate_f(gr_id, 'Version', H5T_STD_I64LE, dataspaceV_id, attrV_id, error)
    call h5awrite_f(attrV_id, H5T_STD_I64LE, VersionData, dimsVersion, error)
    call h5aclose_f(attrV_id, error)
    call h5sclose_f(dataspaceV_id, error)

    ! NumberOfCells
    Nb2DCells=me%M2D_Prop%nb_elements
    dims_nbcel=(/1/)
    maxdims_nbcel=dims_nbcel
    call h5screate_simple_f(1, dims_nbcel, nbcel_space_id, error, maxdims_nbcel)
    call h5dcreate_f(gr_id, "NumberOfCells", H5T_STD_I64LE, nbcel_space_id, nbcel_id, error)
    call h5dwrite_f(nbcel_id, H5T_STD_I64LE, Nb2DCells, dims_nbcel, error)
    call h5sclose_f(nbcel_space_id, error)

    ! NumberOfPoints
    NbPoints=me%M2D_Prop%nb_nodes
    dims_nbpoint=(/1/)
    maxdims_nbpoint=dims_nbpoint
    call h5screate_simple_f(1, dims_nbpoint, nbpoint_space_id, error, maxdims_nbpoint)
    call h5dcreate_f(gr_id, "NumberOfPoints", H5T_STD_I64LE, nbpoint_space_id, nbpoint_id, error)
    call h5dwrite_f(nbpoint_id, H5T_STD_I64LE, NbPoints, dims_nbpoint, error)
    call h5sclose_f(nbpoint_space_id, error)

    ! NumberOfConnectivityIds
    NbConIds=me%M2D_Prop%nb_elements*3 ! 3 points defining a cell
    dims_nbconIds=(/1/)
    maxdims_nbconIds=dims_nbconIds
    call h5screate_simple_f(1, dims_nbconIds, nbconIds_space_id, error, maxdims_nbconIds)
    call h5dcreate_f(gr_id, "NumberOfConnectivityIds", H5T_STD_I64LE, nbconIds_space_id, nbconIds_id, error)
    call h5dwrite_f(nbconIds_id, H5T_STD_I64LE, NbConIds, dims_nbconIds, error)
    call h5sclose_f(nbconIds_space_id, error)

    ! Offsets
    allocate(Offsets(me%M2D_Prop%nb_elements+1))
    Offsets(1)=0
    do nb_e=2,me%M2D_Prop%nb_elements+1
      Offsets(nb_e)=Offsets(nb_e-1)+3 ! 3 points defining a cell
    enddo
    dims_Offsets=(/me%M2D_Prop%nb_elements+1/)
    maxdims_Offsets=dims_Offsets
    call h5screate_simple_f(1, dims_Offsets, Offsets_space_id, error, maxdims_Offsets)
    call h5dcreate_f(gr_id, "Offsets", H5T_STD_I64LE, Offsets_space_id, Offsets_id, error)
    call h5dwrite_f(Offsets_id, H5T_STD_I64LE, Offsets, dims_Offsets, error)
    call h5sclose_f(Offsets_space_id, error)
    deallocate(Offsets)

    ! Types
    allocate(Types(me%M2D_Prop%nb_elements))
    do nb_e=1,me%M2D_Prop%nb_elements
      Types(nb_e)=5
    enddo
    dims_Types=(/me%M2D_Prop%nb_elements/)
    maxdims_Types=dims_Types
    call h5screate_simple_f(1, dims_Types, Types_space_id, error, maxdims_Types)
    call h5dcreate_f(gr_id, "Types", H5T_STD_U8LE, Types_space_id, Types_id, error)
    call h5dwrite_f(Types_id, H5T_STD_U8LE, Types, dims_Types, error)
    call h5sclose_f(Types_space_id, error)
    deallocate(Types)
    
    ! Points
    allocate(Points(3,me%M2D_Prop%nb_nodes))
    do nb_n=1,me%M2D_Prop%nb_nodes
      Points(1,nb_n)=me%M2D_Prop%Node(nb_n)%x_coord
      Points(2,nb_n)=me%M2D_Prop%Node(nb_n)%y_coord
      Points(3,nb_n)=0
    enddo
    dims_Points=(/3, me%M2D_Prop%nb_nodes/)
    maxdims_Points=dims_Points
    call h5screate_simple_f(2, dims_Points, Points_space_id, error, maxdims_Points)
    call h5dcreate_f(gr_id, "Points", H5T_IEEE_F32LE, Points_space_id, Points_id, error)
    call h5dwrite_f(Points_id, H5T_IEEE_F32LE, Points, dims_Points, error)
    call h5sclose_f(Points_space_id, error)
    deallocate(Points)

    ! Connectivity
    allocate(Connectivity(NbConIds))
    idx_offset=0
    do nb_e=1,me%M2D_Prop%nb_elements
      Connectivity(nb_e+idx_offset)=me%M2D_Prop%elem(nb_e)%nod(1)-1
      Connectivity(nb_e+1+idx_offset)=me%M2D_Prop%elem(nb_e)%nod(2)-1
      Connectivity(nb_e+2+idx_offset)=me%M2D_Prop%elem(nb_e)%nod(3)-1
      idx_offset=idx_offset+2
    enddo
    dims_Connectivity=(/NbConIds/)
    maxdims_Connectivity=dims_Connectivity
    call h5screate_simple_f(1, dims_Connectivity, Connectivity_space_id, error, maxdims_Connectivity)
    call h5dcreate_f(gr_id, "Connectivity", H5T_STD_I64LE, Connectivity_space_id, Connectivity_id, error)
    call h5dwrite_f(Connectivity_id, H5T_STD_I64LE, Connectivity, dims_Connectivity, error)
    call h5sclose_f(Connectivity_space_id, error)
    deallocate(Connectivity)

    ! CellData/Temperature
    call h5gcreate_f(gr_id, "CellData", grA_id, error)
      allocate(Temperature(me%M2D_Prop%nb_elements))
      do nb_e=1,me%M2D_Prop%nb_elements
        Temperature(nb_e)=me%StVar%temp(nb_e)
      enddo
      dims_Temperature=(/me%M2D_Prop%nb_elements/)
      maxdims_Temperature=dims_Temperature
      call h5screate_simple_f(1, dims_Temperature, Temperature_space_id, error, maxdims_Temperature)
      call h5dcreate_f(grA_id, "Temperature", H5T_IEEE_F64LE, Temperature_space_id, Temperature_id, error)
      call h5dwrite_f(Temperature_id, H5T_IEEE_F64LE, Temperature, dims_Temperature, error)
      call h5sclose_f(Temperature_space_id, error)
      deallocate(Temperature)

      allocate(TimeArr(me%M2D_Prop%nb_elements))
      do nb_e=1,me%M2D_Prop%nb_elements
        TimeArr(nb_e)=time
      enddo
      dims_Time=(/me%M2D_Prop%nb_elements/)
      maxdims_Time=dims_Time
      call h5screate_simple_f(1, dims_Time, Time_space_id, error, maxdims_Time)
      call h5dcreate_f(grA_id, "Time", H5T_IEEE_F64LE, Time_space_id, Time_id, error)
      call h5dwrite_f(Time_id, H5T_IEEE_F64LE, TimeArr, dims_Time, error)
      call h5sclose_f(Time_space_id, error)
      deallocate(TimeArr)

    call h5gclose_f(grA_id, error)
    call h5gclose_f(gr_id, error)  
    call h5gclose_f(group_id, error)  

end subroutine writing_HDF5_2D_casing



subroutine additional_data_HDF5_2D_casing(me,NbSteps)
  type(mesh2D_t),     intent(inout) :: me
  integer(i_64), intent(in) :: NbSteps

  integer(hsize_t), dimension(1) :: maxdims_nbsteps
  integer(hsize_t), dimension(1) :: dims_nbsteps
  integer(hid_t) :: nbsteps_id,nbsteps_space_id
  integer :: error

  ! NumberOfSteps
  dims_nbsteps=(/1/)
  maxdims_nbsteps=dims_nbsteps
  call h5screate_simple_f(1, dims_nbsteps, nbsteps_space_id, error, maxdims_nbsteps)
  call h5dcreate_f(me%f_id_writing, "NumberOfSteps", H5T_STD_I64LE, nbsteps_space_id, nbsteps_id, error)
  call h5dwrite_f(nbsteps_id, H5T_STD_I64LE, NbSteps+1, dims_nbsteps, error)
  call h5sclose_f(nbsteps_space_id, error)

end subroutine additional_data_HDF5_2D_casing


end module cmp_mesh2D_hdf5_write_m
