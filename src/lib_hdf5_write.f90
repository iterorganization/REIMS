! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module lib_hdf_write_m
    !! Write library for writing hdf5 files use by components 
    use krn_global_tools_m
    use lib_input_m, only: input_t
    use iso_c_binding
    use hdf5
    use h5ds
    use h5lt
    !use h5lt_const

    implicit none
    private
    public hdf5_t, hdf_desc_t, sim_delta

    type hdf_desc_t
        character(:), allocatable :: name      !! name that will appear in HDF5 file
        real(dp),     allocatable :: node_x(:) !! spacial axes (scale) of 1d element
        real(dp),         pointer :: data(:,:) !! data to write into hdf5 data table
    end type hdf_desc_t

    type hdf_desc_p_t
        type(hdf_desc_t), pointer :: p 
    end type hdf_desc_p_t

    type tbl_dsc_t
        character(:), allocatable :: name
        character(:), allocatable :: vars(:)
        type(hdf_desc_p_t), allocatable :: comps(:)
    end type tbl_dsc_t

    type growing_h5_tbl_t
        integer(hid_t)   :: dset        !! Dataset is kept open
        integer          :: rank
        integer(hsize_t) :: dims(8) = 0 !! Value 3 should be enough
        real(dp), allocatable  :: data(:,:)
    end type growing_h5_tbl_t

    type hdf5_t
        logical  :: write_results       !! from cfg: write result to file active
        real(dp) :: min_save_time       !! from cfg: minimum time before saving
        real(dp) :: last_saved_time = 0 !! updated when data is saved
        integer  :: steps_to_revert = 0 !! nb of HDF5 entries to revert after rollback (persistent across rejected steps)
        logical  :: write_results2D     !! from cfg: write 2D results to specific files
        real(dp) :: time_btw_2D_writes  !! from cfg: physical time between 2D writes

        integer(hid_t)                      :: file_id
        type(tbl_dsc_t),        allocatable :: tbl_dsc(:)
        type(growing_h5_tbl_t), allocatable :: tables(:)
        type(growing_h5_tbl_t)              :: time, global !! growing tables
    contains
        procedure :: add_table    => hdf5_add_table
        procedure :: add_to_table => hdf5_add_to_table
        procedure :: prepare_data => hdf5_prepare_data
        procedure :: init         => hdf5_init
        procedure :: close        => hdf5_close
        procedure :: write        => hdf5_write
    end type hdf5_t

    interface
        integer function H5Fstart_swmr_write(f_id) bind(C,name='H5Fstart_swmr_write')
            import :: hid_t
            integer(hid_t), value :: f_id
        end function
    end interface

    interface
        integer function H5Dflush(dset_id) bind(C,name='H5Dflush')
            import :: hid_t
            integer(hid_t), intent(in), value :: dset_id
        end function H5Dflush
    end interface

    real(dp), parameter :: sim_delta = 1.0e-10_dp !! Minimum value for criteria calculation     

contains

subroutine hdf5_add_table(me, table, variables)
    class(hdf5_t),  intent(inout) :: me
    character(*),   intent(in)    :: table, variables(:)
    type(hdf_desc_p_t) :: tmp(0)
    if (allocated(me%tbl_dsc)) then
        me%tbl_dsc = [me%tbl_dsc, tbl_dsc_t(table,variables,tmp)]
    else
        me%tbl_dsc = [tbl_dsc_t(table,variables,tmp)]
    endif
end subroutine hdf5_add_table

subroutine hdf5_add_to_table(me, desc, tbl_name)
    class(hdf5_t),            intent(inout) :: me
    type(hdf_desc_t), target, intent(inout) :: desc
    character(*),             intent(in)    :: tbl_name

    integer :: i
    type(hdf_desc_p_t) :: tmp

    do i = 1, size(me%tbl_dsc)
        if (me%tbl_dsc(i)%name==tbl_name) exit
    enddo
    tmp%p => desc 
    me%tbl_dsc(i)%comps = [me%tbl_dsc(i)%comps, tmp]
end subroutine hdf5_add_to_table

subroutine hdf5_prepare_data(me)
    class(hdf5_t), target, intent(inout) :: me
    integer :: i,j,total,nodes

    allocate(me%tables(size(me%tbl_dsc)))
    do i = 1, size(me%tables); associate(tbl => me%tables(i), dsc => me%tbl_dsc(i))
        total = 0
        do j = 1, size(dsc%comps)
            total = total + size(dsc%comps(j)%p%node_x)
        enddo
        allocate(tbl%data(total,size(dsc%vars)))
        tbl%data = ieee_value(0.0_dp, ieee_quiet_nan)
        total = 0
        do j = 1, size(me%tbl_dsc(i)%comps)
            nodes = size(dsc%comps(j)%p%node_x)
            dsc%comps(j)%p%data => tbl%data(total+1:total+nodes,1:size(dsc%vars))
            total = total + nodes
        enddo
    end associate; enddo
end subroutine hdf5_prepare_data

subroutine hdf5_init(me,cfg)
    class(hdf5_t),  intent(inout) :: me
    class(input_t), intent(in)    :: cfg
    
    integer(hid_t) :: group_id, prop_id, var_id
    real(dp)       :: data_global(2)
    integer        :: i, err

    me%write_results = cfg%bin('active',.true.)
    me%write_results2D = cfg%bin('active2D',.true.)    

    if(.not.me%write_results) return
    print*,'Writing of results ACTIVE'

    me%min_save_time = cfg%dbl('minimum_time',0.0_dp)
    me%time_btw_2D_writes = cfg%dbl('time_between_2D_writes',3600.0_dp)
    
    ! Create hdf file and 'reims' group
    call H5Pcreate_f(H5P_FILE_ACCESS_F, prop_id, err)
    call H5Pset_libver_bounds_f(prop_id, H5F_LIBVER_LATEST_F, H5F_LIBVER_LATEST_F, err)
    call H5Fcreate_f(cfg%str('file'), H5F_ACC_TRUNC_F, me%file_id, err, H5P_DEFAULT_F, prop_id)
    call h5pclose_f(prop_id, err)
    call h5gcreate_f(me%file_id, 'reims', group_id, err)

    call h5ltset_attribute_string_f(me%file_id,'reims','type','reims'//char(0),err)
    call tbl_h5_create_dbl(me%time,group_id,'t',0.)

    data_global = 0.0_dp
    call tbl_h5_create_dbl(me%global,group_id,'global',data_global)
    call write_str_arr(group_id,'global_var',['dt','dt2'],var_id)
    call attach_scale(var_id,       me%global%dset, 1, 'var', 'var')
    call attach_scale(me%time%dset, me%global%dset, 2, 't',   't')
    call h5dclose_f(var_id,err)

    do i = 1, size(me%tbl_dsc)
        call hdf5_data_create(me%tables(i), group_id, me%time%dset, &
            me%tbl_dsc(i)%name, me%tbl_dsc(i)%vars, me%tbl_dsc(i)%comps)
    enddo
    call h5gclose_f(group_id, err)
    err = H5Fstart_swmr_write(me%file_id)

end subroutine hdf5_init

subroutine hdf5_data_create(tbl, group_id, time_id, name, vars, comps)
    type(growing_h5_tbl_t), intent(inout) :: tbl
    integer(hid_t),            intent(in) :: group_id, time_id
    character(*),              intent(in) :: name,vars(:)
    type(hdf_desc_p_t),       intent(inout) :: comps(:)

    type(hdf_desc_t), allocatable :: components(:)
    integer(hid_t)                :: x_id, var_id, name_id, cmp_id
    integer(hsize_t)              :: dims(8)
    real(dp),         allocatable :: x(:)
    integer(hsize_t), allocatable :: size_attr(:), size_cmp(:), size_name(:), size_var(:)
    integer,          allocatable :: int_attr(:)
    character(:),     allocatable :: str_attr(:), str_cmp(:), str_name(:), tmp_str
    character(9)                  :: num_str
    integer                       :: i,j,cell,cell_tot,err,str_trim

    cell_tot = 0
    allocate(components(size(comps)))
    do i = 1, size(comps)
        components(i) = comps(i)%p
        cell_tot = cell_tot + size(components(i)%node_x)
    enddo
    if (cell_tot == 0) return
    j = 0 ! finding longest component name and allocate required arrays
    do i = 1, size(components)
        j = max(j,len(components(i)%name))
    enddo
    allocate(character(j)   :: str_attr(size(components)))
    allocate(                  size_attr(size(components)))
    allocate(                  int_attr(size(components)))
    allocate(                  x(cell_tot))
    allocate(character(j)   :: str_cmp(cell_tot))
    allocate(                  size_cmp(cell_tot))
    allocate(character(j+9) :: str_name(cell_tot))
    allocate(                  size_name(cell_tot))

    size_var = len_trim(vars)
    cell = 0 ! fill arrays: loop over component and later on node
    do i = 1, size(components)
        j = size(components(i)%node_x)
        int_attr(i) = j
        x(cell+1:cell+j) = components(i)%node_x
        str_attr(i) = components(i)%name
        size_attr(i) = len(components(i)%name)
        write(num_str,'(i0)'),j
        str_trim = len_trim(num_str)
        str_cmp(cell+1:cell+j) = components(i)%name
        size_cmp(cell+1:cell+j) = len(components(i)%name)
        size_name(cell+1:cell+j) = len(components(i)%name) + str_trim + 1
        tmp_str = '(i'//to_str(str_trim+1)//'.'//to_str(str_trim)//')'
        do j = 1, size(components(i)%node_x)
            write(num_str,tmp_str),j
            str_name(cell+j) = components(i)%name//num_str
        enddo
        cell = cell + j - 1
    enddo

    call tbl_h5_create_dbl(tbl,group_id,name,tbl%data)
    call write_vstr1d(group_id,name//'_cmp',str_cmp,size_cmp,cmp_id)
    call write_vstr1d(group_id,name//'_name',str_name,size_name,name_id)
    dims = 0
    dims(1) = cell
    call h5ltmake_dataset_f(group_id,name//'_x',1,dims,H5T_NATIVE_DOUBLE,x,err)
    call h5dopen_f(group_id,name//'_x',x_id,err)
    call write_vstr1d(group_id,name//'_var',vars,size_var,var_id)

    call attach_scale(time_id, tbl%dset, 3, 't', 't')
    call attach_scale(var_id,  tbl%dset, 1, 'var','var')
    call attach_scale(x_id,    tbl%dset, 2, 'x')
    call attach_scale(cmp_id,  tbl%dset, 2, 'cmp')
    call attach_scale(name_id, tbl%dset, 2, 'name','name')

    ! Writing attributes
    call write_atribute_int(tbl%dset,'nodes',int_attr)
    call write_attribute_vstr1d(tbl%dset,'names',str_attr,size_attr)

    call h5dclose_f(var_id,err)
    call h5dclose_f(name_id,err)
    call h5dclose_f(x_id,err)
    call h5dclose_f(cmp_id,err)
end subroutine hdf5_data_create


subroutine hdf5_write(me,t,dt,t_final)
    class(hdf5_t), intent(inout) :: me
    real(dp), intent(in) :: t,dt,t_final

    real(dp) :: data_global(2,1), NextWriteTime
    integer :: i
    logical :: LastTimeStep

    if(.not.me%write_results) return
    NextWriteTime = me%last_saved_time + me%min_save_time 
    LastTimeStep = t + sim_delta >= t_final
    if(t < NextWriteTime .and. .not. LastTimeStep) then
        if(me%steps_to_revert > 0) then
            call tbl_h5_revert_dbl(me%time,   me%steps_to_revert)
            call tbl_h5_revert_dbl(me%global, me%steps_to_revert)
            do i = 1, size(me%tables)
                if(me%tables(i)%dims(2) /= 0) &
                    call tbl_h5_revert_dbl(me%tables(i), me%steps_to_revert)
            enddo
            me%steps_to_revert = 0
            me%last_saved_time = t
        endif
        return
    endif
    me%last_saved_time = t
    
    do i = 1, size(me%tables)
        if (me%tables(i)%dims(2) /= 0) &
            call tbl_h5_append_dbl(me%tables(i), me%tables(i)%data)
    enddo

    data_global(:,1) = [dt,dt]
    call tbl_h5_append_dbl(me%global, data_global)

    call tbl_h5_append_dbl(me%time,t)
end subroutine hdf5_write    

subroutine hdf5_close(me)
    class(hdf5_t), intent(inout) :: me
    integer :: i, error
    if (.not.me%write_results) return
    do i = 1, size(me%tables)
        if (me%tables(i)%dims(2)/=0) call h5dclose_f(me%tables(i)%dset,error)
    enddo
    call h5dclose_f(me%time%dset,   error)
    call h5dclose_f(me%global%dset, error)
    call h5fclose_f(me%file_id,     error)
    call h5close_f(error)
end subroutine hdf5_close 

subroutine write_atribute_int(location,attribute_name,attribute_value)
    integer(hid_t), intent(in) :: location            !! hdf5 file or group id
    character(*),   intent(in) :: attribute_name
    integer,        intent(in) :: attribute_value(..) !! from scalar up to rank 3

    integer :: error, rank_in
    integer(hid_t)   :: dspc_id, attr_id
    integer(hsize_t), allocatable :: dims(:)

    rank_in = rank(attribute_value)
    dims = shape(attribute_value)
    call h5screate_simple_f(rank_in, dims, dspc_id, error, dims)
    call h5acreate_f(location, attribute_name, H5T_NATIVE_INTEGER, dspc_id, attr_id, error)
    select rank(attribute_value)
        rank(0)
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attribute_value, dims, error)
        rank(1)
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attribute_value, dims, error)
        rank(2)
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attribute_value, dims, error)
        rank(3)
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attribute_value, dims, error)
        rank default
            error stop 'Integer attribute rank can be from 0 to 3'
    end select
    call h5aclose_f(attr_id, error)
    call h5sclose_f(dspc_id, error)
    deallocate(dims)
end subroutine write_atribute_int

subroutine tbl_h5_revert_dbl(me, n_steps_back)
    class(growing_h5_tbl_t), intent(inout) :: me
    integer, intent(in) :: n_steps_back
    integer :: err
    me%dims(1) = max(1, me%dims(1)-n_steps_back)    ! shrink the dimension counter
    call h5dset_extent_f(me%dset, me%dims, err)     ! truncate the HDF5 dataset
    err = H5Dflush(me%dset)
end subroutine tbl_h5_revert_dbl

subroutine tbl_h5_create_dbl(me,location,name,datas)
    class(growing_h5_tbl_t), intent(inout) :: me
    integer(hid_t),     intent(in) :: location
    character(*),       intent(in) :: name
    real(dp),           intent(in) :: datas(..)

    integer :: err
    integer(hid_t)   :: prop, dspc
    integer(hsize_t) :: max_dims(8), chunk_dims(8)

    me%rank = rank(datas) + 1
    me%dims(1)         = 1
    me%dims(2:me%rank) = shape(datas)
    max_dims           = me%dims
    chunk_dims         = me%dims
    max_dims(1)        = H5S_UNLIMITED_F
    chunk_dims(1)      = 1 ! was 128, setting to 1 to avoid calls to fseek()

    call h5screate_simple_f(me%rank, me%dims, dspc, err, max_dims) 
    call h5pcreate_f(H5P_DATASET_CREATE_F, prop, err)
    call h5pset_chunk_f(prop, me%rank, chunk_dims, err)
    call h5dcreate_f(location, name, H5T_NATIVE_DOUBLE, dspc, me%dset, err, prop)
    call h5pclose_f(prop, err)
    call h5sclose_f(dspc, err)
    select rank(datas)
        rank(0)
            call h5dwrite_f(me%dset, H5T_NATIVE_DOUBLE, datas, me%dims, err)
        rank(1)
            call h5dwrite_f(me%dset, H5T_NATIVE_DOUBLE, datas, me%dims, err)
        rank(2)
            call h5dwrite_f(me%dset, H5T_NATIVE_DOUBLE, datas, me%dims, err)
        rank default
            error stop "growing array can be only rank from 0 to 2"
    end select
end subroutine tbl_h5_create_dbl

subroutine tbl_h5_append_dbl(me,datas)
    class(growing_h5_tbl_t), intent(inout) :: me
    real(dp),target, intent(in)    :: datas(..)

    integer :: err
    integer(hid_t)    :: mem_sp, dset_sp
    integer(hsize_t)  :: offset(8), mem_dim(8)
    real(dp), pointer :: data_2d(:,:), data_3d(:,:,:)
    real(dp)          :: data_1d(1)

    offset     = 0
    offset(1)  = me%dims(1)
    mem_dim    = me%dims
    mem_dim(1) = 1     
    me%dims(1) = 1 + me%dims(1)

    call h5dset_extent_f(me%dset, me%dims, err)             ! increase size of the dataset
    call h5screate_simple_f (me%rank, mem_dim, mem_sp, err) ! create memory space
    call h5dget_space_f(me%dset, dset_sp, err)              ! get dataset space
    call h5sselect_hyperslab_f(dset_sp, H5S_SELECT_SET_F, offset, mem_dim, err)
    select rank(datas)
        rank(0)
            data_1d = [datas]
            call H5dwrite_f(me%dset, H5T_NATIVE_DOUBLE, data_1d, mem_dim, &
                err, mem_sp, dset_sp)
        rank(1)
            data_2d(1:1,1:size(datas)) => datas
            call H5dwrite_f(me%dset, H5T_NATIVE_DOUBLE, data_2d, mem_dim, &
                err, mem_sp, dset_sp)
        rank(2)
            data_3d(1:1,1:size(datas,1),1:size(datas,2)) => datas
            call H5dwrite_f(me%dset, H5T_NATIVE_DOUBLE, data_3d, mem_dim, &
                err, mem_sp, dset_sp)
        rank default
            error stop "growing array can be only rank from 0 to 2"
    end select
    call h5sclose_f(dset_sp, err)
    call h5sclose_f(mem_sp, err)
    err = H5Dflush(me%dset)
end subroutine tbl_h5_append_dbl

subroutine write_vstr1d(location,dset_name,value,length,id_out)
    !! Write 1d variable string array to hdf5 file
    integer(hid_t),  intent(in) :: location  !! hdf5 file or group id
    character(*),    intent(in) :: dset_name
    character(*),    intent(in) :: value(:)  !! from scalar up to rank 3
    integer(size_t), intent(in) :: length(:)
    integer(hid_t), optional, intent(out) :: id_out

    integer(hid_t)   :: dspc_id, dset_id, tstr_id
    integer(hsize_t) :: dims2(2),dims1(1)
    integer          :: err

    dims2(1) = len(value)
    dims2(2) = size(value)
    dims1(1) = dims2(2)
    call h5tcopy_f(H5T_C_S1, tstr_id, err)
    call h5tset_size_f(tstr_id, -1_hsize_t, err)
    call h5screate_simple_f(1, dims1, dspc_id, err)
    call h5dcreate_f(location, dset_name, tstr_id, dspc_id, dset_id, err)
    call h5dwrite_vl_f(dset_id, tstr_id, value, dims2, length, err, dspc_id)
    call h5tclose_f(tstr_id,err)
    call h5sclose_f(dspc_id,err)
    if (present(id_out)) then
        id_out = dset_id
    else
        call h5dclose_f(dset_id,err)
    endif
end subroutine write_vstr1d

subroutine write_attribute_vstr1d(location,attribute_name,value_in,length)
    !! BROKEN: (not variable) Write 1d atribute variable string array to hdf5 file
    integer(hid_t),  intent(in) :: location  !! hdf5 file or group id
    character(*),    intent(in) :: attribute_name
    character(*),    intent(in) :: value_in(:)  !! from scalar up to rank 3
    integer(size_t), intent(in) :: length(:) !! BROKEN: unused right now

    integer(hid_t)   :: dspc_id, dset_id, tstr_id
    integer(hsize_t) :: dims2(2),dims1(1)
    integer          :: i,err
    character(:),allocatable :: val(:)


    dims2(1) = len(value_in) + 1
    dims2(2) = size(value_in)
    dims1(1) = dims2(2)
    allocate(character(dims2(1)) :: val(dims2(2)))
    do i = 1, size(val)
        val(i) = trim(value_in(i))//char(0)
    enddo
    call h5tcopy_f(H5T_C_S1, tstr_id, err)
    call h5tset_size_f(tstr_id, dims2(1), err)
    call h5screate_simple_f(1, dims1, dspc_id, err)
    call h5acreate_f(location, attribute_name, tstr_id, dspc_id, dset_id, err)
    call h5awrite_f(dset_id, tstr_id, val, dims2, err)
    call h5tclose_f(tstr_id,err)
    call h5sclose_f(dspc_id,err)
    call h5aclose_f(dset_id,err)
end subroutine write_attribute_vstr1d

subroutine write_str_arr(location,dset_name,value,id_out)
    integer(hid_t),  intent(in) :: location  !! hdf5 file or group id
    character(*),    intent(in) :: dset_name
    character(*),    intent(in) :: value(..)  !! from rank 1 up to rank 3
    integer(hid_t), optional, intent(out) :: id_out

    integer(hid_t)   :: dspc_id, dset_id, tstr_id
    integer(hsize_t) :: dims(7),str_len
    integer          :: err, rank_in

    rank_in = rank(value)
    dims    = 0
    dims(1:rank_in) = shape(value)
    select rank(value)
        rank(1) 
            str_len = len(value(1))
        rank(2)
            str_len = len(value(1,1))
        rank(3)
            str_len = len(value(1,1,1))
        rank default
            print *, "array can be only rank from 1 to 3"
    end select
    call h5tcopy_f(H5T_FORTRAN_S1, tstr_id, err)
    call h5tset_cset_f(tstr_id, H5T_CSET_ASCII_F, err)
    call h5tset_size_f(tstr_id, str_len, err)
    call h5screate_simple_f(rank_in, dims, dspc_id, err)
    call h5dcreate_f(location, dset_name, tstr_id, dspc_id, dset_id, err)
    select rank(value)
        rank(1) 
            call h5dwrite_f(dset_id, tstr_id, value, dims, err, dspc_id)
        rank(2)
            call h5dwrite_f(dset_id, tstr_id, value, dims, err, dspc_id)
        rank(3)
            call h5dwrite_f(dset_id, tstr_id, value, dims, err, dspc_id)
    end select
    call h5tclose_f(tstr_id,err)
    call h5sclose_f(dspc_id,err)
    if (present(id_out)) then
        id_out = dset_id
    else
        call h5dclose_f(dset_id,err)
    endif
end subroutine write_str_arr

subroutine attach_scale(scale_id, data_id, axies, name, label)
    integer(hid_t),         intent(in) :: scale_id
    integer(hid_t),         intent(in) :: data_id
    integer,                intent(in) :: axies
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: label

    integer :: err

    if(present(label)) call h5dsset_label_f(data_id, axies, label, err)
    call h5dsset_scale_f(scale_id, err, name)
    call h5dsattach_scale_f(data_id, scale_id, axies, err)
end subroutine attach_scale


end module lib_hdf_write_m