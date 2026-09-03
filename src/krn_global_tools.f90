! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module krn_global_tools_m
    !use, intrinsic :: iso_fortran_env ! This is required only by hdf5 write
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp=>real64, sp=>real32, i_64=>int64, i_8=>int8, &
        stdout=>output_unit, stderr=>error_unit
    use iso_c_binding 
    implicit none

    ! PROJECT GLOBALS -----------------------------------------------------------------------
    logical :: R_Correction !! mechanical equilibrium correction by adding new state variable
    character(:), allocatable :: saved_cwd !! current working directory saved 
    ! PROJECT GLOBALS -----------------------------------------------------------------------

    integer,  parameter :: Nb_VarP = 7     !! Number of primitive variables for Helium model
    integer,  parameter :: Nb_VarC = 4     !! Number of conservative variables for Helium model
    real(dp), parameter :: Pi_value = 3.14159265_dp
    real(dp), parameter :: id_4x4(4,4)     = reshape([1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, &
                                                      0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0],[4,4])
    integer,  parameter :: idx_4x4_col(16) = [1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4] 
    integer,  parameter :: idx_4x4_row(16) = [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4] 
    integer,  parameter :: idx_2x1_col(2)  = [1,1] 
    integer,  parameter :: idx_2x1_row(2)  = [3,4] 
    integer,  parameter :: idx_1x4_col(4)  = [1,2,3,4] 
    integer,  parameter :: idx_1x4_row(4)  = [1,1,1,1] 
    integer,  parameter :: idx_2x4_col(8)  = [1,2,3,4,1,2,3,4]
    integer,  parameter :: idx_2x4_row(8)  = [3,3,3,3,4,4,4,4] 

    type(c_ptr), allocatable :: lib_handles(:)

    interface to_str
        module procedure str_from_int,str_from_real
    end interface

    interface is_nan
        module procedure is_nan_int,is_nan_dbl
    end interface

    type timer_point_t
        character(:), pointer :: description
        integer :: count = 0, tic = 0, sum = 0
    end type

    type timer_t(sz)
        integer, len :: sz
        type(timer_point_t) :: list(sz)
        contains
            procedure :: set => timer_set
            procedure :: print => timer_print
    end type timer_t
    
    type str_ptr
        character(:), pointer :: p => null()
    end type

    interface
        function LoadLibraryA(lpLibFileName) bind(C, name="LoadLibraryA")
            import :: c_ptr, c_char
            type(c_ptr) :: LoadLibraryA
            character(kind=c_char), dimension(*) :: lpLibFileName
        end function

        function GetProcAddress(hModule, lpProcName) bind(C, name="GetProcAddress")
            import :: c_funptr, c_ptr, c_char
            type(c_funptr) :: GetProcAddress
            type(c_ptr), value :: hModule
            character(kind=c_char), dimension(*) :: lpProcName
        end function        
    end interface

    abstract interface    
        function itf_1var(c_pt,var1) bind(C) result(outvar)
          import :: c_ptr, c_double
          type(c_ptr) ,   value :: c_pt
          real(c_double), value :: var1
          real(c_double) :: outvar
        end function itf_1var   

        function itf_2var(c_pt,var1,var2) bind(C) result(outvar)
            import :: c_ptr, c_double
            type(c_ptr),    value :: c_pt
            real(c_double), value :: var1,var2
            real(c_double) :: outvar
        end function itf_2var

        function itf_4var(c_pt,var1,var2,var3,var4) bind(C) result(outvar)
          import :: c_ptr, c_double
          type(c_ptr) ,   value :: c_pt
          real(c_double), value :: var1,var2,var3,var4
          real(c_double) :: outvar
        end function itf_4var        
    
        function itf_5var(c_pt,var1,var2,var3,var4,var5) bind(C) result(outvar)
            import :: c_ptr, c_double
            type(c_ptr) ,   value :: c_pt
            real(c_double), value :: var1,var2,var3,var4,var5
            real(c_double) :: outvar
        end function itf_5var
    
        function itf_8var(c_pt,var1,var2,var3,var4,var5,var6,var7,var8) bind(C) result(outvar)
            import :: c_ptr, c_double
            type(c_ptr) ,   value :: c_pt
            real(c_double), value :: var1,var2,var3,var4,var5,var6,var7,var8
            real(c_double) :: outvar
        end function itf_8var
    end interface   
    
contains

subroutine timer_set(me,to,from,description)
    class(timer_t(*)), intent(inout) :: me
    integer,           intent(in)    :: to
    integer, optional, intent(in)    :: from
    character(*), target, optional, intent(in) :: description !! to be printed on summary
    call system_clock(me%list(to)%tic)
    me%list(to)%count = me%list(to)%count + 1
    if (present(from)) me%list(to)%sum = me%list(to)%sum + &
                            (me%list(to)%tic - me%list(from)%tic)
    if (present(description)) me%list(to)%description=>description
end subroutine

subroutine timer_print(me)
    class(timer_t(*)), intent(inout) :: me
    integer :: i, tic, rate
    character(50) :: n_times
    real(dp) :: time_s
    call system_clock(tic, rate)
    do i=1,size(me%list)
        n_times = ''
        if (me%list(i)%count > 1) n_times = ' '//to_str(me%list(i)%count)//' times'
        if (me%list(i)%sum > 0) then
            time_s = real(me%list(i)%sum)/real(rate)
            print '(f9.2,a)', time_s, ' s '//me%list(i)%description//trim(n_times)
        endif
    enddo
end subroutine

pure function str_from_int(value_in) result(str_out)
    integer, intent(in)       :: value_in
    character(:), allocatable :: str_out
    character(20) :: tmp
    write(tmp,'(i0)') value_in
    str_out = trim(tmp)
end function

pure function str_from_real(value_in) result(str_out)
    real(dp), intent(in)      :: value_in
    character(:), allocatable :: str_out
    character(30) :: tmp
    write(tmp,*) value_in
    str_out = trim(adjustl(tmp))
end function

elemental function to_dbl(str_in) result(val_out)
    character(*), intent(in) :: str_in
    real(dp) :: val_out
    integer :: error
    val_out = ieee_value(val_out, ieee_quiet_nan)
    read(str_in,*,ioStat=error) val_out
end function

elemental function to_int(str_in) result(val_out)
    character(*), intent(in) :: str_in
    integer :: val_out
    integer :: error
    val_out = -huge(val_out)
    read(str_in,*,ioStat=error) val_out
end function

function is_nan_dbl(val_in) result(val_out)
    real(dp), intent(in) :: val_in
    logical :: val_out
    val_out = ieee_is_nan(val_in)
end function

function is_nan_int(val_in) result(val_out)
    integer, intent(in) :: val_in
    logical :: val_out
    if(val_in == -huge(val_in)) then
        val_out = .true.
    else
        val_out = .false.
    endif
end function

subroutine check_status(status,error_message)
    !! If status /= 0 abort printing error_message and status
    integer, intent(in) :: status
    character(*), intent(in) :: error_message
    if (status /= 0) then
        write(stderr,*) error_message//to_str(status)
        stop 1
    end if
end subroutine

function reims_config_file() result(path)
    use hdf5
    use iso_c_binding
    character(:),allocatable :: path
    character(1024) :: buff
    integer(c_int) :: err
    integer :: i
    type(c_ptr) :: ptr

    ! C functions interface declaration
    interface
        function c_chDir(new_path) bind(C ,name="chdir")
            import :: c_int
            character :: new_path(*)
            integer(c_int) :: c_chDir 
        end function
        function c_getCwd(buf, ln) bind(C, name="getcwd")
            import :: c_char, c_int, c_ptr
            type(c_ptr) :: c_getCwd
            character(c_char) :: buf(*)
            integer(c_int), value :: ln
        end function c_getCwd
    end interface

    call h5open_f(i)

    ptr = c_getCwd(buff,len(buff))
    i = scan(buff,char(0))
    saved_cwd = buff(:i)

    if (nArgs() == 0) error stop 'Please provide input yaml file name in the command line'
    call getArg(1, buff)
    i = scan(buff,'/\',.true.)
    if(i > 0) then
        err = c_chDir(buff(:i-1)//char(0))
    endif
    path = trim(buff(i+1:))
end function reims_config_file

function convert2int_port(str) result(convert)
    character(len=*), intent(in) :: str
    integer :: convert

    select case (str)
        case ('in')
            convert = 1
        case ('out')
            convert = 2
    end select
end function convert2int_port

integer function count_chars(str, char)
character(len=*), intent(in) :: str
character(len=1), intent(in) :: char
integer :: i

count_chars = 0
do i = 1, len_trim(str)
    if (str(i:i) == char) then
        count_chars = count_chars + 1
    end if
end do
end function count_chars

subroutine check_if_conversion_required(list,name,conversion_required,max_value)
    character(len=*), intent(in) :: list,name
    logical, intent(out) :: conversion_required
    integer, intent(out) :: max_value

    integer :: name_pos, colon_pos, end_pos, num_commas, str_len, i, num_count, pos, start, length
    character(:), allocatable :: value_str
    integer, allocatable :: aa(:)
    
    name_pos = index(list, trim(name))
    colon_pos = index(list(name_pos:), ':') + name_pos - 1
    end_pos = index(list(colon_pos+1:), ' ')
    if (end_pos == 0) then
        end_pos = len(list)+1
    else
        end_pos = colon_pos + end_pos
    end if
    value_str = list(colon_pos+1:end_pos-1)
    conversion_required = (trim(value_str) /= '0')

    max_value=0
    if(conversion_required) then
        start = 1
        pos = 0
        length = len_trim(value_str)
        num_count = count_chars(value_str, ',') + 1
        allocate(aa(num_count))
    
        do i = 1, num_count
            pos = scan(value_str(start:length), ',')
            if (pos == 0) then
                read(value_str(start:length), *) aa(i)
            else
                read(value_str(start:start + pos - 2), *) aa(i)
                start = start + pos
            end if
        end do
        max_value=maxval(aa)
    endif

end subroutine check_if_conversion_required    

subroutine array_conv(list,name,conv)
    character(len=*), intent(in) :: list,name
    integer, intent(out) :: conv(:)

    integer :: name_pos, colon_pos, end_pos, i, idx, position, num, next_comma
    character(:), allocatable :: value_str

    name_pos = index(list, trim(name))
    colon_pos = index(list(name_pos:), ':') + name_pos - 1
    end_pos = index(list(colon_pos+1:), ' ')
    if (end_pos == 0) then
        end_pos = len(list)+1
    else
        end_pos = colon_pos + end_pos
    end if
    value_str = list(colon_pos+1:end_pos-1)   

    conv = 0
    idx = 1
    position = 1
    do
        read(value_str(position:), '(I5)', iostat=i) num
        if (i /= 0) exit
        if (num <= size(conv)) conv(num) = idx
        idx = idx + 1
        next_comma = index(value_str(position:), ',')
        if (next_comma == 0) exit
        position = next_comma + position
        if (position == 1) exit
    end do

end subroutine array_conv    

end module krn_global_tools_m

