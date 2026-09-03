! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module krn_simulation_m
    !! Step management and interpolation in time
    use lib_input_m
    use krn_global_tools_m
    use lib_hdf_write_m
    use fortran_yaml_c, only: YamlFile,  type_key_value_pair, &
                        type_node, type_dictionary, type_error, &
                        type_list, type_list_item,  type_scalar
    implicit none
    private 
    public simulation_t, signal_t, h5, set_error, sim_error, sim_error_msg, MAX_SIM_ERRORS

    abstract interface
      subroutine itf_signal(c_pt, time, signal) bind(C)
        import :: c_ptr, c_double, c_int
        type(c_ptr),   value :: c_pt
        real(c_double), value :: time
        real(c_double)        :: signal(*)
      end subroutine itf_signal   
          
      function init_signal_ext_if(fct_ptr_to_ext,cfg_len,cfg_str,c_pt,sign_ptr,n) bind(C)
        import :: c_int, c_ptr, c_funptr, dp
        type(c_funptr), value :: fct_ptr_to_ext
        integer(c_int), value :: cfg_len
        character(*)    :: cfg_str 
        type(c_ptr)     :: c_pt      
        type(c_funptr)  :: sign_ptr
        integer(c_int), value :: n
        logical         :: init_signal_ext_if
      end function   
    end interface         

    type :: signal_t
        real(dp), allocatable :: time(:)      !! time events (points of interpolation)
        real(dp), allocatable :: offset(:,:)  !! offset values (points of interpolation)
        real(dp), allocatable :: gain(:,:)    !! gain values for 0th order simply 1
        real(dp), allocatable :: x(:)         !! spatial coordinates
        integer :: i = 1     !! last index (starting point for iteration)
        type(simulation_t), pointer :: sim
        type(c_funptr) :: init_signal_ext_ptr = C_NULL_FUNPTR   
        character(:), allocatable :: cfg_str
        type(c_funptr) :: sign_ptr = C_NULL_FUNPTR
        type(c_ptr) :: val_cfg = C_NULL_PTR                
        procedure(itf_signal), pointer, nopass :: val_fpt => null()    
        procedure(init_signal_ext_if), pointer, nopass :: init_signal_ext  => null()             
    contains
        procedure :: init => signal_init
        procedure :: v1d  => signal_v1d
        procedure :: v0d  => signal_v0d
    end type signal_t

    type simulation_t
        ! protected
        real(dp) :: t = 0    !! current simulation time
        real(dp) :: dt       !! current time step
        logical  :: explicit = .true.  !! true for explicit, false for implicit scheme
        logical  :: rejected = .false. !! current step is rejected
        logical  :: rollback = .false.      !! rollback in progress, skip state advance
        logical  :: in_recovery = .false.   !! error recovery explicit steps in progress
        integer  :: recovery_steps_left = 0 !! countdown for recovery explicit steps
        integer  :: consecutive_rollbacks = 0 !! number of consecutive rollbacks without a successful step
        real(dp) :: dtPrev1  !! previous time step
        real(dp) :: dtPrev2  !! time step 2 steps before
        real(dp) :: dtPrev3  !! time step 3 steps before
        real(dp) :: dtPrev4  !! time step 4 steps before
        real(dp) :: dtPrev5  !! time step 5 steps before

        ! private
        real(dp) :: dt_exp !! explicit time step (save for calc)
        type(signal_t) :: impl_tol !! from input file implicit tolerance
        real(dP) :: expl_tol !! from input file explicit tolerance
        real(dP) :: kValue   !! from input file step management filter gain
        real(dP) :: max_step !! from input file maximum step that solver can take
        real(dP) :: t_final  !! from input file time to simulate
        real(dp) :: c0       !! step controller memory
        integer  :: step = 0 !! all steps including rejected
        integer :: expl_steps = 10 !! Minimum number of explicit steps left
        real(dp), allocatable :: events_time(:) !! list of time events forcing step
        logical,  allocatable :: events_expl(:) !! list of time events forcing explicit
        real(dp) :: scheme_diff_fact !! factor to be applied to the acoustic speed estimate
    contains
        procedure :: init => sim_init
        procedure :: update_exp_dt => sim_update_exp_dt
        procedure :: step_management => sim_step_management
    end type simulation_t

    type(hdf5_t) :: h5

    integer                        :: sim_error = 0
    integer, parameter             :: MAX_SIM_ERRORS = 10
    character(512)                 :: sim_error_msg(MAX_SIM_ERRORS)

    type, bind(C) :: comp_desc_c_t
        type(c_ptr)    :: name_ptr
        integer(c_int) :: name_len

        type(c_ptr)    :: node_x_ptr
        integer(c_int) :: n_node     

        type(c_ptr)    :: data_ptr
        integer(c_int) :: n1        
        integer(c_int) :: n2        
    end type comp_desc_c_t

    type, bind(C) :: tbl_dsc_c_t
        type(c_ptr)    :: name_ptr
        integer(c_int) :: name_len

        type(c_ptr)    :: vars_ptr
        integer(c_int) :: n_vars    
        integer(c_int) :: var_len   

        type(c_ptr)    :: comps_ptr     ! -> comp_desc_c_t(:)
        integer(c_int) :: n_comps
    end type tbl_dsc_c_t

    type, bind(C) :: h5_export_c_t
        integer(c_int) :: nb_tables
        type(c_ptr)    :: tables_ptr    ! -> tbl_dsc_c_t(:)
    end type h5_export_c_t

    type(h5_export_c_t), save, target :: H
    type(tbl_dsc_c_t),   allocatable, target, save :: Tbls_c(:)
    type(comp_desc_c_t), allocatable, target, save :: Comps_c(:)

    character(kind=c_char), allocatable, target, save :: vars_buf(:)
    integer(c_int),         allocatable,         save :: vars_offset(:)
    integer(c_int),         allocatable,         save :: n_vars_tbl(:)
    integer(c_int),         allocatable,         save :: var_len_tbl(:)      

contains

subroutine sim_init(me,input)
    class(simulation_t), intent(inout) :: me
    type(input_t), intent(in) :: input
    me%t_final   = input%dbl('simulation_end')
    R_Correction = input%bin('R_correction',.true.)

    call me%impl_tol%init(me, input, 'implicit_tolerance')
    !me%impl_tol = input%dbl('implicit_tolerance',5.0e-4_dp)
    me%expl_tol = input%dbl('explicit_tolerance',1.0e-5_dp)
    me%kValue = input%dbl('step_controller_gain',4.0_dp)
    me%max_step  = input%dbl('max_time_step',50.0_dp)
    me%events_expl = [.false., .false.]
    me%events_time = [me%t_final, huge(0._dp)]
    me%scheme_diff_fact = input%dbl('scheme_diffusion_factor',1.0_dp)
end subroutine sim_init

subroutine sim_update_exp_dt(me,wave_time)
    class(simulation_t), intent(inout) :: me
    real(dp),            intent(in)    :: wave_time
    ! value `dt_exp` is function of `wave_time` and overwritten
    ! for the beginning of the simulation (to have soft start) 
    me%dt_exp = 0.8_dp * wave_time
    if(me%expl_steps > 5) me%dt_exp = 0.5_dp * wave_time
    if(me%expl_steps > 8) me%dt_exp = 0.1_dp * wave_time
    if(me%explicit) me%dt = me%dt_exp
end subroutine sim_update_exp_dt

subroutine sim_step_management(me,err)
    class(simulation_t), intent(inout) :: me
    real(dp), intent(in) :: err
    real(dp) :: c1, ratio, new_dt
    integer  :: i
    me%step = me%step + 1

    if(me%explicit) then
        ! compute pressure criterium to switch from explicit to implicit
        me%expl_steps = me%expl_steps - 1
        me%explicit = err >= me%expl_tol .or. me%expl_steps > 0
        if (me%in_recovery) then
            me%recovery_steps_left = me%recovery_steps_left - 1
            if (me%recovery_steps_left == 0) me%in_recovery = .false.
        end if
        me%rejected = .false.
        new_dt = me%dt_exp
    else
        ! implicit step calculation and test for step rejection
        c1 = me%impl_tol%v0d() / max(err, sim_delta)
        ratio = filter(abs(me%c0 - 1.0_dp) < 1.0e-8_dp, c1, me%c0, me%kValue)
        me%c0 = c1
        me%rejected = ratio < 0.8_dp
        if (.not. me%rejected) me%consecutive_rollbacks = 0
        new_dt = ratio * me%dt
    endif

    if(me%rejected) then
        me%c0 = 1.0_dp
    else
        me%t = me%t + me%dt
        me%dtPrev5 = me%dtPrev4
        me%dtPrev4 = me%dtPrev3
        me%dtPrev3 = me%dtPrev2
        me%dtPrev2 = me%dtPrev1
        me%dtPrev1 = me%dt
    endif
    
    ! search for correct time event indicated by index: `i` shows next event
    i = 1; do while (me%events_time(i) < me%t + sim_delta); i = i + 1; enddo

    ! `dt` calculation by reducing step due to max_step or event
    me%dt = min(new_dt, me%max_step, me%events_time(i) - me%t)

    ! switch to explicit due to too small time step
    if(.not.me%explicit) then
        me%explicit = me%dt < me%dt_exp
        if(me%explicit) me%expl_steps = 3
    endif

    ! switch to explicit if current event requires it
    i = max(1,i-1) ! `i` shows current event
    if(abs(me%t - me%events_time(i)) < sim_delta .and. me%events_expl(i)) then
        me%expl_steps = max(me%expl_steps,3)
        me%explicit = .true. 
    endif

    if(me%explicit) then
        me%c0 = 1.0_dp
        me%dt = me%dt_exp
    endif

    call sim_print(me)
end subroutine sim_step_management

subroutine set_error(msg)
    character(*), intent(in) :: msg
    integer :: slot
    !$omp atomic capture
       slot = sim_error
       sim_error = sim_error + 1
    !$omp end atomic
    if (slot < MAX_SIM_ERRORS) sim_error_msg(slot + 1) = msg
end subroutine set_error

subroutine sim_print(me)
    type(simulation_t), intent(in) :: me
    character(19) :: scheme
    scheme = 'Implicit  new dt = '
    if(me%explicit) scheme = 'Explicit  new dt = '
    if(me%rejected) scheme = 'Rejected  new dt = '
    print*, scheme, s_to_str(me%dt), '  t = ', s_to_str(me%t), '  no. = ', to_str(me%step)
end subroutine sim_print

function s_to_str(val) result(str_out)
    real(dp), intent(in) :: val
    character(5) :: str_out
    character(12) :: buff
    integer :: e
    buff=repeat(' ',12);
    if (val<10000) then
        e = floor(log10(val*1.001))
        select case(e)
            case(1);              write(buff,'(en11.1)') val
            case(2,3);            write(buff,'(f8.0)')  val
            case(0,-3,-6,-9,-12); write(buff,'(en11.2)') val
            case default;         write(buff,'(en11.0)') val
        end select
        select case(e)
            case(-12:-10); buff(7:7)='p'
            case(-9:-7);   buff(7:7)='n'
            case(-6:-4);   buff(7:7)='u'
            case(-3:-1);   buff(7:7)='m'
        end select
        buff(1:4) = buff(4:7)
        buff(5:5)='s'
    else if(val < 3600*100) then
        e = val/60
        write(buff,'(i2":"i2)') e/60, modulo(e,60)
        if (buff(4:4)==' ') buff(4:)='0'
    else if(val < 3600*10000) then
        write(buff,'(f5.0)') val/3600_8
        buff(5:5)='h'
    else
        write(buff,'(f5.0)') val/3600_8/24_8
        buff(5:5)='d'
    endif
    str_out = buff(1:5)
end function

function filter(simple, c1l, c0l, kl) !! Filter/controller with limiter used by step management     
    logical, intent(in) :: simple     !! True: Elementary filter, False: PI42 controller
    real(dp), intent(in) :: c1l, c0l, kl
    real(dp) :: filter
    real(dp), parameter :: kappa = 1., b = 5.
      
    if(simple) then                ! Elementary filter:
        filter = c1l**(1.0_dp/kl) 
    else                           ! PI42 controller:
        filter = c1l**(3.0_dp/b/kl) * c0l**(-1.0_dp/b/kl)
    endif                          ! Limiter:
    filter = 1.0_dp + kappa * atan((filter - 1.0_dp)/kappa);
end function filter

subroutine make_c_string(txt, buf, pos, ptr, len_out)
    character(*), intent(in)                         :: txt
    character(kind=c_char), intent(inout), target    :: buf(:)
    integer,               intent(inout)             :: pos
    type(c_ptr),           intent(out)               :: ptr
    integer(c_int),        intent(out)               :: len_out

    integer :: L, i, start

    L     = len(txt)
    start = pos

    if (start + L > size(buf)) then
        error stop 'make_c_string: buffer too small'
    end if

    do i = 1, L
        buf(start + i - 1) = txt(i:i)
    end do
    buf(start + L) = c_null_char

    ptr     = c_loc(buf(start))
    len_out = L + 1
    pos     = start + L + 1
end subroutine make_c_string

subroutine expose_state_to_ext_library(h5_ptr) bind(C)
    type(c_ptr) :: h5_ptr

    ! REIMS builds an H root structure (C-compatible) and transmits only the C pointer c_loc(H) to the DLL.
    ! H itself contains C pointers to other substructures that describe the entire original HDF hierarchy.      

    character(kind=c_char), allocatable, target, save :: table_names_buf(:)
    character(kind=c_char), allocatable, target, save :: comp_names_buf(:)
    integer :: nb_tables, t, c, h1, h2, total_comps
    integer :: total_name_chars, total_vars_chars, total_comp_name_chars
    integer :: i, pos_name, pos_vars, pos_cname, idx
    integer :: first_idx
    character(:), allocatable :: txt

    nb_tables = size(h5%tbl_dsc)

    allocate(Tbls_c(nb_tables))
    allocate(vars_offset(nb_tables))
    allocate(n_vars_tbl(nb_tables))
    allocate(var_len_tbl(nb_tables))

    total_comps = 0
    do t = 1, nb_tables
        total_comps = total_comps + size(h5%tbl_dsc(t)%comps)
    end do
    allocate(Comps_c(total_comps))

    total_name_chars = 0
    do t = 1, nb_tables
        total_name_chars = total_name_chars + len(h5%tbl_dsc(t)%name) + 1
    end do
    allocate(table_names_buf(total_name_chars))

    total_vars_chars = 0
    do t = 1, nb_tables
        h1 = size(h5%tbl_dsc(t)%vars)
        n_vars_tbl(t) = h1
        var_len_tbl(t) = 0
        do i = 1, h1
            var_len_tbl(t) = max(var_len_tbl(t), len(h5%tbl_dsc(t)%vars(i)))
        end do
        total_vars_chars = total_vars_chars + n_vars_tbl(t)*var_len_tbl(t)
    end do
    allocate(vars_buf(total_vars_chars))

    total_comp_name_chars = 0
    do t = 1, nb_tables
        do c = 1, size(h5%tbl_dsc(t)%comps)
            total_comp_name_chars = total_comp_name_chars + &
                 len(h5%tbl_dsc(t)%comps(c)%p%name) + 1
        end do
    end do
    allocate(comp_names_buf(total_comp_name_chars))

    pos_name  = 1
    pos_vars  = 1
    pos_cname = 1
    idx       = 0

    do t = 1, nb_tables

        txt = h5%tbl_dsc(t)%name
        call make_c_string(txt, table_names_buf, pos_name, &
                           Tbls_c(t)%name_ptr, Tbls_c(t)%name_len)

        h1 = n_vars_tbl(t)
        Tbls_c(t)%n_vars   = h1
        Tbls_c(t)%var_len  = var_len_tbl(t)

        if (h1 > 0) then
            vars_offset(t) = pos_vars
            do i = 1, h1
                txt = h5%tbl_dsc(t)%vars(i)
                do c = 1, var_len_tbl(t)
                    if (c <= len(txt)) then
                        vars_buf(pos_vars + c - 1) = txt(c:c)
                    else
                        vars_buf(pos_vars + c - 1) = ' '
                    end if
                end do
                pos_vars = pos_vars + var_len_tbl(t)
            end do
            Tbls_c(t)%vars_ptr = c_loc(vars_buf(vars_offset(t)))
        else
            vars_offset(t)   = 0
            Tbls_c(t)%vars_ptr = c_null_ptr
        end if

        Tbls_c(t)%n_comps = size(h5%tbl_dsc(t)%comps)
        if (Tbls_c(t)%n_comps > 0) then
            first_idx = idx + 1

            do i = 1, Tbls_c(t)%n_comps
                idx = idx + 1

                txt = h5%tbl_dsc(t)%comps(i)%p%name
                call make_c_string(txt, comp_names_buf, pos_cname, &
                                   Comps_c(idx)%name_ptr, Comps_c(idx)%name_len)

                h2 = size(h5%tbl_dsc(t)%comps(i)%p%node_x)
                Comps_c(idx)%node_x_ptr = c_loc(h5%tbl_dsc(t)%comps(i)%p%node_x(1))
                Comps_c(idx)%n_node     = h2

                Comps_c(idx)%data_ptr = c_loc(h5%tbl_dsc(t)%comps(i)%p%data(1,1))         
                Comps_c(idx)%n1       = h2
                Comps_c(idx)%n2       = h1
            end do

            Tbls_c(t)%comps_ptr = c_loc(Comps_c(first_idx))
        else
            Tbls_c(t)%comps_ptr = c_null_ptr
        end if

    end do

    H%nb_tables  = nb_tables
    H%tables_ptr = c_loc(Tbls_c(1))

    h5_ptr = c_loc(H)
end subroutine expose_state_to_ext_library

! -------------- I N T E R P O L A T I O N   F U N C T I O N S -----------------
subroutine signal_init(me, sim, input_cfg, key, x)
    class(signal_t),             intent(inout) :: me
    type(simulation_t), target, intent(inout) :: sim
    type(input_t),              intent(in)    :: input_cfg
    character(*),               intent(in)    :: key
    real(dp), optional,         intent(in)    :: x(:)

    logical :: explicit
    type(input_t) :: cfg
    integer :: len_int,i,j,repeat_idx,order
    real(dp) :: delta
    real(dp), allocatable :: x_in(:),y_out(:,:)
    character(:), allocatable :: event
    character(:), allocatable :: init_str
    type(str_ptr), allocatable :: keys(:)

    ! if(.not. input_cfg%has_key(key)) ! TODO by Jacek????

    me%sim => sim
    me%x = [0._dp]
    if (present(x)) me%x = x
    me%time   = [-huge(0.0_dp), huge(0.0_dp)]
    me%gain   = spread([0._dp], 1, size(me%x))
    me%offset = spread([0._dp], 1, size(me%x))
    if (.not.input_cfg%has_key(key)) return 

    me%offset = spread([input_cfg%dbl(key,ieee_value(delta, ieee_quiet_nan))], &
                                                                  1, size(me%x))
                                                                 
    if(.not.is_nan(me%offset(1,1))) return

    me%offset(:,1) = input_cfg%dbl1d(key,me%offset(:,1))
    if(.not.is_nan(me%offset(1,1))) return

    cfg = input_cfg%dict(key)
    init_str=cfg%str('external','')
    do i = 1, size(lib_handles)     
        me%init_signal_ext_ptr = GetProcAddress(lib_handles(i),init_str // c_null_char)
        if (c_associated(me%init_signal_ext_ptr)) then
            call c_f_procpointer(me%init_signal_ext_ptr, me%init_signal_ext)
            exit
        endif
    enddo

    if(associated(me%init_signal_ext)) then
      ! serialization
      keys = cfg%keys1d()
      me%cfg_str = ''
      do i = 1, size(keys)
          me%cfg_str = me%cfg_str // keys(i)%p // c_null_char // cfg%str(keys(i)%p) // c_null_char
      end do
      return
    endif

    repeat_idx = cfg%int('repeat',0)
    order = cfg%int('order',0)
    event = cfg%str('event','no')
      
    if(cfg%has_key('time')) then
        me%time = cfg%dbl1d('time')
        if (present(x)) then
            me%offset = cfg%dbl2d('value')
        else
            me%offset = reshape(cfg%dbl1d('value'),[1,size(me%time)])
        endif
    else
        allocate(y_out(size(x),1))
        call interpolate_vector(cfg%dbl1d('x'),cfg%dbl1d('value'),x,y_out(:,1))
        me%offset = y_out
        return
    endif
      
    if(cfg%has_key('x').and.cfg%has_key('time')) then
        x_in = cfg%dbl1d('x')
        allocate(y_out(size(me%x),size(me%offset,2)))
        do i = 1, size(me%offset,2)
            call interpolate_vector(x_in,me%offset(:,i),me%x,y_out(:,i))
        enddo
        me%offset = y_out
    endif
      
    ! Handling repeated time series
    if(repeat_idx>0) then
        do while(me%time(size(me%time))<=me%sim%t_final)
            len_int = size(me%time)
            delta = me%time(len_int) + me%time(repeat_idx)
            me%time = [me%time, me%time(  repeat_idx:len_int) + delta]
            me%offset = reshape([me%offset, me%offset(:,repeat_idx:len_int)], &
                                                    [size(me%x),size(me%time)])
        enddo
    endif
    len_int = size(me%time)
    me%time = [me%time, huge(0._dp)]
    me%gain = spread(me%gain(:,1),2,len_int)
  
    ! Handling the order of interpolation
    if(order == 1) then
        do j = 1,size(me%x)
            me%gain(j,1:len_int-1) = [((me%offset(j,i+1) - &
                me%offset(j,i))/(me%time(i+1) - me%time(i)), i=1, len_int-1)]
            me%gain(j,len_int) = me%gain(j,len_int-1)
            me%offset(j,:) = me%offset(j,:) - (me%gain(j,:) * me%time(1:len_int))
        enddo
    endif
  
    ! Adding events to list
    if     (event == 'explicit') then; explicit = .true.
    elseif (event == 'implicit') then; explicit = .false.
    else; return
    endif
    call merge_events(me%sim%events_time, me%sim%events_expl, me%time, explicit)
end subroutine signal_init

function signal_v1d(me) result(val)
    class(signal_t), intent(inout) :: me
    real(dp) :: val(size(me%x)) !! interpolated value
    real(dp) :: time
    type(c_ptr) :: c_cfg
    time = me%sim%t + me%sim%dt

    if(associated(me%init_signal_ext) .and. .not.c_associated(me%sign_ptr)) then
      c_cfg = C_NULL_PTR
      if(.not. me%init_signal_ext(c_funloc(expose_state_to_ext_library),len(me%cfg_str), me%cfg_str, c_cfg, me%sign_ptr,&
                                  size(me%x))) return
      if(c_associated(me%sign_ptr)) then
          call c_f_procpointer(me%sign_ptr, me%val_fpt)
          me%val_cfg = c_cfg
      endif
    endif

    if (associated(me%val_fpt)) then
        call me%val_fpt(me%val_cfg, time, val)
        return
    endif
    do ! search for correct segment index starting from previous point: me%i
        if (time >= me%time(me%i + 1)) then
            me%i = me%i + 1
        else if (time < me%time(me%i)) then
            if (me%i == 1) exit
            me%i = me%i - 1
        else
            exit
        endif
    enddo
    val = time * me%gain(:,me%i) + me%offset(:,me%i)
end function signal_v1d

function signal_v0d(me) result(val)
    class(signal_t), intent(inout) :: me
    real(dp) :: val !! interpolated value
    real(dp) :: tmp(1)
    tmp = signal_v1d(me)
    val = tmp(1)
end function signal_v0d

subroutine interpolate(x, y, x_out, y_out) ! TODO this should be a function not subroutine
    real(dp), intent(in) :: x(:), y(:), x_out
    real(dp), intent(out) :: y_out
    integer :: i

    ! Check if x_out is outside the range of x
    if (x_out <= x(1)) then
        y_out = y(1)
    elseif (x_out >= x(size(x))) then
        y_out = y(size(y))
    else
        ! Find the interval [x(i), x(i+1)] that contains x_out
        do i = 1, size(x) - 1
            if (x_out >= x(i) .and. x_out <= x(i+1)) then
                ! Perform linear interpolation
                y_out = y(i) + (y(i+1) - y(i)) * (x_out - x(i)) / (x(i+1) - x(i))
                exit
            endif
        enddo
    endif
end subroutine interpolate

subroutine interpolate_vector(x, y, x_out, y_out) ! TODO this should be a function not subroutine
    real(dp), intent(in) :: x(:), y(:), x_out(:)
    real(dp), intent(out) :: y_out(:)
    integer :: i
    do i = 1, size(x_out)
        call interpolate(x,y,x_out(i),y_out(i))
    enddo
end subroutine interpolate_vector

subroutine merge_events(time_inout, expl_inout, time, explicit)
    !! Merging sorted events of simulation
    !!
    !! Both input arrays should be sorted and the  last elements
    !! should alvays be huge() and they will be not sorted
    real(dp), allocatable, intent(inout) :: time_inout(:)
    logical,  allocatable, intent(inout) :: expl_inout(:)
    real(dp), intent(in) :: time(:) !! sorted, without duplicates input arrays ending with huge
    logical,  intent(in) :: explicit

    real(dp), allocatable :: tmp(:) !! temporary time event array
    logical,  allocatable :: msk(:) !! temporary explicit event array
    integer :: i1, i2, o !! input, output indexes

    allocate(tmp(size(time_inout) + size(time)))
    allocate(msk(size(tmp)))
    o = 0; i1 = 1; i2 = 1
    do while (i1 <= size(time_inout) .and. i2 <= size(time))
        o = o + 1
        if (abs(time_inout(i1)-time(i2))<=sim_delta) then
            tmp(o) = time_inout(i1)
            msk(o) = expl_inout(i1) .or. explicit
            i1 = i1 + 1
            i2 = i2 + 1
        else if (time_inout(i1)<time(i2)) then
            tmp(o) = time_inout(i1)
            msk(o) = expl_inout(i1)
            i1 = i1 + 1
        else
            tmp(o) = time(i2)
            msk(o) = explicit
            i2 = i2 + 1
        endif
    enddo
    time_inout = tmp(1:o)
    expl_inout = msk(1:o)
end subroutine merge_events

end module krn_simulation_m