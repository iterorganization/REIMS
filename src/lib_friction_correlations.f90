! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module lib_friction_correlations_m
  use krn_global_tools_m
  use iso_c_binding
  use fortran_yaml_c  
  use lib_input_m
  implicit none

  type friction_t
      character(:), allocatable :: friction_corr
      type(c_ptr) :: frict_cfg = C_NULL_PTR
      procedure(itf_1var), pointer, nopass :: frict => friction_correlation_not_implemented
      procedure(itf_1var), pointer, nopass :: frict_der => null()
  contains
      procedure :: friction     => friction_correlation
      procedure :: friction_der => friction_correlation_der  
  end type friction_t
  type(friction_t), allocatable, target :: frictions(:)     

  abstract interface
    function init_friction_ext_if(cfg_len,cfg_str,c_pt,frict_ptr) bind(C)
      import :: c_int, c_ptr, c_funptr, dp
      integer(c_int), value :: cfg_len
      character(*)    :: cfg_str 
      type(c_ptr)     :: c_pt      
      type(c_funptr)  :: frict_ptr
      logical         :: init_friction_ext_if      
    end function 
  end interface

  procedure(init_friction_ext_if), pointer :: init_friction_ext  => null()  

  type blasius_cfg_t  
     real(dp) :: alpha, beta
  end type blasius_cfg_t  

  type katheder_cfg_t
     real(dp) :: alpha, beta, VoidFr, Re_min
  end type katheder_cfg_t
  
  type central_spiral_cfg_t  
     real(dp) :: alpha, beta
  end type central_spiral_cfg_t  
            
  contains

    subroutine load_frictions(cfgs)
        class(input_t), intent(in) :: cfgs(:)
        
        character(20)  :: names(3) = ['blasius','katheder','central_spiral']
        integer        :: i
        type(c_funptr) :: init_friction_ext_ptr

        do i = 1, size(lib_handles)     
            init_friction_ext_ptr = GetProcAddress(lib_handles(i), 'init_friction_ext' // c_null_char)
            if (c_associated(init_friction_ext_ptr)) then
                call c_f_procpointer(init_friction_ext_ptr, init_friction_ext)
                exit
            endif
        enddo

        allocate(frictions(size(cfgs) + size(names)))
        do i = 1, size(cfgs)
            call set_friction_correlation(frictions(i), cfgs(i), '')
        enddo 
        do i = 1, size(names)
            call set_friction_correlation(frictions(size(cfgs)+i), empty_cfg(), names(i))
        enddo 
    end subroutine load_frictions


    subroutine set_friction_correlation(fct, cfg, built_in_name)
        type(friction_t), intent(inout) :: fct
        class(input_t),   intent(in)    :: cfg
        character(*),     intent(in)    :: built_in_name
        
        type(c_funptr) :: frict_ptr
        integer  :: i

        character(:), allocatable :: cfg_str
        type(str_ptr), allocatable :: keys(:)
        type(c_ptr) :: c_cfg
        c_cfg = C_NULL_PTR

        fct%friction_corr = cfg%str('friction',trim(built_in_name))
        select case (cfg%str('base',fct%friction_corr))
        case('blasius')
            call friction_factor_blasius_init(cfg, c_cfg)
            fct%frict_cfg      =  c_cfg
            fct%frict          => friction_factor_blasius
            fct%frict_der      => friction_factor_blasius_der
        case('katheder')
            call friction_factor_katheder_init(cfg, c_cfg)
            fct%frict_cfg      =  c_cfg
            fct%frict          => friction_factor_katheder
            fct%frict_der      => friction_factor_katheder_der
        case('central_spiral')
            call friction_factor_central_spiral_init(cfg, c_cfg)
            fct%frict_cfg      =  c_cfg
            fct%frict          => friction_factor_central_spiral
            fct%frict_der      => friction_factor_central_spiral_der
        end select

        if(.not. associated(init_friction_ext)) return

        ! serialization
        keys = cfg%keys1d()
        cfg_str = ''
        do i = 1, size(keys)
            cfg_str =  cfg_str // keys(i)%p // c_null_char // cfg%str(keys(i)%p) // c_null_char
        end do

        frict_ptr = C_NULL_FUNPTR
        if(.not. init_friction_ext(len(cfg_str), cfg_str, c_cfg, frict_ptr)) return  

        if(c_associated(frict_ptr)) then
            call c_f_procpointer(frict_ptr, fct%frict)
            fct%frict_cfg = c_cfg
        endif        
    end subroutine set_friction_correlation


   subroutine friction_init(me,cfg)
        class(friction_t), pointer, intent(inout) :: me
        class(input_t), target, intent(in) :: cfg
    
        integer :: i
        character(:), allocatable :: fct   
        
        fct=cfg%str('friction')
        do i=1,size(frictions)
            if(fct==frictions(i)%friction_corr) me => frictions(i)
        enddo
        if(.NOT.associated(me)) error stop fct//'friction does not exist'
   end subroutine friction_init   


  subroutine friction_factor_blasius_init(cfg,c_pt)  
      class(input_t), target, intent(in) :: cfg
      type(c_ptr), intent(out) :: c_pt
      type(blasius_cfg_t), pointer :: blasius_ptr ! fortran pointer

      allocate(blasius_ptr)
      blasius_ptr%alpha = cfg%dbl('alpha',0.079_dp)      
      blasius_ptr%beta  = cfg%dbl('beta',0.250_dp)     

      c_pt = c_loc(blasius_ptr) ! translates fortran pointer to C pointer
  end subroutine friction_factor_blasius_init      

   
  function friction_factor_blasius(c_pt,Re) result(friction) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: Re
      real(dp) :: friction  
      real(dp) :: Remin
      type(blasius_cfg_t), pointer :: blasius_ptr
      call c_f_pointer(c_pt,blasius_ptr)      
  
      Remin=10.0_dp  
      if(Re>Remin) then
          friction=blasius_ptr%alpha*(Re**(-blasius_ptr%beta))
      else
          friction=blasius_ptr%alpha*(Remin**(-blasius_ptr%beta))
      endif
  end function friction_factor_blasius   
  
  
  function friction_factor_blasius_der(c_pt,Re) result(friction_der) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: Re
      real(dp) :: friction_der
      real(dp) :: Remin
      type(blasius_cfg_t), pointer :: blasius_ptr
      call c_f_pointer(c_pt,blasius_ptr)      
  
      Remin=10.0_dp  
      if(Re>Remin) then
          friction_der=-blasius_ptr%alpha*blasius_ptr%beta*(Re**(-blasius_ptr%beta-1.0_dp))
      else
          friction_der=0.0_dp
      endif      
  end function friction_factor_blasius_der     


  subroutine friction_factor_katheder_init(cfg,c_pt)  
      class(input_t), target, intent(in) :: cfg
      type(c_ptr), intent(out) :: c_pt
      type(katheder_cfg_t), pointer :: katheder_ptr ! fortran pointer

      allocate(katheder_ptr)
      katheder_ptr%alpha  = cfg%dbl('alpha',0.0231_dp)
      katheder_ptr%beta   = cfg%dbl('beta',0.7953_dp)
      katheder_ptr%VoidFr = cfg%dbl('VoidFr',0.297_dp)
      katheder_ptr%Re_min = cfg%dbl('Re_min',1000.0_dp)

      c_pt = c_loc(katheder_ptr) ! translates fortran pointer to C pointer
  end subroutine friction_factor_katheder_init 


  function friction_factor_katheder(c_pt,Re) result(friction) bind(C)
    type(c_ptr), value :: c_pt  
    real(dp), value :: Re
    real(dp) :: friction
    type(katheder_cfg_t), pointer :: katheder_ptr
    call c_f_pointer(c_pt,katheder_ptr)     

    if(Re > katheder_ptr%Re_min) then
      friction=0.25_dp*(19.5_dp/(Re**katheder_ptr%beta)+katheder_ptr%alpha)/(katheder_ptr%VoidFr**0.742_dp)
    else
      friction=0.0_dp
    endif
  end function friction_factor_katheder


  function friction_factor_katheder_der(c_pt,Re) result(friction_der) bind(C)
    type(c_ptr), value :: c_pt
    real(dp), value :: Re
    real(dp) :: friction_der
    type(katheder_cfg_t), pointer :: katheder_ptr
    call c_f_pointer(c_pt,katheder_ptr)

    if(Re > katheder_ptr%Re_min) then
      friction_der=-0.25_dp*19.5_dp*katheder_ptr%beta*(Re**(-katheder_ptr%beta-1.0_dp))/(katheder_ptr%VoidFr**0.742_dp)
    else
      friction_der=0.0_dp
    endif
  end function friction_factor_katheder_der 


  subroutine friction_factor_central_spiral_init(cfg,c_pt)  
      class(input_t), target, intent(in) :: cfg
      type(c_ptr), intent(out) :: c_pt
      type(central_spiral_cfg_t), pointer :: central_spiral_ptr ! fortran pointer

      allocate(central_spiral_ptr)
      central_spiral_ptr%alpha  = cfg%dbl('alpha',0.36_dp)      
      central_spiral_ptr%beta   = cfg%dbl('beta',0.038_dp)

      c_pt = c_loc(central_spiral_ptr) ! translates fortran pointer to C pointer
  end subroutine friction_factor_central_spiral_init     


  function friction_factor_central_spiral(c_pt,Re) result(friction) bind(C)
    type(c_ptr), value :: c_pt
    real(dp), value :: Re
    real(dp) :: friction
    type(central_spiral_cfg_t), pointer :: central_spiral_ptr
    call c_f_pointer(c_pt,central_spiral_ptr)    

    if(Re>0.0_dp) then
      friction=0.25_dp*central_spiral_ptr%alpha/Re**central_spiral_ptr%beta
    else
      friction=0.0_dp
    endif
  end function friction_factor_central_spiral      


  function friction_factor_central_spiral_der(c_pt,Re) result(friction_der) bind(C)
    type(c_ptr), value :: c_pt
    real(dp), value :: Re
    real(dp) :: friction_der
    type(central_spiral_cfg_t), pointer :: central_spiral_ptr
    call c_f_pointer(c_pt,central_spiral_ptr)    

    if(Re>0.0_dp) then
      friction_der=-0.25_dp*central_spiral_ptr%alpha*central_spiral_ptr%beta*(Re**(-central_spiral_ptr%beta-1.0_dp))
    else
      friction_der=0.0_dp
    endif    
  end function friction_factor_central_spiral_der      
 

  ! Member functions implementations --- calls to lower-level function pointers  
  function friction_correlation(me,Re)
    class(friction_t),  intent(inout) :: me
    real(dp), intent(in) :: Re
    real(dp) :: friction_correlation

    friction_correlation=me%frict(me%frict_cfg,Re)
  end function friction_correlation


  function friction_correlation_der(me,Re)
    class(friction_t),  intent(inout) :: me
    real(dp), intent(in) :: Re
    real(dp) :: friction_correlation_der
    real(dp) :: eps
    if(associated(me%frict_der)) then
        friction_correlation_der=me%frict_der(me%frict_cfg,Re)
    else
        eps=2.0e-1_dp
        friction_correlation_der=(me%frict(me%frict_cfg,Re+0.5_dp*eps) &
                                 -me%frict(me%frict_cfg,Re-0.5_dp*eps))/eps
    endif
  end function friction_correlation_der

  
  ! Dummy function - not implemented, to handle errors when friction does not set or override them
  function friction_correlation_not_implemented(c_pt,Re) result(ret_val) bind(C)
    type(c_ptr), value :: c_pt
    real(dp), value :: Re
    real(dp) :: ret_val
    ret_val=0.0_dp
    if(.false.) print*,c_associated(c_pt),Re
    print*,'friction_correlation - not implemented'
    stop
  end function friction_correlation_not_implemented   


  function friction_factor_laminar(Re)  
    real(dp), intent(in) :: Re
    real(dp) :: friction_factor_laminar
    real(dp) :: a,Remin

    a=16.0_dp
    Remin=1.0e-8_dp

    if(Re>Remin) then
      friction_factor_laminar=a/Re
    else
      friction_factor_laminar=a/Remin
    endif
  end function friction_factor_laminar
  

  function der_friction_factor_laminar(Re)
    real(dp), intent(in) :: Re
    real(dp) :: der_friction_factor_laminar
    real(dp) :: a,Remin

    a=16.0_dp
    Remin=1.0e-8_dp

    if(Re>Remin) then
      der_friction_factor_laminar=-a/(Re**2)
    else
      der_friction_factor_laminar=0.0_dp
    endif
  end function der_friction_factor_laminar    
    
end module lib_friction_correlations_m