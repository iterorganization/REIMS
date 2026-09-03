! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module lib_nusselt_correlations_m
  use krn_global_tools_m
  use iso_c_binding
  use fortran_yaml_c  
  use lib_input_m
  implicit none

  type nusselt_t
      character(:), allocatable :: nusselt_corr
      type(c_ptr) :: nuss_cfg = C_NULL_PTR
      procedure(itf_4var), pointer, nopass :: nuss => nusselt_correlation_not_implemented
      procedure(itf_4var), pointer, nopass :: nuss_der_Re   => null()
      procedure(itf_4var), pointer, nopass :: nuss_der_Pra  => null()
      procedure(itf_4var), pointer, nopass :: nuss_der_Tother => null()
      procedure(itf_4var), pointer, nopass :: nuss_der_T    => null()
  contains
      procedure :: nusselt            => nusselt_correlation
      procedure :: nusselt_der_Re     => nusselt_correlation_der_Re
      procedure :: nusselt_der_Pra    => nusselt_correlation_der_Pra
      procedure :: nusselt_der_Tother => nusselt_correlation_der_Tother
      procedure :: nusselt_der_T      => nusselt_correlation_der_T
  end type nusselt_t
  type(nusselt_t), allocatable, target :: nusselts(:)  

  abstract interface
    function init_nusselt_ext_if(cfg_len,cfg_str,c_pt,nuss_ptr) bind(C)
      import :: c_int, c_ptr, c_funptr, dp
      integer(c_int), value :: cfg_len
      character(*)    :: cfg_str 
      type(c_ptr)     :: c_pt      
      type(c_funptr)  :: nuss_ptr
      logical         :: init_nusselt_ext_if      
    end function 
  end interface

  procedure(init_nusselt_ext_if), pointer :: init_nusselt_ext  => null()  

  type pipe_cfg_t  
     real(dp) :: a, b, c
  end type pipe_cfg_t  

  type DBG_cfg_t  
     real(dp) :: NuL, a, b, c, d
  end type DBG_cfg_t  

  contains

    subroutine load_nusselts(cfgs)
        class(input_t), intent(in) :: cfgs(:)
        
        character(20)  :: names(2) = ['pipe','DBG']
        integer        :: i
        type(c_funptr) :: init_nusselt_ext_ptr

        do i = 1, size(lib_handles)     
            init_nusselt_ext_ptr = GetProcAddress(lib_handles(i), 'init_nusselt_ext' // c_null_char)
            if (c_associated(init_nusselt_ext_ptr)) then
                call c_f_procpointer(init_nusselt_ext_ptr, init_nusselt_ext)
                exit
            endif
        enddo

        allocate(nusselts(size(cfgs) + size(names)))
        do i = 1, size(cfgs)
            call set_nusselt_correlation(nusselts(i), cfgs(i), '')
        enddo 
        do i = 1, size(names)
            call set_nusselt_correlation(nusselts(size(cfgs)+i), empty_cfg(), names(i))
        enddo 
    end subroutine load_nusselts    


    subroutine set_nusselt_correlation(nss, cfg, built_in_name)
        type(nusselt_t), intent(inout) :: nss
        class(input_t),   intent(in)    :: cfg
        character(*),     intent(in)    :: built_in_name
        
        type(c_funptr) :: nusselt_ptr
        integer  :: i

        character(:), allocatable :: cfg_str
        type(str_ptr), allocatable :: keys(:)
        type(c_ptr) :: c_cfg
        c_cfg = C_NULL_PTR

        nss%nusselt_corr = cfg%str('nusselt',trim(built_in_name))
        select case (cfg%str('base',nss%nusselt_corr))
        case('pipe')
            call nusselt_pipe_init(cfg, c_cfg)
            nss%nuss_cfg        =  c_cfg
            nss%nuss            => nusselt_pipe
        case('DBG')
            call nusselt_DBG_init(cfg, c_cfg)
            nss%nuss_cfg        =  c_cfg
            nss%nuss            => nusselt_DBG
            nss%nuss_der_Re     => nusselt_DBG_der_Re
            nss%nuss_der_Pra    => nusselt_DBG_der_Pra
            nss%nuss_der_Tother => nusselt_DBG_der_Tother
            nss%nuss_der_T      => nusselt_DBG_der_T
        end select

        if(.not. associated(init_nusselt_ext)) return

        ! serialization
        keys = cfg%keys1d()
        cfg_str = ''
        do i = 1, size(keys)
            cfg_str =  cfg_str // keys(i)%p // c_null_char // cfg%str(keys(i)%p) // c_null_char
        end do

        nusselt_ptr = C_NULL_FUNPTR
        if(.not. init_nusselt_ext(len(cfg_str), cfg_str, c_cfg, nusselt_ptr)) return  

        if(c_associated(nusselt_ptr)) then
            call c_f_procpointer(nusselt_ptr, nss%nuss)
            nss%nuss_cfg = c_cfg
        endif        
    end subroutine set_nusselt_correlation


    subroutine nusselt_init(me,cfg)
        class(nusselt_t), pointer, intent(inout) :: me
        class(input_t), target, intent(in) :: cfg
    
        integer :: i
        character(:), allocatable :: nss   
        
        nss=cfg%str('nusselt','')
        do i=1,size(nusselts)
            if(nss==nusselts(i)%nusselt_corr) me => nusselts(i)
        enddo
    end subroutine nusselt_init


    subroutine nusselt_pipe_init(cfg,c_pt)  
        class(input_t), target, intent(in) :: cfg
        type(c_ptr), intent(out) :: c_pt
        type(pipe_cfg_t), pointer :: pipe_ptr ! fortran pointer
  
        allocate(pipe_ptr)
        pipe_ptr%a = cfg%dbl('a',0.023_dp)
        pipe_ptr%b = cfg%dbl('b',0.8_dp)
        pipe_ptr%c = cfg%dbl('c',0.4_dp)

        c_pt = c_loc(pipe_ptr) ! translates fortran pointer to C pointer
    end subroutine nusselt_pipe_init      
  
     
    function nusselt_pipe(c_pt,Re,Pra,Tother,T) result(nusselt) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: Re,Pra,Tother,T
        real(dp) :: nusselt
        type(pipe_cfg_t), pointer :: pipe_ptr
        call c_f_pointer(c_pt,pipe_ptr)      
    
        nusselt=pipe_ptr%a*(Re**pipe_ptr%b)*(Pra**pipe_ptr%c)
    end function nusselt_pipe


    subroutine nusselt_DBG_init(cfg,c_pt)  
        class(input_t), target, intent(in) :: cfg
        type(c_ptr), intent(out) :: c_pt
        type(DBG_cfg_t), pointer :: DBG_ptr ! fortran pointer
  
        allocate(DBG_ptr)
        DBG_ptr%NuL = cfg%dbl('NuL',8.235_dp)
        DBG_ptr%a   = cfg%dbl('a',0.0259_dp)
        DBG_ptr%b   = cfg%dbl('b',0.8_dp)
        DBG_ptr%c   = cfg%dbl('c',0.4_dp)
        DBG_ptr%d   = cfg%dbl('d',-0.716_dp)

        c_pt = c_loc(DBG_ptr) ! translates fortran pointer to C pointer
    end subroutine nusselt_DBG_init      


    function nusselt_DBG(c_pt,Re,Pra,Tother,T) result(nusselt) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: Re,Pra,Tother,T
        real(dp) :: nusselt,NuT
        type(DBG_cfg_t), pointer :: DBG_ptr
        call c_f_pointer(c_pt,DBG_ptr)      
     
        NuT=DBG_ptr%a*(Re**DBG_ptr%b)*(Pra**DBG_ptr%c)*((Tother/T)**DBG_ptr%d)
        if(Re<1.0e-12_dp .or. Pra<0.0_dp) NuT=0.0_dp

        nusselt=max(DBG_ptr%NuL,NuT)
    end function nusselt_DBG   


    function nusselt_DBG_der_Re(c_pt,Re,Pra,Tother,T) result(nusselt_der_Re) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: Re,Pra,Tother,T
        real(dp) :: nusselt_der_Re,NuT
        type(DBG_cfg_t), pointer :: DBG_ptr
        call c_f_pointer(c_pt,DBG_ptr)      
     
        NuT=DBG_ptr%a*(Re**DBG_ptr%b)*(Pra**DBG_ptr%c)*((Tother/T)**DBG_ptr%d)
        if(Re<1.0e-12_dp) NuT=0.0_dp
  
        if(DBG_ptr%NuL>NuT) then
          nusselt_der_Re=0.0_dp
        else
          nusselt_der_Re=DBG_ptr%a*DBG_ptr%b*(Re**(DBG_ptr%b-1.0_dp))*(Pra**DBG_ptr%c)*((Tother/T)**DBG_ptr%d)
        endif
    end function nusselt_DBG_der_Re    
    
    
    function nusselt_DBG_der_Pra(c_pt,Re,Pra,Tother,T) result(nusselt_der_Pra) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: Re,Pra,Tother,T
        real(dp) :: nusselt_der_Pra,NuT
        type(DBG_cfg_t), pointer :: DBG_ptr
        call c_f_pointer(c_pt,DBG_ptr)      
     
        NuT=DBG_ptr%a*(Re**DBG_ptr%b)*(Pra**DBG_ptr%c)*((Tother/T)**DBG_ptr%d)
        if(Re<1.0e-12_dp) NuT=0.0_dp
  
        if(DBG_ptr%NuL>NuT) then
          nusselt_der_Pra=0.0_dp
        else
          nusselt_der_Pra=DBG_ptr%a*DBG_ptr%c*(Re**DBG_ptr%b)*(Pra**(DBG_ptr%c-1.0_dp))*((Tother/T)**DBG_ptr%d)
        endif
    end function nusselt_DBG_der_Pra    
    
    
    function nusselt_DBG_der_Tother(c_pt,Re,Pra,Tother,T) result(nusselt_der_Tother) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: Re,Pra,Tother,T
        real(dp) :: nusselt_der_Tother,NuT
        type(DBG_cfg_t), pointer :: DBG_ptr
        call c_f_pointer(c_pt,DBG_ptr)      

        NuT=DBG_ptr%a*(Re**DBG_ptr%b)*(Pra**DBG_ptr%c)*((Tother/T)**DBG_ptr%d)
        if(Re<1.0e-12_dp) NuT=0.0_dp
    
        if(DBG_ptr%NuL>NuT) then
          nusselt_der_Tother=0.0_dp
        else
          nusselt_der_Tother=DBG_ptr%a*(Re**DBG_ptr%b)*(Pra**DBG_ptr%c)*DBG_ptr%d*(1.0_dp/T)*&
                                 ((Tother/T)**(DBG_ptr%d-1.0_dp))
        endif        
    end function nusselt_DBG_der_Tother    
    
    
    function nusselt_DBG_der_T(c_pt,Re,Pra,Tother,T) result(nusselt_der_T) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: Re,Pra,Tother,T
        real(dp) :: nusselt_der_T,NuT
        type(DBG_cfg_t), pointer :: DBG_ptr
        call c_f_pointer(c_pt,DBG_ptr)      
     
        NuT=DBG_ptr%a*(Re**DBG_ptr%b)*(Pra**DBG_ptr%c)*((Tother/T)**DBG_ptr%d)
        if(Re<1.0e-12_dp) NuT=0.0_dp
    
        if(DBG_ptr%NuL>NuT) then
          nusselt_der_T=0.0_dp
        else
          nusselt_der_T=DBG_ptr%a*(Re**DBG_ptr%b)*(Pra**DBG_ptr%c)*DBG_ptr%d*(-Tother/(T**2))*&
                             ((Tother/T)**(DBG_ptr%d-1.0_dp))
        endif
    end function nusselt_DBG_der_T        


    ! Member functions implementations --- calls to lower-level function pointers  
    function nusselt_correlation(me,Re,Pra,Tother,T)
      class(nusselt_t),  intent(inout) :: me
      real(dp), intent(in) :: Re,Pra
      real(dp), intent(in), optional :: Tother,T
      real(dp) :: nusselt_correlation,Tother_loc,T_loc

      if(present(Tother) .and. present(T)) then
          Tother_loc=Tother; T_loc=T
      else
          Tother_loc=0.0_dp; T_loc=0.0_dp
      endif

      nusselt_correlation=me%nuss(me%nuss_cfg,Re,Pra,Tother_loc,T_loc)
    end function nusselt_correlation
  
  
    function nusselt_correlation_der_Re(me,Re,Pra,Tother,T)
      class(nusselt_t),  intent(inout) :: me
      real(dp), intent(in) :: Re,Pra
      real(dp), intent(in), optional :: Tother,T
      real(dp) :: nusselt_correlation_der_Re,Tother_loc,T_loc
      real(dp) :: eps

      if(present(Tother) .and. present(T)) then
          Tother_loc=Tother; T_loc=T
      else
          Tother_loc=0.0_dp; T_loc=0.0_dp
      endif

      if(associated(me%nuss_der_Re)) then
          nusselt_correlation_der_Re=me%nuss_der_Re(me%nuss_cfg,Re,Pra,Tother_loc,T_loc)
      else
          eps=2.0e-1_dp
          if(Re >= 0.5_dp*eps) then
            nusselt_correlation_der_Re=(me%nuss(me%nuss_cfg,Re+0.5_dp*eps,Pra,Tother_loc,T_loc) &
                                       -me%nuss(me%nuss_cfg,Re-0.5_dp*eps,Pra,Tother_loc,T_loc))/eps
          else
            nusselt_correlation_der_Re=(me%nuss(me%nuss_cfg,Re+eps,Pra,Tother_loc,T_loc) &
                                       -me%nuss(me%nuss_cfg,Re,Pra,Tother_loc,T_loc))/eps
          endif
      endif
    end function nusselt_correlation_der_Re
  
  
    function nusselt_correlation_der_Pra(me,Re,Pra,Tother,T)
      class(nusselt_t),  intent(inout) :: me
      real(dp), intent(in) :: Re,Pra
      real(dp), intent(in), optional :: Tother,T
      real(dp) :: nusselt_correlation_der_Pra,Tother_loc,T_loc
      real(dp) :: eps
      
      if(present(Tother) .and. present(T)) then
          Tother_loc=Tother; T_loc=T
      else
          Tother_loc=0.0_dp; T_loc=0.0_dp
      endif

      if(associated(me%nuss_der_Pra)) then
          nusselt_correlation_der_Pra=me%nuss_der_Pra(me%nuss_cfg,Re,Pra,Tother_loc,T_loc)
      else
          eps=2.0e-4_dp
          if(Pra >= 0.5_dp*eps) then
            nusselt_correlation_der_Pra=(me%nuss(me%nuss_cfg,Re,Pra+0.5_dp*eps,Tother_loc,T_loc) &
                                        -me%nuss(me%nuss_cfg,Re,Pra-0.5_dp*eps,Tother_loc,T_loc))/eps
          else
            nusselt_correlation_der_Pra=(me%nuss(me%nuss_cfg,Re,Pra+eps,Tother_loc,T_loc) &
                                        -me%nuss(me%nuss_cfg,Re,Pra,Tother_loc,T_loc))/eps
          endif
      endif
    end function nusselt_correlation_der_Pra
  
  
    function nusselt_correlation_der_Tother(me,Re,Pra,Tother,T)
      class(nusselt_t),  intent(inout) :: me
      real(dp), intent(in) :: Re,Pra
      real(dp), intent(in), optional :: Tother,T
      real(dp) :: nusselt_correlation_der_Tother,Tother_loc,T_loc
      real(dp) :: eps
      
      if(present(Tother) .and. present(T)) then
          Tother_loc=Tother; T_loc=T
      else
          Tother_loc=0.0_dp; T_loc=0.0_dp
      endif
  
      if(associated(me%nuss_der_Tother)) then
          nusselt_correlation_der_Tother=me%nuss_der_Tother(me%nuss_cfg,Re,Pra,Tother_loc,T_loc)
      else        
          eps=2.0e-1_dp
          nusselt_correlation_der_Tother=(me%nuss(me%nuss_cfg,Re,Pra,Tother_loc+0.5_dp*eps,T_loc) &
                                       -me%nuss(me%nuss_cfg,Re,Pra,Tother_loc-0.5_dp*eps,T_loc))/eps
      endif
    end function nusselt_correlation_der_Tother
  
  
    function nusselt_correlation_der_T(me,Re,Pra,Tother,T)
      class(nusselt_t),  intent(inout) :: me
      real(dp), intent(in) :: Re,Pra
      real(dp), intent(in), optional :: Tother,T
      real(dp) :: nusselt_correlation_der_T,Tother_loc,T_loc
      real(dp) :: eps
      
      if(present(Tother) .and. present(T)) then
          Tother_loc=Tother; T_loc=T
      else
          Tother_loc=0.0_dp; T_loc=0.0_dp
      endif

      if(associated(me%nuss_der_T)) then
          nusselt_correlation_der_T=me%nuss_der_T(me%nuss_cfg,Re,Pra,Tother_loc,T_loc)
      else        
          eps=2.0e-1_dp
          nusselt_correlation_der_T=(me%nuss(me%nuss_cfg,Re,Pra,Tother_loc,T_loc+0.5_dp*eps) &
                                    -me%nuss(me%nuss_cfg,Re,Pra,Tother_loc,T_loc-0.5_dp*eps))/eps
      endif
    end function nusselt_correlation_der_T
  
    
    ! Dummy function - not implemented, to handle errors when nusselt does not set or override them
    function nusselt_correlation_not_implemented(c_pt,Re,Pra,Tother,T) result(ret_val) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: Re,Pra,Tother,T
      real(dp) :: ret_val
      ret_val=0.0_dp
      if(.false.) print*,c_associated(c_pt),Re,Pra,Tother,T
      print*,'nusselt_correlation - not implemented'
      stop
    end function nusselt_correlation_not_implemented      
    
end module lib_nusselt_correlations_m
