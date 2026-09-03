! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module lib_material_insulation_m
  use krn_global_tools_m
  use iso_c_binding
  use lib_input_m, only: input_t

  implicit none

  type glass_epoxy_cfg_t
  end type glass_epoxy_cfg_t    

  type glass_kapton_glass_cfg_t
  end type glass_kapton_glass_cfg_t

  contains

  subroutine material_glass_epoxy_init(cfg,density,c_pt)  
      class(input_t), target, intent(in) :: cfg
      real(dp), intent(out) :: density
      type(c_ptr), intent(out) :: c_pt

      type(glass_epoxy_cfg_t), pointer :: glass_epoxy_ptr ! fortran pointer

      density=cfg%dbl('density',1948.0_dp)

      allocate(glass_epoxy_ptr)
      c_pt = c_loc(glass_epoxy_ptr) ! translates fortran pointer to C pointer
  end subroutine material_glass_epoxy_init   


  function thermal_conductivity_glass_epoxy(c_pt,T,B) result(thermal_conductivity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B
      real(dp) :: thermal_conductivity    
      real(dp) :: Tmin,Tmax,TT,cc,dd,c,d,nc,nd,t0,aa,a,na,bb,bbb,nb
      type(glass_epoxy_cfg_t), pointer :: glass_epoxy_ptr
      call c_f_pointer(c_pt,glass_epoxy_ptr)

      Tmin = 1.0_dp
      Tmax = 300.0_dp
      TT=min(T,Tmax)
      TT=max(TT,Tmin)

      thermal_conductivity = 3.81925332e-8_dp*TT**3 - 2.24517344e-5_dp*TT**2 + 0.00591967022_dp*TT + 0.0645106161_dp
  end function thermal_conductivity_glass_epoxy  


  function heat_capacity_glass_epoxy(c_pt,T,B,TC,TCS,TC0) result(heat_capacity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B,TC,TCS,TC0
      real(dp) :: heat_capacity    
      real(dp) :: Tmin,Tmax,TT,u,w,Rb,Rh
      real(dp), parameter :: T1 = 20.0_dp, T2 = 30.0_dp
      type(glass_epoxy_cfg_t), pointer :: glass_epoxy_ptr
      call c_f_pointer(c_pt,glass_epoxy_ptr)

      Tmin = 1.0_dp
      Tmax = 300.0_dp
      TT=min(T,Tmax)
      TT=max(TT,Tmin)

      u  = (TT - T1) / (T2 - T1)
      u  = max(0.0_dp, min(1.0_dp, u))
      w  = u**3 * (6.0_dp*u**2 - 15.0_dp*u + 10.0_dp)
      Rb = -0.00630328_dp*TT**3 + 0.440052_dp*TT**2 - 1.719_dp*TT + 1.92937_dp
      Rh = 5.58686404_dp * TT**0.98980198_dp
      heat_capacity = (1.0_dp - w)*Rb + w*Rh
  end function heat_capacity_glass_epoxy  
  
  
  subroutine material_glass_kapton_glass_init(cfg,density,c_pt)  
      class(input_t), target, intent(in) :: cfg
      real(dp), intent(out) :: density
      type(c_ptr), intent(out) :: c_pt

      type(glass_kapton_glass_cfg_t), pointer :: glass_kapton_glass_ptr ! fortran pointer

      density=cfg%dbl('density',1800.0_dp)

      allocate(glass_kapton_glass_ptr)
      c_pt = c_loc(glass_kapton_glass_ptr) ! translates fortran pointer to C pointer
  end subroutine material_glass_kapton_glass_init  
    

  function thermal_conductivity_glass_kapton_glass(c_pt,T,B) result(thermal_conductivity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T, B
      real(dp) :: thermal_conductivity    
      real(dp) :: Tmin,Tmax,TT
      type(glass_kapton_glass_cfg_t), pointer :: glass_kapton_glass_ptr
      call c_f_pointer(c_pt,glass_kapton_glass_ptr)

      Tmin = 1.0_dp
      Tmax = 300.0_dp
      TT=min(T,Tmax)
      TT=max(TT,Tmin)

      thermal_conductivity=exp(-4.5664911_dp+0.7191848_dp*log(TT)-0.0128252_dp*(log(TT)**2))      
  end function thermal_conductivity_glass_kapton_glass  


  function heat_capacity_glass_kapton_glass(c_pt,T,B,TC,TCS,TC0) result(heat_capacity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B,TC,TCS,TC0
      real(dp) :: heat_capacity    
      real(dp) :: Tmin,Tmax,TT,u,w,Rb,Rh
      real(dp), parameter :: T1 = 20.0_dp, T2 = 30.0_dp
      type(glass_kapton_glass_cfg_t), pointer :: glass_kapton_glass_ptr
      call c_f_pointer(c_pt,glass_kapton_glass_ptr)

      Tmin = 1.0_dp
      Tmax = 300.0_dp
      TT=min(T,Tmax)
      TT=max(TT,Tmin)

      u  = (TT - T1) / (T2 - T1)
      u  = max(0.0_dp, min(1.0_dp, u))
      w  = u**3 * (6.0_dp*u**2 - 15.0_dp*u + 10.0_dp)
      Rb = -0.00630328_dp*TT**3 + 0.440052_dp*TT**2 - 1.719_dp*TT + 1.92937_dp
      Rh = 5.58686404_dp * TT**0.98980198_dp
      heat_capacity = (1.0_dp - w)*Rb + w*Rh
  end function heat_capacity_glass_kapton_glass   

end module lib_material_insulation_m


