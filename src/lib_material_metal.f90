! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module lib_material_metal_m
  use krn_global_tools_m
  use iso_c_binding
  use lib_input_m, only: input_t

  implicit none

  type copper_cfg_t
    real(dp) :: RRR
  end type copper_cfg_t   

  type stainless_steel_cfg_t
  end type stainless_steel_cfg_t  
  
  contains


  subroutine material_copper_init(cfg,density,c_pt)  
      class(input_t), target, intent(in) :: cfg
      real(dp), intent(out) :: density
      type(c_ptr), intent(out) :: c_pt

      type(copper_cfg_t), pointer :: copper_ptr ! fortran pointer

      density=cfg%dbl('density',8960.0_dp)

      allocate(copper_ptr)
      copper_ptr%RRR=cfg%dbl('RRR',100.0_dp)      

      c_pt = c_loc(copper_ptr) ! translates fortran pointer to C pointer
  end subroutine material_copper_init


  function resistivity_copper(c_pt,T,B) result(resistivity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B
      real(dp) :: resistivity    
      real(dp) :: Tmin,Tmax,TT,arg,rho1,A,aa,bb
      type(copper_cfg_t), pointer :: copper_ptr
      call c_f_pointer(c_pt,copper_ptr)

      Tmin   = 0.1_dp ; Tmax   = 1000.0_dp
      TT=T
      TT=max(TT,Tmin)
      TT=min(TT,Tmax)

      arg     = (50.0_dp/TT)**6.428_dp
      rho1    = 1.171e-17_dp*TT**4.49_dp/(1.0_dp+4.5e-7_dp*TT**3.35_dp*exp(-arg))

      bb=1.69e-8_dp
      resistivity = bb/copper_ptr%RRR + rho1 + 0.4531_dp*(bb*rho1/(copper_ptr%RRR*rho1+bb))

      if(abs(B)<1.0e-15_dp) return

      A = log10(1.553e-8_dp*B/resistivity)
      aa = -2.662_dp + 0.3168_dp*A + 0.6229_dp*A**2 - 0.1839_dp*A**3 + 0.01827_dp*A**4

      resistivity = resistivity*(1.0_dp+10.0_dp**aa)   
  end function resistivity_copper


  function thermal_conductivity_copper(c_pt,T,B) result(thermal_conductivity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B
      real(dp) :: thermal_conductivity    
      real(dp) :: Tmin,Tmax,TT,BETA,BETAR,W1,W10,w0,p1,p2,p3,p4,p5,p6,p7,alpha,L0,BB      
      type(copper_cfg_t), pointer :: copper_ptr
      call c_f_pointer(c_pt,copper_ptr)

      Tmin=1.0_dp;    Tmax=3.00e2_dp
      TT=T
      TT=min(TT,Tmax)
      TT=max(TT,Tmin)

      BETA = 0.634_dp / copper_ptr%RRR
      BETAR = BETA / (0.0003_dp)
      p1 = 0.00000001754_dp
      p2 = 2.763_dp
      p3 = 1102_dp
      p4 = -0.165_dp
      p5 = 70_dp
      p6 = 1.756_dp
      p7 = 0.838_dp / BETAR**0.1661_dp
      alpha=5e-11_dp ;  L0=2.44e-8_dp

      W0 = BETA / TT
      W0 = W0+alpha*B/TT/L0
      W1 = p1 * TT**p2 / (1 + p1 * p3 * TT**(p2 + p4) * Exp(-(p5 / TT)**p6))
      W10 = p7 * W0 * W1 / (W0 + W1)
      thermal_conductivity = 1.0_dp / (W0 + W1 + W10)
  end function thermal_conductivity_copper  


  function heat_capacity_copper(c_pt,T,B,TC,TCS,TC0) result(heat_capacity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B,TC,TCS,TC0
      real(dp) :: heat_capacity    
      real(dp) :: Tmin,Tmax,TT,Cp300,Cplow
      real(dp) :: beta,gamma

      type(copper_cfg_t), pointer :: copper_ptr
      call c_f_pointer(c_pt,copper_ptr)

      gamma=1.1e-2_dp
      beta=1.1e-3_dp
      cp300=385.49_dp
      Tmin=1.0_dp
      Tmax=3.00e2_dp
  
      TT=T
      TT=min(TT,Tmax)
      TT=max(TT,Tmin)
  
      Cplow=gamma*TT+beta*TT**3
      heat_capacity=1_dp/(1_dp/Cp300+1_dp/Cplow)
  end function heat_capacity_copper


  subroutine material_stainless_steel_init(cfg,density,c_pt)  
      class(input_t), target, intent(in) :: cfg
      real(dp), intent(out) :: density
      type(c_ptr), intent(out) :: c_pt

      type(stainless_steel_cfg_t), pointer :: stainless_steel_ptr ! fortran pointer

      density=cfg%dbl('density',7900.0_dp)

      allocate(stainless_steel_ptr)
      c_pt = c_loc(stainless_steel_ptr) ! translates fortran pointer to C pointer
  end subroutine material_stainless_steel_init  
    

  function thermal_conductivity_stainless_steel(c_pt,T,B) result(thermal_conductivity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B
      real(dp) :: thermal_conductivity    
      real(dp) :: tmin,tmax,TT,a0,a1,a2,a3,a4,a5
      type(stainless_steel_cfg_t), pointer :: stainless_steel_ptr
      call c_f_pointer(c_pt,stainless_steel_ptr)

      tmin=1.0d0
      tmax=1000.0d0
    
      TT=min(T,tmax)
      TT=max(TT,tmin)

      a0=-2.7522797_dp
      a1=1.47636_dp
      a2=-0.5109092_dp
      a3=0.2846192_dp
      a4=-0.063919_dp
      a5=0.0047142_dp

      thermal_conductivity = exp(a0 + a1*log(TT) + a2*log(TT)**2 + a3*log(TT)**3 + a4*log(TT)**4 + a5*log(TT)**5)
  end function thermal_conductivity_stainless_steel  


  function heat_capacity_stainless_steel(c_pt,T,B,TC,TCS,TC0) result(heat_capacity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B,TC,TCS,TC0
      real(dp) :: heat_capacity    
      real(dp) :: tmin,tmax,TT,CpHigh,Cplow,gamma,beta
      type(stainless_steel_cfg_t), pointer :: stainless_steel_ptr
      call c_f_pointer(c_pt,stainless_steel_ptr)

      gamma=0.48_dp
      beta=0.00075_dp
      CpHigh=500.0_dp
      tmin=1.0_dp
      tmax=3.00e2_dp
  
      TT=T
      TT=min(TT,tmax)
      TT=max(TT,tmin)
  
      Cplow=gamma*TT+beta*TT**3
      heat_capacity=1_dp/(1.0_dp/CpHigh+1.0_dp/Cplow)
  end function heat_capacity_stainless_steel

end module lib_material_metal_m


