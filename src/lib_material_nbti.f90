! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module lib_material_nbti_m
  use krn_global_tools_m
  use iso_c_binding
  use lib_input_m, only: input_t

  implicit none

  type nbti_cfg_t
    real(dp) :: E0,nPow,str1,str2,Bc20,Tc0_p,nu,CC0,p,q,n
  end type nbti_cfg_t
  
  contains


  subroutine material_nbti_init(cfg,density,c_pt)  
      class(input_t), target, intent(in) :: cfg
      real(dp), intent(out) :: density
      type(c_ptr), intent(out) :: c_pt

      type(nbti_cfg_t), pointer :: nbti_ptr ! fortran pointer

      density=cfg%dbl('density',6000.0_dp)

      allocate(nbti_ptr)
      nbti_ptr%E0    = cfg%dbl('E0',1.0e-5_dp)
      nbti_ptr%nPow  = cfg%int('nPow',5)      
      nbti_ptr%str1  = cfg%dbl('str1',0.00742_dp)
      nbti_ptr%str2  = cfg%dbl('str2',1.301e-9_dp)
      nbti_ptr%Bc20  = cfg%dbl('Bc20',13.72_dp)
      nbti_ptr%Tc0_p = cfg%dbl('Tc0_p',8.79_dp)
      nbti_ptr%nu    = cfg%dbl('nu',1.7_dp)      
      nbti_ptr%CC0   = cfg%dbl('CC0',8.92534e11_dp)
      nbti_ptr%p     = cfg%dbl('p',0.98_dp)
      nbti_ptr%q     = cfg%dbl('q',0.98_dp)
      nbti_ptr%n     = cfg%dbl('n',1.96_dp)

      c_pt = c_loc(nbti_ptr) ! translates fortran pointer to C pointer
  end subroutine material_nbti_init      


   subroutine get_power_law_parameters_nbti(c_pt,E0,nPow) 
      type(c_ptr), intent(in) :: c_pt
      real(dp), intent(out) :: E0
      integer, intent(out) :: nPow

      type(nbti_cfg_t), pointer :: nbti_ptr
      call c_f_pointer(c_pt,nbti_ptr)
      E0 = nbti_ptr%E0
      nPow = nbti_ptr%nPow
  end subroutine get_power_law_parameters_nbti


  function thermal_conductivity_nbti(c_pt,T,B) result(thermal_conductivity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B
      real(dp) :: thermal_conductivity
      real(dp) :: tmin,tmax,TT,a0,a1,a2,a3,a4,a5,a6
      real(dp) :: T_CUT,k_cut,slope_cut
      type(nbti_cfg_t), pointer :: nbti_ptr
      call c_f_pointer(c_pt,nbti_ptr)

      tmin=1.0d0
      tmax=1000.0d0
    
      TT=min(T,tmax)
      TT=max(TT,tmin)

      a0=6.60e-2_dp
      a1=4.56e-2_dp
      a2=3.00e-4_dp
      a3=-3.00e-6_dp
      a4=6.00e-9_dp
      a5=1.5e-11_dp
      a6=-5.0e-14_dp

      T_CUT     = 200.0_dp
      k_cut     = a0 + a1*T_CUT + a2*T_CUT**2 + a3*T_CUT**3 &
                + a4*T_CUT**4 + a5*T_CUT**5 + a6*T_CUT**6
      slope_cut = a1 + 2.0_dp*a2*T_CUT + 3.0_dp*a3*T_CUT**2 + 4.0_dp*a4*T_CUT**3 &
                + 5.0_dp*a5*T_CUT**4 + 6.0_dp*a6*T_CUT**5

      if (TT <= T_CUT) then
          thermal_conductivity=a0+a1*TT+a2*TT**2+a3*TT**3+a4*TT**4+a5*TT**5+a6*TT**6
      else
          thermal_conductivity = k_cut + slope_cut * (TT - T_CUT)
      end if
  end function thermal_conductivity_nbti  


  function heat_capacity_nbti(c_pt,T,B,TC,TCS,TC0) result(heat_capacity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B,TC,TCS,TC0
      real(dp) :: heat_capacity
      real(dp) :: cpHigh, Bc20_loc, gama, Beta, CpLow  

      type(nbti_cfg_t), pointer :: nbti_ptr
      call c_f_pointer(c_pt,nbti_ptr)

      cpHigh=400.0_dp

      Bc20_loc=14.0_dp
      gama=0.145_dp
      Beta=0.0023_dp

      if(T<Tc) then
          CpLow=(Beta+3.0_dp*gama/Tc**2)*T**3+gama*B*T/Bc20_loc
      else
          CpLow=Beta*T**3+gama*T
      endif

      heat_capacity=1.0_dp/(1.0_dp/cpHigh+1.0_dp/CpLow)          
  end function heat_capacity_nbti


  function strain_nbti(c_pt,B,I) result(strain) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: B,I
      real(dp) :: strain

      type(nbti_cfg_t), pointer :: nbti_ptr
      call c_f_pointer(c_pt,nbti_ptr)

      strain = -nbti_ptr%str1-nbti_ptr%str2*I*B
  end function strain_nbti


  function critical_temperature_nbti(c_pt,B,St) result(critical_temperature) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: B,St
      real(dp) :: critical_temperature
      real(dp) :: Blim,bb0,tt,Tc,Blow

      type(nbti_cfg_t), pointer :: nbti_ptr
      call c_f_pointer(c_pt,nbti_ptr)

      ! M.S.Lubell, Scaling formulas for critical current and critical field
      ! for commercial NbTi, IEEE Trans. Mag. ,19, (1983).      

      Blow=1.0e-3_dp
      Blim=max(abs(B),Blow)
      bb0 = Blim/nbti_ptr%Bc20

      if(bb0>0.0_dp .and. bb0<1.0_dp) then
         tt = (1.0_dp-bb0)**(1.0_dp/nbti_ptr%nu)
      else if(abs(bb0)<=1.0e-12_dp) then
         tt = 1.0_dp
      else if(bb0>1.0_dp) then
         tt = 0.0_dp
      endif
      Tc = nbti_ptr%Tc0_p * tt
      critical_temperature = Tc
  end function critical_temperature_nbti


  function critical_field_nbti(c_pt,T,St) result(critical_field) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,St
      real(dp) :: critical_field
      real(dp) :: Tlim,tt,bb,Tlow

      type(nbti_cfg_t), pointer :: nbti_ptr
      call c_f_pointer(c_pt,nbti_ptr)

      ! M.S.Lubell, Scaling formulas for critical current and critical field
      ! for commercial NbTi, IEEE Trans. Mag. ,19, (1983).      

      Tlow=0.0_dp
      Tlim=max(T,Tlow)
      tt = Tlim/nbti_ptr%Tc0_p

      if(tt>0.0_dp .and. tt<1.0_dp) then
        bb = 1.0_dp-tt**nbti_ptr%nu
      else if(tt<=0.0_dp) then
        bb = 1.0_dp
      else if(tt>1.0_dp) then
        bb = 0.0_dp
      endif
 
      critical_field = nbti_ptr%Bc20 * bb
  end function critical_field_nbti


  real(dp) function hNbTi(t,n)
    real(dp), intent(in) :: t,n
    real(dp) :: h

    if(t>0.0_dp .and. t<1.0_dp) then
       h = 1.0_dp-t**n
    else if(t<=0.0_dp) then
       h = 1.0_dp
    else if(t>=1.0_dp) then
       h = 0.0_dp
    endif

    hNbTi = h
    return
  end function hNbTi


  real(dp) function fpNbTi(b,p,q)
    real(dp), intent(in) :: b,p,q
    real(dp) :: fp

    if(b>0.0_dp .and. b<1.0_dp) then
       fp = b**p * (1.0_dp-b)**q
    else
       fp = 0.0_dp
    endif

    fpNbTi=fp
    return
  end function fpNbTi    


  function critical_current_density_nbti(c_pt,T,B,St,Tc0,Bc) result(critical_current_density) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B,St,Tc0,Bc
      real(dp) :: critical_current_density
      real(dp) :: Blim,Tlim,tt,bb,h,fp,Blow,Tlow

      type(nbti_cfg_t), pointer :: nbti_ptr
      call c_f_pointer(c_pt,nbti_ptr)

      ! M.A. Green, Calculating the Jc, B, T Surface for Niobium Titanium
      ! Using a Reduced State Model, IEEE Trans. Mag., 25, 2, (1989).
  
      ! G. Morgan, A Comparison of Two Analytic Forms for the Jc(B,T)
      ! surface, SSC Magnet Division Notes, 310-1 (SSC-MD-218), (1989)
  
      ! L. Bottura, B. Bordini, Jc(B,T,e) Parameterization for the ITER
      ! Nb3Sn Production, IEEE Trans. Appl. Sup., 19(2), 1477-1480, 2009
              
      Blow=1.0e-3_dp
      Tlow=0.0_dp
      Blim=max(abs(B),Blow)
      Tlim=max(T,Tlow)
      tt = Tlim/nbti_ptr%Tc0_p

      if(tt>1.0_dp) then
         critical_current_density=0.0_dp
         return
      endif

      bb = Blim/Bc
      if(bb>1.0_dp) then
         critical_current_density=0.0_dp
         return
      endif

      h = hNbTi(tt,nbti_ptr%nu)
      fp = fpNbTi(bb,nbti_ptr%p,nbti_ptr%q)
      critical_current_density = nbti_ptr%CC0/Blim * h**nbti_ptr%n * fp
  end function critical_current_density_nbti
  

  function current_sharing_temperature_nbti(c_pt,B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc) result(current_sharing_temperature) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc
      real(dp) :: current_sharing_temperature
      real(dp) :: T,tt,ttlow,ttup,error,tolerance,Jc
      logical :: converged

      type(nbti_cfg_t), pointer :: nbti_ptr
      call c_f_pointer(c_pt,nbti_ptr)

      tolerance=1.0e-5_dp

      if(B>=nbti_ptr%Bc20) then
        current_sharing_temperature=0.0_dp
        return
      endif

      if(Jop>=Jc0) then
        current_sharing_temperature=0.0_dp
        return
      endif

      if(Jop<=0.0_dp) then
        current_sharing_temperature = nbti_ptr%Tc0_p
        return
      endif

      ttup =1.0_dp
      ttlow=0.0_dp

      converged=.false.
      do while(.not.converged) 
        tt=0.5_dp*(ttlow+ttup)
        T  = tt * nbti_ptr%Tc0_p
        Jc = critical_current_density_nbti(c_pt,T,B,St,Tc0,Bc)

        if(Jc>Jop) then
           ttlow = tt
        elseif(Jc<=Jop) then
           ttup  = tt
        elseif(abs(Jc-Jop)<=1.0e-12_dp) then
           ttup  = tt
           ttlow = tt
        endif

        error     = abs(ttup-ttlow)
        converged = error<=tolerance
      enddo
      current_sharing_temperature = tt*nbti_ptr%Tc0_p
  end function current_sharing_temperature_nbti
    
end module lib_material_nbti_m


