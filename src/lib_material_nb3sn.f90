! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module lib_material_nb3sn_m
  use krn_global_tools_m
  use iso_c_binding
  use lib_input_m, only: input_t

  implicit none

  type nb3sn_cfg_t
    real(dp) :: E0,nPow,str1,str2,Bc20m,Tc0m,nu,Ca1,Ca2,e0a,emax,C0,p,q
  end type nb3sn_cfg_t  
  
  contains

 
  subroutine material_nb3sn_init(cfg,density,c_pt)  
      class(input_t), target, intent(in) :: cfg
      real(dp), intent(out) :: density
      type(c_ptr), intent(out) :: c_pt

      type(nb3sn_cfg_t), pointer :: nb3sn_ptr ! fortran pointer

      density=cfg%dbl('density',8040.0_dp)

      allocate(nb3sn_ptr)
      nb3sn_ptr%E0    = cfg%dbl('E0',1.0e-5_dp)
      nb3sn_ptr%nPow  = cfg%int('nPow',5)  
      nb3sn_ptr%str1  = cfg%dbl('str1',0.0060942_dp)
      nb3sn_ptr%str2  = cfg%dbl('str2',1.0777e-9_dp)      
      nb3sn_ptr%Bc20m = cfg%dbl('Bc20m',29.39_dp)
      nb3sn_ptr%Tc0m  = cfg%dbl('Tc0m',16.48_dp)
      nb3sn_ptr%nu    = cfg%dbl('nu',1.52_dp)
      nb3sn_ptr%Ca1   = cfg%dbl('Ca1',45.74_dp)
      nb3sn_ptr%Ca2   = cfg%dbl('Ca2',4.431_dp)
      nb3sn_ptr%e0a   = cfg%dbl('e0a',0.00232_dp)
      nb3sn_ptr%emax  = cfg%dbl('emax',0.00_dp)      
      nb3sn_ptr%C0    = cfg%dbl('C0',8.0771e10_dp)
      nb3sn_ptr%p     = cfg%dbl('p',0.556_dp)
      nb3sn_ptr%q     = cfg%dbl('q',1.698_dp)   

      c_pt = c_loc(nb3sn_ptr) ! translates fortran pointer to C pointer
  end subroutine material_nb3sn_init    


  subroutine get_power_law_parameters_nb3sn(c_pt,E0,nPow) 
      type(c_ptr), intent(in) :: c_pt
      real(dp), intent(out) :: E0
      integer, intent(out) :: nPow

      type(nb3sn_cfg_t), pointer :: nb3sn_ptr
      call c_f_pointer(c_pt,nb3sn_ptr)
      E0 = nb3sn_ptr%E0
      nPow = nb3sn_ptr%nPow
  end subroutine get_power_law_parameters_nb3sn


  function thermal_conductivity_nb3sn(c_pt,T,B) result(thermal_conductivity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B
      real(dp) :: thermal_conductivity    
      real(dp) :: tmin,tmax,TT,T0,AA1,BB1,n1,m1,AA2,BB2,n2,m2  
      type(nb3sn_cfg_t), pointer :: nb3sn_ptr
      call c_f_pointer(c_pt,nb3sn_ptr)

      tmin=1.0d0
      tmax=1000.0d0
    
      TT=min(T,tmax)
      TT=max(TT,tmin)

      AA1=1.39049E+14_dp  ;   BB1=76.61272374_dp
      n1=3.663439606_dp   ;   m1=9.322084303_dp
      AA2=35.77281992_dp  ;   BB2=5.967307791_dp
      n2=2.766424235_dp   ;   m2=3.335767669_dp
      T0=20.92684715_dp
    
      if(TT <= T0) then
         thermal_conductivity = AA1*TT**n1/(BB1+TT)**m1
      else
         thermal_conductivity = AA2*TT**n2/(BB2+TT)**m2
      endif     

  end function thermal_conductivity_nb3sn


  function heat_capacity_nb3sn(c_pt,T,B,TC,TCS,TC0) result(heat_capacity) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B,TC,TCS,TC0
      real(dp) :: heat_capacity
      real(dp) :: CP,CPN,CPS,AAA,BBB,CCC,DDD,F,CPN2
      real(dp) :: AA,BB,CC,DD,a,bbbb,c,d,na,nb,nc,nd,TMAX,TT

      type(nb3sn_cfg_t), pointer :: nb3sn_ptr
      call c_f_pointer(c_pt,nb3sn_ptr)

      AA = 38.2226876_dp  ; BB = -848.36422_dp
      CC = 1415.13807_dp  ; DD = -346.83796_dp
      a  = 6.804586085_dp ; bbbb  = 59.92091818_dp
      c  = 25.82863336_dp ; d  = 8.779183354_dp
      na = 1.0_dp         ; nb = 2.0_dp
      nc = 3.0_dp         ; nd = 4.0_dp
      TMAX = 400.0_dp

      if(T<=TC) then
        TT=TC/TC0
        CCC=(-0.46306_dp) - (0.067830_dp)*TC
        DDD=27.2_dp/(1.0_dp+(0.34_dp*TT))**2
        AAA=1500.0_dp*(CCC**2)/(2.0_dp*DDD-1.0_dp)
        if (TC<=10.0_dp) then
           BBB=(7.5475E-3_dp)*TC**2
        else if(TC>10.0_dp .and. TC<=20.0_dp) then
           BBB=(-0.3_dp+0.00375_dp*TC**2)/0.09937_dp
        endif
        CPS=(AAA+BBB)*(T/TC)**3
      else
        CPS=0.0_dp
      endif

      if (T<=10.0_dp .and. Tc<=10.0_dp) then
        CPN=(7.5475e-3_dp)*T**2
      else if (T>10.0_dp .and. T<=20.0_dp .and. Tc>10.0_dp .and. Tc<=20.0_dp) then
         CPN=(-0.3_dp + 0.00375_dp*T**2)/(0.09937_dp)
      else
         TT=T
         TT=min(TT,TMAX)
         CPN = AA*TT /(a+TT)**na + BB*TT**2/(bbbb+TT)**nb + CC*TT**3/(c+TT)**nc + DD*TT**4/(d+TT)**nd
      endif

      if (T<=10.0_dp) then
        CPN2=(7.5475e-3_dp)*T**2
      else if (T>10.0_dp .and. T<=20.0_dp) then
         CPN2=(-0.3_dp + 0.00375_dp*T**2)/(0.09937_dp)
      else
         TT=T
         TT=min(TT,TMAX)
         CPN2 = AA*TT /(a+TT)**na + BB*TT**2/(bbbb+TT)**nb + CC*TT**3/(c+TT)**nc + DD*TT**4/(d+TT)**nd
      endif   

      if(T<=TCS) then
        CP=CPS
      else if (T>TCS .and. T<=TC) then
        if(TCS<TC) then
           F= (T-TCS)/(TC-TCS)
        else
           F= 1.0_dp
        endif
        CP= F*CPN + (1.0_dp-F)*CPS
      else if (T>TC) then
        CP=CPN2
      endif
      heat_capacity=CP      
  end function heat_capacity_nb3sn

  
  function strain_nb3sn(c_pt,B,I) result(strain) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: B,I
      real(dp) :: strain

      type(nb3sn_cfg_t), pointer :: nb3sn_ptr
      call c_f_pointer(c_pt,nb3sn_ptr)

      strain = -nb3sn_ptr%str1-nb3sn_ptr%str2*I*B
  end function strain_nb3sn


  real(dp) function sNb3Sn(e,Ca1,Ca2,e0a)
      real(dp), intent(in) :: e,Ca1,Ca2,e0a
      real(dp) :: esh,s,esh2,e0a2,ed,ed2

      esh  = ca2 * e0a / sqrt(ca1*ca1+ca2*ca2)
      ed   = e-esh
      esh2 = esh*esh
      e0a2 = e0a*e0a
      ed2  = ed*ed

      s = 1.0_dp+(ca1*(sqrt(esh2+e0a2)-sqrt(ed2+e0a2))-Ca2*e)/(1.0_dp-Ca1*e0a)
      if(s<0.0) s = 0.0_dp
      sNb3Sn = s
      return
  end function sNb3Sn    


  real(dp) function hNb3Sn(t,nu)
    real(dp), intent(in) :: t,nu
    real(dp) :: h

    if(t>0.0_dp .and. t<1.0_dp) then
       h = (1.0_dp-t**nu) * (1.0_dp-t*t)
    else if(t<=0.0_dp) then
       h = 1.0_dp
    else if(t>=1.0_dp) then
       h = 0.0_dp
    endif

    hNb3Sn = h
    return
  end function hNb3Sn     


  real(dp) function fpNb3Sn(b,p,q)
    real(dp), intent(in) :: b,p,q
    real(dp) :: fp

    if(b>0.0_dp .and. b<1.0_dp) then
       fp = b**p * (1.0_dp-b)**q
    else
       fp = 0.0_dp
    endif

    fpNb3Sn = fp
    return
  end function fpNb3Sn      
  

  function critical_temperature_nb3sn(c_pt,B,St) result(critical_temperature) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: B,St
      real(dp) :: critical_temperature
      real(dp) :: Blim,ee,s,Bc20,bb0,tt,Tc,Blow

      type(nb3sn_cfg_t), pointer :: nb3sn_ptr
      call c_f_pointer(c_pt,nb3sn_ptr)

      ! L. Bottura, B. Bordini, Jc(B,T,e) Parameterization for the ITER
      ! Nb3Sn Production, IEEE Trans. Appl. Sup., 19(2), 1477-1480, 2009      

      Blow=1.0e-3_dp
      Blim=max(abs(B),Blow)

      ee = St-nb3sn_ptr%emax
      s  = sNb3Sn(ee,nb3sn_ptr%Ca1,nb3sn_ptr%Ca2,nb3sn_ptr%e0a)
      Bc20 = nb3sn_ptr%Bc20m * s
      bb0 = Blim/Bc20

      if(bb0>0.0_dp .and. bb0<1.0_dp) then
         tt = s**0.3333333_dp * (1.0_dp-bb0)**(1.0_dp/nb3sn_ptr%nu)
      else if(abs(bb0)<=1.0e-12_dp) then
         tt = 1.0_dp
      else if(bb0>1.0_dp) then
         tt = 0.0_dp
      endif

      Tc = nb3sn_ptr%Tc0m * tt
      critical_temperature = Tc
  end function critical_temperature_nb3sn


  function critical_field_nb3sn(c_pt,T,St) result(critical_field) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,St
      real(dp) :: critical_field
      real(dp) :: Tlim,ee,s,Tc0,tt,bb,Tlow

      type(nb3sn_cfg_t), pointer :: nb3sn_ptr
      call c_f_pointer(c_pt,nb3sn_ptr)

      ! L. Bottura, B. Bordini, Jc(B,T,e) Parameterization for the ITER
      ! Nb3Sn Production, IEEE Trans. Appl. Sup., 19(2), 1477-1480, 2009      

      Tlow=0.0_dp
      Tlim=max(T,Tlow)

      ee = St-nb3sn_ptr%emax
      s  = sNb3Sn(ee,nb3sn_ptr%Ca1,nb3sn_ptr%Ca2,nb3sn_ptr%e0a)

      Tc0 = nb3sn_ptr%Tc0m * s**0.333333_dp
      tt = Tlim/Tc0
 
      if(tt>0.0_dp .and. tt<1.0_dp) then
         bb = s * (1.0_dp-tt**nb3sn_ptr%nu)
      else if(tt<=0.0_dp) then
         bb = 1.0_dp
      else if(tt>1.0_dp) then
         bb = 0.0_dp
      endif

      critical_field = nb3sn_ptr%Bc20m * bb
  end function critical_field_nb3sn


  function critical_current_density_nb3sn(c_pt,T,B,St,Tc0,Bc) result(critical_current_density) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: T,B,St,Tc0,Bc
      real(dp) :: critical_current_density
      real(dp) :: ee,s,tt,bb,h,fp,Blim,Tlim,Blow,Tlow

      type(nb3sn_cfg_t), pointer :: nb3sn_ptr
      call c_f_pointer(c_pt,nb3sn_ptr)

      ! L. Bottura, B. Bordini, Jc(B,T,e) Parameterization for the ITER
      ! Nb3Sn Production, IEEE Trans. Appl. Sup., 19(2), 1477-1480, 2009      

      Blow=1.0e-3_dp
      Tlow=0.0_dp
      Blim=max(abs(B),Blow)
      Tlim=max(T,Tlow)
      ee = St-nb3sn_ptr%emax
      s  = sNb3Sn(ee,nb3sn_ptr%Ca1,nb3sn_ptr%Ca2,nb3sn_ptr%e0a)

      tt = Tlim/Tc0
      if(tt>1.0_dp) then
         critical_current_density = 0.0_dp
         return
      endif

      bb = Blim/Bc
      if(bb>1.0_dp) then
         critical_current_density = 0.0_dp
         return
      endif

      h = hNb3Sn(tt,nb3sn_ptr%nu)
      fp = fpNb3Sn(bb,nb3sn_ptr%p,nb3sn_ptr%q)
      critical_current_density = nb3sn_ptr%C0/Blim * s * h * fp
  end function critical_current_density_nb3sn


  function current_sharing_temperature_nb3sn(c_pt,B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc) result(current_sharing_temperature) bind(C)
      type(c_ptr), value :: c_pt
      real(dp), value :: B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc
      real(dp) :: current_sharing_temperature
      real(dp) :: T,tt,ttlow,ttup,error,tolerance,Jc
      logical :: converged      

      type(nb3sn_cfg_t), pointer :: nb3sn_ptr
      call c_f_pointer(c_pt,nb3sn_ptr)

      tolerance=1.0e-5_dp
    
      if(B>=Bc0) then
        current_sharing_temperature=0.0_dp
        return
      endif

      if(Jop>=Jc0) then
         current_sharing_temperature=0.0_dp
         return
      endif

      if(Jop<=0.0_dp) then
         current_sharing_temperature = Tc
         return
      endif

      ttup =Tc/nb3sn_ptr%Tc0m
      ttlow=0.0_dp

      converged=.false.
      do while(.not.converged)
        tt=0.5_dp*(ttlow+ttup)
        T  = tt * nb3sn_ptr%Tc0m
        Jc = critical_current_density_nb3sn(c_pt,T,B,St,Tc0,Bc)
        if(Jc>Jop) then
           ttlow = tt
        else if(Jc<=Jop) then
           ttup  = tt
        else if(abs(Jc-Jop)<=1.0e-12_dp) then
           ttup  = tt
           ttlow = tt
        endif
        error     = abs(ttup-ttlow)
        converged = error<=tolerance
      enddo

      current_sharing_temperature = tt*nb3sn_ptr%Tc0m
  end function current_sharing_temperature_nb3sn 

end module lib_material_nb3sn_m