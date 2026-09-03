! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_boundary_calc_m
    use cmp_boundary_init_m
  use lib_He_thermo_m
    use lib_ext_math_m, only: brent_t, zero
    use ieee_arithmetic, only: ieee_is_nan
    implicit none

    type :: hug_star_t
        real(dp) :: vstar, ustar, W, K
        real(dp) :: ev_v
        real(dp) :: dustar_dPstar
        real(dp) :: dvstar_dPstar, dvstar_dv_LR, dvstar_dp_LR
        logical  :: ok
        integer  :: n_iter
    end type hug_star_t

    type, extends(brent_t) :: hug_brent_t
        real(dp) :: Pstar, p_LR, v_LR, e_LR
    contains
        procedure :: f => hug_brent_f
    end type hug_brent_t

contains


subroutine jacobian_star(Cons_star,Jacob_star)
  real(dp), dimension(:), intent(in) :: Cons_star
  real(dp), dimension(:,:), intent(out) :: Jacob_star

  real(dp), dimension(Nb_VarC) :: U

  U(:)=Cons_star(:)

  Jacob_star(1,1)=0.0_dp
  Jacob_star(2,1)=1.0_dp
  Jacob_star(3,1)=0.0_dp
  Jacob_star(4,1)=0.0_dp

  Jacob_star(1,2)=-(3.0_dp/2.0_dp)*((U(2)**2)/(U(1)**2))
  Jacob_star(2,2)=3.0_dp*U(2)/U(1)
  Jacob_star(3,2)=-1.0_dp
  Jacob_star(4,2)=1.0_dp

  Jacob_star(1,3)=-U(4)*U(2)/(U(1)**2)-((U(2)**3)/(U(1)**3))
  Jacob_star(2,3)=U(4)/U(1)+(3.0_dp/2.0_dp)*((U(2)**2)/(U(1)**2))
  Jacob_star(3,3)=0.0_dp
  Jacob_star(4,3)=U(2)/U(1)
  
  Jacob_star(1,4)=-U(4)*U(2)/(U(1)**2)
  Jacob_star(2,4)=U(4)/U(1)
  Jacob_star(3,4)=0.0_dp
  Jacob_star(4,4)=U(2)/U(1)
end subroutine jacobian_star          


subroutine jacobian_star_Req_rem(PrimStar,Jacob_star)
  real(dp), dimension(:), intent(in) :: PrimStar
  real(dp), dimension(:,:), intent(out) :: Jacob_star

  real(dp) :: HH,kk,KKK,dPde,dPdT,dedT

  HH=PrimStar(Pri_e)+0.5_dp*PrimStar(Pri_u)**2+PrimStar(Pri_p)/PrimStar(Pri_ro)
  call jacobian_roT(PrimStar(Pri_ro), PrimStar(Pri_T), dedT, dPdT)
  dPde=dPdT/dedT
  kk=dPde/PrimStar(Pri_ro)
  KKK=PrimStar(Pri_c)**2+kk*(PrimStar(Pri_u)**2-HH)

  Jacob_star(1,1)=0.0_dp
  Jacob_star(2,1)=1.0_dp
  Jacob_star(3,1)=0.0_dp
  Jacob_star(4,1)=0.0_dp

  Jacob_star(1,2)=KKK-PrimStar(Pri_u)**2
  Jacob_star(2,2)=PrimStar(Pri_u)*(2.0_dp-kk)
  Jacob_star(3,2)=kk
  Jacob_star(4,2)=0.0_dp

  Jacob_star(1,3)=(KKK-HH)*PrimStar(Pri_u)
  Jacob_star(2,3)=HH-kk*PrimStar(Pri_u)**2
  Jacob_star(3,3)=PrimStar(Pri_u)*(1.0_dp+kk)
  Jacob_star(4,3)=0.0_dp
  
  Jacob_star(1,4)=0.0_dp
  Jacob_star(2,4)=0.0_dp
  Jacob_star(3,4)=0.0_dp
  Jacob_star(4,4)=0.0_dp
end subroutine jacobian_star_Req_rem     



function function_of_pstar_PT(sgnBC,Pstar,Rotank,Ctank,Ptank,Einttank,p_LR,u_LR,z_LR)
        
  real(dp), intent(in) :: sgnBC,Pstar,Rotank,Ctank,Ptank,Einttank,p_LR,u_LR,z_LR
  real(dp) :: function_of_pstar_PT, ustar, rostar_LR, Tstar_LR, R_Corr_star_LR, estar_LR, Cstar_LR

  rostar_LR=Rotank+(Pstar-Ptank)/(Ctank**2)
  Ustar=u_LR+sgnBC*(Pstar-p_LR)/z_LR
  call state_roP(rostar_LR, Pstar, R_Corr_star_LR, estar_LR, Tstar_LR, Cstar_LR)

  function_of_pstar_PT=estar_LR+Pstar/rostar_LR+0.5_dp*(ustar**2)-Einttank-Ptank/Rotank      

end function function_of_pstar_PT



function derivative_function_of_pstar_PT(sgnBC,Pstar,Rotank,Ctank,Ptank,p_LR,u_LR,z_LR)
  
  real(dp), intent(in) :: sgnBC,Pstar,Rotank,Ctank,Ptank,p_LR,u_LR,z_LR
  real(dp) :: derivative_function_of_pstar_PT
  real(dp) :: Tstar_LR,rostar_LR,dTdp_Ro,dTdRo_p,dRodP,dUdP,Ustar,dRstardRo,dRstardT,dRdP,Rstar
  real(dp) :: dPdT_local, cv_local

  rostar_LR=Rotank+(Pstar-Ptank)/(Ctank**2)
  Tstar_LR=T_roP(rostar_LR,Pstar)
  call jacobian_roT(rostar_LR, Tstar_LR, cv_local, dPdT_local, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
  Ustar=u_LR+sgnBC*(Pstar-p_LR)/z_LR
  dRodP=1.0_dp/(Ctank**2)
  dUdP=sgnBC/z_LR
  dRdP=(dRstardRo+dRstardT*dTdRo_p)/(Ctank**2)+dRstardT*dTdp_Ro

  derivative_function_of_pstar_PT=(dRdP*rostar_LR-Rstar*dRodP)/(rostar_LR**2)+Ustar*dUdP

end function derivative_function_of_pstar_PT



subroutine solve_boundary_PT(me,ss)
  type(boundary_t), intent(inout) :: me
  real(dp), intent(out) :: ss
 
  real(dp) :: kk, KKK, dedT, dPdT, HH, dPde, sgnBC
  real(dp), dimension(Nb_VarP) :: PrimStar
  real(dp) :: Ptank, Ttank, Rotank, Ctank, Ro_LR, p_LR, u_LR, Einttank, c_LR, z_LR, R_Correction_Star
  real(dp) :: Pstar, Tstar, Ustar, Rostar, r_Corr_LR, estar, Cstar, Funct, DerFunct
  real(dp) :: DerhdPstar,dcdro,dcdp,BB,dTdp_Ro,dTdRo_p,dRstardRo,dRstardT,Rstar,dRdP,T_temp
  real(dp) :: Rout,Pout,eout,cout,dR_dRo,R_Phys_Star
  real(dp), dimension(Nb_VarP) :: Pr_cell
  real(dp), dimension(Nb_VarC) :: Cs_cell,Cons_star,DerBBdU_LR,DerUstardU_LR_Glob
  real(dp), dimension(Nb_VarC) :: DerPstardU_LR_Glob,DerRostardU_LR_Glob,DerRcorrstardU_LR_Glob,DerRphysstardU_LR_Glob
  real(dp), dimension(Nb_VarC) :: DerUstardU_LR,DerhdU_LR,dc,DerRostardU_LR,DerRcorrstardU_LR,DerRphystardU_LR
  real(dp), dimension(Nb_VarC,Nb_VarC) :: Jacob_star,DerConstardU_LR
  integer :: NbreMaxIte, ite
  logical :: Inflow, use_hugoniot
  type(hug_star_t) :: st
  real(dp), dimension(Nb_VarC) :: dvstar_dU, dPsi_dU
  real(dp) :: v_LR, DeltaP_hug, Deltav_hug

  st%ok = .false.
  Ptank = me%imposed_PorM%v0d()
  Ttank = me%imposed_T%v0d()
  Rotank = ro_pT(Ptank,Ttank)
  call state_roT(Rotank, Ttank, Rout, Pout, eout, cout)
  Ctank = cout
  Einttank = eout

  Pr_cell(:)=me%port%Prim(:)
  Cs_cell(:)=me%port%Cons(:)
  sgnBC=-me%port%Sgn4j  ! sgnBC=1 if 1-port, -1 if 2-port

  Ro_LR=Pr_cell(Pri_ro)
  p_LR=Pr_cell(Pri_p)
  u_LR=Pr_cell(Pri_u)
  c_LR=Pr_cell(Pri_c)
  z_LR=Ro_LR*c_LR
  r_Corr_LR=Pr_cell(Pri_R)
  use_hugoniot = me%hugoniot_bc .and. &
                 abs(Ptank - p_LR) >= 1.0e-7_dp * Ro_LR * c_LR**2

  ! Inflow/outflow criterion: Hugoniot-consistent when use_hugoniot, acoustic otherwise
  if (use_hugoniot) then
    call hugoniot_star_state(sgnBC, Ptank, Ro_LR, p_LR, u_LR, c_LR, st)
    if (.not. st%ok) then; call set_error('hugoniot criterion: solver failed'); return; end if
    Inflow = sgnBC * st%ustar > 0.0_dp
  else
    Inflow = Ptank > p_LR - sgnBC*z_LR*u_LR
  end if

  if (Inflow) then
    Pstar=0.5_dp*(Ptank+p_LR)
    NbreMaxIte=200

    ! Newton method for Pstar calculation
    do ite=1,NbreMaxIte

      if (use_hugoniot) then
        call hugoniot_star_state(sgnBC, Pstar, Ro_LR, p_LR, u_LR, c_LR, st)
        if (.not. st%ok) then; call set_error('hugoniot inflow: solver failed'); return; end if
        Rostar=Rotank+(Pstar-Ptank)/(Ctank**2)
        call state_roP(Rostar, Pstar, R_Phys_Star, estar, Tstar, Cstar)
        if (ieee_is_nan(Tstar)) then; call set_error('convergence failure in trop (hugoniot inflow)'); return; end if
        Funct=estar+Pstar/Rostar+0.5_dp*(st%ustar**2)-Einttank-Ptank/Rotank
        call jacobian_roT(Rostar, Tstar, dedT, dPdT, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
        dRdP=(dRstardRo+dRstardT*dTdRo_p)/(Ctank**2)+dRstardT*dTdp_Ro
        DerFunct=(dRdP*Rostar-Rstar/(Ctank**2))/(Rostar**2)+st%ustar*st%dustar_dPstar
      else
        Funct=function_of_pstar_PT(sgnBC,Pstar,Rotank,Ctank,Ptank,Einttank,p_LR,u_LR,z_LR)
        DerFunct=derivative_function_of_pstar_PT(sgnBC,Pstar,Rotank,Ctank,Ptank,p_LR,u_LR,z_LR)
      end if

      if(abs(Funct)>1.0e-8_dp) then
        Pstar=Pstar-Funct/DerFunct
        Pstar=max(Pstar, 1.0e-6_dp*p_LR)
      else
        if (use_hugoniot) then
          Ustar=st%ustar
        else
          Ustar=u_LR+sgnBC*(Pstar-p_LR)/z_LR
        end if
        Rostar=Rotank+(Pstar-Ptank)/(Ctank**2)
        call state_roP(Rostar, Pstar, R_Phys_Star, estar, Tstar, Cstar)
        R_Correction_Star=1.0_dp
        if(R_Correction) R_Correction_Star=R_Phys_Star
        exit
      endif

      if(ite==NbreMaxIte) then
        call set_error('convergence failure in solve_boundary_condition_pipe_left')
        return
      endif

    enddo
    ! End of Newton method

  else

  ! Outflow
    Pstar=Ptank
    if(R_Correction) then
      if (use_hugoniot) then
        Rostar = 1.0_dp / st%vstar
        Ustar  = st%ustar
      else
        Ustar=u_LR+sgnBC*(Pstar-p_LR)/z_LR
        Rostar=(Pstar-p_LR)/(c_LR**2)+Ro_LR
      end if
      R_Correction_Star=r_Corr_LR*Rostar/Ro_LR+(c_LR**2)*(Rostar-Ro_LR)
      call state_roP_withR(Rostar, Pstar, R_Correction_Star, estar, Tstar, Cstar)
    else
      if (use_hugoniot) then
        Rostar = 1.0_dp / st%vstar
        Ustar  = st%ustar
      else
        Ustar=u_LR+sgnBC*(Pstar-p_LR)/z_LR
        Rostar=(Pstar-p_LR)/(c_LR**2)+Ro_LR
      end if
      call state_roP(Rostar, Pstar, R_Phys_Star, estar, Tstar, Cstar)
      if (ieee_is_nan(Tstar)) then; call set_error('convergence failure in T_roP (boundary)'); return; end if
      R_Correction_Star=r_Corr_LR
    endif
  endif
   
  me%port%flx(Con_Mas)=Rostar*Ustar
  me%port%flx(Con_Qdm)=Rostar*Ustar*Ustar+Pstar
  me%port%flx(Con_Ene)=(Rostar*(estar+0.5_dp*Ustar*Ustar)+Pstar)*Ustar
  me%port%flx(Con_R)=0.0_dp
  if(R_Correction) me%port%flx(Con_R)=R_Correction_Star*Ustar

  me%port%vit=Ustar

  call dc2_roT(Ro_LR, Pr_cell(Pri_T), dcdro, dcdp)
  dcdro = dcdro/(2.0_dp*c_LR)
  dcdp  = dcdp /(2.0_dp*c_LR)

  dc(Con_Mas)=dcdro-0.5_dp*dcdp*(Cs_cell(Con_Qdm)**2)/(Cs_cell(Con_Mas)**2)
  dc(Con_Qdm)=dcdp*Cs_cell(Con_Qdm)/Cs_cell(Con_Mas)
  dc(Con_Ene)=-dcdp
  dc(Con_R)=dcdp      

  if(Inflow) then

    if(R_Correction) then

      Cons_star(Con_Mas)=Rostar
      Cons_star(Con_Qdm)=Rostar*Ustar
      Cons_star(Con_Ene)=Rostar*(estar+0.5_dp*Ustar*Ustar)
      Cons_star(Con_R)=R_Correction_Star
      call jacobian_star(Cons_star,Jacob_star)

      if (use_hugoniot) then
        DerhdPstar = DerFunct
      else
        DerhdPstar=derivative_function_of_pstar_PT(sgnBC,Pstar,Rotank,Ctank,Ptank,p_LR,u_LR,z_LR)
      end if

      if (use_hugoniot) then
        v_LR       = 1.0_dp / Ro_LR
        DeltaP_hug = Pstar - p_LR
        Deltav_hug = v_LR - st%vstar

        dvstar_dU(1) = st%dvstar_dv_LR*(-v_LR**2) + st%dvstar_dp_LR*(-u_LR**2/2.0_dp)
        dvstar_dU(2) = st%dvstar_dp_LR*u_LR
        dvstar_dU(3) = st%dvstar_dp_LR*(-1.0_dp)
        dvstar_dU(4) = st%dvstar_dp_LR*1.0_dp

        dPsi_dU(1) = ((-v_LR**2 - dvstar_dU(1))*DeltaP_hug + Deltav_hug*(-u_LR**2/2.0_dp)) / DeltaP_hug**2
        dPsi_dU(2) = (            -dvstar_dU(2) *DeltaP_hug + Deltav_hug*u_LR)              / DeltaP_hug**2
        dPsi_dU(3) = (            -dvstar_dU(3) *DeltaP_hug + Deltav_hug*(-1.0_dp))         / DeltaP_hug**2
        dPsi_dU(4) = (            -dvstar_dU(4) *DeltaP_hug + Deltav_hug*1.0_dp)            / DeltaP_hug**2

        DerUstardU_LR(1) = -u_LR*v_LR + sgnBC*st%K*(u_LR**2/2.0_dp) + sgnBC*DeltaP_hug*dPsi_dU(1)/(2.0_dp*st%K)
        DerUstardU_LR(2) =  v_LR      - sgnBC*st%K*u_LR              + sgnBC*DeltaP_hug*dPsi_dU(2)/(2.0_dp*st%K)
        DerUstardU_LR(3) =              sgnBC*st%K                    + sgnBC*DeltaP_hug*dPsi_dU(3)/(2.0_dp*st%K)
        DerUstardU_LR(4) =            - sgnBC*st%K                    + sgnBC*DeltaP_hug*dPsi_dU(4)/(2.0_dp*st%K)

        DerhdU_LR(:) = Ustar * DerUstardU_LR(:)
        DerPstardU_LR_Glob(:) = -DerhdU_LR(:)/DerhdPstar
        DerRostardU_LR_Glob(:) = DerPstardU_LR_Glob(:)/(Ctank**2)
        DerUstardU_LR_Glob(:) = DerUstardU_LR(:) + st%dustar_dPstar * DerPstardU_LR_Glob(:)

      else

        BB=pstar-(Cs_cell(Con_R)-Cs_cell(Con_Ene)+0.5_dp*(Cs_cell(Con_Qdm)**2)/Cs_cell(Con_Mas))
        DerBBdU_LR(1)=0.5_dp*(Cs_cell(Con_Qdm)**2)/(Cs_cell(Con_Mas)**2)
        DerBBdU_LR(2)=-Cs_cell(Con_Qdm)/Cs_cell(Con_Mas)
        DerBBdU_LR(3)=1.0_dp
        DerBBdU_LR(4)=-1.0_dp

        call derUstar_acoustic(BB, DerBBdU_LR, dc, sgnBC, c_LR, Ro_LR, u_LR, DerUstardU_LR)
        DerhdU_LR(:)=Ustar*DerUstardU_LR(:)

        DerPstardU_LR_Glob(:)=-DerhdU_LR(:)/DerhdPstar

        DerRostardU_LR_Glob(:)=DerPstardU_LR_Glob(:)/(Ctank**2)

        BB=pstar-(Cs_cell(Con_R)-Cs_cell(Con_Ene)+0.5_dp*(Cs_cell(Con_Qdm)**2)/Cs_cell(Con_Mas))
        DerBBdU_LR(1)=DerPstardU_LR_Glob(1)+0.5_dp*(Cs_cell(Con_Qdm)**2)/(Cs_cell(Con_Mas)**2)
        DerBBdU_LR(2)=DerPstardU_LR_Glob(2)-Cs_cell(Con_Qdm)/Cs_cell(Con_Mas)
        DerBBdU_LR(3)=DerPstardU_LR_Glob(3)+1.0_dp
        DerBBdU_LR(4)=DerPstardU_LR_Glob(4)-1.0_dp

        call derUstar_acoustic(BB, DerBBdU_LR, dc, sgnBC, c_LR, Ro_LR, u_LR, DerUstardU_LR_Glob)

      end if

      call jacobian_roT(rostar, Tstar, dedT, dPdT, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
      dRdP=(dRstardRo+dRstardT*dTdRo_p)/(Ctank**2)+dRstardT*dTdp_Ro
      DerRcorrstardU_LR_Glob(:)=dRdP*DerPstardU_LR_Glob(:)

      DerConstardU_LR(:,1)=DerRostardU_LR_Glob(:)
      DerConstardU_LR(:,2)=Rostar*DerUstardU_LR_Glob(:)+DerRostardU_LR_Glob(:)*Ustar
      DerConstardU_LR(:,3)=DerRcorrstardU_LR_Glob(:)+0.5_dp*(DerRostardU_LR_Glob(:)*Ustar*Ustar+&
                         2.0_dp*Rostar*Ustar*DerUstardU_LR_Glob(:))-DerPstardU_LR_Glob(:)
      DerConstardU_LR(:,4)=DerRcorrstardU_LR_Glob(:)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      me%port%derFlx_derCon(:,:)=matmul(DerConstardU_LR(:,:),Jacob_star(:,:))

      me%port%derVit_derCon(:)=DerUstardU_LR_Glob(:)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else

      PrimStar(Pri_ro)=Rostar
      PrimStar(Pri_u)=Ustar
      PrimStar(Pri_p)=pstar
      PrimStar(Pri_e)=estar
      PrimStar(Pri_T)=Tstar
      PrimStar(Pri_c)=Cstar
      PrimStar(Pri_R)=R_Correction_Star
      call jacobian_star_Req_rem(PrimStar,Jacob_star)

      if (use_hugoniot) then
        DerhdPstar=DerFunct
      else
        DerhdPstar=derivative_function_of_pstar_PT(sgnBC,Pstar,Rotank,Ctank,Ptank,p_LR,u_LR,z_LR)
      end if

      HH=Pr_cell(Pri_e)+0.5_dp*Pr_cell(Pri_u)**2+Pr_cell(Pri_p)/Pr_cell(Pri_ro)
      call jacobian_roT(Pr_cell(Pri_ro), Pr_cell(Pri_T), dedT, dPdT)
      dPde=dPdT/dedT
      kk=dPde/Pr_cell(Pri_ro)
      KKK=Pr_cell(Pri_c)**2+kk*(Pr_cell(Pri_u)**2-HH)

      if (use_hugoniot) then
        v_LR       = 1.0_dp / Ro_LR
        DeltaP_hug = Pstar - p_LR
        Deltav_hug = v_LR - st%vstar

        dvstar_dU(1) = st%dvstar_dv_LR*(-v_LR**2) + st%dvstar_dp_LR*KKK
        dvstar_dU(2) = st%dvstar_dp_LR*(-kk*u_LR)
        dvstar_dU(3) = st%dvstar_dp_LR*kk
        dvstar_dU(4) = 0.0_dp

        dPsi_dU(1) = ((-v_LR**2 - dvstar_dU(1))*DeltaP_hug + Deltav_hug*KKK)       / DeltaP_hug**2
        dPsi_dU(2) = (           -dvstar_dU(2) *DeltaP_hug + Deltav_hug*(-kk*u_LR)) / DeltaP_hug**2
        dPsi_dU(3) = (           -dvstar_dU(3) *DeltaP_hug + Deltav_hug*kk)         / DeltaP_hug**2
        dPsi_dU(4) = 0.0_dp

        DerUstardU_LR(1) = -u_LR*v_LR - sgnBC*st%K*KKK     + sgnBC*DeltaP_hug*dPsi_dU(1)/(2.0_dp*st%K)
        DerUstardU_LR(2) =  v_LR      + sgnBC*st%K*kk*u_LR + sgnBC*DeltaP_hug*dPsi_dU(2)/(2.0_dp*st%K)
        DerUstardU_LR(3) =            - sgnBC*st%K*kk       + sgnBC*DeltaP_hug*dPsi_dU(3)/(2.0_dp*st%K)
        DerUstardU_LR(4) = 0.0_dp

        DerhdU_LR(:) = Ustar*DerUstardU_LR(:)
        DerPstardU_LR_Glob(:) = -DerhdU_LR(:)/DerhdPstar
        DerRostardU_LR_Glob(:) = DerPstardU_LR_Glob(:)/(Ctank**2)

        DerUstardU_LR_Glob(:) = DerUstardU_LR(:) + st%dustar_dPstar*DerPstardU_LR_Glob(:)

      else
        BB=Pstar-p_LR
        DerBBdU_LR(1)=-KKK
        DerBBdU_LR(2)=kk*Pr_cell(Pri_u)
        DerBBdU_LR(3)=-kk
        DerBBdU_LR(4)=0.0_dp

        dc(Con_Mas)=dcdro+dcdp*KKK
        dc(Con_Qdm)=-dcdp*kk*Pr_cell(Pri_u)
        dc(Con_Ene)=dcdp*kk
        dc(Con_R)=0.0_dp

        call derUstar_acoustic(BB, DerBBdU_LR, dc, sgnBC, c_LR, Ro_LR, u_LR, DerUstardU_LR)
        DerhdU_LR(:)=Ustar*DerUstardU_LR(:)

        DerPstardU_LR_Glob(:)=-DerhdU_LR(:)/DerhdPstar

        DerRostardU_LR_Glob(:)=DerPstardU_LR_Glob(:)/(Ctank**2)

        BB=Pstar-p_LR
        DerBBdU_LR(1)=DerPstardU_LR_Glob(1)-KKK
        DerBBdU_LR(2)=DerPstardU_LR_Glob(2)+kk*Pr_cell(Pri_u)
        DerBBdU_LR(3)=DerPstardU_LR_Glob(3)-kk
        DerBBdU_LR(4)=0.0_dp

        call derUstar_acoustic(BB, DerBBdU_LR, dc, sgnBC, c_LR, Ro_LR, u_LR, DerUstardU_LR_Glob)
      end if

      call jacobian_roT(rostar, Tstar, dedT, dPdT, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
      dRdP=(dRstardRo+dRstardT*dTdRo_p)/(Ctank**2)+dRstardT*dTdp_Ro
      DerRphysstardU_LR_Glob(:)=dRdP*DerPstardU_LR_Glob(:)

      DerRcorrstardU_LR_Glob(:)=0.0_dp

      DerConstardU_LR(:,1)=DerRostardU_LR_Glob(:)
      DerConstardU_LR(:,2)=Rostar*DerUstardU_LR_Glob(:)+DerRostardU_LR_Glob(:)*Ustar
      DerConstardU_LR(:,3)=DerRphysstardU_LR_Glob(:)+0.5_dp*(DerRostardU_LR_Glob(:)*Ustar*Ustar+&
                         2.0_dp*Rostar*Ustar*DerUstardU_LR_Glob(:))-DerPstardU_LR_Glob(:)
      DerConstardU_LR(:,4)=0.0_dp! DerRcorrstardU_LR_Glob(:)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      me%port%derFlx_derCon(:,:)=matmul(DerConstardU_LR(:,:),Jacob_star(:,:))
      me%port%derFlx_derCon(:,4)=0.0_dp

      me%port%derVit_derCon(:)=DerUstardU_LR_Glob(:)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    endif

  else

    if(R_Correction) then

      if (use_hugoniot) then
        v_LR       = 1.0_dp / Ro_LR
        DeltaP_hug = Pstar - p_LR
        Deltav_hug = v_LR - st%vstar

        dvstar_dU(1) = st%dvstar_dv_LR*(-v_LR**2) + st%dvstar_dp_LR*(-u_LR**2/2.0_dp)
        dvstar_dU(2) = st%dvstar_dp_LR*u_LR
        dvstar_dU(3) = st%dvstar_dp_LR*(-1.0_dp)
        dvstar_dU(4) = st%dvstar_dp_LR*1.0_dp

        dPsi_dU(1) = ((-v_LR**2 - dvstar_dU(1))*DeltaP_hug + Deltav_hug*(-u_LR**2/2.0_dp)) / DeltaP_hug**2
        dPsi_dU(2) = (            -dvstar_dU(2) *DeltaP_hug + Deltav_hug*u_LR)              / DeltaP_hug**2
        dPsi_dU(3) = (            -dvstar_dU(3) *DeltaP_hug + Deltav_hug*(-1.0_dp))         / DeltaP_hug**2
        dPsi_dU(4) = (            -dvstar_dU(4) *DeltaP_hug + Deltav_hug*1.0_dp)            / DeltaP_hug**2

        DerUstardU_LR(1) = -u_LR*v_LR + sgnBC*st%K*(u_LR**2/2.0_dp) + sgnBC*DeltaP_hug*dPsi_dU(1)/(2.0_dp*st%K)
        DerUstardU_LR(2) =  v_LR      - sgnBC*st%K*u_LR              + sgnBC*DeltaP_hug*dPsi_dU(2)/(2.0_dp*st%K)
        DerUstardU_LR(3) =              sgnBC*st%K                    + sgnBC*DeltaP_hug*dPsi_dU(3)/(2.0_dp*st%K)
        DerUstardU_LR(4) =            - sgnBC*st%K                    + sgnBC*DeltaP_hug*dPsi_dU(4)/(2.0_dp*st%K)

        DerRostardU_LR(:)    = -dvstar_dU(:) * Rostar**2
        DerRcorrstardU_LR(:) = (estar - st%ev_v*st%vstar) * DerRostardU_LR(:)

      else

        BB=pstar-(Cs_cell(Con_R)-Cs_cell(Con_Ene)+0.5_dp*(Cs_cell(Con_Qdm)**2)/Cs_cell(Con_Mas))
        DerBBdU_LR(1)=0.5_dp*(Cs_cell(Con_Qdm)**2)/(Cs_cell(Con_Mas)**2)
        DerBBdU_LR(2)=-Cs_cell(Con_Qdm)/Cs_cell(Con_Mas)
        DerBBdU_LR(3)=1.0_dp
        DerBBdU_LR(4)=-1.0_dp

        call derUstar_acoustic(BB, DerBBdU_LR, dc, sgnBC, c_LR, Ro_LR, u_LR, DerUstardU_LR)

        DerRostardU_LR(1)=(DerBBdU_LR(1)*(c_LR**2)-2.0_dp*BB*c_LR*dc(1))/(c_LR**4)+1.0_dp
        DerRostardU_LR(2:4)=(DerBBdU_LR(2:4)*(c_LR**2)-2.0_dp*BB*c_LR*dc(2:4))/(c_LR**4)

        DerRcorrstardU_LR(1)=-Cs_cell(Con_R)*Rostar/(Cs_cell(Con_Mas)**2)+Cs_cell(Con_R)*DerRostardU_LR(1)/Cs_cell(Con_Mas)+&
                           2.0_dp*c_LR*dc(1)*(Rostar-Cs_cell(Con_Mas))+(c_LR**2)*(DerRostardU_LR(1)-1.0_dp)
        DerRcorrstardU_LR(2:3)=Cs_cell(Con_R)*DerRostardU_LR(2:3)/Cs_cell(Con_Mas)+2.0_dp*c_LR*dc(2:3)*(Rostar-Cs_cell(Con_Mas))+&
                            (c_LR**2)*DerRostardU_LR(2:3)
        DerRcorrstardU_LR(4)=Rostar/Cs_cell(Con_Mas)+Cs_cell(Con_R)*DerRostardU_LR(4)/Cs_cell(Con_Mas)+2.0_dp*c_LR*dc(4)*&
                           (Rostar-Cs_cell(Con_Mas))+(c_LR**2)*DerRostardU_LR(4)

      end if

      DerConstardU_LR(:,1)=DerRostardU_LR(:)
      DerConstardU_LR(:,2)=Rostar*DerUstardU_LR(:)+DerRostardU_LR(:)*Ustar
      DerConstardU_LR(:,3)=DerRcorrstardU_LR(:)+0.5_dp*(DerRostardU_LR(:)*Ustar*Ustar+2.0_dp*Rostar*Ustar*DerUstardU_LR(:))
      DerConstardU_LR(:,4)=DerRcorrstardU_LR(:)
    
      me%port%derFlx_derCon(:,1)=DerConstardU_LR(:,2)
      me%port%derFlx_derCon(:,2)=DerConstardU_LR(:,2)*Ustar+Rostar*Ustar*DerUstardU_LR(:)
      me%port%derFlx_derCon(:,3)=DerConstardU_LR(:,3)*Ustar+(Rostar*(estar+0.5_dp*Ustar*Ustar)+pstar)*DerUstardU_LR(:)
      me%port%derFlx_derCon(:,4)=DerRcorrstardU_LR(:)*Ustar+R_Correction_Star*DerUstardU_LR(:)

      me%port%derVit_derCon(:)=DerUstardU_LR(:)

    else

      HH=Pr_cell(Pri_e)+0.5_dp*Pr_cell(Pri_u)**2+Pr_cell(Pri_p)/Pr_cell(Pri_ro)
      call jacobian_roT(Pr_cell(Pri_ro), Pr_cell(Pri_T), dedT, dPdT)
      dPde=dPdT/dedT
      kk=dPde/Pr_cell(Pri_ro)
      KKK=Pr_cell(Pri_c)**2+kk*(Pr_cell(Pri_u)**2-HH)

      if (use_hugoniot) then

        v_LR       = 1.0_dp / Ro_LR
        DeltaP_hug = Pstar - p_LR
        Deltav_hug = v_LR - st%vstar

        dvstar_dU(1) = st%dvstar_dv_LR*(-v_LR**2) + st%dvstar_dp_LR*KKK
        dvstar_dU(2) = st%dvstar_dp_LR*(-kk*u_LR)
        dvstar_dU(3) = st%dvstar_dp_LR*kk
        dvstar_dU(4) = 0.0_dp

        dPsi_dU(1) = ((-v_LR**2 - dvstar_dU(1))*DeltaP_hug + Deltav_hug*KKK)       / DeltaP_hug**2
        dPsi_dU(2) = (           -dvstar_dU(2) *DeltaP_hug + Deltav_hug*(-kk*u_LR)) / DeltaP_hug**2
        dPsi_dU(3) = (           -dvstar_dU(3) *DeltaP_hug + Deltav_hug*kk)         / DeltaP_hug**2
        dPsi_dU(4) = 0.0_dp

        DerUstardU_LR(1) = -u_LR*v_LR - sgnBC*st%K*KKK     + sgnBC*DeltaP_hug*dPsi_dU(1)/(2.0_dp*st%K)
        DerUstardU_LR(2) =  v_LR      + sgnBC*st%K*kk*u_LR + sgnBC*DeltaP_hug*dPsi_dU(2)/(2.0_dp*st%K)
        DerUstardU_LR(3) =            - sgnBC*st%K*kk       + sgnBC*DeltaP_hug*dPsi_dU(3)/(2.0_dp*st%K)
        DerUstardU_LR(4) = 0.0_dp

        DerRostardU_LR(:) = -dvstar_dU(:) * Rostar**2

      else

        BB=Pstar-p_LR
        DerBBdU_LR(1)=-KKK
        DerBBdU_LR(2)=kk*Pr_cell(Pri_u)
        DerBBdU_LR(3)=-kk
        DerBBdU_LR(4)=0.0_dp

        dc(Con_Mas)=dcdro+dcdp*KKK
        dc(Con_Qdm)=-dcdp*kk*Pr_cell(Pri_u)
        dc(Con_Ene)=dcdp*kk
        dc(Con_R)=0.0_dp

        call derUstar_acoustic(BB, DerBBdU_LR, dc, sgnBC, c_LR, Ro_LR, u_LR, DerUstardU_LR)

        DerRostardU_LR(1)=(DerBBdU_LR(1)*(c_LR**2)-2.0_dp*BB*c_LR*dc(1))/(c_LR**4)+1.0_dp
        DerRostardU_LR(2)=(DerBBdU_LR(2)*(c_LR**2)-2.0_dp*BB*c_LR*dc(2))/(c_LR**4)
        DerRostardU_LR(3)=(DerBBdU_LR(3)*(c_LR**2)-2.0_dp*BB*c_LR*dc(3))/(c_LR**4)
        DerRostardU_LR(4)=0.0_dp

      end if

      call jacobian_roT(rostar, Tstar, dedT, dPdT, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
      dR_dRo=dRstardRo+dRstardT*dTdRo_p
      DerRphystardU_LR(:)=dR_dRo*DerRostardU_LR(:)

      DerRcorrstardU_LR(:)=0.0_dp
      DerRcorrstardU_LR(4)=1.0_dp

      DerConstardU_LR(:,1)=DerRostardU_LR(:)
      DerConstardU_LR(:,2)=Rostar*DerUstardU_LR(:)+DerRostardU_LR(:)*Ustar
      DerConstardU_LR(:,3)=DerRphystardU_LR(:)+0.5_dp*(DerRostardU_LR(:)*Ustar*Ustar+2.0_dp*Rostar*Ustar*DerUstardU_LR(:))
      DerConstardU_LR(:,4)=DerRcorrstardU_LR(:)

      me%port%derFlx_derCon(:,1)=DerConstardU_LR(:,2)
      me%port%derFlx_derCon(:,2)=DerConstardU_LR(:,2)*Ustar+Rostar*Ustar*DerUstardU_LR(:)
      me%port%derFlx_derCon(:,3)=DerConstardU_LR(:,3)*Ustar+(Rostar*(estar+0.5_dp*Ustar*Ustar)+pstar)*DerUstardU_LR(:)
      me%port%derFlx_derCon(:,4)=0.0_dp

      me%port%derVit_derCon(:)=DerUstardU_LR(:)

    endif

  endif

  if (use_hugoniot) then
    ss = abs(u_LR) + max(c_LR, st%W / Ro_LR)
  else
    ss = abs(u_LR) + c_LR
  end if
end subroutine solve_boundary_PT


pure subroutine derUstar_acoustic(BB, DerBBdU, dc, sgnBC, c_LR, ro_LR, u_LR, DerUstar)
  real(dp), intent(in)  :: BB, DerBBdU(Nb_VarC), dc(Nb_VarC), sgnBC, c_LR, ro_LR, u_LR
  real(dp), intent(out) :: DerUstar(Nb_VarC)
  real(dp) :: inv_roc2
  inv_roc2 = sgnBC / (c_LR**2 * ro_LR)
  DerUstar(:) = inv_roc2 * (DerBBdU(:)*c_LR - BB*dc(:))
  DerUstar(1) = DerUstar(1) - u_LR/ro_LR - sgnBC*BB/(ro_LR**2 * c_LR)
  DerUstar(2) = DerUstar(2) + 1.0_dp/ro_LR
end subroutine derUstar_acoustic


function hug_brent_f(me, x) result(H)
  class(hug_brent_t), intent(in) :: me
  real(dp),           intent(in) :: x
  real(dp) :: H, e_v, ev_v, ep_v, droeint_dro
  droeint_dro = droeint_droP(1.0_dp/x, me%Pstar, e=e_v, dedv=ev_v, dedp=ep_v)
  H = (me%Pstar + me%p_LR)*(x - me%v_LR) + 2.0_dp*(e_v - me%e_LR)
end function hug_brent_f


subroutine hugoniot_star_state(sgnBC, Pstar, Ro_LR, p_LR, u_LR, c_LR, st)
  real(dp),         intent(in)  :: sgnBC, Pstar, Ro_LR, p_LR, u_LR, c_LR
  type(hug_star_t), intent(out) :: st

  integer,  parameter :: NbreMaxIte = 200
  real(dp), parameter :: MachError = 1.0e-15_dp, Tol_brent = 1.0e-10_dp

  real(dp) :: DeltaP, v_LR
  real(dp) :: e_LR, ev_LR, ep_LR
  real(dp) :: v, v_pred, v_lo, v_hi, Psi
  real(dp) :: e_v, ev_v, ep_v
  real(dp) :: Funct, DerFunct, D, tol_H, droeint_dro
  integer  :: ite
  logical  :: have_lo, have_hi
  type(hug_brent_t) :: hug_obj

  st%ok     = .false.
  st%n_iter = 0

  DeltaP = Pstar - p_LR
  v_LR   = 1.0_dp / Ro_LR

  droeint_dro = droeint_droP(Ro_LR, p_LR, e=e_LR, dedv=ev_LR, dedp=ep_LR)

  v_pred = v_LR - DeltaP * v_LR**2 / c_LR**2
  if (v_pred <= 0.0_dp) v_pred = v_LR * 0.5_dp
  v = v_pred

  hug_obj%Pstar = Pstar
  hug_obj%p_LR  = p_LR
  hug_obj%v_LR  = v_LR
  hug_obj%e_LR  = e_LR

  v_lo = 0.0_dp; v_hi = 0.0_dp
  have_lo = .false.; have_hi = .false.

  tol_H = 1.0e-8_dp * (Pstar + p_LR) * v_LR

  do ite = 1, NbreMaxIte
    droeint_dro = droeint_droP(1.0_dp/v, Pstar, e=e_v, dedv=ev_v, dedp=ep_v)
    Funct   = (Pstar + p_LR)*(v - v_LR) + 2.0_dp*(e_v - e_LR)
    DerFunct = (Pstar + p_LR) + 2.0_dp*ev_v
    if (Funct < 0.0_dp) then
      if (.not. have_lo .or. v > v_lo) then; v_lo = v; have_lo = .true.; endif
    else
      if (.not. have_hi .or. v < v_hi) then; v_hi = v; have_hi = .true.; endif
    endif
    if (abs(Funct) > tol_H) then
      v = v - Funct/DerFunct
      if (v <= 0.0_dp) exit
    else
      st%n_iter = ite; exit
    endif
    if (ite == NbreMaxIte) then
      if (have_lo .and. have_hi) then
        v = zero(hug_obj, v_lo, v_hi, MachError, Tol_brent)
        droeint_dro = droeint_droP(1.0_dp/v, Pstar, e=e_v, dedv=ev_v, dedp=ep_v)
      else
        call set_error('hugoniot_star_state: Newton did not converge, no bracket')
        return
      endif
      exit
    endif
  enddo

  if (v <= 0.0_dp) then
    if (.not. (have_lo .and. have_hi)) then
      call set_error('hugoniot_star_state: cannot bracket root')
      return
    endif
    v = zero(hug_obj, v_lo, v_hi, MachError, Tol_brent)
    droeint_dro = droeint_droP(1.0_dp/v, Pstar, e=e_v, dedv=ev_v, dedp=ep_v)
  endif

  Psi = (v_LR - v) / DeltaP
  if (Psi <= 0.0_dp) then
    call set_error('hugoniot_star_state: Psi = (v_LR - v*)/DeltaP <= 0')
    return
  endif

  st%vstar = v
  st%K     = sqrt(Psi)
  st%W     = 1.0_dp / st%K
  st%ustar = u_LR + sgnBC * DeltaP * st%K
  st%ok    = .true.

  D = (Pstar + p_LR) + 2.0_dp * ev_v
  st%dvstar_dPstar = -((v - v_LR) + 2.0_dp * ep_v) / D
  st%dustar_dPstar = sgnBC * (st%K + (-st%dvstar_dPstar - Psi) / (2.0_dp * st%K))
  st%dvstar_dv_LR  = ((Pstar + p_LR) + 2.0_dp * ev_LR) / D
  st%dvstar_dp_LR  = -((v - v_LR) - 2.0_dp * ep_LR) / D
  st%ev_v          = ev_v

end subroutine hugoniot_star_state


function function_of_pstar_MT(sgnBC,Pstar,Ttank,mdot,AreaPP,p_LR,u_LR,z_LR)

  real(dp), intent(in) :: sgnBC,Pstar,Ttank,mdot,AreaPP,p_LR,u_LR,z_LR
  real(dp) :: function_of_pstar_MT, Ustar, rostar_LR, T_temp

  Ustar=u_LR+sgnBC*(Pstar-p_LR)/z_LR
  rostar_LR=mdot/(Ustar*AreaPP)
  T_temp=T_roP(rostar_LR,Pstar)
  function_of_pstar_MT=T_temp-Ttank

end function function_of_pstar_MT



subroutine solve_boundary_MT(me,SS)     
  type(boundary_t), intent(inout) :: me
  real(dp), intent(out) :: ss
 
  real(dp) :: Ttank, ro_LR, p_LR, u_LR, c_LR, z_LR, R_Correction_Star,R_Phys_Star
  real(dp) :: Pstar, Tstar, Ustar, Rostar, estar, Cstar
  real(dp) :: DerhdPstar,dcdro,dcdp,BB,dTdp_Ro,dTdRo_p,dRstardRo,dRstardT,Rstar,T_temp,T_temp_plus
  real(dp), dimension(Nb_VarC) :: DerBBdU_LR,DerUstardU_LR_Glob,Cs_cell,Cons_star
  real(dp), dimension(Nb_VarC) :: DerPstardU_LR_Glob,DerRostardU_LR_Glob,DerRcorrstardU_LR_Glob,DerRphysstardU_LR_Glob
  real(dp), dimension(Nb_VarC) :: DerUstardU_LR,DerhdU_LR,dc,DerRostardU_LR
  real(dp), dimension(Nb_VarC,Nb_VarC) :: Jacob_star,DerConstardU_LR
  real(dp) :: kk, KKK, dedT, dPdT, HH, dPde
  real(dp), dimension(Nb_VarP) :: PrimStar
  real(dp), dimension(Nb_VarP) :: Pr_cell
  integer :: NbreMaxIte, ite, i
  real(dp) :: pmin,pmid,pmax,f_pmin,f_pmid,f_pmax,AreaPP,mdot,sgnBC,Nbcel,dpp,rostarLoc,FctLoc

  Pr_cell(:)=me%port%Prim(:)
  Cs_cell(:)=me%port%Cons(:)
  AreaPP=me%port%Area
  sgnBC=-me%port%Sgn4j  ! sgnBC=1 if 1-port, -1 if 2-port

  Ttank=me%imposed_T%v0d()
  mdot=me%imposed_PorM%v0d()

  ro_LR=Pr_cell(Pri_ro)
  p_LR=Pr_cell(Pri_P)
  u_LR=Pr_cell(Pri_u)
  c_LR=Pr_cell(Pri_c)
  z_LR=ro_LR*c_LR

  ! Analysis of increasing function
  pmin=p_LR/2.0_dp
  pmax=p_LR*2.0_dp
  Nbcel=100000.0_dp
  dpp=(pmax-pmin)/Nbcel
  do i=1,Nbcel
    rostarLoc=mdot/((u_LR+sgnBC*(pmin+dpp*(i-1)-p_LR)/z_LR)*AreaPP)
    if(rostarLoc<=0.0_dp .or. rostarLoc>200.0_dp) cycle
    FctLoc=function_of_pstar_MT(sgnBC,pmin+dpp*(i-1),Ttank,mdot,AreaPP,p_LR,u_LR,z_LR)/1.0e3_dp
    pmin=pmin+dpp*(i-1)
    exit
  enddo

  ! Dichotomy
  f_pmin=function_of_pstar_MT(sgnBC,pmin,Ttank,mdot,AreaPP,p_LR,u_LR,z_LR)/1.0e3_dp

  f_pmax=function_of_pstar_MT(sgnBC,pmax,Ttank,mdot,AreaPP,p_LR,u_LR,z_LR)/1.0e3_dp  

  if(f_pmin*f_pmax>0.0_dp) then
    call set_error('initial range does not contain solution in solve_boundary_MT')
    return
  else
    NbreMaxIte=200
    do ite=1,NbreMaxIte
      pmid=0.5_dp*(pmin+pmax)
      f_pmid=function_of_pstar_MT(sgnBC,pmid,Ttank,mdot,AreaPP,p_LR,u_LR,z_LR)/1.0e3_dp
      if(f_pmin*f_pmid<0.0_dp) then
        pmax=pmid
      else
        pmin=pmid
        f_pmin=f_pmid
      endif
      if(abs(pmax-pmin)<1.0e-2_dp) exit
      if(ite==NbreMaxIte) then
        call set_error('convergence failure in solve_boundary_MT')
        return
      endif
    enddo  
  endif
  Pstar=0.5_dp*(pmin+pmax)

  Ustar=u_LR+sgnBC*(Pstar-p_LR)/z_LR
  Rostar=mdot/(Ustar*AreaPP)    ! Ustar and mdot always share the same sign
  call state_roP(Rostar, Pstar, R_Phys_Star, estar, Tstar, Cstar)
  R_Correction_Star=1.0_dp
  if(R_Correction) R_Correction_Star=R_Phys_Star

  me%port%flx(Con_Mas)=Rostar*Ustar
  me%port%flx(Con_Qdm)=Rostar*Ustar*Ustar+Pstar
  me%port%flx(Con_Ene)=(Rostar*(estar+0.5_dp*Ustar*Ustar)+Pstar)*Ustar
  me%port%flx(Con_R)=0.0_dp
  if(R_Correction) me%port%flx(Con_R)=R_Correction_Star*Ustar

  me%port%vit=Ustar

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                            Derivatives                           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call dc2_roT(ro_LR, Pr_cell(Pri_T), dcdro, dcdp)
  dcdro = dcdro/(2.0_dp*c_LR)
  dcdp  = dcdp /(2.0_dp*c_LR)

  dc(Con_Mas)=dcdro-0.5_dp*dcdp*(Cs_cell(Con_Qdm)**2)/(Cs_cell(Con_Mas)**2)
  dc(Con_Qdm)=dcdp*Cs_cell(Con_Qdm)/Cs_cell(Con_Mas)
  dc(Con_Ene)=-dcdp
  dc(Con_R)=dcdp   

  if(R_Correction) then

    Cons_star(Con_Mas)=Rostar
    Cons_star(Con_Qdm)=Rostar*Ustar
    Cons_star(Con_Ene)=Rostar*(estar+0.5_dp*Ustar*Ustar)
    Cons_star(Con_R)=R_Correction_Star
    call jacobian_star(Cons_star,Jacob_star)

    DerhdPstar=(function_of_pstar_MT(sgnBC,Pstar+1.0e-3_dp,Ttank,mdot,AreaPP,p_LR,u_LR,z_LR)-&
                function_of_pstar_MT(sgnBC,Pstar,Ttank,mdot,AreaPP,p_LR,u_LR,z_LR))    

    BB=pstar-(Cs_cell(Con_R)-Cs_cell(Con_Ene)+0.5_dp*(Cs_cell(Con_Qdm)**2)/Cs_cell(Con_Mas))
    DerBBdU_LR(1)=0.5_dp*(Cs_cell(Con_Qdm)**2)/(Cs_cell(Con_Mas)**2)
    DerBBdU_LR(2)=-Cs_cell(Con_Qdm)/Cs_cell(Con_Mas)
    DerBBdU_LR(3)=1.0_dp
    DerBBdU_LR(4)=-1.0_dp            
    
    DerUstardU_LR(1)=-1.0_dp/(Cs_cell(Con_Mas)**2)*(Cs_cell(Con_Qdm)+sgnBC*BB/c_LR)+&
                   +(sgnBC*(DerBBdU_LR(1)*c_LR-BB*dc(1))/(c_LR**2))/Cs_cell(Con_Mas)
                   DerUstardU_LR(2)=(1.0_dp+sgnBC*(DerBBdU_LR(2)*c_LR-BB*dc(2))/(c_LR**2))/Cs_cell(Con_Mas)
    DerUstardU_LR(3)=(sgnBC*(DerBBdU_LR(3)*c_LR-BB*dc(3))/(c_LR**2))/Cs_cell(Con_Mas)
    DerUstardU_LR(4)=(sgnBC*(DerBBdU_LR(4)*c_LR-BB*dc(4))/(c_LR**2))/Cs_cell(Con_Mas)

    DerRostardU_LR(:)=-(mdot/AreaPP)*DerUstardU_LR(:)/(Ustar**2)

    T_temp=T_roP(Rostar,Pstar)
    T_temp_plus=T_roP(Rostar+1.0e-3_dp,Pstar)
    DerhdU_LR(:)=DerRostardU_LR(:)*(T_temp_plus-T_temp)
    DerPstardU_LR_Glob(:)=-DerhdU_LR(:)/DerhdPstar    

    BB=pstar-(Cs_cell(Con_R)-Cs_cell(Con_Ene)+0.5_dp*(Cs_cell(Con_Qdm)**2)/Cs_cell(Con_Mas))
    DerBBdU_LR(1)=DerPstardU_LR_Glob(1)+0.5_dp*(Cs_cell(Con_Qdm)**2)/(Cs_cell(Con_Mas)**2)
    DerBBdU_LR(2)=DerPstardU_LR_Glob(2)-Cs_cell(Con_Qdm)/Cs_cell(Con_Mas)
    DerBBdU_LR(3)=DerPstardU_LR_Glob(3)+1.0_dp
    DerBBdU_LR(4)=DerPstardU_LR_Glob(4)-1.0_dp    

    DerUstardU_LR_Glob(1)=-1.0_dp/(Cs_cell(Con_Mas)**2)*(Cs_cell(Con_Qdm)+sgnBC*BB/c_LR)+&
                        +(sgnBC*(DerBBdU_LR(1)*c_LR-BB*dc(1))/(c_LR**2))/Cs_cell(Con_Mas)
    DerUstardU_LR_Glob(2)=(1.0_dp+sgnBC*(DerBBdU_LR(2)*c_LR-BB*dc(2))/(c_LR**2))/Cs_cell(Con_Mas)
    DerUstardU_LR_Glob(3)=(sgnBC*(DerBBdU_LR(3)*c_LR-BB*dc(3))/(c_LR**2))/Cs_cell(Con_Mas)
    DerUstardU_LR_Glob(4)=(sgnBC*(DerBBdU_LR(4)*c_LR-BB*dc(4))/(c_LR**2))/Cs_cell(Con_Mas)    

    DerRostardU_LR_Glob(:)=-(mdot/AreaPP)*DerUstardU_LR_Glob(:)/(Ustar**2)

    call jacobian_roT(rostar, Tstar, dedT, dPdT, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
    DerRcorrstardU_LR_Glob(:)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardU_LR_Glob(:)+dRstardT*dTdp_Ro*DerPstardU_LR_Glob(:) 

    DerConstardU_LR(:,1)=DerRostardU_LR_Glob(:)
    DerConstardU_LR(:,2)=Rostar*DerUstardU_LR_Glob(:)+DerRostardU_LR_Glob(:)*Ustar
    DerConstardU_LR(:,3)=DerRcorrstardU_LR_Glob(:)+0.5_dp*(DerRostardU_LR_Glob(:)*Ustar*Ustar+&
                         2.0_dp*Rostar*Ustar*DerUstardU_LR_Glob(:))-DerPstardU_LR_Glob(:)
    DerConstardU_LR(:,4)=DerRcorrstardU_LR_Glob(:)

    me%port%derFlx_derCon(:,:)=matmul(DerConstardU_LR(:,:),Jacob_star(:,:))
    
    me%port%derVit_derCon(:)=DerUstardU_LR_Glob(:)    

  else

    PrimStar(Pri_ro)=Rostar
    PrimStar(Pri_u)=Ustar
    PrimStar(Pri_p)=pstar
    PrimStar(Pri_e)=estar
    PrimStar(Pri_T)=Tstar
    PrimStar(Pri_c)=Cstar
    PrimStar(Pri_R)=R_Correction_Star
    call jacobian_star_Req_rem(PrimStar,Jacob_star)

    DerhdPstar=(function_of_pstar_MT(sgnBC,Pstar+1.0e-3_dp,Ttank,mdot,AreaPP,p_LR,u_LR,z_LR)-&
                function_of_pstar_MT(sgnBC,Pstar,Ttank,mdot,AreaPP,p_LR,u_LR,z_LR))

    HH=Pr_cell(Pri_e)+0.5_dp*Pr_cell(Pri_u)**2+Pr_cell(Pri_p)/Pr_cell(Pri_ro)
    call jacobian_roT(Pr_cell(Pri_ro), Pr_cell(Pri_T), dedT, dPdT)
    dPde=dPdT/dedT
    kk=dPde/Pr_cell(Pri_ro)
    KKK=Pr_cell(Pri_c)**2+kk*(Pr_cell(Pri_u)**2-HH)

    BB=Pstar-p_LR
    DerBBdU_LR(1)=-KKK
    DerBBdU_LR(2)=kk*Pr_cell(Pri_u)
    DerBBdU_LR(3)=-kk
    DerBBdU_LR(4)=0.0_dp

    dc(Con_Mas)=dcdro+dcdp*KKK
    dc(Con_Qdm)=-dcdp*kk*Pr_cell(Pri_u)
    dc(Con_Ene)=dcdp*kk
    dc(Con_R)=0.0_dp

    DerUstardU_LR(1)=-Pr_cell(Pri_u)/Pr_cell(Pri_ro)+sgnBC*(DerBBdU_LR(1)*Pr_cell(Pri_ro)*c_LR-BB*(c_LR+Pr_cell(Pri_ro)*dc(1)))/&
                    ((Pr_cell(Pri_ro)*c_LR)**2)
    DerUstardU_LR(2)=(1.0_dp/Pr_cell(Pri_ro))*(1.0_dp+sgnBC*(DerBBdU_LR(2)*c_LR-BB*dc(2))/(c_LR**2))
    DerUstardU_LR(3)=sgnBC*(1.0_dp/Pr_cell(Pri_ro))*(DerBBdU_LR(3)*c_LR-BB*dc(3))/(c_LR**2)
    DerUstardU_LR(4)=0.0_dp

    DerRostardU_LR(:)=-(mdot/AreaPP)*DerUstardU_LR(:)/(Ustar**2)
     
    T_temp=T_roP(Rostar,Pstar)
    T_temp_plus=T_roP(Rostar+1.0e-3_dp,Pstar)
    DerhdU_LR(:)=DerRostardU_LR(:)*(T_temp_plus-T_temp)
    DerPstardU_LR_Glob(:)=-DerhdU_LR(:)/DerhdPstar

    BB=Pstar-p_LR
    DerBBdU_LR(1)=DerPstardU_LR_Glob(1)-KKK
    DerBBdU_LR(2)=DerPstardU_LR_Glob(2)+kk*Pr_cell(Pri_u)
    DerBBdU_LR(3)=DerPstardU_LR_Glob(3)-kk
    DerBBdU_LR(4)=0.0_dp

    DerUstardU_LR_Glob(1)=-Pr_cell(Pri_u)/Pr_cell(Pri_ro)+sgnBC*(DerBBdU_LR(1)*Pr_cell(Pri_ro)*c_LR-BB*(c_LR+Pr_cell(Pri_ro)*&
                          dc(1)))/((Pr_cell(Pri_ro)*c_LR)**2)
    DerUstardU_LR_Glob(2)=(1.0_dp/Pr_cell(Pri_ro))*(1.0_dp+sgnBC*(DerBBdU_LR(2)*c_LR-BB*dc(2))/(c_LR**2))
    DerUstardU_LR_Glob(3)=sgnBC*(1.0_dp/Pr_cell(Pri_ro))*(DerBBdU_LR(3)*c_LR-BB*dc(3))/(c_LR**2)
    DerUstardU_LR_Glob(4)=0.0_dp

    DerRostardU_LR_Glob(:)=-(mdot/AreaPP)*DerUstardU_LR_Glob(:)/(Ustar**2)

    call jacobian_roT(rostar, Tstar, dedT, dPdT, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
    DerRphysstardU_LR_Glob(:)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardU_LR_Glob(:)+dRstardT*dTdp_Ro*DerPstardU_LR_Glob(:)

    DerConstardU_LR(:,1)=DerRostardU_LR_Glob(:)
    DerConstardU_LR(:,2)=Rostar*DerUstardU_LR_Glob(:)+DerRostardU_LR_Glob(:)*Ustar
    DerConstardU_LR(:,3)=DerRphysstardU_LR_Glob(:)+0.5_dp*(DerRostardU_LR_Glob(:)*Ustar*Ustar+&
                       2.0_dp*Rostar*Ustar*DerUstardU_LR_Glob(:))-DerPstardU_LR_Glob(:)
    DerConstardU_LR(:,4)=0.0_dp     

    me%port%derFlx_derCon(:,:)=matmul(DerConstardU_LR(:,:),Jacob_star(:,:))
    me%port%derFlx_derCon(:,4)=0.0_dp

    me%port%derVit_derCon(:)=DerUstardU_LR_Glob(:)

  endif

  ss=abs(u_LR)+c_LR
end subroutine solve_boundary_MT



subroutine boundary_resolution_from_and_to_ports(me)
    type(boundary_t), intent(inout) :: me
    real(dp) :: wave_time

    if(me%is_PT) then
        call solve_boundary_PT(me,wave_time)
    else
        call solve_boundary_MT(me,wave_time)
    endif
    if (sim_error > 0) return
    me%wave_time = me%port%dxLoc/wave_time
end subroutine boundary_resolution_from_and_to_ports

end module cmp_boundary_calc_m
