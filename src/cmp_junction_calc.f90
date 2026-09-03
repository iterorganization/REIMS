! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_junction_calc_m
    use cmp_junction_init_m
    use lib_He_thermo_m
    use lib_ext_math_m
    use ieee_arithmetic, only: ieee_is_nan
    implicit none

    type der_cons_vars_t
        real(dp) :: dU1 = 0, dU2 = 0, dU3 = 0, dU4 = 0
    end type der_cons_vars_t
    
    type conservative_variables_t
        real(dp) :: VarC(Nb_VarC)
    end type conservative_variables_t

    type conservative_variable_derivatives_t
        real(dp) :: DerVarC(Nb_VarC,Nb_VarC)
    end type conservative_variable_derivatives_t    
  
    type primitive_variables_t
        real(dp) :: VarP(Nb_VarP)
    end type primitive_variables_t

contains


subroutine solve_junction(dyn,FlJ)
    type(junction_dynamic_parameters_t), intent(inout) :: dyn
    type(flux_and_derivatives_t), intent(inout) :: FlJ(:)

    type(conservative_variable_derivatives_t),dimension(dyn%NbTotBr,dyn%NbTotBr) :: DerConStar
    type(primitive_variables_t),              dimension(dyn%NbTotBr) :: PrimStar
    type(conservative_variables_t),           dimension(dyn%NbTotBr) :: ConStar
    type(der_cons_vars_t), dimension(dyn%NbTotBr,dyn%NbTotBr+dyn%NbOut) :: DerFct,DerXStar
    type(der_cons_vars_t), dimension(dyn%NbTotBr,dyn%NbOut) :: DerCoef,DerAlphaGlob
    type(der_cons_vars_t), dimension(dyn%NbTotBr,dyn%NbTotBr) :: DerPStarGlob,DerVelStarGlob
    type(der_cons_vars_t), dimension(dyn%NbTotBr,dyn%NbTotBr) :: DerRPhysStarGlob,DerRCorStarGlob
    type(der_cons_vars_t), dimension(dyn%NbTotBr,dyn%NbTotBr) :: DerRhoStarGlob
    real(dp), dimension(dyn%NbTotBr) :: p_star,rhoStar,vitStar,ETotStar,rCorrStar,EnthalpyStar,eStar
    real(dp), dimension(dyn%NbTotBr) :: rhoStar_Int,QDMStar_Int,EnergyStar_Int,QDMStar,EnergyStar
    real(dp), dimension(dyn%NbTotBr) :: DerVitStar_dU1,DerVitStar_dU2,DerVitStar_dU3,DerVitStar_dU4
    real(dp), dimension(dyn%NbTotBr) :: DerRhoStar_dU1,DerRhoStar_dU2,DerRhoStar_dU3,DerRhoStar_dU4
    real(dp), dimension(dyn%NbTotBr) :: DerETotStar_Int_dU1,DerETotStar_Int_dU2,DerETotStar_Int_dU3
    real(dp), dimension(dyn%NbTotBr) :: DerHStar_dU1,DerHStar_dU2,DerHStar_dU3,DerHStar_dU4
    real(dp), dimension(dyn%NbTotBr) :: DerETotStar_dU1,DerETotStar_dU2,DerETotStar_dU3
    real(dp), dimension(dyn%NbTotBr) :: DerRCorrStar_dU1,DerRCorrStar_dU2,DerRCorrStar_dU3
    real(dp), dimension(dyn%NbTotBr) :: DerEntMixStar_dU1,DerEntMixStar_dU2,DerEntMixStar_dU3
    real(dp), dimension(dyn%NbTotBr) :: DerVitStar_dX,DerRhoStar_dX,DerETotStar_Int_dX
    real(dp), dimension(dyn%NbTotBr) :: DerRCorrStar_dX,DerEntMixStar_dX,DerHStar_dAlpha
    real(dp), dimension(dyn%NbTotBr) :: DerRhoETotStar_dAlpha,DerRhoETotStar_dX,DerHStar_dX
    real(dp), dimension(dyn%NbTotBr) :: DerETotStar_dX,DerETotStar_Int_dU4,DerEntMixStar_dU4
    real(dp), dimension(dyn%NbTotBr) :: DerRCorrStar_dU4,DerETotStar_dU4
    real(dp), dimension(Nb_VarC) :: DerBBdUk,DerDDdUk,DerRhoETotStar_dUk,Der_qj_dUk,Der_qRef_dUk
    real(dp), dimension(Nb_VarC,Nb_VarC) :: dFdFJ,dUdUjj
    real(dp), dimension(Nb_VarC,Nb_VarC,dyn%NbTotBr) :: Jacob_star
    real(dp), dimension(dyn%NbTotBr+dyn%NbOut) :: x, fVec
    real(dp), dimension(dyn%NbTotBr+dyn%NbOut,dyn%NbTotBr+dyn%NbOut) :: fJac, invJac
    real(dp), dimension(dyn%NbOut+1) :: DerCoef_dX
    real(dp), dimension(dyn%NbOut) :: Der_dFdRo_dP,ETotStar_Int
    real(dp), dimension(dyn%NbOut) :: alpha,dFdRo,CoefPLoss,Der_dFdRo_dRhoInt
    real(dp), dimension(dyn%NbTotBr) :: T_temp
    real(dp) :: BB,DD,dTdp_Ro,dTdRo_p,dRStar_dRo,dRStar_dT,RStar
    real(dp) :: denominator,EnthalpyMixStar,EnthalpyMixStarNum,qj,ps_ij
    real(dp) :: Der_qjdX,Der_qRef_dX,Der_qjdAlpha,DerCoef_dAlpha
    real(dp) :: kk, KKK, dedT, dPdT, HH, dPde
    integer :: info, NbOut, j, jj, NbBr, NbIn, i(dyn%NbIn), o(dyn%NbOut)
    logical, dimension(dyn%NbTotBr) :: om, im  
    real(dp), parameter :: gg = 9.81_dp

    ! Preparing short name for indexing, masking and sizes
    NbBr  = dyn%NbTotBr
    NbIn  = dyn%NbIn
    NbOut = dyn%NbOut
    om = .false.; im(1:NbOut)      = .true.
    im = .false.; im(NbOut+1:NbBr) = .true.
    i = [(j, j = NbOut+1, NbBr)]
    o = [(j, j = 1, dyn%NbOut)]

    ! Calling fSolve to provide solution of `p_star` and `alpha`
    x(1:NbBr) = dyn%br%p0
    x(NbBr+1:NbBr+NbOut) = 1.0e-1_dp
    call fSolve(dyn,NbBr+NbOut,x,fVec,1.0e-8_dp,info)
    if (info /= 1) then
        x(1:NbBr) = dyn%p_star_prev
        x(NbBr+1:NbBr+NbOut) = dyn%alpha_prev(1:NbOut)
        call fSolve(dyn,NbBr+NbOut,x,fVec,1.0e-8_dp,info)
        if (info /= 1) then
          call set_error('fSolve did not converge in solve_junction')
          return
        endif          
    end if
    p_star = x(1:NbBr)
    alpha  = x(NbBr+1:NbBr+NbOut)
    dyn%p_star_prev = p_star
    dyn%alpha_prev(1:NbOut) = alpha

    ! Incoming quantities as a function of p_star(j) 

    ! Approximate jump relations through the rarefaction wave
    vitStar(i)=dyn%br(i)%vit0-(p_star(i)-dyn%br(i)%p0-dyn%br(i)%rho0*(dyn%br(i)%q0 &
        -dyn%br(i)%qBar) + dyn%br(i)%dx0*dyn%br(i)%Frc/2.0_dp )/(dyn%br(i)%rho0 &
        *(dyn%br(i)%vit0-dyn%br(i)%SpeedS0))
    rhoStar(i)=dyn%br(i)%rho0*(dyn%br(i)%vit0-dyn%br(i)%SpeedS0)/(vitStar(i)-dyn%br(i)%SpeedS0)

    if(R_Correction) then
        rCorrStar(i)=(dyn%br(i)%rCorr0*(dyn%br(i)%vit0-dyn%br(i)%SpeedS0) &
            -dyn%br(i)%dx0*dyn%br(i)%FrcR/2.0_dp )/(vitStar(i)-dyn%br(i)%SpeedS0)
        ETotStar(i)=(rCorrStar(i)-p_star(i))/rhoStar(i)+0.5_dp*vitStar(i)**2 &
            +(dyn%br(i)%q0-dyn%br(i)%qBar)
    else
        ETotStar(i)=dyn%br(i)%ETot0+(dyn%br(i)%p0*dyn%br(i)%vit0-p_star(i) &
            *vitStar(i))/(dyn%br(i)%rho0*(dyn%br(i)%vit0-dyn%br(i)%SpeedS0)) &
            +(dyn%br(i)%q0-dyn%br(i)%qBar)
    endif
    eStar(i)=ETotStar(i)-0.5_dp*vitStar(i)**2
    EnthalpyStar(i)=ETotStar(i)+p_star(i)/rhoStar(i)

    EnthalpyMixStarNum = sum(dyn%br%AreaJGB * rhoStar * vitStar * EnthalpyStar, im)
    denominator = sum(dyn%br%AreaJGB * rhoStar * vitStar, im)

    if(abs(denominator)<1.0e-10_dp) denominator = 1.0e-10_dp
    EnthalpyMixStar = EnthalpyMixStarNum / denominator

    ! Outgoing quantities as a function of p_star(j) and correction term alpha(j)
    ! through the contact discontinuity

    ! Approximate jump relations through the shock wave
    vitStar(o)=dyn%br(o)%vit0-(p_star(o)-dyn%br(o)%p0-dyn%br(o)%rho0*(dyn%br(o)%q0 &
        -dyn%br(o)%qBar) + dyn%br(o)%dx0*dyn%br(o)%Frc/2.0_dp)/(dyn%br(o)%rho0 &
        *(dyn%br(o)%vit0-dyn%br(o)%SpeedS0))
    rhoStar_Int(o)=dyn%br(o)%rho0*(dyn%br(o)%vit0-dyn%br(o)%SpeedS0)/(vitStar(o) &
        -dyn%br(o)%SpeedS0)
    ETotStar_Int(o)=dyn%br(o)%ETot0+(dyn%br(o)%p0*dyn%br(o)%vit0-p_star(o)*vitStar(o)) &
        /(dyn%br(o)%rho0*(dyn%br(o)%vit0-dyn%br(o)%SpeedS0))+(dyn%br(o)%q0-dyn%br(o)%qBar)

    QDMStar_Int(o)=rhoStar_Int(o)*vitStar(o)
    EnergyStar_Int(o)=rhoStar_Int(o)*ETotStar_Int(o)

    ! Correction through the contact discontinuity
    ! Additional term related to the real gas EOS
    do jj = 1, NbOut
        dFdRo(o(jj)) = droeint_droP(rhoStar_Int(o(jj)), p_star(o(jj)))
    enddo
    rhoStar(o)=rhoStar_Int(o)+alpha(o)
    QDMStar(o)=QDMStar_Int(o)+alpha(o)*vitStar(o)
    EnergyStar(o)=EnergyStar_Int(o)+alpha(o)*(dFdRo(o)+0.5_dp*vitStar(o)*vitStar(o))
    ETotStar(o)=EnergyStar(o)/rhoStar(o)

    EnthalpyStar(o)=(EnergyStar(o)+p_star(o))/rhoStar(o)

    eStar(o)=ETotStar(o)-0.5_dp*vitStar(o)**2
    !do j=1,NbOut
        if(R_Correction) then
            do jj = 1, NbOut
                T_temp(o(jj)) = T_roE(rhoStar(o(jj)), eStar(o(jj)))
                rCorrStar(o(jj)) = r_roT(rhoStar(o(jj)), T_temp(o(jj)))
            enddo
        endif
    !enddo

    !            Derivative computation for implicit scheme
    !..........................................................................!
    ! Jacobian matrix evaluation for the computation of derivatives of p_star
    ! and alpha with respect to 0-variables (all branches involved)
    do j=1,NbBr
        DerVitStar_dX(j)=-1.0_dp/(dyn%br(j)%rho0*(dyn%br(j)%vit0-dyn%br(j)%SpeedS0))
        DerRhoStar_dX(j)=-(dyn%br(j)%rho0*(dyn%br(j)%vit0-dyn%br(j)%SpeedS0)) &
            *DerVitStar_dX(j)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)

        if(j<=NbOut) then
            DerETotStar_Int_dX(j)=-(vitStar(j)+p_star(j)*DerVitStar_dX(j)) &
                /(dyn%br(j)%rho0*(dyn%br(j)%vit0-dyn%br(j)%SpeedS0))
            Der_dFdRo_dRhoInt(j)=(droeint_droP(rhoStar_Int(j)+1.0e-2_dp,p_star(j)) &
                -droeint_droP(rhoStar_Int(j),p_star(j)))/1.0e-2_dp         
            Der_dFdRo_dP(j)=(droeint_droP(rhoStar_Int(j),p_star(j)+1.0e-1_dp) &
                -droeint_droP(rhoStar_Int(j),p_star(j)))/1.0e-1_dp
            DerRhoETotStar_dX(j)=DerRhoStar_dX(j)*ETotStar_Int(j)+rhoStar_Int(j) &
                *DerETotStar_Int_dX(j)+alpha(j)*(Der_dFdRo_dRhoInt(j)*DerRhoStar_dX(j) &
                +Der_dFdRo_dP(j)+vitStar(j)*DerVitStar_dX(j))
            DerHStar_dX(j)=((DerRhoETotStar_dX(j)+1.0_dp)*rhoStar(j)-(rhoStar(j)*ETotStar(j) &
                +p_star(j))*DerRhoStar_dX(j))/(rhoStar(j)**2)
            DerRhoETotStar_dAlpha(j)=(dFdRo(j)+0.5_dp*vitStar(j)*vitStar(j))
            DerHStar_dAlpha(j)=(DerRhoETotStar_dAlpha(j)*rhoStar(j)-(rhoStar(j)*ETotStar(j) &
                +p_star(j)))/(rhoStar(j)**2)
        else
            if(R_Correction) then
                DerRCorrStar_dX(j)=-(dyn%br(j)%rCorr0*(dyn%br(j)%vit0-dyn%br(j)%SpeedS0) &
                    -dyn%br(j)%dx0*dyn%br(j)%FrcR/2.0_dp )*DerVitStar_dX(j)/((vitStar(j) &
                    -dyn%br(j)%SpeedS0)**2)
                DerETotStar_dX(j)=((DerRCorrStar_dX(j)-1.0_dp)*rhoStar(j)-(rCorrStar(j) &
                    -p_star(j))*DerRhoStar_dX(j))/(rhoStar(j)**2)+vitStar(j)*DerVitStar_dX(j)
            else
                DerETotStar_dX(j)=-(vitStar(j)+p_star(j)*DerVitStar_dX(j))/(dyn%br(j)%rho0 &
                    *(dyn%br(j)%vit0-dyn%br(j)%SpeedS0))
            endif
            DerHStar_dX(j)=DerETotStar_dX(j)+(rhoStar(j)-p_star(j)*DerRhoStar_dX(j))/(rhoStar(j)**2)
        endif
    enddo
    
    DerEntMixStar_dX(1:NbOut)=0.0_dp
    do j=NbOut+1,NbBr
        DerEntMixStar_dX(j)=((dyn%br(j)%AreaJGB*(DerRhoStar_dX(j)*vitStar(j)+rhoStar(j) &
            *DerVitStar_dX(j))*EnthalpyStar(j)+dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j) &
            *DerHStar_dX(j))*denominator-EnthalpyMixStarNum*dyn%br(j)%AreaJGB*(DerRhoStar_dX(j) &
            *vitStar(j)+rhoStar(j)*DerVitStar_dX(j)))/(denominator**2)
    enddo
  
    fJac = 0.0_dp
    do j=1,NbBr+NbOut 
        if(j<=NbBr) then
            fJac(j,1)=dyn%br(j)%AreaJGB*(DerRhoStar_dX(j)*vitStar(j)+rhoStar(j)*DerVitStar_dX(j))
        else
            fJac(j,1)=dyn%br(j-NbBr)%AreaJGB*vitStar(j-NbBr)
        endif
    enddo

    do j=1,NbOut  ! Outgoing pipes
        if(abs(rhoStar(j)*vitStar(j))<epsCoef.or.abs(rhoStar(NbOut+1)*vitStar(NbOut+1))<epsCoef)then
            fJac(j,j+1)=-1.0_dp
            fJac(NbOut+1,j+1)=1.0_dp
            fJac(NbBr+j,j+1)=0.0_dp
        else
            qj=-dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j)/(dyn%br(NbOut+1)%AreaJGB &
                *rhoStar(NbOut+1)*vitStar(NbOut+1))
            ps_ij=dyn%br(NbOut+1)%AreaJGB/dyn%br(j)%AreaJGB

            Der_qjdX=-(dyn%br(j)%AreaJGB/(dyn%br(NbOut+1)%AreaJGB*rhoStar(NbOut+1) &
                *vitStar(NbOut+1)))*(DerRhoStar_dX(j)*vitStar(j)+rhoStar(j)*DerVitStar_dX(j))        
            Der_qRef_dX=(dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j)/dyn%br(NbOut+1)%AreaJGB) &
                *(DerRhoStar_dX(NbOut+1)*vitStar(NbOut+1)+rhoStar(NbOut+1)*DerVitStar_dX(NbOut+1)) &
                /((rhoStar(NbOut+1)*vitStar(NbOut+1))**2)                   
            CoefPLoss(j)=1.0_dp-cos((3.0_dp/4.0_dp)*(3.14159_dp-dyn%br(j)%ThetaJGB))/(qj*ps_ij)
            DerCoef_dX(j)=Der_qjdX*cos((3.0_dp/4.0_dp)*(3.14159_dp-dyn%br(j)%ThetaJGB)) &
                /(ps_ij*qj**2)
            DerCoef_dX(NbOut+1)=Der_qRef_dX*cos((3.0_dp/4.0_dp)*(3.14159_dp-dyn%br(j)%ThetaJGB)) &
                /(ps_ij*qj**2)
            Der_qjdAlpha=-dyn%br(j)%AreaJGB*vitStar(j)/(dyn%br(NbOut+1)%AreaJGB*rhoStar(NbOut+1) &
                *vitStar(NbOut+1))
            DerCoef_dAlpha=Der_qjdAlpha*cos((3.0_dp/4.0_dp)*(3.14159_dp-dyn%br(j)%ThetaJGB)) &
                /(ps_ij*qj**2)
            if(NbOut+1==dyn%idxKappaRef) then
                CoefPLoss(j)=CoefPLoss(j)*dyn%kappa
                DerCoef_dX(j)=DerCoef_dX(j)*dyn%kappa
                DerCoef_dX(NbOut+1)=DerCoef_dX(NbOut+1)*dyn%kappa
                DerCoef_dAlpha=DerCoef_dAlpha*dyn%kappa
            else
                if(j==dyn%idxKappaRef) then
                    CoefPLoss(j)=CoefPLoss(j)*dyn%kappa
                    DerCoef_dX(j)=DerCoef_dX(j)*dyn%kappa
                    DerCoef_dX(NbOut+1)=DerCoef_dX(NbOut+1)*dyn%kappa
                    DerCoef_dAlpha=DerCoef_dAlpha*dyn%kappa     
                endif
            endif

            fJac(j,j+1)=-1.0_dp-DerCoef_dX(j)*rhoStar(j)*vitStar(j)**2-CoefPLoss(j) &
                *(DerRhoStar_dX(j)*vitStar(j)**2+2.0_dp*rhoStar(j)*vitStar(j)*DerVitStar_dX(j))
            fJac(NbOut+1,j+1)=1.0_dp-rhoStar(j)*(vitStar(j)**2)*DerCoef_dX(NbOut+1)
            fJac(NbBr+j,j+1)=-vitStar(j)**2*(CoefPLoss(j)+rhoStar(j)*DerCoef_dAlpha)
        endif
    enddo

    do j=NbOut+2,NbBr
        fJac(j,j)=-1.0_dp
        fJac(NbOut+1,j)=1.0_dp
    enddo

    do j=1,NbOut
        fJac(j,NbBr+j)=DerHStar_dX(j)

        fJac(NbBr+j,NbBr+j)=DerHStar_dAlpha(j)

        do jj=NbOut+1,NbBr
            fJac(jj,NbBr+j)=-DerEntMixStar_dX(jj) 
        enddo
    enddo    
   
    !..........................................................................!

    ! Jacobian inverse
    invJac(:,:)=inv(fJac)
    if (sim_error > 0) return

    do j=1,NbOut
        DerCoef(:,j)%dU1=0.0_dp
        DerCoef(:,j)%dU2=0.0_dp
        DerCoef(:,j)%dU3=0.0_dp
        DerCoef(:,j)%dU4=0.0_dp
    enddo

    do j=1,NbBr
        if(R_Correction) then
            BB=p_star(j)-(dyn%br(j)%U4-dyn%br(j)%U3+0.5_dp*(dyn%br(j)%U2**2)/dyn%br(j)%U1)
            BB=BB-dyn%br(j)%U1*(dyn%br(j)%q0-dyn%br(j)%qBar)    + dyn%br(j)%dx0*dyn%br(j)%Frc/2.0_dp 
            DerBBdUk(1)=0.5_dp*(dyn%br(j)%U2**2)/(dyn%br(j)%U1**2)-(dyn%br(j)%q0-dyn%br(j)%qBar) &
                                                           + dyn%br(j)%dx0*dyn%br(j)%DerFric1/2.0_dp 
            DerBBdUk(2)=-dyn%br(j)%U2/dyn%br(j)%U1         + dyn%br(j)%dx0*dyn%br(j)%DerFric2/2.0_dp
            DerBBdUk(3)=1.0_dp                             + dyn%br(j)%dx0*dyn%br(j)%DerFric3/2.0_dp
            DerBBdUk(4)=-1.0_dp                            + dyn%br(j)%dx0*dyn%br(j)%DerFric4/2.0_dp
        else
            HH=dyn%br(j)%ETot0+dyn%br(j)%p0/dyn%br(j)%rho0
            call jacobian_roT(dyn%br(j)%rho0, dyn%br(j)%T0, dedT, dPdT)
            dPde=dPdT/dedT
            kk=dPde/dyn%br(j)%rho0
            KKK=dyn%br(j)%CSound0**2+kk*(dyn%br(j)%vit0**2-HH)

            BB=p_star(j)-dyn%br(j)%p0
            BB=BB-dyn%br(j)%U1*(dyn%br(j)%q0-dyn%br(j)%qBar)    + dyn%br(j)%dx0*dyn%br(j)%Frc/2.0_dp 
            DerBBdUk(1)=-KKK-(dyn%br(j)%q0-dyn%br(j)%qBar) + dyn%br(j)%dx0*dyn%br(j)%DerFric1/2.0_dp 
            DerBBdUk(2)=kk*dyn%br(j)%vit0                  + dyn%br(j)%dx0*dyn%br(j)%DerFric2/2.0_dp  
            DerBBdUk(3)=-kk                                + dyn%br(j)%dx0*dyn%br(j)%DerFric3/2.0_dp  
            DerBBdUk(4)=0.0_dp                             + dyn%br(j)%dx0*dyn%br(j)%DerFric4/2.0_dp  
        endif

        DerVitStar_dU1(j)=(-dyn%br(j)%U2/(dyn%br(j)%U1**2)-(DerBBdUk(1)*(dyn%br(j)%U2 &
            -dyn%br(j)%U1*dyn%br(j)%SpeedS0)+BB*dyn%br(j)%SpeedS0)/((dyn%br(j)%U2 &
            -dyn%br(j)%U1*dyn%br(j)%SpeedS0)**2))
        DerVitStar_dU2(j)=(1.0_dp/dyn%br(j)%U1-(DerBBdUk(2)*(dyn%br(j)%U2-dyn%br(j)%U1&
            *dyn%br(j)%SpeedS0)-BB)/((dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)**2))
        DerVitStar_dU3(j)=(-DerBBdUk(3)/(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0))
        DerVitStar_dU4(j)=(-DerBBdUk(4)/(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0))

        DerRhoStar_dU1(j)=(-dyn%br(j)%SpeedS0*(vitStar(j)-dyn%br(j)%SpeedS0)-&
            (dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)*DerVitStar_dU1(j)) &
            /((vitStar(j)-dyn%br(j)%SpeedS0)**2)
        DerRhoStar_dU2(j)=((vitStar(j)-dyn%br(j)%SpeedS0)-(dyn%br(j)%U2-dyn%br(j)%U1*&
            dyn%br(j)%SpeedS0)*DerVitStar_dU2(j))/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
        DerRhoStar_dU3(j)=-(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0) &
            *DerVitStar_dU3(j)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
        DerRhoStar_dU4(j)=-(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0) &
            *DerVitStar_dU4(j)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)

        if(j<=NbOut) then
            if(R_Correction) then
                DD=dyn%br(j)%U4-dyn%br(j)%U3+0.5_dp*(dyn%br(j)%U2**2)/dyn%br(j)%U1
                DerDDdUk(1)=-0.5_dp*(dyn%br(j)%U2**2)/(dyn%br(j)%U1**2)
                DerDDdUk(2)=dyn%br(j)%U2/dyn%br(j)%U1
                DerDDdUk(3)=-1.0_dp
                DerDDdUk(4)=1.0_dp
            else
                HH=dyn%br(j)%ETot0+dyn%br(j)%p0/dyn%br(j)%rho0
                call jacobian_roT(dyn%br(j)%rho0, dyn%br(j)%T0, dedT, dPdT)
                dPde=dPdT/dedT
                kk=dPde/dyn%br(j)%rho0
                KKK=dyn%br(j)%CSound0**2+kk*(dyn%br(j)%vit0**2-HH)
                
                DD=dyn%br(j)%p0
                DerDDdUk(1)=KKK
                DerDDdUk(2)=-kk*dyn%br(j)%vit0
                DerDDdUk(3)=kk
                DerDDdUk(4)=0.0_dp
            endif
    
            DerETotStar_Int_dU1(j)=-dyn%br(j)%U3/(dyn%br(j)%U1**2)+(((DerDDdUk(1) &
                *dyn%br(j)%U2/dyn%br(j)%U1-DD*dyn%br(j)%U2/(dyn%br(j)%U1**2)) &
                -p_star(j)*DerVitStar_dU1(j))*(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0) &
                +(DD*dyn%br(j)%U2/dyn%br(j)%U1-p_star(j)*vitStar(j))*dyn%br(j)%SpeedS0) &
                /((dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)**2)
            DerETotStar_Int_dU2(j)=(((DerDDdUk(2)*dyn%br(j)%U2/dyn%br(j)%U1+DD/dyn%br(j)%U1) &
                -p_star(j)*DerVitStar_dU2(j))*(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0) &
                -(DD*dyn%br(j)%U2/dyn%br(j)%U1-p_star(j)*vitStar(j)))/((dyn%br(j)%U2 &
                -dyn%br(j)%U1*dyn%br(j)%SpeedS0)**2)
            DerETotStar_Int_dU3(j)=(1.0_dp/dyn%br(j)%U1)+(DerDDdUk(3)*dyn%br(j)%U2/dyn%br(j)%U1 &
                -p_star(j)*DerVitStar_dU3(j))/(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)
            DerETotStar_Int_dU4(j)=(DerDDdUk(4)*dyn%br(j)%U2/dyn%br(j)%U1-p_star(j) &
                *DerVitStar_dU4(j))/(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)

            Der_dFdRo_dRhoInt(j)=(droeint_droP(rhoStar_Int(j)+1.0e-2_dp,p_star(j))-&
                                droeint_droP(rhoStar_Int(j),p_star(j)))/1.0e-2_dp                        

            DerRhoETotStar_dUk(1)=DerRhoStar_dU1(j)*ETotStar_Int(j)+rhoStar_Int(j) &
                *DerETotStar_Int_dU1(j)+alpha(j)*(Der_dFdRo_dRhoInt(j)*DerRhoStar_dU1(j) &
                +vitStar(j)*DerVitStar_dU1(j))
            DerRhoETotStar_dUk(2)=DerRhoStar_dU2(j)*ETotStar_Int(j)+rhoStar_Int(j) &
                *DerETotStar_Int_dU2(j)+alpha(j)*(Der_dFdRo_dRhoInt(j)*DerRhoStar_dU2(j) &
                +vitStar(j)*DerVitStar_dU2(j))
            DerRhoETotStar_dUk(3)=DerRhoStar_dU3(j)*ETotStar_Int(j)+rhoStar_Int(j) &
                *DerETotStar_Int_dU3(j)+alpha(j)*(Der_dFdRo_dRhoInt(j)*DerRhoStar_dU3(j) &
                +vitStar(j)*DerVitStar_dU3(j))
            DerRhoETotStar_dUk(4)=DerRhoStar_dU4(j)*ETotStar_Int(j)+rhoStar_Int(j) &
                *DerETotStar_Int_dU4(j)+alpha(j)*(Der_dFdRo_dRhoInt(j)*DerRhoStar_dU4(j) &
                +vitStar(j)*DerVitStar_dU4(j))

            DerHStar_dU1(j)=(DerRhoETotStar_dUk(1)*rhoStar(j)-(rhoStar(j)*ETotStar(j) &
                +p_star(j))*DerRhoStar_dU1(j))/(rhoStar(j)**2)
            DerHStar_dU2(j)=(DerRhoETotStar_dUk(2)*rhoStar(j)-(rhoStar(j)*ETotStar(j) &
                +p_star(j))*DerRhoStar_dU2(j))/(rhoStar(j)**2)
            DerHStar_dU3(j)=(DerRhoETotStar_dUk(3)*rhoStar(j)-(rhoStar(j)*ETotStar(j) &
                +p_star(j))*DerRhoStar_dU3(j))/(rhoStar(j)**2)
            DerHStar_dU4(j)=(DerRhoETotStar_dUk(4)*rhoStar(j)-(rhoStar(j)*ETotStar(j) &
                +p_star(j))*DerRhoStar_dU4(j))/(rhoStar(j)**2)

        else
            if(R_Correction) then
                DerRCorrStar_dU1(j)=((-dyn%br(j)%U4*dyn%br(j)%U2/(dyn%br(j)%U1**2) &
                    -dyn%br(j)%dx0*dyn%br(j)%DerFricR1/2.0_dp)*(vitStar(j)-dyn%br(j)%SpeedS0) &
                    -(dyn%br(j)%U4*(dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0) &
                    -dyn%br(j)%dx0*dyn%br(j)%FrcR/2.0_dp)*DerVitStar_dU1(j)) &
                    /((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                DerRCorrStar_dU2(j)=((dyn%br(j)%U4/dyn%br(j)%U1-dyn%br(j)%dx0 &
                    *dyn%br(j)%DerFricR2/2.0_dp)*(vitStar(j)-dyn%br(j)%SpeedS0) &
                    -(dyn%br(j)%U4*(dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0) &
                    -dyn%br(j)%dx0*dyn%br(j)%FrcR/2.0_dp)*DerVitStar_dU2(j)) &
                    /((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                DerRCorrStar_dU3(j)=((-dyn%br(j)%dx0*dyn%br(j)%DerFricR3/2.0_dp)*(vitStar(j)&
                    -dyn%br(j)%SpeedS0)-(dyn%br(j)%U4*(dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0) &
                    -dyn%br(j)%dx0*dyn%br(j)%FrcR/2.0_dp)*DerVitStar_dU3(j))/((vitStar(j) &
                    -dyn%br(j)%SpeedS0)**2)
                DerRCorrStar_dU4(j)=((dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0 &
                    -dyn%br(j)%dx0*dyn%br(j)%DerFricR4/2.0_dp)*(vitStar(j)-dyn%br(j)%SpeedS0) &
                    -(dyn%br(j)%U4*(dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0) &
                    -dyn%br(j)%dx0*dyn%br(j)%FrcR/2.0_dp )*DerVitStar_dU4(j))/((vitStar(j) &
                    -dyn%br(j)%SpeedS0)**2)
            
                DerETotStar_dU1(j)=(DerRCorrStar_dU1(j)*rhoStar(j)-(rCorrStar(j)-p_star(j)) &
                    *DerRhoStar_dU1(j))/(rhoStar(j)**2)+vitStar(j)*DerVitStar_dU1(j)
                DerETotStar_dU2(j)=(DerRCorrStar_dU2(j)*rhoStar(j)-(rCorrStar(j)-p_star(j)) &
                    *DerRhoStar_dU2(j))/(rhoStar(j)**2)+vitStar(j)*DerVitStar_dU2(j)
                DerETotStar_dU3(j)=(DerRCorrStar_dU3(j)*rhoStar(j)-(rCorrStar(j)-p_star(j)) &
                    *DerRhoStar_dU3(j))/(rhoStar(j)**2)+vitStar(j)*DerVitStar_dU3(j)
                DerETotStar_dU4(j)=(DerRCorrStar_dU4(j)*rhoStar(j)-(rCorrStar(j)-p_star(j)) &
                    *DerRhoStar_dU4(j))/(rhoStar(j)**2)+vitStar(j)*DerVitStar_dU4(j)

                DerHStar_dU1(j)=DerETotStar_dU1(j)-p_star(j)*DerRhoStar_dU1(j)/(rhoStar(j)**2)
                DerHStar_dU2(j)=DerETotStar_dU2(j)-p_star(j)*DerRhoStar_dU2(j)/(rhoStar(j)**2)
                DerHStar_dU3(j)=DerETotStar_dU3(j)-p_star(j)*DerRhoStar_dU3(j)/(rhoStar(j)**2)
                DerHStar_dU4(j)=DerETotStar_dU4(j)-p_star(j)*DerRhoStar_dU4(j)/(rhoStar(j)**2)          
            else
                HH=dyn%br(j)%ETot0+dyn%br(j)%p0/dyn%br(j)%rho0
                call jacobian_roT(dyn%br(j)%rho0, dyn%br(j)%T0, dedT, dPdT)
                dPde=dPdT/dedT
                kk=dPde/dyn%br(j)%rho0
                KKK=dyn%br(j)%CSound0**2+kk*(dyn%br(j)%vit0**2-HH)
                
                DD=dyn%br(j)%p0
                DerDDdUk(1)=KKK
                DerDDdUk(2)=-kk*dyn%br(j)%vit0
                DerDDdUk(3)=kk
                DerDDdUk(4)=0.0_dp

                DerETotStar_dU1(j)=-dyn%br(j)%U3/(dyn%br(j)%U1**2)+(((DerDDdUk(1)*dyn%br(j)%U2 &
                    /dyn%br(j)%U1-DD*dyn%br(j)%U2/(dyn%br(j)%U1**2))-p_star(j)*DerVitStar_dU1(j)) &
                    *(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)+(DD*dyn%br(j)%U2/dyn%br(j)%U1 &
                    -p_star(j)*vitStar(j))*dyn%br(j)%SpeedS0)/((dyn%br(j)%U2-dyn%br(j)%U1 &
                    *dyn%br(j)%SpeedS0)**2)
                DerETotStar_dU2(j)=(((DerDDdUk(2)*dyn%br(j)%U2/dyn%br(j)%U1+DD/dyn%br(j)%U1) &
                    -p_star(j)*DerVitStar_dU2(j))*(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0) &
                    -(DD*dyn%br(j)%U2/dyn%br(j)%U1-p_star(j)*vitStar(j)))/((dyn%br(j)%U2 &
                    -dyn%br(j)%U1*dyn%br(j)%SpeedS0)**2)
                DerETotStar_dU3(j)=(1.0_dp/dyn%br(j)%U1)+(DerDDdUk(3)*dyn%br(j)%U2/dyn%br(j)%U1 &
                    -p_star(j)*DerVitStar_dU3(j))/(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)
                DerETotStar_dU4(j)=(DerDDdUk(4)*dyn%br(j)%U2/dyn%br(j)%U1-p_star(j) &
                    *DerVitStar_dU4(j))/(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)           

                DerHStar_dU1(j)=DerETotStar_dU1(j)-p_star(j)*DerRhoStar_dU1(j)/(rhoStar(j)**2)
                DerHStar_dU2(j)=DerETotStar_dU2(j)-p_star(j)*DerRhoStar_dU2(j)/(rhoStar(j)**2)
                DerHStar_dU3(j)=DerETotStar_dU3(j)-p_star(j)*DerRhoStar_dU3(j)/(rhoStar(j)**2)
                DerHStar_dU4(j)=DerETotStar_dU4(j)-p_star(j)*DerRhoStar_dU4(j)/(rhoStar(j)**2)
            endif
        endif
    enddo

    DerEntMixStar_dU1(1:NbOut)=0.0_dp
    DerEntMixStar_dU2(1:NbOut)=0.0_dp
    DerEntMixStar_dU3(1:NbOut)=0.0_dp
    DerEntMixStar_dU4(1:NbOut)=0.0_dp
    do j=NbOut+1,NbBr
      DerEntMixStar_dU1(j)=((dyn%br(j)%AreaJGB*(DerRhoStar_dU1(j)*vitStar(j) &
        +rhoStar(j)*DerVitStar_dU1(j))*EnthalpyStar(j)+dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j) &
        *DerHStar_dU1(j))*denominator-EnthalpyMixStarNum*dyn%br(j)%AreaJGB*(DerRhoStar_dU1(j) &
        *vitStar(j)+rhoStar(j)*DerVitStar_dU1(j)))/(denominator**2)
      DerEntMixStar_dU2(j)=((dyn%br(j)%AreaJGB*(DerRhoStar_dU2(j)*vitStar(j) &
        +rhoStar(j)*DerVitStar_dU2(j))*EnthalpyStar(j)+dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j) &
        *DerHStar_dU2(j))*denominator-EnthalpyMixStarNum*dyn%br(j)%AreaJGB*(DerRhoStar_dU2(j) &
        *vitStar(j)+rhoStar(j)*DerVitStar_dU2(j)))/(denominator**2)
      DerEntMixStar_dU3(j)=((dyn%br(j)%AreaJGB*(DerRhoStar_dU3(j)*vitStar(j) &
        +rhoStar(j)*DerVitStar_dU3(j))*EnthalpyStar(j)+dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j) &
        *DerHStar_dU3(j))*denominator-EnthalpyMixStarNum*dyn%br(j)%AreaJGB*(DerRhoStar_dU3(j) &
        *vitStar(j)+rhoStar(j)*DerVitStar_dU3(j)))/(denominator**2)
      DerEntMixStar_dU4(j)=((dyn%br(j)%AreaJGB*(DerRhoStar_dU4(j)*vitStar(j) &
        +rhoStar(j)*DerVitStar_dU4(j))*EnthalpyStar(j)+dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j) &
        *DerHStar_dU4(j))*denominator-EnthalpyMixStarNum*dyn%br(j)%AreaJGB*(DerRhoStar_dU4(j) &
        *vitStar(j)+rhoStar(j)*DerVitStar_dU4(j)))/(denominator**2)
    enddo

    do j=1,NbOut
      if(abs(rhoStar(j)*vitStar(j))<epsCoef.or.abs(rhoStar(NbOut+1)*vitStar(NbOut+1))<epsCoef) then
        CoefPLoss(j)=0.0_dp
        DerCoef(j,j)%dU1=0.0_dp
        DerCoef(j,j)%dU2=0.0_dp
        DerCoef(j,j)%dU3=0.0_dp
        DerCoef(j,j)%dU4=0.0_dp
        DerCoef(NbOut+1,j)%dU1=0.0_dp
        DerCoef(NbOut+1,j)%dU2=0.0_dp
        DerCoef(NbOut+1,j)%dU3=0.0_dp
        DerCoef(NbOut+1,j)%dU4=0.0_dp
      else
        qj=-dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j)/(dyn%br(NbOut+1)%AreaJGB*rhoStar(NbOut+1) &
            *vitStar(NbOut+1))
        ps_ij=dyn%br(NbOut+1)%AreaJGB/dyn%br(j)%AreaJGB

        Der_qj_dUk(1)=-(dyn%br(j)%AreaJGB/(dyn%br(NbOut+1)%AreaJGB*rhoStar(NbOut+1) &
            *vitStar(NbOut+1)))*(DerRhoStar_dU1(j)*vitStar(j)+rhoStar(j)*DerVitStar_dU1(j)) 
        Der_qj_dUk(2)=-(dyn%br(j)%AreaJGB/(dyn%br(NbOut+1)%AreaJGB*rhoStar(NbOut+1) &
            *vitStar(NbOut+1)))*(DerRhoStar_dU2(j)*vitStar(j)+rhoStar(j)*DerVitStar_dU2(j)) 
        Der_qj_dUk(3)=-(dyn%br(j)%AreaJGB/(dyn%br(NbOut+1)%AreaJGB*rhoStar(NbOut+1) &
            *vitStar(NbOut+1)))*(DerRhoStar_dU3(j)*vitStar(j)+rhoStar(j)*DerVitStar_dU3(j)) 
        Der_qj_dUk(4)=-(dyn%br(j)%AreaJGB/(dyn%br(NbOut+1)%AreaJGB*rhoStar(NbOut+1) &
            *vitStar(NbOut+1)))*(DerRhoStar_dU4(j)*vitStar(j)+rhoStar(j)*DerVitStar_dU4(j))           

        Der_qRef_dUk(1)=(dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j)/dyn%br(NbOut+1)%AreaJGB) &
            *(DerRhoStar_dU1(NbOut+1)*vitStar(NbOut+1)+rhoStar(NbOut+1)*DerVitStar_dU1(NbOut+1)) &
            /((rhoStar(NbOut+1)*vitStar(NbOut+1))**2)
        Der_qRef_dUk(2)=(dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j)/dyn%br(NbOut+1)%AreaJGB) &
            *(DerRhoStar_dU2(NbOut+1)*vitStar(NbOut+1)+rhoStar(NbOut+1)*DerVitStar_dU2(NbOut+1)) &
            /((rhoStar(NbOut+1)*vitStar(NbOut+1))**2)
        Der_qRef_dUk(3)=(dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j)/dyn%br(NbOut+1)%AreaJGB) &
            *(DerRhoStar_dU3(NbOut+1)*vitStar(NbOut+1)+rhoStar(NbOut+1)*DerVitStar_dU3(NbOut+1)) &
            /((rhoStar(NbOut+1)*vitStar(NbOut+1))**2)
        Der_qRef_dUk(4)=(dyn%br(j)%AreaJGB*rhoStar(j)*vitStar(j)/dyn%br(NbOut+1)%AreaJGB) &
            *(DerRhoStar_dU4(NbOut+1)*vitStar(NbOut+1)+rhoStar(NbOut+1)*DerVitStar_dU4(NbOut+1)) &
            /((rhoStar(NbOut+1)*vitStar(NbOut+1))**2)                      

        CoefPLoss(j)=1.0_dp-cos((3.0_dp/4.0_dp)*(3.14159_dp-dyn%br(j)%ThetaJGB))/(qj*ps_ij)
        DerCoef(j,j)%dU1=Der_qj_dUk(1)*cos((3.0_dp/4.0_dp)*(3.14159_dp &
            -dyn%br(j)%ThetaJGB))/(ps_ij*qj**2)
        DerCoef(j,j)%dU2=Der_qj_dUk(2)*cos((3.0_dp/4.0_dp)*(3.14159_dp &
            -dyn%br(j)%ThetaJGB))/(ps_ij*qj**2)
        DerCoef(j,j)%dU3=Der_qj_dUk(3)*cos((3.0_dp/4.0_dp)*(3.14159_dp &
            -dyn%br(j)%ThetaJGB))/(ps_ij*qj**2)
        DerCoef(j,j)%dU4=Der_qj_dUk(4)*cos((3.0_dp/4.0_dp)*(3.14159_dp &
            -dyn%br(j)%ThetaJGB))/(ps_ij*qj**2)
        DerCoef(NbOut+1,j)%dU1=Der_qRef_dUk(1)*cos((3.0_dp/4.0_dp) &
            *(3.14159_dp-dyn%br(j)%ThetaJGB))/(ps_ij*qj**2)
        DerCoef(NbOut+1,j)%dU2=Der_qRef_dUk(2)*cos((3.0_dp/4.0_dp) &
            *(3.14159_dp-dyn%br(j)%ThetaJGB))/(ps_ij*qj**2)
        DerCoef(NbOut+1,j)%dU3=Der_qRef_dUk(3)*cos((3.0_dp/4.0_dp) &
            *(3.14159_dp-dyn%br(j)%ThetaJGB))/(ps_ij*qj**2)
        DerCoef(NbOut+1,j)%dU4=Der_qRef_dUk(4)*cos((3.0_dp/4.0_dp) &
            *(3.14159_dp-dyn%br(j)%ThetaJGB))/(ps_ij*qj**2)

        if(NbOut+1==dyn%idxKappaRef) then
          CoefPLoss(j)=CoefPLoss(j)*dyn%kappa
          DerCoef(j,j)%dU1=DerCoef(j,j)%dU1*dyn%kappa
          DerCoef(j,j)%dU2=DerCoef(j,j)%dU2*dyn%kappa
          DerCoef(j,j)%dU3=DerCoef(j,j)%dU3*dyn%kappa
          DerCoef(j,j)%dU4=DerCoef(j,j)%dU4*dyn%kappa
          DerCoef(NbOut+1,j)%dU1=DerCoef(NbOut+1,j)%dU1*dyn%kappa
          DerCoef(NbOut+1,j)%dU2=DerCoef(NbOut+1,j)%dU2*dyn%kappa
          DerCoef(NbOut+1,j)%dU3=DerCoef(NbOut+1,j)%dU3*dyn%kappa
          DerCoef(NbOut+1,j)%dU4=DerCoef(NbOut+1,j)%dU4*dyn%kappa
        else
          if(j==dyn%idxKappaRef) then
            CoefPLoss(j)=CoefPLoss(j)*dyn%kappa
            DerCoef(j,j)%dU1=DerCoef(j,j)%dU1*dyn%kappa
            DerCoef(j,j)%dU2=DerCoef(j,j)%dU2*dyn%kappa
            DerCoef(j,j)%dU3=DerCoef(j,j)%dU3*dyn%kappa
            DerCoef(j,j)%dU4=DerCoef(j,j)%dU4*dyn%kappa
            DerCoef(NbOut+1,j)%dU1=DerCoef(NbOut+1,j)%dU1*dyn%kappa
            DerCoef(NbOut+1,j)%dU2=DerCoef(NbOut+1,j)%dU2*dyn%kappa
            DerCoef(NbOut+1,j)%dU3=DerCoef(NbOut+1,j)%dU3*dyn%kappa
            DerCoef(NbOut+1,j)%dU4=DerCoef(NbOut+1,j)%dU4*dyn%kappa              
          endif
        endif

      endif
    enddo

    do j=1,NbBr  
        DerFct(j,1)%dU1=dyn%br(j)%AreaJGB*(DerRhoStar_dU1(j)*vitStar(j)+rhoStar(j)*DerVitStar_dU1(j))   
        DerFct(j,1)%dU2=dyn%br(j)%AreaJGB*(DerRhoStar_dU2(j)*vitStar(j)+rhoStar(j)*DerVitStar_dU2(j))
        DerFct(j,1)%dU3=dyn%br(j)%AreaJGB*(DerRhoStar_dU3(j)*vitStar(j)+rhoStar(j)*DerVitStar_dU3(j))
        DerFct(j,1)%dU4=dyn%br(j)%AreaJGB*(DerRhoStar_dU4(j)*vitStar(j)+rhoStar(j)*DerVitStar_dU4(j))      
    enddo

    do j=1,NbOut  ! Outgoing pipes
        DerFct(j,j+1)%dU1=-DerCoef(j,j)%dU1*rhoStar(j)*vitStar(j)**2-CoefPLoss(j)*&
                            (DerRhoStar_dU1(j)*vitStar(j)**2+2.0_dp*rhoStar(j)*vitStar(j)*&
                            DerVitStar_dU1(j))
        DerFct(j,j+1)%dU2=-DerCoef(j,j)%dU2*rhoStar(j)*vitStar(j)**2-CoefPLoss(j)*&
                            (DerRhoStar_dU2(j)*vitStar(j)**2+2.0_dp*rhoStar(j)*vitStar(j)*&
                            DerVitStar_dU2(j))
        DerFct(j,j+1)%dU3=-DerCoef(j,j)%dU3*rhoStar(j)*vitStar(j)**2-CoefPLoss(j)*&
                            (DerRhoStar_dU3(j)*vitStar(j)**2+2.0_dp*rhoStar(j)*vitStar(j)*&
                            DerVitStar_dU3(j))
        DerFct(j,j+1)%dU4=-DerCoef(j,j)%dU4*rhoStar(j)*vitStar(j)**2-CoefPLoss(j)*&
                            (DerRhoStar_dU4(j)*vitStar(j)**2+2.0_dp*rhoStar(j)*vitStar(j)*&
                            DerVitStar_dU4(j))

        DerFct(NbOut+1,j+1)%dU1=-rhoStar(j)*(vitStar(j)**2)*DerCoef(NbOut+1,j)%dU1
        DerFct(NbOut+1,j+1)%dU2=-rhoStar(j)*(vitStar(j)**2)*DerCoef(NbOut+1,j)%dU2
        DerFct(NbOut+1,j+1)%dU3=-rhoStar(j)*(vitStar(j)**2)*DerCoef(NbOut+1,j)%dU3
        DerFct(NbOut+1,j+1)%dU4=-rhoStar(j)*(vitStar(j)**2)*DerCoef(NbOut+1,j)%dU4
    enddo

    do j=NbOut+2,NbBr  ! Incoming pipes
        DerFct(j,j)%dU1=0.0_dp
        DerFct(j,j)%dU2=0.0_dp
        DerFct(j,j)%dU3=0.0_dp
        DerFct(j,j)%dU4=0.0_dp

        DerFct(NbOut+1,j)%dU1=0.0_dp
        DerFct(NbOut+1,j)%dU2=0.0_dp
        DerFct(NbOut+1,j)%dU3=0.0_dp
        DerFct(NbOut+1,j)%dU4=0.0_dp
    enddo

    do j=1,NbOut
        DerFct(j,NbBr+j)%dU1=DerHStar_dU1(j)
        DerFct(j,NbBr+j)%dU2=DerHStar_dU2(j)
        DerFct(j,NbBr+j)%dU3=DerHStar_dU3(j)
        DerFct(j,NbBr+j)%dU4=DerHStar_dU4(j)

        do jj=NbOut+1,NbBr
            DerFct(jj,NbBr+j)%dU1=-DerEntMixStar_dU1(jj)
            DerFct(jj,NbBr+j)%dU2=-DerEntMixStar_dU2(jj)
            DerFct(jj,NbBr+j)%dU3=-DerEntMixStar_dU3(jj)
            DerFct(jj,NbBr+j)%dU4=-DerEntMixStar_dU4(jj)
        enddo
    enddo

    do j=1,NbBr
        DerXStar(j,:)%dU1=-matmul(DerFct(j,:)%dU1,invJac(:,:))
        DerXStar(j,:)%dU2=-matmul(DerFct(j,:)%dU2,invJac(:,:))
        DerXStar(j,:)%dU3=-matmul(DerFct(j,:)%dU3,invJac(:,:))
        DerXStar(j,:)%dU4=-matmul(DerFct(j,:)%dU4,invJac(:,:))
    enddo

    do j=1,NbBr    !! Variables written in the junction referential
        DerPStarGlob(:,j)%dU1=DerXStar(:,j)%dU1
        DerPStarGlob(:,j)%dU2=DerXStar(:,j)%dU2
        DerPStarGlob(:,j)%dU3=DerXStar(:,j)%dU3
        DerPStarGlob(:,j)%dU4=DerXStar(:,j)%dU4
    enddo

    do j=1,NbOut    !! Variables written in the junction referential
        DerAlphaGlob(:,j)%dU1=DerXStar(:,NbBr+j)%dU1
        DerAlphaGlob(:,j)%dU2=DerXStar(:,NbBr+j)%dU2
        DerAlphaGlob(:,j)%dU3=DerXStar(:,NbBr+j)%dU3
        DerAlphaGlob(:,j)%dU4=DerXStar(:,NbBr+j)%dU4
    enddo

    do j=1,NbBr
        do jj=1,NbBr
            if(jj==j) then
                if(R_Correction) then
                    BB=p_star(j)-(dyn%br(j)%U4-dyn%br(j)%U3+0.5_dp*(dyn%br(j)%U2**2)/dyn%br(j)%U1)
                    BB=BB-dyn%br(j)%U1*(dyn%br(j)%q0-dyn%br(j)%qBar) &
                                                                + dyn%br(j)%dx0*dyn%br(j)%Frc/2.0_dp
                    DerBBdUk(1)=DerPStarGlob(j,j)%dU1+0.5_dp*(dyn%br(j)%U2**2)/(dyn%br(j)%U1**2) &
                        -(dyn%br(j)%q0-dyn%br(j)%qBar) +     dyn%br(j)%dx0*dyn%br(j)%DerFric1/2.0_dp
                    DerBBdUk(2)=DerPStarGlob(j,j)%dU2-dyn%br(j)%U2/dyn%br(j)%U1 &
                                                           + dyn%br(j)%dx0*dyn%br(j)%DerFric2/2.0_dp
                    DerBBdUk(3)=DerPStarGlob(j,j)%dU3+1.0_dp &
                                                           + dyn%br(j)%dx0*dyn%br(j)%DerFric3/2.0_dp
                    DerBBdUk(4)=DerPStarGlob(j,j)%dU4-1.0_dp &
                                                           + dyn%br(j)%dx0*dyn%br(j)%DerFric4/2.0_dp
                else
                    HH=dyn%br(j)%ETot0+dyn%br(j)%p0/dyn%br(j)%rho0
                    call jacobian_roT(dyn%br(j)%rho0, dyn%br(j)%T0, dedT, dPdT)
                    dPde=dPdT/dedT
                    kk=dPde/dyn%br(j)%rho0
                    KKK=dyn%br(j)%CSound0**2+kk*(dyn%br(j)%vit0**2-HH)

                    BB=p_star(j)-dyn%br(j)%p0
                    BB=BB-dyn%br(j)%U1*(dyn%br(j)%q0-dyn%br(j)%qBar) &
                                                                + dyn%br(j)%dx0*dyn%br(j)%Frc/2.0_dp
                    DerBBdUk(1)=DerPStarGlob(j,j)%dU1-KKK-(dyn%br(j)%q0-dyn%br(j)%qBar) &
                                                           + dyn%br(j)%dx0*dyn%br(j)%DerFric1/2.0_dp
                    DerBBdUk(2)=DerPStarGlob(j,j)%dU2+kk*dyn%br(j)%vit0 &
                                                           + dyn%br(j)%dx0*dyn%br(j)%DerFric2/2.0_dp
                    DerBBdUk(3)=DerPStarGlob(j,j)%dU3-kk   + dyn%br(j)%dx0*dyn%br(j)%DerFric3/2.0_dp
                    DerBBdUk(4)=DerPStarGlob(j,j)%dU4      + dyn%br(j)%dx0*dyn%br(j)%DerFric4/2.0_dp
                endif
            
                DerVelStarGlob(jj,j)%dU1=(-dyn%br(j)%U2/(dyn%br(j)%U1**2)-(DerBBdUk(1) &
                    *(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)+BB*dyn%br(j)%SpeedS0) &
                    /((dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)**2))
                DerVelStarGlob(jj,j)%dU2=(1.0_dp/dyn%br(j)%U1-(DerBBdUk(2)*(dyn%br(j)%U2 &
                    -dyn%br(j)%U1*dyn%br(j)%SpeedS0)-BB)/((dyn%br(j)%U2 &
                    -dyn%br(j)%U1*dyn%br(j)%SpeedS0)**2))
                DerVelStarGlob(jj,j)%dU3=(-DerBBdUk(3)/(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0))
                DerVelStarGlob(jj,j)%dU4=(-DerBBdUk(4)/(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0))

                DerRhoStarGlob(jj,j)%dU1=(-dyn%br(j)%SpeedS0*(vitStar(j)-dyn%br(j)%SpeedS0) &
                    -(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0)*DerVelStarGlob(jj,j)%dU1) &
                    /((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                DerRhoStarGlob(jj,j)%dU2=((vitStar(j)-dyn%br(j)%SpeedS0)-(dyn%br(j)%U2 &
                    -dyn%br(j)%U1*dyn%br(j)%SpeedS0)*DerVelStarGlob(jj,j)%dU2)/((vitStar(j) &
                    -dyn%br(j)%SpeedS0)**2)
                DerRhoStarGlob(jj,j)%dU3=-(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0) &
                    *DerVelStarGlob(jj,j)%dU3/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                DerRhoStarGlob(jj,j)%dU4=-(dyn%br(j)%U2-dyn%br(j)%U1*dyn%br(j)%SpeedS0) &
                    *DerVelStarGlob(jj,j)%dU4/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
            else
                DerVelStarGlob(jj,j)%dU1=-DerPStarGlob(jj,j)%dU1/(dyn%br(j)%U2 &
                    -dyn%br(j)%U1*dyn%br(j)%SpeedS0)
                DerVelStarGlob(jj,j)%dU2=-DerPStarGlob(jj,j)%dU2/(dyn%br(j)%U2 &
                    -dyn%br(j)%U1*dyn%br(j)%SpeedS0)
                DerVelStarGlob(jj,j)%dU3=-DerPStarGlob(jj,j)%dU3/(dyn%br(j)%U2 &
                    -dyn%br(j)%U1*dyn%br(j)%SpeedS0)
                DerVelStarGlob(jj,j)%dU4=-DerPStarGlob(jj,j)%dU4/(dyn%br(j)%U2 &
                    -dyn%br(j)%U1*dyn%br(j)%SpeedS0)

                DerRhoStarGlob(jj,j)%dU1=-DerVelStarGlob(jj,j)%dU1*(dyn%br(j)%U2-dyn%br(j)%U1 &
                    *dyn%br(j)%SpeedS0)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                DerRhoStarGlob(jj,j)%dU2=-DerVelStarGlob(jj,j)%dU2*(dyn%br(j)%U2-dyn%br(j)%U1 &
                    *dyn%br(j)%SpeedS0)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                DerRhoStarGlob(jj,j)%dU3=-DerVelStarGlob(jj,j)%dU3*(dyn%br(j)%U2-dyn%br(j)%U1 &
                    *dyn%br(j)%SpeedS0)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                DerRhoStarGlob(jj,j)%dU4=-DerVelStarGlob(jj,j)%dU4*(dyn%br(j)%U2-dyn%br(j)%U1 &
                    *dyn%br(j)%SpeedS0)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)           
                                                    
            endif
            if(j<=NbOut) then
                DerRhoStarGlob(jj,j)%dU1=DerRhoStarGlob(jj,j)%dU1+DerAlphaGlob(jj,j)%dU1
                DerRhoStarGlob(jj,j)%dU2=DerRhoStarGlob(jj,j)%dU2+DerAlphaGlob(jj,j)%dU2
                DerRhoStarGlob(jj,j)%dU3=DerRhoStarGlob(jj,j)%dU3+DerAlphaGlob(jj,j)%dU3
                DerRhoStarGlob(jj,j)%dU4=DerRhoStarGlob(jj,j)%dU4+DerAlphaGlob(jj,j)%dU4
            endif   
        enddo
    enddo

    if(R_Correction) then
        do j=1,NbBr
            do jj=1,NbBr
                if(j<=NbOut) then
                    T_temp(j)=T_roE(rhoStar(j),eStar(j))
                    call jacobian_roT(rhoStar(j), T_temp(j), dedT, dPdT, dTdp_Ro, dTdRo_p, dRStar_dRo, dRStar_dT, RStar)
                    DerRCorStarGlob(jj,j)%dU1=(dRStar_dRo+dRStar_dT*dTdRo_p) &
                        *DerRhoStarGlob(jj,j)%dU1+dRStar_dT*dTdp_Ro*DerPStarGlob(jj,j)%dU1
                    DerRCorStarGlob(jj,j)%dU2=(dRStar_dRo+dRStar_dT*dTdRo_p) &
                        *DerRhoStarGlob(jj,j)%dU2+dRStar_dT*dTdp_Ro*DerPStarGlob(jj,j)%dU2  
                    DerRCorStarGlob(jj,j)%dU3=(dRStar_dRo+dRStar_dT*dTdRo_p) &
                        *DerRhoStarGlob(jj,j)%dU3+dRStar_dT*dTdp_Ro*DerPStarGlob(jj,j)%dU3
                    DerRCorStarGlob(jj,j)%dU4=(dRStar_dRo+dRStar_dT*dTdRo_p) &
                        *DerRhoStarGlob(jj,j)%dU4+dRStar_dT*dTdp_Ro*DerPStarGlob(jj,j)%dU4
                else
                    if(j==jj) then
                        DerRCorStarGlob(jj,j)%dU1=((-dyn%br(j)%U4*dyn%br(j)%U2/(dyn%br(j)%U1**2) &
                            -dyn%br(j)%dx0*dyn%br(j)%DerFricR1/2.0_dp)*(vitStar(j) &
                            -dyn%br(j)%SpeedS0)-(dyn%br(j)%U4*(dyn%br(j)%U2/dyn%br(j)%U1 &
                            -dyn%br(j)%SpeedS0)-dyn%br(j)%dx0*dyn%br(j)%FrcR/2.0_dp) &
                            *DerVelStarGlob(jj,j)%dU1)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                        DerRCorStarGlob(jj,j)%dU2=((dyn%br(j)%U4/dyn%br(j)%U1-dyn%br(j)%dx0 &
                            *dyn%br(j)%DerFricR2/2.0_dp)*(vitStar(j)-dyn%br(j)%SpeedS0) &
                            -(dyn%br(j)%U4*(dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0) &
                            -dyn%br(j)%dx0*dyn%br(j)%FrcR/2.0_dp)*DerVelStarGlob(jj,j)%dU2) &
                            /((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                        DerRCorStarGlob(jj,j)%dU3=((-dyn%br(j)%dx0*dyn%br(j)%DerFricR3/2.0_dp) &
                            *(vitStar(j)-dyn%br(j)%SpeedS0)-(dyn%br(j)%U4*(dyn%br(j)%U2 &
                            /dyn%br(j)%U1-dyn%br(j)%SpeedS0)-dyn%br(j)%dx0*dyn%br(j)%FrcR/2.0_dp) &
                            *DerVelStarGlob(jj,j)%dU3)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                        DerRCorStarGlob(jj,j)%dU4=((dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0 &
                            -dyn%br(j)%dx0*dyn%br(j)%DerFricR4/2.0_dp)*(vitStar(j) &
                            -dyn%br(j)%SpeedS0)-(dyn%br(j)%U4*(dyn%br(j)%U2/dyn%br(j)%U1 &
                            -dyn%br(j)%SpeedS0) - dyn%br(j)%dx0*dyn%br(j)%FrcR/2.0_dp) &
                            *DerVelStarGlob(jj,j)%dU4)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                    else
                        DerRCorStarGlob(jj,j)%dU1=-DerVelStarGlob(jj,j)%dU1*(dyn%br(j)%U4 &
                            *(dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0)-dyn%br(j)%dx0 &
                            *dyn%br(j)%FrcR/2.0_dp)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                        DerRCorStarGlob(jj,j)%dU2=-DerVelStarGlob(jj,j)%dU2*(dyn%br(j)%U4 &
                            *(dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0)-dyn%br(j)%dx0 &
                            *dyn%br(j)%FrcR/2.0_dp)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                        DerRCorStarGlob(jj,j)%dU3=-DerVelStarGlob(jj,j)%dU3*(dyn%br(j)%U4 &
                            *(dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0)-dyn%br(j)%dx0 &
                            *dyn%br(j)%FrcR/2.0_dp)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)
                        DerRCorStarGlob(jj,j)%dU4=-DerVelStarGlob(jj,j)%dU4*(dyn%br(j)%U4 &
                            *(dyn%br(j)%U2/dyn%br(j)%U1-dyn%br(j)%SpeedS0)-dyn%br(j)%dx0 &
                            *dyn%br(j)%FrcR/2.0_dp)/((vitStar(j)-dyn%br(j)%SpeedS0)**2)  
                    endif
                endif        
            enddo
        enddo
    else
        do j=1,NbBr
            do jj=1,NbBr
                T_temp(j)=T_roE(rhoStar(j),eStar(j))
                call jacobian_roT(rhoStar(j), T_temp(j), dedT, dPdT, dTdp_Ro, dTdRo_p, dRStar_dRo, dRStar_dT, RStar)
                DerRPhysStarGlob(jj,j)%dU1=(dRStar_dRo+dRStar_dT*dTdRo_p) &
                    *DerRhoStarGlob(jj,j)%dU1+dRStar_dT*dTdp_Ro*DerPStarGlob(jj,j)%dU1
                DerRPhysStarGlob(jj,j)%dU2=(dRStar_dRo+dRStar_dT*dTdRo_p) &
                    *DerRhoStarGlob(jj,j)%dU2+dRStar_dT*dTdp_Ro*DerPStarGlob(jj,j)%dU2  
                DerRPhysStarGlob(jj,j)%dU3=(dRStar_dRo+dRStar_dT*dTdRo_p) &
                    *DerRhoStarGlob(jj,j)%dU3+dRStar_dT*dTdp_Ro*DerPStarGlob(jj,j)%dU3
                DerRPhysStarGlob(jj,j)%dU4=(dRStar_dRo+dRStar_dT*dTdRo_p) &
                    *DerRhoStarGlob(jj,j)%dU4+dRStar_dT*dTdp_Ro*DerPStarGlob(jj,j)%dU4
            enddo
        enddo
    endif

    do j=1,NbBr
        ConStar(j)%VarC(Con_Mas)=rhoStar(j)
        ConStar(j)%VarC(Con_Qdm)=rhoStar(j)*vitStar(j)
        ConStar(j)%VarC(Con_Ene)=rhoStar(j)*ETotStar(j)
        ConStar(j)%VarC(Con_R)=1.0_dp
        if(R_Correction) ConStar(j)%VarC(Con_R)=rCorrStar(j)

        if(R_Correction) then
            call jacobian_star_loc(ConStar(j)%VarC,Jacob_star(:,:,j))
        else
            PrimStar(j)%VarP(Pri_ro)=rhoStar(j)
            PrimStar(j)%VarP(Pri_u)=vitStar(j)
            PrimStar(j)%VarP(Pri_p)=p_star(j)
            call state_roP(rhoStar(j), p_star(j), PrimStar(j)%VarP(Pri_R), PrimStar(j)%VarP(Pri_e), &
                PrimStar(j)%VarP(Pri_T), PrimStar(j)%VarP(Pri_c))
            call jacobian_star_Req_rem_loc(PrimStar(j)%VarP,Jacob_star(:,:,j))
        endif

        dFdFJ(:,:)=0.0_dp
        dFdFJ(1,1)=dyn%br(j)%SgnJGB
        dFdFJ(2,2)=1.0_dp
        dFdFJ(3,3)=dyn%br(j)%SgnJGB
        dFdFJ(4,4)=dyn%br(j)%SgnJGB
      
        do jj=1,NbBr
            DerConStar(jj,j)%DerVarC(:,1)=[DerRhoStarGlob(jj,j)%dU1,&
                                    DerRhoStarGlob(jj,j)%dU2,&
                                    DerRhoStarGlob(jj,j)%dU3,&
                                    DerRhoStarGlob(jj,j)%dU4]

            DerConStar(jj,j)%DerVarC(:,2)=[&
                rhoStar(j)*DerVelStarGlob(jj,j)%dU1+DerRhoStarGlob(jj,j)%dU1*vitStar(j),&
                rhoStar(j)*DerVelStarGlob(jj,j)%dU2+DerRhoStarGlob(jj,j)%dU2*vitStar(j),&
                rhoStar(j)*DerVelStarGlob(jj,j)%dU3+DerRhoStarGlob(jj,j)%dU3*vitStar(j),&
                rhoStar(j)*DerVelStarGlob(jj,j)%dU4+DerRhoStarGlob(jj,j)%dU4*vitStar(j)]
        
            if(R_Correction) then
                DerConStar(jj,j)%DerVarC(:,3)=[&
                    DerRCorStarGlob(jj,j)%dU1+0.5_dp*(DerRhoStarGlob(jj,j)%dU1*vitStar(j)*vitStar(j)+&
                        2.0_dp*rhoStar(j)*vitStar(j)*DerVelStarGlob(jj,j)%dU1)-DerPStarGlob(jj,j)%dU1,&
                    DerRCorStarGlob(jj,j)%dU2+0.5_dp*(DerRhoStarGlob(jj,j)%dU2*vitStar(j)*vitStar(j)+&
                        2.0_dp*rhoStar(j)*vitStar(j)*DerVelStarGlob(jj,j)%dU2)-DerPStarGlob(jj,j)%dU2,&
                    DerRCorStarGlob(jj,j)%dU3+0.5_dp*(DerRhoStarGlob(jj,j)%dU3*vitStar(j)*vitStar(j)+&
                        2.0_dp*rhoStar(j)*vitStar(j)*DerVelStarGlob(jj,j)%dU3)-DerPStarGlob(jj,j)%dU3,&
                    DerRCorStarGlob(jj,j)%dU4+0.5_dp*(DerRhoStarGlob(jj,j)%dU4*vitStar(j)*vitStar(j)+&
                        2.0_dp*rhoStar(j)*vitStar(j)*DerVelStarGlob(jj,j)%dU4)-DerPStarGlob(jj,j)%dU4]                        
            else
                DerConStar(jj,j)%DerVarC(:,3)=[&
                    DerRPhysStarGlob(jj,j)%dU1+0.5_dp*(DerRhoStarGlob(jj,j)%dU1*vitStar(j)*vitStar(j)+&
                        2.0_dp*rhoStar(j)*vitStar(j)*DerVelStarGlob(jj,j)%dU1)-DerPStarGlob(jj,j)%dU1,&
                    DerRPhysStarGlob(jj,j)%dU2+0.5_dp*(DerRhoStarGlob(jj,j)%dU2*vitStar(j)*vitStar(j)+&
                        2.0_dp*rhoStar(j)*vitStar(j)*DerVelStarGlob(jj,j)%dU2)-DerPStarGlob(jj,j)%dU2,&
                    DerRPhysStarGlob(jj,j)%dU3+0.5_dp*(DerRhoStarGlob(jj,j)%dU3*vitStar(j)*vitStar(j)+&
                        2.0_dp*rhoStar(j)*vitStar(j)*DerVelStarGlob(jj,j)%dU3)-DerPStarGlob(jj,j)%dU3,&
                    DerRPhysStarGlob(jj,j)%dU4+0.5_dp*(DerRhoStarGlob(jj,j)%dU4*vitStar(j)*vitStar(j)+&
                        2.0_dp*rhoStar(j)*vitStar(j)*DerVelStarGlob(jj,j)%dU4)-DerPStarGlob(jj,j)%dU4]
            endif

            DerConStar(jj,j)%DerVarC(:,4)=0.0_dp
            if(R_Correction) DerConStar(jj,j)%DerVarC(:,4) = [&
                    DerRCorStarGlob(jj,j)%dU1,DerRCorStarGlob(jj,j)%dU2,&
                    DerRCorStarGlob(jj,j)%dU3,DerRCorStarGlob(jj,j)%dU4]          
        enddo

        ! Velocity star rewritten in the referential intrinsic to each branch
        vitStar(j)=dyn%br(j)%SgnJGB*vitStar(j)
        FlJ(j)%Vit=vitStar(j)
        FlJ(j)%Flx(Con_Mas)=rhoStar(j)*vitStar(j)
        FlJ(j)%Flx(Con_Qdm)=rhoStar(j)*vitStar(j)*vitStar(j)+p_star(j)
        FlJ(j)%Flx(Con_Ene)=(rhoStar(j)*ETotStar(j)+p_star(j))*vitStar(j)
        FlJ(j)%Flx(Con_R)=0.0_dp
        if(R_Correction) FlJ(j)%Flx(Con_R)=rCorrStar(j)*vitStar(j)

        do jj=1,NbBr
            dUdUjj(:,:)=0.0_dp
            dUdUjj(1,1)=1.0_dp
            dUdUjj(2,2)=dyn%br(jj)%SgnJGB
            dUdUjj(3,3)=1.0_dp
            dUdUjj(4,4)=1.0_dp

            FlJ(j)%DerFlx(:,:,jj) = matmul(dUdUjj(:,:), matmul(DerConStar(jj,j)%DerVarC(:,:),&
                    matmul(Jacob_star(:,:,j),dFdFJ(:,:))))
            if(.not. R_Correction) FlJ(j)%DerFlx(:,4,jj)=0.0_dp

            FlJ(j)%DerVit(1,jj)=dyn%br(j)%SgnJGB*DerVelStarGlob(jj,j)%dU1*dUdUjj(1,1)
            FlJ(j)%DerVit(2,jj)=dyn%br(j)%SgnJGB*DerVelStarGlob(jj,j)%dU2*dUdUjj(2,2)
            FlJ(j)%DerVit(3,jj)=dyn%br(j)%SgnJGB*DerVelStarGlob(jj,j)%dU3*dUdUjj(3,3)
            FlJ(j)%DerVit(4,jj)=dyn%br(j)%SgnJGB*DerVelStarGlob(jj,j)%dU4*dUdUjj(4,4)       
        enddo
    enddo

end subroutine solve_junction



function inv(A) result(A_inv)
    !! private called by junction solve
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(size(A,1),size(A,2)) :: A_inv

    real(dp), dimension(size(A,1)) :: work  ! work array for LaPACK
    integer, dimension(size(A,1)) :: i_piv  ! pivot indices
    integer :: n, info

    ! External procedures defined in LaPACK
    external DGETRF
    external DGETRI

    ! Store A in A_inv to prevent it from being overwritten by LaPACK
    A_inv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, A_inv, n, i_piv, info)

    if (info /= 0 ) then
        call set_error('matrix is numerically singular in junction inv(A)')
        return
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, A_inv, n, i_piv, work, n, info)

    if (info /= 0) then
        call set_error('matrix inversion failed in junction inv(A)')
        return
    end if
end function inv

  
subroutine jacobian_star_loc(U,Jacob_star)
    !! private called by junction solve
    real(dp), intent(in) :: U(:)
    real(dp), intent(out) :: Jacob_star(:,:)

    Jacob_star      = 0.0_dp
    Jacob_star(2,1) = 1.0_dp

    Jacob_star(1,2) = -3.0_dp / 2.0_dp * U(2)**2 / U(1)**2
    Jacob_star(2,2) =  3.0_dp * U(2) / U(1)
    Jacob_star(3,2) = -1.0_dp
    Jacob_star(4,2) =  1.0_dp

    Jacob_star(1,3) = -U(4) * U(2) / U(1)**2 - U(2)**3 / U(1)**3
    Jacob_star(2,3) =  U(4) / U(1) + 3.0_dp / 2.0_dp * U(2)**2 / U(1)**2
    Jacob_star(4,3) =  U(2) / U(1)
    
    Jacob_star(1,4) = -U(4) * U(2) / U(1)**2
    Jacob_star(2,4) =  U(4) / U(1)
    Jacob_star(4,4) =  U(2) / U(1)  
end subroutine jacobian_star_loc



subroutine jacobian_star_Req_rem_loc(PPP, Jacob_star)
    !! private called by junction solve
    real(dp), intent(in)  :: PPP(Nb_VarP)
    real(dp), intent(out) :: Jacob_star(Nb_VarC,Nb_VarC)

    real(dp) :: HH,kk,KKK,dPde,dPdT,dedT

    HH = PPP(Pri_e) + 0.5_dp * PPP(Pri_u)**2 + PPP(Pri_p) / PPP(Pri_ro)
    call jacobian_roT(PPP(Pri_ro), PPP(Pri_T), dedT, dPdT)
    dPde = dPdT / dedT
    kk   = dPde / PPP(Pri_ro)
    KKK  = PPP(Pri_c)**2 + kk * (PPP(Pri_u)**2 - HH)
      
    Jacob_star      = 0.0_dp
    Jacob_star(2,1) = 1.0_dp
      
    Jacob_star(1,2) = KKK - PPP(Pri_u)**2
    Jacob_star(2,2) = PPP(Pri_u) * (2.0_dp - kk)
    Jacob_star(3,2) = kk
      
    Jacob_star(1,3) = (KKK-HH) * PPP(Pri_u)
    Jacob_star(2,3) = HH - kk * PPP(Pri_u)**2
    Jacob_star(3,3) = PPP(Pri_u) * (1.0_dp + kk)          
end subroutine jacobian_star_Req_rem_loc


subroutine Junction_resolution_from_and_to_ports(me)
    !! Called at the beginning of main loop
    type(junction_t), intent(inout) :: me

    real(dp), parameter :: epsM = 1.0e-12_dp
    real(dp) :: m_dot(me%NbTotBr), angle(me%NbTotBr)
    integer :: i, i_ref, i_conv(me%NbTotBr), i_range(me%NbTotBr)
    logical :: in(me%NbTotBr), out(me%NbTotBr)

    ! calculation of the mass flow rate per branch to determine inlets, outlets
    in = .false.
    do i = 1, size(me%br); associate(prt => me%br(i)%p) 
        m_dot(i) = prt%Sgn4j * prt%Prim(Pri_u) * prt%Prim(Pri_ro) * prt%Area
        if (m_dot(i) > epsM) in(i) = .true. ! Identification of the incoming branches
    end associate; enddo
    if (any(ieee_is_nan(m_dot))) then; call set_error('NaN in junction m_dot (junction)'); return; end if
    i_ref = maxLoc(m_dot, 1) ! reference branch should be with highest mass flow rate
    if(m_dot(i_ref) <= epsM) in = .true. ! no inlet branches so make them all inlet

    ! calculate conversion index array
    out = .not.in
    in(i_ref) = .false. ! exclude reference branch from incoming
    i_range = [(i, i = 1, size(me%br))]
    i_conv = [pack(i_range, out), i_ref, pack(i_range, in)]

    angle = abs(me%br(i_ref)%angle - me%br%angle) / 180.0_dp * Pi_value
    
    me%dynamic%NbTotBr = size(me%br)
    me%dynamic%NbIn = count(in)+1
    me%dynamic%NbOut = count(out)
     me%dynamic%kappa = me%kappa
    me%dynamic%idxKappaRef = findLoc(i_conv, 1, 1)   
    do i = 1, size(me%br); associate(dyn => me%dynamic%br(i), prt => me%br(i_conv(i))%p)
        dyn%rho0 = prt%Prim(Pri_ro)
        dyn%vit0 = -prt%Sgn4j * prt%Prim(Pri_u)
        dyn%p0 = prt%Prim(Pri_p)
        dyn%ETot0 = prt%Prim(Pri_e) + 0.5_dp * prt%Prim(Pri_u) * prt%Prim(Pri_u)
        dyn%rCorr0 = prt%Prim(Pri_R)
        dyn%SpeedS0 = abs(prt%Prim(Pri_u)) + prt%Prim(Pri_c)
        dyn%CSound0 = prt%Prim(Pri_c)
        dyn%T0 = prt%Prim(Pri_T)
        dyn%AreaJGB = prt%Area
        dyn%ThetaJGB = angle(i_conv(i))
        dyn%SgnJGB = -prt%Sgn4j
        dyn%U1 = prt%Cons(Con_Mas)
        dyn%U2 = -prt%Sgn4j * prt%Cons(Con_Qdm)
        dyn%U3 = prt%Cons(Con_Ene)
        dyn%U4 = prt%Cons(Con_R)
        dyn%q0 = 0.0_dp
        dyn%qBar = 0.0_dp
        dyn%dx0 = prt%dxLoc
        dyn%Frc      = -prt%Sgn4j * prt%Su_fric
        dyn%DerFric1 = -prt%Sgn4j * prt%DerSu_fric(Con_Mas)
        dyn%DerFric2 =              prt%DerSu_fric(Con_Qdm)
        dyn%DerFric3 = -prt%Sgn4j * prt%DerSu_fric(Con_Ene)
        dyn%DerFric4 = -prt%Sgn4j * prt%DerSu_fric(Con_R)
        if(R_Correction) then
            dyn%FrcR      = prt%Sr_fric
            dyn%DerFricR1 = prt%DerSr_fric(Con_Mas)
            dyn%DerFricR2 = -prt%Sgn4j * prt%DerSr_fric(Con_Qdm)
            dyn%DerFricR3 = prt%DerSr_fric(Con_Ene)
            dyn%DerFricR4 = prt%DerSr_fric(Con_R)
        endif
    end associate; enddo

    call solve_junction(me%dynamic,me%FlJ)
    if (sim_error > 0) return

    do i = 1, size(me%br)
        me%DerFlx_dBr(:,:,i_conv,i_conv(i)) = me%FlJ(i)%DerFlx
        me%DerVel_dBr(:,  i_conv,i_conv(i)) = me%FlJ(i)%DerVit
    enddo

    do i = 1, size(me%br)
        me%br(i_conv(i))%p%flx = me%FlJ(i)%Flx
        me%br(i_conv(i))%p%vit = me%FlJ(i)%Vit
        me%br(i)%p%derFlx_derCon = me%DerFlx_dBr(:,:,i,i)
        me%br(i)%p%derVit_derCon = me%DerVel_dBr(:,i,i)
    enddo

    me%wave_time = min(minVal(me%dynamic%br%dx0 / me%dynamic%br%SpeedS0), 1e30_dp)
end subroutine Junction_resolution_from_and_to_ports


subroutine links_junction_update(me,dt)
    !! Called by Reims for during implicit step. It updates `Matrix` with new values from port
    type(junction_t), intent(inout) :: me
    real(dp), intent(in) :: dt

    real(dp) :: VarNCCenter
    integer :: j

    do j=1,me%NbTotBr ! loop required due to the pointer array: me%br
        VarNCCenter = 0.0_dp
        if(R_Correction) VarNCCenter = me%br(j)%p%Prim(Pri_ro) * me%br(j)%p%Prim(Pri_c)**2
        me%Matrix(:,:,:,j) = me%br(j)%p%Sgn4j * (dt/me%br(j)%p%dxLoc) * me%DerFlx_dBr(:,:,:,j)
        me%Matrix(:,Con_R,:,j) = me%Matrix(:,Con_R,:,j) + &
                me%br(j)%p%Sgn4j * (dt/me%br(j)%p%dxLoc) * VarNCCenter * me%DerVel_dBr(:,:,j)
        me%Matrix(:,:,j,j) = 0.0_dp
    enddo
end subroutine links_junction_update

end module cmp_junction_calc_m

! TODO: Check if sign of the flow velocity in a given branch is in agreement with
! the sign expected before the flux computation. If not, reiterate the junction
! treatment with a changed status of the branch under consideration (incoming or outgoing)
 