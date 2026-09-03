! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_channel_source_terms_m
    use lib_He_thermo_m
    use cmp_channel_init_m
    implicit none
 
contains 
    

subroutine friction_source_term(me)
    type(channel_t), intent(inout) :: me
  
      integer :: i
      real(dp) :: eps, ccv, ccp, mmu, llambda, dPdT_Ro, dEdp_Ro, Re
    real(dp) :: fffL, fffT, fff, dfdRe
      real(dp) :: Vit, Press, e, csound
      real(dp) :: dTdp_Ro,dTdRo_p,Der_dTdP_Ro_dRo_T,Der_dTdP_Ro_dT_Ro
      real(dp) :: ccv_plusRo,ccp_plusRo,mmu_plusRo,llambda_plusRo,dPdT_Ro_plusRo
      real(dp) :: ccv_plusT,ccp_plusT,mmu_plusT,llambda_plusT,dPdT_Ro_plusT
      real(dp) :: DermmudRo_T, DermmudT_Ro, DerccvdRo_T, DerccvdT_Ro
    real(dp) :: Rok,Tempk,Phik
      real(dp), dimension(Nb_VarC) :: dRedUi, dfdUi, Ui, DerpdUi, DerVitdUi
      real(dp), dimension(Nb_VarC) :: DermmudUi, DerccvdUi, DerRok, DerTempk
      real(dp), dimension(Nb_VarC) :: DerPhik, Der_dEdp_Ro_dUi, Der_dTdp_Ro_dUi
    real(dp) :: kk, KKK, dedT, dPdT, HH, dPde
  
      me%srcPP%Su(:)=0.0_dp
      me%srcPP%DerSu(:,:)=0.0_dp
      me%srcPP%Sr(:)=0.0_dp
      me%srcPP%DerSr(:,:)=0.0_dp

      do i=1,me%HeProp%NbCells
    
        eps=1.0e-2_dp
        
        Ui(:)=me%StVar%He_cs(:,i)
        Rok=Ui(Con_Mas)
        Vit=Ui(Con_Qdm)/Rok
        e=Ui(Con_Ene)/Ui(Con_Mas)-0.5_dp*Vit**2
        if(R_Correction) then
            call state_roE_withR(Rok, e, Ui(Con_R), Press, Tempk, csound)
        else
            call state_roE(Rok, e, Press, Tempk, csound)
        endif
     
        call he_prop(Rok,Tempk,ccv,ccp,mmu,llambda,dPdT_Ro)
        if (sim_error > 0) return
        dEdp_Ro=ccv/dPdT_Ro
      
        Re=Rok*abs(Vit)*me%HeProp%Diam/mmu
        fffL=friction_factor_laminar(Re)    ! Considers potential laminar configuration even if Blasius is imposed
        fffT=me%fct%friction(Re)
        fff=max(fffL,fffT)
        if(fffL>fffT) then
          dfdRe=der_friction_factor_laminar(Re)
        else
          dfdRe=me%fct%friction_der(Re)
        endif
      
        DerRok(:)=0.0_dp
        DerRok(Con_Mas)=1.0_dp
                  
        call jacobian_roT(Rok, Tempk, dedT, dPdT, dTdp_Ro, dTdRo_p, &
            d2TdP_dT=Der_dTdP_Ro_dT_Ro, d2TdP_dRo=Der_dTdP_Ro_dRo_T)
        
        call he_prop(Rok+eps,Tempk,ccv_plusRo,ccp_plusRo,mmu_plusRo,llambda_plusRo,dPdT_Ro_plusRo)
        if (sim_error > 0) return
        call he_prop(Rok,Tempk+eps,ccv_plusT,ccp_plusT,mmu_plusT,llambda_plusT,dPdT_Ro_plusT)
        if (sim_error > 0) return
        DerccvdRo_T=(ccv_plusRo-ccv)/eps
        DerccvdT_Ro=(ccv_plusT-ccv)/eps
        DermmudRo_T=(mmu_plusRo-mmu)/eps
        DermmudT_Ro=(mmu_plusT-mmu)/eps
        
        if(R_Correction) then
            DerpdUi(1)=-0.5_dp*(Ui(2)**2)/(Ui(1)**2)
            DerpdUi(2)=Ui(2)/Ui(1)
            DerpdUi(3)=-1.0_dp
            DerpdUi(4)=1.0_dp
        else
            HH=e+0.5_dp*Vit**2+Press/Rok
            call jacobian_roT(Rok, Tempk, dedT, dPdT)
            dPde=dPdT/dedT
            kk=dPde/Rok
            KKK=csound**2+kk*(Vit**2-HH)
    
            DerpdUi(1)=KKK
            DerpdUi(2)=-kk*Vit
            DerpdUi(3)=kk
            DerpdUi(4)=0.0_dp
        endif
  
        DerVitdUi(:)=0.0_dp
        if(Vit>1.0e-6_dp) then
            DerVitdUi(1)=-Ui(2)/(Ui(1)**2)
            DerVitdUi(2)=1.0_dp/Ui(1)          
        endif      
              
        DerTempk(:)=dTdRo_p*DerRok(:)+dTdp_Ro*DerPdUi(:)
     
        DerccvdUi(:)=DerccvdRo_T*DerRok(:)+DerccvdT_Ro*DerTempk(:)
        DermmudUi(:)=DermmudRo_T*DerRok(:)+DermmudT_Ro*DerTempk(:)
        
        Der_dTdp_Ro_dUi(:)=Der_dTdP_Ro_dRo_T*DerRok(:)+Der_dTdP_Ro_dT_Ro*DerTempk(:)
        Der_dEdp_Ro_dUi(:)=DerccvdUi(:)*dTdp_Ro+ccv*Der_dTdp_Ro_dUi(:)
        
        if(Re>1.0e-8_dp) then
            dRedUi(:)=-me%HeProp%Diam*abs(Ui(2))*DermmudUi(:)/(mmu**2)
            dRedUi(2)=dRedUi(2)+me%HeProp%Diam*Ui(2)/(mmu*abs(Ui(2)))
        else
            dRedUi(:)=0.0_dp
        endif          
      
        dfdUi(:)=dfdRe*dRedUi(:)
                
        me%srcPP%Su(i)=-2.0_dp*Rok*fff*Vit*abs(Vit)/me%HeProp%Diam
    
        if(Re>1.0e-8_dp) then
            me%srcPP%DerSu(1,i)=-2.0_dp*Ui(2)*abs(Ui(2))*((dfdUi(1)*Ui(1)-fff)/(Ui(1)**2))/me%HeProp%Diam
            me%srcPP%DerSu(2,i)=-2.0_dp*(dfdUi(2)*Ui(2)*abs(Ui(2))+2.0_dp*fff*(Ui(2)**2)/abs(Ui(2)))/(Ui(1)*me%HeProp%Diam)
            me%srcPP%DerSu(3,i)=-2.0_dp*Ui(2)*abs(Ui(2))*dfdUi(3)/(Ui(1)*me%HeProp%Diam)
            me%srcPP%DerSu(4,i)=-2.0_dp*Ui(2)*abs(Ui(2))*dfdUi(4)/(Ui(1)*me%HeProp%Diam)
        endif 
  
        Phik=1.0_dp/(Rok*dEdp_Ro)
        DerPhik(:)=-(DerRok(:)*dEdp_Ro+Rok*Der_dEdp_Ro_dUi(:))/((Rok*dEdp_Ro)**2)
  
        if(R_Correction) then
          me%srcPP%Sr(i)=-Vit*(1.0_dp+Phik)*me%srcPP%Su(i)
          me%srcPP%DerSr(:,i)=(-DerVitdUi(:)*me%srcPP%Su(i)-Vit*me%srcPP%DerSu(:,i))*(1.0_dp+Phik)-Vit*me%srcPP%Su(i)*DerPhik(:)
        else
          me%srcPP%Sr(i)=0.0_dp
          me%srcPP%DerSr(:,i)=0.0_dp
        endif
        
      enddo
  
end subroutine friction_source_term   


subroutine SourceTerms_NOTfriction_channel(me,i,QorT)
    type(channel_t), intent(inout) :: me
    integer, intent(in) :: i
    real(dp), intent(in) :: QorT

    real(dp) :: eps, ccv, ccp, mmu, llambda, dPdT_Ro, dEdp_Ro, Re
    real(dp) :: Vit, Press, e, csound
    real(dp) :: dTdp_Ro,dTdRo_p,Der_dTdP_Ro_dRo_T,Der_dTdP_Ro_dT_Ro
    real(dp) :: ccv_plusRo,ccp_plusRo,mmu_plusRo,llambda_plusRo,dPdT_Ro_plusRo
    real(dp) :: ccv_plusT,ccp_plusT,mmu_plusT,llambda_plusT,dPdT_Ro_plusT
    real(dp) :: DerccvdRo_T, DerccvdT_Ro, DerccpdRo_T, DerccpdT_Ro, DermmudRo_T, DermmudT_Ro
    real(dp) :: DerllambdadRo_T, DerllambdadT_Ro, Pra, Nu, dNudRe_Pr, dNudPr_Re, dNudPr_Tempk
    real(dp) :: Rok,Tempk,Ht,Phik,WetPerim
    real(dp), dimension(Nb_VarC) :: dRedUi, Ui, DerpdUi, Der_dEdp_Ro_dUi
    real(dp), dimension(Nb_VarC) :: DerccvdUi, DerccpdUi, DermmudUi, DerllambdadUi, Der_dTdp_Ro_dUi
    real(dp), dimension(Nb_VarC) :: dPrdUi, DerNudUi
    real(dp), dimension(Nb_VarC) :: DerRok, DerHt
    real(dp), dimension(Nb_VarC) :: DerPhik, DerTempk
    real(dp) :: kk, KKK, dedT, dPdT, HH, dPde

    eps=1.0e-2_dp
    
    Ui(:)=me%StVar%He_cs(:,i)
    Rok=Ui(Con_Mas)
    Vit=Ui(Con_Qdm)/Rok
    e=Ui(Con_Ene)/Ui(Con_Mas)-0.5_dp*Vit**2
    if(R_Correction) then
        call state_roE_withR(Rok, e, Ui(Con_R), Press, Tempk, csound)
    else
        call state_roE(Rok, e, Press, Tempk, csound)
    endif
 
    call he_prop(Rok,Tempk,ccv,ccp,mmu,llambda,dPdT_Ro)
    if (sim_error > 0) return
    dEdp_Ro=ccv/dPdT_Ro
  
    Re=Rok*abs(Vit)*me%HeProp%Diam/mmu
  
    DerRok(:)=0.0_dp
    DerRok(Con_Mas)=1.0_dp
              
    call jacobian_roT(Rok, Tempk, dedT, dPdT, dTdp_Ro, dTdRo_p, &
        d2TdP_dT=Der_dTdP_Ro_dT_Ro, d2TdP_dRo=Der_dTdP_Ro_dRo_T)
    
    call he_prop(Rok+eps,Tempk,ccv_plusRo,ccp_plusRo,mmu_plusRo,llambda_plusRo,dPdT_Ro_plusRo)
    if (sim_error > 0) return
    call he_prop(Rok,Tempk+eps,ccv_plusT,ccp_plusT,mmu_plusT,llambda_plusT,dPdT_Ro_plusT)
    if (sim_error > 0) return
    DerccvdRo_T=(ccv_plusRo-ccv)/eps
    DerccvdT_Ro=(ccv_plusT-ccv)/eps
    DerccpdRo_T=(ccp_plusRo-ccp)/eps
    DerccpdT_Ro=(ccp_plusT-ccp)/eps
    DermmudRo_T=(mmu_plusRo-mmu)/eps
    DermmudT_Ro=(mmu_plusT-mmu)/eps
    DerllambdadRo_T=(llambda_plusRo-llambda)/eps
    DerllambdadT_Ro=(llambda_plusT-llambda)/eps
    
    if(R_Correction) then
        DerpdUi(1)=-0.5_dp*(Ui(2)**2)/(Ui(1)**2)
        DerpdUi(2)=Ui(2)/Ui(1)
        DerpdUi(3)=-1.0_dp
        DerpdUi(4)=1.0_dp
    else
        HH=e+0.5_dp*Vit**2+Press/Rok
        call jacobian_roT(Rok, Tempk, dedT, dPdT)
        dPde=dPdT/dedT
        kk=dPde/Rok
        KKK=csound**2+kk*(Vit**2-HH)

        DerpdUi(1)=KKK
        DerpdUi(2)=-kk*Vit
        DerpdUi(3)=kk
        DerpdUi(4)=0.0_dp
    endif
           
    DerTempk(:)=dTdRo_p*DerRok(:)+dTdp_Ro*DerPdUi(:)
   
    DerccvdUi(:)=DerccvdRo_T*DerRok(:)+DerccvdT_Ro*DerTempk(:)
    DerccpdUi(:)=DerccpdRo_T*DerRok(:)+DerccpdT_Ro*DerTempk(:)
    DermmudUi(:)=DermmudRo_T*DerRok(:)+DermmudT_Ro*DerTempk(:)
    DerllambdadUi(:)=DerllambdadRo_T*DerRok(:)+DerllambdadT_Ro*DerTempk(:)
    
    Der_dTdp_Ro_dUi(:)=Der_dTdP_Ro_dRo_T*DerRok(:)+Der_dTdP_Ro_dT_Ro*DerTempk(:)
    Der_dEdp_Ro_dUi(:)=DerccvdUi(:)*dTdp_Ro+ccv*Der_dTdp_Ro_dUi(:)
    
    if(Re>1.0e-8_dp) then
        dRedUi(:)=-me%HeProp%Diam*abs(Ui(2))*DermmudUi(:)/(mmu**2)
        dRedUi(2)=dRedUi(2)+me%HeProp%Diam*Ui(2)/(mmu*abs(Ui(2)))
    else
        dRedUi(:)=0.0_dp
    endif          
  
    Pra=mmu*ccp/llambda
    dPrdUi(:)=(llambda*(DermmudUi(:)*ccp+mmu*DerccpdUi(:))-mmu*ccp*DerllambdadUi(:))/(llambda**2)
        
    if(me%HeProp%ExtHeating=='temp') then ! used in general for the heat exchanger treatment
        if(associated(me%nuss)) then
           Nu=me%nuss%nusselt(Re,Pra,QorT,Tempk)
           dNudRe_Pr=me%nuss%nusselt_der_Re(Re,Pra,QorT,Tempk)
           dNudPr_Re=me%nuss%nusselt_der_Pra(Re,Pra,QorT,Tempk)
           dNudPr_Tempk=me%nuss%nusselt_der_T(Re,Pra,QorT,Tempk)
           DerNudUi(:)=dNudRe_Pr*dRedUi(:)+dNudPr_Re*dPrdUi(:)+dNudPr_Tempk*DerTempk(:)
           Ht=Nu*llambda/me%HeProp%Diam
           DerHt(:)=(1.0_dp/me%HeProp%Diam)*(DerNudUi(:)*llambda+Nu*DerllambdadUi(:))   
        else        
           Ht=me%HeProp%Ht
           DerHt(:)=0.0_dp            
        endif      
        WetPerim=4.0_dp*me%HeProp%Area/me%HeProp%Diam          
        me%srcPP%Se=WetPerim*Ht*(QorT-Tempk)/me%HeProp%Area
        me%srcPP%DerSe(:)=(WetPerim/me%HeProp%Area)*(DerHt(:)*(QorT-Tempk)-Ht*DerTempk(:))
    else if(me%HeProp%ExtHeating=='flux') then
        me%srcPP%Se=QorT/me%HeProp%Area
        me%srcPP%DerSe(:)=0.0_dp   
    else
        me%srcPP%Se=0.0_dp
        me%srcPP%DerSe(:)=0.0_dp
    endif
    
    Phik=1.0_dp/(Rok*dEdp_Ro)
    DerPhik(:)=-(DerRok(:)*dEdp_Ro+Rok*Der_dEdp_Ro_dUi(:))/((Rok*dEdp_Ro)**2)
    
    if(R_Correction) then
        me%srcPP%Sr(i)=me%srcPP%Sr(i)+me%srcPP%Se*(1.0_dp+Phik)  
        me%srcPP%DerSr(:,i)=me%srcPP%DerSr(:,i)+me%srcPP%DerSe(:)*(1.0_dp+Phik)+me%srcPP%Se*DerPhik(:)
    endif

    if(me%HeProp%ExtHeating=='link') then
      ! Communication to FS port associated to channel node i
      me%HeProp%thermP(i)%p%Pra=Pra
      me%HeProp%thermP(i)%p%Re=Re
      me%HeProp%thermP(i)%p%dPrdUi(:)=dPrdUi(:)
      me%HeProp%thermP(i)%p%dRedUi(:)=dRedUi(:)
      me%HeProp%thermP(i)%p%Phi=Phik
      me%HeProp%thermP(i)%p%DerPhi(:)=DerPhik(:)
      me%HeProp%thermP(i)%p%llambda=llambda
      me%HeProp%thermP(i)%p%DerllambdadUi(:)=DerllambdadUi(:)
      me%HeProp%thermP(i)%p%TempP=Tempk
      me%HeProp%thermP(i)%p%DerTempPdUi(:)=DerTempk(:)
    endif

    if(me%HeProp%FFsrcLink) then
      ! Communication to FF_src port associated to channel node i
      me%HeProp%thermF(i)%p%Pra=Pra
      me%HeProp%thermF(i)%p%Re=Re
      me%HeProp%thermF(i)%p%dPrdUi(:)=dPrdUi(:)
      me%HeProp%thermF(i)%p%dRedUi(:)=dRedUi(:)
      me%HeProp%thermF(i)%p%Phi=Phik
      me%HeProp%thermF(i)%p%DerPhi(:)=DerPhik(:)
      me%HeProp%thermF(i)%p%llambda=llambda
      me%HeProp%thermF(i)%p%DerllambdadUi(:)=DerllambdadUi(:)
      me%HeProp%thermF(i)%p%TempP=Tempk
      me%HeProp%thermF(i)%p%DerTempPdUi(:)=DerTempk(:)
    endif

    end subroutine SourceTerms_NOTfriction_channel    
     
end module cmp_channel_source_terms_m