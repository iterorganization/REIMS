! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_channel_flux_m
    use lib_He_thermo_m
    use cmp_channel_init_m
    implicit none

    type cell_state_t
      real(dp) :: Prim(Nb_VarP)
      real(dp) :: Cons(Nb_VarC)
      real(dp) :: Jcb(Nb_VarC,Nb_VarC)
      real(dp) :: SuFrct,SrFrct,dx
      real(dp) :: DerSuFrct(Nb_VarC),DerSrFrct(Nb_VarC)
    end type cell_state_t
    
contains

subroutine states_L_and_R_channel(me,i,Cell_L,Cell_R)  

    type(channel_t), intent(in) :: me
    type(cell_state_t), intent(out) :: Cell_L,Cell_R
    integer, intent(in) :: i

    ! Internal boundary cells
    Cell_L%Prim(:)=me%StVar%He_pr(:,i)
    Cell_R%Prim(:)=me%StVar%He_pr(:,i+1)

    Cell_L%Cons(:)=me%StVar%He_cs(:,i)
    Cell_R%Cons(:)=me%StVar%He_cs(:,i+1)

    Cell_L%Jcb(:,:)=me%big%Jacob(:,:,i)
    Cell_R%Jcb(:,:)=me%big%Jacob(:,:,i+1)

    Cell_L%SuFrct=me%srcPP%Su(i)
    Cell_R%SuFrct=me%srcPP%Su(i+1)
    Cell_L%DerSuFrct(:)=me%srcPP%DerSu(:,i)
    Cell_R%DerSuFrct(:)=me%srcPP%DerSu(:,i+1)

    Cell_L%SrFrct=me%srcPP%Sr(i)
    Cell_R%SrFrct=me%srcPP%Sr(i+1)
    Cell_L%DerSrFrct(:)=me%srcPP%DerSr(:,i)
    Cell_R%DerSrFrct(:)=me%srcPP%DerSr(:,i+1)

    Cell_L%dx=me%HeProp%dxLoc(i)
    Cell_R%dx=me%HeProp%dxLoc(i+1)
end subroutine states_L_and_R_channel
    
  
subroutine wave_speed_calculation_channel(Cell_L,Cell_R,dUldUl,dUrdUr,sm,sl,sr,dSmdUl,dSmdUr,sim)
    type(cell_state_t), intent(in) :: Cell_L,Cell_R
    real(dp), dimension(Nb_VarC,Nb_VarC) :: dUldUl,dUrdUr
    real(dp), intent(out) :: sm,sl,sr
    real(dp), dimension(:), intent(out) :: dSmdUl,dSmdUr
    type(simulation_t), intent(in) :: sim

    real(dp), dimension(Nb_VarP) :: StatL,StatR
    real(dp), dimension(Nb_VarC) :: Ul,Ur,SrcL,SrcR
    real(dp), dimension(Nb_VarC,Nb_VarC) :: JacobL,JacobR,DerSrcL,DerSrcR
    real(dp) :: SS,num,den,dxL,dxR
    real(dp), dimension(2) :: Fl,Fr
    
    SrcL=0.0_dp;SrcR=0.0_dp;DerSrcL=0.0_dp;DerSrcR=0.0_dp

    StatL=Cell_L%Prim
    StatR=Cell_R%Prim
    Ul=Cell_L%Cons
    Ur=Cell_R%Cons
    JacobL=Cell_L%Jcb
    JacobR=Cell_R%Jcb
    SrcL(2)=Cell_L%SuFrct
    SrcR(2)=Cell_R%SuFrct
    DerSrcL(:,2)=Cell_L%DerSuFrct(:)
    DerSrcR(:,2)=Cell_R%DerSuFrct(:)
    dxL=Cell_L%dx
    dxR=Cell_R%dx

    SS=max(abs(StatL(Pri_u))+StatL(Pri_c),abs(StatR(Pri_u))+StatR(Pri_c))*sim%scheme_diff_fact
    sr=SS
    sl=-SS
    
    Fl(1)=Ul(2)
    Fl(2)=StatL(Pri_ro)*StatL(Pri_u)**2+StatL(Pri_P)
    if(R_Correction) Fl(2)=(3.0_dp/2.0_dp)*(Ul(2)**2)/Ul(1)+Ul(4)-Ul(3)

    Fr(1)=Ur(2)
    Fr(2)=StatR(Pri_ro)*StatR(Pri_u)**2+StatR(Pri_P)
    if(R_Correction) Fr(2)=(3.0_dp/2.0_dp)*(Ur(2)**2)/Ur(1)+Ur(4)-Ur(3)

    num=Fr(2)-Fl(2)-sr*Ur(2)+sl*Ul(2)-(dxL*SrcL(2)+dxR*SrcR(2))/2.0_dp
    den=Fr(1)-Fl(1)-sr*Ur(1)+sl*Ul(1)

    sm=num/den

    dSmdUl(:)=((-JacobL(:,2)+sl*dUldUl(:,2)-0.5_dp*dxL*DerSrcL(:,2))*den-num*(-JacobL(:,1)+sl*dUldUl(:,1)))/(den**2)
    dSmdUr(:)=((JacobR(:,2)-sr*dUrdUr(:,2)-0.5_dp*dxR*DerSrcR(:,2))*den-num*(JacobR(:,1)-sr*dUrdUr(:,1)))/(den**2)
end subroutine wave_speed_calculation_channel    
  
  
subroutine flux_internal_calculation_channel(i,HepropPP,Cell_L,Cell_R,dtMaxLoc,flux,sim)    
    type(cell_state_t), intent(in) :: Cell_L,Cell_R
    real(dp), intent(out) :: dtMaxLoc
    type(flux_He_channel_t), intent(inout) :: flux
    integer, intent(in) :: i
    type(He_prop_t), intent(in) :: HePropPP
    type(simulation_t), intent(in) :: sim

    real(dp), dimension(Nb_VarC,Nb_VarC) :: dUldUl,dUrdUr
    real(dp), dimension(Nb_VarP) :: StatLstar,StatRstar
    real(dp), dimension(Nb_VarC) :: Fl,Fr,dSmdUl,dSmdUr,dpstardUl,dpstardUr,ConsLstar,ConsRstar
    real(dp), dimension(Nb_VarC) :: Fstar,LastTermL,LastTermR
    real(dp), dimension(Nb_VarC,Nb_VarC) :: derUstarL_dUL,derUstarL_dUR,derUstarR_dUL,derUstarR_dUR
    real(dp), dimension(Nb_VarC,Nb_VarC) :: dFstardUl,dFstardUr,LastMatrix,AddThermL,AddThermR
    real(dp), dimension(Nb_VarP) :: StatL,StatR
    real(dp), dimension(Nb_VarC) :: Ul,Ur,SrcL,SrcR
    real(dp), dimension(Nb_VarC,Nb_VarC) :: JacobL,JacobR,DerSrcL,DerSrcR
    real(dp) :: Etstar,pstar,dxL,dxR
    real(dp) :: sm,sl,sr
    integer :: n,imPP

    StatLstar=0.0_dp;StatRstar=0.0_dp;Fl=0.0_dp;Fr=0.0_dp

    imPP=HepropPP%NbCells

    dUldUl(:,:)=id_4x4(:,:)
    dUrdUr(:,:)=id_4x4(:,:)

    SrcL=0.0_dp;SrcR=0.0_dp;DerSrcL=0.0_dp;DerSrcR=0.0_dp

    call wave_speed_calculation_channel(Cell_L,Cell_R,dUldUl,dUrdUr,sm,sl,sr,dSmdUl(:),dSmdUr(:),sim)

    StatL=Cell_L%Prim
    StatR=Cell_R%Prim
    Ul=Cell_L%Cons
    Ur=Cell_R%Cons
    JacobL=Cell_L%Jcb
    JacobR=Cell_R%Jcb
    SrcL(2)=Cell_L%SuFrct
    SrcR(2)=Cell_R%SuFrct
    DerSrcL(:,2)=Cell_L%DerSuFrct(:)
    DerSrcR(:,2)=Cell_R%DerSuFrct(:)

    if(R_Correction) then
      SrcL(4)=Cell_L%SrFrct
      SrcR(4)=Cell_R%SrFrct
      DerSrcL(:,4)=Cell_L%DerSrFrct(:)
      DerSrcR(:,4)=Cell_R%DerSrFrct(:)
    endif
    dxL=Cell_L%dx
    dxR=Cell_R%dx
    
    if(max(abs(sl),abs(sr))<=1.0e-8_dp) then
        call set_error('flux denominator zero - CFL criterion for non-uniform mesh channels')
        return
    endif
    if(i==0) then
        dtMaxLoc=Cell_R%dx/max(abs(sl),abs(sr))
    else if(i==imPP) then
        dtMaxLoc=Cell_L%dx/max(abs(sl),abs(sr))
    else
        dtMaxLoc=min(Cell_L%dx/max(abs(sl),abs(sr)),Cell_R%dx/max(abs(sl),abs(sr)))   ! Facteur 1/2 ?
    endif
    
    Fl(Con_Mas)=StatL(Pri_ro)*StatL(Pri_u)
    Fl(Con_Qdm)=StatL(Pri_ro)*StatL(Pri_u)*StatL(Pri_u)+StatL(Pri_p)
    Fl(Con_Ene)=(StatL(Pri_ro)*(StatL(Pri_e)+0.5_dp*StatL(Pri_u)*StatL(Pri_u))+StatL(Pri_p))*StatL(Pri_u)
    Fl(Con_R)=0.0_dp
    if(R_Correction) Fl(Con_R)=StatL(Pri_R)*StatL(Pri_u)

    Fr(Con_Mas)=StatR(Pri_ro)*StatR(Pri_u)
    Fr(Con_Qdm)=StatR(Pri_ro)*StatR(Pri_u)*StatR(Pri_u)+StatR(Pri_p)
    Fr(Con_Ene)=(StatR(Pri_ro)*(StatR(Pri_e)+0.5_dp*StatR(Pri_u)*StatR(Pri_u))+StatR(Pri_p))*StatR(Pri_u)
    Fr(Con_R)=0.0_dp
    if(R_Correction) Fr(Con_R)=StatR(Pri_R)*StatR(Pri_u)

    ! Wave sampling
    pstar=0.5_dp*(Fr(2)-sr*Ur(2)-sm*(Fr(1)-sr*Ur(1))+Fl(2)-sl*Ul(2)-sm*(Fl(1)-sl*Ul(1)))

    if(sm>=0.0_dp) then
        ! Left star state
        StatLstar(Pri_ro)=StatL(Pri_ro)*(StatL(Pri_u)-sl)/(sm-sl)
        Etstar=(StatL(Pri_e)+0.5_dp*StatL(Pri_u)*StatL(Pri_u))+(StatL(Pri_p)*&
               StatL(Pri_u)-pstar*sm)/(StatL(Pri_ro)*(StatL(Pri_u)-sl))            
        StatLstar(Pri_R)=(StatL(Pri_R)*(StatL(Pri_u)-sl)+dxL*SrcL(4)/2.0_dp)/(sm-sl)
  
        flux%Cons(Con_Mas,i)=Fl(Con_Mas)+sl*(StatLstar(Pri_ro)-StatL(Pri_ro))
        flux%Cons(Con_Qdm,i)=Fl(Con_Qdm)+sl*(StatLstar(Pri_ro)*sm-StatL(Pri_ro)*StatL(Pri_u))+dxL*SrcL(2)/2.0_dp
        flux%Cons(Con_Ene,i)=Fl(Con_Ene)+sl*(StatLstar(Pri_ro)*Etstar-StatL(Pri_ro)*(StatL(Pri_e)+&
                                  0.5_dp*StatL(Pri_u)*StatL(Pri_u)))        
        flux%Cons(Con_R,i)=0.0_dp                                     
        if(R_Correction) flux%Cons(Con_R,i)=Fl(Con_R)+sl*(StatLstar(Pri_R)-StatL(Pri_R))+dxL*SrcL(4)/2.0_dp
        flux%VitTNC(i)=sm
    else
        ! Right star state
        StatRstar(Pri_ro)=StatR(Pri_ro)*(StatR(Pri_u)-sr)/(sm-sr)
        Etstar=(StatR(Pri_e)+0.5_dp*StatR(Pri_u)*StatR(Pri_u))+(StatR(Pri_p)*&
               StatR(Pri_u)-pstar*sm)/(StatR(Pri_ro)*(StatR(Pri_u)-sr))
        StatRstar(Pri_R)=(StatR(Pri_R)*(StatR(Pri_u)-sr)-dxR*SrcR(4)/2.0_dp)/(sm-sr)
  
        flux%Cons(Con_Mas,i)=Fr(Con_Mas)+sr*(StatRstar(Pri_ro)-StatR(Pri_ro))
        flux%Cons(Con_Qdm,i)=Fr(Con_Qdm)+sr*(StatRstar(Pri_ro)*sm-StatR(Pri_ro)*StatR(Pri_u))-dxR*SrcR(2)/2.0_dp
        flux%Cons(Con_Ene,i)=Fr(Con_Ene)+sr*(StatRstar(Pri_ro)*Etstar-StatR(Pri_ro)*(StatR(Pri_e)+&
                                  0.5_dp*StatR(Pri_u)*StatR(Pri_u)))
        flux%Cons(Con_R,i)=0.0_dp
        if(R_Correction) flux%Cons(Con_R,i)=Fr(Con_R)+sr*(StatRstar(Pri_R)-StatR(Pri_R))-dxR*SrcR(4)/2.0_dp
        flux%VitTNC(i)=sm
    endif
  
    ! Flux derivative calculation
    Fstar(1)=0.0_dp
    Fstar(2)=pstar
    Fstar(3)=pstar*sm
    Fstar(4)=0.0_dp

    LastTermL(:)=Fl(:)-sl*Ul(:)-Fstar(:)+dxL*SrcL(:)/2.0_dp
    LastTermR(:)=Fr(:)-sr*Ur(:)-Fstar(:)-dxR*SrcR(:)/2.0_dp
    if(.not. R_Correction) LastTermL(4)=0.0_dp; LastTermR(4)=0.0_dp

    ConsLstar(:)=LastTermL(:)/(sm-sl)
    ConsRstar(:)=LastTermR(:)/(sm-sr)
    do n=1,Nb_VarC
        AddThermL(:,n)=(ConsRstar(n)-ConsLstar(n))*dSmdUl(:)
        AddThermR(:,n)=(ConsRstar(n)-ConsLstar(n))*dSmdUr(:)
    enddo
  
    ! DerFlux_dUL calculation
    dpstardUl(:)=0.5_dp*(JacobL(:,2)-sl*dUldUl(:,2)-sm*(JacobL(:,1)-sl*dUldUl(:,1))-&
                (Fl(1)-sl*Ul(1))*dSmdUl(:)-(Fr(1)-sr*Ur(1))*dSmdUl(:)) 

    dFstardUl(:,:)=0.0_dp
    dFstardUl(:,2)=dpstardUl(:)
    dFstardUl(:,3)=sm*dpstardUl(:)+pstar*dSmdUl(:)

    do n=1,4
        LastMatrix(:,n)=LastTermL(n)*dSmdUl(:)
    enddo

    derUstarL_dUL(:,:)=((JacobL(:,:)-sl*id_4x4(:,:)-dFstardUl(:,:)+dxL*DerSrcL(:,:)/2.0_dp)*(sm-sl)-LastMatrix(:,:))/((sm-sl)**2)

    do n=1,4
        LastMatrix(:,n)=LastTermR(n)*dSmdUl(:)
    enddo    
    derUstarR_dUL(:,:)=((-dFstardUl(:,:))*(sm-sr)-LastMatrix(:,:))/((sm-sr)**2)

    if(i/=0) THEN
        flux%Cons_DerdUL(:,:,i)=0.5_dp*JacobL(:,:)-&
                        sign(1.0_dp,sl)*(sl/2.0_dp)*(derUstarL_dUL(:,:)-id_4x4(:,:))-&
                        sign(1.0_dp,sm)*(sm/2.0_dp)*(derUstarR_dUL(:,:)-derUstarL_dUL(:,:))+&
                        sign(1.0_dp,sr)*(sr/2.0_dp)*derUstarR_dUL(:,:)-&
                        0.5_dp*sign(1.0_dp,sm)*AddThermL(:,:)  
        flux%VitTNC_DerdUL(:,i)=dSmdUl(:)
    endif
                          
    ! DerFlux_dUR calculation   
    dpstardUr(:)=0.5_dp*(JacobR(:,2)-sr*dUrdUr(:,2)-sm*(JacobR(:,1)-sr*dUrdUr(:,1))-&
                (Fr(1)-sr*Ur(1))*dSmdUr(:)-(Fl(1)-sl*Ul(1))*dSmdUr(:))

    dFstardUr(:,:)=0.0_dp
    dFstardUr(:,2)=dpstardUr(:)
    dFstardUr(:,3)=sm*dpstardUr(:)+pstar*dSmdUr(:)

    do n=1,4
    LastMatrix(:,n)=LastTermL(n)*dSmdUr(:)
    enddo    
    derUstarL_dUR(:,:)=((-dFstardUr(:,:))*(sm-sl)-LastMatrix(:,:))/((sm-sl)**2)

    do n=1,4
    LastMatrix(:,n)=LastTermR(n)*dSmdUr(:)
    enddo    

    derUstarR_dUR(:,:)=((JacobR(:,:)-sr*id_4x4(:,:)-dFstardUr(:,:)-dxR*DerSrcR(:,:)/2.0_dp)*(sm-sr)-LastMatrix(:,:))/((sm-sr)**2)    

    if(i/=imPP) then
        flux%Cons_DerdUR(:,:,i)=0.5_dp*JacobR(:,:)-&
                        sign(1.0_dp,sl)*(sl/2.0_dp)*derUstarL_dUR(:,:)-&
                        sign(1.0_dp,sm)*(sm/2.0_dp)*(derUstarR_dUR(:,:)-derUstarL_dUR(:,:))-&
                        sign(1.0_dp,sr)*(sr/2.0_dp)*(id_4x4(:,:)-derUstarR_dUR(:,:))-&
                        0.5_dp*sign(1.0_dp,sm)*AddThermR(:,:)
        
        flux%VitTNC_DerdUR(:,i)=dSmdUr(:)
    endif 
end subroutine flux_internal_calculation_channel    

  
subroutine jacobian_channel(i,Cons,Jacob)
    integer, intent(in) :: i
    real(dp), dimension(:,:), intent(in) :: Cons
    real(dp), dimension(:,:,:), intent(inout) :: Jacob

    real(dp), dimension(Nb_VarC) :: U

    U(:)=Cons(:,i)

    Jacob(1,1,i)=0.0_dp
    Jacob(2,1,i)=1.0_dp
    Jacob(3,1,i)=0.0_dp
    Jacob(4,1,i)=0.0_dp

    Jacob(1,2,i)=-(3.0_dp/2.0_dp)*((U(2)**2)/(U(1)**2))
    Jacob(2,2,i)=3.0_dp*U(2)/U(1)
    Jacob(3,2,i)=-1.0_dp
    Jacob(4,2,i)=1.0_dp

    Jacob(1,3,i)=-U(4)*U(2)/(U(1)**2)-((U(2)**3)/(U(1)**3))
    Jacob(2,3,i)=U(4)/U(1)+(3.0_dp/2.0_dp)*((U(2)**2)/(U(1)**2))
    Jacob(3,3,i)=0.0_dp
    Jacob(4,3,i)=U(2)/U(1)
    
    Jacob(1,4,i)=-U(4)*U(2)/(U(1)**2)
    Jacob(2,4,i)=U(4)/U(1)
    Jacob(3,4,i)=0.0_dp
    Jacob(4,4,i)=U(2)/U(1)
      
end subroutine jacobian_channel        

  
subroutine jacobian_channel_Req_rem(i,Prim,Jacob)
    integer, intent(in) :: i
    real(dp), dimension(:,:), intent(in) :: Prim
    real(dp), dimension(:,:,:), intent(inout) :: Jacob

    real(dp), dimension(Nb_VarP) :: PPP
    real(dp) :: HH,kk,KKK,dPde,dPdT,dedT

    PPP(:)=Prim(:,i)
    HH=PPP(Pri_e)+0.5_dp*PPP(Pri_u)**2+PPP(Pri_p)/PPP(Pri_ro)
    call jacobian_roT(PPP(Pri_ro), PPP(Pri_T), dedT, dPdT)
    dPde=dPdT/dedT
    kk=dPde/PPP(Pri_ro) ! checked at room temperatures --> k+1=1.667==gamma
    KKK=PPP(Pri_c)**2+kk*(PPP(Pri_u)**2-HH)

    Jacob(1,1,i)=0.0_dp
    Jacob(2,1,i)=1.0_dp
    Jacob(3,1,i)=0.0_dp
    Jacob(4,1,i)=0.0_dp

    Jacob(1,2,i)=KKK-PPP(Pri_u)**2
    Jacob(2,2,i)=PPP(Pri_u)*(2.0_dp-kk)
    Jacob(3,2,i)=kk
    Jacob(4,2,i)=0.0_dp

    Jacob(1,3,i)=(KKK-HH)*PPP(Pri_u)
    Jacob(2,3,i)=HH-kk*PPP(Pri_u)**2
    Jacob(3,3,i)=PPP(Pri_u)*(1.0_dp+kk)
    Jacob(4,3,i)=0.0_dp
    
    Jacob(1,4,i)=0.0_dp
    Jacob(2,4,i)=0.0_dp
    Jacob(3,4,i)=0.0_dp
    Jacob(4,4,i)=0.0_dp
end subroutine jacobian_channel_Req_rem      
  
       
subroutine He_Riemann_solver_channel(me,sim)
    type(channel_t), intent(inout) :: me
    type(simulation_t), intent(in) :: sim
  
    type(cell_state_t) :: Cell_L,Cell_R
    real(dp) :: dtMaxLoc
    integer :: i,imPP
      
    me%big%Jacob=0.0_dp
    me%wave_time = 1.0e30_dp

    imPP=me%HeProp%NbCells
      
    do i=1,imPP     
        call jacobian_channel_Req_rem(i,me%StVar%He_pr,me%big%Jacob)
        if(R_Correction) call jacobian_channel(i,me%StVar%He_cs,me%big%Jacob)
    enddo

    do i=1,imPP-1
        call states_L_and_R_channel(me,i,Cell_L,Cell_R)
        call flux_internal_calculation_channel(i,me%HeProp,Cell_L,Cell_R,dtMaxLoc,me%flxHe,sim)
        if (sim_error > 0) return

        me%wave_time = min(abs(dtMaxLoc),me%wave_time)    
    enddo
      
end subroutine He_Riemann_solver_channel      

end module cmp_channel_flux_m