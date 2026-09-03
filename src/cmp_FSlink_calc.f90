! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_FSlink_calc_m
    use cmp_FSlink_init_m
    implicit none
  
contains


subroutine FSlink_src_reinitialization(me)
  type(FSlink_t), intent(inout) :: me

  integer :: i

  do i = 1, me%sizeNodes
    me%lk(1,i)%p%rhs_SrcP_Se   = 0.0_dp
    me%lk(1,i)%p%rhs_SrcP_Sr   = 0.0_dp
    me%lk(1,i)%p%DerSrcP_Se(:) = 0.0_dp
    me%lk(1,i)%p%DerSrcP_Sr(:) = 0.0_dp

    me%lk(2,i)%p%rhs_SrcS  = 0.0_dp
    me%lk(2,i)%p%DerSrcSdT = 0.0_dp
  enddo

end subroutine FSlink_src_reinitialization



subroutine FSlink_resolution_from_and_to_ports(me)
    type(FSlink_t), intent(inout) :: me

    real(dp) :: Pra,Re,Phi,llambda,TempP,AreaP,DiamP,TempS,AreaMS,rhoMS,cpMS,dcpMSdT
    real(dp), dimension(Nb_VarC) :: dPrdUi,dRedUi,DerPhi,DerllambdadUi,DerTempPdUi,DerFlx2DdUi
    real(dp) :: Nu,dNudRe,dNudPr,dNudTemp,dNudTsolid,Ht,WetPerwire,SrcP_Se,SrcP_Sr,SrcS,DerHtdTsolid,DeltaT,DerSrcSdT
    real(dp) :: DerSrcP_SedT,DerSrcP_SrdT
    real(dp), dimension(Nb_VarC) :: DerNudUi,DerSrcP_Se,DerSrcP_Sr,DerSrcSdUi,DerHt,DerTempInt_dUi
    real(dp) :: SurfCont,Dist4Grad,lambdaS,dlambdaSdT,dxLocP,VolMS,TempInt,DerTempInt_dTm,LengthCont
    real(dp) :: Flx2D,derFlx2D_derCon
    integer :: i, ii

    do i = 1, me%sizeNodes

      ! From pipe node
      Pra=me%lk(1,i)%p%Pra
      Re=me%lk(1,i)%p%Re
      dPrdUi(:)=me%lk(1,i)%p%dPrdUi(:)
      dRedUi(:)=me%lk(1,i)%p%dRedUi(:)
      Phi=me%lk(1,i)%p%Phi
      DerPhi(:)=me%lk(1,i)%p%DerPhi(:)
      llambda=me%lk(1,i)%p%llambda
      DerllambdadUi(:)=me%lk(1,i)%p%DerllambdadUi(:)
      TempP=me%lk(1,i)%p%TempP
      DerTempPdUi(:)=me%lk(1,i)%p%DerTempPdUi(:)
      AreaP=me%lk(1,i)%p%AreaP
      DiamP=me%lk(1,i)%p%DiamP
      dxLocP=me%lk(1,i)%p%dxLoc

      if(me%lk(2,i)%p%typeS==trim("strand")) then

        WetPerwire = me%WetPer

        ! From strand node
        TempS   = me%lk(2,i)%p%TempS
        AreaMS  = me%lk(2,i)%p%AreaS
        rhoMS   = me%lk(2,i)%p%rhoMS
        cpMS    = me%lk(2,i)%p%cpMS
        dcpMSdT = me%lk(2,i)%p%dcpMSdT
  
        Nu=me%nuss%nusselt(Re,Pra,TempS,TempP)
        dNudRe=me%nuss%nusselt_der_Re(Re,Pra,TempS,TempP)
        dNudPr=me%nuss%nusselt_der_Pra(Re,Pra,TempS,TempP)
        dNudTemp=me%nuss%nusselt_der_T(Re,Pra,TempS,TempP)
        dNudTsolid=me%nuss%nusselt_der_Tother(Re,Pra,TempS,TempP)
        DerNudUi(:)=dNudRe*dRedUi(:)+dNudPr*dPrdUi(:)+dNudTemp*DerTempPdUi(:)

        Ht=Nu*llambda/DiamP
        DerHt(:)=(1.0_dp/DiamP)*(DerNudUi(:)*llambda+Nu*DerllambdadUi(:))

        SrcP_Se=WetPerwire*Ht*(TempS-TempP)/AreaP
        DerSrcP_Se(:)=(WetPerwire/AreaP)*(DerHt(:)*(TempS-TempP)-Ht*DerTempPdUi(:))
      
        DerHtdTsolid=dNudTsolid*llambda/DiamP

        if(R_Correction) then
          SrcP_Sr=SrcP_Se*(1.0_dp+Phi)  
          DerSrcP_Sr(:)=DerSrcP_Se(:)*(1.0_dp+Phi)+SrcP_Se*DerPhi
        else
          SrcP_Sr=0.0_dp
          DerSrcP_Sr(:)=0.0_dp
        endif

        SrcS=-AreaP*SrcP_Se/(AreaMS*rhoMS*cpMS)

        DeltaT=TempS-TempP

        DerSrcSdUi(:)=-AreaP*DerSrcP_Se(:)/(AreaMS*rhoMS*cpMS)
        DerSrcSdT=-(WetPerwire/(AreaMS*rhoMS))*((DerHtdTsolid*DeltaT+Ht)*cpMS-&
                     Ht*DeltaT*dcpMSdT)/(cpMS**2)

        DerSrcP_SedT=(WetPerwire/AreaP)*(DerHtdTsolid*DeltaT+Ht)

        if(R_Correction) then
          DerSrcP_SrdT=(1.0_dp+Phi)*DerSrcP_SedT
        else
          DerSrcP_SrdT=0.0_dp
        endif

        ! To strand node --> considering FS accumulation
        !$omp critical(FS_state_accum)
        me%lk(2,i)%p%rhs_SrcS  = me%lk(2,i)%p%rhs_SrcS+SrcS
        me%lk(2,i)%p%DerSrcSdT = me%lk(2,i)%p%DerSrcSdT+DerSrcSdT     
        !$omp end critical(FS_state_accum)  
        
        ! Link part --> FS accumulation not required because different locations for links
        me%DerSrcP_SedT(i)=DerSrcP_SedT
        me%DerSrcP_SrdT(i)=DerSrcP_SrdT
        me%DerSrcSdUi(:,i)=DerSrcSdUi(:)        

      else if(me%lk(2,i)%p%typeS==trim("solid")) then

        SurfCont  = me%SurfCont(i) 

        ! From solid node
        rhoMS     = me%lk(2,i)%p%rhoMS
        VolMS     = me%lk(2,i)%p%VolMS
        Dist4Grad = me%lk(2,i)%p%Dist4Grad

        TempS      = me%lk(2,i)%p%TempS
        lambdaS    = me%lk(2,i)%p%lambdS
        dlambdaSdT = me%lk(2,i)%p%DerlambdS
        cpMs       = me%lk(2,i)%p%cpS
        dcpMsdT    = me%lk(2,i)%p%dcpSdT
       
        Nu=me%nuss%nusselt(Re,Pra,TempS,TempP)
        dNudRe=me%nuss%nusselt_der_Re(Re,Pra,TempS,TempP)
        dNudPr=me%nuss%nusselt_der_Pra(Re,Pra,TempS,TempP)
        dNudTemp=me%nuss%nusselt_der_T(Re,Pra,TempS,TempP)
        dNudTsolid=me%nuss%nusselt_der_Tother(Re,Pra,TempS,TempP)
        DerNudUi(:)=dNudRe*dRedUi(:)+dNudPr*dPrdUi(:)+dNudTemp*DerTempPdUi(:)

        Ht=Nu*llambda/DiamP
        DerHt(:)=(1.0_dp/DiamP)*(DerNudUi(:)*llambda+Nu*DerllambdadUi(:))      
        DerHtdTsolid=dNudTsolid*llambda/DiamP

        TempInt=(Ht*TempP+lambdaS*TempS/Dist4Grad)/(Ht+lambdaS/Dist4Grad)

        DerTempInt_dUi(:)=((DerHt(:)*TempP+Ht*DerTempPdUi(:))*(Ht+lambdaS/Dist4Grad)-&
                          (Ht*TempP+lambdaS*TempS/Dist4Grad)*DerHt(:))/((Ht+lambdaS/Dist4Grad)**2)

        DerTempInt_dTm=((DerHtdTsolid*TempP+dlambdaSdT*TempS/Dist4Grad+lambdaS/Dist4Grad)*(Ht+lambdaS/Dist4Grad)-&
                        (Ht*TempP+lambdaS*TempS/Dist4Grad)*(DerHtdTsolid+dlambdaSdT/Dist4Grad))/((Ht+lambdaS/Dist4Grad)**2)

        ! FH       
        SrcP_Se=Ht*SurfCont*(TempInt-TempP)/(AreaP*dxLocP)
        DerSrcP_Se(:)=(SurfCont*(DerHt(:)*(TempInt-TempP)+Ht*(DerTempInt_dUi(:)-DerTempPdUi(:))))/(AreaP*dxLocP)

        SrcP_Sr=0.0_dp
        DerSrcP_Sr(:)=0.0_dp
        if(R_Correction) then    
          SrcP_Sr=SrcP_Se*(1.0_dp+Phi)
          DerSrcP_Sr(:)=DerSrcP_Se(:)*(1.0_dp+Phi)+SrcP_Se*DerPhi(:)
        endif

        DerSrcP_SedT=(DerHtdTsolid*(TempInt-TempP)+Ht*DerTempInt_dTm)*SurfCont/(AreaP*dxLocP)
        DerSrcP_SrdT=0.0_dp
        if(R_Correction) DerSrcP_SrdT=(1.0_dp+Phi)*(DerHtdTsolid*(TempInt-TempP)+Ht*DerTempInt_dTm)*SurfCont/(AreaP*dxLocP)

        ! HF
        SrcS=-SurfCont*lambdaS*(TempS-TempInt)/(Dist4Grad*rhoMS*cpMs*VolMS)
        DerSrcSdUi(:)=lambdaS*DerTempInt_dUi(:)*SurfCont/(rhoMS*cpMs*VolMS*Dist4Grad)

        DerSrcSdT=((dlambdaSdT*(TempS-TempInt)+lambdaS*(1.0_dp-DerTempInt_dTm))*cpMs-(lambdaS*(TempS-TempInt))*dcpMsdT)/(cpMs**2)
        DerSrcSdT=-DerSrcSdT*SurfCont/(rhoMS*VolMS*Dist4Grad) 

        ! To solid node --> considering FS accumulation
        !$omp critical(FS_state_accum)
        me%lk(2,i)%p%rhs_SrcS  = me%lk(2,i)%p%rhs_SrcS+SrcS
        me%lk(2,i)%p%DerSrcSdT = me%lk(2,i)%p%DerSrcSdT+DerSrcSdT       
        !$omp end critical(FS_state_accum)      

        ! Link part --> FS accumulation not required because different locations for links
        me%DerSrcP_SedT(i)=DerSrcP_SedT
        me%DerSrcP_SrdT(i)=DerSrcP_SrdT
        me%DerSrcSdUi(:,i)=DerSrcSdUi(:)        

      else if(me%lk(2,i)%p%typeS==trim("mesh2D")) then

        SrcP_Se=0.0_dp
        SrcP_Sr=0.0_dp
        DerSrcP_Se(:)=0.0_dp
        DerSrcP_Sr(:)=0.0_dp

        do ii=1,me%nb_2D_Ports(i)
          Dist4Grad = me%lk(1+ii,i)%p%Dist4Grad

          TempS     = me%lk(1+ii,i)%p%TempS
          lambdaS     = me%lk(1+ii,i)%p%lambdS
          dlambdaSdT     = me%lk(1+ii,i)%p%DerlambdS
           
          LengthCont = me%lk(1+ii,i)%p%LengthCont
          SurfCont = LengthCont*dxLocP
            
          Ht=me%Ht
          DerHt(:)=0.0_dp
          DerHtdTsolid=0.0_dp

          TempInt=(Ht*TempP+lambdaS*TempS/Dist4Grad)/(Ht+lambdaS/Dist4Grad)

          DerTempInt_dUi(:)=((DerHt(:)*TempP+Ht*DerTempPdUi(:))*(Ht+lambdaS/Dist4Grad)-&
                            (Ht*TempP+lambdaS*TempS/Dist4Grad)*DerHt(:))/((Ht+lambdaS/Dist4Grad)**2)
  
          DerTempInt_dTm=((DerHtdTsolid*TempP+dlambdaSdT*TempS/Dist4Grad+lambdaS/Dist4Grad)*(Ht+lambdaS/Dist4Grad)-&
                          (Ht*TempP+lambdaS*TempS/Dist4Grad)*(DerHtdTsolid+dlambdaSdT/Dist4Grad))/((Ht+lambdaS/Dist4Grad)**2)          

          SrcP_Se=SrcP_Se+Ht*SurfCont*(TempInt-TempP)/(AreaP*dxLocP)
          DerSrcP_Se(:)=DerSrcP_Se(:)+(SurfCont*(DerHt(:)*(TempInt-TempP)+Ht*(DerTempInt_dUi(:)-DerTempPdUi(:))))/(AreaP*dxLocP)
                  
          if(R_Correction) then    
            SrcP_Sr=SrcP_Se*(1.0_dp+Phi)
            DerSrcP_Sr(:)=DerSrcP_Se(:)*(1.0_dp+Phi)+SrcP_Se*DerPhi(:)
          endif      

          DerSrcP_SedT=(DerHtdTsolid*(TempInt-TempP)+Ht*DerTempInt_dTm)*SurfCont/(AreaP*dxLocP)
          DerSrcP_SrdT=0.0_dp
          if(R_Correction) DerSrcP_SrdT=(1.0_dp+Phi)*(DerHtdTsolid*(TempInt-TempP)+Ht*DerTempInt_dTm)*SurfCont/(AreaP*dxLocP)
  
          Flx2D=-LengthCont*lambdaS*(TempInt-TempS)/Dist4Grad
          DerFlx2DdUi(:)=-lambdaS*DerTempInt_dUi(:)*LengthCont/Dist4Grad
  
          derFlx2D_derCon=(dlambdaSdT*(TempS-TempInt)+lambdaS*(1.0_dp-DerTempInt_dTm))*LengthCont/Dist4Grad

          me%lk(1+ii,i)%p%Flx = Flx2D
          me%lk(1+ii,i)%p%derFlx_derCon = derFlx2D_derCon
    
          ! Link part --> SS_src accumulation not required because different locations for links
          me%DerSrcP_SedT2D(ii,i)=DerSrcP_SedT
          me%DerSrcP_SrdT2D(ii,i)=DerSrcP_SrdT          
          me%derFlx2DdUi(:,ii,i)=DerFlx2DdUi(:)
        enddo  

      endif

      ! To pipe node --> considering FS accumulation
      !$omp critical(FS_pipe_accum)
      me%lk(1,i)%p%rhs_SrcP_Se   = me%lk(1,i)%p%rhs_SrcP_Se+SrcP_Se
      me%lk(1,i)%p%rhs_SrcP_Sr   = me%lk(1,i)%p%rhs_SrcP_Sr+SrcP_Sr
      me%lk(1,i)%p%DerSrcP_Se(:) = me%lk(1,i)%p%DerSrcP_Se(:)+DerSrcP_Se(:)
      me%lk(1,i)%p%DerSrcP_Sr(:) = me%lk(1,i)%p%DerSrcP_Sr(:)+DerSrcP_Sr(:)
      !$omp end critical(FS_pipe_accum)

    enddo      
    
end subroutine FSlink_resolution_from_and_to_ports



subroutine links_FSlink_update(me,dt)
    type(FSlink_t), intent(inout) :: me
    real(dp), intent(in) :: dt

    integer :: i, ii
    real(dp) :: VarcentS

    do i = 1, me%sizeNodes

      if(me%lk(2,i)%p%typeS==trim("mesh2D")) then
        do ii=1,me%nb_2D_Ports(i)
          me%MatrixPS_2D(1,ii,i)=-dt*me%DerSrcP_SedT2D(ii,i)
          me%MatrixPS_2D(2,ii,i)=-dt*me%DerSrcP_SrdT2D(ii,i)
  
          VarcentS=1.0_dp/(me%lk(1+ii,i)%p%rhoS*me%lk(1+ii,i)%p%cpS)
          me%MatrixSP_2D(:,ii,i) = dt/me%lk(1+ii,i)%p%SurfS*VarcentS*me%derFlx2DdUi(:,ii,i)       
        enddo

      else

        me%MatrixPS(1,i)=-dt*me%DerSrcP_SedT(i)
        me%MatrixPS(2,i)=-dt*me%DerSrcP_SrdT(i)
        
        me%MatrixSP(:,i)=-dt*me%DerSrcSdUi(:,i)
      endif
    enddo

end subroutine links_FSlink_update

end module cmp_FSlink_calc_m

