! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_SSsrcLink_calc_m
    use cmp_SSsrcLink_init_m

    implicit none
  
contains


subroutine SSsrcLink_src_reinitialization(me)
  type(SSsrcLink_t), intent(inout) :: me

  integer :: i

  do i = 1, me%sizeNodes
    me%lk(1,i)%p%rhs_SrcS  = 0.0_dp
    me%lk(1,i)%p%DerSrcSdT = 0.0_dp
    me%lk(2,i)%p%rhs_SrcS  = 0.0_dp
    me%lk(2,i)%p%DerSrcSdT = 0.0_dp
  enddo

  ! TODO if mesh2D involves accumulation (case of a slice angle for example)

end subroutine SSsrcLink_src_reinitialization



subroutine SSsrcLink_resolution_from_and_to_ports(me)
  type(SSsrcLink_t), intent(inout) :: me

  real(dp), dimension(2) :: TempS,rhoS,VolS,Dist4Grad,cpS,dcpSdT,lambdaS,dlambdaSdT,SrcS,TempInt,DerSrcMC_dT
  real(dp), dimension(2,2) :: DerTempInt_dT,DerSrcS_dT
  real(dp) :: SurfCont,LengthCont,Fact,lambdaIns,Lins,Num,Den,SrcMC,Flx2D
  real(dp) :: derFlx2D_derMC, derFlx2D_derCon, Height
  integer :: ck, oth, ii, i
  integer, parameter :: MC  = 1
  integer, parameter :: M2D = 2

  do i = 1, me%sizeNodes

    if(me%lk(2,i)%p%typeS==trim("solid")) then

      SurfCont = me%SurfCont(i)
      Lins     = me%Lins(i)

      do ck=1,2
        rhoS(ck)       = me%lk(ck,i)%p%rhoS
        VolS(ck)       = me%lk(ck,i)%p%VolS
        Dist4Grad(ck)  = me%lk(ck,i)%p%Dist4Grad

        TempS(ck)      = me%lk(ck,i)%p%TempS
        lambdaS(ck)    = me%lk(ck,i)%p%lambdS
        dlambdaSdT(ck) = me%lk(ck,i)%p%DerlambdS
        cpS(ck)        = me%lk(ck,i)%p%cpS
        dcpSdT(ck)     = me%lk(ck,i)%p%dcpSdT        
      enddo

      lambdaIns=me%mat_ins%thermal_conductivity(0.5_dp*(TempS(1)+TempS(2)))
      ! TODO : How to deal with missing insulation mass?

      do ck=1,2
       if(ck==1) oth=2
       if(ck==2) oth=1
       Num=lambdaIns*(lambdaS(ck)*Dist4Grad(oth)*TempS(ck)+lambdaS(oth)*Dist4Grad(ck)*TempS(oth))+&
           lambdaS(ck)*lambdaS(oth)*Lins*TempS(ck)
       Den=lambdaIns*(lambdaS(ck)*Dist4Grad(oth)+lambdaS(oth)*Dist4Grad(ck))+lambdaS(ck)*lambdaS(oth)*Lins
       TempInt(ck)=Num/Den
       DerTempInt_dT(ck,ck)=((lambdaIns*(dlambdaSdT(ck)*Dist4Grad(oth)*TempS(ck)+lambdaS(ck)*Dist4Grad(oth))+lambdaS(oth)*Lins*&
                       (dlambdaSdT(ck)*TempS(ck)+lambdaS(ck)))*Den-Num*&
                       (lambdaIns*dlambdaSdT(ck)*Dist4Grad(oth)+lambdaS(oth)*Lins*dlambdaSdT(ck)))/(Den**2)   
       DerTempInt_dT(ck,oth)=((lambdaIns*(dlambdaSdT(oth)*Dist4Grad(ck)*TempS(oth)+lambdaS(oth)*Dist4Grad(ck))+&
                       lambdaS(ck)*Lins*TempS(ck)*dlambdaSdT(oth))*Den-Num*&
                       (lambdaIns*dlambdaSdT(oth)*Dist4Grad(ck)+lambdaS(ck)*Lins*dlambdaSdT(oth)))/(Den**2) 
      
       SrcS(ck)=SurfCont*lambdaS(ck)*(TempInt(ck)-TempS(ck))/(rhoS(ck)*cpS(ck)*VolS(ck)*Dist4Grad(ck))
       Fact=SurfCont/(rhoS(ck)*VolS(ck)*Dist4Grad(ck))
       DerSrcS_dT(ck,ck)=Fact*((dlambdaSdT(ck)*(TempInt(ck)-TempS(ck))+lambdaS(ck)*(DerTempInt_dT(ck,ck)-1.0_dp))*cpS(ck)-&
                        (lambdaS(ck)*(TempInt(ck)-TempS(ck)))*dcpSdT(ck))/(cpS(ck)**2)
       DerSrcS_dT(ck,oth)=Fact*lambdaS(ck)*DerTempInt_dT(ck,oth)/cpS(ck)
      enddo
  
      ! To solid node --> considering SS_src accumulation
      !$omp critical(SS_src_accum)
      me%lk(1,i)%p%rhs_SrcS  = me%lk(1,i)%p%rhs_SrcS+SrcS(1)
      me%lk(1,i)%p%DerSrcSdT = me%lk(1,i)%p%DerSrcSdT+DerSrcS_dT(1,1)
      me%lk(2,i)%p%rhs_SrcS  = me%lk(2,i)%p%rhs_SrcS+SrcS(2)
      me%lk(2,i)%p%DerSrcSdT = me%lk(2,i)%p%DerSrcSdT+DerSrcS_dT(2,2)    
      !$omp end critical(SS_src_accum)   

      ! Link part --> SS_src accumulation not required because different locations for links
      me%DerSrcS1dT2(i)=DerSrcS_dT(1,2)
      me%DerSrcS2dT1(i)=DerSrcS_dT(2,1)

    else if(me%lk(2,i)%p%typeS==trim("mesh2D")) then

      Lins = me%Lins(i)

      SrcMC=0.0_dp	
      DerSrcMC_dT=0.0_dp

      ! solid side
      rhoS(MC)       = me%lk(1,i)%p%rhoS
      VolS(MC)       = me%lk(1,i)%p%VolS
      Dist4Grad(MC)  = me%lk(1,i)%p%Dist4Grad
      Height         = me%lk(1,i)%p%Height

      TempS(MC)      = me%lk(1,i)%p%TempS
      lambdaS(MC)    = me%lk(1,i)%p%lambdS
      dlambdaSdT(MC) = me%lk(1,i)%p%DerlambdS   
      cpS(MC)        = me%lk(1,i)%p%cpS
      dcpSdT(MC)     = me%lk(1,i)%p%dcpSdT      

      do ii=1,me%nb_2D_Ports(i)

        ! 2D cell side
        Dist4Grad(M2D) = me%lk(1+ii,i)%p%Dist4Grad

        TempS(M2D)      = me%lk(1+ii,i)%p%TempS
        lambdaS(M2D)    = me%lk(1+ii,i)%p%lambdS
        dlambdaSdT(M2D) = me%lk(1+ii,i)%p%DerlambdS         

        lambdaIns=me%mat_ins%thermal_conductivity(0.5_dp*(TempS(MC)+TempS(M2D)))
  
        LengthCont = me%lk(1+ii,i)%p%LengthCont ! face of 1 cell
     
        do ck=1,2
          if(ck==1) oth=2
          if(ck==2) oth=1
          Num=lambdaIns*(lambdaS(ck)*Dist4Grad(oth)*TempS(ck)+lambdaS(oth)*Dist4Grad(ck)*TempS(oth))+&
              lambdaS(ck)*lambdaS(oth)*Lins*TempS(ck)
          Den=lambdaIns*(lambdaS(ck)*Dist4Grad(oth)+lambdaS(oth)*Dist4Grad(ck))+lambdaS(ck)*lambdaS(oth)*Lins
          TempInt(ck)=Num/Den
          DerTempInt_dT(ck,ck)=((lambdaIns*(dlambdaSdT(ck)*Dist4Grad(oth)*TempS(ck)+lambdaS(ck)*Dist4Grad(oth))+lambdaS(oth)*Lins*&
                          (dlambdaSdT(ck)*TempS(ck)+lambdaS(ck)))*Den-Num*&
                          (lambdaIns*dlambdaSdT(ck)*Dist4Grad(oth)+lambdaS(oth)*Lins*dlambdaSdT(ck)))/(Den**2)   
          DerTempInt_dT(ck,oth)=((lambdaIns*(dlambdaSdT(oth)*Dist4Grad(ck)*TempS(oth)+lambdaS(oth)*Dist4Grad(ck))+&
                          lambdaS(ck)*Lins*TempS(ck)*dlambdaSdT(oth))*Den-Num*&
                          (lambdaIns*dlambdaSdT(oth)*Dist4Grad(ck)+lambdaS(ck)*Lins*dlambdaSdT(oth)))/(Den**2) 
        enddo

        SrcMC = SrcMC + LengthCont*Height*lambdaS(MC)*(TempInt(MC)-TempS(MC))/(rhoS(MC)*cpS(MC)*VolS(MC)*Dist4Grad(MC))
        Flx2D = -LengthCont*lambdaS(M2D)*(TempInt(M2D)-TempS(M2D))/Dist4Grad(M2D)
  
        Fact=LengthCont*Height/(rhoS(MC)*VolS(MC)*Dist4Grad(MC))
        DerSrcMC_dT(MC)  = DerSrcMC_dT(MC) + Fact*((dlambdaSdT(MC)*(TempInt(MC)-TempS(MC))+lambdaS(MC)*&
                             (DerTempInt_dT(MC,MC)-1.0_dp))*cpS(MC)-(lambdaS(MC)*(TempInt(MC)-TempS(MC)))*dcpSdT(MC))/(cpS(MC)**2)
        DerSrcMC_dT(M2D) = Fact*lambdaS(MC)*DerTempInt_dT(MC,M2D)/cpS(MC)  
  
        Fact=-LengthCont/Dist4Grad(M2D)
        derFlx2D_derCon=Fact*(dlambdaSdT(M2D)*(TempInt(M2D)-TempS(M2D))+lambdaS(M2D)*(DerTempInt_dT(M2D,M2D)-1.0_dp))
        derFlx2D_derMC=Fact*lambdaS(M2D)*DerTempInt_dT(M2D,MC)
  
        me%lk(1+ii,i)%p%Flx = Flx2D
        me%lk(1+ii,i)%p%derFlx_derCon = derFlx2D_derCon
  
        ! Link part --> SS_src accumulation not required because different locations for links
        me%Der4LkArr(1,1+ii,i)=DerSrcMC_dT(M2D)
        me%Der4LkArr(1+ii,1,i)=derFlx2D_derMC

      enddo

      ! To solid node --> considering SS_src accumulation
      !$omp critical(SS_src_accum)
      me%lk(1,i)%p%rhs_SrcS  = me%lk(1,i)%p%rhs_SrcS+SrcMC
      me%lk(1,i)%p%DerSrcSdT = me%lk(1,i)%p%DerSrcSdT+DerSrcMC_dT(MC)
      !$omp end critical(SS_src_accum)

    endif
  enddo
    
end subroutine SSsrcLink_resolution_from_and_to_ports



subroutine links_SSsrcLink_update(me,dt)
    type(SSsrcLink_t), intent(inout) :: me
    real(dp), intent(in) :: dt

    integer  :: ii, i
    real(dp) :: VarcentS

    do i = 1, me%sizeNodes
      if(me%lk(2,i)%p%typeS==trim("solid")) then
  
        me%Lk12(i) = -dt*me%DerSrcS1dT2(i)
        me%Lk21(i) = -dt*me%DerSrcS2dT1(i)
  
      else if(me%lk(2,i)%p%typeS==trim("mesh2D")) then
        
        do ii=1,me%nb_2D_Ports(i)
          me%LkArr(1,1+ii,i) = -dt*me%Der4LkArr(1,1+ii,i) 
  
          VarcentS=1.0_dp/(me%lk(1+ii,i)%p%rhoS*me%lk(1+ii,i)%p%cpS)
          me%LkArr(1+ii,1,i) = dt/me%lk(1+ii,i)%p%SurfS*VarcentS*me%Der4LkArr(1+ii,1,i)
        enddo
  
      endif
    enddo

end subroutine links_SSsrcLink_update

end module cmp_SSsrcLink_calc_m