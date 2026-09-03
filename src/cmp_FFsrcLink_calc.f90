! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_FFsrcLink_calc_m
  use cmp_FFsrcLink_init_m
  use lib_He_thermo_m
    implicit none

  
contains


subroutine FFsrcLink_resolution_from_and_to_ports(me)
    type(FFsrcLink_t), intent(inout) :: me

    real(dp) :: Pra,Re,llambda,Nu,dNudRe,dNudPr,dNudTemp,dNudTother,HtExch,WetPerChan
    real(dp), dimension(Nb_VarC) :: dPrdUi,dRedUi,DerllambdadUi,DerNudUi,DerNudUiOther
    real(dp), dimension(2) :: TempP,Phi,AreaP,DiamP,Ht,TempOther
    real(dp), dimension(Nb_VarC,2) :: DerPhi,DerTempPdUi,DerHt,DerTempOtherdUiOther,DerHtDerOther
    real(dp), dimension(2) :: Sexch_Ene,Sexch_R
    real(dp), dimension(Nb_VarC,2) :: DerHtExchdUi,DerHtExchdUiOther,DerSexch_EnedUi,DerSexch_EnedOther
    real(dp), dimension(Nb_VarC,2) :: DerSexch_RdUi,DerSexch_RdOther
    integer :: k,i

    ! Only thermal exchanges between fluids are present in this routine

    WetPerChan = me%WetPerChan

    do i = 1, me%sizeNodes
      do k = 1, 2
        Pra              = me%lk(k,i)%p%Pra
        Re               = me%lk(k,i)%p%Re
        dPrdUi(:)        = me%lk(k,i)%p%dPrdUi(:)
        dRedUi(:)        = me%lk(k,i)%p%dRedUi(:)
        llambda          = me%lk(k,i)%p%llambda
        DerllambdadUi(:) = me%lk(k,i)%p%DerllambdadUi(:)
        TempP(k)         = me%lk(k,i)%p%TempP
        DerTempPdUi(:,k) = me%lk(k,i)%p%DerTempPdUi(:)
        Phi(k)           = me%lk(k,i)%p%Phi
        DerPhi(:,k)      = me%lk(k,i)%p%DerPhi(:)
        AreaP(k)         = me%lk(k,i)%p%AreaP
        DiamP(k)         = me%lk(k,i)%p%DiamP

        if(k==1) then
          TempOther(k)=me%lk(2,i)%p%TempP
          DerTempOtherdUiOther(:,k)=me%lk(2,i)%p%DerTempPdUi(:)
        else
          TempOther(k)=me%lk(1,i)%p%TempP
          DerTempOtherdUiOther(:,k)=me%lk(1,i)%p%DerTempPdUi(:)
        endif
        Nu=me%nuss%nusselt(Re,Pra,TempOther(k),TempP(k))
        Ht(k)=Nu*llambda/DiamP(k)

        dNudRe=me%nuss%nusselt_der_Re(Re,Pra,TempOther(k),TempP(k))
        dNudPr=me%nuss%nusselt_der_Pra(Re,Pra,TempOther(k),TempP(k))
        dNudTemp=me%nuss%nusselt_der_T(Re,Pra,TempOther(k),TempP(k))
        DerNudUi(:)=dNudRe*dRedUi(:)+dNudPr*dPrdUi(:)+dNudTemp*DerTempPdUi(:,k)
        DerHt(:,k)=(1.0_dp/DiamP(k))*(DerNudUi(:)*llambda+Nu*DerllambdadUi(:))

        dNudTother=me%nuss%nusselt_der_Tother(Re,Pra,TempOther(k),TempP(k))
        DerNudUiOther(:)=dNudTother*DerTempOtherdUiOther(:,k)
        DerHtDerOther(:,k)=DerNudUiOther(:)*llambda/DiamP(k)     
      enddo

      if (min(Ht(1), Ht(2)) < 1.0e-15_dp) then
        HtExch = 0.0_dp
        DerHtExchdUi = 0.0_dp
        DerHtExchdUiOther = 0.0_dp
      else
        HtExch=1.0_dp/(1.0_dp/Ht(1)+1.0_dp/Ht(2))
        do k = 1, 2
          DerHtExchdUi(:,k)=(DerHt(:,k)/(Ht(k)**2))/((1.0_dp/Ht(1)+1.0_dp/Ht(2))**2)
          DerHtExchdUiOther(:,k)=(DerHtDerOther(:,k)/(Ht(k)**2))/((1.0_dp/Ht(1)+1.0_dp/Ht(2))**2)
        enddo
      end if
      do k = 1, 2
        Sexch_Ene(k)=WetPerChan*HtExch*(TempOther(k)-TempP(k))/AreaP(k)

        DerSexch_EnedUi(:,k)=(WetPerChan/AreaP(k))*(DerHtExchdUi(:,k)*(TempOther(k)-TempP(k))-HtExch*DerTempPdUi(:,k))
        DerSexch_EnedOther(:,k)=(WetPerChan/AreaP(k))*(DerHtExchdUiOther(:,k)*(TempOther(k)-TempP(k))+&
                                HtExch*DerTempOtherdUiOther(:,k))

        if(R_Correction) then
          Sexch_R(k)=Sexch_Ene(k)*(1.0_dp+Phi(k))
          DerSexch_RdUi(:,k)=DerSexch_EnedUi(:,k)*(1.0_dp+Phi(k))+Sexch_Ene(k)*DerPhi(:,k)
          DerSexch_RdOther(:,k)=DerSexch_EnedOther(:,k)*(1.0_dp+Phi(k))
        else
          Sexch_R(k)=0.0_dp
          DerSexch_RdUi(:,k)=0.0_dp
          DerSexch_RdOther(:,k)=0.0_dp
        endif

        ! To pipe nodes
        me%lk(k,i)%p%rhs_Se=Sexch_Ene(k)
        me%lk(k,i)%p%rhs_Sr=Sexch_R(k)
        me%lk(k,i)%p%DerSe(:)=DerSexch_EnedUi(:,k)
        me%lk(k,i)%p%DerSr(:)=DerSexch_RdUi(:,k)
        ! Link part
        me%DerSedOther(:,k,i)=DerSexch_EnedOther(:,k)
        me%DerSrdOther(:,k,i)=DerSexch_RdOther(:,k)      
      enddo
    enddo
   
end subroutine FFsrcLink_resolution_from_and_to_ports



subroutine links_FFsrcLink_update(me,dt)
    type(FFsrcLink_t), intent(inout) :: me
    real(dp), intent(in) :: dt

    integer :: k,i

    do i = 1, me%sizeNodes
      do k = 1, 2
        me%Matrix(:,1,k,i)=-dt*me%DerSedOther(:,k,i)
        me%Matrix(:,2,k,i)=-dt*me%DerSrdOther(:,k,i)
      enddo

      if (associated(me%MxUp(1,i)%p)) then
        call accumulate(me%Matrix(:,:,1,i),me%MxUp(:,i))
        call accumulate(me%Matrix(:,:,2,i),me%MxDown(:,i))
      endif
    enddo

end subroutine links_FFsrcLink_update



subroutine FFsrcLink_Prelax_resolution_from_and_to_ports(me)
    type(FFsrcLink_t), intent(inout) :: me

    real(dp), dimension(Nb_VarP,2) :: Pr
    real(dp), dimension(Nb_VarC,2) :: Cs
    real(dp), dimension(2) :: Kappak, rok, pk, rokstar, tkstar, rCorrkstar, ekstar, ckstar
    real(dp), dimension(2) :: ukstar, Vitk, Area
    real(dp) :: VitBar, hBar, pstar
    real(dp) :: ccv,ccp,mmu,llambda,dPdT_Ro,hk,Phik,T_temp
    integer :: k,i

    do i = 1, me%sizeNodes
      do k = 1, 2
        Pr(:,k)=me%lk(k,i)%p%Prim(:)
        Cs(:,k)=me%lk(k,i)%p%Cons(:)
        Area(k)=me%lk(k,i)%p%AreaP
      enddo

      if(Pr(Pri_p,1)>=Pr(Pri_p,2)) then
        VitBar=Pr(Pri_u,1)
        if(R_Correction) then
          hBar=Pr(Pri_R,1)/Pr(Pri_ro,1)
        else
          hBar=Pr(Pri_e,1)+Pr(Pri_p,1)/Pr(Pri_ro,1)
        endif
      else
        VitBar=Pr(Pri_u,2)
        if(R_Correction) then
          hBar=Pr(Pri_R,2)/Pr(Pri_ro,2)
        else
          hBar=Pr(Pri_e,2)+Pr(Pri_p,2)/Pr(Pri_ro,2)
        endif
      endif

      do k = 1, 2
        call he_prop(Pr(Pri_ro,k),Pr(Pri_T,k),ccv,ccp,mmu,llambda,dPdT_Ro)
        if (sim_error > 0) return
        rok(k)=Pr(Pri_ro,k)
        pk(k)=Pr(Pri_p,k)
        if(R_Correction) then
          hk=Pr(Pri_R,k)/Pr(Pri_ro,k)
        else
          hk=Pr(Pri_e,k)+Pr(Pri_p,k)/Pr(Pri_ro,k)
        endif
        Vitk(k)=Pr(Pri_u,k)
        Phik=1.0_dp/(Pr(Pri_ro,k)*(ccv/dPdT_Ro))
        Kappak(k)=Phik*((hBar-hk)+0.5_dp*((VitBar-Vitk(k))**2))/Area(k)+(Pr(Pri_c,k)**2)/Area(k)
      enddo  

      pstar=0.5_dp*(pk(1)+pk(2))+0.5_dp*(pk(1)-pk(2))*(Kappak(2)-Kappak(1))/(Kappak(1)+Kappak(2))

      do k = 1, 2
        rokstar(k)=rok(k)+(pstar-pk(k))/(Area(k)*Kappak(k))
        T_temp=T_roP(rokstar(k),pstar)
        tkstar(k)=T_temp
        call state_roT(rokstar(k), tkstar(k), rCorrkstar(k), Pstar, ekstar(k), ckstar(k))
        ukstar(k)=Vitk(k)+(VitBar-Vitk(k))*(rokstar(k)-rok(k))/rok(k)
  
        Pr(Pri_ro,k)=rokstar(k)
        Pr(Pri_u,k)=ukstar(k)
        Pr(Pri_p,k)=Pstar
        Pr(Pri_e,k)=ekstar(k)
        Pr(Pri_T,k)=tkstar(k)
        Pr(Pri_c,k)=ckstar(k)
        if(R_Correction) Pr(Pri_R,k)=rCorrkstar(k)
  
        Cs(Con_Mas,k)=rokstar(k)
        Cs(Con_Qdm,k)=rokstar(k)*ukstar(k)
        Cs(Con_Ene,k)=rokstar(k)*(ekstar(k)+0.5_dp*(ukstar(k)**2))
        if(R_Correction) Cs(Con_R,k)=rCorrkstar(k)
  
        me%lk(k,i)%p%Prim(:)=Pr(:,k)
        me%lk(k,i)%p%Cons(:)=Cs(:,k)
      enddo
    enddo

end subroutine FFsrcLink_Prelax_resolution_from_and_to_ports

end module cmp_FFsrcLink_calc_m
