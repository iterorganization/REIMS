! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_channel_calc_m
    use cmp_channel_init_m
    use krn_global_tools_m
    use cmp_channel_source_terms_m
    use krn_simulation_m
    implicit none

contains

subroutine flux_from_FF_flux_ports_self(me)
    type(channel_t), intent(inout) :: me

    ! in component
    me%flxHe%Cons(:,0)=me%HeProp%in%flx(:)
    me%flxHe%VitTNC(0)=me%HeProp%in%vit
    me%flxHe%Cons_DerdUR(:,:,0)=me%HeProp%in%derFlx_derCon(:,:)
    me%flxHe%VitTNC_DerdUR(:,0)=me%HeProp%in%derVit_derCon(:)

    ! out component
    me%flxHe%Cons(:,me%HeProp%NbCells)=me%HeProp%out%flx(:)
    me%flxHe%VitTNC(me%HeProp%NbCells)=me%HeProp%out%vit
    me%flxHe%Cons_DerdUL(:,:,me%HeProp%NbCells)=me%HeProp%out%derFlx_derCon(:,:)
    me%flxHe%VitTNC_DerdUL(:,me%HeProp%NbCells)=me%HeProp%out%derVit_derCon(:)
        
end subroutine flux_from_FF_flux_ports_self



subroutine channels_FF_flux_port_comm(me)
    type(channel_t), intent(inout) :: me

    me%HeProp%in%Prim(:)         = me%StVar%He_pr(:,1)
    me%HeProp%in%Cons(:)         = me%StVar%He_cs(:,1)
    me%HeProp%in%Su_fric        = me%srcPP%Su(1)
    me%HeProp%in%DerSu_fric(:)  = me%srcPP%DerSu(:,1)
    me%HeProp%in%Sr_fric        = me%srcPP%Sr(1)
    me%HeProp%in%DerSr_fric(:)  = me%srcPP%DerSr(:,1)
        
    me%HeProp%out%Prim(:)        = me%StVar%He_pr(:,me%HeProp%NbCells)
    me%HeProp%out%Cons(:)        = me%StVar%He_cs(:,me%HeProp%NbCells)
    me%HeProp%out%Su_fric       = me%srcPP%Su(me%HeProp%NbCells)
    me%HeProp%out%DerSu_fric(:) = me%srcPP%DerSu(:,me%HeProp%NbCells)
    me%HeProp%out%Sr_fric       = me%srcPP%Sr(me%HeProp%NbCells)
    me%HeProp%out%DerSr_fric(:) = me%srcPP%DerSr(:,me%HeProp%NbCells)
        
end subroutine channels_FF_flux_port_comm



subroutine source_from_FS_ports_self_channel(me,dt)
    type(channel_t), intent(inout) :: me
    real(dp), intent(in) :: dt

    integer :: i
    
    if(me%HeProp%ExtHeating=='link') then
      do i=1,me%HeProp%NbCells
        me%big%bmx(1:Nb_VarC,Con_Ene,i)=me%big%bmx(1:Nb_VarC,Con_Ene,i)-dt*me%HeProp%thermP(i)%p%DerSrcP_Se(:)
        me%big%bmx(1:Nb_VarC,Con_R,i)=me%big%bmx(1:Nb_VarC,Con_R,i)-dt*me%HeProp%thermP(i)%p%DerSrcP_Sr(:)
    
        me%big%bvx(Con_Ene,i)=me%big%bvx(Con_Ene,i)+dt*me%HeProp%thermP(i)%p%rhs_SrcP_Se
        me%big%bvx(Con_R,i)=me%big%bvx(Con_R,i)+dt*me%HeProp%thermP(i)%p%rhs_SrcP_Sr
      enddo
    endif

end subroutine source_from_FS_ports_self_channel



subroutine source_from_FF_src_ports_self_channel(me,dt)
    type(channel_t), intent(inout) :: me
    real(dp), intent(in) :: dt

    integer :: i
    
    if(me%HeProp%FFsrcLink) then
      do i=1,me%HeProp%NbCells
        me%big%bmx(1:Nb_VarC,Con_Ene,i)=me%big%bmx(1:Nb_VarC,Con_Ene,i)-dt*me%HeProp%thermF(i)%p%DerSe(:)
        me%big%bmx(1:Nb_VarC,Con_R,i)=me%big%bmx(1:Nb_VarC,Con_R,i)-dt*me%HeProp%thermF(i)%p%DerSr(:)
    
        me%big%bvx(Con_Ene,i)=me%big%bvx(Con_Ene,i)+dt*me%HeProp%thermF(i)%p%rhs_Se
        me%big%bvx(Con_R,i)=me%big%bvx(Con_R,i)+dt*me%HeProp%thermF(i)%p%rhs_Sr
      enddo
    endif

end subroutine source_from_FF_src_ports_self_channel



subroutine from_sol_to_prim_channel(me)    
    type(channel_t), intent(inout) :: me

    real(dp) :: PoutCP,ToutCP,coutCP,rocheck
    integer :: jj

    do jj=1,me%HeProp%NbCells
        me%StVar%He_cs(:,jj)=me%StVar%He_cs(:,jj)+me%big%bvx(1:Nb_VarC,jj) ! update of cons by the solution
      
        me%StVar%He_pr(Pri_ro,jj)=me%StVar%He_cs(Con_Mas,jj)

        rocheck=me%StVar%He_pr(Pri_ro,jj)
        if(ieee_is_nan(rocheck) .or. rocheck<0.001_dp .or. rocheck>160.0_dp) then
            call set_error('unphysical density in channel '//trim(me%HeProp%name))
            return
        endif

        me%StVar%He_pr(Pri_u,jj)=me%StVar%He_cs(Con_Qdm,jj)/me%StVar%He_cs(Con_Mas,jj)
        me%StVar%He_pr(Pri_e,jj)=me%StVar%He_cs(Con_Ene,jj)/me%StVar%He_cs(Con_Mas,jj)-&
                                       0.5_dp*me%StVar%He_pr(Pri_u,jj)**2

        me%StVar%He_pr(Pri_R,jj)=me%StVar%He_cs(Con_R,jj)

        if(R_Correction) then
          call state_roE_withR(me%StVar%He_pr(Pri_ro,jj), me%StVar%He_pr(Pri_e,jj), &
                            me%StVar%He_pr(Pri_R,jj), PoutCP, ToutCP, coutCP)
        else
          call state_roE(me%StVar%He_pr(Pri_ro,jj), me%StVar%He_pr(Pri_e,jj), PoutCP, ToutCP, coutCP)
        endif
        if (ToutCP >= T_He_max .or. ToutCP <= T_He_min) then
            call set_error('temperature out of EOS range in channel '//trim(me%HeProp%name))
            return
        end if
        me%StVar%He_pr(Pri_p,jj)=PoutCP
        me%StVar%He_pr(Pri_T,jj)=ToutCP
        me%StVar%He_pr(Pri_c,jj)=coutCP
    enddo
      
end subroutine from_sol_to_prim_channel



subroutine channels_to_FF_src_port_comm(me)
    type(channel_t), intent(inout) :: me
    
    integer :: jj

    if(me%HeProp%FFsrcLink) then
      do jj=1,me%HeProp%NbCells
        me%HeProp%thermF(jj)%p%Prim(:) = me%StVar%He_pr(:,jj)
        me%HeProp%thermF(jj)%p%Cons(:) = me%StVar%He_cs(:,jj)
      enddo
    endif
        
end subroutine channels_to_FF_src_port_comm



subroutine FF_src_port_to_channels_comm(me)
    type(channel_t), intent(inout) :: me
    
    integer :: jj

    if(me%HeProp%FFsrcLink) then
      do jj=1,me%HeProp%NbCells
        me%StVar%He_pr(:,jj) = me%HeProp%thermF(jj)%p%Prim(:)
        me%StVar%He_cs(:,jj) = me%HeProp%thermF(jj)%p%Cons(:)
      enddo
    endif
        
end subroutine FF_src_port_to_channels_comm



subroutine last_tasks_for_channels(me,sim)    
    type(channel_t), intent(inout) :: me
    type(simulation_t), intent(in) :: sim

    real(dp) :: mdot,HHH,Rout,Pout,eout,cout
    integer :: jj

    me%err = 0_dp
    me%err_den = 0_dp

    do jj=1,me%HeProp%NbCells

        ! Write to HDF5
        me%hdf%data(jj,1) = me%StVar%He_pr(Pri_p,jj)/1.0e5_dp                               ! He Pressure (Bar)
        me%hdf%data(jj,2) = me%StVar%He_pr(Pri_ro,jj)                                       ! He Density (kg/m3)
        me%hdf%data(jj,3) = me%StVar%He_pr(Pri_u,jj)                                        ! He Velocity (m/s)
        me%hdf%data(jj,4) = me%HeProp%Area*me%StVar%He_pr(Pri_u,jj)*me%StVar%He_pr(Pri_ro,jj) ! He MassFlowRate (kg/s)
        me%hdf%data(jj,5) = me%StVar%He_pr(Pri_T,jj)                                        ! He Temperature (K)
        call state_roT(me%StVar%He_pr(Pri_ro,jj), me%StVar%He_pr(Pri_T,jj), Rout, Pout, eout, cout)
        mdot=me%HeProp%Area*me%StVar%He_pr(Pri_u,jj)*me%StVar%He_pr(Pri_ro,jj)
        HHH=eout+Pout/me%StVar%He_pr(Pri_ro,jj)+0.5_dp*(me%StVar%He_pr(Pri_u,jj)**2)
        me%hdf%data(jj,6)=mdot*HHH                                                               ! He Mdot*H
        me%hdf%data(jj,7)=me%StVar%He_pr(Pri_R,jj)                                          ! He R eq
        me%hdf%data(jj,8)=me%srcPP%Su(jj)

        ! For step management
        if(sim%explicit) then
            me%err = me%err + (me%StVar%He_pr(Pri_p,jj) - me%StVarOld%He_pr(Pri_p,jj))**2
            me%err_den = me%err_den + me%StVarOld%He_pr(Pri_p,jj)**2
        else
            me%err = me%err + ( (1.0_dp/6.0_dp) * (sim%dtPrev1+sim%dt) * ( &
                (me%StVar    %He_cs(2,jj)-me%StVarOld %He_cs(2,jj))               / sim%dt - &
                (me%StVarOld %He_cs(2,jj)-me%StVarOld2%He_cs(2,jj)) * (sim%dt+sim%dtPrev1) / sim%dtPrev1**2 + &
                (me%StVarOld2%He_cs(2,jj)-me%StVarOld3%He_cs(2,jj)) *  sim%dt/sim%dtPrev1  / sim%dtPrev2 ))**2
            me%err_den = me%err_den + (me%StVarOld%He_cs(2,jj))**2                      
        endif

    enddo
      
end subroutine last_tasks_for_channels

subroutine explicit_scheme_for_fluxes_channel(me,dt)          
    type(channel_t), intent(inout) :: me
    real(dp),     intent(in)    :: dt
    
    real(dp) :: Varcentree, VitL, VitR
    real(dp) :: FluxL(Nb_VarC), FluxR(Nb_VarC)
    integer :: i

    do i=1,me%HeProp%NbCells
        Varcentree=0.0_dp
        if(R_Correction) Varcentree=me%StVar%He_pr(Pri_ro,i)*(me%StVar%He_pr(Pri_c,i)**2)

        FluxL(:)=me%flxHe%Cons(:,i-1)
        VitL=me%flxHe%VitTNC(i-1)
        FluxR(:)=me%flxHe%Cons(:,i)
        VitR=me%flxHe%VitTNC(i)

        me%StVar%He_cs(:,i)=me%StVar%He_cs(:,i)-(dt/me%HeProp%dxloc(i))*(FluxR(:)-FluxL(:))
        me%StVar%He_cs(Con_R,i)=me%StVar%He_cs(Con_R,i)-(dt/me%HeProp%dxloc(i))*Varcentree*(VitR-VitL)
    enddo
     
end subroutine explicit_scheme_for_fluxes_channel       

subroutine full_physics_definition_channel(me,sim)
    type(channel_t),       intent(inout) :: me
    type(simulation_t), intent(in)    :: sim

      real(dp) :: wwn    
      real(dp) :: Varcentree_i,Varcentree_iplus !! Centered variable ro c2
      real(dp), dimension(Nb_VarC) :: DerVarcentree
      real(dp), dimension(Nb_VarC) :: dcSQ
      real(dp), dimension(Nb_VarC) :: DerRodUi
      real(dp) :: dcSQdro,dcSQdp,dxloc,dxlocA
      integer :: i,NbCells

      wwn = sim%dt/sim%dtPrev1
      me%big%amx=0.0_dp;me%big%bmx=0.0_dp;me%big%cmx=0.0_dp;me%big%bvx=0.0_dp
  
      DerRodUi(:)=0.0_dp
      DerRodUi(Con_Mas)=1.0_dp

      NbCells=me%HeProp%NbCells
  
      do i=1,NbCells

        dxloc=me%HeProp%dxloc(i)
        if(i/=NbCells) dxlocA=me%HeProp%dxloc(i+1)

        if(i/=NbCells) me%big%amx(:,:,i)=-(sim%dt/dxlocA)*(me%flxHe%Cons_DerdUL(:,:,i)) ! 1st element --> cell 2
        me%big%bmx(:,:,i)=id_4x4(:,:)+(sim%dt/dxloc)*(me%flxHe%Cons_DerdUL(:,:,i)-&
                                     me%flxHe%Cons_DerdUR(:,:,i-1))
        if(i/=NbCells) me%big%cmx(:,:,i)=(sim%dt/dxloc)*me%flxHe%Cons_DerdUR(:,:,i)   

        me%big%bvx(1:Nb_VarC,i)=-(sim%dt/dxloc)*(me%flxHe%Cons(:,i)-me%flxHe%Cons(:,i-1))
        
        ! Correction with non-conservative terms      
        if(R_Correction) then
          call dc2_roT(me%StVar%He_pr(Pri_ro,i), me%StVar%He_pr(Pri_T,i), dcSQdro, dcSQdp)

          dcSQ(Con_Mas)=dcSQdro-0.5_dp*dcSQdp*(me%StVar%He_cs(Con_Qdm,i)**2)/(me%StVar%He_cs(Con_Mas,i)**2)
          dcSQ(Con_Qdm)=dcSQdp*me%StVar%He_cs(Con_Qdm,i)/me%StVar%He_cs(Con_Mas,i)
          dcSQ(Con_Ene)=-dcSQdp
          dcSQ(Con_R)=dcSQdp

          Varcentree_i=me%StVar%He_pr(Pri_ro,i)*(me%StVar%He_pr(Pri_c,i)**2)
          if(i/=NbCells) Varcentree_iplus=me%StVar%He_pr(Pri_ro,i+1)*(me%StVar%He_pr(Pri_c,i+1)**2)
          DerVarcentree(:)=DerRodUi(:)*(me%StVar%He_pr(Pri_c,i)**2)+me%StVar%He_pr(Pri_ro,i)*dcSQ(:)
        else    
          Varcentree_i=0.0_dp
          if(i/=NbCells) Varcentree_iplus=0.0_dp
          DerVarcentree(:)=0.0_dp
        endif
    
        if(i/=NbCells) me%big%amx(1:Nb_VarC,Con_R,i)=me%big%amx(1:Nb_VarC,Con_R,i)-(sim%dt/dxlocA)*Varcentree_iplus*&
                                                         me%flxHe%VitTNC_DerdUL(:,i)
        me%big%bmx(1:Nb_VarC,Con_R,i)=me%big%bmx(1:Nb_VarC,Con_R,i)+(sim%dt/dxloc)*(Varcentree_i*&
                                           (me%flxHe%VitTNC_DerdUL(:,i)-me%flxHe%VitTNC_DerdUR(:,i-1))+&
                                           (me%flxHe%VitTNC(i)-me%flxHe%VitTNC(i-1))*DerVarcentree(:))
        if(i/=NbCells) me%big%cmx(1:Nb_VarC,Con_R,i)=me%big%cmx(1:Nb_VarC,Con_R,i)+(sim%dt/dxloc)*Varcentree_i*&
                                                          me%flxHe%VitTNC_DerdUR(:,i)
        me%big%bvx(Con_R,i)=me%big%bvx(Con_R,i)-(sim%dt/dxloc)*Varcentree_i*(me%flxHe%VitTNC(i)-&
                                 me%flxHe%VitTNC(i-1))                 

      enddo

      call source_term_definition_channel(me,sim)
    
      do i=1,NbCells
        ! Modification of the diagonal terms 
          me%big%bmx(1:Nb_VarC,1:Nb_VarC,i)=me%big%bmx(1:Nb_VarC,1:Nb_VarC,i)-id_4x4(:,:)+&
                                                 ((1.0_dp+2.0_dp*wwn)/(1.0_dp+wwn))*id_4x4(:,:)
          
        ! Modification of the vector B by adding time step n-1 contribution
          me%big%bvx(1:Nb_VarC,i)=me%big%bvx(1:Nb_VarC,i)+((wwn**2)/(1.0_dp+wwn))*&
                                       (me%StVar%He_cs(:,i)-me%StVarOld2%He_cs(:,i))
      enddo

end subroutine full_physics_definition_channel

  

subroutine source_term_definition_channel(me,sim)
    type(channel_t), intent(inout) :: me
    type(simulation_t), intent(in) :: sim

    real(dp) :: QorT(me%HeProp%NbCells)
  
        integer :: i
    
        if(sim%explicit) then
          me%big%bmx=0.0_dp; me%big%bvx=0.0_dp
          
          do i=1,me%HeProp%NbCells
            me%big%bmx(1:Nb_VarC,1:Nb_VarC,i)=id_4x4(:,:)            
          enddo
        endif

        if(me%HeProp%ExtHeating=='flux' .or. me%HeProp%ExtHeating=='temp') QorT=me%HeProp%Temp_or_Q_Wall%v1d()

        do i=1,me%HeProp%NbCells
          ! Adding of source term contributions
            call SourceTerms_NOTfriction_channel(me,i,QorT(i))
      
            me%big%bmx(1:Nb_VarC,Con_Qdm,i)=me%big%bmx(1:Nb_VarC,Con_Qdm,i)-sim%dt*(me%srcPP%DerSu(:,i))
            me%big%bmx(1:Nb_VarC,Con_Ene,i)=me%big%bmx(1:Nb_VarC,Con_Ene,i)-sim%dt*(me%srcPP%DerSe(:))
            me%big%bmx(1:Nb_VarC,Con_R,i)=me%big%bmx(1:Nb_VarC,Con_R,i)-sim%dt*(me%srcPP%DerSr(:,i))
      
            me%big%bvx(Con_Qdm,i)=me%big%bvx(Con_Qdm,i)+sim%dt*(me%srcPP%Su(i))
            me%big%bvx(Con_Ene,i)=me%big%bvx(Con_Ene,i)+sim%dt*(me%srcPP%Se)
            me%big%bvx(Con_R,i)=me%big%bvx(Con_R,i)+sim%dt*(me%srcPP%Sr(i))    
        enddo
        
end subroutine source_term_definition_channel   

end module cmp_channel_calc_m
