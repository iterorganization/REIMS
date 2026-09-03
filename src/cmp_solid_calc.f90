! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_solid_calc_m
    use cmp_solid_init_m
    use krn_simulation_m
    implicit none

contains



subroutine solid_scenario_update(me)
    type(solid_t), intent(inout) :: me

    me%Q_ext = me%Q_ext_load%v1d()

end subroutine solid_scenario_update



subroutine source_from_FS_ports_self_solid(me,dt)
    type(solid_t), intent(inout) :: me
    real(dp), intent(in) :: dt

    integer :: i

    do i=1,me%MC_Prop%NbCells  
      if(me%MC_Prop%FSlink(i)) then        
        me%big%bmS(i)=me%big%bmS(i)-dt*me%MC_Prop%thermP(me%MC_Prop%idxFSlk(i))%p%DerSrcSdT
        me%big%bvS(i)=me%big%bvS(i)+dt*me%MC_Prop%thermP(me%MC_Prop%idxFSlk(i))%p%rhs_SrcS
      endif
    enddo

end subroutine source_from_FS_ports_self_solid



subroutine source_from_SS_src_ports_self_solid(me,dt)
  type(solid_t), intent(inout) :: me
  real(dp), intent(in) :: dt

  integer :: i

  if(me%MC_Prop%SSsrcLink) then
    do i=1,me%MC_Prop%NbCells          
      me%big%bmS(i)=me%big%bmS(i)-dt*me%MC_Prop%thermS(i)%p%DerSrcSdT
      me%big%bvS(i)=me%big%bvS(i)+dt*me%MC_Prop%thermS(i)%p%rhs_SrcS
    enddo   
  endif

end subroutine source_from_SS_src_ports_self_solid



subroutine from_sol_to_temp_solid(me,sim)    
    type(solid_t), intent(inout) :: me
    type(simulation_t), intent(in) :: sim

    integer :: jj

    me%err = 0_dp
    me%err_den = 0_dp

    do jj=1,me%MC_Prop%NbCells
        me%StVar%MCtemp(jj)=me%StVar%MCtemp(jj)+me%big%bvS(jj) ! update of solid temperature by the solution
        
        ! Write to HDF5
        me%hdf%data(jj,1) = me%StVar%MCtemp(jj)

        ! For step management
        if(.not. sim%explicit) then
            me%err = me%err + ( (1.0_dp/6.0_dp) * (sim%dtPrev1+sim%dt) * ( &
                (me%StVar%MCtemp(jj)-me%StVarOld%MCtemp(jj))               / sim%dt - &
                (me%StVarOld%MCtemp(jj)-me%StVarOld2%MCtemp(jj))  * (sim%dt+sim%dtPrev1) / sim%dtPrev1**2 + &
                (me%StVarOld2%MCtemp(jj)-me%StVarOld3%MCtemp(jj)) *  sim%dt/sim%dtPrev1  / sim%dtPrev2 ))**2
            me%err_den = me%err_den + (me%StVarOld%MCtemp(jj))**2
        endif
    enddo
      
end subroutine from_sol_to_temp_solid



subroutine explicit_scheme_for_fluxes_solid(me,dt)
    type(solid_t), intent(inout) :: me
    real(dp),     intent(in)    :: dt

    real(dp) :: VarcentreeS, cpMS
    integer :: i, ii

    if(me%MC_Prop%cond_btw_nodes) then
        do i=1,me%MC_Prop%NbCells
            cpMS=me%mat%heat_capacity(me%StVar%MCtemp(i))
            VarcentreeS=1.0_dp/(me%mat%density*cpMS)
            me%StVar%MCtemp(i)=me%StVar%MCtemp(i)-(dt/me%MC_Prop%dxLoc(i))*VarcentreeS*(me%flxS%sd(i)-me%flxS%sd(i-1))
        enddo
    endif

end subroutine explicit_scheme_for_fluxes_solid



subroutine full_physics_definition_solid(me,sim)
    type(solid_t),      intent(inout) :: me
    type(simulation_t), intent(in)    :: sim
    
    real(dp) :: wwn,ro_M,dx,cpMS,dcpMSdT,cpMSplus,VarcentreeS_i,VarcentreeS_iplus,DerVarcentreeS
    integer :: i,NbCells
    wwn = sim%dt/sim%dtPrev1
    me%big%amS=0.0_dp;me%big%bmS=0.0_dp;me%big%cmS=0.0_dp;me%big%bvS=0.0_dp

    NbCells=me%MC_Prop%NbCells

    if(me%MC_Prop%cond_btw_nodes) then
      do i=1,NbCells
        dx=me%MC_Prop%dxLoc(i)
         
        cpMS=me%mat%heat_capacity(me%StVar%MCtemp(i))
        dcpMSdT=me%mat%heat_capacity_der(me%StVar%MCtemp(i))
        if(i/=NbCells) then
          cpMSplus=me%mat%heat_capacity(me%StVar%MCtemp(i+1))
        endif
      
        ro_M=me%mat%density
        VarcentreeS_i=1.0_dp/(ro_M*cpMS)
        if(i/=NbCells) VarcentreeS_iplus=1.0_dp/(ro_M*cpMSplus)
        DerVarcentreeS=-dcpMSdT/(ro_M*cpMS*cpMS)
      
        if(i/=NbCells) me%big%amS(i)=me%big%amS(i)-(sim%dt/me%MC_Prop%dxLoc(i+1))*VarcentreeS_iplus*me%flxS%sd_DerdTL(i)
        me%big%bmS(i)=me%big%bmS(i)+1.0_dp-sim%dt/dx*(VarcentreeS_i*(me%flxS%sd_DerdTR(i-1)-me%flxS%sd_DerdTL(i))+&
                                                  (me%flxS%sd(i-1)-me%flxS%sd(i))*DerVarcentreeS)
        if(i/=NbCells) me%big%cmS(i)=me%big%cmS(i)+(sim%dt/dx)*VarcentreeS_i*me%flxS%sd_DerdTR(i)
        me%big%bvS(i)=me%big%bvS(i)-(sim%dt/dx)*VarcentreeS_i*(me%flxS%sd(i)-me%flxS%sd(i-1))
      enddo   
    else
      me%big%bmS(:)=1.0_dp
    endif 

    call source_term_definition_solid(me,sim)
    
    do i=1,NbCells
      ! Modification of the diagonal terms 
        me%big%bmS(i)=me%big%bmS(i)-1.0_dp+((1.0_dp+2.0_dp*wwn)/(1.0_dp+wwn))
          
      ! Modification of the vector B by adding time step n-1 contribution
        me%big%bvS(i)=me%big%bvS(i)+((wwn**2)/(1.0_dp+wwn))*(me%StVar%MCtemp(i)-me%StVarOld2%MCtemp(i))
    enddo

end subroutine full_physics_definition_solid



subroutine source_term_definition_solid(me,sim)
  type(solid_t), intent(inout) :: me
  type(simulation_t), intent(in) :: sim

  integer :: i
  real(dp) :: cp, dcpdT, ro_M, Spoint, DerSpoint
  
    if(sim%explicit) then
      me%big%bmS=1.0_dp; me%big%bvS=0.0_dp
    endif

    do i=1,me%MC_Prop%NbCells
      cp    = me%mat%heat_capacity(me%StVar%MCtemp(i))
      dcpdT = me%mat%heat_capacity_der(me%StVar%MCtemp(i))
      ro_M  = me%mat%density

      Spoint    =  me%Q_ext(i) / (me%MC_Prop%volLoc(i)/me%MC_Prop%dxLoc(i) * ro_M * cp)
      DerSpoint = -me%Q_ext(i) * dcpdT / (me%MC_Prop%volLoc(i)/me%MC_Prop%dxLoc(i) * ro_M * cp**2)
      
      me%big%bmS(i) = me%big%bmS(i) - sim%dt * DerSpoint
      me%big%bvS(i) = me%big%bvS(i) + sim%dt * Spoint

      if(me%MC_Prop%FSlink(i)) then
        ! Communication to FS port associated to solid node i
        me%MC_Prop%thermP(me%MC_Prop%idxFSlk(i))%p%TempS=me%StVar%MCtemp(i)
        me%MC_Prop%thermP(me%MC_Prop%idxFSlk(i))%p%lambdS=me%mat%thermal_conductivity(me%StVar%MCtemp(i))
        me%MC_Prop%thermP(me%MC_Prop%idxFSlk(i))%p%DerlambdS=me%mat%thermal_conductivity_der(me%StVar%MCtemp(i))
        me%MC_Prop%thermP(me%MC_Prop%idxFSlk(i))%p%cpS=me%mat%heat_capacity(me%StVar%MCtemp(i))
        me%MC_Prop%thermP(me%MC_Prop%idxFSlk(i))%p%dcpSdT=me%mat%heat_capacity_der(me%StVar%MCtemp(i))        
      endif

      if(me%MC_Prop%SSsrcLink) then
        ! Communication to SS_src port associated to solid node i
        me%MC_Prop%thermS(i)%p%TempS=me%StVar%MCtemp(i)
        me%MC_Prop%thermS(i)%p%lambdS=me%mat%thermal_conductivity(me%StVar%MCtemp(i))
        me%MC_Prop%thermS(i)%p%DerlambdS=me%mat%thermal_conductivity_der(me%StVar%MCtemp(i))
        me%MC_Prop%thermS(i)%p%cpS=me%mat%heat_capacity(me%StVar%MCtemp(i))
        me%MC_Prop%thermS(i)%p%dcpSdT=me%mat%heat_capacity_der(me%StVar%MCtemp(i))        
      endif
    enddo      
end subroutine source_term_definition_solid 

end module cmp_solid_calc_m
