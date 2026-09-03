! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_mesh2D_calc_m
    use cmp_mesh2D_init_m
    use krn_simulation_m
    implicit none

contains


subroutine from_sol_to_temp_mesh2D(me,sim)    
    type(mesh2D_t), intent(inout) :: me
    type(simulation_t), intent(in) :: sim

    integer :: nb_e

    me%err = 0_dp
    me%err_den = 0_dp

    do nb_e=1,me%M2D_Prop%nb_elements
        me%StVar%temp(nb_e)=me%StVar%temp(nb_e)+me%big%bvx(nb_e) ! update of mesh2D temperature by the solution
        
        ! For step management
        if(.not. sim%explicit) then
            me%err = me%err + ( (1.0_dp/6.0_dp) * (sim%dtPrev1+sim%dt) * ( &
                (me%StVar    %temp(nb_e)-me%StVarOld %temp(nb_e))               / sim%dt - &
                (me%StVarOld %temp(nb_e)-me%StVarOld2%temp(nb_e)) * (sim%dt+sim%dtPrev1) / sim%dtPrev1**2 + &
                (me%StVarOld2%temp(nb_e)-me%StVarOld3%temp(nb_e)) *  sim%dt/sim%dtPrev1  / sim%dtPrev2 ))**2
            me%err_den = me%err_den + (me%StVarOld%temp(nb_e))**2
        endif
    enddo
      
end subroutine from_sol_to_temp_mesh2D



subroutine full_physics_definition_mesh2D(me,sim,Qext_Check_case_loc)
    type(mesh2D_t),     intent(inout) :: me
    type(simulation_t), intent(in)    :: sim
    real(dp), intent(out), optional :: Qext_Check_case_loc

    real(dp) :: wwn
    real(dp) :: ro,cp,dcpdT,VarcentS,DerVarcentS,dt_Vi,Q_case_sum
    integer  :: nb_e, NbCells, jj, i
  
      wwn = sim%dt/sim%dtPrev1
      me%big%bmx=0.0_dp;me%big%bvx=0.0_dp;me%big%dmx=0.0_dp

      NbCells=me%M2D_Prop%nb_elements

      do nb_e=1,NbCells
        do i = 1, size(me%label%regionArr)
            if(trim(me%M2D_Prop%elem(nb_e)%PhysE)==trim(me%label%regionArr(i)%id)) then
                cp=me%label%regionArr(i)%mat%heat_capacity(me%StVar%temp(nb_e))
                dcpdT=me%label%regionArr(i)%mat%heat_capacity_der(me%StVar%temp(nb_e))
                ro=me%label%regionArr(i)%mat%density
                Q_case_sum = me%label%regionArr(i)%value%v0d()/(me%label%regionArr(i)%SurfaceRegion*me%extrusion_length)
                exit
            endif
        enddo

        VarcentS=1.0_dp/(ro*cp)
        DerVarcentS=-dcpdT/(ro*cp*cp)

        dt_Vi=sim%dt/me%M2D_Prop%elem(nb_e)%surface

        me%big%bmx(nb_e)=1.0_dp+dt_Vi*DerVarcentS*me%flxS%SomNormQT(nb_e)+dt_Vi*VarcentS*me%flxS%DerSomNormQT(nb_e)
           
        me%big%bvx(nb_e)=-dt_Vi*VarcentS*me%flxS%SomNormQT(nb_e)
        
        do jj=1,3
          me%big%dmx(nb_e,jj)=dt_Vi*VarcentS*me%flxS%DerNormQTOth2D(nb_e,jj)
        enddo

        ! Source terms
        me%big%bmx(nb_e) = me%big%bmx(nb_e) - sim%dt * DerVarcentS * Q_case_sum
        me%big%bvx(nb_e) = me%big%bvx(nb_e) + sim%dt * VarcentS    * Q_case_sum

        ! if(present(Qext_Check_case_loc) .and. trim(me%M2D_Prop%elem(nb_e)%PhysE)=='SS_Outer') then
        !   Qext_Check_case_loc = Q_case_sum * sim%dt
        ! endif
      enddo
    
      if(.not.sim%explicit) then
        do nb_e=1,NbCells
          ! Modification of the diagonal terms 
            me%big%bmx(nb_e)=me%big%bmx(nb_e)-1.0_dp+((1.0_dp+2.0_dp*wwn)/(1.0_dp+wwn))

          ! Modification of the vector B by adding time step n-1 contribution
            me%big%bvx(nb_e)=me%big%bvx(nb_e)+((wwn**2)/(1.0_dp+wwn))*(me%StVar%temp(nb_e)-me%StVarOld2%temp(nb_e))
        enddo
      endif

end subroutine full_physics_definition_mesh2D



subroutine mesh2Ds_SS_src_port_comm(me)
  type(mesh2D_t), intent(inout) :: me

  integer :: ii,nb_e,idx_prev,elem, i
  real(dp) :: cp, ro, dcpdT, lambdS, DerlambdS

  if(me%M2D_Prop%Nb_labelMC>0) then
    idx_prev=0
    do ii=1,me%M2D_Prop%Nb_labelMC
      do nb_e=1,me%label%MC_Arr(ii)%NbElem
        elem=me%label%MC_Arr(ii)%array(nb_e)
        me%M2D_Prop%thermS(idx_prev+nb_e)%p%TempS=me%StVar%temp(elem)
        do i = 1, size(me%label%regionArr)
            if(trim(me%M2D_Prop%elem(elem)%PhysE)==trim(me%label%regionArr(i)%id)) then
                cp=me%label%regionArr(i)%mat%heat_capacity(me%StVar%temp(elem))
                dcpdT=me%label%regionArr(i)%mat%heat_capacity_der(me%StVar%temp(elem))
                ro=me%label%regionArr(i)%mat%density
                lambdS = me%label%regionArr(i)%mat%thermal_conductivity(me%StVar%temp(elem))
                DerlambdS = me%label%regionArr(i)%mat%thermal_conductivity_der(me%StVar%temp(elem))
                exit
            endif
        enddo
        me%M2D_Prop%thermS(idx_prev+nb_e)%p%cpS        = cp
        me%M2D_Prop%thermS(idx_prev+nb_e)%p%dcpSdT     = dcpdT
        me%M2D_Prop%thermS(idx_prev+nb_e)%p%lambdS     = lambdS
        me%M2D_Prop%thermS(idx_prev+nb_e)%p%DerlambdS  = DerlambdS        
        me%M2D_Prop%thermS(idx_prev+nb_e)%p%rhoS       = ro
        me%M2D_Prop%thermS(idx_prev+nb_e)%p%Dist4Grad  = me%M2D_Prop%elem(elem)%Delta(3) ! 3 forced here because 2 direct neighbors expected per 2D cell - Check is done during preprocessing
        me%M2D_Prop%thermS(idx_prev+nb_e)%p%LengthCont = me%M2D_Prop%elem(elem)%face(3)  ! 3 forced here because 2 direct neighbors expected per 2D cell - Check is done during preprocessing
        me%M2D_Prop%thermS(idx_prev+nb_e)%p%SurfS      = me%M2D_Prop%elem(elem)%surface
      enddo
      idx_prev=idx_prev+me%label%MC_Arr(ii)%NbElem
    enddo
  endif
  
end subroutine mesh2Ds_SS_src_port_comm



subroutine mesh2Ds_FS_port_comm(me)
  type(mesh2D_t), intent(inout) :: me

  integer :: ii,nb_e,idx_prev,elem, i
  real(dp) :: cp, ro, dcpdT, lambdS, DerlambdS

  if(me%M2D_Prop%Nb_labelCh>0) then
    idx_prev=0
    do ii=1,me%M2D_Prop%Nb_labelCh
      do nb_e=1,me%label%Ch_Arr(ii)%NbElem
        elem=me%label%Ch_Arr(ii)%array(nb_e)
        me%M2D_Prop%thermP(idx_prev+nb_e)%p%TempS=me%StVar%temp(elem)
        do i = 1, size(me%label%regionArr)
            if(trim(me%M2D_Prop%elem(elem)%PhysE)==trim(me%label%regionArr(i)%id)) then
                cp=me%label%regionArr(i)%mat%heat_capacity(me%StVar%temp(elem))
                dcpdT=me%label%regionArr(i)%mat%heat_capacity_der(me%StVar%temp(elem))
                ro=me%label%regionArr(i)%mat%density
                lambdS = me%label%regionArr(i)%mat%thermal_conductivity(me%StVar%temp(elem))
                DerlambdS = me%label%regionArr(i)%mat%thermal_conductivity_der(me%StVar%temp(elem))
                exit
            endif
        enddo
        me%M2D_Prop%thermP(idx_prev+nb_e)%p%cpS        = cp
        me%M2D_Prop%thermP(idx_prev+nb_e)%p%dcpSdT     = dcpdT
        me%M2D_Prop%thermP(idx_prev+nb_e)%p%lambdS     = lambdS
        me%M2D_Prop%thermP(idx_prev+nb_e)%p%DerlambdS  = DerlambdS
        me%M2D_Prop%thermP(idx_prev+nb_e)%p%rhoS       = ro
        me%M2D_Prop%thermP(idx_prev+nb_e)%p%Dist4Grad  = me%M2D_Prop%elem(elem)%Delta(3) ! 3 forced here because 2 direct neighbors expected per 2D cell - Check is done during preprocessing
        me%M2D_Prop%thermP(idx_prev+nb_e)%p%LengthCont = me%M2D_Prop%elem(elem)%face(3)  ! 3 forced here because 2 direct neighbors expected per 2D cell - Check is done during preprocessing
        me%M2D_Prop%thermP(idx_prev+nb_e)%p%SurfS      = me%M2D_Prop%elem(elem)%surface
      enddo
      idx_prev=idx_prev+me%label%Ch_Arr(ii)%NbElem
    enddo
  endif
  
end subroutine mesh2Ds_FS_port_comm



subroutine source_from_SS_src_ports_self_mesh2D(me,dt)
  type(mesh2D_t), intent(inout) :: me
  real(dp), intent(in) :: dt

  real(dp) :: VarcentS,DerVarcentS,dt_Vi
  integer :: ii,idx_prev,nb_e,elem

  if(me%M2D_Prop%Nb_labelMC>0) then
    idx_prev=0
    do ii=1,me%M2D_Prop%Nb_labelMC
      do nb_e=1,me%label%MC_Arr(ii)%NbElem
        elem=me%label%MC_Arr(ii)%array(nb_e)

        dt_Vi=dt/me%M2D_Prop%elem(elem)%surface
        VarcentS=1.0_dp/(me%M2D_Prop%thermS(idx_prev+nb_e)%p%rhoS*me%M2D_Prop%thermS(idx_prev+nb_e)%p%cpS)
        DerVarcentS=-me%M2D_Prop%thermS(idx_prev+nb_e)%p%dcpSdT/&
                    (me%M2D_Prop%thermS(idx_prev+nb_e)%p%rhoS*me%M2D_Prop%thermS(idx_prev+nb_e)%p%cpS*&
                    me%M2D_Prop%thermS(idx_prev+nb_e)%p%cpS)

        me%big%bmx(elem)=me%big%bmx(elem)+dt_Vi*DerVarcentS*me%M2D_Prop%thermS(idx_prev+nb_e)%p%Flx+ &
                         dt_Vi*VarcentS*me%M2D_Prop%thermS(idx_prev+nb_e)%p%derFlx_derCon

        me%big%bvx(elem)=me%big%bvx(elem)-dt_Vi*VarcentS*me%M2D_Prop%thermS(idx_prev+nb_e)%p%Flx
      enddo
      idx_prev=idx_prev+me%label%MC_Arr(ii)%NbElem
    enddo
  endif

end subroutine source_from_SS_src_ports_self_mesh2D



subroutine source_from_FS_ports_self_mesh2D(me,dt)
  type(mesh2D_t), intent(inout) :: me
  real(dp), intent(in) :: dt

  real(dp) :: VarcentS,DerVarcentS,dt_Vi
  integer :: ii,idx_prev,nb_e,elem

  if(me%M2D_Prop%Nb_labelCh>0) then
    idx_prev=0
    do ii=1,me%M2D_Prop%Nb_labelCh
      do nb_e=1,me%label%Ch_Arr(ii)%NbElem
        elem=me%label%Ch_Arr(ii)%array(nb_e)

        dt_Vi=dt/me%M2D_Prop%elem(elem)%surface
        VarcentS=1.0_dp/(me%M2D_Prop%thermP(idx_prev+nb_e)%p%rhoS*me%M2D_Prop%thermP(idx_prev+nb_e)%p%cpS)
        DerVarcentS=-me%M2D_Prop%thermP(idx_prev+nb_e)%p%dcpSdT/&
                    (me%M2D_Prop%thermP(idx_prev+nb_e)%p%rhoS*me%M2D_Prop%thermP(idx_prev+nb_e)%p%cpS*&
                    me%M2D_Prop%thermP(idx_prev+nb_e)%p%cpS)

        me%big%bmx(elem)=me%big%bmx(elem)+dt_Vi*DerVarcentS*me%M2D_Prop%thermP(idx_prev+nb_e)%p%Flx+ &
                         dt_Vi*VarcentS*me%M2D_Prop%thermP(idx_prev+nb_e)%p%derFlx_derCon

        me%big%bvx(elem)=me%big%bvx(elem)-dt_Vi*VarcentS*me%M2D_Prop%thermP(idx_prev+nb_e)%p%Flx
      enddo
      idx_prev=idx_prev+me%label%Ch_Arr(ii)%NbElem
    enddo
  endif

end subroutine source_from_FS_ports_self_mesh2D

end module cmp_mesh2D_calc_m
