! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_strand_calc_m
    use cmp_strand_init_m
    use cmp_strand_source_terms_m
    implicit none

contains


subroutine compute_effective_field(me)
  type(strand_t), intent(inout) :: me

  real(kind=dp) :: dBdx,dBdy,T,Jop,E0,n,Area,Jc,B0,eps
  real(kind=dp) :: fIntegral,dTheta,dR,Theta,R,xR,yR,f,A,E,JcEq
  real(kind=dp) :: tolerance,Bmin,Bmax,BError,BEff
  real(kind=dp) :: r0,emax,St,Tc0,Bc
  integer :: i,nTheta,nR,iTheta,iR,Iterat
  logical :: converged

  E0   = me%mat_supc%E0
  n    = me%mat_supc%nPow
  Area = me%ct(SUPC)%area
  eps=1.0e-8_dp

  do i=1,me%SC_Prop%NbCells
    T=me%StVar%SCtemp(i) ! Solid strand temperature
    Jop=abs(me%scen%ElCur)/Area !  Current density
    dBdx=me%scen%dBfield(i)
    dBdy=0.0_dp
    B0=me%scen%Bfield(i)

    if(Jop<=0.0_dp) cycle

    St  = me%mat_supc%strain(me%scen%Bfield(i),me%scen%ElCur)
    Tc0 = me%mat_supc%critical_temperature(0.0_dp,St)
    Bc  = me%mat_supc%critical_field(T,St)
    Jc  = me%mat_supc%critical_current_density(T,me%scen%Bfield(i),St,Tc0,Bc)
    if(Jc<=0.0_dp) cycle

    ! Integrate the electric field
    fIntegral = 0.0_dp
    nTheta    = 5
    nR        = 5
    dTheta    = 2.0_dp*Pi_value/dfloat(nTheta)
    dR        = (me%SC_Prop%outer_rad-me%SC_Prop%inner_rad)/dfloat(nR)
    r0=1.0e7_dp
	  emax=r0*Jop
    do iTheta = 1,nTheta
      ! compute angle
      Theta  = 2.0_dp*Pi_value*dfloat(iTheta-1)/dfloat(nTheta-1)
      do iR  = 1,nR
        ! compute radius
        R = dfloat(iR-1)/dfloat(nR-1)
        R = me%SC_Prop%inner_rad*(1.0_dp-R) + me%SC_Prop%outer_rad*R
        ! compute x and y locations in the cross section
        xR = R*cos(Theta)
        yR = R*sin(Theta)
        ! compute B at the location
        me%scen%Bfield(i) = B0 + dBdx*xR + dBdy*yR
        ! compute Jc
        St  = me%mat_supc%strain(me%scen%Bfield(i),me%scen%ElCur)
        Tc0 = me%mat_supc%critical_temperature(0.0_dp,St)
        Bc  = me%mat_supc%critical_field(T,St)
        Jc  = me%mat_supc%critical_current_density(T,me%scen%Bfield(i),St,Tc0,Bc)        
        ! compute the local electric field checking for normal state
        if(Jc>0.0_dp) then
           f = E0*(Jop/Jc)**n
        else
           f=emax
          ! print*, 'critical current calculation failed'
          ! print*, iTheta,iR
          ! print*, Jc,i,T,me%scen%Bfield(i)
          ! print*, 'failed in compute_effective_field'
          ! read(*,*)
        endif
        ! integrate adding contributions of single area differentials
        fIntegral = fIntegral+f*dR*R*dTheta
      enddo
    enddo
    ! average electric field from normalised integral
    A = Pi_value*(me%SC_Prop%outer_rad**2-me%SC_Prop%inner_rad**2)
    E = fIntegral/A

    ! trap limiting case of zero electric field
    if(E<=eps*E0) then
      me%scen%Bfield(i)=B0
      cycle
    endif

    ! find equivalent Jc, corresponding to the computed electric field
    JcEq = Jop/(E/E0)**(1.0/n)

    ! find iteratively equivalent magnetic field
    converged = .false.
    Tolerance = eps
    Bmin      = B0
    Bmax      = B0 + sqrt(dBdx**2+dBdy**2)*me%SC_Prop%outer_rad
    Iterat = 0
    do while(.not.converged)
      Iterat = Iterat + 1
      BEff = 0.5*(BMin+BMax)
      if(Iterat>50) then
        call set_error('convergence failure in compute_effective_field')
        return
      endif
      me%scen%Bfield(i)=BEff
      St  = me%mat_supc%strain(me%scen%Bfield(i),me%scen%ElCur)
      Tc0 = me%mat_supc%critical_temperature(0.0_dp,St)
      Bc  = me%mat_supc%critical_field(T,St)
      Jc  = me%mat_supc%critical_current_density(T,me%scen%Bfield(i),St,Tc0,Bc)
      if(Jc>JcEq) then
         Bmin=BEff
      elseif(Jc<JcEq) then
         Bmax=BEff
      else
         Bmin=BEff
         Bmax=BEff
      endif
      BError    = (Bmax-Bmin)/BEff
      converged = BError<=Tolerance
    enddo
    me%scen%Bfield(i)=BEff
  enddo
end subroutine compute_effective_field



subroutine scenario_update(me,sim)
  type(strand_t), intent(inout) :: me
  type(simulation_t), intent(in) :: sim

  me%scen%Q_ext   = me%scen%Q_ext_load%v1d()
  me%scen%Bfield  = me%scen%Bfield_load%v1d()
  me%scen%dBfield = me%scen%dBfield_load%v1d()
  me%scen%ElCur   = me%scen%ElCur_load%v0d()

  call compute_effective_field(me)
  if (sim_error > 0) return
end subroutine scenario_update



subroutine source_from_FS_ports_self_strand(me,dt)
    type(strand_t), intent(inout) :: me
    real(dp), intent(in) :: dt

    integer :: i

    if(me%SC_Prop%FSlink) then
      do i=1,me%SC_Prop%NbCells
          me%big%bmS(i)=me%big%bmS(i)-dt*me%SC_Prop%thermP(i)%p%DerSrcSdT
          me%big%bvS(i)=me%big%bvS(i)+dt*me%SC_Prop%thermP(i)%p%rhs_SrcS
      enddo
    endif

end subroutine source_from_FS_ports_self_strand



subroutine strands_SS_flux_port_comm(me)
  type(strand_t), intent(inout) :: me
  integer :: ii
  real(dp) :: lambdS, DerlambdS, cpMS
  real(dp) :: B,St,Tc,Tc0,Bc,Bc0,Jc0,Jop,TcS

  lambdS=0.0_dp
  DerlambdS=0.0_dp
  cpMS=0.0_dp

  lambdS=lambdS+me%ct(STAB)%area*me%mat_stab%thermal_conductivity(me%StVar%SCtemp(1),me%scen%Bfield(1))      
  lambdS=lambdS+me%ct(SUPC)%area*me%mat_supc%thermal_conductivity(me%StVar%SCtemp(1))      

  DerlambdS=DerlambdS+me%ct(STAB)%area*me%mat_stab%thermal_conductivity_der(me%StVar%SCtemp(1),me%scen%Bfield(1))
  DerlambdS=DerlambdS+me%ct(SUPC)%area*me%mat_supc%thermal_conductivity_der(me%StVar%SCtemp(1))

  cpMS=cpMS+me%ct(STAB)%area*me%mat_stab%density*me%mat_stab%heat_capacity(me%StVar%SCtemp(1))

  B   = me%scen%Bfield(1)
  St  = me%mat_supc%strain(B,me%scen%ElCur)
  TC  = me%mat_supc%critical_temperature(B,St)
  Tc0 = me%mat_supc%critical_temperature(0.0_dp,St)
  Bc  = me%mat_supc%critical_field(me%StVar%SCtemp(1),St)
  Bc0 = me%mat_supc%critical_field(0.0_dp,St)
  Jc0 = me%mat_supc%critical_current_density(0.0_dp,B,St,Tc0,Bc)
  Jop = abs(me%scen%ElCur)/me%ct(SUPC)%area
  TcS=me%mat_supc%current_sharing_temperature(B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc)
  cpMS=cpMS+me%ct(SUPC)%area*me%mat_supc%density*me%mat_supc%heat_capacity(me%StVar%SCtemp(1),B,TC,TcS,Tc0)

  lambdS=lambdS/sum(me%ct(:)%area)
  DerlambdS=DerlambdS/sum(me%ct(:)%area)
  cpMS=cpMS/(me%ct(STAB)%area*me%mat_stab%density+me%ct(SUPC)%area*me%mat_supc%density)

  me%SC_Prop%in%lambdS    = lambdS
  me%SC_Prop%in%DerlambdS = DerlambdS
  me%SC_Prop%in%cpMS      = cpMS
  me%SC_Prop%in%TempS     = me%StVar%SCtemp(1)


  lambdS=0.0_dp
  DerlambdS=0.0_dp
  cpMS=0.0_dp

  lambdS=lambdS+me%ct(STAB)%area* &
         me%mat_stab%thermal_conductivity(me%StVar%SCtemp(me%SC_Prop%NbCells),me%scen%Bfield(me%SC_Prop%NbCells))
  lambdS=lambdS+me%ct(SUPC)%area* &
         me%mat_supc%thermal_conductivity(me%StVar%SCtemp(me%SC_Prop%NbCells))             

  DerlambdS=DerlambdS+me%ct(STAB)%area*&
            me%mat_stab%thermal_conductivity_der(me%StVar%SCtemp(me%SC_Prop%NbCells),me%scen%Bfield(me%SC_Prop%NbCells))
  DerlambdS=DerlambdS+me%ct(SUPC)%area*&
            me%mat_supc%thermal_conductivity_der(me%StVar%SCtemp(me%SC_Prop%NbCells))                

  cpMS=cpMS+me%ct(STAB)%area*me%mat_stab%density*me%mat_stab%heat_capacity(me%StVar%SCtemp(me%SC_Prop%NbCells))

  B   = me%scen%Bfield(me%SC_Prop%NbCells)
  St  = me%mat_supc%strain(B,me%scen%ElCur)
  TC  = me%mat_supc%critical_temperature(B,St)
  Tc0 = me%mat_supc%critical_temperature(0.0_dp,St)
  Bc  = me%mat_supc%critical_field(me%StVar%SCtemp(me%SC_Prop%NbCells),St)
  Bc0 = me%mat_supc%critical_field(0.0_dp,St)
  Jc0 = me%mat_supc%critical_current_density(0.0_dp,B,St,Tc0,Bc)
  Jop = abs(me%scen%ElCur)/me%ct(SUPC)%area
  TcS=me%mat_supc%current_sharing_temperature(B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc)      
  cpMS=cpMS+me%ct(SUPC)%area*me%mat_supc%density*me%mat_supc%heat_capacity(me%StVar%SCtemp(me%SC_Prop%NbCells),B,TC,TcS,Tc0)

  lambdS=lambdS/sum(me%ct(:)%area)
  DerlambdS=DerlambdS/sum(me%ct(:)%area)
  cpMS=cpMS/(me%ct(STAB)%area*me%mat_stab%density+me%ct(SUPC)%area*me%mat_supc%density)

  me%SC_Prop%out%lambdS    = lambdS
  me%SC_Prop%out%DerlambdS = DerlambdS
  me%SC_Prop%out%cpMS      = cpMS
  me%SC_Prop%out%TempS     = me%StVar%SCtemp(me%SC_Prop%NbCells)

end subroutine strands_SS_flux_port_comm



subroutine flux_from_SS_flux_ports_self(me)
  type(strand_t), intent(inout) :: me

  ! in component
  me%flxS%wire(0)=me%SC_Prop%in%Flx
  me%flxS%wire_DerdTR(0)=me%SC_Prop%in%derFlx_derCon

  ! out component
  me%flxS%wire(me%SC_Prop%NbCells)=me%SC_Prop%out%Flx
  me%flxS%wire_DerdTL(me%SC_Prop%NbCells)=me%SC_Prop%out%derFlx_derCon

end subroutine flux_from_SS_flux_ports_self



subroutine from_sol_to_temp_strand(me,sim)
    type(strand_t), intent(inout) :: me
    type(simulation_t), intent(in) :: sim

    integer :: jj
    real(dp) :: B,St,Tc,Tc0,Bc,Bc0,Jc0,Jop,TcS

    me%err = 0_dp
    me%err_den = 0_dp

    do jj=1,me%SC_Prop%NbCells
        me%StVar%SCtemp(jj)=me%StVar%SCtemp(jj)+me%big%bvS(jj) ! update of strand temperature by the solution

        B   = me%scen%Bfield(jj)
        St  = me%mat_supc%strain(B,me%scen%ElCur)
        TC  = me%mat_supc%critical_temperature(B,St)
        Tc0 = me%mat_supc%critical_temperature(0.0_dp,St)
        Bc  = me%mat_supc%critical_field(me%StVar%SCtemp(jj),St)
        Bc0 = me%mat_supc%critical_field(0.0_dp,St)
        Jc0 = me%mat_supc%critical_current_density(0.0_dp,B,St,Tc0,Bc)
        Jop = abs(me%scen%ElCur)/me%ct(SUPC)%area
        TcS=me%mat_supc%current_sharing_temperature(B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc)        

        ! Write to HDF5
        me%hdf%data(jj,1) = me%StVar%SCtemp(jj)
        me%hdf%data(jj,2) = me%scen%Bfield(jj)
        me%hdf%data(jj,3) = St
        me%hdf%data(jj,4) = TcS
        me%hdf%data(jj,5) = TcS-me%StVar%SCtemp(jj)

        ! For step management
        if(.not. sim%explicit) then
            me%err = me%err + ( (1.0_dp/6.0_dp) * (sim%dtPrev1+sim%dt) * ( &
                (me%StVar    %SCtemp(jj)-me%StVarOld %SCtemp(jj))               / sim%dt - &
                (me%StVarOld %SCtemp(jj)-me%StVarOld2%SCtemp(jj)) * (sim%dt+sim%dtPrev1) / sim%dtPrev1**2 + &
                (me%StVarOld2%SCtemp(jj)-me%StVarOld3%SCtemp(jj)) *  sim%dt/sim%dtPrev1  / sim%dtPrev2 ))**2
            me%err_den = me%err_den + (me%StVarOld%SCtemp(jj))**2
        endif
    enddo

end subroutine from_sol_to_temp_strand



subroutine explicit_scheme_for_fluxes_strand(me,dt)
    type(strand_t), intent(inout) :: me
    real(dp),     intent(in)    :: dt

    real(dp) :: VarcentreeS, cpMS
    integer :: i, ii
    real(dp) :: B,St,Tc,Tc0,Bc,Bc0,Jc0,Jop,TcS

    do i=1,me%SC_Prop%NbCells
        cpMS=0.0_dp

        cpMS=cpMS+me%ct(STAB)%area*me%mat_stab%density*me%mat_stab%heat_capacity(me%StVar%SCtemp(i))
      
        B   = me%scen%Bfield(i)
        St  = me%mat_supc%strain(B,me%scen%ElCur)
        TC  = me%mat_supc%critical_temperature(B,St)
        Tc0 = me%mat_supc%critical_temperature(0.0_dp,St)
        Bc  = me%mat_supc%critical_field(me%StVar%SCtemp(i),St)
        Bc0 = me%mat_supc%critical_field(0.0_dp,St)
        Jc0 = me%mat_supc%critical_current_density(0.0_dp,B,St,Tc0,Bc)
        Jop = abs(me%scen%ElCur)/me%ct(SUPC)%area
        TcS = me%mat_supc%current_sharing_temperature(B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc)
        cpMS=cpMS+me%ct(SUPC)%area*me%mat_supc%density*me%mat_supc%heat_capacity(me%StVar%SCtemp(i),B,TC,TcS,Tc0)

        cpMS=cpMS/(me%ct(STAB)%area*me%mat_stab%density+me%ct(SUPC)%area*me%mat_supc%density)

        VarcentreeS=1.0_dp/(me%ro_M*cpMS)
        me%StVar%SCtemp(i)=me%StVar%SCtemp(i)-(dt/me%SC_Prop%dxLoc(i))*VarcentreeS*(me%flxS%wire(i)-me%flxS%wire(i-1))
    enddo

end subroutine explicit_scheme_for_fluxes_strand



subroutine full_physics_definition_strand(me,sim,Qext_Check)
    type(strand_t),     intent(inout) :: me
    type(simulation_t), intent(in)    :: sim
    real(dp), intent(out), optional :: Qext_Check

    real(dp) :: wwn
    real(dp) :: dx,cpMS,dcpMSdT,cpMSplus,VarcentreeS_i,VarcentreeS_iplus,DerVarcentreeS
    real(dp) :: B,St,Tc,Tc0,Bc,Bc0,Jc0,Jop,TcS
    integer :: i,ii,NbCells

      wwn = sim%dt/sim%dtPrev1
      me%big%amS=0.0_dp;me%big%bmS=0.0_dp;me%big%cmS=0.0_dp;me%big%bvS=0.0_dp

      NbCells=me%SC_Prop%NbCells

      do i=1,NbCells

        dx=me%SC_Prop%dxLoc(i)

        cpMS=0.0_dp
        dcpMSdT=0.0_dp
        cpMSplus=0.0_dp
        
        cpMS=cpMS+me%ct(STAB)%area*me%mat_stab%density*me%mat_stab%heat_capacity(me%StVar%SCtemp(i))
        dcpMSdT=dcpMSdT+me%ct(STAB)%area*me%mat_stab%density*me%mat_stab%heat_capacity_der(me%StVar%SCtemp(i))

        B   = me%scen%Bfield(i)
        St  = me%mat_supc%strain(B,me%scen%ElCur)
        TC  = me%mat_supc%critical_temperature(B,St)
        Tc0 = me%mat_supc%critical_temperature(0.0_dp,St)
        Bc  = me%mat_supc%critical_field(me%StVar%SCtemp(i),St)
        Bc0 = me%mat_supc%critical_field(0.0_dp,St)
        Jc0 = me%mat_supc%critical_current_density(0.0_dp,B,St,Tc0,Bc)
        Jop = abs(me%scen%ElCur)/me%ct(SUPC)%area
        TcS=me%mat_supc%current_sharing_temperature(B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc)      
        cpMS=cpMS+me%ct(SUPC)%area*me%mat_supc%density*me%mat_supc%heat_capacity(me%StVar%SCtemp(i),B,TC,TcS,Tc0)
        dcpMSdT=dcpMSdT+me%ct(SUPC)%area*me%mat_supc%density*me%mat_supc%heat_capacity_der(me%StVar%SCtemp(i),B,TC,TcS,Tc0)

        if(i/=NbCells) then
          cpMSplus=cpMSplus+me%ct(STAB)%area*me%mat_stab%density*me%mat_stab%heat_capacity(me%StVar%SCtemp(i+1))
  
          B   = me%scen%Bfield(i+1)
          St  = me%mat_supc%strain(B,me%scen%ElCur)
          TC  = me%mat_supc%critical_temperature(B,St)
          Tc0 = me%mat_supc%critical_temperature(0.0_dp,St)
          Bc  = me%mat_supc%critical_field(me%StVar%SCtemp(i+1),St)
          Bc0 = me%mat_supc%critical_field(0.0_dp,St)
          Jc0 = me%mat_supc%critical_current_density(0.0_dp,B,St,Tc0,Bc)
          Jop = abs(me%scen%ElCur)/me%ct(SUPC)%area
          TcS=me%mat_supc%current_sharing_temperature(B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc)      
          cpMSplus=cpMSplus+me%ct(SUPC)%area*me%mat_supc%density*me%mat_supc%heat_capacity(me%StVar%SCtemp(i+1),B,TC,TcS,Tc0)
        endif

        cpMS=cpMS/(me%ct(STAB)%area*me%mat_stab%density+me%ct(SUPC)%area*me%mat_supc%density)
        dcpMSdT=dcpMSdT/(me%ct(STAB)%area*me%mat_stab%density+me%ct(SUPC)%area*me%mat_supc%density)
        cpMSplus=cpMSplus/(me%ct(STAB)%area*me%mat_stab%density+me%ct(SUPC)%area*me%mat_supc%density)

        VarcentreeS_i=1.0_dp/(me%ro_M*cpMS)
        if(i/=NbCells) VarcentreeS_iplus=1.0_dp/(me%ro_M*cpMSplus)
        DerVarcentreeS=-dcpMSdT/(me%ro_M*cpMS*cpMS)

        if(i/=NbCells) me%big%amS(i)=me%big%amS(i)-(sim%dt/me%SC_Prop%dxLoc(i+1))*VarcentreeS_iplus*me%flxS%wire_DerdTL(i)
        me%big%bmS(i)=me%big%bmS(i)+1.0_dp-sim%dt/dx*(VarcentreeS_i*(me%flxS%wire_DerdTR(i-1)-me%flxS%wire_DerdTL(i))+&
                                                  (me%flxS%wire(i-1)-me%flxS%wire(i))*DerVarcentreeS)
        if(i/=NbCells) me%big%cmS(i)=me%big%cmS(i)+(sim%dt/dx)*VarcentreeS_i*me%flxS%wire_DerdTR(i)
        me%big%bvS(i)=me%big%bvS(i)-(sim%dt/dx)*VarcentreeS_i*(me%flxS%wire(i)-me%flxS%wire(i-1))

      enddo

      if(present(Qext_Check)) then
        call source_term_definition_strand(me,sim,Qext_Check)
      else
        call source_term_definition_strand(me,sim)
      endif

      do i=1,NbCells
        ! Modification of the diagonal terms
          me%big%bmS(i)=me%big%bmS(i)-1.0_dp+((1.0_dp+2.0_dp*wwn)/(1.0_dp+wwn))

        ! Modification of the vector B by adding time step n-1 contribution
          me%big%bvS(i)=me%big%bvS(i)+((wwn**2)/(1.0_dp+wwn))*(me%StVar%SCtemp(i)-me%StVarOld2%SCtemp(i))
      enddo

end subroutine full_physics_definition_strand



subroutine source_term_definition_strand(me,sim,Qext_Check)
  type(strand_t), intent(inout) :: me
  type(simulation_t), intent(in) :: sim
  real(dp), intent(out), optional :: Qext_Check

  real(dp) :: Spoint, DerSpoint, SJoule, DerSJoule, Joule
  integer :: i

      if(sim%explicit) then
        me%big%bmS=1.0_dp; me%big%bvS=0.0_dp
      endif

      do i=1,me%SC_Prop%NbCells
        ! Adding of source term contributions (external heating + Joule effect)
        call solid_source_terms(me,i,me%scen%Q_ext(i),Spoint,DerSpoint,SJoule,DerSJoule,Joule,sim)
        if (sim_error > 0) return

        if(present(Qext_Check)) Qext_Check=Qext_Check+me%scen%Q_ext(i)*me%SC_Prop%dxLoc(i)*sim%dt

        me%big%bmS(i)=me%big%bmS(i)-sim%dt*(DerSpoint+DerSJoule)
        me%big%bvS(i)=me%big%bvS(i)+sim%dt*(Spoint+SJoule)
        me%hdf%data(i,6) = Joule
      enddo

end subroutine source_term_definition_strand


end module cmp_strand_calc_m
