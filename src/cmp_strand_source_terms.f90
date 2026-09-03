! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_strand_source_terms_m 
    use cmp_strand_init_m
    implicit none
 
contains


subroutine supercond_current(Curr,CritCurr,n,V0,Area,Cond,ratio,derImpFunctdRatio)

    real(dp), intent(in) :: Curr,CritCurr,V0,Area,Cond
    integer, intent(in) :: n
    real(dp), intent(out) :: ratio,derImpFunctdRatio
    real(dp) :: aa,bb,funct
    integer :: NbreMaxIte,ite

    aa=CritCurr/(Area*Cond*V0)
    bb=Curr/(Area*Cond*V0)

    ratio=min(bb/aa,bb**(1.0_dp/float(n)))
    NbreMaxIte=200
    ! Newton method
    do ite=1,NbreMaxIte
      funct=ratio**n+aa*ratio-bb
      derImpFunctdRatio=n*(ratio**(n-1))+aa
      if(abs(funct)>1.0e-8_dp) then
        ratio=ratio-funct/derImpFunctdRatio
      else
        exit
      endif
      if(ite==NbreMaxIte) then
        call set_error('convergence failure in supercond_current')
        return
      endif
    enddo

end subroutine supercond_current



subroutine solid_source_terms(me,i,Qext,Spoint,DerSpoint,SJoule,DerSJoule,Joule,sim)
    
    type(strand_t), intent(inout) :: me
    integer, intent(in) :: i
    real(dp), intent(in) :: Qext
    real(dp), intent(out) :: Spoint, DerSpoint
    real(dp), intent(out) :: SJoule, DerSJoule, Joule
    type(simulation_t), intent(in) :: sim

    integer :: ii,idxSC,idxST
    real(dp) :: AreaMS,rhoMS,cpMS,dcpMSdT
    real(dp) :: CritCurr,ConductMet,rapp
    real(dp) :: Derivee,DerImpFuntdRapp
    real(dp) :: DerCritCurr,DerConductMet
    real(dp) :: denum,DRappdCritCurr,DRappdConductMet,DRappdT
    real(dp) :: B,TC,TcS,Tc0,St,Bc,Bc0,Jc0,Jop

      rhoMS=me%ro_M
      AreaMS=sum(me%ct(:)%area)

      cpMS=0.0_dp
      dcpMSdT=0.0_dp

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
      TcS = me%mat_supc%current_sharing_temperature(B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc)
      
      cpMS=cpMS+me%ct(SUPC)%area*me%mat_supc%density*me%mat_supc%heat_capacity(me%StVar%SCtemp(i),B,TC,TcS,Tc0)      
      dcpMSdT=dcpMSdT+me%ct(SUPC)%area*me%mat_supc%density*me%mat_supc%heat_capacity_der(me%StVar%SCtemp(i),B,TC,TcS,Tc0)

      cpMS=cpMS/(me%ct(STAB)%area*me%mat_stab%density+me%ct(SUPC)%area*me%mat_supc%density)
      dcpMSdT=dcpMSdT/(me%ct(STAB)%area*me%mat_stab%density+me%ct(SUPC)%area*me%mat_supc%density)

    !--------------------------

    ! External heating
      Spoint=Qext/(AreaMS*rhoMS*cpMS)
      DerSpoint=-Qext*dcpMSdT/(AreaMS*rhoMS*(cpMS**2))

    !--------------------------

      ! Joule effect      
      do ii = 1, me%NbConst  
        if(me%ct(ii)%type=='superconductor') then 
          CritCurr=me%mat_supc%critical_current_density(me%StVar%SCtemp(i),B,St,Tc0,Bc)*me%ct(ii)%area
          idxSC=ii
          if(me%mat_supc%nPow>250) then
            call set_error('infinite nPower in solid_source_terms')
            return
          endif
        else if(me%ct(ii)%type=='stabilizer') then 
          ConductMet=1.0_dp/me%mat_stab%resistivity(me%StVar%SCtemp(i),B)
          idxST=ii
        endif
      enddo            
        
      if(abs(me%scen%ElCur)<1.0e-15_dp) then
        Joule =0.0_dp
        SJoule=0.0_dp
        DerSJoule=0.0_dp
      else if(CritCurr<=0.0_dp) then
        Joule = (me%scen%ElCur**2)/(me%ct(idxST)%area*ConductMet)
        SJoule=((me%scen%ElCur**2)/(me%ct(idxST)%area*AreaMS*rhoMS))/(ConductMet*cpMS)

        DerConductMet=-me%mat_stab%resistivity_der(me%StVar%SCtemp(i),B)/&
                      ((me%mat_stab%resistivity(me%StVar%SCtemp(i),B))**2)
        DerSJoule=-((me%scen%ElCur**2)/(me%ct(idxST)%area*AreaMS*rhoMS))*(DerConductMet*cpMS+&
                      ConductMet*dcpMSdT)/((ConductMet*cpMS)**2)
      else

        call supercond_current(abs(me%scen%ElCur),CritCurr,me%mat_supc%nPow,me%mat_supc%E0,me%ct(idxST)%area,&
                               ConductMet,rapp,DerImpFuntdRapp)
        if (sim_error > 0) return
        Joule =me%mat_supc%E0*(rapp**me%mat_supc%nPow)*abs(me%scen%ElCur)
        SJoule=me%mat_supc%E0*(rapp**me%mat_supc%nPow)*abs(me%scen%ElCur)/(AreaMS*rhoMS*cpMS)
              
        DerCritCurr=me%mat_supc%critical_current_density_der(me%StVar%SCtemp(i),B,St,Tc0,Bc)*me%ct(idxSC)%area
        DerConductMet=-me%mat_stab%resistivity_der(me%StVar%SCtemp(i),B)/&
                      ((me%mat_stab%resistivity(me%StVar%SCtemp(i),B))**2)

        denum=me%ct(idxST)%area*ConductMet*me%mat_supc%E0
        DRappdCritCurr=-(rapp/denum)/DerImpFuntdRapp
        DRappdConductMet=-(-CritCurr*rapp/(denum*ConductMet)+abs(me%scen%ElCur)/(denum*ConductMet))/DerImpFuntdRapp
        
        DRappdT=DRappdCritCurr*DerCritCurr+DRappdConductMet*DerConductMet
        
        Derivee=(me%mat_supc%nPow*(rapp**(me%mat_supc%nPow-1))*DRappdT*cpMS-(rapp**me%mat_supc%nPow)*dcpMSdT)/&
                (cpMS**2)
        
        DerSJoule=me%mat_supc%E0*abs(me%scen%ElCur)*Derivee/(AreaMS*rhoMS)
      endif

      if(me%SC_Prop%FSlink) then
        ! Communication to FS port associated to strand node i
        me%SC_Prop%thermP(i)%p%TempS=me%StVar%SCtemp(i)
        me%SC_Prop%thermP(i)%p%cpMS=cpMS
        me%SC_Prop%thermP(i)%p%dcpMSdT=dcpMSdT
      endif

end subroutine solid_source_terms
     
end module cmp_strand_source_terms_m