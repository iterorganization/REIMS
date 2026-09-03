! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_circulator_calc_m
    use cmp_circulator_init_m   
    implicit none

    type flux_and_derivatives_for_cold_circulator_t
      real(dp), dimension(Nb_VarC) :: flx,DerVitdBr_a,DerVitdBr_b
      real(dp), dimension(Nb_VarC,Nb_VarC) :: DerFlxdBr_a,DerFlxdBr_b
      real(dp) :: vit      
    end type flux_and_derivatives_for_cold_circulator_t

    type derivatives_involving_both_branches_t
      real(dp), dimension(Nb_VarC) :: DerVarIndUin,DerVarIndUout,DerMdotDerUin,DerMdotDerUout
      real(dp), dimension(Nb_VarC) :: DerHHHdUin,DerHHHdUout
    end type derivatives_involving_both_branches_t

contains

function imposed_mass_flow_rate(me,dPres,mdot0,dp0)
    type(circulator_t), intent(inout) :: me
    real(dp), intent(in) :: dPres,mdot0
    real(dp), intent(in), optional :: dp0
    real(dp) :: imposed_mass_flow_rate,eps
    
    if(me%SubType=="pump") then
        imposed_mass_flow_rate=mdot0
    else if(me%SubType=="compressor") then
        eps=1.0e-3_dp  
        if(dPres<=0.0_dp) then
            imposed_mass_flow_rate=mdot0
        else if(dPres>=dp0) then
            imposed_mass_flow_rate=eps*mdot0 !0.0_dp
        else
            imposed_mass_flow_rate=mdot0*(1.0_dp-(dPres/dp0)**2)
        endif
    endif
end function imposed_mass_flow_rate


function imposed_compression_energy(me,rostar,g,dPres,dp0)
    type(circulator_t), intent(inout) :: me
    real(dp), intent(in) :: rostar,g,dPres,dp0
    real(dp) :: imposed_compression_energy,eps,Ro_a,Ro_b
    
    if(me%SubType=="pump") then
        Ro_a=me%br(1)%p%Prim(Pri_ro)
        Ro_b=me%br(2)%p%Prim(Pri_ro)
        imposed_compression_energy=0.5_dp*(1.0_dp/Ro_a+1.0_dp/Ro_b)*dPres
    else if(me%SubType=="compressor") then
        eps=1.0e-3_dp  
        if(dPres<=0.0_dp) then
            imposed_compression_energy=dp0/(rostar*g)
        else if(dPres>=dp0) then
            imposed_compression_energy=eps*dp0/(rostar*g) !0.0_dp
        else
            imposed_compression_energy=dp0*(1.0_dp-(dPres/dp0)**2)/(rostar*g)
        endif
    endif
end function imposed_compression_energy    


subroutine mdot_derivatives(me,deriv_wr_both_br)
    type(circulator_t), intent(inout) :: me
    type(derivatives_involving_both_branches_t), intent(out) :: deriv_wr_both_br

    real(dp) :: dPres,dMdotdPin,dMdotdPout
    real(dp), dimension(Nb_VarC) :: dPindUin,dPoutdUout
    real(dp) :: kk, KKK, dedT, dPdT, HH, dPde

    if(me%SubType=="pump") then
        deriv_wr_both_br%DerMdotDerUin(:)=0.0_dp
        deriv_wr_both_br%DerMdotDerUout(:)=0.0_dp
    else if(me%SubType=="compressor") then
        dPres=me%br(2)%p%Prim(Pri_P)-me%br(1)%p%Prim(Pri_P)

        if(dPres<=0.0_dp) then
            deriv_wr_both_br%DerMdotDerUin(:)=0.0_dp
            deriv_wr_both_br%DerMdotDerUout(:)=0.0_dp
        else if(dPres>=me%dp0) then
            deriv_wr_both_br%DerMdotDerUin(:)=0.0_dp
            deriv_wr_both_br%DerMdotDerUout(:)=0.0_dp
        else
            dMdotdPin=2.0_dp*me%mdot0*dPres/(me%dp0**2)
            dMdotdPout=-2.0_dp*me%mdot0*dPres/(me%dp0**2)

            if(R_Correction) then
                dPindUin(1)=-0.5_dp*(me%br(1)%p%Cons(Con_Qdm)**2)/(me%br(1)%p%Cons(Con_Mas)**2)
                dPindUin(2)=me%br(1)%p%Cons(Con_Qdm)/me%br(1)%p%Cons(Con_Mas)
                dPindUin(3)=-1.0_dp
                dPindUin(4)=1.0_dp
                
                dPoutdUout(1)=-0.5_dp*(me%br(2)%p%Cons(Con_Qdm)**2)/(me%br(2)%p%Cons(Con_Mas)**2)
                dPoutdUout(2)=me%br(2)%p%Cons(Con_Qdm)/me%br(2)%p%Cons(Con_Mas)
                dPoutdUout(3)=-1.0_dp
                dPoutdUout(4)=1.0_dp         
            else
                HH=me%br(1)%p%Prim(Pri_e)+me%br(1)%p%Prim(Pri_P)/me%br(1)%p%Prim(Pri_ro)+0.5_dp*me%br(1)%p%Prim(Pri_u)**2
                call jacobian_roT(me%br(1)%p%Prim(Pri_ro), me%br(1)%p%Prim(Pri_T), dedT, dPdT)
                dPde=dPdT/dedT
                kk=dPde/me%br(1)%p%Prim(Pri_ro)
                KKK=me%br(1)%p%Prim(Pri_c)**2+kk*(me%br(1)%p%Prim(Pri_u)**2-HH)
                dPindUin(1)=KKK
                dPindUin(2)=-kk*me%br(1)%p%Prim(Pri_u)
                dPindUin(3)=kk
                dPindUin(4)=0.0_dp
                
                HH=me%br(2)%p%Prim(Pri_e)+me%br(2)%p%Prim(Pri_P)/me%br(2)%p%Prim(Pri_ro)+0.5_dp*me%br(2)%p%Prim(Pri_u)**2
                call jacobian_roT(me%br(2)%p%Prim(Pri_ro), me%br(2)%p%Prim(Pri_T), dedT, dPdT)
                dPde=dPdT/dedT
                kk=dPde/me%br(2)%p%Prim(Pri_ro)
                KKK=me%br(2)%p%Prim(Pri_c)**2+kk*(me%br(2)%p%Prim(Pri_u)**2-HH)
                dPoutdUout(1)=KKK
                dPoutdUout(2)=-kk*me%br(2)%p%Prim(Pri_u)
                dPoutdUout(3)=kk
                dPoutdUout(4)=0.0_dp 
            endif
            deriv_wr_both_br%DerMdotDerUin(:)=dMdotdPin*dPindUin(:)
            deriv_wr_both_br%DerMdotDerUout(:)=dMdotdPout*dPoutdUout(:)
        endif
    endif
end subroutine mdot_derivatives


subroutine HHH_derivatives(me,rostarIn,g,dRostardUin,deriv_wr_both_br)
      type(circulator_t), intent(inout) :: me
      real(dp), intent(in) :: rostarIn,g
      real(dp), dimension(:), intent(in) :: dRostardUin
      type(derivatives_involving_both_branches_t), intent(inout) :: deriv_wr_both_br

      real(dp) :: dPres,dHHHdPin,dHHHdPout,dHHHdRostarIn,eps,f_ro,Ro_a,Ro_b
      real(dp), dimension(Nb_VarC) :: dPindUin,dPoutdUout,df_ro_DUin,df_ro_DUout
      real(dp) :: kk, KKK, dedT, dPdT, HH, dPde

      if(me%SubType=="pump") then

        Ro_a=me%br(1)%p%Prim(Pri_ro)
        Ro_b=me%br(2)%p%Prim(Pri_ro)
        f_ro=(1.0_dp/Ro_a+1.0_dp/Ro_b)
        dPres=me%br(2)%p%Prim(Pri_P)-me%br(1)%p%Prim(Pri_P)

        df_ro_DUin(1)=-1.0_dp/(Ro_a**2)
        df_ro_DUin(2:4)=0.0_dp

        df_ro_DUout(1)=-1.0_dp/(Ro_b**2)
        df_ro_DUout(2:4)=0.0_dp

        if(R_Correction) then
          dPindUin(1)=-0.5_dp*(me%br(1)%p%Cons(Con_Qdm)**2)/(me%br(1)%p%Cons(Con_Mas)**2)
          dPindUin(2)=me%br(1)%p%Cons(Con_Qdm)/me%br(1)%p%Cons(Con_Mas)
          dPindUin(3)=-1.0_dp
          dPindUin(4)=1.0_dp
  
          dPoutdUout(1)=-0.5_dp*(me%br(2)%p%Cons(Con_Qdm)**2)/(me%br(2)%p%Cons(Con_Mas)**2)
          dPoutdUout(2)=me%br(2)%p%Cons(Con_Qdm)/me%br(2)%p%Cons(Con_Mas)
          dPoutdUout(3)=-1.0_dp
          dPoutdUout(4)=1.0_dp        
        else
          HH=me%br(1)%p%Prim(Pri_e)+me%br(1)%p%Prim(Pri_P)/me%br(1)%p%Prim(Pri_ro)+0.5_dp*me%br(1)%p%Prim(Pri_u)**2
          call jacobian_roT(me%br(1)%p%Prim(Pri_ro), me%br(1)%p%Prim(Pri_T), dedT, dPdT)
          dPde=dPdT/dedT
          kk=dPde/me%br(1)%p%Prim(Pri_ro)
          KKK=me%br(1)%p%Prim(Pri_c)**2+kk*(me%br(1)%p%Prim(Pri_u)**2-HH)
          dPindUin(1)=KKK
          dPindUin(2)=-kk*me%br(1)%p%Prim(Pri_u)
          dPindUin(3)=kk
          dPindUin(4)=0.0_dp
  
          HH=me%br(2)%p%Prim(Pri_e)+me%br(2)%p%Prim(Pri_P)/me%br(2)%p%Prim(Pri_ro)+0.5_dp*me%br(2)%p%Prim(Pri_u)**2
          call jacobian_roT(me%br(2)%p%Prim(Pri_ro), me%br(2)%p%Prim(Pri_T), dedT, dPdT)
          dPde=dPdT/dedT
          kk=dPde/me%br(2)%p%Prim(Pri_ro)
          KKK=me%br(2)%p%Prim(Pri_c)**2+kk*(me%br(2)%p%Prim(Pri_u)**2-HH)
          dPoutdUout(1)=KKK
          dPoutdUout(2)=-kk*me%br(2)%p%Prim(Pri_u)
          dPoutdUout(3)=kk
          dPoutdUout(4)=0.0_dp  
        endif

        deriv_wr_both_br%DerHHHdUin(:)=0.5_dp*(df_ro_DUin(:)*dPres+f_ro*(-dPindUin(:)))
        deriv_wr_both_br%DerHHHdUout(:)=0.5_dp*(df_ro_DUout(:)*dPres+f_ro*dPoutdUout(:))
        
      else if(me%SubType=="compressor") then

        eps=1.0e-3_dp
        dPres=me%br(2)%p%Prim(Pri_P)-me%br(1)%p%Prim(Pri_P)
  
        if(dPres<=0.0_dp) then
          deriv_wr_both_br%DerHHHdUin(:)=-(me%dp0/g)*dRostardUin(:)/(rostarIn**2)
          deriv_wr_both_br%DerHHHdUout(:)=0.0_dp
        else if(dPres>=me%dp0) then
          deriv_wr_both_br%DerHHHdUin(:)=-(eps*me%dp0/g)*dRostardUin(:)/(rostarIn**2)
          deriv_wr_both_br%DerHHHdUout(:)=0.0_dp
        else
          dHHHdPin=2.0_dp*(me%dp0/(rostarIn*g))*dPres/(me%dp0**2)
          dHHHdRostarIn=-(me%dp0/g)*(1.0_dp-(dPres/me%dp0)**2)/(rostarIn**2)
  
          dHHHdPout=-2.0_dp*(me%dp0/(rostarIn*g))*dPres/(me%dp0**2)
  
          if(R_Correction) then
            dPindUin(1)=-0.5_dp*(me%br(1)%p%Cons(Con_Qdm)**2)/(me%br(1)%p%Cons(Con_Mas)**2)
            dPindUin(2)=me%br(1)%p%Cons(Con_Qdm)/me%br(1)%p%Cons(Con_Mas)
            dPindUin(3)=-1.0_dp
            dPindUin(4)=1.0_dp
    
            dPoutdUout(1)=-0.5_dp*(me%br(2)%p%Cons(Con_Qdm)**2)/(me%br(2)%p%Cons(Con_Mas)**2)
            dPoutdUout(2)=me%br(2)%p%Cons(Con_Qdm)/me%br(2)%p%Cons(Con_Mas)
            dPoutdUout(3)=-1.0_dp
            dPoutdUout(4)=1.0_dp        
          else
            HH=me%br(1)%p%Prim(Pri_e)+me%br(1)%p%Prim(Pri_P)/me%br(1)%p%Prim(Pri_ro)+0.5_dp*me%br(1)%p%Prim(Pri_u)**2
            call jacobian_roT(me%br(1)%p%Prim(Pri_ro), me%br(1)%p%Prim(Pri_T), dedT, dPdT)
            dPde=dPdT/dedT
            kk=dPde/me%br(1)%p%Prim(Pri_ro)
            KKK=me%br(1)%p%Prim(Pri_c)**2+kk*(me%br(1)%p%Prim(Pri_u)**2-HH)
            dPindUin(1)=KKK
            dPindUin(2)=-kk*me%br(1)%p%Prim(Pri_u)
            dPindUin(3)=kk
            dPindUin(4)=0.0_dp
    
            HH=me%br(2)%p%Prim(Pri_e)+me%br(2)%p%Prim(Pri_P)/me%br(2)%p%Prim(Pri_ro)+0.5_dp*me%br(2)%p%Prim(Pri_u)**2
            call jacobian_roT(me%br(2)%p%Prim(Pri_ro), me%br(2)%p%Prim(Pri_T), dedT, dPdT)
            dPde=dPdT/dedT
            kk=dPde/me%br(2)%p%Prim(Pri_ro)
            KKK=me%br(2)%p%Prim(Pri_c)**2+kk*(me%br(2)%p%Prim(Pri_u)**2-HH)
            dPoutdUout(1)=KKK
            dPoutdUout(2)=-kk*me%br(2)%p%Prim(Pri_u)
            dPoutdUout(3)=kk
            dPoutdUout(4)=0.0_dp  
          endif
  
          deriv_wr_both_br%DerHHHdUin(:)=dHHHdPin*dPindUin(:)+dHHHdRostarIn*dRostardUin(:)
          deriv_wr_both_br%DerHHHdUout(:)=dHHHdPout*dPoutdUout(:)
        endif
      endif

end subroutine HHH_derivatives


subroutine incoming_branch_cold_circulator(me,mdot,dPres,VarIn,deriv_wr_both_br,ssIn,FluxIn)
      type(circulator_t), intent(inout) :: me
      real(dp), intent(in) :: mdot,dPres
      real(dp), intent(out) :: VarIn,ssIn
      type(derivatives_involving_both_branches_t), intent(inout) :: deriv_wr_both_br
      type(flux_and_derivatives_for_cold_circulator_t), intent(out) :: FluxIn

      real(dp) :: Rol,uL,pL,rL,cL,ML,Zl,Area,Sl,g
      real(dp) :: Rostar,pstar,Ustar,R_Correction_Star,eintStar,Tstar,Cstar
      real(dp) :: BB,HHH,EtotStar
      real(dp), dimension(Nb_VarC) :: U,Cons_star,DerUstardUl,DerRostardUl,DerRcorrstardUl,DerPstardUl,DerBBdUl
      real(dp), dimension(Nb_VarC) :: DerUstardUoth,DerRostardUoth,DerRcorrstardUoth,DerPstardUoth
      real(dp), dimension(Nb_VarC) :: dMdotdUl,dMdotdUoth
      real(dp), dimension(Nb_VarC,Nb_VarC) :: Jacob_star,DerConstardUl,DerConstardUoth
      real(dp) :: kk, KKK, dedT, dPdT, HH, dPde
      real(dp), dimension(Nb_VarP) :: PPP
      real(dp), dimension(Nb_VarP) :: PrimStar
      real(dp) :: dTdp_Ro,dTdRo_p,dRstardRo,dRstardT,Rstar

      Rol=me%br(1)%p%Prim(Pri_ro)
      uL=me%br(1)%p%Prim(Pri_u)
      pL=me%br(1)%p%Prim(Pri_P)
      rL=me%br(1)%p%Prim(Pri_R)
      cL=me%br(1)%p%Prim(Pri_c)
      ML=uL/cL
      Zl=Rol*cL
      Area=me%br(1)%p%Area
      U(:)=me%br(1)%p%Cons(:)

      ! Approximate RH relation
      Sl=-(abs(uL)+cL)
      Ustar=Sl*mdot/(mdot-Rol*(uL-Sl)*Area)                                          
      Rostar=Rol*(uL-Sl)/(Ustar-Sl)
      pstar=pL+Rol*(uL-Sl)*(uL-Ustar)
      call state_roP(Rostar, pstar, R_Correction_Star, eintStar, Tstar, Cstar)
      if(R_Correction) R_Correction_Star=rL*(uL-Sl)/(Ustar-Sl)
      eintStar=(R_Correction_Star-pstar)/Rostar
      EtotStar=eintStar+0.5_dp*Ustar**2

      FluxIn%flx(Con_Mas)=Rostar*Ustar
      FluxIn%flx(Con_Qdm)=Rostar*Ustar*Ustar+pstar
      FluxIn%flx(Con_Ene)=(Rostar*EtotStar+pstar)*Ustar
      FluxIn%flx(Con_R)=0.0_dp
      if(R_Correction) FluxIn%flx(Con_R)=R_Correction_Star*Ustar

      FluxIn%vit=Ustar

      !------------------------------------------------------------------!
      !                            Derivatives                           !
      dMdotdUl(:)=deriv_wr_both_br%DerMdotDerUin(:)
      DerUstardUl(1)=(Sl*dMdotdUl(1)*(mdot-(U(2)-U(1)*Sl)*Area)-Sl*mdot*(dMdotdUl(1)+Sl*Area))/((mdot-(U(2)-U(1)*Sl)*Area)**2)
      DerUstardUl(2)=(Sl*dMdotdUl(2)*(mdot-(U(2)-U(1)*Sl)*Area)-Sl*mdot*(dMdotdUl(2)-Area))/((mdot-(U(2)-U(1)*Sl)*Area)**2)
      DerUstardUl(3)=(Sl*dMdotdUl(3)*(mdot-(U(2)-U(1)*Sl)*Area)-Sl*mdot*dMdotdUl(3))/((mdot-(U(2)-U(1)*Sl)*Area)**2)
      DerUstardUl(4)=(Sl*dMdotdUl(4)*(mdot-(U(2)-U(1)*Sl)*Area)-Sl*mdot*dMdotdUl(4))/((mdot-(U(2)-U(1)*Sl)*Area)**2)
      
      DerRostardUl(1)=(-Sl*(Ustar-Sl)-(U(2)-U(1)*Sl)*DerUstardUl(1))/((Ustar-Sl)**2)
      DerRostardUl(2)=((Ustar-Sl)-(U(2)-U(1)*Sl)*DerUstardUl(2))/((Ustar-Sl)**2)
      DerRostardUl(3)=-(U(2)-U(1)*Sl)*DerUstardUl(3)/((Ustar-Sl)**2)
      DerRostardUl(4)=-(U(2)-U(1)*Sl)*DerUstardUl(4)/((Ustar-Sl)**2)

      if(R_Correction) then
        Cons_star(Con_Mas)=Rostar
        Cons_star(Con_Qdm)=Rostar*Ustar
        Cons_star(Con_Ene)=Rostar*EtotStar
        Cons_star(Con_R)=R_Correction_Star
        call jacobian_star(Cons_star,Jacob_star)
  
        BB=U(4)-U(3)+0.5_dp*(U(2)**2)/U(1)
        DerBBdUl(1)=-0.5_dp*(U(2)**2)/(U(1)**2)
        DerBBdUl(2)=U(2)/U(1)
        DerBBdUl(3)=-1.0_dp
        DerBBdUl(4)=1.0_dp

        DerPstardUl(1)=DerBBdUl(1)-Sl*(U(2)/U(1)-Ustar)+(U(2)-U(1)*Sl)*(-U(2)/(U(1)**2)-DerUstardUl(1))
        DerPstardUl(2)=DerBBdUl(2)+(U(2)/U(1)-Ustar)+(U(2)-U(1)*Sl)*(1.0_dp/U(1)-DerUstardUl(2))
        DerPstardUl(3)=DerBBdUl(3)-(U(2)-U(1)*Sl)*DerUstardUl(3)
        DerPstardUl(4)=DerBBdUl(4)-(U(2)-U(1)*Sl)*DerUstardUl(4)

        DerRcorrstardUl(1)=-U(4)*(U(2)*(Ustar-Sl)/(U(1)**2)+(U(2)/U(1)-Sl)*DerUstardUl(1))/((Ustar-Sl)**2)
        DerRcorrstardUl(2)=U(4)*((Ustar-Sl)/U(1)-(U(2)/U(1)-Sl)*DerUstardUl(2))/((Ustar-Sl)**2)
        DerRcorrstardUl(3)=-U(4)*(U(2)/U(1)-Sl)*DerUstardUl(3)/((Ustar-Sl)**2)
        DerRcorrstardUl(4)=(U(2)/U(1)-Sl)*(Ustar-Sl-U(4)*DerUstardUl(4))/((Ustar-Sl)**2)  

        DerConstardUl(:,1)=DerRostardUl(:)
        DerConstardUl(:,2)=Rostar*DerUstardUl(:)+DerRostardUl(:)*Ustar
        DerConstardUl(:,3)=DerRcorrstardUl(:)+0.5_dp*(DerRostardUl(:)*Ustar*Ustar+2.0_dp*Rostar*Ustar*DerUstardUl(:))-DerPstardUl(:)
        DerConstardUl(:,4)=DerRcorrstardUl(:)
  
        FluxIn%DerFlxdBr_a(:,:)=matmul(DerConstardUl(:,:),Jacob_star(:,:))
        FluxIn%DerVitdBr_a(:)=DerUstardUl(:)            

      else

        PrimStar(Pri_ro)=Rostar
        PrimStar(Pri_u)=Ustar
        PrimStar(Pri_p)=pstar
        PrimStar(Pri_e)=eintStar
        PrimStar(Pri_T)=Tstar
        PrimStar(Pri_c)=Cstar
        call jacobian_star_Req_rem(PrimStar,Jacob_star)
  
        PPP(:)=me%br(1)%p%Prim(:)
        HH=PPP(Pri_e)+0.5_dp*PPP(Pri_u)**2+PPP(Pri_p)/PPP(Pri_ro)
        call jacobian_roT(PPP(Pri_ro), PPP(Pri_T), dedT, dPdT)
        dPde=dPdT/dedT
        kk=dPde/PPP(Pri_ro)
        KKK=PPP(Pri_c)**2+kk*(PPP(Pri_u)**2-HH)

        BB=pL
        DerBBdUl(1)=KKK
        DerBBdUl(2)=-kk*PPP(Pri_u)
        DerBBdUl(3)=kk
        DerBBdUl(4)=0.0_dp

        DerPstardUl(1)=DerBBdUl(1)-Sl*(U(2)/U(1)-Ustar)+(U(2)-U(1)*Sl)*(-U(2)/(U(1)**2)-DerUstardUl(1))
        DerPstardUl(2)=DerBBdUl(2)+(U(2)/U(1)-Ustar)+(U(2)-U(1)*Sl)*(1.0_dp/U(1)-DerUstardUl(2))
        DerPstardUl(3)=DerBBdUl(3)-(U(2)-U(1)*Sl)*DerUstardUl(3)
        DerPstardUl(4)=DerBBdUl(4)-(U(2)-U(1)*Sl)*DerUstardUl(4)

        call jacobian_roT(Rostar, Tstar, dedT, dPdT, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
        DerRcorrstardUl(1)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUl(1)+dRstardT*dTdp_Ro*DerPstardUl(1)
        DerRcorrstardUl(2)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUl(2)+dRstardT*dTdp_Ro*DerPstardUl(2)
        DerRcorrstardUl(3)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUl(3)+dRstardT*dTdp_Ro*DerPstardUl(3)
        DerRcorrstardUl(4)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUl(4)+dRstardT*dTdp_Ro*DerPstardUl(4)

        DerConstardUl(:,1)=DerRostardUl(:)
        DerConstardUl(:,2)=Rostar*DerUstardUl(:)+DerRostardUl(:)*Ustar
        DerConstardUl(:,3)=DerRcorrstardUl(:)+0.5_dp*(DerRostardUl(:)*Ustar*Ustar+2.0_dp*Rostar*Ustar*DerUstardUl(:))-DerPstardUl(:)
        DerConstardUl(:,4)=0.0_dp
  
        FluxIn%DerFlxdBr_a(:,:)=matmul(DerConstardUl(:,:),Jacob_star(:,:))
        FluxIn%DerFlxdBr_a(:,4)=0.0_dp
        FluxIn%DerVitdBr_a(:)=DerUstardUl(:)   

      endif
      
      !----------------------------------------------------------------------------------------------------------
      ! Contributions of the CC outlet
      dMdotdUoth(:)=deriv_wr_both_br%DerMdotDerUout(:)
      DerUstardUoth(:)=(Sl*dMdotdUoth(:)*(mdot-(U(2)-U(1)*Sl)*Area)-Sl*mdot*dMdotdUoth(:))/((mdot-(U(2)-U(1)*Sl)*Area)**2)
      DerRostardUoth(:)=-(U(2)-U(1)*Sl)*DerUstardUoth(:)/((Ustar-Sl)**2)
      DerPstardUoth(:)=-(U(2)-U(1)*Sl)*DerUstardUoth(:)

      if(R_Correction) then
        DerRcorrstardUoth(:)=-U(4)*(U(2)/U(1)-Sl)*DerUstardUoth(:)/((Ustar-Sl)**2)
      else
        call jacobian_roT(Rostar, Tstar, dedT, dPdT, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
        DerRcorrstardUoth(:)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUoth(:)+dRstardT*dTdp_Ro*DerPstardUoth(:)  ! DerRcorrstardUoth(4) already equals to 0
      endif

      DerConstardUoth(:,1)=DerRostardUoth(:)
      DerConstardUoth(:,2)=Rostar*DerUstardUoth(:)+DerRostardUoth(:)*Ustar
      DerConstardUoth(:,3)=DerRcorrstardUoth(:)+0.5_dp*(DerRostardUoth(:)*Ustar*Ustar+2.0_dp*Rostar*Ustar*DerUstardUoth(:))-&
                           DerPstardUoth(:)
      DerConstardUoth(:,4)=DerRcorrstardUoth(:)

      FluxIn%DerFlxdBr_b(:,:)=matmul(DerConstardUoth(:,:),Jacob_star(:,:))
      if(.not. R_Correction) FluxIn%DerFlxdBr_b(:,4)=0.0_dp
      FluxIn%DerVitdBr_b(:)=DerUstardUoth(:)
      !------------------------------------------------------------------------------------------------------

      g=9.81_dp
      HHH=imposed_compression_energy(me,Rostar,g,dPres,me%dp0)

      if(me%SubType=="pump") then
        VarIn=R_Correction_Star/Rostar+0.5_dp*Ustar**2+HHH
  
        call HHH_derivatives(me,rostar,g,DerRostardUl,deriv_wr_both_br)
  
        deriv_wr_both_br%DerVarIndUin(:)=(DerRcorrstardUl(:)*Rostar-R_Correction_Star*DerRostardUl(:))/&
                                         (Rostar**2)+Ustar*DerUstardUl(:)+deriv_wr_both_br%DerHHHdUin(:)
  
        deriv_wr_both_br%DerVarIndUout(:)=(DerRcorrstardUoth(:)*Rostar-R_Correction_Star*DerRostardUoth(:))/&
                                          (Rostar**2)+Ustar*DerUstardUoth(:)+deriv_wr_both_br%DerHHHdUout(:)

      else if(me%SubType=="compressor") then
        VarIn=R_Correction_Star/Rostar+0.5_dp*Ustar**2+g*HHH
  
        call HHH_derivatives(me,rostar,g,DerRostardUl,deriv_wr_both_br)
  
        deriv_wr_both_br%DerVarIndUin(:)=(DerRcorrstardUl(:)*Rostar-R_Correction_Star*DerRostardUl(:))/&
                                         (Rostar**2)+Ustar*DerUstardUl(:)+g*deriv_wr_both_br%DerHHHdUin(:)
  
        deriv_wr_both_br%DerVarIndUout(:)=(DerRcorrstardUoth(:)*Rostar-R_Correction_Star*DerRostardUoth(:))/&
                                          (Rostar**2)+Ustar*DerUstardUoth(:)+g*deriv_wr_both_br%DerHHHdUout(:)

      endif

      ssIn=abs(uL)+cL

end subroutine incoming_branch_cold_circulator

      
subroutine jacobian_star(Cons_star,Jacob_star)
      
        real(dp), dimension(:), intent(in) :: Cons_star
        real(dp), dimension(:,:), intent(out) :: Jacob_star
    
        real(dp), dimension(Nb_VarC) :: U 
    
        U(:)=Cons_star(:)
    
        Jacob_star(1,1)=0.0_dp
        Jacob_star(2,1)=1.0_dp
        Jacob_star(3,1)=0.0_dp
        Jacob_star(4,1)=0.0_dp
    
        Jacob_star(1,2)=-(3.0_dp/2.0_dp)*((U(2)**2)/(U(1)**2))
        Jacob_star(2,2)=3.0_dp*U(2)/U(1)
        Jacob_star(3,2)=-1.0_dp
        Jacob_star(4,2)=1.0_dp
    
        Jacob_star(1,3)=-U(4)*U(2)/(U(1)**2)-((U(2)**3)/(U(1)**3))
        Jacob_star(2,3)=U(4)/U(1)+(3.0_dp/2.0_dp)*((U(2)**2)/(U(1)**2))
        Jacob_star(3,3)=0.0_dp
        Jacob_star(4,3)=U(2)/U(1)
        
        Jacob_star(1,4)=-U(4)*U(2)/(U(1)**2)
        Jacob_star(2,4)=U(4)/U(1)
        Jacob_star(3,4)=0.0_dp
        Jacob_star(4,4)=U(2)/U(1)
        
end subroutine jacobian_star        

      
subroutine jacobian_star_Req_rem(PrimStar,Jacob_star)
      
          real(dp), dimension(:), intent(in) :: PrimStar
          real(dp), dimension(:,:), intent(out) :: Jacob_star

          real(dp) :: HH,kk,KKK,dPde,dPdT,dedT

          HH=PrimStar(Pri_e)+0.5_dp*PrimStar(Pri_u)**2+PrimStar(Pri_p)/PrimStar(Pri_ro)
          call jacobian_roT(PrimStar(Pri_ro), PrimStar(Pri_T), dedT, dPdT)
          dPde=dPdT/dedT
          kk=dPde/PrimStar(Pri_ro)
          KKK=PrimStar(Pri_c)**2+kk*(PrimStar(Pri_u)**2-HH)
      
          Jacob_star(1,1)=0.0_dp
          Jacob_star(2,1)=1.0_dp
          Jacob_star(3,1)=0.0_dp
          Jacob_star(4,1)=0.0_dp
      
          Jacob_star(1,2)=KKK-PrimStar(Pri_u)**2
          Jacob_star(2,2)=PrimStar(Pri_u)*(2.0_dp-kk)
          Jacob_star(3,2)=kk
          Jacob_star(4,2)=0.0_dp
      
          Jacob_star(1,3)=(KKK-HH)*PrimStar(Pri_u)
          Jacob_star(2,3)=HH-kk*PrimStar(Pri_u)**2
          Jacob_star(3,3)=PrimStar(Pri_u)*(1.0_dp+kk)
          Jacob_star(4,3)=0.0_dp
          
          Jacob_star(1,4)=0.0_dp
          Jacob_star(2,4)=0.0_dp
          Jacob_star(3,4)=0.0_dp
          Jacob_star(4,4)=0.0_dp
          
end subroutine jacobian_star_Req_rem            


subroutine outgoing_branch_cold_circulator(me,mdot,VarIn,deriv_wr_both_br,ssOut,FluxOut)

      type(circulator_t), intent(inout) :: me
      real(dp), intent(in) :: mdot,VarIn
      type(derivatives_involving_both_branches_t), intent(in) :: deriv_wr_both_br
      real(dp), intent(out) :: ssOut
      type(flux_and_derivatives_for_cold_circulator_t), intent(out) :: FluxOut

      real(dp) :: Ror,Pr,Ur,Cr,Zr,pmin,pmax,f_pmin,f_pmax,dpp
      real(dp) :: fct,Utest,Rotest,Ptest,Ttest
      real(dp) :: Pstar,Ustar,Rostar,R_Correction_Star,estar,Tstar,Cstar
      real(dp) :: dcdro,dcdp,DerhdPstar,BB,dTdp_Ro,dTdRo_p,dRstardRo,dRstardT,Rstar
      real(dp), dimension(Nb_VarC) :: dc,Cons_star,DerBBdUr,DerhDur,DerRcorrstardUr,DerRostardUr,DerUstardUr,DerhdUoth
      real(dp), dimension(Nb_VarC) :: DerRcorrstardUr_Glob,DerRostardUr_Glob,DerUstardUr_Glob,DerPstardUr_Glob
      real(dp), dimension(Nb_VarC) :: dRodUoth,dRdUoth
      real(dp), dimension(Nb_VarC) :: DerRcorrstardUoth_Glob,DerRostardUoth_Glob,DerUstardUoth_Glob,DerPstardUoth_Glob
      real(dp), dimension(Nb_VarC,Nb_VarC) :: Jacob_star,DerConstardUr,DerConstardUoth
      real(dp) :: kk, KKK, dedT, dPdT, HH, dPde
      real(dp), dimension(Nb_VarP) :: PPP
      real(dp), dimension(Nb_VarP) :: PrimStar
      integer :: i,nbpoints

      Ror=me%br(2)%p%Prim(Pri_ro)
      Pr=me%br(2)%p%Prim(Pri_P)
      Ur=me%br(2)%p%Prim(Pri_u)
      Cr=me%br(2)%p%Prim(Pri_c)
      Zr=Ror*Cr

      me%dynCC%mdot=mdot
      me%dynCC%VarIn=VarIn
      me%dynCC%area=me%br(2)%p%Area
      me%dynCC%Pr=Pr
      me%dynCC%Ur=Ur
      me%dynCC%Zr=Zr

      ! Preparatory work for the solution of F to be found, knowing the fact that F is diminishing
      fct=-1.0_dp
      nbpoints=1000
      pmin=Pr/1.5_dp
      pmax=Pr*1.5_dp
      dpp=(pmax-pmin)/(nbpoints-1)
      do i=1,nbpoints
        Ptest=Pmin+dpp*(i-1)
        Utest=Ur+(Ptest-Pr)/Zr
        Rotest=mdot/(Utest*me%br(2)%p%Area)
        if(Rotest>1.0e-3_dp .and. Rotest<290.0_dp) then
          Ttest=T_roP(Rotest,Ptest)
          if(Ttest>2.5_dp .and. Ttest<500.0_dp) then
            fct=me%dynCC%f(Ptest)
            if(fct>0.0_dp) exit
          endif
        endif
      enddo

      if(fct<=0.0_dp .or. i==nbpoints) then
        fct=-1.0_dp
        nbpoints=1000000
        pmin=Pr/1.5_dp
        pmax=Pr*1.5_dp
        dpp=(pmax-pmin)/(nbpoints-1)
        do i=1,nbpoints
          Ptest=Pmin+dpp*(i-1)
          Utest=Ur+(Ptest-Pr)/Zr
          Rotest=mdot/(Utest*me%br(2)%p%Area)
          if(Rotest>1.0e-3_dp .and. Rotest<290.0_dp) then
            Ttest=T_roP(Rotest,Ptest)
            if(Ttest>2.5_dp .and. Ttest<500.0_dp) then
              fct=me%dynCC%f(Ptest)
              if(fct>0.0_dp) exit
            endif
          endif
        enddo
        if(fct<=0.0_dp .or. i==nbpoints) then
            call set_error('initial bracket search failed in outgoing_branch_cold_circulator')
            return
        endif     
      endif

      ! Dichotomy
      pmin=Ptest
      f_pmin=me%dynCC%f(pmin)

      pmax=Ptest*1.1_dp
      f_pmax=me%dynCC%f(pmax)

      if(f_pmin*f_pmax>0.0_dp) then
        call set_error('initial range does not contain solution in outgoing_branch_cold_circulator')
        return
      else
        !Brent's method
        Pstar=zero(me%dynCC,pmin,pmax,1.0e-15_dp,1.0e-2_dp)        
      endif

      Ustar=Ur+(Pstar-Pr)/Zr   
      Rostar=mdot/(Ustar*me%br(2)%p%Area)
      call state_roP(Rostar, Pstar, R_Correction_Star, estar, Tstar, Cstar)
      if(Tstar<2.0_dp .or. Tstar>50.0_dp) then
        call set_error('unphysical temperature in outgoing_branch_cold_circulator')
        return
      endif

      FluxOut%flx(Con_Mas)=Rostar*Ustar
      FluxOut%flx(Con_Qdm)=Rostar*Ustar*Ustar+Pstar
      FluxOut%flx(Con_Ene)=(Rostar*(estar+0.5_dp*Ustar*Ustar)+Pstar)*Ustar
      FluxOut%flx(Con_R)=0.0_dp
      if(R_Correction) FluxOut%flx(Con_R)=R_Correction_Star*Ustar
  
      FluxOut%vit=Ustar

      !------------------------------------------------------------------!
      !                            Derivatives                           !
      !------------------------------------------------------------------!
      call dc2_roT(Ror, me%br(2)%p%Prim(Pri_T), dcdro, dcdp)
      dcdro = dcdro/(2.0_dp*Cr)
      dcdp  = dcdp /(2.0_dp*Cr)

      if(R_Correction) then

        dc(Con_Mas)=dcdro-0.5_dp*dcdp*(me%br(2)%p%Cons(Con_Qdm)**2)/(me%br(2)%p%Cons(Con_Mas)**2)
        dc(Con_Qdm)=dcdp*me%br(2)%p%Cons(Con_Qdm)/me%br(2)%p%Cons(Con_Mas)
        dc(Con_Ene)=-dcdp
        dc(Con_R)=dcdp

        Cons_star(Con_Mas)=Rostar
        Cons_star(Con_Qdm)=Rostar*Ustar
        Cons_star(Con_Ene)=Rostar*(estar+0.5_dp*Ustar*Ustar)
        Cons_star(Con_R)=R_Correction_Star
        call jacobian_star(Cons_star,Jacob_star)
    
        DerhdPstar=me%dynCC%df(Pstar)
  
        BB=pstar-(me%br(2)%p%Cons(Con_R)-me%br(2)%p%Cons(Con_Ene)+0.5_dp*(me%br(2)%p%Cons(Con_Qdm)**2)/me%br(2)%p%Cons(Con_Mas))
        DerBBdUr(1)=0.5_dp*(me%br(2)%p%Cons(Con_Qdm)**2)/(me%br(2)%p%Cons(Con_Mas)**2)
        DerBBdUr(2)=-me%br(2)%p%Cons(Con_Qdm)/me%br(2)%p%Cons(Con_Mas)
        DerBBdUr(3)=1.0_dp
        DerBBdUr(4)=-1.0_dp
  
        DerUstardUr(1)=-1.0_dp/(me%br(2)%p%Cons(Con_Mas)**2)*(me%br(2)%p%Cons(Con_Qdm)+BB/Cr)+&
                       +((DerBBdUr(1)*Cr-BB*dc(1))/(Cr**2))/me%br(2)%p%Cons(Con_Mas)
        DerUstardUr(2)=(1.0_dp+(DerBBdUr(2)*Cr-BB*dc(2))/(Cr**2))/me%br(2)%p%Cons(Con_Mas)
        DerUstardUr(3)=((DerBBdUr(3)*Cr-BB*dc(3))/(Cr**2))/me%br(2)%p%Cons(Con_Mas)
        DerUstardUr(4)=((DerBBdUr(4)*Cr-BB*dc(4))/(Cr**2))/me%br(2)%p%Cons(Con_Mas)
        
        DerRostardUr(:)=(deriv_wr_both_br%DerMdotDerUout(:)*Ustar*me%br(2)%p%Area-mdot*DerUstardUr(:)*me%br(2)%p%Area)/&
                        ((Ustar*me%br(2)%p%Area)**2)
  
        call jacobian_roT(rostar, Tstar, dedT, dPdT, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
        DerRcorrstardUr(:)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUr(:)
        
        DerhDur(:)=-((DerRcorrstardUr(:)*rostar-R_Correction_Star*DerRostardUr(:))/(rostar**2)+Ustar*DerUstardUr(:)-&
                     deriv_wr_both_br%DerVarIndUout(:))/1.0e3_dp
        DerPstardUr_Glob(:)=-DerhDur(:)/DerhdPstar
  
        BB=pstar-(me%br(2)%p%Cons(Con_R)-me%br(2)%p%Cons(Con_Ene)+0.5_dp*(me%br(2)%p%Cons(Con_Qdm)**2)/me%br(2)%p%Cons(Con_Mas))
        DerBBdUr(1)=DerPstardUr_Glob(1)+0.5_dp*(me%br(2)%p%Cons(Con_Qdm)**2)/(me%br(2)%p%Cons(Con_Mas)**2)
        DerBBdUr(2)=DerPstardUr_Glob(2)-me%br(2)%p%Cons(Con_Qdm)/me%br(2)%p%Cons(Con_Mas)
        DerBBdUr(3)=DerPstardUr_Glob(3)+1.0_dp
        DerBBdUr(4)=DerPstardUr_Glob(4)-1.0_dp
    
        DerUstardUr_Glob(1)=-1.0_dp/(me%br(2)%p%Cons(Con_Mas)**2)*(me%br(2)%p%Cons(Con_Qdm)+BB/Cr)+&
                            +((DerBBdUr(1)*Cr-BB*dc(1))/(Cr**2))/me%br(2)%p%Cons(Con_Mas)
        DerUstardUr_Glob(2)=(1.0_dp+(DerBBdUr(2)*Cr-BB*dc(2))/(Cr**2))/me%br(2)%p%Cons(Con_Mas)
        DerUstardUr_Glob(3)=((DerBBdUr(3)*Cr-BB*dc(3))/(Cr**2))/me%br(2)%p%Cons(Con_Mas)
        DerUstardUr_Glob(4)=((DerBBdUr(4)*Cr-BB*dc(4))/(Cr**2))/me%br(2)%p%Cons(Con_Mas)
  
        DerRostardUr_Glob(:)=(deriv_wr_both_br%DerMdotDerUout(:)*Ustar*me%br(2)%p%Area-mdot*DerUstardUr_Glob(:)*me%br(2)%p%Area)/&
                             ((Ustar*me%br(2)%p%Area)**2)
  
        DerRcorrstardUr_Glob(:)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUr_Glob(:)+dRstardT*dTdp_Ro*DerPstardUr_Glob(:)
  
        DerConstardUr(:,1)=DerRostardUr_Glob(:)
        DerConstardUr(:,2)=Rostar*DerUstardUr_Glob(:)+DerRostardUr_Glob(:)*Ustar
        DerConstardUr(:,3)=DerRcorrstardUr_Glob(:)+0.5_dp*(DerRostardUr_Glob(:)*Ustar*Ustar+&
                           2.0_dp*Rostar*Ustar*DerUstardUr_Glob(:))-DerPstardUr_Glob(:)
        DerConstardUr(:,4)=DerRcorrstardUr_Glob(:)
  
        FluxOut%DerFlxdBr_b(:,:)=matmul(DerConstardUr(:,:),Jacob_star(:,:))
        FluxOut%DerVitdBr_b(:)=DerUstardUr_Glob(:)
  
        ! Contribution of the cold circulator inlet
        dRodUoth(:)=deriv_wr_both_br%DerMdotDerUin(:)/(Ustar*me%br(2)%p%Area)
        dRdUoth(:)=(dRstardRo+dRstardT*dTdRo_p)*dRodUoth(:)
        DerhdUoth(:)=-((dRdUoth(:)*Rostar-R_Correction_Star*dRodUoth(:))/(Rostar**2)-deriv_wr_both_br%DerVarIndUin(:))/1.0e3_dp
  
        DerPstardUoth_Glob(:)=-DerhDuOth(:)/DerhdPstar
        DerUstardUoth_Glob(:)=DerPstardUoth_Glob(:)/Zr
        DerRostardUoth_Glob(:)=(deriv_wr_both_br%DerMdotDerUin(:)*Ustar*me%br(2)%p%Area-mdot*DerUstardUoth_Glob(:)*me%br(2)%p%Area)/&
                               ((Ustar*me%br(2)%p%Area)**2)
        DerRcorrstardUoth_Glob(:)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUoth_Glob(:)+dRstardT*dTdp_Ro*DerPstardUoth_Glob(:)
  
        DerConstardUoth(:,1)=DerRostardUoth_Glob(:)
        DerConstardUoth(:,2)=Rostar*DerUstardUoth_Glob(:)+DerRostardUoth_Glob(:)*Ustar
        DerConstardUoth(:,3)=DerRcorrstardUoth_Glob(:)+0.5_dp*(DerRostardUoth_Glob(:)*Ustar*Ustar+&
                           2.0_dp*Rostar*Ustar*DerUstardUoth_Glob(:))-DerPstardUoth_Glob(:)
        DerConstardUoth(:,4)=DerRcorrstardUoth_Glob(:)
  
        FluxOut%DerFlxdBr_a(:,:)=matmul(DerConstardUoth(:,:),Jacob_star(:,:))
        FluxOut%DerVitdBr_a(:)=DerUstardUoth_Glob(:)             

      else
        
        PrimStar(Pri_ro)=Rostar
        PrimStar(Pri_u)=Ustar
        PrimStar(Pri_p)=pstar
        PrimStar(Pri_e)=estar
        PrimStar(Pri_T)=Tstar
        PrimStar(Pri_c)=Cstar
        call jacobian_star_Req_rem(PrimStar,Jacob_star)

        DerhdPstar=me%dynCC%df(Pstar)

        PPP(:)=me%br(2)%p%Prim(:)
        HH=PPP(Pri_e)+0.5_dp*PPP(Pri_u)**2+PPP(Pri_p)/PPP(Pri_ro)
        call jacobian_roT(PPP(Pri_ro), PPP(Pri_T), dedT, dPdT)
        dPde=dPdT/dedT
        kk=dPde/PPP(Pri_ro)
        KKK=PPP(Pri_c)**2+kk*(PPP(Pri_u)**2-HH)

        BB=Pstar-Pr
        DerBBdUr(1)=-KKK
        DerBBdUr(2)=kk*PPP(Pri_u)
        DerBBdUr(3)=-kk
        DerBBdUr(4)=0.0_dp

        dc(Con_Mas)=dcdro+dcdp*KKK
        dc(Con_Qdm)=-dcdp*kk*PPP(Pri_u)
        dc(Con_Ene)=dcdp*kk
        dc(Con_R)=0.0_dp         

        DerUstardUr(1)=-PPP(Pri_u)/PPP(Pri_ro)+(DerBBdUr(1)*PPP(Pri_ro)*Cr-BB*(Cr+PPP(Pri_ro)*dc(1)))/((PPP(Pri_ro)*Cr)**2)
        DerUstardUr(2)=(1.0_dp/PPP(Pri_ro))*(1.0_dp+(DerBBdUr(2)*Cr-BB*dc(2))/(Cr**2))
        DerUstardUr(3)=(1.0_dp/PPP(Pri_ro))*(DerBBdUr(3)*Cr-BB*dc(3))/(Cr**2)
        DerUstardUr(4)=0.0_dp        

        DerRostardUr(:)=(deriv_wr_both_br%DerMdotDerUout(:)*Ustar*me%br(2)%p%Area-mdot*DerUstardUr(:)*me%br(2)%p%Area)/&
                ((Ustar*me%br(2)%p%Area)**2)

        call jacobian_roT(rostar, Tstar, dedT, dPdT, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
        DerRcorrstardUr(:)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUr(:)
        
        DerhDur(:)=-((DerRcorrstardUr(:)*rostar-R_Correction_Star*DerRostardUr(:))/(rostar**2)+Ustar*DerUstardUr(:)-&
                     deriv_wr_both_br%DerVarIndUout(:))/1.0e3_dp
        DerPstardUr_Glob(:)=-DerhDur(:)/DerhdPstar                

        BB=Pstar-Pr
        DerBBdUr(1)=DerPstardUr_Glob(1)-KKK
        DerBBdUr(2)=DerPstardUr_Glob(2)+kk*PPP(Pri_u)
        DerBBdUr(3)=DerPstardUr_Glob(3)-kk
        DerBBdUr(4)=0.0_dp
    
        DerUstardUr_Glob(1)=-PPP(Pri_u)/PPP(Pri_ro)+(DerBBdUr(1)*PPP(Pri_ro)*Cr-BB*(Cr+PPP(Pri_ro)*dc(1)))/((PPP(Pri_ro)*Cr)**2)
        DerUstardUr_Glob(2)=(1.0_dp/PPP(Pri_ro))*(1.0_dp+(DerBBdUr(2)*Cr-BB*dc(2))/(Cr**2))
        DerUstardUr_Glob(3)=(1.0_dp/PPP(Pri_ro))*(DerBBdUr(3)*Cr-BB*dc(3))/(Cr**2)
        DerUstardUr_Glob(4)=0.0_dp
  
        DerRostardUr_Glob(:)=(deriv_wr_both_br%DerMdotDerUout(:)*Ustar*me%br(2)%p%Area-mdot*DerUstardUr_Glob(:)*me%br(2)%p%Area)/&
                             ((Ustar*me%br(2)%p%Area)**2)
  
        DerRcorrstardUr_Glob(:)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUr_Glob(:)+dRstardT*dTdp_Ro*DerPstardUr_Glob(:)

        DerConstardUr(:,1)=DerRostardUr_Glob(:)
        DerConstardUr(:,2)=Rostar*DerUstardUr_Glob(:)+DerRostardUr_Glob(:)*Ustar
        DerConstardUr(:,3)=DerRcorrstardUr_Glob(:)+0.5_dp*(DerRostardUr_Glob(:)*Ustar*Ustar+&
                           2.0_dp*Rostar*Ustar*DerUstardUr_Glob(:))-DerPstardUr_Glob(:)
        DerConstardUr(:,4)=0.0_dp             
 
        FluxOut%DerFlxdBr_b(:,:)=matmul(DerConstardUr(:,:),Jacob_star(:,:))
        FluxOut%DerFlxdBr_b(:,4)=0.0_dp

        FluxOut%DerVitdBr_b(:)=DerUstardUr_Glob(:)

        ! Contribution of the cold circulator inlet
        dRodUoth(:)=deriv_wr_both_br%DerMdotDerUin(:)/(Ustar*me%br(2)%p%Area)
        dRdUoth(:)=(dRstardRo+dRstardT*dTdRo_p)*dRodUoth(:)
        DerhdUoth(:)=-((dRdUoth(:)*Rostar-R_Correction_Star*dRodUoth(:))/(Rostar**2)-deriv_wr_both_br%DerVarIndUin(:))/1.0e3_dp
  
        DerPstardUoth_Glob(:)=-DerhDuOth(:)/DerhdPstar
        DerUstardUoth_Glob(:)=DerPstardUoth_Glob(:)/Zr
        DerRostardUoth_Glob(:)=(deriv_wr_both_br%DerMdotDerUin(:)*Ustar*me%br(2)%p%Area-mdot*DerUstardUoth_Glob(:)*me%br(2)%p%Area)/&
                               ((Ustar*me%br(2)%p%Area)**2)
        DerRcorrstardUoth_Glob(:)=(dRstardRo+dRstardT*dTdRo_p)*DerRostardUoth_Glob(:)+dRstardT*dTdp_Ro*DerPstardUoth_Glob(:)
  
        DerConstardUoth(:,1)=DerRostardUoth_Glob(:)
        DerConstardUoth(:,2)=Rostar*DerUstardUoth_Glob(:)+DerRostardUoth_Glob(:)*Ustar
        DerConstardUoth(:,3)=DerRcorrstardUoth_Glob(:)+0.5_dp*(DerRostardUoth_Glob(:)*Ustar*Ustar+&
                           2.0_dp*Rostar*Ustar*DerUstardUoth_Glob(:))-DerPstardUoth_Glob(:)
        DerConstardUoth(:,4)=0.0_dp
  
        FluxOut%DerFlxdBr_a(:,:)=matmul(DerConstardUoth(:,:),Jacob_star(:,:))
        FluxOut%DerFlxdBr_a(:,4)=0.0_dp

        FluxOut%DerVitdBr_a(:)=DerUstardUoth_Glob(:)   

      endif
 
      ssOut=abs(Ur)+Cr

end subroutine outgoing_branch_cold_circulator


subroutine circulator_resolution_from_and_to_ports(me)
    type(circulator_t), intent(inout) :: me

    real(dp) :: Pin,Pout,dPres,mdot,VarIn,ssOut,ssIn
    type(flux_and_derivatives_for_cold_circulator_t) :: FluxIn,Fluxout
    type(derivatives_involving_both_branches_t) :: deriv_wr_both_br

    Pin=me%br(1)%p%Prim(Pri_P)
    Pout=me%br(2)%p%Prim(Pri_P)
    dPres=Pout-Pin

    if(me%SubType=="pump") then
        mdot=imposed_mass_flow_rate(me,dPres,me%mdot0)
    else
        mdot=imposed_mass_flow_rate(me,dPres,me%mdot0,me%dp0)
    endif

    call mdot_derivatives(me,deriv_wr_both_br)

    call incoming_branch_cold_circulator(me,mdot,dPres,VarIn,deriv_wr_both_br,ssIn,FluxIn)

    call outgoing_branch_cold_circulator(me,mdot,VarIn,deriv_wr_both_br,ssOut,Fluxout)
    if (sim_error > 0) return

    me%wave_time = min(me%br(1)%p%dxLoc/ssIn,me%br(2)%p%dxLoc/ssOut)

    me%br(1)%p%flx(:)=FluxIn%flx(:)
    me%br(1)%p%vit=FluxIn%vit
    me%br(1)%p%derFlx_derCon(:,:)=FluxIn%DerFlxdBr_a(:,:)
    me%br(1)%p%derVit_derCon(:)=FluxIn%DerVitdBr_a(:)

    me%in_DerFlxdBr_out(:,:)=FluxIn%DerFlxdBr_b(:,:)
    me%in_DerVitdBr_out(:)=FluxIn%DerVitdBr_b(:)

    me%br(2)%p%flx(:)=Fluxout%flx(:)
    me%br(2)%p%vit=Fluxout%vit
    me%br(2)%p%derFlx_derCon(:,:)=FluxOut%DerFlxdBr_b(:,:)
    me%br(2)%p%derVit_derCon(:)=FluxOut%DerVitdBr_b(:)

    me%out_DerFlxdBr_in(:,:)=FluxOut%DerFlxdBr_a(:,:)
    me%out_DerVitdBr_in(:)=FluxOut%DerVitdBr_a(:)
end subroutine circulator_resolution_from_and_to_ports    


subroutine links_circulator_update(me,dt)
    real(dp), intent(in) :: dt
    type(circulator_t), intent(inout) :: me
  
    real(dp) :: VarNCcenter
    
    ! Link IN-out
    VarNCcenter=0.0_dp
    if(R_Correction) VarNCcenter=me%br(1)%p%Prim(Pri_ro)*(me%br(1)%p%Prim(Pri_c)**2)
    me%MxInOut(:,:)=me%br(1)%p%Sgn4j*(dt/me%br(1)%p%dxLoc)*me%in_DerFlxdBr_out(:,:)
    me%MxInOut(:,Con_R)=me%MxInOut(:,Con_R)+me%br(1)%p%Sgn4j*(dt/me%br(1)%p%dxLoc)*VarNCcenter*me%in_DerVitdBr_out(:)
  
    ! Link OUT-in
    VarNCcenter=0.0_dp
    if(R_Correction) VarNCcenter=me%br(2)%p%Prim(Pri_ro)*(me%br(2)%p%Prim(Pri_c)**2)
    me%MxOutIn(:,:)=me%br(2)%p%Sgn4j*(dt/me%br(2)%p%dxLoc)*me%out_DerFlxdBr_in(:,:)
    me%MxOutIn(:,Con_R)=me%MxOutIn(:,Con_R)+me%br(2)%p%Sgn4j*(dt/me%br(2)%p%dxLoc)*VarNCcenter*me%out_DerVitdBr_in(:)    
end subroutine links_circulator_update    

end module cmp_circulator_calc_m