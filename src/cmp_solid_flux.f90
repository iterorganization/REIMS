! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_solid_flux_m
    use cmp_solid_init_m
    implicit none

contains

subroutine S_heat_flux_internal_calculation(me,i)
     
    type(solid_t), intent(inout) :: me
    integer, intent(in) :: i

    real(dp) :: S_TempL,S_TempR,lambdSL,lambdSR,DerlambdSLdTL,DerlambdSRdTR
    real(dp) :: dxL,dxR,Tint,DerTintL,DerTintR,num,den
    integer :: ii

    S_TempL=me%StVar%MCtemp(i)
    S_TempR=me%StVar%MCtemp(i+1)

    lambdSL=me%mat%thermal_conductivity(S_TempL)
    DerlambdSLdTL=me%mat%thermal_conductivity_der(S_TempL)
    lambdSR=me%mat%thermal_conductivity(S_TempR)
    DerlambdSRdTR=me%mat%thermal_conductivity_der(S_TempR)  

    dxL=me%MC_Prop%dxLoc(i)
    dxR=me%MC_Prop%dxLoc(i+1)
    Tint=(lambdSL*dxR*S_TempL+lambdSR*dxL*S_TempR)/(lambdSL*dxR+lambdSR*dxL)
    num=lambdSL*dxR*S_TempL+lambdSR*dxL*S_TempR
    den=lambdSL*dxR+lambdSR*dxL
    DerTintL=dxR*((DerlambdSLdTL*S_TempL+lambdSL)*den-num*DerlambdSLdTL)/(den**2)
    DerTintR=dxL*((DerlambdSRdTR*S_TempR+lambdSR)*den-num*DerlambdSRdTR)/(den**2)

    me%flxS%sd(i)=-2.0_dp*lambdSL*(Tint-S_TempL)/dxL
    me%flxS%sd_DerdTL(i)=-(2.0_dp/dxL)*(DerlambdSLdTL*(Tint-S_TempL)+lambdSL*(DerTintL-1.0_dp))
    me%flxS%sd_DerdTR(i)=-(2.0_dp/dxL)*lambdSL*DerTintR

end subroutine S_heat_flux_internal_calculation



subroutine heat_diffusion_solid(me)
    type(solid_t), intent(inout) :: me
    
    integer :: i
    
    if(me%MC_Prop%cond_btw_nodes) then
      do i=1,me%MC_Prop%NbCells-1
        call S_heat_flux_internal_calculation(me,i)
      enddo
  
      ! Adiabatic boundary conditions
      me%flxS%sd(0)=0.0_dp
      me%flxS%sd_DerdTR(0)=0.0_dp
  
      me%flxS%sd(me%MC_Prop%NbCells)=0.0_dp
      me%flxS%sd_DerdTL(me%MC_Prop%NbCells)=0.0_dp
    endif
        
end subroutine heat_diffusion_solid 

end module cmp_solid_flux_m