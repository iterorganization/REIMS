! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_strand_flux_m
    use cmp_strand_init_m
    implicit none

contains

subroutine S_heat_flux_internal_calculation(me,i)
     
    type(strand_t), intent(inout) :: me
    integer, intent(in) :: i

    real(dp) :: S_TempL,S_TempR,lambdSL,lambdSR,DerlambdSLdTL,DerlambdSRdTR
    real(dp) :: dxL,dxR,Tint,DerTintL,DerTintR,num,den
    integer :: ii

    S_TempL=me%StVar%SCtemp(i)
    S_TempR=me%StVar%SCtemp(i+1)

    lambdSL=0.0_dp
    DerlambdSLdTL=0.0_dp
    lambdSR=0.0_dp
    DerlambdSRdTR=0.0_dp

    lambdSL=lambdSL+me%ct(STAB)%area*me%mat_stab%thermal_conductivity(S_TempL,me%scen%Bfield(i))
    lambdSL=lambdSL+me%ct(SUPC)%area*me%mat_supc%thermal_conductivity(S_TempL)
    
    DerlambdSLdTL=DerlambdSLdTL+me%ct(STAB)%area*me%mat_stab%thermal_conductivity_der(S_TempL,me%scen%Bfield(i))
    DerlambdSLdTL=DerlambdSLdTL+me%ct(SUPC)%area*me%mat_supc%thermal_conductivity_der(S_TempL)

    lambdSR=lambdSR+me%ct(STAB)%area*me%mat_stab%thermal_conductivity(S_TempR,me%scen%Bfield(i+1))
    lambdSR=lambdSR+me%ct(SUPC)%area*me%mat_supc%thermal_conductivity(S_TempR)
    
    DerlambdSRdTR=DerlambdSRdTR+me%ct(STAB)%area*me%mat_stab%thermal_conductivity_der(S_TempR,me%scen%Bfield(i+1))
    DerlambdSRdTR=DerlambdSRdTR+me%ct(SUPC)%area*me%mat_supc%thermal_conductivity_der(S_TempR)

    lambdSL=lambdSL/sum(me%ct(:)%area)
    DerlambdSLdTL=DerlambdSLdTL/sum(me%ct(:)%area)
    lambdSR=lambdSR/sum(me%ct(:)%area)
    DerlambdSRdTR=DerlambdSRdTR/sum(me%ct(:)%area)

    dxL=me%SC_Prop%dxLoc(i)
    dxR=me%SC_Prop%dxLoc(i+1)
    Tint=(lambdSL*dxR*S_TempL+lambdSR*dxL*S_TempR)/(lambdSL*dxR+lambdSR*dxL)
    num=lambdSL*dxR*S_TempL+lambdSR*dxL*S_TempR
    den=lambdSL*dxR+lambdSR*dxL
    DerTintL=dxR*((DerlambdSLdTL*S_TempL+lambdSL)*den-num*DerlambdSLdTL)/(den**2)
    DerTintR=dxL*((DerlambdSRdTR*S_TempR+lambdSR)*den-num*DerlambdSRdTR)/(den**2)

    me%flxS%wire(i)=-2.0_dp*lambdSL*(Tint-S_TempL)/dxL
    me%flxS%wire_DerdTL(i)=-(2.0_dp/dxL)*(DerlambdSLdTL*(Tint-S_TempL)+lambdSL*(DerTintL-1.0_dp))
    me%flxS%wire_DerdTR(i)=-(2.0_dp/dxL)*lambdSL*DerTintR

end subroutine S_heat_flux_internal_calculation



subroutine heat_diffusion_strand(me)
    type(strand_t), intent(inout) :: me
    
    integer :: i
    
    do i=1,me%SC_Prop%NbCells-1
      call S_heat_flux_internal_calculation(me,i)
    enddo
        
end subroutine heat_diffusion_strand 

end module cmp_strand_flux_m