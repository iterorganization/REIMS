! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_circulator_init_m
    use krn_interface_m
    use krn_simulation_m
    use lib_input_m, only: input_t
    use krn_global_tools_m
    use lib_ext_math_m
    use lib_He_thermo_m
    implicit none

    type, extends(brent_t) :: circulator_dynamic_parameters_t
        real(dp) :: mdot,VarIn,area,Pr,Ur,Zr
        contains
            procedure :: f => circulator_dynamic_parameters_f      
            procedure :: df => circulator_dynamic_parameters_df
    end type circulator_dynamic_parameters_t 

    type FF_flux_port_pointer_CC_t
        type(FF_flux_port_t), pointer :: p
    end type FF_flux_port_pointer_CC_t    

    type circulator_t
        character(:), allocatable :: SubType
        type(FF_flux_port_pointer_CC_t) :: br(2)

        real(dp) :: MxInOut(Nb_VarC,Nb_VarC)
        real(dp) :: MxOutIn(Nb_VarC,Nb_VarC)
       
        real(dp) :: mdot0,dp0
        real(dp) :: a_pr(Nb_VarP),b_pr(Nb_VarP)
        real(dp) :: a_cs(Nb_VarC),b_cs(Nb_VarC)
        real(dp) :: out_DerFlxdBr_in(Nb_VarC,Nb_VarC),in_DerFlxdBr_out(Nb_VarC,Nb_VarC)
        real(dp) :: out_DerVitdBr_in(Nb_VarC),in_DerVitdBr_out(Nb_VarC)
        real(dp) :: wave_time !! Minimum time for wave propagation
        type(circulator_dynamic_parameters_t) :: dynCC
    end type circulator_t

contains

subroutine circulator_init_part1(me,krn,cfg)
    type(circulator_t),   intent(out)   :: me
    type(krn_t),          intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg

    integer :: nb_non_zeros_exp, nb_non_zeros_imp

    me%SubType = cfg%path
    me%mdot0 = cfg%dbl('m0')

    if(me%SubType=="compressor") me%dp0 = cfg%dbl('dp0')
     
    nb_non_zeros_exp=0
    nb_non_zeros_imp=2*Nb_VarC*Nb_VarC

    call krn%add('',0,nb_non_zeros_exp,nb_non_zeros_imp,0,0,0,0,0)
    
end subroutine circulator_init_part1

subroutine circulator_init_part2(me,krn,cfg)
    type(circulator_t),   intent(inout) :: me
    type(krn_t), target,  intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg

    type(input_t), allocatable :: branches(:)

    branches = cfg%dict1d('link')

    me%br(1)%p => krn%FF_flux_ports(krn%FF_flux_list%find(branches(1)%str('id'),&
                                    convert2int_port(branches(1)%str('node'))))
    if(me%br(1)%p%connected) then
        print*,'Port ',branches(1)%int('node'),' belonging to ',branches(1)%str('id'),' already connected to a link'
        stop
    endif
    if(convert2int_port(branches(1)%str('node'))/=2) then
        print*,'The node of the 1st link must be out'
        stop
    endif
    me%br(1)%p%connected=.true.

    me%br(2)%p => krn%FF_flux_ports(krn%FF_flux_list%find(branches(2)%str('id'),&
                                    convert2int_port(branches(2)%str('node'))))
    if(me%br(2)%p%connected) then
        print*,'Port ',branches(2)%int('node'),' belonging to ',branches(2)%str('id'),' already connected to a link'
        stop
    endif
    me%br(2)%p%connected=.true.    

    call krn%coo_add_link(me%br(1)%p, me%br(2)%p, idx_4x4_col, idx_4x4_row, to_ptr_arr(me%MxInOut))
    call krn%coo_add_link(me%br(2)%p, me%br(1)%p, idx_4x4_col, idx_4x4_row, to_ptr_arr(me%MxOutin))
end subroutine circulator_init_part2

function circulator_dynamic_parameters_f(me,x)
    class(circulator_dynamic_parameters_t), intent(in) :: me
    real(dp),                               intent(in) :: x
        
    real(dp) :: circulator_dynamic_parameters_f    
    real(dp) :: ustar, rostarL, TstarL, R_Corr_StarL, estarL, CstarL

    Ustar=me%Ur+(x-me%Pr)/me%Zr
    rostarL=me%mdot/(Ustar*me%area)
    call state_roP(rostarL, x, R_Corr_StarL, estarL, TstarL, CstarL)
    circulator_dynamic_parameters_f=-(estarL+x/rostarL+0.5_dp*(ustar**2)-me%VarIn)/1.0e3_dp
end function circulator_dynamic_parameters_f
    
function circulator_dynamic_parameters_df(dynamic,x)
    class(circulator_dynamic_parameters_t), intent(inout) :: dynamic
    real(dp), intent(inout) :: x

    real(dp) :: circulator_dynamic_parameters_df
    real(dp) :: TstarL,rostarL,dTdp_Ro,dTdRo_p,dRodP,dUdP,Ustar,dRstardRo,dRstardT,dRdP,Rstar
    real(dp) :: cv_local, dPdT_local

    Ustar=dynamic%Ur+(x-dynamic%Pr)/dynamic%Zr
    rostarL=dynamic%mdot/(Ustar*dynamic%area)
    TstarL=T_roP(rostarL,x)
    call jacobian_roT(rostarL, TstarL, cv_local, dPdT_local, dTdp_Ro, dTdRo_p, dRstardRo, dRstardT, Rstar)
    dUdP=1.0_dp/dynamic%Zr
    dRodP=-(dynamic%mdot/dynamic%area)*dUdP/(Ustar**2)
    dRdP=(dRstardRo+dRstardT*dTdRo_p)*dRodP+dRstardT*dTdp_Ro
    circulator_dynamic_parameters_df=-((dRdP*rostarL-Rstar*dRodP)/(rostarL**2)+Ustar*dUdP)/1.0e3_dp
end function circulator_dynamic_parameters_df

end module cmp_circulator_init_m