! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_SSfluxLink_m
    use krn_interface_m
    use lib_input_m, only: input_t
    implicit none

    type SSfluxLink_t
        integer :: NbSolLk        
        type(SS_flux_port_pointer_t), allocatable :: solLk(:)

        real(dp) :: derFlx_derOth12, derFlx_derOth21
        real(dp) :: Lk12, Lk21
    end type SSfluxLink_t

contains

subroutine all_SSfluxLink_allocation(me)
    type(SSfluxLink_t), intent(inout) :: me

    allocate(me%solLk(me%NbSolLk))

end subroutine all_SSfluxLink_allocation


subroutine SSfluxLink_init_part1(me,krn,cfg)
    type(SSfluxLink_t), intent(out) :: me
    type(krn_t), intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg

    integer :: nb_non_zeros_exp, nb_non_zeros_imp

    nb_non_zeros_exp = 0
    nb_non_zeros_imp = 2

    call krn%add('',0,nb_non_zeros_exp,nb_non_zeros_imp,0,0,0,0,0)
end subroutine SSfluxLink_init_part1


subroutine SSfluxLink_init_part2(me,krn,cfg)
    type(SSfluxLink_t), intent(inout) :: me
    type(krn_t), target, intent(inout) :: krn
    class(input_t), pointer, intent(in)    :: cfg

    type(input_t), allocatable :: links(:)
    integer :: i

    me%NbSolLk=size(cfg%dict1d('link'))
    call all_SSfluxLink_allocation(me)

    links = cfg%dict1d('link')
    do i = 1, me%NbSolLk
        me%solLk(i)%p => krn%SS_flux_ports(krn%SS_flux_list%find(links(i)%str('id'),convert2int_port(links(i)%str('node'))))  
    enddo

    call krn%coo_add_link(me%solLk(1)%p, me%solLk(2)%p, [1], [1], to_ptr_arr(me%Lk12))
    call krn%coo_add_link(me%solLk(2)%p, me%solLk(1)%p, [1], [1], to_ptr_arr(me%Lk21))

end subroutine SSfluxLink_init_part2


subroutine SSfluxLink_resolution_from_and_to_ports(me)
    type(SSfluxLink_t), intent(inout) :: me

    real(dp) :: TempS(me%NbSolLk), lambdS(me%NbSolLk), DerlambdS(me%NbSolLk), AreaMS(me%NbSolLk), sgn(me%NbSolLk), dxLoc(me%NbSolLk)
    real(dp) :: FlxS, DerFlxSdT1, DerFlxSdT2, Tint, num, den, DerTint(me%NbSolLk)
    integer :: NbSolLk, ibr

    NbSolLk=me%NbSolLk
    do ibr=1,NbSolLk
      TempS(ibr)=me%solLk(ibr)%p%TempS
      lambdS(ibr)=me%solLk(ibr)%p%lambdS
      DerlambdS(ibr)=me%solLk(ibr)%p%DerlambdS
      AreaMS(ibr)=me%solLk(ibr)%p%AreaS
      sgn(ibr)=me%solLk(ibr)%p%Sgn4lk 
      dxLoc(ibr)=me%solLk(ibr)%p%dxLoc
    enddo

    Tint=(lambdS(1)*dxLoc(2)*TempS(1)+lambdS(2)*dxLoc(1)*TempS(2))/(lambdS(1)*dxLoc(2)+lambdS(2)*dxLoc(1))
    num=lambdS(1)*dxLoc(2)*TempS(1)+lambdS(2)*dxLoc(1)*TempS(2)
    den=lambdS(1)*dxLoc(2)+lambdS(2)*dxLoc(1)
    DerTint(1)=dxLoc(2)*((DerlambdS(1)*TempS(1)+lambdS(1))*den-num*DerlambdS(1))/(den**2)
    DerTint(2)=dxLoc(1)*((DerlambdS(2)*TempS(2)+lambdS(2))*den-num*DerlambdS(2))/(den**2)

    FlxS=-2.0_dp*lambdS(1)*(Tint-TempS(1))/dxLoc(1)
    DerFlxSdT1=-(2.0_dp/dxLoc(1))*(DerlambdS(1)*(Tint-TempS(1))+lambdS(1)*(DerTint(1)-1.0_dp))
    DerFlxSdT2=-(2.0_dp/dxLoc(1))*lambdS(1)*DerTint(2)

    me%solLk(1)%p%Flx=sgn(1)*FlxS
    me%solLk(1)%p%derFlx_derCon=sgn(1)*DerFlxSdT1
    me%solLk(2)%p%Flx=-sgn(2)*FlxS
    me%solLk(2)%p%derFlx_derCon=-sgn(2)*DerFlxSdT2

    me%derFlx_derOth12=sgn(1)*DerFlxSdT2
    me%derFlx_derOth21=-sgn(2)*DerFlxSdT1
   
end subroutine SSfluxLink_resolution_from_and_to_ports


subroutine links_SSfluxLink_update(me,dt)
    type(SSfluxLink_t), intent(inout) :: me
    real(dp), intent(in) :: dt
    real(dp) :: VarcentreeS_1, VarcentreeS_2

    VarcentreeS_1=1.0_dp/(me%solLk(1)%p%rhoMS*me%solLk(1)%p%cpMs)
    me%Lk12=me%solLk(1)%p%Sgn4lk*(dt/me%solLk(1)%p%dxLoc)*VarcentreeS_1*me%derFlx_derOth12

    VarcentreeS_2=1.0_dp/(me%solLk(2)%p%rhoMS*me%solLk(2)%p%cpMs)
    me%Lk21=me%solLk(2)%p%Sgn4lk*(dt/me%solLk(2)%p%dxLoc)*VarcentreeS_2*me%derFlx_derOth21

end subroutine links_SSfluxLink_update

end module cmp_SSfluxLink_m