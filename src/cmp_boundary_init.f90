! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_boundary_init_m
    use krn_interface_m
    use krn_simulation_m
    use lib_input_m, only: input_t
    use krn_global_tools_m
    use lib_He_thermo_m
    implicit none
 
    type boundary_t
        logical :: is_PT
        logical :: hugoniot_bc  !! exact Rankine-Hugoniot relations at this boundary_PT
        type(signal_t) :: imposed_PorM, imposed_T
        type(FF_flux_port_t), pointer :: port
        real(dp) :: wave_time !! Minimum time for wave propagation
    end type boundary_t

contains


subroutine boundary_init_part1(me,krn,cfg,sim)
    type(boundary_t), intent(out) :: me
    type(krn_t), intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg
    type(simulation_t), intent(inout) :: sim

    call me%imposed_T%init(sim, cfg,'t')
    if(cfg%path=="boundary_PT") then
        me%is_PT = .true.
        call me%imposed_PorM%init(sim, cfg,'p')
        me%hugoniot_bc = cfg%bin('hugoniot_boundary', .false.)
    else
        me%is_PT = .false.
        me%hugoniot_bc = .false.
        call me%imposed_PorM%init(sim, cfg,'mdot')
    endif
    call krn%add('',0,0,0,0,0,0,0,0)
end subroutine boundary_init_part1


subroutine boundary_init_part2(me,krn,cfg)
    type(boundary_t), intent(inout) :: me
    type(krn_t), target, intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg

    me%port => krn%FF_flux_ports(krn%FF_flux_list%find(cfg%str('link/id'), &
                                                       convert2int_port(cfg%str('link/node'))))

    if(me%port%connected) then
        print*,'Port ',cfg%str('link/node'),' belonging to ',cfg%str('link/id'),' already connected to a link'
        stop
    endif
    me%port%connected=.true.                                                       
end subroutine boundary_init_part2

end module cmp_boundary_init_m