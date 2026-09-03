! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

program reims_p
    use krn_linear_system_m
    use krn_global_tools_m
    use krn_simulation_m
    use lib_input_m, only: input_t
    use lib_material_m, only: load_materials
    use cmp_channel_flux_m
    use cmp_channel_calc_m
    use cmp_strand_flux_m
    use cmp_strand_calc_m
    use cmp_solid_flux_m
    use cmp_solid_calc_m
    use cmp_mesh2D_flux_m
    use cmp_mesh2D_calc_m
    use cmp_mesh2D_hdf5_write_m
    use cmp_junction_calc_m
    use cmp_circulator_calc_m
    use cmp_boundary_calc_m
    use cmp_FSlink_calc_m
    use cmp_SSfluxLink_m
    use cmp_SSsrcLink_calc_m
    use cmp_FFsrcLink_calc_m
    implicit none

    class(input_t), pointer :: channels_input(:), junctions_input(:), circulators_input(:)
    class(input_t), pointer :: boundaries_input(:), FSlinks_input(:), strands_input(:), mesh2Ds_input(:)
    class(input_t), pointer :: SSfluxLinks_input(:), SSsrcLinks_input(:), FFsrcLinks_input(:)
    class(input_t), pointer :: solids_input(:)
    type(channel_t),       allocatable :: channels(:)
    type(strand_t),     allocatable :: strands(:)
    type(solid_t),      allocatable :: solids(:)
    type(mesh2D_t),     allocatable :: mesh2Ds(:)
    type(junction_t),   allocatable :: junctions(:)
    type(circulator_t), allocatable :: circulators(:)
    type(boundary_t),   allocatable :: boundaries(:)
    type(FSlink_t),     allocatable :: FSlinks(:)
    type(FFsrcLink_t),  allocatable :: FFsrcLinks(:)
    type(SSfluxLink_t), allocatable :: SSfluxLinks(:)
    type(SSsrcLink_t),  allocatable :: SSsrcLinks(:)
    type(simulation_t) :: sim
    type(timer_t(20)) :: timer
    type(input_t) :: file
    type(krn_t) :: krn

    integer :: i, error
    integer(i_64) :: step2D
    real(dp) :: next2DWriteTime
    type(str_ptr), allocatable :: external_libs(:)    

    ! Variable declaration for some checks
    ! real(dp) :: FullEnergy_MC, FullEnergy_2D, mDot, HHH, mDotHOut, mDotHin, checkMDotHDt, MassSS_MC, MassSS_2D
    ! integer :: ii
    ! real(dp) :: Qext_Check,Qext_Check_step
    ! real(dp) :: Qext_Check_case,Qext_Check_case_step,SurfReg,Qext_Check_case_loc
    ! integer :: nb_e
    ! real(dp) :: Q_MC, Q_MC1, Q_MC2, Q_2D

    call timer%set(1,description='REIMS start')
    print*,'                                                                   '
    print*,'           _____________________________________________________   '
    print*,'                                                                   '
    print*,'                                    REIMS                          '
    print*,'                (Riemann Explicit Implicit Magnet Simulator)       '
    print*,'                                    v2.1a1                         '
    print*,'           _____________________________________________________   '
    print*,'                                                                   '
    print*,'                                                                   '

    call file%open(reims_config_file())
    print*, 'Read config file completed'
    call sim%init(file%dict('simulation'))

    external_libs=file%str1d('external_libs',empty=.true.)
    allocate(lib_handles(size(external_libs)))
    do i = 1, size(external_libs)
        lib_handles(i) = LoadLibraryA(external_libs(i)%p // '.dll' // c_null_char)
        if (.not. c_associated(lib_handles(i))) error stop 'one of the external libraries is not associated'
    enddo
    call load_materials(file%dict1d('materials',empty=.true.))
    call load_frictions(file%dict1d('friction_correlations',empty=.true.))
    call load_nusselts(file%dict1d('nusselt_correlations',empty=.true.))
    
    channels_input    => file%cmp1d(['channel'])
    junctions_input   => file%cmp1d(['junction'])
    circulators_input => file%cmp1d(['pump','compressor'])
    boundaries_input  => file%cmp1d(['boundary_PT','boundary_mT'])
    FSlinks_input     => file%cmp1d(['thermalink'])
    strands_input     => file%cmp1d(['strand'])
    solids_input      => file%cmp1d(['solid'])
    mesh2Ds_input     => file%cmp1d(['mesh2D'])
    FFsrcLinks_input  => file%cmp1d(['fluidlink'])
    SSfluxLinks_input => file%cmp1d(['strandlink'])
    SSsrcLinks_input  => file%cmp1d(['solidlink'])
    call krn%init_part1(file%nb_cmp())
    call h5%add_table('channel',['P','rho','u','m','T','mH','R','Su'])
    call h5%add_table('strand',['T','B','st','Tcs','marg','Joule'])
    call h5%add_table('solid',['T'])

    allocate(channels(size(channels_input)))
    allocate(strands(size(strands_input)))
    allocate(solids(size(solids_input)))
    allocate(mesh2Ds(size(mesh2Ds_input)))
    allocate(junctions(size(junctions_input)))
    allocate(circulators(size(circulators_input)))
    allocate(boundaries(size(boundaries_input)))
    allocate(FSlinks(size(FSlinks_input)))
    allocate(FFsrcLinks(size(FFsrcLinks_input)))
    allocate(SSfluxLinks(size(SSfluxLinks_input)))
    allocate(SSsrcLinks(size(SSsrcLinks_input)))

    do i = 1, size(channels)
        call channel_init_part1(channels(i), krn, channels_input(i), h5, sim)
    enddo
    do i = 1, size(strands)
        call strand_init_part1(strands(i), krn, strands_input(i), h5, sim)
    enddo
    do i = 1, size(solids)
        call solid_init_part1(solids(i), krn, solids_input(i), h5, sim)
    enddo
    do i = 1, size(mesh2Ds)
        call mesh2D_init_part1(mesh2Ds(i), krn, mesh2Ds_input(i), sim)
    enddo
    do i = 1, size(junctions)
        call junction_init_part1(junctions(i), krn, junctions_input(i))
    enddo
    do i = 1, size(circulators)
        call circulator_init_part1(circulators(i), krn, circulators_input(i))
    enddo
    do i = 1, size(boundaries)
        call boundary_init_part1(boundaries(i), krn, boundaries_input(i), sim)
    enddo
    do i = 1, size(FSlinks)
        call FSlink_init_part1(FSlinks(i), krn, FSlinks_input(i))
    enddo
    do i = 1, size(FFsrcLinks)
        call FFsrcLink_init_part1(FFsrcLinks(i), krn, FFsrcLinks_input(i))
    enddo
    do i = 1, size(SSfluxLinks)
        call SSfluxLink_init_part1(SSfluxLinks(i), krn, SSfluxLinks_input(i))
    enddo
    do i = 1, size(SSsrcLinks)
        call SSsrcLink_init_part1(SSsrcLinks(i), krn, SSsrcLinks_input(i))
    enddo

    call h5%prepare_data()
    call krn%init_part2()
    do i = 1, size(channels)
        call channel_init_part2(channels(i), krn, channels_input(i))
    enddo
    do i = 1, size(strands)
        call strand_init_part2(strands(i), krn, strands_input(i))
    enddo
    do i = 1, size(solids)
        call solid_init_part2(solids(i), krn, solids_input(i))
    enddo
    do i = 1, size(mesh2Ds)
        call mesh2D_init_part2(mesh2Ds(i), krn, mesh2Ds_input(i))
    enddo
    do i = 1, size(junctions)
        call junction_init_part2(junctions(i), krn, junctions_input(i))
    enddo
    do i = 1, size(circulators)
        call circulator_init_part2(circulators(i), krn, circulators_input(i))
    enddo
    do i = 1, size(boundaries)
        call boundary_init_part2(boundaries(i), krn, boundaries_input(i))
    enddo
    do i = 1, size(FSlinks)
        call FSlink_init_part2(FSlinks(i), krn, FSlinks_input(i))
    enddo
    do i = 1, size(FFsrcLinks)
        call FFsrcLink_init_part2(FFsrcLinks(i), krn, FFsrcLinks_input(i))
    enddo
    do i = 1, size(SSfluxLinks)
        call SSfluxLink_init_part2(SSfluxLinks(i), krn, SSfluxLinks_input(i))
    enddo
    do i = 1, size(SSsrcLinks)
        call SSsrcLink_init_part2(SSsrcLinks(i), krn, SSsrcLinks_input(i))
    enddo
    call consistency_ports(channels)
    call krn%coo_to_csr()

    call h5%init(file%dict('write_results'))
    step2D=0
    next2DWriteTime = h5%time_btw_2D_writes
    do i = 1, size(mesh2Ds)
        if(h5%write_results2D) call writing_HDF5_2D_casing(mesh2Ds(i),sim%t,step2D)
    enddo

    !Heat balance check
    ! open(unit=132,file='CheckBalance.out',status='replace')
    ! checkMDotHDt=0.0_dp
    ! MassSS_MC=0.0_dp
    ! do i=1,size(solids)
    !     do ii=1,solids(i)%MC_Prop%NbCells
    !         MassSS_MC = MassSS_MC + solids(i)%MC_Prop%ro_M * solids(i)%MC_Prop%volLoc(ii)
    !     enddo
    ! enddo
    ! MassSS_2D=0.0_dp
    ! do i=1,size(mesh2Ds)
    !     do ii=1,mesh2Ds(i)%M2D_Prop%nb_elements
    !         MassSS_2D = MassSS_2D + 7900.0_dp * mesh2Ds(i)%M2D_Prop%elem(ii)%surface * 0.5_dp
    !     enddo
    ! enddo
    ! ! Energy to be lost by the  solid part to be cooled down from 10K to 4.3K
    ! FullEnergy_MC=MassSS_MC*(CpSS_Integral(4.3_dp)-CpSS_Integral(10.0_dp))
    ! FullEnergy_2D=MassSS_2D*(CpSS_Integral(4.3_dp)-CpSS_Integral(6.0_dp))

    ! Heat balance check
    ! open(unit=134,file='CheckBalanceSS.out',status='replace')
    ! Q_MC=0.0_dp
    ! Q_MC1=0.0_dp
    ! Q_MC2=0.0_dp
    ! Q_2D=0.0_dp

    ! Load check
    ! open(unit=133,file='CheckHeatLoad.out',status='replace')
    ! Qext_Check=0.0_dp
    ! Qext_Check_case=0.0_dp    
    ! print*, 'Volume of solid:',sum(chunks%MC_Prop%VolM * chunks%MC_Prop%NbCells)
    call timer%set(2,1,'initialisation time')

    ! ---------------------------- M A I N   L O O P ---------------------------
    do while(sim%t < sim%t_final - sim_delta); call main_loop(); enddo
    ! --------------------------------------------------------------------------

    ! write(132,*)
    ! flush(132)
    ! close(132)

    ! write(133,*)
    ! flush(133)
    ! close(133)

    ! write(134,*)
    ! flush(134)
    ! close(134)

    call close_hdf5_files()
    print*,''
    print*,'   -----------------------------------'
    call timer%set(3,1,'total computation time')
    call timer%print()
    print*,'   -----------------------------------'
    print*,''

    print*,''
    print*, 'D. Furfaro, J. Kosek, A. Ovcharov, T. Schioler, R. Rotella, T. Luce,'
    print*, 'A new fast and robust thermo-hydraulic code for ITER superconducting'
    print*, 'magnet simulation, Cryogenics, Volume 144, 2024, 103978'
    print*,''

contains

subroutine main_loop() !! Main loop
    real(dp) :: error_val

    call timer%set(9,description='Main loop start')
    if(sim%rollback) then ! Rollback in progress: state already set, skip state advance
        sim%rollback = .false.
        h5%steps_to_revert = 3
    else if(sim%rejected) then ! Retrieve previous state
        channels%StVar       = channels%StVarOld
        strands%StVar     = strands%StVarOld
        solids%StVar      = solids%StVarOld
        mesh2Ds%StVar     = mesh2Ds%StVarOld
    else                  ! Advance state variables
        channels%StVarOld6   = channels%StVarOld5
        channels%StVarOld5   = channels%StVarOld4
        channels%StVarOld4   = channels%StVarOld3
        channels%StVarOld3   = channels%StVarOld2
        channels%StVarOld2   = channels%StVarOld
        channels%StVarOld    = channels%StVar
        strands%StVarOld6 = strands%StVarOld5
        strands%StVarOld5 = strands%StVarOld4
        strands%StVarOld4 = strands%StVarOld3
        strands%StVarOld3 = strands%StVarOld2
        strands%StVarOld2 = strands%StVarOld
        strands%StVarOld  = strands%StVar
        solids%StVarOld6  = solids%StVarOld5
        solids%StVarOld5  = solids%StVarOld4
        solids%StVarOld4  = solids%StVarOld3
        solids%StVarOld3  = solids%StVarOld2
        solids%StVarOld2  = solids%StVarOld
        solids%StVarOld   = solids%StVar
        mesh2Ds%StVarOld6 = mesh2Ds%StVarOld5
        mesh2Ds%StVarOld5 = mesh2Ds%StVarOld4
        mesh2Ds%StVarOld4 = mesh2Ds%StVarOld3
        mesh2Ds%StVarOld3 = mesh2Ds%StVarOld2
        mesh2Ds%StVarOld2 = mesh2Ds%StVarOld
        mesh2Ds%StVarOld  = mesh2Ds%StVar
    endif

    ! Qext_Check_step=0.0_dp

    !$omp parallel do private(i)
        do i=1, size(channels)
            call friction_source_term(channels(i))
            call channels_FF_flux_port_comm(channels(i))
        enddo
    !$omp end parallel do
    if (handle_numerical_error()) return
    !$omp parallel do private(i)
        do i=1, size(junctions)
            call junction_resolution_from_and_to_ports(junctions(i))
        enddo
    !$omp end parallel do
    if (handle_numerical_error()) return
    do i=1, size(circulators)
        call circulator_resolution_from_and_to_ports(circulators(i))
    enddo
    if (handle_numerical_error()) return
    do i=1, size(boundaries)
        call boundary_resolution_from_and_to_ports(boundaries(i))
    enddo
    if (handle_numerical_error()) return
    !$omp parallel do private(i)
        do i=1, size(channels)
            call He_Riemann_solver_channel(channels(i),sim)
            call flux_from_FF_flux_ports_self(channels(i))
        enddo
    !$omp end parallel do
    if (handle_numerical_error()) return
    !$omp parallel do private(i)
        do i=1, size(strands)
            call scenario_update(strands(i),sim)
            call strands_SS_flux_port_comm(strands(i))
        enddo
    !$omp end parallel do
    if (handle_numerical_error()) return
    !$omp parallel do private(i)
        do i=1, size(solids)
            call solid_scenario_update(solids(i))
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(SSfluxLinks)
            call SSfluxLink_resolution_from_and_to_ports(SSfluxLinks(i))
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(strands)
            call heat_diffusion_strand(strands(i))
            call flux_from_SS_flux_ports_self(strands(i))
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(solids)
            call heat_diffusion_solid(solids(i))
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(mesh2Ds)
            call heat_diffusion_mesh2D(mesh2Ds(i),sim)
            call mesh2Ds_FS_port_comm(mesh2Ds(i))
            call mesh2Ds_SS_src_port_comm(mesh2Ds(i))
        enddo
    !$omp end parallel do

    ! Calculation of minimal time for wave propagation through the node
    ! and updating explicit step
    call sim%update_exp_dt( min( &
        minVal(channels%wave_time), &
        minVal(junctions%wave_time), &
        minVal(circulators%wave_time), &
        minVal(boundaries%wave_time)))

    call timer%set(10,9,'Fluxes + Riemann')
    call timer%set(11)
    if(sim%explicit) then
        ! explicit fluxes / implicit source terms ------------------------------
        !print*,'Explicit scheme at physical time',sim%t,'(in s)'

        !$omp parallel do private(i)
            do i=1,size(channels)
                call explicit_scheme_for_fluxes_channel(channels(i),sim%dt)
                call friction_source_term(channels(i))
                call source_term_definition_channel(channels(i),sim)
            enddo
        !$omp end parallel do
        if (handle_numerical_error()) return
        !$omp parallel do private(i)
            do i=1,size(strands)
                call explicit_scheme_for_fluxes_strand(strands(i),sim%dt)
                call source_term_definition_strand(strands(i),sim) !,Qext_Check_step)   ! Attention in case of parallelization --> Qext_Check
            enddo
        !$omp end parallel do
        if (handle_numerical_error()) return
        !$omp parallel do private(i)
            do i=1,size(solids)
                call explicit_scheme_for_fluxes_solid(solids(i),sim%dt)
                call source_term_definition_solid(solids(i),sim)
            enddo
        !$omp end parallel do
    else ! full implicit -------------------------------------------------------
        !$omp parallel do private(i)
            do i=1,size(channels)
                call full_physics_definition_channel(channels(i),sim)
            enddo
        !$omp end parallel do
        if (handle_numerical_error()) return
        !$omp parallel do private(i)
            do i=1,size(strands)
                call full_physics_definition_strand(strands(i),sim) !,Qext_Check_step)
            enddo
        !$omp end parallel do
        if (handle_numerical_error()) return  
        !$omp parallel do private(i)                      
            do i=1,size(solids)
                call full_physics_definition_solid(solids(i),sim)
            enddo
        !$omp end parallel do
        !$omp parallel do private(i)
            do i=1,size(junctions)
                call links_junction_update(junctions(i),sim%dt)
            enddo
        !$omp end parallel do
        do i=1,size(circulators)
            call links_circulator_update(circulators(i),sim%dt)
        enddo
        !$omp parallel do private(i)
            do i=1,size(SSfluxLinks)
                call links_SSfluxLink_update(SSfluxLinks(i),sim%dt)
            enddo
        !$omp end parallel do
    endif ! --------------------------------------------------------------------

    ! Qext_Check_case_step=0.0_dp
    !$omp parallel do private(i)
        do i=1,size(mesh2Ds)
            ! SurfReg=0.0_dp
            ! do nb_e=1,mesh2Ds(i)%M2D_Prop%nb_elements
            !     if(trim(mesh2Ds(i)%M2D_Prop%elem(nb_e)%PhysE)=='SS_Outer') then
            !         SurfReg=SurfReg+mesh2Ds(i)%M2D_Prop%elem(nb_e)%surface
            !     endif
            ! enddo
            call full_physics_definition_mesh2D(mesh2Ds(i),sim) !,Qext_Check_case_loc)       
            ! Qext_Check_case_step = Qext_Check_case_step + Qext_Check_case_loc * SurfReg * 1.062_dp
        enddo
    !$omp end parallel do

    !$omp parallel do private(i)
        do i=1, size(FSlinks)
            call FSlink_src_reinitialization(FSlinks(i))
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(FSlinks)
            call FSlink_resolution_from_and_to_ports(FSlinks(i))
            call links_FSlink_update(FSlinks(i),sim%dt)
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(channels)
            call source_from_FS_ports_self_channel(channels(i),sim%dt)
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(strands)
            call source_from_FS_ports_self_strand(strands(i),sim%dt)
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(solids)
            call source_from_FS_ports_self_solid(solids(i),sim%dt)
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(mesh2Ds)
            call source_from_FS_ports_self_mesh2D(mesh2Ds(i),sim%dt)
        enddo
    !$omp end parallel do

    !$omp parallel do private(i)
        do i=1, size(SSsrcLinks)
            call SSsrcLink_src_reinitialization(SSsrcLinks(i))
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(SSsrcLinks)
            call SSsrcLink_resolution_from_and_to_ports(SSsrcLinks(i))
            call links_SSsrcLink_update(SSsrcLinks(i),sim%dt)
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(solids)
            call source_from_SS_src_ports_self_solid(solids(i),sim%dt)
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(mesh2Ds)
            call source_from_SS_src_ports_self_mesh2D(mesh2Ds(i),sim%dt)
        enddo
    !$omp end parallel do

    !$omp parallel do private(i)
        do i=1, size(FFsrcLinks)
            call FFsrcLink_resolution_from_and_to_ports(FFsrcLinks(i))
            call links_FFsrcLink_update(FFsrcLinks(i),sim%dt)
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(channels)
            call source_from_FF_src_ports_self_channel(channels(i),sim%dt)
        enddo
    !$omp end parallel do
    if (handle_numerical_error()) return

    ! ---------------------- L I N E A R   S O L V E R -------------------------
    call Pardiso_fact_solve(krn,sim%explicit,timer)
    if (handle_numerical_error()) return
    ! --------------------------------------------------------------------------
    call timer%set(12)

    !$omp parallel do private(i)
        do i=1,size(channels)
            call from_sol_to_prim_channel(channels(i))
        enddo
    !$omp end parallel do
    if (handle_numerical_error()) return
    !$omp parallel do private(i)
        do i=1,size(strands)
            call from_sol_to_temp_strand(strands(i),sim)
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1,size(solids)
            call from_sol_to_temp_solid(solids(i),sim)
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1,size(mesh2Ds)
            call from_sol_to_temp_mesh2D(mesh2Ds(i),sim)
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1,size(channels)
            call channels_to_FF_src_port_comm(channels(i))
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1, size(FFsrcLinks)
            call FFsrcLink_Prelax_resolution_from_and_to_ports(FFsrcLinks(i))
        enddo
    !$omp end parallel do
    !$omp parallel do private(i)
        do i=1,size(channels)
            call FF_src_port_to_channels_comm(channels(i))
            call last_tasks_for_channels(channels(i),sim)
        enddo
    !$omp end parallel do

    if (handle_numerical_error()) return

    if(sim%explicit) then
        ! pressure criteria for switching to implicit
        error_val = maxVal(sqrt(channels%err) / max(1.0e-6_dp, sqrt(channels%err_den)))
    else; error_val = &
        ! compute local truncation error (lte) - energy calculation as
        ! criteria for implicit step
        sqrt(sum(channels%err))   / max(sqrt(sum(channels%err_den)),   sim_delta) + &
        sqrt(sum(strands%err)) / max(sqrt(sum(strands%err_den)), sim_delta) + &
        sqrt(sum(solids%err))  / max(sqrt(sum(solids%err_den)),  sim_delta) + &
        sqrt(sum(mesh2Ds%err)) / max(sqrt(sum(mesh2Ds%err_den)), sim_delta)
    endif
    call timer%set(13,12,'Solution distribution')

    call sim%step_management(error_val)
    if(.not.sim%rejected) call h5%write(sim%t,sim%dt,sim%t_final)
    if(h5%write_results2D .and. .not.sim%rejected) then
        if(sim%t >= next2DWriteTime) then
            step2D = step2D + 1
            do i=1,size(mesh2Ds)
                call writing_HDF5_2D_casing(mesh2Ds(i),sim%t,step2D)
            enddo
            next2DWriteTime = next2DWriteTime + h5%time_btw_2D_writes
        endif
    endif

    ! if(.not.sim%rejected) then
    !   mDotHOut = channels(1)%flxHe%Cons(Con_Ene,15) * channels(1)%HeProp%Area
    !   mDotHin  = channels(1)%flxHe%Cons(Con_Ene,5)*channels(1)%HeProp%Area
    !   checkMDotHDt = checkMDotHDt + (mDotHOut-mDotHin)*sim%dt
    !   write(132,'(4(e13.6,1x))') sim%t, FullEnergy_MC, FullEnergy_2D, -checkMDotHDt

    !   Qext_Check = Qext_Check + Qext_Check_step
    !   write(133,'(2(e13.6,1x))') sim%t, Qext_Check

    !   Qext_Check_case = Qext_Check_case + Qext_Check_case_step
    !   write(133,'(2(e13.6,1x))') sim%t, Qext_Check_case

    !   do i=1,size(solids)
    !       do ii=1,solids(i)%MC_Prop%NbCells
    !           if(abs(solids(i)%StVarOld%MCtemp(ii))>1.0e-10_dp) then
    !               Q_MC = Q_MC + solids(i)%MC_Prop%ro_M * solids(i)%MC_Prop%volLoc(ii) * &
    !                             (CpSS_Integral(solids(i)%StVar%MCtemp(ii)) - CpSS_Integral(solids(i)%StVarOld%MCtemp(ii)))
    !               if(i==1) then
    !                   Q_MC1 = Q_MC1 + solids(i)%MC_Prop%ro_M * solids(i)%MC_Prop%volLoc(ii) * &
    !                             (CpSS_Integral(solids(i)%StVar%MCtemp(ii)) - CpSS_Integral(solids(i)%StVarOld%MCtemp(ii)))
    !               else if(i==2) then
    !                   Q_MC2 = Q_MC2 + solids(i)%MC_Prop%ro_M * solids(i)%MC_Prop%volLoc(ii) * &
    !                   (CpSS_Integral(solids(i)%StVar%MCtemp(ii)) - CpSS_Integral(solids(i)%StVarOld%MCtemp(ii)))
    !               endif
    !           endif
    !       enddo
    !   enddo
    !   do i=1,size(mesh2Ds)
    !       do ii=1,mesh2Ds(i)%M2D_Prop%nb_elements
    !           if(abs(mesh2Ds(i)%StVarOld%temp(ii))>1.0e-10_dp) then
    !               Q_2D = Q_2D + 7900.0_dp * mesh2Ds(i)%M2D_Prop%elem(ii)%surface * 0.5_dp * &   ! 0.5_dp for 0.5m of contact with MC
    !                             (CpSS_Integral(mesh2Ds(i)%StVar%temp(ii)) - CpSS_Integral(mesh2Ds(i)%StVarOld%temp(ii)))
    !           endif
    !       enddo
    !   enddo
    !   write(134,'(9(e13.6,1x))') sim%t, Q_MC, Q_MC1, Q_MC2, Q_2D
    ! endif

end subroutine main_loop

subroutine close_hdf5_files()
    integer :: j
    if(h5%write_results2D) then
        do j = 1, size(mesh2Ds)
            call additional_data_HDF5_2D_casing(mesh2Ds(j),step2D)
            call h5fclose_f(mesh2Ds(j)%f_id_writing, error)
        enddo
    endif
    call h5%close()
end subroutine close_hdf5_files

logical function handle_numerical_error()
    integer :: i
    handle_numerical_error = (sim_error > 0)
    if (sim_error == 0) return
    do i = 1, min(sim_error, MAX_SIM_ERRORS)
        print*, trim(sim_error_msg(i))
    enddo
    if (sim%in_recovery) then
        call close_hdf5_files()
        error stop trim(sim_error_msg(1))
    end if
    sim%consecutive_rollbacks = sim%consecutive_rollbacks + 1
    if (sim%consecutive_rollbacks > 3) then
        print*, 'Too many consecutive rollbacks (', sim%consecutive_rollbacks, ')'
        call close_hdf5_files()
        error stop trim(sim_error_msg(1))
    end if
    ! 3-step rollback: restart from n-3 (= Old4)
    ! new Old=n-3, Old2=n-4, Old3=n-5 (BDF2 + LTE history)
    channels%StVar     = channels%StVarOld4
    channels%StVarOld  = channels%StVarOld4
    channels%StVarOld2 = channels%StVarOld5
    channels%StVarOld3 = channels%StVarOld6
    strands%StVar      = strands%StVarOld4
    strands%StVarOld   = strands%StVarOld4
    strands%StVarOld2  = strands%StVarOld5
    strands%StVarOld3  = strands%StVarOld6
    solids%StVar       = solids%StVarOld4
    solids%StVarOld    = solids%StVarOld4
    solids%StVarOld2   = solids%StVarOld5
    solids%StVarOld3   = solids%StVarOld6
    mesh2Ds%StVar      = mesh2Ds%StVarOld4
    mesh2Ds%StVarOld   = mesh2Ds%StVarOld4
    mesh2Ds%StVarOld2  = mesh2Ds%StVarOld5
    mesh2Ds%StVarOld3  = mesh2Ds%StVarOld6
    sim%t          = sim%t - sim%dtPrev1 - sim%dtPrev2 - sim%dtPrev3
    sim%dtPrev1    = sim%dtPrev4
    sim%dtPrev2    = sim%dtPrev5
    sim%explicit            = .true.
    sim%expl_steps          = 3
    sim%rollback            = .true.
    sim%in_recovery         = .true.
    sim%recovery_steps_left = 3
    sim_error               = 0
end function handle_numerical_error

end program reims_p
