! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_junction_init_m
    use krn_interface_m
    use krn_simulation_m
    use lib_input_m, only: input_t
    use krn_global_tools_m
    use lib_ext_math_m
    use lib_He_thermo_m
    implicit none

    type :: junction_dynamic_per_branch_t
        !! Parameters stored per branch used in fSolve
        !! I've put large values just for test 
        real(dp) :: rho0, vit0, p0, ETot0, rCorr0
        real(dp) :: SpeedS0 , CSound0, T0, AreaJGB, ThetaJGB
        real(dp) :: SgnJGB, U1 , U2, U3, U4, q0
        real(dp) :: qBar, Frc = 0, DerFric1 = 0, DerFric2 = 0, DerFric3 = 0
        real(dp) :: DerFric4 = 0, dx0, FrcR, DerFricR1
        real(dp) :: DerFricR2, DerFricR3, DerFricR4
    end type junction_dynamic_per_branch_t
    type :: junction_dynamic_star_t
        !! Used locally by fSolve function as variable noted with '*' and '**'
        real(dp) :: rho, ETot, QDM, Energy, vit, rCorr, q, Enthalpy
        real(dp) :: rho_int, ETot_int, QDM_int, Energy_int
    end type junction_dynamic_star_t

    type, extends(fSolve_t) :: junction_dynamic_parameters_t
        !! Object given to non-linear equation solver (fSolve)
        integer :: NbTotBr, NbIn, NbOut, idxKappaRef
        real(dp) :: kappa
        type(junction_dynamic_per_branch_t), allocatable :: br(:)
        real(dp), allocatable :: p_star_prev(:)
        real(dp), allocatable :: alpha_prev(:)
    contains
        procedure :: fcn => junction_dynamic_parameters_fcn
            !! Function given to non-linear equation solver (fSolve)
    end type junction_dynamic_parameters_t

    type FF_flux_port_pointer_t
        type(FF_flux_port_t), pointer :: p
        real(dp) :: angle
    end type FF_flux_port_pointer_t

    type flux_and_derivatives_t
        real(dp) :: Flx(Nb_VarC)
        real(dp), allocatable :: DerFlx(:,:,:)
        real(dp) :: Vit 
        real(dp), allocatable :: DerVit(:,:)
    end type flux_and_derivatives_t    

    type junction_t
        integer :: NbTotBr
        
        type(FF_flux_port_pointer_t), allocatable :: br(:)

        real(dp), allocatable :: DerFlx_dBr(:,:,:,:)
        real(dp), allocatable :: DerVel_dBr(:,:,:)

        type(flux_and_derivatives_t), allocatable :: FlJ(:)
        real(dp), allocatable :: Matrix(:,:,:,:)

        real(dp) :: wave_time !! Minimum time for wave propagation
        type(junction_dynamic_parameters_t) :: dynamic

        real(dp) :: kappa
        logical, allocatable :: idxKappaRef(:)
    end type junction_t

    real(dp), parameter :: epsCoef = 1.0e-8_dp

contains

subroutine all_junction_allocation(me)
    type(junction_t), intent(inout) :: me
    integer :: j

    allocate(me%br(me%NbTotBr))
    allocate(me%FlJ(me%NbTotBr))
    allocate(me%Matrix(Nb_VarC,Nb_VarC,me%NbTotBr,me%NbTotBr))
    me%Matrix = 0.0_dp
    allocate(me%DerFlx_dBr(Nb_VarC,Nb_VarC,me%NbTotBr,me%NbTotBr))
    me%DerFlx_dBr = 0.0_dp
    allocate(me%idxKappaRef(me%NbTotBr))
    me%idxKappaRef = .false.; me%idxKappaRef(1) = .true.
    allocate(me%DerVel_dBr(Nb_VarC,me%NbTotBr,me%NbTotBr))
    me%DerVel_dBr=0.0_dp
    
    allocate(me%dynamic%br(me%NbTotBr))
    allocate(me%dynamic%p_star_prev(me%NbTotBr))
    me%dynamic%p_star_prev = 0.0_dp
    allocate(me%dynamic%alpha_prev(me%NbTotBr))
    me%dynamic%alpha_prev = 1.0e-1_dp

    do j=1,me%NbTotBr
        allocate(me%FlJ(j)%DerFlx(Nb_VarC,Nb_VarC,me%NbTotBr))
        me%FlJ(j)%DerFlx = 0.0_dp
        allocate(me%FlJ(j)%DerVit(Nb_VarC,me%NbTotBr))
        me%FlJ(j)%DerVit = 0.0_dp
    enddo
end subroutine all_junction_allocation


subroutine junction_init_part1(me,krn,cfg)
    type(junction_t), intent(out) :: me
    type(krn_t), intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg

    integer :: nb_non_zeros_exp, nb_non_zeros_imp

    me%NbTotBr=size(cfg%dict1d('link'))    

    nb_non_zeros_exp = 0
    nb_non_zeros_imp = me%NbTotBr * (me%NbTotBr-1) * Nb_VarC * Nb_VarC

    call krn%add('',0,nb_non_zeros_exp,nb_non_zeros_imp,0,0,0,0,0)
end subroutine junction_init_part1


subroutine junction_init_part2(me,krn,cfg)
    type(junction_t), intent(inout) :: me
    type(krn_t), target, intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg

    type(input_t), allocatable :: branches(:)
    integer :: i,j

    call all_junction_allocation(me)

    branches = cfg%dict1d('link')
    do i = 1, me%NbTotBr
        me%br(i)%p => krn%FF_flux_ports(krn%FF_flux_list%find(branches(i)%str('id'),&
                                                              convert2int_port(branches(i)%str('node'))))
        if(me%br(i)%p%connected) then
            print*,'Port ',branches(i)%str('node'),' belonging to ',branches(i)%str('id'),' already connected to a link'
            stop
        endif
        me%br(i)%p%connected=.true.
        me%br(i)%angle = branches(i)%dbl('angle',0.0_dp)
    enddo

    me%kappa=branches(1)%dbl('kappa',1.0_dp)

    do j=1,me%NbTotBr
      do i=1,me%NbTotBr
        if(j/=i) then
            call krn%coo_add_link(me%br(j)%p, me%br(i)%p, idx_4x4_col, &
                    idx_4x4_row, to_ptr_arr(me%Matrix(:,:,i,j)))
        endif
      enddo
    enddo
end subroutine junction_init_part2


subroutine junction_dynamic_parameters_fcn(me, n, x, fVec)

    class(junction_dynamic_parameters_t), intent(in) :: me
    integer, intent(in) :: n !! number of incoming branches + 2 * outgoing branches
    real(dp), intent(inout) :: x(n)
    real(dp), intent(inout) :: fVec(n)

    real(dp), dimension(me%NbOut) :: dFdRo, psi, q, Coef
    real(dp) :: EnthalpyMixStar, Denominator
    type(junction_dynamic_star_t) :: star(me%NbTotBr)
    integer :: i, i_out, i_in, i_ref, i_tot, i_kap
    i_tot = me%NbTotBr !! incoming branch stop index - total number of branches
    i_out = me%NbOut !! outgoing branch stop index - number of outgoing branches
    i_ref = i_out + 1  !! reference branch index
    i_in  = i_out + 2  !! incoming branch start index
    i_kap = me%idxKappaRef

    ! Incoming quantities as a function of p(j)
    associate(br => me%br(i_ref:i_tot), st => star(i_ref:i_tot), xx => x(i_ref:i_tot))

        ! Approximate jump relations through the rarefaction wave
        st%vit = br%vit0 - (xx - br%p0 - br%rho0*(br%q0 - br%qBar) + &
                br%dx0*br%Frc/2.0_dp) / (br%rho0 * (br%vit0 - br%SpeedS0))
        st%rho = br%rho0 * (br%vit0 - br%SpeedS0)/(st%vit - br%SpeedS0)

        if(R_Correction) then
            st%rCorr = (br%rCorr0*(br%vit0 - br%SpeedS0) - br%dx0*br%FrcR/2.0_dp) &
                / (st%vit-br%SpeedS0)
            st%ETot = (st%rCorr - xx)/st%rho + 0.5_dp*st%vit**2 + (br%q0 - br%qBar)        
        else
            st%ETot = br%ETot0 + (br%p0*br%vit0 - xx*st%vit) &
                / (br%rho0*(br%vit0 - br%SpeedS0)) + (br%q0-br%qBar)
        endif

        ! Variables involved in the system to solve
        st%q = st%rho * st%vit
        st%Enthalpy = st%ETot + xx/st%rho  ! Total enthalpy

        EnthalpyMixStar = sum(br%AreaJGB * st%q * st%Enthalpy)
        Denominator = sum(br%AreaJGB * st%q)
    end associate

    if(abs(Denominator) < 1.0e-10_dp) Denominator = 1.0e-10_dp
    EnthalpyMixStar = EnthalpyMixStar / Denominator

    ! Outgoing quantities as a function of p(j) and correction term alpha(j)
    ! through the contact discontinuity
    associate(br => me%br(1:i_out), st => star(1:i_out), xx => x(1:i_out))

        ! Approximate jump relations through the shock wave
        st%vit = br%vit0 - (xx - br%p0 - br%rho0*(br%q0 - br%qBar) + &
                br%dx0 * br%Frc / 2.0_dp ) / (br%rho0*(br%vit0 - br%SpeedS0))
        st%rho_int = br%rho0*(br%vit0 - br%SpeedS0) / (st%vit - br%SpeedS0)
        st%ETot_int = br%ETot0 + (br%p0*br%vit0 - xx*st%vit) / (br%rho0 * &
                (br%vit0 - br%SpeedS0)) + (br%q0 - br%qBar)

        st%QDM_int = st%rho_int * st%vit
        st%Energy_int = st%rho_int * st%ETot_int

        ! Correction through the contact discontinuity
        do i = 1, size(dFdRo)
            dFdRo(i) = droeint_droP(st(i)%rho_int, xx(i))
        enddo

        ! Additional term related to the real gas EOS  
        st%rho = st%rho_int + x(i_tot+1:n)
        st%QDM = st%QDM_int + x(i_tot+1:n) * st%vit
        st%Energy = st%Energy_int + x(i_tot+1:n) * (dFdRo + 0.5_dp*st%vit**2)

        ! Variables involved in the system to solve
        st%q = st%QDM
        st%Enthalpy = (st%Energy + xx)/st%rho
    
        ! SYSTEM TO SOLVE : fVec(:)=0.0_dp
        fVec(1) = sum(me%br(1:i_tot)%AreaJGB * star(1:i_tot)%q)

        ! The dat branch is supposed to be the first given incoming pipe
        Coef = 0.0_dp
        where (abs(st%q)>=epsCoef .and. abs(star(i_ref)%q)>=epsCoef)
            q = -br%AreaJGB * st%q / (me%br(i_ref)%AreaJGB * star(i_ref)%q)
            psi = me%br(i_ref)%AreaJGB / br%AreaJGB
            Coef = 1.0_dp - cos((3.0_dp/4.0_dp)*(3.14159_dp-br%ThetaJGB))/(q*psi)
        end where

        fVec(2:i_ref) = x(i_ref) - xx - Coef * st%rho * st%vit**2
        if(me%idxKappaRef == i_ref) &
            fVec(2:i_ref) = x(i_ref) - xx - Coef * me%kappa * st%rho * st%vit**2
        if(me%idxKappaRef < i_ref) &
            fVec(i_kap+1) = x(i_ref) - x(i_kap) - &
                Coef(i_kap) * me%kappa * star(i_kap)%rho * star(i_kap)%vit**2
    end associate
    ! Incoming pipes
    ! Coef_j=0.0_dp ! see Bassett (2003)
    ! - Coef_j*star%rho(i_in:i_tot)*vitStar(i_in:i_tot)**2    
    fVec(i_in:i_tot) = x(i_ref) - x(i_in:i_tot)
    fVec(i_tot+1:n) = star(1:i_out)%Enthalpy - EnthalpyMixStar
end subroutine junction_dynamic_parameters_fcn

end module cmp_junction_init_m