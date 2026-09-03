! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_channel_init_m
    use krn_interface_m
    use krn_simulation_m
    use lib_input_m, only: input_t
    use lib_hdf_write_m
    use lib_He_thermo_m
    use lib_friction_correlations_m
    use lib_nusselt_correlations_m
    implicit none

    type source_terms_channels_t
        real(dp) :: Se 
        real(dp) :: DerSe(Nb_VarC)
        real(dp), allocatable :: Su(:),Sr(:)
        real(dp), allocatable :: DerSu(:,:),DerSr(:,:)
    end type source_terms_channels_t 
 
    type He_prop_t
        integer :: NbCells
        real(dp), allocatable :: dxLoc(:)
        real(dp) :: Length  !! channel length
        real(dp) :: Diam    !! channel diameter
        real(dp) :: Area    !! channel cross section area
        logical :: FFsrcLink !! Fluid connexion with another channel
        character(:), allocatable :: ExtHeating !! Heat exchange activation - T or Q imposed      
        type(signal_t) :: Temp_or_Q_Wall !! Temperature of the wall OR Flux
        real(dp) :: Ht  !! Heat transfer coefficient in case of imposed value        
        type(FF_flux_port_t), pointer :: in  !! inlet fluid port
        type(FF_flux_port_t), pointer :: out !! outlet fluid port
        type(FS_port_pointer_t), allocatable :: thermP(:)
        type(FF_src_port_pointer_t), allocatable :: thermF(:)
        character(:), allocatable :: name !! channel name
    end type He_prop_t

    type StateVariable_channel_t
        real(dp), allocatable :: He_pr(:,:) !! He primitive variables
        real(dp), allocatable :: He_cs(:,:) !! He conservative variables
    end type StateVariable_channel_t

    type big_arrays_channel_t
        real(dp), allocatable :: amx(:,:,:) !! lower diagonal subMatrix
        real(dp), allocatable :: bmx(:,:,:) !! main diagonal subMatrix
        real(dp), allocatable :: cmx(:,:,:) !! upper diagonal subMatrix
        
        real(dp), pointer :: bvx(:,:)    
        real(dp), allocatable :: Jacob(:,:,:)  !! Jacobian matrix of fluid model
    end type big_arrays_channel_t  

    type arrays_linear_system_channel_t
        type(dbl_pointer_t), allocatable :: ValExp(:), ValImp(:) !! Non zero values in COO format for 1 channel
        integer, allocatable :: RowExp(:), RowImp(:)  !! row index array in COO format for 1 channel
        integer, allocatable :: ColExp(:), ColImp(:)  !! col index array in COO format for 1 channel
    end type arrays_linear_system_channel_t     

    type flux_He_channel_t
        real(dp), allocatable :: Cons(:,:)   !! He flux (Conservative terms)
        real(dp), allocatable :: VitTNC(:)   !! He flux (Velocity - Non Conservative Term)
        real(dp), allocatable :: Cons_DerdUL(:,:,:), Cons_DerdUR(:,:,:)
        real(dp), allocatable :: VitTNC_DerdUL(:,:), VitTNC_DerdUR(:,:)
    end type flux_He_channel_t  

    type channel_t
        type(hdf_desc_t)                     :: hdf
        type(StateVariable_channel_t)        :: StVar     !! StateVariable updated by BDF 2 scheme -> time n+1
        type(StateVariable_channel_t)        :: StVarOld  !! StateVariable -> time n
        type(StateVariable_channel_t)        :: StVarOld2 !! StateVariable -> time n-1
        type(StateVariable_channel_t)        :: StVarOld3 !! StateVariable -> time n-2
        type(StateVariable_channel_t)        :: StVarOld4 !! StateVariable -> time n-3
        type(StateVariable_channel_t)        :: StVarOld5 !! StateVariable -> time n-4
        type(StateVariable_channel_t)        :: StVarOld6 !! StateVariable -> time n-5
        type(He_prop_t)                      :: HeProp    !! Variables relative to Helium channel global properties
        type(flux_He_channel_t)              :: flxHe     !! Fluxes - He flow
        type(big_arrays_channel_t)           :: big       !! Big arrays
        type(arrays_linear_system_channel_t) :: bigLS       !! Big arrays for Linear system
        type(source_terms_channels_t)        :: srcPP
        real(dp)                             :: wave_time   !! Minimum time for wave propagation
        real(dp)                             :: err
        real(dp)                             :: err_den
        class(friction_t), pointer           :: fct         !! Friction correlation
        class(nusselt_t),  pointer           :: nuss        !! Nusselt correlation
    end type channel_t

contains

subroutine all_channel_allocation_nb(me)
    type(channel_t), intent(inout) :: me
    integer :: NbCells
      
    NbCells=me%HeProp%NbCells
    allocate(me%StVar%He_pr(Nb_VarP,NbCells))
    allocate(me%StVar%He_cs(Nb_VarC,NbCells))
    allocate(me%StVarOld%He_pr(Nb_VarP,NbCells))
    allocate(me%StVarOld%He_cs(Nb_VarC,NbCells))
    allocate(me%StVarOld2%He_pr(Nb_VarP,NbCells))
    allocate(me%StVarOld2%He_cs(Nb_VarC,NbCells))
    allocate(me%StVarOld3%He_pr(Nb_VarP,NbCells))
    allocate(me%StVarOld3%He_cs(Nb_VarC,NbCells))
    allocate(me%StVarOld4%He_pr(Nb_VarP,NbCells))
    allocate(me%StVarOld4%He_cs(Nb_VarC,NbCells))
    allocate(me%StVarOld5%He_pr(Nb_VarP,NbCells))
    allocate(me%StVarOld5%He_cs(Nb_VarC,NbCells))
    allocate(me%StVarOld6%He_pr(Nb_VarP,NbCells))
    allocate(me%StVarOld6%He_cs(Nb_VarC,NbCells))
    allocate(me%flxHe%Cons(Nb_VarC,0:NbCells))
    allocate(me%flxHe%VitTNC(0:NbCells))
    allocate(me%flxHe%Cons_DerdUL(Nb_VarC,Nb_VarC,0:NbCells))
    allocate(me%flxHe%Cons_DerdUR(Nb_VarC,Nb_VarC,0:NbCells))
    allocate(me%flxHe%VitTNC_DerdUL(Nb_VarC,0:NbCells))
    allocate(me%flxHe%VitTNC_DerdUR(Nb_VarC,0:NbCells))
    allocate(me%big%amx(Nb_VarC,Nb_VarC,NbCells-1))
    allocate(me%big%cmx(Nb_VarC,Nb_VarC,NbCells-1))
    allocate(me%big%bmx(Nb_VarC,Nb_VarC,NbCells))
    allocate(me%big%Jacob(Nb_VarC,Nb_VarC,NbCells))
    allocate(me%srcPP%Su(NbCells))
    allocate(me%srcPP%DerSu(Nb_VarC,NbCells))
    allocate(me%srcPP%Sr(NbCells))
    allocate(me%srcPP%DerSr(Nb_VarC,NbCells))
      
    me%StVar%He_pr=0.0_dp; me%StVar%He_cs=0.0_dp
    me%StVarOld%He_pr=0.0_dp; me%StVarOld%He_cs=0.0_dp
    me%StVarOld2%He_pr=0.0_dp; me%StVarOld2%He_cs=0.0_dp
    me%StVarOld3%He_pr=0.0_dp; me%StVarOld3%He_cs=0.0_dp
    me%StVarOld4%He_pr=0.0_dp; me%StVarOld4%He_cs=0.0_dp
    me%StVarOld5%He_pr=0.0_dp; me%StVarOld5%He_cs=0.0_dp
    me%StVarOld6%He_pr=0.0_dp; me%StVarOld6%He_cs=0.0_dp
    me%flxHe%Cons=0.0_dp;me%flxHe%VitTNC=0.0_dp
    me%flxHe%Cons_DerdUL=0.0_dp;me%flxHe%Cons_DerdUR=0.0_dp
    me%flxHe%VitTNC_DerdUL=0.0_dp;me%flxHe%VitTNC_DerdUR=0.0_dp
    me%big%amx=0.0_dp;me%big%cmx=0.0_dp;me%big%bmx=0.0_dp
    me%big%Jacob=0.0_dp;me%srcPP%Su=0.0_dp;me%srcPP%DerSu=0.0_dp
    me%srcPP%Sr=0.0_dp;me%srcPP%DerSr=0.0_dp
end subroutine all_channel_allocation_nb

subroutine v_prim_to_v_Conss(i,Prim,Cons)  
    integer, intent(in) :: i
    real(dp), intent(in) :: Prim(:,:)
    real(dp), intent(inout) :: Cons(:,:)
    
    Cons(Con_Mas,i)=Prim(Pri_ro,i)
    Cons(Con_Qdm,i)=Prim(Pri_ro,i)*Prim(Pri_u,i)
    Cons(Con_Ene,i)=Prim(Pri_ro,i)*(Prim(Pri_e,i)+0.5_dp*Prim(Pri_u,i)**2)
    Cons(Con_R,i)=Prim(Pri_R,i)
end subroutine v_prim_to_v_Conss

subroutine channel_init_part1(me, krn, cfg, h5, sim)
    type(channel_t),         intent(out)   :: me
    type(krn_t),          intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg
    class(hdf5_t),        intent(inout) :: h5
    type(simulation_t), intent(inout) :: sim

    real(dp) :: Rout,Pout,eout,cout, u_init, p_init, T_init
    integer :: i, nb_non_zeros_exp, nb_non_zeros_imp, nb_FS_ports, nb_FF_src_ports
            
    if(cfg%has_key('nodes')) then ! means uniform mesh
        me%HeProp%NbCells = cfg%int('nodes')
        me%HeProp%Length = cfg%dbl('length')
        allocate(me%HeProp%dxLoc(me%HeProp%NbCells))
        me%HeProp%dxLoc(:)=me%HeProp%Length/me%HeProp%NbCells
    else ! variable mesh
        me%HeProp%NbCells=size(cfg%dbl1d('length'))
        allocate(me%HeProp%dxLoc(me%HeProp%NbCells))
        me%HeProp%dxLoc(:)=cfg%dbl1d('length')
        me%HeProp%Length = sum(me%HeProp%dxLoc(:))
    endif

    ! For HDF5
    me%hdf%name = cfg%str('id')
    me%hdf%node_x = [(sum(me%HeProp%dxLoc(1:i-1)) + me%HeProp%dxLoc(i)/2, &
                    i = 1, me%HeProp%NbCells)]
    call h5%add_to_table(me%hdf,'channel')

    me%HeProp%FFsrcLink = cfg%bin('channel_link',.false.)

    me%HeProp%ExtHeating = cfg%str('thermal','flux')
    if(me%HeProp%ExtHeating=='flux') then
        call me%HeProp%Temp_or_Q_Wall%init(sim,cfg,'flux',me%hdf%node_x)
    else if(me%HeProp%ExtHeating=='temp') then
        call me%HeProp%Temp_or_Q_Wall%init(sim,cfg,'temp',me%hdf%node_x)
        call nusselt_init(me%nuss,cfg)        
        if(.not. associated(me%nuss)) me%HeProp%Ht = cfg%dbl('heat_transfert_coef')
    endif

    call friction_init(me%fct,cfg)     
    if(cfg%has_key('diameter') .and. cfg%has_key('area')) then  ! expects the 2 values for a Katheder type friction factor
        me%HeProp%Diam = cfg%dbl('diameter')
        me%HeProp%Area = cfg%dbl('area')
    else
        if(cfg%has_key('diameter')) then
            me%HeProp%Diam = cfg%dbl('diameter')
            me%HeProp%Area = Pi_value*me%HeProp%Diam**2/4.0_dp
        else if(cfg%has_key('area')) then
            me%HeProp%Area = cfg%dbl('area')
            me%HeProp%Diam = sqrt(4.0_dp*me%HeProp%Area/Pi_value)
        else
            error stop 'channel: missing required key "diameter" or "area"'
        endif
    endif

    call all_channel_allocation_nb(me)
          
    u_init=cfg%dbl('initial/u',0.0_dp)
    p_init=cfg%dbl('initial/p')
    T_init=cfg%dbl('initial/t')
    do i=1,me%HeProp%NbCells
        me%StVar%He_pr(Pri_u,i)=u_init
        me%StVar%He_pr(Pri_ro,i)=ro_pT(p_init,T_init)
        me%StVar%He_pr(Pri_T,i)=T_init
        call state_roT(me%StVar%He_pr(Pri_ro,i), me%StVar%He_pr(Pri_T,i), Rout, Pout, eout, cout)

        me%StVar%He_pr(Pri_R,i) = 1.0_dp ! Value by default
        if(R_Correction) me%StVar%He_pr(Pri_R,i) = Rout
        me%StVar%He_pr(Pri_p,i)=Pout
        me%StVar%He_pr(Pri_e,i)=eout
        me%StVar%He_pr(Pri_c,i)=cout
        call v_prim_to_v_Conss(i,me%StVar%He_pr,me%StVar%He_cs)
    enddo    
        
    nb_non_zeros_exp=Nb_VarC*Nb_VarC*me%HeProp%NbCells

    nb_non_zeros_imp=2*Nb_VarC*Nb_VarC
    do i=2,me%HeProp%NbCells-1
        nb_non_zeros_imp=nb_non_zeros_imp+3*Nb_VarC*Nb_VarC
    enddo
    nb_non_zeros_imp=nb_non_zeros_imp+2*Nb_VarC*Nb_VarC

    if(me%HeProp%ExtHeating=='link') then
        nb_FS_ports=me%HeProp%NbCells
    else
        nb_FS_ports=0
    endif
    
    if(me%HeProp%FFsrcLink) then
        nb_FF_src_ports=me%HeProp%NbCells
    else
        nb_FF_src_ports=0
    endif    

    me%HeProp%name = cfg%str('id')

    call krn%add(me%HeProp%name,Nb_VarC*me%HeProp%NbCells,nb_non_zeros_exp,nb_non_zeros_imp,2,nb_FS_ports,nb_FF_src_ports,0,0)
    
end subroutine channel_init_part1


subroutine channel_init_part2(me, krn, cfg)
    !! Fluid port initialistion for channels
    type(channel_t), target, intent(inout) :: me
    type(krn_t),          intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg  

    type(FF_flux_port_t), pointer :: FF_flux_ports(:)
    type(FF_src_port_t), pointer :: FF_src_ports(:)
    type(FS_port_t), pointer :: FS_ports(:)
    integer :: NbSubM,i,idx,idxI,idxJ,PrevRowIdx,PrevColIdx
    real(dp), pointer :: rhs(:)
    integer, allocatable :: list_thp_loc(:),list_FF_src_p_loc(:)

    ! FF_flux ports
    call krn%update(rhs_or_solution_view=rhs,FF_flux_p_loc=[1,(me%HeProp%NbCells-1)*Nb_VarC+1],FF_flux_p_view=FF_flux_ports)
    me%HeProp%in  => FF_flux_ports(1)
    me%HeProp%out => FF_flux_ports(2)
    me%HeProp%in%Area       = me%HeProp%Area
    me%HeProp%in%dxLoc      = me%HeProp%dxLoc(1)
    me%HeProp%in%Sgn4j      = -1.0_dp   
    me%HeProp%in%connected  = .false.  
    me%HeProp%out%Area      = me%HeProp%Area
    me%HeProp%out%dxLoc     = me%HeProp%dxLoc(me%HeProp%NbCells)
    me%HeProp%out%Sgn4j     = 1.0_dp   
    me%HeProp%out%connected = .false.   

    ! FS ports
    if(me%HeProp%ExtHeating=='link') then
      allocate(list_thp_loc(me%HeProp%NbCells),me%HeProp%thermP(me%HeProp%NbCells))
      list_thp_loc=[(1+(i-1)*Nb_VarC, i=1,me%HeProp%NbCells)]
      call krn%update(FS_p_loc=list_thp_loc,FS_p_view=FS_ports)
      do i=1,me%HeProp%NbCells
        me%HeProp%thermP(i)%p => FS_ports(i)
        me%HeProp%thermP(i)%p%AreaP = me%HeProp%Area
        me%HeProp%thermP(i)%p%DiamP = me%HeProp%Diam
        me%HeProp%thermP(i)%p%dxLoc = me%HeProp%dxLoc(i)
      enddo
    endif

    ! FF_src ports
    if(me%HeProp%FFsrcLink) then
      allocate(list_FF_src_p_loc(me%HeProp%NbCells),me%HeProp%thermF(me%HeProp%NbCells))
      list_FF_src_p_loc=[(1+(i-1)*Nb_VarC, i=1,me%HeProp%NbCells)]
      call krn%update(FF_src_p_loc=list_FF_src_p_loc,FF_src_p_view=FF_src_ports)
      do i=1,me%HeProp%NbCells
        me%HeProp%thermF(i)%p => FF_src_ports(i)
        me%HeProp%thermF(i)%p%AreaP = me%HeProp%Area
        me%HeProp%thermF(i)%p%DiamP = me%HeProp%Diam
      enddo
      me%HeProp%thermF(1)%p%can_be_shared = .true.
      me%HeProp%thermF(me%HeProp%NbCells)%p%can_be_shared = .true.
    endif       

    me%big%bvx(1:Nb_VarC,1:me%HeProp%NbCells) => rhs

    NbSubM=me%HeProp%NbCells
    allocate(me%bigLS%ValExp(NbSubM*Nb_VarC*Nb_VarC),me%bigLS%RowExp(NbSubM*Nb_VarC*Nb_VarC))
    allocate(me%bigLS%ColExp(NbSubM*Nb_VarC*Nb_VarC))

    NbSubM=2*2+(me%HeProp%NbCells-2)*3
    allocate(me%bigLS%ValImp(NbSubM*Nb_VarC*Nb_VarC),me%bigLS%RowImp(NbSubM*Nb_VarC*Nb_VarC))
    allocate(me%bigLS%ColImp(NbSubM*Nb_VarC*Nb_VarC))

    idx=1

    PrevColIdx=0
    PrevRowIdx=0
    do i=1,me%HeProp%NbCells
        ! PPt%bigLS%ColExp((i-1)*16+1:i*16) = PrevColIdx+idx_4x4_col
        ! PPt%bigLS%RowExp((i-1)*16+1:i*16) = PrevRowIdx+idx_4x4_row
        ! PPt%bigLS%ValExp((i-1)*16+1:i*16) = to_ptr_arr(PPt%bigPP%bmxPP(:,:,i)) ! TODO adapt this part
        do idxI=1,Nb_VarC
            do idxJ=1,Nb_VarC
                me%bigLS%ValExp(idx)%p => me%big%bmx(idxJ,idxI,i)
                me%bigLS%ColExp(idx)=PrevColIdx+idxJ
                me%bigLS%RowExp(idx)=PrevRowIdx+idxI

                me%bigLS%ValImp(idx)%p => me%big%bmx(idxJ,idxI,i)
                me%bigLS%ColImp(idx)=PrevColIdx+idxJ
                me%bigLS%RowImp(idx)=PrevRowIdx+idxI
                idx=idx+1
            enddo  
        enddo
        PrevColIdx=PrevColIdx+Nb_VarC
        PrevRowIdx=PrevRowIdx+Nb_VarC
    enddo
    call krn%coo_add(.true.,me%bigLS%ColExp,me%bigLS%RowExp,me%bigLS%ValExp)

    PrevColIdx=Nb_VarC
    PrevRowIdx=0
    do i=1,me%HeProp%NbCells-1
        do idxI=1,Nb_VarC
            do idxJ=1,Nb_VarC
                me%bigLS%ValImp(idx)%p => me%big%cmx(idxJ,idxI,i)
                me%bigLS%ColImp(idx)=PrevColIdx+idxJ
                me%bigLS%RowImp(idx)=PrevRowIdx+idxI          
                idx=idx+1
            enddo  
        enddo
        PrevColIdx=PrevColIdx+Nb_VarC
        PrevRowIdx=PrevRowIdx+Nb_VarC      
    enddo
    PrevColIdx=0
    PrevRowIdx=Nb_VarC
    do i=1,me%HeProp%NbCells-1
        do idxI=1,Nb_VarC
            do idxJ=1,Nb_VarC
                me%bigLS%ValImp(idx)%p => me%big%amx(idxJ,idxI,i)
                me%bigLS%ColImp(idx)=PrevColIdx+idxJ
                me%bigLS%RowImp(idx)=PrevRowIdx+idxI 
                idx=idx+1
            enddo  
        enddo
        PrevColIdx=PrevColIdx+Nb_VarC
        PrevRowIdx=PrevRowIdx+Nb_VarC
    enddo
    call krn%coo_add(.false.,me%bigLS%ColImp,me%bigLS%RowImp,me%bigLS%ValImp)

    deallocate(me%bigLS%ValExp,me%bigLS%RowExp,me%bigLS%ColExp)
    deallocate(me%bigLS%ValImp,me%bigLS%RowImp,me%bigLS%ColImp)

end subroutine channel_init_part2

subroutine consistency_ports(PP)
    type(channel_t), intent(in) :: PP(:)
    integer :: i

    do i = 1, size(PP)
        if(PP(i)%HeProp%in%connected==.false.) then
            print*,'FF_flux_port 1 of ',PP(i)%hdf%name,' is not connected to any link'
            read(*,*)
        endif

        if(PP(i)%HeProp%out%connected==.false.) then
            print*,'FF_flux_port 2 of ',PP(i)%hdf%name,' is not connected to any link'
            read(*,*)
        endif   
    enddo
    flush(stdout)
end subroutine consistency_ports

end module cmp_channel_init_m