! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_strand_init_m
    use krn_interface_m
    use krn_simulation_m
    use lib_input_m, only: input_t
    use lib_hdf_write_m
    use lib_material_m
    implicit none

    type scenario_t
      type(signal_t) :: Q_ext_load,Bfield_load,dBfield_load,ElCur_load
      real(dp), allocatable :: Q_ext(:),Bfield(:),dBfield(:)
      real(dp) :: ElCur
    end type scenario_t    

    type constituent_t
        character(:), allocatable ::  type
        real(dp) :: area
    end type constituent_t

    type strand_prop_t
        integer :: NbCells
        real(dp), allocatable :: dxLoc(:)        
        real(dp) :: Length !! Strand length
        logical  :: FSlink
        real(dp) :: inner_rad  !! CICC inner radius 
        real(dp) :: outer_rad  !! CICC outer radius
        type(FS_port_pointer_t), allocatable :: thermP(:)
        type(SS_flux_port_t), pointer       :: in             !! inlet SS_flux port
        type(SS_flux_port_t), pointer       :: out            !! outlet SS_flux port     
    end type strand_prop_t 

    type StateVariable_strand_t
        real(dp), allocatable :: SCtemp(:) !! Temperature of composite strand
    end type StateVariable_strand_t    

    type big_arrays_strand_t
        real(dp), allocatable :: amS(:) !! lower diagonal subMatrix
        real(dp), allocatable :: bmS(:) !! main diagonal subMatrix    
        real(dp), allocatable :: cmS(:) !! upper diagonal subMatrix

        real(dp), pointer :: bvS(:)    
    end type big_arrays_strand_t     

    type arrays_linear_system_strand_t
        type(dbl_pointer_t), allocatable :: ValExp(:), ValImp(:) !! Non zero values in COO format for 1 strand
        integer, allocatable :: RowExp(:), RowImp(:)  !! row index array in COO format for 1 strand
        integer, allocatable :: ColExp(:), ColImp(:)  !! col index array in COO format for 1 strand
    end type arrays_linear_system_strand_t      

    type flux_strand_t
        real(dp), allocatable :: wire(:)   !! heat diffusion flux
        real(dp), allocatable :: wire_DerdTL(:), wire_DerdTR(:)
    end type flux_strand_t

    type strand_t
        type(hdf_desc_t)                           ::  hdf      
        type(StateVariable_strand_t)               ::  StVar          !! StateVariable updated by BDF 2 scheme -> time n+1
        type(StateVariable_strand_t)               ::  StVarOld       !! StateVariable -> time n
        type(StateVariable_strand_t)               ::  StVarOld2      !! StateVariable -> time n-1
        type(StateVariable_strand_t)               ::  StVarOld3      !! StateVariable -> time n-2
        type(StateVariable_strand_t)               ::  StVarOld4      !! StateVariable -> time n-3
        type(StateVariable_strand_t)               ::  StVarOld5      !! StateVariable -> time n-4
        type(StateVariable_strand_t)               ::  StVarOld6      !! StateVariable -> time n-5
        type(strand_prop_t)                        ::  SC_Prop        !! Variables relative to strand global properties
        integer                                    ::  NbConst        !! Number of constituents
        type(constituent_t), allocatable           ::  ct(:)          !! Constituent fixed properties
        real(dp)                                   ::  ro_M           !! Density of the composite
        type(scenario_t)                           ::  scen           !! Electro-magnetic + Heat scenario
        type(big_arrays_strand_t)                  ::  big            !! Big arrays
        type(arrays_linear_system_strand_t)        ::  bigLS          !! Big arrays for Linear system
        type(flux_strand_t)                        ::  flxS           !! Heat diffusion fluxes
        real(dp)                                   ::  err
        real(dp)                                   ::  err_den
        class(material_t), pointer                 ::  mat_stab, mat_supc  !! Material dependent properties
    end type strand_t

    ! Index definitions
    integer,parameter :: STAB = 1  !! Stabilizer
    integer,parameter :: SUPC = 2  !! Superconductor

contains

subroutine all_strand_allocation(me)
    type(strand_t), intent(inout) :: me
    integer :: NbCells

    allocate(me%ct(me%NbConst))
      
    NbCells=me%SC_Prop%NbCells
    allocate(me%StVar%SCtemp(NbCells))
    allocate(me%StVarOld%SCtemp(NbCells))
    allocate(me%StVarOld2%SCtemp(NbCells))
    allocate(me%StVarOld3%SCtemp(NbCells))
    allocate(me%StVarOld4%SCtemp(NbCells))
    allocate(me%StVarOld5%SCtemp(NbCells))
    allocate(me%StVarOld6%SCtemp(NbCells))
    allocate(me%flxS%wire(0:NbCells))
    allocate(me%flxS%wire_DerdTR(0:NbCells))
    allocate(me%flxS%wire_DerdTL(0:NbCells))
    allocate(me%scen%Q_ext(NbCells))
    allocate(me%scen%Bfield(NbCells))
    allocate(me%scen%dBfield(NbCells))
    allocate(me%big%amS(NbCells-1))
    allocate(me%big%cmS(NbCells-1))
    allocate(me%big%bmS(NbCells))
      
    me%StVar%SCtemp=0.0_dp; me%StVarOld%SCtemp=0.0_dp
    me%StVarOld2%SCtemp=0.0_dp; me%StVarOld3%SCtemp=0.0_dp
    me%StVarOld4%SCtemp=0.0_dp; me%StVarOld5%SCtemp=0.0_dp; me%StVarOld6%SCtemp=0.0_dp
    me%flxS%wire=0.0_dp
    me%flxS%wire_DerdTR=0.0_dp; me%flxS%wire_DerdTL=0.0_dp
    me%scen%Q_ext=0.0_dp;me%scen%Bfield=0.0_dp;me%scen%dBfield=0.0_dp
    me%big%amS=0.0_dp;me%big%bmS=0.0_dp;me%big%cmS=0.0_dp
end subroutine all_strand_allocation


subroutine strand_init_part1(me, krn, cfg, h5, sim)
    type(strand_t),         intent(out) :: me
    type(krn_t),          intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg
    class(hdf5_t),        intent(inout) :: h5
    type(simulation_t), intent(inout) :: sim

    real(dp) :: T_init, num, den
    integer :: i, ii, nb_non_zeros_exp, nb_non_zeros_imp, nb_FS_ports

    if(cfg%has_key('nodes')) then ! means uniform mesh
        me%SC_Prop%NbCells = cfg%int('nodes')
        me%SC_Prop%Length = cfg%dbl('length')
        allocate(me%SC_Prop%dxLoc(me%SC_Prop%NbCells))
        me%SC_Prop%dxLoc(:)=me%SC_Prop%Length/me%SC_Prop%NbCells
    else ! variable mesh
        me%SC_Prop%NbCells=size(cfg%dbl1d('length'))
        allocate(me%SC_Prop%dxLoc(me%SC_Prop%NbCells))
        me%SC_Prop%dxLoc(:)=cfg%dbl1d('length')
        me%SC_Prop%Length = sum(me%SC_Prop%dxLoc(:))
    endif    

    ! For HDF5
    me%hdf%name = cfg%str('id')
    me%hdf%node_x = [(sum(me%SC_Prop%dxLoc(1:i-1)) + me%SC_Prop%dxLoc(i)/2, i = 1, me%SC_Prop%NbCells)]    
    call h5%add_to_table(me%hdf,'strand')

    me%NbConst=2
    call all_strand_allocation(me)     

    T_init=cfg%dbl('initial/t')
    do ii=1,me%SC_Prop%NbCells
        me%StVar%SCtemp(ii)=T_init
    enddo

    call material_init(me%mat_stab,cfg%dict('stabilizer'))
    call material_init(me%mat_supc,cfg%dict('superconductor'))

    me%ct(1)%area = cfg%dbl('stabilizer/area')
    me%ct(1)%type = 'stabilizer'

    me%ct(2)%area = cfg%dbl('superconductor/area')
    me%ct(2)%type = 'superconductor'

    num=me%ct(STAB)%area*me%mat_stab%density+me%ct(SUPC)%area*me%mat_supc%density
    den=me%ct(STAB)%area+me%ct(SUPC)%area
    me%ro_M=num/den

    me%SC_Prop%inner_rad=cfg%dbl('CICC_inner_radius')
    me%SC_Prop%outer_rad=cfg%dbl('CICC_outer_radius')

    me%SC_Prop%FSlink = cfg%bin('channel_link')

    ! default values for 1d_signal
    call me%scen%Q_ext_load%init(sim,cfg,'flux',me%hdf%node_x)
    call me%scen%Bfield_load%init(sim,cfg,'field',me%hdf%node_x)
    call me%scen%dBfield_load%init(sim,cfg,'field_gradient',me%hdf%node_x)
    call me%scen%ElCur_load%init(sim,cfg,'current')
        
    nb_non_zeros_exp=me%SC_Prop%NbCells

    nb_non_zeros_imp=2
    do i=2,me%SC_Prop%NbCells-1
        nb_non_zeros_imp=nb_non_zeros_imp+3
    enddo
    nb_non_zeros_imp=nb_non_zeros_imp+2

    if(me%SC_Prop%FSlink) then
        nb_FS_ports=me%SC_Prop%NbCells
    else
        nb_FS_ports=0
    endif

    call krn%add(cfg%str('id'),me%SC_Prop%NbCells,nb_non_zeros_exp,nb_non_zeros_imp,0,nb_FS_ports,0,2,0)
    
end subroutine strand_init_part1


subroutine strand_init_part2(me, krn, cfg)
    !! Fluid port initialisation for strands
    type(strand_t), target, intent(inout) :: me
    type(krn_t),          intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg  

    integer :: NbSubM,i,idx,PrevRowIdx,PrevColIdx
    real(dp), pointer :: rhs(:)
    integer, allocatable :: list_thp_loc(:)
    type(FS_port_t), pointer :: FS_ports(:)
    type(SS_flux_port_t), pointer :: SS_flux_ports(:)

    ! SS_flux ports --> always 2 ports implemented
    call krn%update(rhs_or_solution_view=rhs,SS_flux_p_loc=[1,me%SC_Prop%NbCells],SS_flux_p_view=SS_flux_ports)
    me%SC_Prop%in  => SS_flux_ports(1)
    me%SC_Prop%out => SS_flux_ports(2)
    me%SC_Prop%in%AreaS   = sum(me%ct(:)%area)
    me%SC_Prop%in%rhoMS   = me%ro_M
    me%SC_Prop%in%Sgn4lk  = -1.0_dp   
    me%SC_Prop%in%dxLoc   = me%SC_Prop%dxLoc(1)
    me%SC_Prop%in%Flx           = 0.0_dp
    me%SC_Prop%in%derFlx_derCon = 0.0_dp
    me%SC_Prop%out%AreaS  = sum(me%ct(:)%area) 
    me%SC_Prop%out%rhoMS  = me%ro_M
    me%SC_Prop%out%Sgn4lk = 1.0_dp    
    me%SC_Prop%out%dxLoc  = me%SC_Prop%dxLoc(me%SC_Prop%NbCells)
    me%SC_Prop%out%Flx           = 0.0_dp
    me%SC_Prop%out%derFlx_derCon = 0.0_dp      
    
    me%big%bvS(1:me%SC_Prop%NbCells) => rhs

    if(me%SC_Prop%FSlink) then
      allocate(list_thp_loc(me%SC_Prop%NbCells),me%SC_Prop%thermP(me%SC_Prop%NbCells))
      list_thp_loc=[(i,i=1,me%SC_Prop%NbCells)]
      call krn%update(FS_p_loc=list_thp_loc,CompName=cfg%str('id'),FS_p_view=FS_ports)
      do i=1,me%SC_Prop%NbCells
          me%SC_Prop%thermP(i)%p => FS_ports(i)
          me%SC_Prop%thermP(i)%p%typeS = 'strand'
          me%SC_Prop%thermP(i)%p%AreaS = sum(me%ct(:)%area)
          me%SC_Prop%thermP(i)%p%rhoMS = me%ro_M
      enddo
    endif

    NbSubM=me%SC_Prop%NbCells
    allocate(me%bigLS%ValExp(NbSubM),me%bigLS%RowExp(NbSubM),me%bigLS%ColExp(NbSubM))

    NbSubM=2*2+(me%SC_Prop%NbCells-2)*3
    allocate(me%bigLS%ValImp(NbSubM),me%bigLS%RowImp(NbSubM),me%bigLS%ColImp(NbSubM))

    idx=1

    PrevColIdx=0
    PrevRowIdx=0
    do i=1,me%SC_Prop%NbCells
        me%bigLS%ValExp(idx)%p => me%big%bmS(i)
        me%bigLS%ColExp(idx)=PrevColIdx+1
        me%bigLS%RowExp(idx)=PrevRowIdx+1

        me%bigLS%ValImp(idx)%p => me%big%bmS(i)
        me%bigLS%ColImp(idx)=PrevColIdx+1
        me%bigLS%RowImp(idx)=PrevRowIdx+1

        idx=idx+1
        PrevColIdx=PrevColIdx+1
        PrevRowIdx=PrevRowIdx+1
    enddo
    call krn%coo_add(.true.,me%bigLS%ColExp,me%bigLS%RowExp,me%bigLS%ValExp)

    PrevColIdx=1
    PrevRowIdx=0
    do i=1,me%SC_Prop%NbCells-1
        me%bigLS%ValImp(idx)%p => me%big%cmS(i)
        me%bigLS%ColImp(idx)=PrevColIdx+1
        me%bigLS%RowImp(idx)=PrevRowIdx+1          
        
        idx=idx+1
        PrevColIdx=PrevColIdx+1
        PrevRowIdx=PrevRowIdx+1      
    enddo
    PrevColIdx=0
    PrevRowIdx=1
    do i=1,me%SC_Prop%NbCells-1
        me%bigLS%ValImp(idx)%p => me%big%amS(i)
        me%bigLS%ColImp(idx)=PrevColIdx+1
        me%bigLS%RowImp(idx)=PrevRowIdx+1 
        
        idx=idx+1
        PrevColIdx=PrevColIdx+1
        PrevRowIdx=PrevRowIdx+1
    enddo
    call krn%coo_add(.false.,me%bigLS%ColImp,me%bigLS%RowImp,me%bigLS%ValImp)

    deallocate(me%bigLS%ValExp,me%bigLS%RowExp,me%bigLS%ColExp)
    deallocate(me%bigLS%ValImp,me%bigLS%RowImp,me%bigLS%ColImp)

end subroutine strand_init_part2

end module cmp_strand_init_m