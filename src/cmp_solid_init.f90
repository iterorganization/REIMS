! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_solid_init_m
    use krn_interface_m
    use krn_simulation_m, only: simulation_t, signal_t
    use lib_input_m, only: input_t
    use lib_hdf_write_m
    use lib_material_m
    implicit none

    type solid_prop_t
        integer  :: NbCells, nb_FS_ports
        real(dp), allocatable :: dxLoc(:), volLoc(:)
        real(dp) :: Area        !! Cross-section area
        logical  :: SSsrcLink, cond_btw_nodes
        logical, allocatable :: FSlink(:)
        type(FS_port_pointer_t), allocatable :: thermP(:)
        type(SS_src_port_pointer_t), allocatable :: thermS(:)
        integer, allocatable :: idxFSlk(:) 
    end type solid_prop_t 

    type StateVariable_solid_t
        real(dp), allocatable :: MCtemp(:) !! Temperature of solid
    end type StateVariable_solid_t    

    type big_arrays_solid_t
        real(dp), allocatable :: amS(:) !! lower diagonal subMatrix
        real(dp), allocatable :: bmS(:) !! main diagonal subMatrix    
        real(dp), allocatable :: cmS(:) !! upper diagonal subMatrix

        real(dp), pointer     :: bvS(:)    
    end type big_arrays_solid_t     

    type arrays_linear_system_solid_t
        type(dbl_pointer_t), allocatable :: ValExp(:), ValImp(:) !! Non zero values in COO format for 1 solid
        integer, allocatable :: RowExp(:), RowImp(:)  !! row index array in COO format for 1 solid
        integer, allocatable :: ColExp(:), ColImp(:)  !! col index array in COO format for 1 solid
    end type arrays_linear_system_solid_t      

    type flux_solid_t
        real(dp), allocatable :: sd(:)   !! heat diffusion flux
        real(dp), allocatable :: sd_DerdTL(:), sd_DerdTR(:)
    end type flux_solid_t    

    type solid_t
        type(hdf_desc_t)                   :: hdf      
        type(StateVariable_solid_t)        :: StVar          !! StateVariable updated by BDF 2 scheme -> time n+1
        type(StateVariable_solid_t)        :: StVarOld       !! StateVariable -> time n
        type(StateVariable_solid_t)        :: StVarOld2      !! StateVariable -> time n-1
        type(StateVariable_solid_t)        :: StVarOld3      !! StateVariable -> time n-2
        type(StateVariable_solid_t)        :: StVarOld4      !! StateVariable -> time n-3
        type(StateVariable_solid_t)        :: StVarOld5      !! StateVariable -> time n-4
        type(StateVariable_solid_t)        :: StVarOld6      !! StateVariable -> time n-5
        type(solid_prop_t)                 :: MC_Prop        !! Variables relative to solid global properties
        type(big_arrays_solid_t)           :: big            !! Big arrays
        type(arrays_linear_system_solid_t) :: bigLS          !! Big arrays for Linear system
        type(signal_t)                     :: Q_ext_load     !! External heat load signal (W/m)
        real(dp), allocatable              :: Q_ext(:)       !! External heat load per cell (W/m)        
        type(flux_solid_t)                 :: flxS           !! Heat diffusion fluxes along spiral
        real(dp)                           :: err
        real(dp)                           :: err_den
        class(material_t), pointer         :: mat            !! Material dependent properties
    end type solid_t

contains

subroutine all_solid_allocation(me)
    type(solid_t), intent(inout) :: me
    integer :: NbCells
      
    NbCells=me%MC_Prop%NbCells
    allocate(me%StVar%MCtemp(NbCells))
    allocate(me%StVarOld%MCtemp(NbCells))
    allocate(me%StVarOld2%MCtemp(NbCells))
    allocate(me%StVarOld3%MCtemp(NbCells))
    allocate(me%StVarOld4%MCtemp(NbCells))
    allocate(me%StVarOld5%MCtemp(NbCells))
    allocate(me%StVarOld6%MCtemp(NbCells))
    allocate(me%flxS%sd(0:NbCells))
    allocate(me%flxS%sd_DerdTR(0:NbCells))
    allocate(me%flxS%sd_DerdTL(0:NbCells))    
    allocate(me%big%amS(NbCells-1))
    allocate(me%big%cmS(NbCells-1))
    allocate(me%big%bmS(NbCells))
    allocate(me%MC_Prop%idxFSlk(me%MC_Prop%NbCells))
    allocate(me%Q_ext(NbCells))
      
    me%StVar%MCtemp=0.0_dp; me%StVarOld%MCtemp=0.0_dp
    me%StVarOld2%MCtemp=0.0_dp; me%StVarOld3%MCtemp=0.0_dp
    me%StVarOld4%MCtemp=0.0_dp; me%StVarOld5%MCtemp=0.0_dp; me%StVarOld6%MCtemp=0.0_dp
    me%flxS%sd=0.0_dp
    me%flxS%sd_DerdTR=0.0_dp; me%flxS%sd_DerdTL=0.0_dp    
    me%big%amS=0.0_dp;me%big%bmS=0.0_dp;me%big%cmS=0.0_dp
    me%MC_Prop%idxFSlk=0; me%Q_ext=0.0_dp
end subroutine all_solid_allocation



subroutine solid_init_part1(me, krn, cfg, h5, sim)
    type(solid_t),         intent(out) :: me
    type(krn_t),          intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg
    class(hdf5_t),        intent(inout) :: h5
    type(simulation_t),   intent(inout) :: sim

    real(dp) :: T_init
    integer :: i, ii, nb_non_zeros_exp, nb_non_zeros_imp, nb_SS_src_ports
     
    if(cfg%has_key('nodes')) then ! means uniform mesh
        me%MC_Prop%NbCells = cfg%int('nodes')
        allocate(me%MC_Prop%dxLoc(me%MC_Prop%NbCells),me%MC_Prop%volLoc(me%MC_Prop%NbCells))
        me%MC_Prop%dxLoc(:)  = cfg%dbl('length')/me%MC_Prop%NbCells
        me%MC_Prop%volLoc(:) = cfg%dbl('volume')/me%MC_Prop%NbCells
    else ! variable mesh
        me%MC_Prop%NbCells=size(cfg%dbl1d('length'))
        allocate(me%MC_Prop%dxLoc(me%MC_Prop%NbCells),me%MC_Prop%volLoc(me%MC_Prop%NbCells))
        me%MC_Prop%dxLoc(:)  = cfg%dbl1d('length')
        me%MC_Prop%volLoc(:) = cfg%dbl1d('volume', n=me%MC_Prop%NbCells)
    endif

    me%MC_Prop%SSsrcLink = cfg%bin('solid_link')

    call material_init(me%mat,cfg)

    allocate(me%MC_Prop%FSlink(me%MC_Prop%NbCells))
    me%MC_Prop%FSlink(:) = cfg%bin1D('channel_link', n=me%MC_Prop%NbCells)
    
    ! For HDF5
    me%hdf%name   = cfg%str('id')
    me%hdf%node_x = [(sum(me%MC_Prop%dxLoc(1:i-1)) + me%MC_Prop%dxLoc(i)/2, &
                     i = 1, me%MC_Prop%NbCells)]
    call h5%add_to_table(me%hdf,'solid')

    call all_solid_allocation(me)
    
    call me%Q_ext_load%init(sim,cfg,'flux',me%hdf%node_x)

    T_init=cfg%dbl('initial/t')
    do ii=1,me%MC_Prop%NbCells
        me%StVar%MCtemp(ii)=T_init
    enddo

    me%MC_Prop%cond_btw_nodes=cfg%bin('conduction_between_nodes')
        
    nb_non_zeros_exp=me%MC_Prop%NbCells

    if(me%MC_Prop%cond_btw_nodes) then
        nb_non_zeros_imp=2
        do i=2,me%MC_Prop%NbCells-1
            nb_non_zeros_imp=nb_non_zeros_imp+3
        enddo
        nb_non_zeros_imp=nb_non_zeros_imp+2
    else
        nb_non_zeros_imp=me%MC_Prop%NbCells
    endif

    me%MC_Prop%nb_FS_ports=0
    do ii=1,me%MC_Prop%NbCells
        if(me%MC_Prop%FSlink(ii)) me%MC_Prop%nb_FS_ports=me%MC_Prop%nb_FS_ports+1
    enddo 

    if(me%MC_Prop%SSsrcLink) then
        nb_SS_src_ports=me%MC_Prop%NbCells
    else
        nb_SS_src_ports=0
    endif

    call krn%add(cfg%str('id'),me%MC_Prop%NbCells,nb_non_zeros_exp,nb_non_zeros_imp,0,me%MC_Prop%nb_FS_ports,0,0,nb_SS_src_ports)
    
end subroutine solid_init_part1


subroutine solid_init_part2(me, krn, cfg)
    !! Fluid port initialisation for solids
    type(solid_t), target, intent(inout) :: me
    type(krn_t),          intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg  

    integer :: NbSubM,i,idx,PrevRowIdx,PrevColIdx,j
    real(dp), pointer :: rhs(:)
    integer, allocatable :: list_thp_loc(:),list_thS_loc(:)
    type(FS_port_t), pointer :: FS_ports(:)
    type(SS_src_port_t), pointer :: SS_src_ports(:)

    call krn%update(rhs_or_solution_view=rhs)
    me%big%bvS(1:me%MC_Prop%NbCells) => rhs

    if(me%MC_Prop%nb_FS_ports>0) then
      allocate(list_thp_loc(me%MC_Prop%nb_FS_ports),me%MC_Prop%thermP(me%MC_Prop%nb_FS_ports))

      list_thp_loc = pack([(i, i=1, me%MC_Prop%NbCells)], me%MC_Prop%FSlink)
      call krn%update(FS_p_loc=list_thp_loc,CompName=cfg%str('id'),FS_p_view=FS_ports)
      
      j=1
      do i=1,me%MC_Prop%NbCells
          if(me%MC_Prop%FSlink(i)) then
              me%MC_Prop%idxFSlk(i)=j
              me%MC_Prop%thermP(j)%p => FS_ports(j)
              me%MC_Prop%thermP(j)%p%typeS   = 'solid'
              me%MC_Prop%thermP(j)%p%rhoMS   = me%mat%density
              me%MC_Prop%thermP(j)%p%VolMS   = me%MC_Prop%volLoc(i)
              j=j+1
          endif
      enddo
    endif

    if(me%MC_Prop%SSsrcLink) then
      allocate(list_thS_loc(me%MC_Prop%NbCells),me%MC_Prop%thermS(me%MC_Prop%NbCells))
      list_thS_loc=[(i,i=1,me%MC_Prop%NbCells)]
      call krn%update(SS_src_p_loc=list_thS_loc,SS_src_p_view=SS_src_ports)
      do i=1,me%MC_Prop%NbCells
          me%MC_Prop%thermS(i)%p => SS_src_ports(i)
          me%MC_Prop%thermS(i)%p%typeS   = 'solid'
          me%MC_Prop%thermS(i)%p%rhoS    = me%mat%density
          me%MC_Prop%thermS(i)%p%VolS    = me%MC_Prop%volLoc(i)
          me%MC_Prop%thermS(i)%p%Height  = me%MC_Prop%dxLoc(i)
      enddo
    endif    

    NbSubM=me%MC_Prop%NbCells
    allocate(me%bigLS%ValExp(NbSubM),me%bigLS%RowExp(NbSubM),me%bigLS%ColExp(NbSubM))

    if(me%MC_Prop%cond_btw_nodes) then
        NbSubM=2*2+(me%MC_Prop%NbCells-2)*3
    else
        NbSubM=me%MC_Prop%NbCells
    endif
    allocate(me%bigLS%ValImp(NbSubM),me%bigLS%RowImp(NbSubM),me%bigLS%ColImp(NbSubM))

    idx=1

    PrevColIdx=0
    PrevRowIdx=0
    do i=1,me%MC_Prop%NbCells
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

    if(me%MC_Prop%cond_btw_nodes) then
        PrevColIdx=1
        PrevRowIdx=0
        do i=1,me%MC_Prop%NbCells-1
            me%bigLS%ValImp(idx)%p => me%big%cmS(i)
            me%bigLS%ColImp(idx)=PrevColIdx+1
            me%bigLS%RowImp(idx)=PrevRowIdx+1          
            
            idx=idx+1
            PrevColIdx=PrevColIdx+1
            PrevRowIdx=PrevRowIdx+1      
        enddo
        PrevColIdx=0
        PrevRowIdx=1
        do i=1,me%MC_Prop%NbCells-1
            me%bigLS%ValImp(idx)%p => me%big%amS(i)
            me%bigLS%ColImp(idx)=PrevColIdx+1
            me%bigLS%RowImp(idx)=PrevRowIdx+1 
            
            idx=idx+1
            PrevColIdx=PrevColIdx+1
            PrevRowIdx=PrevRowIdx+1
        enddo
    endif
    call krn%coo_add(.false.,me%bigLS%ColImp,me%bigLS%RowImp,me%bigLS%ValImp)

    deallocate(me%bigLS%ValExp,me%bigLS%RowExp,me%bigLS%ColExp)
    deallocate(me%bigLS%ValImp,me%bigLS%RowImp,me%bigLS%ColImp)

end subroutine solid_init_part2

end module cmp_solid_init_m