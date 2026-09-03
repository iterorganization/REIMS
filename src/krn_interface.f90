! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

include 'mkl_spblas.f90'
include 'mkl_pardiso.f90'
module krn_interface_m
    !! Module used by kernel and components to assure communication between them
    !!
    !! This module should be used (include) by the components at the initiation phase to:
    !!
    !!  1. Set global matrix sizes required by each component (state owners and links)
    !!  2. Expose state owners existing ports and find them by links
    !!  3. Set positions and variables in global matrices and converting them to pardiso format
    !!
    !! This module should be used (include) by kernel to:
    !!
    !!  1. Allocate array of components and full matrices for linear solver at initiation phase
    !!  2. Provide arrays to linear solver including translation array of pointers to array of doubles
    !!
    !! Sequence of execution in 2 stages of initialisation:
    !!
    !!  1. Stage 1 -- Kernel -> Initial allocation with number of components
    !!  2. Stage 1 -- Components -> Reports required sizes and ports
    !!  3. Stage 2 -- Kernel -> Allocation of the main matrices
    !!  4. Stage 2 -- Components (state owners) -> Set positions of their ports and receives back
    !!                   allocated ports and arrays like rhs or solution
    !!  5. Stage 2 -- Components (links) -> Get requested ports
    !!  6. Stage 2 -- Components -> Set non zero elements using they own coordinate system
    !!  7. Stage 2 -- Kernel -> Translate coo matrices to csr matrices
    use krn_global_tools_m
    use mkl_pardiso
    implicit none

    type dbl_pointer_t
        !! encapsulate pointer in order to create array of pointers
        real(dp), pointer :: p !! pointer to double
    end type dbl_pointer_t

    type coo_csr_t(nnz)
        integer, len :: nnz     !! number of non-zeros (major size of the type)
        integer :: idx = 1      !! tracing last filled index of global array
        integer :: coo_col(nnz) !! column idx of non zero array in coo
        integer :: coo_row(nnz) !! row idx of non zero array in coo
        type(dbl_pointer_t) :: coo_ptr(nnz) !! non zero length array of pointers to physical value in coo
        integer, pointer :: csr_col(:)      !! for pardiso final column idx of non zero array in csr
        integer, pointer :: csr_row(:)      !! for pardiso final row idx of non zero array in csr
        type(dbl_pointer_t), pointer :: csr_ptr(:) !! for pardiso final non zero length array of pointers to physical value in csr
        real(dp) ::  csr_val(nnz)           !! non zero value array for pardiso
        type(mkl_pardiso_handle) :: pt(64)  !! internal pointers of pardiso
        integer :: iparm(64) = 0            !! parameters for pardiso
        integer :: perm(nnz)                !! non zero size after analysis for pardiso
        logical :: analysed = .false.       !! is pardiso analysed?
    end type coo_csr_t

    type, abstract :: port_t
        integer :: loc !! top left location of port matrix
        logical :: can_be_shared = .false.
    end type port_t

    type :: port_list_t
        integer :: idx = 1                !! Tracking new entries
        integer, allocatable :: nb(:)     !! Number of ports per component
        character(:), allocatable :: list !! contains entries like this: <space><component_name>:<number_of_ports>,<starting_idx>
    contains
        procedure :: add          => port_list_add           !! Used by State components exposing ports here
        procedure :: update       => port_list_update        !! Used by krn_update to calculate offsets etc...
        procedure :: find         => port_list_find          !! Finding index of the port in the ports list by links components
        procedure :: find_nbPorts => port_list_find_nbPorts  !! Finding number of ports for a given label in the ports list by links components
    end type port_list_t

    type krn_t
        !! Structure which will be used to size matrices like: coo,csr and ports 
        integer :: idx = 1
        integer, allocatable :: nb_state_vars(:)
        integer, allocatable :: nb_non_zeros_exp(:)
        integer, allocatable :: nb_non_zeros_imp(:)

        type(port_list_t) :: FF_flux_list, FF_src_list, FS_list, SS_flux_list, SS_src_list
        character(:), allocatable :: PortIdxForFSLink_List

        type(coo_csr_t(:)), allocatable :: imp      !! Implicit matrices in both format coo and csr
        type(coo_csr_t(:)), allocatable :: exp      !! Explicit matrices in both format coo and csr

        real(dp), allocatable :: rhs_or_solution(:) !! Right hand side of linear equation or it's solution
        real(dp), allocatable :: solution(:)        !! solution vector of the linear system
        type(FF_flux_port_t),    allocatable :: FF_flux_ports(:)      !! Fluid/Fluid port array --> flux
        type(FF_src_port_t),     allocatable :: FF_src_ports(:)       !! Fluid/Fluid port array --> source        
        type(FS_port_t),         allocatable :: FS_ports(:)           !! Fluid/Solid port array
        type(SS_flux_port_t),    allocatable :: SS_flux_ports(:)      !! Solid/Solid port array --> flux
        type(SS_src_port_t),     allocatable :: SS_src_ports(:)       !! Solid/Solid port array --> source
        ! TODO should be here also state space 3 vectors?

    contains
        procedure :: init_part1   => krn_init_part1   !! Allocating kernel arrays with number components
        procedure :: add          => krn_add          !! Used by all components add their array and port requirements 
        procedure :: add_label    => krn_add_label    !! Add label of 2d mesh to expose ports
        procedure :: init_part2   => krn_init_part2   !! Allocating full matrices
                
        procedure :: update       => krn_update       !! Used by state owners to set ports locations and get back link to them
        procedure :: coo_add      => krn_coo_add      !! Used by state owners to add coordinate format array to the full array
        procedure :: coo_add_link => krn_coo_add_link !! Used by links to add coordinate format array to the full array
        ! before it was also get function now is used by find index
        procedure :: coo_to_csr   => krn_coo_to_csr   !! Conversion full coordinate format array to sparse compressed row array
    end type krn_t



    type, extends(port_t) :: FF_flux_port_t
        real(dp) :: Prim(Nb_VarP)                  !! Channel->link, primitive variables
        real(dp) :: Cons(Nb_VarC)                  !! Channel->link, conservative variables
        real(dp) :: Su_fric                        !! Channel->link, friction source term
        real(dp) :: DerSu_fric(Nb_VarC)            !! Channel->link, derivatives of friction source term  
        real(dp) :: Sr_fric                        !! Channel->link, friction source term (4th eq)
        real(dp) :: DerSr_fric(Nb_VarC)            !! Channel->link, derivatives of friction source term (4th eq)                                            
        real(dp) :: Area                           !! Channel->link, cross-section area of the pipe
        real(dp) :: dxLoc                          !! Channel->link, space step of the cell close to the boundary
        real(dp) :: Sgn4j                          !! Channel->link, sign used by junctions to identify incoming branches
        real(dp) :: flx(Nb_VarC)                   !! Link->channel, flux calculated by link 
        real(dp) :: vit                            !! Link->channel, fluid velocity
        real(dp) :: derFlx_derCon(Nb_VarC,Nb_VarC) !! Link->channel, partial derivative of flux calculated by link
        real(dp) :: derVit_derCon(Nb_VarC)         !! Link->channel, derivative of velocity for every variable
        logical  :: connected                      !! Link connected to a port or not
    end type FF_flux_port_t   

    type, extends(port_t) :: FF_src_port_t
        ! Pipe -> Link
        real(dp) :: Pra, Re, Phi, llambda, TempP, AreaP, DiamP
        real(dp) :: dPrdUi(Nb_VarC), dRedUi(Nb_VarC), DerPhi(Nb_VarC)
        real(dp) :: DerllambdadUi(Nb_VarC), DerTempPdUi(Nb_VarC)
        ! Link -> Pipe
        real(dp) :: rhs_Se,rhs_Sr
        real(dp) :: DerSe(Nb_VarC),DerSr(Nb_VarC)
        ! Pipe -> Link -> Pipe
        real(dp) :: Prim(Nb_VarP)
        real(dp) :: Cons(Nb_VarC)
    end type FF_src_port_t  

    type FF_src_port_pointer_t
        type(FF_src_port_t), pointer :: p
    end type FF_src_port_pointer_t     

    type, extends(port_t) :: FS_port_t
        character(6) :: typeS ! strand or solid
        ! Pipe -> Link
        real(dp) :: Pra, Re, Phi, llambda, TempP, AreaP, DiamP, dxLoc
        real(dp) :: dPrdUi(Nb_VarC), dRedUi(Nb_VarC), DerPhi(Nb_VarC)
        real(dp) :: DerllambdadUi(Nb_VarC), DerTempPdUi(Nb_VarC)
        ! Strand -> Link
        real(dp) :: TempS, AreaS, rhoMS, cpMS, dcpMSdT, Dist4Grad
        ! solid -> Link (same as the ones from strand + the following)
        real(dp) :: VolMS,lambdS,DerlambdS
        ! mesh2D -> Link (same as the previous ones + the following)
        real(dp) ::  LengthCont, cpS, dcpSdT, rhoS, SurfS
        ! Link -> Pipe/(Strand, solid or mesh2D)
        real(dp) :: rhs_SrcP_Se, rhs_SrcP_Sr, rhs_SrcS, DerSrcSdT
        real(dp) :: DerSrcP_Se(Nb_VarC), DerSrcP_Sr(Nb_VarC)  
        real(dp) :: Flx, derFlx_derCon
    end type FS_port_t
       
    type FS_port_pointer_t
        type(FS_port_t), pointer :: p
    end type FS_port_pointer_t       
    
    type, extends(port_t) :: SS_flux_port_t
      ! Strand -> Link
      real(dp) :: TempS, lambdS, DerlambdS, AreaS, Sgn4lk, dxLoc, rhoMS, cpMs
      ! Link -> Strand
      real(dp) :: Flx, derFlx_derCon
    end type SS_flux_port_t

    type SS_flux_port_pointer_t
        type(SS_flux_port_t), pointer :: p
    end type SS_flux_port_pointer_t      
    
    type, extends(port_t) :: SS_src_port_t
      character(6) :: typeS ! solid or mesh2D
      ! Solid -> Link
      real(dp) :: TempS, rhoS, VolS, Dist4Grad, LengthCont, Height, SurfS, cpS, dcpSdT, lambdS, DerlambdS
      ! Link -> Solid
      real(dp) :: rhs_SrcS, DerSrcSdT
      real(dp) :: Flx, derFlx_derCon
    end type SS_src_port_t  

    type SS_src_port_pointer_t
        type(SS_src_port_t), pointer :: p
    end type SS_src_port_pointer_t       

    interface to_ptr_arr !! Translation array 0d, 1d or 2d of doubles to array of pointers
        module procedure to_ptr_arr0d, to_ptr_arr1d, to_ptr_arr2d
    end interface
contains

subroutine port_list_add(me,idx,name,nb_ports)
    class(port_list_t), intent(inout) :: me
    integer, intent(in) :: idx        !!
    integer, intent(in) :: nb_ports   
    character(*), intent(in) :: name  !! Component (port) name

    me%nb(idx) = me%nb(idx) + nb_ports ! WARNING why to add if no of ports is equal 0 or the name is not given?
    if (name/='' .and. nb_ports/=0) then
        me%list = me%list//' '//name//':'//to_str(nb_ports)//','//to_str(me%idx)
        me%idx = me%idx + nb_ports
    endif
end subroutine port_list_add

function port_list_update(me, idx, loc, start_state_idx, ports_array) result(view)
    class(port_list_t), intent(inout) :: me
    integer, intent(in) :: idx !! 
    integer, intent(in) :: loc(:) !!
    integer, intent(in) :: start_state_idx !!  
    class(port_t), target, intent(in) :: ports_array(:) !! which will have updated loc (should be inout?)
    class(port_t), pointer :: view(:) !! Returned view on the ports array with updated values of: loc

    integer :: start_idx

    start_idx = sum(me%nb(1:idx)) - me%nb(idx)
    view => ports_array(start_idx+1:start_idx+me%nb(idx))
    view%loc = loc + start_state_idx
end function port_list_update

function port_list_find(me,name,idx) result(idx_out)
    class(port_list_t), target, intent(in) :: me
    character(*), intent(in) :: name !! Name of the component
    integer, intent(in) :: idx       !! Port number
    integer :: idx_out               !! Index of the port on ports array
    integer :: nbPorts, offset, i, pos

    pos = index(me%list,' '//name//':')
    if (pos == 0) then
        print*, 'ERROR: component "'//name//'" not found in port list.'
        error stop 'Aborting: unknown component ID in connection.'
    endif
    i = pos + len(name) + 2
    nbPorts = to_int(me%list(i:))
    if (idx > nbPorts) then
        print*, 'ERROR: port index '//to_str(idx)//' exceeds nb ports ('//to_str(nbPorts)//') for "'//name//'"'
        error stop 'Aborting: port index out of range.'
    endif
    offset = to_int(me%list(i+index(me%list(i:),','):))
    idx_out = idx + offset - 1
end function port_list_find

function port_list_find_nbPorts(me,name) result(idx_out)
    class(port_list_t), target, intent(in) :: me
    character(*), intent(in) :: name !! Name of the component
    integer :: idx_out               !! Index of the port on ports array
    integer :: i, pos

    pos = index(me%list,' '//name//':')
    if (pos == 0) then
        print*, 'ERROR: component "'//name//'" not found in port list.'
        error stop 'Aborting: unknown component ID in connection.'
    endif
    i = pos + len(name) + 2
    idx_out = to_int(me%list(i:))

end function port_list_find_nbPorts

subroutine krn_init_part1(me,size)
    class(krn_t), intent(inout) :: me
    integer, intent(in) :: size
    allocate(me%nb_state_vars(size))
    allocate(me%nb_non_zeros_exp(size))
    allocate(me%nb_non_zeros_imp(size))

    me%FF_flux_list%list = ''
    me%FF_src_list%list = ''
    me%FS_list%list = ''
    me%SS_flux_list%list = ''
    me%SS_src_list%list = ''
    allocate(me%FF_flux_list%nb(size))
    allocate(me%FF_src_list%nb(size))
    allocate(me%FS_list%nb(size))
    allocate(me%SS_flux_list%nb(size))
    allocate(me%SS_src_list%nb(size))
    me%FF_flux_list%nb = 0
    me%FF_src_list%nb = 0
    me%FS_list%nb = 0
    me%SS_flux_list%nb = 0
    me%SS_src_list%nb = 0
    me%PortIdxForFSLink_List = ''

end subroutine krn_init_part1

subroutine krn_init_part2(me)
    class(krn_t), intent(inout) :: me
    integer :: imp_size,exp_size,rhs_size
    imp_size = sum(me%nb_non_zeros_imp)
    exp_size = sum(me%nb_non_zeros_exp)
    rhs_size = sum(me%nb_state_vars)
    allocate(me%FF_flux_ports(sum(me%FF_flux_list%nb)))
    allocate(me%FS_ports(sum(me%FS_list%nb)))
    allocate(me%FF_src_ports(sum(me%FF_src_list%nb)))
    allocate(me%SS_flux_ports(sum(me%SS_flux_list%nb)))
    allocate(me%SS_src_ports(sum(me%SS_src_list%nb)))
    allocate(me%rhs_or_solution(rhs_size))
    allocate(me%solution(rhs_size))

    allocate(coo_csr_t(imp_size)::me%imp)
    allocate(coo_csr_t(exp_size)::me%exp)
    me%imp%perm = 0
    me%exp%perm = 0
    me%idx = 0
end subroutine krn_init_part2

subroutine krn_update(me, rhs_or_solution_view, FF_flux_p_loc, FS_p_loc, &
        FF_src_p_loc, SS_flux_p_loc, SS_src_p_loc, CompName, FF_flux_p_view, FS_p_view, &
        FF_src_p_view, SS_flux_p_view, SS_src_p_view)
    class(krn_t), target, intent(inout) :: me
    real(dp), pointer, intent(out), optional :: rhs_or_solution_view(:) !! Returns pointer with rhs or solution array
    integer, intent(in), optional :: FF_flux_p_loc(:) !! Position (state index) of the FF_flux ports local to the component
    integer, intent(in), optional :: FS_p_loc(:) !! Position (state index) of the FS ports local to the component
    integer, intent(in), optional :: FF_src_p_loc(:) !! Position (state index) of the FF_src ports local to the component
    integer, intent(in), optional :: SS_flux_p_loc(:) !! Position (state index) of the SS_flux ports local to the component   
    integer, intent(in), optional :: SS_src_p_loc(:) !! Position (state index) of the SS_src ports local to the component  
    character(*), intent(in), optional :: CompName 
    type(FF_flux_port_t), pointer, intent(out), optional :: FF_flux_p_view(:) !! Returns pointer with array of FF_flux ports for this component
    type(FS_port_t), pointer, intent(out), optional :: FS_p_view(:) !! Returns pointer with array of FS ports for this component
    type(FF_src_port_t), pointer, intent(out), optional :: FF_src_p_view(:) !! Returns pointer with array of FF_src ports for this component
    type(SS_flux_port_t), pointer, intent(out), optional :: SS_flux_p_view(:) !! Returns pointer with array of SS_flux ports for this component
    type(SS_src_port_t), pointer, intent(out), optional :: SS_src_p_view(:) !! Returns pointer with array of SS_src ports for this component

    integer :: start_state_idx,ii
    logical :: all_equal
    class(port_t), pointer :: port_ptr(:)

    if(present(rhs_or_solution_view)) me%idx = me%idx + 1
    start_state_idx = sum(me%nb_state_vars(1:me%idx)) - me%nb_state_vars(me%idx)
    if(present(rhs_or_solution_view)) then
        rhs_or_solution_view => me%rhs_or_solution(start_state_idx+1:start_state_idx+me%nb_state_vars(me%idx))
    endif
    if(present(FF_flux_p_loc) .and. present(FF_flux_p_view)) then
        port_ptr => me%FF_flux_list%update(me%idx,FF_flux_p_loc,start_state_idx,me%FF_flux_ports)
        select type (port_ptr); class is (FF_flux_port_t)
            FF_flux_p_view => port_ptr
        end select
    endif
    if(present(FS_p_loc) .and. present(FS_p_view)) then
        port_ptr => me%FS_list%update(me%idx,FS_p_loc,start_state_idx,me%FS_ports)
        select type (port_ptr); class is (FS_port_t)
            FS_p_view => port_ptr
        end select
        if(present(CompName)) then
            all_equal = all(FS_p_loc == [(ii, ii = 1, size(FS_p_loc))])
            if(all_equal) then ! case of strands or solid if uniform channel_link array
                me%PortIdxForFSLink_List = me%PortIdxForFSLink_List//' '//CompName//':0'
            else ! case of solids only
                me%PortIdxForFSLink_List = me%PortIdxForFSLink_List//' '//CompName//':'
                do ii=1,size(FS_p_loc)-1
                    me%PortIdxForFSLink_List = me%PortIdxForFSLink_List//to_str(FS_p_loc(ii))//','
                enddo
                me%PortIdxForFSLink_List = me%PortIdxForFSLink_List//to_str(FS_p_loc(size(FS_p_loc)))
            endif
        endif
    endif
    if(present(FF_src_p_loc) .and. present(FF_src_p_view)) then
        port_ptr => me%FF_src_list%update(me%idx,FF_src_p_loc,start_state_idx,me%FF_src_ports)
        select type (port_ptr); class is (FF_src_port_t)
            FF_src_p_view => port_ptr
        end select
    endif     
    if(present(SS_flux_p_loc) .and. present(SS_flux_p_view)) then
        port_ptr => me%SS_flux_list%update(me%idx,SS_flux_p_loc,start_state_idx,me%SS_flux_ports)
        select type (port_ptr); class is (SS_flux_port_t)
            SS_flux_p_view => port_ptr
        end select
    endif         
    if(present(SS_src_p_loc) .and. present(SS_src_p_view)) then
        port_ptr => me%SS_src_list%update(me%idx,SS_src_p_loc,start_state_idx,me%SS_src_ports)
        select type (port_ptr); class is (SS_src_port_t)
            SS_src_p_view => port_ptr
        end select
    endif     
end subroutine krn_update

subroutine krn_add(me, name, nb_state_vars, nb_non_zeros_exp, nb_non_zeros_imp, &
        nb_FF_flux_ports, nb_FS_ports, nb_FF_src_ports, nb_SS_flux_ports, &
        nb_SS_src_ports)
    class(krn_t), intent(inout) :: me
    character(*), intent(in) :: name
    integer, intent(in) :: nb_state_vars
    integer, intent(in) :: nb_non_zeros_exp
    integer, intent(in) :: nb_non_zeros_imp
    integer, intent(in) :: nb_FF_flux_ports
    integer, intent(in) :: nb_FS_ports
    integer, intent(in) :: nb_FF_src_ports
    integer, intent(in) :: nb_SS_flux_ports
    integer, intent(in) :: nb_SS_src_ports
    
    me%nb_state_vars(me%idx) = nb_state_vars
    me%nb_non_zeros_exp(me%idx) = nb_non_zeros_exp
    me%nb_non_zeros_imp(me%idx) = nb_non_zeros_imp
 
    call me%FF_flux_list%add(me%idx,name,nb_FF_flux_ports)
    call me%FS_list%add(me%idx,name,nb_FS_ports)
    call me%FF_src_list%add(me%idx,name,nb_FF_src_ports)
    call me%SS_flux_list%add(me%idx,name,nb_SS_flux_ports)
    call me%SS_src_list%add(me%idx,name,nb_SS_src_ports)
    me%idx = me%idx + 1
end subroutine krn_add

subroutine krn_add_label(me, name_plus_label, nb_nodes_FS, nb_nodes_SS)
    class(krn_t), intent(inout) :: me
    character(*), intent(in) :: name_plus_label
    integer, intent(in), optional :: nb_nodes_FS
    integer, intent(in), optional :: nb_nodes_SS

    if(present(nb_nodes_FS)) call me%FS_list%add(me%idx-1,name_plus_label,nb_nodes_FS)
    if(present(nb_nodes_SS)) call me%SS_src_list%add(me%idx-1,name_plus_label,nb_nodes_SS)
end subroutine krn_add_label


subroutine krn_coo_add(me,for_explicit,col,row,val)
    class(krn_t), target, intent(inout) :: me
    logical, intent(in) :: for_explicit !! if true add to explicit arrays otherwise to implicit
    integer, intent(in) :: col(:), row(:) !! relative indexes of the variables placed into global matrix
    type(dbl_pointer_t), intent(in) :: val(:) !! pointer to the value which should be placed into global matrix
    integer :: loc !! starting location for the first state variable
    type(coo_csr_t(:)), pointer :: arr

    loc = sum(me%nb_state_vars(1:me%idx))-me%nb_state_vars(me%idx)

    arr => me%imp
    if(for_explicit) arr => me%exp
    arr%coo_col(arr%idx:arr%idx+size(col)-1) = col + loc
    arr%coo_row(arr%idx:arr%idx+size(col)-1) = row + loc
    arr%coo_ptr(arr%idx:arr%idx+size(col)-1) = val
    arr%idx = arr%idx + size(val)
end subroutine krn_coo_add

subroutine krn_coo_add_link(me,port1,port2,col,row,val,val_out)
    class(krn_t), intent(inout) :: me
    class(port_t), intent(in) :: port1, port2
    integer, intent(in) :: col(:), row(:) !! relative indexes of the variables placed into global matrix
    type(dbl_pointer_t), intent(in) :: val(:) !! pointer to the value which should be placed into global matrix
    type(dbl_pointer_t), optional, intent(out) :: val_out(:) !! in case of overlapping matrix
 
    integer :: i
    type(FF_flux_port_t) :: FF_flux_port
    type(SS_flux_port_t) :: SS_flux_port

    if (.not.same_type_as(port1,FF_flux_port) .and. .not.same_type_as(port1,SS_flux_port)) then
        me%exp%coo_col(me%exp%idx:me%exp%idx+size(col)-1) = col + (port2%loc - 1)
        me%exp%coo_row(me%exp%idx:me%exp%idx+size(col)-1) = row + (port1%loc - 1)
        me%exp%coo_ptr(me%exp%idx:me%exp%idx+size(col)-1) = val
        me%exp%idx = me%exp%idx + size(val)
    end if

    if (present(val_out)) then
        val_out(1)%p => null()
        if (port1%can_be_shared .and. port2%can_be_shared) then
            do i=1,me%imp%idx
                if( me%imp%coo_col(i)==col(1)+(port2%loc-1) .and. &
                    me%imp%coo_row(i)==row(1)+(port1%loc-1)) then
                    val_out(:) = me%imp%coo_ptr(i:i+size(col)-1)
                    return
                endif
            enddo
        endif
    endif
    me%imp%coo_col(me%imp%idx:me%imp%idx+size(col)-1) = col + (port2%loc - 1)
    me%imp%coo_row(me%imp%idx:me%imp%idx+size(col)-1) = row + (port1%loc - 1)
    me%imp%coo_ptr(me%imp%idx:me%imp%idx+size(col)-1) = val
    me%imp%idx = me%imp%idx + size(val)
end subroutine krn_coo_add_link

subroutine accumulate(Mat2add,Mx)
    real(dp), intent(in) :: Mat2add(:,:)
    type(dbl_pointer_t), intent(inout) :: Mx(:)
    integer :: i,j,k

    k=1
    do i = 1,size(Mat2add,2)
        do j = 1,size(Mat2add,1)
            Mx(k)%p = Mx(k)%p + Mat2add(j,i)
            k=k+1
        enddo
    enddo
end subroutine accumulate

subroutine krn_coo_to_csr(me)
    class(krn_t), intent(inout) :: me
    integer :: nb_state_vars,i
    type(coo_csr_t(:)), allocatable :: temp
    nb_state_vars = sum(me%nb_state_vars)
    call coo_to_csr(nb_state_vars, me%exp)

    if(size(me%imp%coo_col)>me%imp%idx-1) then
        allocate(coo_csr_t(me%imp%idx-1)::temp) 
        do i = 1, me%imp%idx-1
          temp%coo_col(i) = me%imp%coo_col(i)
          temp%coo_row(i) = me%imp%coo_row(i)
          temp%coo_ptr(i) = me%imp%coo_ptr(i)

          temp%csr_val(i) = me%imp%csr_val(i)
          temp%perm(i)    = me%imp%perm(i)
        enddo          

        call move_alloc(temp,me%imp)
    endif

    call coo_to_csr(nb_state_vars, me%imp)
end subroutine krn_coo_to_csr

subroutine coo_to_csr(nb_state_vars,arr)
    !! Using Intel sblas functions to translate coo to csr matrix
    !!
    !! If you want to replace this you can start by checking sorting algorithms:
    !! https://www.mjr19.org.uk/IT/sorts/
    use mkl_spblas
    use iso_c_binding
    integer, intent(in) :: nb_state_vars
    type(coo_csr_t(:)), allocatable, target, intent(inout) :: arr

    real(dp), pointer :: coo_ptr_mold(:)
    type(sparse_matrix_t) :: coo_mkl, csr_mkl
    integer :: stat, nbRow, nbCol, index_base, nbNonZeros
    type(c_ptr) :: row_cptr, col_cptr, val_cptr, dummy_cptr, coo_cptr

    nbNonZeros = size(arr%coo_col)
    ! pretend that array of pointers is array of doubles
    coo_cptr = c_loc(arr%coo_ptr)
    call c_f_pointer(coo_cptr,coo_ptr_mold,[nbNonZeros])

    stat = mkl_sparse_d_create_coo(coo_mkl, sparse_index_base_one, nb_state_vars, &
        nb_state_vars, nbNonZeros, arr%coo_row, arr%coo_col, coo_ptr_mold)
    call check_status(stat,'In mkl_sparse_d_create_coo error returned: ')
    stat = mkl_sparse_convert_csr(coo_mkl, sparse_operation_non_transpose, csr_mkl)
    call check_status(stat,'In mkl_sparse_convert_csr error returned: ')
    stat = mkl_sparse_d_export_csr(csr_mkl, index_base, nbRow, nbCol, row_cptr, &
        dummy_cptr, col_cptr, val_cptr)
    call check_status(stat,'In mkl_sparse_d_export_csr error returned: ')

    call c_f_pointer(row_cptr,arr%csr_row,[nb_state_vars+1]) ! last value is number of non zeros + 1
    call c_f_pointer(col_cptr,arr%csr_col,[nbNonZeros])
    call c_f_pointer(val_cptr,arr%csr_ptr,[nbNonZeros])
end subroutine coo_to_csr

function to_ptr_arr0d(arr) result(ptr_arr)
    real(dp), target, intent(in) :: arr
    type(dbl_pointer_t) :: ptr_arr(1)
    ptr_arr(1)%p => arr
end function to_ptr_arr0d

function to_ptr_arr1d(arr) result(ptr_arr)
    real(dp), target, intent(in) :: arr(:)
    type(dbl_pointer_t) :: ptr_arr(size(arr))
    integer :: i
    do i = 1,size(arr)
        ptr_arr(i)%p => arr(i)
    enddo
end function to_ptr_arr1d

function to_ptr_arr2d(arr) result(ptr_arr)
    real(dp), target, intent(in) :: arr(:,:)
    type(dbl_pointer_t) :: ptr_arr(size(arr))
    !real(dp), pointer :: arr_ptr(:) ! TODO after updating compiler by Damien
    integer :: i,j,k
    k=1
    !arr_ptr(1:size(arr)) => arr
    do i = 1,size(arr,2)
        do j = 1,size(arr,1)
            ptr_arr(k)%p => arr(j,i)
            k=k+1
        enddo
    enddo
end function to_ptr_arr2d

function from_ptr_arr(ptr_arr) result(arr)
    type(dbl_pointer_t), intent(inout) :: ptr_arr(:)
    real(dp), allocatable :: arr(:)
    integer :: i

    allocate(arr(size(ptr_arr)))
    do i = 1,size(ptr_arr)
        arr(i) = ptr_arr(i)%p
    enddo
end function from_ptr_arr

end module krn_interface_m