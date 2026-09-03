! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_FFsrcLink_init_m
    use krn_interface_m
    use lib_input_m, only: input_t
    use lib_nusselt_correlations_m
    use krn_simulation_m
    implicit none

    type FFsrcLink_t
        integer :: sizeNodes        
        type(FF_src_port_pointer_t), allocatable :: lk(:,:)

        real(dp) :: WetPerChan
        class(nusselt_t), pointer :: nuss       !! Nusselt correlation 

        real(dp), allocatable :: DerSedOther(:,:,:), DerSrdOther(:,:,:)

        real(dp), allocatable :: Matrix(:,:,:,:)
        type(dbl_pointer_t), allocatable :: MxUp(:,:), MxDown(:,:)
    end type FFsrcLink_t

contains

subroutine all_FFsrcLink_allocation(me)
    type(FFsrcLink_t), intent(inout) :: me

    allocate(me%lk(2,me%sizeNodes))
    
    allocate(me%DerSedOther(Nb_VarC,2,me%sizeNodes),me%DerSrdOther(Nb_VarC,2,me%sizeNodes))
    allocate(me%Matrix(Nb_VarC,2,2,me%sizeNodes))
    allocate(me%MxUp(8,me%sizeNodes),me%MxDown(8,me%sizeNodes))
    me%DerSedOther=0.0_dp; me%DerSrdOther=0.0_dp; me%Matrix=0.0_dp

end subroutine all_FFsrcLink_allocation


subroutine FFsrcLink_init_part1(me,krn,cfg)
    type(FFsrcLink_t), intent(out) :: me
    type(krn_t), intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg

    type(input_t), allocatable :: links(:)
    integer :: nb_non_zeros_exp, nb_non_zeros_imp

    links = cfg%dict1d('link')
    me%sizeNodes=size(links(1)%int1d('node'))

    nb_non_zeros_exp = 2*(2*Nb_VarC)*me%sizeNodes ! Only thermal exchanges between pipes inside the matrix - Only remains pressure relaxation after Pardiso call
    nb_non_zeros_imp = 2*(2*Nb_VarC)*me%sizeNodes ! Only thermal exchanges between pipes inside the matrix - Only remains pressure relaxation after Pardiso call

    call krn%add('',0,nb_non_zeros_exp,nb_non_zeros_imp,0,0,0,0,0)
end subroutine FFsrcLink_init_part1


subroutine FFsrcLink_init_part2(me,krn,cfg)
    type(FFsrcLink_t), intent(inout) :: me
    type(krn_t), target, intent(inout) :: krn
    class(input_t), pointer, intent(in)    :: cfg

    type(input_t), allocatable :: links(:)
    integer, allocatable :: nodes1(:), nodes2(:)
    integer :: i,offset1,offset2

    call all_FFsrcLink_allocation(me)

    links = cfg%dict1d('link')
    nodes1=links(1)%int1d('node')
    nodes2=links(2)%int1d('node')
    offset1 = krn%FF_src_list%find(links(1)%str('id'), 1) - 1
    offset2 = krn%FF_src_list%find(links(2)%str('id'), 1) - 1
    do i = 1, me%sizeNodes
        me%lk(1,i)%p => krn%FF_src_ports(offset1 + nodes1(i))
        me%lk(2,i)%p => krn%FF_src_ports(offset2 + nodes2(i))
    enddo

    me%WetPerChan = cfg%dbl('wetted_perim_channels')

    call nusselt_init(me%nuss,cfg) ! nusselt required in the link       
    if(.NOT.associated(me%nuss)) error stop 'Nusselt does not exist in FluidLink'     

    do i = 1, me%sizeNodes
      call krn%coo_add_link(me%lk(1,i)%p, me%lk(2,i)%p, idx_2x4_col, idx_2x4_row, to_ptr_arr(me%Matrix(:,:,1,i)), me%MxUp(:,i))
      call krn%coo_add_link(me%lk(2,i)%p, me%lk(1,i)%p, idx_2x4_col, idx_2x4_row, to_ptr_arr(me%Matrix(:,:,2,i)), me%MxDown(:,i))
    enddo
    
end subroutine FFsrcLink_init_part2

end module cmp_FFsrcLink_init_m