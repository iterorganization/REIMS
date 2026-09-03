! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_FSlink_init_m
    use krn_interface_m
    use lib_input_m, only: input_t
    use lib_material_m    
    use lib_nusselt_correlations_m
    implicit none

    type FSlink_t   
        integer :: sizeNodes 
        type(FS_port_pointer_t), allocatable :: lk(:,:)

        integer, allocatable  :: nb_2D_Ports(:)

        real(dp) :: WetPer, Ht
        class(nusselt_t), pointer :: nuss       !! Nusselt correlation        

        real(dp), allocatable :: DerSrcP_SedT(:), DerSrcP_SrdT(:), SurfCont(:), DerSrcP_SedT2D(:,:), DerSrcP_SrdT2D(:,:)
        real(dp), allocatable :: DerSrcSdUi(:,:), derFlx2DdUi(:,:,:)

        real(dp), allocatable :: MatrixPS(:,:), MatrixSP(:,:), MatrixPS_2D(:,:,:), MatrixSP_2D(:,:,:)
    end type FSlink_t

contains


subroutine FSlink_init_part1(me,krn,cfg)
    type(FSlink_t), intent(out) :: me
    type(krn_t), intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg

    type(input_t), allocatable :: links(:)
    type(str_ptr), allocatable :: labels(:)
    integer :: nb_non_zeros_exp, nb_non_zeros_imp, i

    links = cfg%dict1d('link')
    me%sizeNodes=size(links(1)%int1d('node'))

    if(links(2)%has_key('node')) then 
        nb_non_zeros_exp = (Nb_VarC+2)*me%sizeNodes
        nb_non_zeros_imp = (Nb_VarC+2)*me%sizeNodes
    else if(links(2)%has_key('label')) then
        allocate(me%nb_2D_Ports(me%sizeNodes))
        labels=links(2)%str1d('label')
        do i = 1, me%sizeNodes
          me%nb_2D_Ports(i) = krn%FS_list%find_nbPorts(links(2)%str('id')//"+"//trim(labels(i)%p))
        enddo
        nb_non_zeros_exp = (Nb_VarC+2)*sum(me%nb_2D_Ports(:))
        nb_non_zeros_imp = (Nb_VarC+2)*sum(me%nb_2D_Ports(:))
    else
        error stop 'thermalink: link[2] requires "node" or "label"'
    endif

    call krn%add('',0,nb_non_zeros_exp,nb_non_zeros_imp,0,0,0,0,0)
end subroutine FSlink_init_part1


subroutine FSlink_init_part2(me,krn,cfg)
    type(FSlink_t), intent(inout) :: me
    type(krn_t), target, intent(inout) :: krn
    class(input_t), pointer, intent(in)    :: cfg

    type(input_t), allocatable :: links(:)
    real(dp), allocatable :: distances(:), surfaces(:)
    integer, allocatable :: nodes1(:), nodes2(:), conv(:)
    type(str_ptr), allocatable :: labels(:)
    integer :: i,offset1,offset2,max_value,j
    character(:), allocatable :: name2
    logical :: conversion_required

    links = cfg%dict1d('link')
    nodes1=links(1)%int1d('node')
    offset1 = krn%FS_list%find(links(1)%str('id'), 1) - 1

    if(links(2)%has_key('node')) then
        allocate(me%lk(2,me%sizeNodes))
        allocate(me%DerSrcP_SedT(me%sizeNodes),me%DerSrcP_SrdT(me%sizeNodes),me%DerSrcSdUi(Nb_VarC,me%sizeNodes))
        allocate(me%MatrixPS(2,me%sizeNodes),me%MatrixSP(Nb_VarC,me%sizeNodes),me%SurfCont(me%sizeNodes))
        me%MatrixPS=0.0_dp; me%MatrixSP=0.0_dp; me%SurfCont=0.0_dp; me%DerSrcP_SedT=0.0_dp; me%DerSrcP_SrdT=0.0_dp
        me%DerSrcSdUi=0.0_dp

        do i = 1, me%sizeNodes
            me%lk(1,i)%p => krn%FS_ports(offset1 + nodes1(i))
        enddo

        nodes2=links(2)%int1d('node')
        name2=links(2)%str('id')
        offset2 = krn%FS_list%find(name2, 1) - 1
    
        call check_if_conversion_required(krn%PortIdxForFSLink_List,name2,conversion_required,max_value)
        if(conversion_required) then
            allocate(conv(max_value))
            call array_conv(krn%PortIdxForFSLink_List,name2,conv)
            do i = 1, me%sizeNodes
                me%lk(2,i)%p => krn%FS_ports(offset2 + conv(nodes2(i)))
            enddo
        else
            do i = 1, me%sizeNodes
                me%lk(2,i)%p => krn%FS_ports(offset2 + nodes2(i))
            enddo
        endif
        
        if(me%lk(2,1)%p%typeS=="solid") then
            distances=links(2)%dbl1d('distance', n=me%sizeNodes)
            surfaces=cfg%dbl1d('contact_surface', n=me%sizeNodes)
            do i = 1, me%sizeNodes
              me%lk(2,i)%p%Dist4Grad=distances(i)
              me%SurfCont(i)=surfaces(i)
            enddo
        else if(me%lk(2,1)%p%typeS=="strand") then
            me%WetPer=cfg%dbl('wetted_perimeter')
        endif
        
        call nusselt_init(me%nuss,cfg) ! nusselt required in the link  
        if(.NOT.associated(me%nuss)) error stop 'Nusselt does not exist in thermalink'       
              
        do i = 1, me%sizeNodes
            call krn%coo_add_link(me%lk(2,i)%p, me%lk(1,i)%p, idx_1x4_col, idx_1x4_row, to_ptr_arr(me%MatrixSP(:,i)))
            call krn%coo_add_link(me%lk(1,i)%p, me%lk(2,i)%p, idx_2x1_col, idx_2x1_row, to_ptr_arr(me%MatrixPS(:,i)))
        enddo

    else if(links(2)%has_key('label')) then

        allocate(me%lk(maxval(me%nb_2D_Ports(:))+1,me%sizeNodes))
        allocate(me%DerSrcP_SedT2D(maxval(me%nb_2D_Ports(:)),me%sizeNodes))
        allocate(me%DerSrcP_SrdT2D(maxval(me%nb_2D_Ports(:)),me%sizeNodes))
        allocate(me%derFlx2DdUi(Nb_VarC,maxval(me%nb_2D_Ports(:)),me%sizeNodes))
        allocate(me%MatrixPS_2D(2,maxval(me%nb_2D_Ports(:)),me%sizeNodes))
        allocate(me%MatrixSP_2D(Nb_VarC,maxval(me%nb_2D_Ports(:)),me%sizeNodes))

        labels=links(2)%str1d('label')
        do i = 1, me%sizeNodes
            me%lk(1,i)%p => krn%FS_ports(offset1 + nodes1(i))
            do j=1,me%nb_2D_Ports(i)
              me%lk(1+j,i)%p => krn%FS_ports(krn%FS_list%find(links(2)%str('id')//"+"//trim(labels(i)%p),j))
            enddo
        enddo

        me%Ht=cfg%dbl('heat_coefficient')    

        do i = 1, me%sizeNodes
            do j=1,me%nb_2D_Ports(i)
                call krn%coo_add_link(me%lk(1+j,i)%p, me%lk(1,i)%p, idx_1x4_col, idx_1x4_row, to_ptr_arr(me%MatrixSP_2D(:,j,i)))
                call krn%coo_add_link(me%lk(1,i)%p, me%lk(1+j,i)%p, idx_2x1_col, idx_2x1_row, to_ptr_arr(me%MatrixPS_2D(:,j,i)))
            enddo
        enddo        

    endif

end subroutine FSlink_init_part2

end module cmp_FSlink_init_m