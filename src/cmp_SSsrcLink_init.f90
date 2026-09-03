! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module cmp_SSsrcLink_init_m
    use krn_interface_m
    use lib_input_m, only: input_t
    use lib_material_m
    implicit none

    type SSsrcLink_t   
        integer :: sizeNodes 
        type(SS_src_port_pointer_t), allocatable :: lk(:,:)
        class(material_t), pointer :: mat_ins

        integer, allocatable  :: nb_2D_Ports(:)

        real(dp), allocatable :: SurfCont(:), Lins(:), Lk12(:), Lk21(:),  DerSrcS1dT2(:), DerSrcS2dT1(:)
        real(dp), allocatable :: LkArr(:,:,:)
        real(dp), allocatable :: Der4LkArr(:,:,:)
    end type SSsrcLink_t

contains

subroutine SSsrcLink_init_part1(me,krn,cfg)
    type(SSsrcLink_t), intent(out) :: me
    type(krn_t), intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg

    integer :: nb_non_zeros_exp, nb_non_zeros_imp, nbPorts, i
    type(input_t), allocatable :: links(:)
    type(str_ptr), allocatable :: labels(:)

    links = cfg%dict1d('link')
    me%sizeNodes=size(links(1)%int1d('node'))

    if(links(2)%has_key('node')) then 
      nb_non_zeros_exp = 2*me%sizeNodes
      nb_non_zeros_imp = 2*me%sizeNodes
    else if(links(2)%has_key('label')) then 
      allocate(me%nb_2D_Ports(me%sizeNodes))
      labels=links(2)%str1d('label')
      do i = 1, me%sizeNodes
        me%nb_2D_Ports(i) = krn%SS_src_list%find_nbPorts(links(2)%str('id')//"+"//trim(labels(i)%p))
      enddo
      nb_non_zeros_exp = 2*sum(me%nb_2D_Ports(:))
      nb_non_zeros_imp = 2*sum(me%nb_2D_Ports(:))
    else
        error stop 'solidlink: link[2] requires "node" or "label"'
    endif

    call krn%add('',0,nb_non_zeros_exp,nb_non_zeros_imp,0,0,0,0,0)
end subroutine SSsrcLink_init_part1


subroutine SSsrcLink_init_part2(me,krn,cfg)
    type(SSsrcLink_t), intent(inout) :: me
    type(krn_t), target, intent(inout) :: krn
    class(input_t), pointer, intent(in)    :: cfg

    type(input_t), allocatable :: links(:)
    real(dp), allocatable :: distances1(:), distances2(:), surfaces(:), thicks(:)
    integer, allocatable :: nodes1(:), nodes2(:)
    type(str_ptr), allocatable :: labels(:)
    integer :: i,offset1,offset2,j

    ! TODO : Check that chunck is always index 1, mesh2D/solid index 2
    
    links = cfg%dict1d('link')

    nodes1=links(1)%int1d('node')
    offset1 = krn%SS_src_list%find(links(1)%str('id'), 1) - 1
    distances1=links(1)%dbl1d('distance', n=me%sizeNodes)

    call material_init(me%mat_ins,cfg)
    
    if(links(2)%has_key('node')) then 
      allocate(me%lk(2,me%sizeNodes),me%SurfCont(me%sizeNodes),me%Lins(me%sizeNodes))
      allocate(me%DerSrcS1dT2(me%sizeNodes),me%DerSrcS2dT1(me%sizeNodes))
      allocate(me%Lk12(me%sizeNodes),me%Lk21(me%sizeNodes))
      nodes2=links(2)%int1d('node')
      offset2 = krn%SS_src_list%find(links(2)%str('id'), 1) - 1
      distances2=links(2)%dbl1d('distance', n=me%sizeNodes)
      surfaces=cfg%dbl1d('contact_surface', n=me%sizeNodes)
      thicks=cfg%dbl1d('thickness', n=me%sizeNodes)
      do i = 1, me%sizeNodes
          me%lk(1,i)%p => krn%SS_src_ports(offset1 + nodes1(i))
          me%lk(1,i)%p%Dist4Grad = distances1(i)
          me%lk(2,i)%p => krn%SS_src_ports(offset2 + nodes2(i))    
          me%lk(2,i)%p%Dist4Grad = distances2(i)
          me%SurfCont(i)=surfaces(i)
          me%Lins(i)=thicks(i)
      enddo

      do i = 1, me%sizeNodes
        call krn%coo_add_link(me%lk(1,i)%p, me%lk(2,i)%p, [1], [1], to_ptr_arr(me%Lk12(i)))
        call krn%coo_add_link(me%lk(2,i)%p, me%lk(1,i)%p, [1], [1], to_ptr_arr(me%Lk21(i)))    
      enddo

    else if(links(2)%has_key('label')) then 
      ! Concept is one flag per chunk.
      ! In old version, one chunk could be associated to 2 flags (contacts _H and _V) --> requires here accumulation
      allocate(me%lk(maxval(me%nb_2D_Ports(:))+1,me%sizeNodes))
      allocate(me%Lins(me%sizeNodes),me%LkArr(maxval(me%nb_2D_Ports(:))+1,maxval(me%nb_2D_Ports(:))+1,me%sizeNodes))
      allocate(me%Der4LkArr(maxval(me%nb_2D_Ports(:))+1,maxval(me%nb_2D_Ports(:))+1,me%sizeNodes))
      labels=links(2)%str1d('label')
      thicks=cfg%dbl1d('thickness', n=me%sizeNodes)
      do i = 1, me%sizeNodes
        me%lk(1,i)%p => krn%SS_src_ports(offset1 + nodes1(i))
        me%lk(1,i)%p%Dist4Grad = distances1(i)
        do j=1,me%nb_2D_Ports(i)
          me%lk(1+j,i)%p => krn%SS_src_ports(krn%SS_src_list%find(links(2)%str('id')//"+"//trim(labels(i)%p),j))
        enddo
        me%Lins(i)=thicks(i)
      enddo

      do i = 1, me%sizeNodes      
        do j=1,me%nb_2D_Ports(i)        
          call krn%coo_add_link(me%lk(1,i)%p, me%lk(1+j,i)%p, [1], [1], to_ptr_arr(me%LkArr(1,1+j,i)))
          call krn%coo_add_link(me%lk(1+j,i)%p, me%lk(1,i)%p, [1], [1], to_ptr_arr(me%LkArr(1+j,1,i)))    
        enddo
      enddo

    endif

end subroutine SSsrcLink_init_part2

end module cmp_SSsrcLink_init_m