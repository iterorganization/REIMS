! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_mesh2D_init_m
    use krn_interface_m
    use krn_simulation_m, only: signal_t, simulation_t
    use lib_input_m, only: input_t
    use hdf5    ! TODO : To be removed
    use lib_material_m
    implicit none

    type coord_t
      real(kind=dp) :: x_coord,y_coord
    end type coord_t

    type elements_1D_t
      integer :: nod(2)                ! Node indices defining a given 1D element
      integer :: PhysIdx                          ! Physical label for each edge
    end type elements_1D_t

    type elements_t
      integer :: nod(3)                !! Node indices defining a given element
      integer :: Neigh(3)              !! Neighbor indices of a given element (3rd element = zero if only 2 neighbors)
      real(kind=dp) :: face(3)         !! Lengths of the 3 faces of the element (consistent with the sorting of neighbors)
      real(kind=dp) :: Delta(3)        !! Distance normal to the face from the cell center (consistent with the sorting of neighbors)
      real(kind=dp) :: DeltaNeigh(3)   !! Distance normal to the face from the cell center of the associated neighbor (consistent with the sorting of neighbors)
      type(coord_t) :: norm(3)         !! Unit outgoing normales at the level of the 3 faces (consistent with the sorting of neighbors)
      type(coord_t) :: NeighCenter(3)  !! Coordinates of the center of each neighbor cell (consistent with the sorting of neighbors)
      type(coord_t) :: centTOface(3)   !! Vector from the center of the element to the face center (consistent with the sorting of neighbors)
      type(coord_t) :: centTOfaceN(3)  !! Vector from the center of the neighbor to the face center (consistent with the sorting of neighbors)
      character(len=50) :: PhysN(3)    !! Physical group associated to each of the 2 or 3 neighbors (domain) OR 1 boundary (consistent with the sorting of neighbors)
      character(len=50) :: PhysE                  !! Physical group associated to the element (as a domain)
      integer :: PhysIdxE                         !! Physical group index associated to the element (as a domain)
      type(coord_t) :: center                     !! Coordinates of the center of the cell
      integer :: NbNeigh                          !! Number of neighbors for a given element
      real(kind=dp) :: surface                    !! Surface of a given element
    end type elements_t    

    type mesh2D_prop_t
        integer :: nb_nodes                          !! Number of nodes
        type(coord_t), allocatable :: node(:)        !! Coordinates of nodes
        integer :: nb_elements                       !! Number of 2D elements
        type(elements_t), allocatable :: elem(:)     !! 2D element of mesh
        integer :: nnz                               !! Number of non zeros useful at different locations of init file
        integer  :: Nb_SSports, Nb_FSports, Nb_label, Nb_labelMC, Nb_labelCh
        character(len=50), allocatable :: LabelName(:)
        type(FS_port_pointer_t), allocatable :: thermP(:)
        type(SS_src_port_pointer_t), allocatable :: thermS(:)
    end type mesh2D_prop_t 

    type StateVariable_mesh2D_t
        real(dp), allocatable :: temp(:) !! Temperature of composite mesh2D
    end type StateVariable_mesh2D_t    

    type big_arrays_mesh2D_t
        real(dp), allocatable :: bmx(:) !! main diagonal subMatrix    
        real(dp), allocatable :: dmx(:,:) !! regarding connection with direct 2D neighbors

        real(dp), pointer :: bvx(:)    
    end type big_arrays_mesh2D_t     

    type arrays_linear_system_mesh2D_t
        type(dbl_pointer_t), allocatable :: ValExp(:), ValImp(:) !! Non zero values in COO format for 1 mesh2D
        integer, allocatable :: RowExp(:), RowImp(:)  !! row index array in COO format for 1 mesh2D
        integer, allocatable :: ColExp(:), ColImp(:)  !! col index array in COO format for 1 mesh2D
    end type arrays_linear_system_mesh2D_t      

    type flux_mesh2D_t
        real(dp), allocatable :: SomNormQT(:), DerSomNormQT(:)
        real(dp), allocatable :: DerNormQTOth2D(:,:)
    end type flux_mesh2D_t

    type lab_prop_t
        character(:), allocatable :: id, type
        type(signal_t) :: value
        real(kind=dp) :: SurfaceRegion, LengthBC
        class(material_t), pointer :: mat
    end type lab_prop_t

    type labelLk_t
      integer :: NbElem
      integer :: array(100)
      character(:), allocatable :: name
    end type labelLk_t    

    type label_t
        type(lab_prop_t), allocatable :: regionArr(:), edgeArr(:)
        type(labelLk_t), allocatable  :: MC_Arr(:), Ch_Arr(:)
    end type label_t

    type mesh2D_t    
        type(StateVariable_mesh2D_t)        :: StVar            !! StateVariable updated by BDF 2 scheme -> time n+1
        type(StateVariable_mesh2D_t)        :: StVarOld         !! StateVariable -> time n
        type(StateVariable_mesh2D_t)        :: StVarOld2        !! StateVariable -> time n-1
        type(StateVariable_mesh2D_t)        :: StVarOld3        !! StateVariable -> time n-2
        type(StateVariable_mesh2D_t)        :: StVarOld4        !! StateVariable -> time n-3
        type(StateVariable_mesh2D_t)        :: StVarOld5        !! StateVariable -> time n-4
        type(StateVariable_mesh2D_t)        :: StVarOld6        !! StateVariable -> time n-5
        type(mesh2D_prop_t)                 :: M2D_Prop         !! Variables relative to mesh2D global properties
        ! type(scenario_mesh2D_t)           :: scen             !! Electro-magnetic + Heat scenario
        type(big_arrays_mesh2D_t)           :: big              !! Big arrays
        type(arrays_linear_system_mesh2D_t) :: bigLS            !! Big arrays for Linear system
        type(flux_mesh2D_t)                 :: flxS             !! Heat diffusion fluxes
        character(:), allocatable           :: filename_extless !! filename without extension
        type(label_t)                       :: label
        real(dp)                            :: err
        real(dp)                            :: err_den   
        integer(hid_t) :: f_id_writing   ! TODO : to be removed
        real(dp)                            :: extrusion_length  !! extrusion length of the 2D slice (m)
    end type mesh2D_t    

contains



subroutine preprocessing_2Delements(me,mshFile)  
    type(mesh2D_prop_t), intent(out) :: me
    character(*), intent(in) :: mshFile

    type(elements_1D_t), allocatable :: elem1D(:) ! Type for 1D element (only nodes are required)
    
    real(dp) :: den, x_face, y_face, check
    real(dp) :: nn(2),cf(2)
 
    integer, allocatable :: cnt(:,:), nod1(:,:)
    integer :: nod_shar(2)
    integer :: cnt_nod_shar,cnt_nod_shar_1D,nod_alone,ii
    integer :: Nb_1D_elements
    integer :: line,nb_n,nb_e,nb_p,nb_1De
    integer :: j_e,nb_v,j_v,j_1De
    integer :: pt1,pt2,pt3,pt4,pt5

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(unit=1,file=mshFile) ! msh file generated by GMSH     -->    export - File.msh - Version 2 ASCII  (without any additional option)
    do line=1,3 ! Mesh Format
      read(1,*)
    enddo
    read(1,*)  ! Physical names
    read(1,*) me%Nb_label
    do line=1,me%Nb_label
      read(1,*)
    enddo
    read(1,*)
    read(1,*)  ! nodes
    read(1,*) me%nb_nodes
    do line=1,me%nb_nodes
      read(1,*)
    enddo
    read(1,*)
    read(1,*)
    read(1,*) me%nb_elements  ! Elements (multi-D at this stage)
    do line=1,me%nb_elements
      read(1,*) pt1, pt2
      if(pt2==2) then
        Nb_1D_elements=pt1-1
        exit
      endif
    enddo
    me%nb_elements=me%nb_elements-Nb_1D_elements ! to keep only the number of 2D elements
    close(1)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(me%LabelName(me%Nb_label))
    allocate(me%node(me%nb_nodes),me%elem(me%nb_elements),elem1D(Nb_1D_elements))
    allocate(cnt(me%nb_elements,me%nb_elements),nod1(me%nb_elements,me%nb_elements))
  
    do nb_n=1,me%nb_nodes
      me%node(nb_n)%x_coord=0.0_dp;me%node(nb_n)%y_coord=0.0_dp
    enddo
    do nb_e=1,me%nb_elements
      me%elem(nb_e)%nod(:)=0
    enddo
    do nb_e=1,Nb_1D_elements
      elem1D(nb_e)%nod(:)=0
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(unit=1,file=mshFile)
    do line=1,3 ! Mesh Format
      read(1,*)
    enddo
    read(1,*)  ! Physical names
    read(1,*) 
    do line=1,me%Nb_label
      read(1,*) pt1, pt2, me%LabelName(line)  ! In practice, pt2=line
    enddo
    read(1,*)
    read(1,*)  ! nodes
    read(1,*)
    do line=1,me%nb_nodes
      read(1,*) pt1, me%node(line)%x_coord, me%node(line)%y_coord, check
    enddo
    read(1,*)
    read(1,*)  ! Elements
    read(1,*)
    do line=1,Nb_1D_elements
      read(1,*) pt1, pt2, pt3, elem1D(line)%PhysIdx, pt4, elem1D(line)%nod(:)
    enddo
    do line=1,me%nb_elements
      read(1,*) pt1, pt2, pt3, me%elem(line)%PhysIdxE, pt5, me%elem(line)%nod(:)
    enddo
    close(1) ! End of reading
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    cnt(:,:)=0
    me%nnz=0
    do nb_e=1,me%nb_elements

      me%elem(nb_e)%PhysE=trim(me%LabelName(me%elem(nb_e)%PhysIdxE))
  
      me%elem(nb_e)%surface=0.5_dp*&
                  abs(me%node(me%elem(nb_e)%nod(1))%x_coord*(me%node(me%elem(nb_e)%nod(2))%y_coord-&
                      me%node(me%elem(nb_e)%nod(3))%y_coord)+&
                      me%node(me%elem(nb_e)%nod(2))%x_coord*(me%node(me%elem(nb_e)%nod(3))%y_coord-&
                      me%node(me%elem(nb_e)%nod(1))%y_coord)+&
                      me%node(me%elem(nb_e)%nod(3))%x_coord*(me%node(me%elem(nb_e)%nod(1))%y_coord-&
                      me%node(me%elem(nb_e)%nod(2))%y_coord))
  
      me%elem(nb_e)%center%x_coord=(me%node(me%elem(nb_e)%nod(1))%x_coord+me%node(me%&
                                           elem(nb_e)%nod(2))%x_coord+&
                                           me%node(me%elem(nb_e)%nod(3))%x_coord)/3.0_dp
      me%elem(nb_e)%center%y_coord=(me%node(me%elem(nb_e)%nod(1))%y_coord+me%node(me%&
                                           elem(nb_e)%nod(2))%y_coord+&
                                           me%node(me%elem(nb_e)%nod(3))%y_coord)/3.0_dp
  
      me%elem(nb_e)%NbNeigh=0
      do nb_v=1,me%nb_elements
        do j_e=1,3
          do j_v=1,3    
            if((me%elem(nb_e)%nod(j_e)==me%elem(nb_v)%nod(j_v)) .and. (nb_v/=nb_e)) then
              cnt(nb_e,nb_v)=cnt(nb_e,nb_v)+1
              if(cnt(nb_e,nb_v)/=2) nod1(nb_e,nb_v)=me%elem(nb_e)%nod(j_e)
            endif
          enddo
  
          if(cnt(nb_e,nb_v)==2) then
            me%elem(nb_e)%NbNeigh=me%elem(nb_e)%NbNeigh+1            
            me%elem(nb_e)%Neigh(me%elem(nb_e)%NbNeigh)=nb_v
  
            me%elem(nb_e)%NeighCenter(me%elem(nb_e)%NbNeigh)%x_coord=(me%node(me%elem(nb_v)%&
                                                                      nod(1))%x_coord+&
                                                                      me%node(me%elem(nb_v)%nod(2))%x_coord+&
                                                                      me%node(me%elem(nb_v)%nod(3))%x_coord)/3.0_dp
            me%elem(nb_e)%NeighCenter(me%elem(nb_e)%NbNeigh)%y_coord=(me%node(me%elem(nb_v)%&
                                                                      nod(1))%y_coord+&
                                                                      me%node(me%elem(nb_v)%nod(2))%y_coord+&
                                                                      me%node(me%elem(nb_v)%nod(3))%y_coord)/3.0_dp
  
            me%elem(nb_e)%face(me%elem(nb_e)%NbNeigh)=sqrt((me%node(me%elem(nb_e)%nod(j_e))%&
                                               x_coord-me%node(nod1(nb_e,nb_v))%x_coord)**2+&
                                               (me%node(me%elem(nb_e)%nod(j_e))%y_coord-&
                                               me%node(nod1(nb_e,nb_v))%y_coord)**2)
  

            x_face=0.5_dp*(me%node(me%elem(nb_e)%nod(j_e))%x_coord+me%node(nod1(nb_e,nb_v))%x_coord)
            y_face=0.5_dp*(me%node(me%elem(nb_e)%nod(j_e))%y_coord+me%node(nod1(nb_e,nb_v))%y_coord)

            me%elem(nb_e)%centTOface(me%elem(nb_e)%NbNeigh)%x_coord=x_face-me%elem(nb_e)%center%x_coord
            me%elem(nb_e)%centTOface(me%elem(nb_e)%NbNeigh)%y_coord=y_face-me%elem(nb_e)%center%y_coord

            me%elem(nb_e)%centTOfaceN(me%elem(nb_e)%NbNeigh)%x_coord=x_face-&
                                                         me%elem(nb_e)%NeighCenter(me%elem(nb_e)%NbNeigh)%x_coord
            me%elem(nb_e)%centTOfaceN(me%elem(nb_e)%NbNeigh)%y_coord=y_face-&
                                                         me%elem(nb_e)%NeighCenter(me%elem(nb_e)%NbNeigh)%y_coord
  

            den=me%elem(nb_e)%face(me%elem(nb_e)%NbNeigh)
            me%elem(nb_e)%norm(me%elem(nb_e)%NbNeigh)%x_coord=-(me%node(me%&
                                                                elem(nb_e)%nod(j_e))%y_coord-&
                                                                me%node(nod1(nb_e,nb_v))%y_coord)/den
            me%elem(nb_e)%norm(me%elem(nb_e)%NbNeigh)%y_coord= (me%node(me%&
                                                                elem(nb_e)%nod(j_e))%x_coord-&
                                                                me%node(nod1(nb_e,nb_v))%x_coord)/den
  
            ! check for unit normales to be outgoing
            nn(1)=me%elem(nb_e)%norm(me%elem(nb_e)%NbNeigh)%x_coord
            nn(2)=me%elem(nb_e)%norm(me%elem(nb_e)%NbNeigh)%y_coord
            cf(1)=me%elem(nb_e)%centTOface(me%elem(nb_e)%NbNeigh)%x_coord
            cf(2)=me%elem(nb_e)%centTOface(me%elem(nb_e)%NbNeigh)%y_coord
            if(dot_product(nn,cf)<=0.0_dp) then
              me%elem(nb_e)%norm(me%elem(nb_e)%NbNeigh)%x_coord=&
                                               -me%elem(nb_e)%norm(me%elem(nb_e)%NbNeigh)%x_coord
              me%elem(nb_e)%norm(me%elem(nb_e)%NbNeigh)%y_coord=&
                                               -me%elem(nb_e)%norm(me%elem(nb_e)%NbNeigh)%y_coord
            endif
            nn(1)=me%elem(nb_e)%norm(me%elem(nb_e)%NbNeigh)%x_coord
            nn(2)=me%elem(nb_e)%norm(me%elem(nb_e)%NbNeigh)%y_coord
            me%elem(nb_e)%Delta(me%elem(nb_e)%NbNeigh)=abs(dot_product(nn,cf))

            cf(1)=me%elem(nb_e)%centTOfaceN(me%elem(nb_e)%NbNeigh)%x_coord
            cf(2)=me%elem(nb_e)%centTOfaceN(me%elem(nb_e)%NbNeigh)%y_coord
            me%elem(nb_e)%DeltaNeigh(me%elem(nb_e)%NbNeigh)=abs(dot_product(-nn,cf))
  
            cnt(nb_e,nb_v)=-10 ! to avoid duplication
  
          endif
        enddo
      enddo
  
      if(me%elem(nb_e)%NbNeigh==1) then

        me%elem(nb_e)%Neigh(2)=0
        me%elem(nb_e)%Neigh(3)=0

        nod_shar(:)=0
        do j_e=1,3
          cnt_nod_shar=0
          ! Neigh 1 (only one neighbour)
          do j_v=1,3
            if(me%elem(nb_e)%nod(j_e)==me%elem(me%elem(nb_e)%Neigh(1))%nod(j_v)) then
              cnt_nod_shar=cnt_nod_shar+1
            endif
          enddo
  
          if(cnt_nod_shar==1) then
            if(nod_shar(1)==0) then
              nod_shar(1)=me%elem(nb_e)%nod(j_e)
            else
              nod_shar(2)=me%elem(nb_e)%nod(j_e)
            endif
          else if(cnt_nod_shar==0) then
              nod_alone=me%elem(nb_e)%nod(j_e)
          endif
        enddo

        do ii=2,3
          me%elem(nb_e)%face(ii)=sqrt((me%node(nod_shar(ii-1))%x_coord-me%node(nod_alone)%x_coord)**2+&
                                    (me%node(nod_shar(ii-1))%y_coord-&
                                    me%node(nod_alone)%y_coord)**2)
    
          x_face=0.5_dp*(me%node(nod_shar(ii-1))%x_coord+me%node(nod_alone)%x_coord)
          y_face=0.5_dp*(me%node(nod_shar(ii-1))%y_coord+me%node(nod_alone)%y_coord)
          me%elem(nb_e)%centTOface(ii)%x_coord=x_face-me%elem(nb_e)%center%x_coord
          me%elem(nb_e)%centTOface(ii)%y_coord=y_face-me%elem(nb_e)%center%y_coord
    
          den=me%elem(nb_e)%face(ii)
          me%elem(nb_e)%norm(ii)%x_coord=-(me%node(nod_shar(ii-1))%y_coord-me%node(nod_alone)%y_coord)/den
          me%elem(nb_e)%norm(ii)%y_coord= (me%node(nod_shar(ii-1))%x_coord-me%node(nod_alone)%x_coord)/den
    
          ! check for unit normales to be outgoing
          nn(1)=me%elem(nb_e)%norm(ii)%x_coord
          nn(2)=me%elem(nb_e)%norm(ii)%y_coord
          cf(1)=me%elem(nb_e)%centTOface(ii)%x_coord
          cf(2)=me%elem(nb_e)%centTOface(ii)%y_coord
          if(dot_product(nn,cf)<=0.0_dp) then
            me%elem(nb_e)%norm(ii)%x_coord=-me%elem(nb_e)%norm(ii)%x_coord
            me%elem(nb_e)%norm(ii)%y_coord=-me%elem(nb_e)%norm(ii)%y_coord
          endif
          nn(1)=me%elem(nb_e)%norm(ii)%x_coord
          nn(2)=me%elem(nb_e)%norm(ii)%y_coord
          me%elem(nb_e)%Delta(ii)=abs(dot_product(nn,cf))
        enddo

        me%elem(nb_e)%PhysN(1)=trim(me%LabelName(me%elem(me%elem(nb_e)%Neigh(1))%PhysIdxE))

        do ii=2,3
          do nb_1De=1,Nb_1D_elements
            cnt_nod_shar_1D=0
            do j_1De=1,2
              if(nod_shar(ii-1)==elem1D(nb_1De)%nod(j_1De) .or. nod_alone==elem1D(nb_1De)%nod(j_1De)) then
                cnt_nod_shar_1D=cnt_nod_shar_1D+1
              endif
            enddo        
            if(cnt_nod_shar_1D==2) then
              me%elem(nb_e)%PhysN(ii)=trim(me%LabelName(elem1D(nb_1De)%PhysIdx))
            endif
          enddo
        enddo
  
      else if(me%elem(nb_e)%NbNeigh==3) then
        do j_v=1,3
          me%elem(nb_e)%PhysN(j_v)=trim(me%LabelName(me%elem(me%elem(nb_e)%Neigh(j_v))%PhysIdxE))
        enddo

      else if(me%elem(nb_e)%NbNeigh==2) then
        me%elem(nb_e)%Neigh(3)=0
  
        nod_shar(:)=0
        do j_e=1,3
          cnt_nod_shar=0
          ! Neigh 1
          do j_v=1,3
            if(me%elem(nb_e)%nod(j_e)==me%elem(me%elem(nb_e)%Neigh(1))%nod(j_v)) then
              cnt_nod_shar=cnt_nod_shar+1
            endif
          enddo
  
          ! Neigh 2
          do j_v=1,3
            if(me%elem(nb_e)%nod(j_e)==me%elem(me%elem(nb_e)%Neigh(2))%nod(j_v)) then
              cnt_nod_shar=cnt_nod_shar+1
            endif
          enddo
  
          if(cnt_nod_shar==1) then
            if(nod_shar(1)==0) then
              nod_shar(1)=me%elem(nb_e)%nod(j_e)
            else
              nod_shar(2)=me%elem(nb_e)%nod(j_e)
            endif
          endif
        enddo
  
        me%elem(nb_e)%face(3)=sqrt((me%node(nod_shar(2))%x_coord-me%node(nod_shar(1))%x_coord)**2+&
                                  (me%node(nod_shar(2))%y_coord-&
                                  me%node(nod_shar(1))%y_coord)**2)
  
        x_face=0.5_dp*(me%node(nod_shar(2))%x_coord+me%node(nod_shar(1))%x_coord)
        y_face=0.5_dp*(me%node(nod_shar(2))%y_coord+me%node(nod_shar(1))%y_coord)
        me%elem(nb_e)%centTOface(3)%x_coord=x_face-me%elem(nb_e)%center%x_coord
        me%elem(nb_e)%centTOface(3)%y_coord=y_face-me%elem(nb_e)%center%y_coord
  
        den=me%elem(nb_e)%face(3)
        me%elem(nb_e)%norm(3)%x_coord=-(me%node(nod_shar(2))%y_coord-me%node(nod_shar(1))%y_coord)/den
        me%elem(nb_e)%norm(3)%y_coord= (me%node(nod_shar(2))%x_coord-me%node(nod_shar(1))%x_coord)/den
  
        ! check for unit normales to be outgoing
        nn(1)=me%elem(nb_e)%norm(3)%x_coord
        nn(2)=me%elem(nb_e)%norm(3)%y_coord
        cf(1)=me%elem(nb_e)%centTOface(3)%x_coord
        cf(2)=me%elem(nb_e)%centTOface(3)%y_coord
        if(dot_product(nn,cf)<=0.0_dp) then
          me%elem(nb_e)%norm(3)%x_coord=-me%elem(nb_e)%norm(3)%x_coord
          me%elem(nb_e)%norm(3)%y_coord=-me%elem(nb_e)%norm(3)%y_coord
        endif
        nn(1)=me%elem(nb_e)%norm(3)%x_coord
        nn(2)=me%elem(nb_e)%norm(3)%y_coord
        me%elem(nb_e)%Delta(3)=abs(dot_product(nn,cf))
  
        do j_v=1,2
          me%elem(nb_e)%PhysN(j_v)=trim(me%LabelName(me%elem(me%elem(nb_e)%Neigh(j_v))%PhysIdxE))
        enddo

        do nb_1De=1,Nb_1D_elements
          cnt_nod_shar_1D=0
          do j_1De=1,2
            if(nod_shar(1)==elem1D(nb_1De)%nod(j_1De) .or. nod_shar(2)==elem1D(nb_1De)%nod(j_1De)) then
              cnt_nod_shar_1D=cnt_nod_shar_1D+1
            endif
          enddo        
          if(cnt_nod_shar_1D==2) then
            me%elem(nb_e)%PhysN(3)=trim(me%LabelName(elem1D(nb_1De)%PhysIdx))
          endif
        enddo
  
      endif

      me%nnz=me%nnz+1+me%elem(nb_e)%NbNeigh

      if(me%elem(nb_e)%NbNeigh==1) then
        print*,"ACCUMULATION TO BE CONSIDERED FOR 2D CELLS THAT ARE CONNECTED TO 2 DIFFERENT LABELS ?"
        print*," check if it is the same for the 2D cells connected to labels"
        read(*,*)
      endif
 
    enddo

    print*,mshFile,' --->  Preprocessing performed'

end subroutine preprocessing_2Delements



subroutine all_mesh2D_allocation(me)
    type(mesh2D_t), intent(inout) :: me
    integer :: NbCells
    
    NbCells=me%M2D_Prop%nb_elements
    allocate(me%StVar%temp(NbCells))
    allocate(me%StVarOld%temp(NbCells))
    allocate(me%StVarOld2%temp(NbCells))
    allocate(me%StVarOld3%temp(NbCells))
    allocate(me%StVarOld4%temp(NbCells))
    allocate(me%StVarOld5%temp(NbCells))
    allocate(me%StVarOld6%temp(NbCells))
    allocate(me%flxS%SomNormQT(NbCells))
    allocate(me%flxS%DerSomNormQT(NbCells))
    allocate(me%flxS%DerNormQTOth2D(NbCells,3))
    ! allocate(me%scen%Q_ext(NbCells))
    allocate(me%big%bmx(NbCells))
    allocate(me%big%dmx(NbCells,3))
      
    me%StVar%temp=0.0_dp; me%StVarOld%temp=0.0_dp
    me%StVarOld2%temp=0.0_dp; me%StVarOld3%temp=0.0_dp
    me%StVarOld4%temp=0.0_dp; me%StVarOld5%temp=0.0_dp; me%StVarOld6%temp=0.0_dp
    me%flxS%SomNormQT=0.0_dp; me%flxS%DerSomNormQT=0.0_dp; me%flxS%DerNormQTOth2D=0.0_dp
    ! me%scen%Q_ext=0.0_dp
    me%big%bmx=0.0_dp;me%big%dmx=0.0_dp
end subroutine all_mesh2D_allocation



subroutine mesh2D_init_part1(me, krn, cfg,sim)
    type(mesh2D_t),         intent(out) :: me
    type(krn_t),          intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg
    type(simulation_t), intent(inout) :: sim

    type(input_t), allocatable :: constit(:)
    real(dp) :: T_init
    integer :: nb_e, nb_SS_src_ports, ii, jj, i, j, count
    type(input_t), allocatable :: regions(:), edges(:)
    type(str_ptr), allocatable :: labelsMC(:), labelsCh(:)

    call preprocessing_2Delements(me%M2D_Prop,cfg%str('mesh'))
    call all_mesh2D_allocation(me)     

    me%filename_extless = cfg%str('mesh')
    me%filename_extless = me%filename_extless(9:index(me%filename_extless, '.msh')-1)

    me%extrusion_length = cfg%dbl('extrusion_length',1.062_dp)

    T_init=cfg%dbl('initial/t')
    do nb_e=1,me%M2D_Prop%nb_elements
        me%StVar%temp(nb_e)=T_init
    enddo

    regions = cfg%dict1d('regions')
    allocate(me%label%regionArr(size(regions)))
    do ii=1, size(regions)
        call material_init(me%label%regionArr(ii)%mat, regions(ii))
        me%label%regionArr(ii)%id=regions(ii)%str('label')
        call me%label%regionArr(ii)%value%init(sim, regions(ii),'source')  ! value zero by default if it is not found
        me%label%regionArr(ii)%SurfaceRegion=0.0_dp
        do nb_e=1,me%M2D_Prop%nb_elements
            if(trim(me%M2D_Prop%elem(nb_e)%PhysE)==me%label%regionArr(ii)%id) then
                me%label%regionArr(ii)%SurfaceRegion=me%label%regionArr(ii)%SurfaceRegion+me%M2D_Prop%elem(nb_e)%surface
            endif
        enddo
        if (me%label%regionArr(ii)%SurfaceRegion == 0.0_dp) then
            print*, 'ERROR: regions label "'// &
                trim(me%label%regionArr(ii)%id)//'" not found in mesh.'
            print*, 'Available labels in mesh:'
            do j = 1, me%M2D_Prop%Nb_label
                print*, '  - '//trim(me%M2D_Prop%LabelName(j))
            enddo
            error stop 'Aborting: regions/mesh mismatch.'
        endif
    enddo

    edges = cfg%dict1d('edges')
    allocate(me%label%edgeArr(size(edges)))
    do ii=1, size(edges)
        me%label%edgeArr(ii)%id=edges(ii)%str('label')
        me%label%edgeArr(ii)%type=edges(ii)%str('type')
        call me%label%edgeArr(ii)%value%init(sim, edges(ii),'value')
        me%label%edgeArr(ii)%LengthBC=0.0_dp
        do nb_e=1,me%M2D_Prop%nb_elements
          do jj=1,3
            if((jj==3 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==2) .or. &
               (jj==3 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==1) .or. &
               (jj==2 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==1)) then       ! Boundary condition
                 if(trim(me%M2D_Prop%elem(nb_e)%PhysN(jj))==me%label%edgeArr(ii)%id) then
                   me%label%edgeArr(ii)%LengthBC = me%label%edgeArr(ii)%LengthBC + me%M2D_Prop%elem(nb_e)%face(jj)
                 endif
            endif
          enddo
        enddo
        if (me%label%edgeArr(ii)%LengthBC == 0.0_dp) then
            print*, 'ERROR: edges label "'// &
                trim(me%label%edgeArr(ii)%id)//'" not found in mesh.'
            print*, 'Available labels in mesh:'
            do j = 1, me%M2D_Prop%Nb_label
                print*, '  - '//trim(me%M2D_Prop%LabelName(j))
            enddo
            error stop 'Aborting: edges/mesh mismatch.'
        endif
    enddo

    if(cfg%has_key('channel_link')) then
        labelsCh=cfg%str1d('channel_link')
        me%M2D_Prop%Nb_labelCh=size(labelsCh)
        if(me%M2D_Prop%Nb_labelCh>0) then
            allocate(me%label%Ch_Arr(me%M2D_Prop%Nb_labelCh))
            do i= 1, me%M2D_Prop%Nb_labelCh
                me%label%Ch_Arr(i)%name = labelsCh(i)%p
            end do    
            me%M2D_Prop%Nb_FSports = 0
            me%label%Ch_Arr(:)%NbElem=0
            do nb_e=1,me%M2D_Prop%nb_elements
                do jj=1,3
                    if((jj==3 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==2) .or. &
                       (jj==3 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==1) .or. &
                       (jj==2 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==1)) then       ! Boundary condition
                        do ii=1,me%M2D_Prop%Nb_labelCh
                            if((trim(me%M2D_Prop%elem(nb_e)%PhysN(jj)) == me%label%Ch_Arr(ii)%name) ) then
                                me%M2D_Prop%Nb_FSports = me%M2D_Prop%Nb_FSports + 1
                                me%label%Ch_Arr(ii)%NbElem = me%label%Ch_Arr(ii)%NbElem + 1
                                me%label%Ch_Arr(ii)%array(me%label%Ch_Arr(ii)%NbElem) = nb_e
                            endif
                        enddo
                    endif
                enddo
            enddo
            do i = 1, me%M2D_Prop%Nb_labelCh
                if (me%label%Ch_Arr(i)%NbElem == 0) then
                    print*, 'ERROR: channel_link label "'// &
                        trim(me%label%Ch_Arr(i)%name)//'" not found in mesh.'
                    print*, 'Available labels in mesh:'
                    do ii = 1, me%M2D_Prop%Nb_label
                        print*, '  - '//trim(me%M2D_Prop%LabelName(ii))
                    enddo
                    error stop 'Aborting: channel_link/mesh mismatch.'
                endif
            enddo
        endif
    endif

    if(cfg%has_key('solid_link')) then
        labelsMC=cfg%str1d('solid_link')
        me%M2D_Prop%Nb_labelMC=size(labelsMC)
        if(me%M2D_Prop%Nb_labelMC>0) then
            allocate(me%label%MC_Arr(me%M2D_Prop%Nb_labelMC))
            do i= 1, me%M2D_Prop%Nb_labelMC
                me%label%MC_Arr(i)%name = labelsMC(i)%p
            end do
            me%M2D_Prop%Nb_SSports = 0
            me%label%MC_Arr(:)%NbElem=0
            do nb_e=1,me%M2D_Prop%nb_elements
                do jj=1,3
                    if((jj==3 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==2) .or. &
                       (jj==3 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==1) .or. &
                       (jj==2 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==1)) then       ! Boundary condition
                        do ii=1,me%M2D_Prop%Nb_labelMC
                            if((trim(me%M2D_Prop%elem(nb_e)%PhysN(jj)) == me%label%MC_Arr(ii)%name) ) then
                                me%M2D_Prop%Nb_SSports = me%M2D_Prop%Nb_SSports + 1
                                me%label%MC_Arr(ii)%NbElem = me%label%MC_Arr(ii)%NbElem + 1
                                me%label%MC_Arr(ii)%array(me%label%MC_Arr(ii)%NbElem) = nb_e
                            endif
                        enddo
                    endif
                enddo
            enddo
            do i = 1, me%M2D_Prop%Nb_labelMC
                if (me%label%MC_Arr(i)%NbElem == 0) then
                    print*, 'ERROR: solid_link label "'// &
                        trim(me%label%MC_Arr(i)%name)//'" not found in mesh.'
                    print*, 'Available labels in mesh:'
                    do ii = 1, me%M2D_Prop%Nb_label
                        print*, '  - '//trim(me%M2D_Prop%LabelName(ii))
                    enddo
                    error stop 'Aborting: solid_link/mesh mismatch.'
                endif
            enddo
        endif
    endif

    call krn%add(cfg%str('id'),me%M2D_Prop%nb_elements,me%M2D_Prop%nnz,me%M2D_Prop%nnz,0,0,0,0,0)

    if(me%M2D_Prop%Nb_labelCh>0) then
        if(me%M2D_Prop%Nb_labelCh>99) print*,"Nb_labelCh>99"
        do ii=1,me%M2D_Prop%Nb_labelCh
            call krn%add_label(cfg%str('id')//"+"//trim(me%label%Ch_Arr(ii)%name), nb_nodes_FS=me%label%Ch_Arr(ii)%NbElem)
        enddo
    endif

    if(me%M2D_Prop%Nb_labelMC>0) then
        if(me%M2D_Prop%Nb_labelMC>99) print*,"Nb_labelMC>99"
        do ii=1,me%M2D_Prop%Nb_labelMC
            call krn%add_label(cfg%str('id')//"+"//trim(me%label%MC_Arr(ii)%name), nb_nodes_SS=me%label%MC_Arr(ii)%NbElem)
        enddo
    endif    
    
end subroutine mesh2D_init_part1


subroutine mesh2D_init_part2(me, krn, cfg)
    !! Fluid port initialisation for mesh2Ds
    type(mesh2D_t), target, intent(inout) :: me
    type(krn_t),          intent(inout) :: krn
    class(input_t), pointer, intent(in) :: cfg  

    integer :: nb_e,idx,PrevRowIdx,PrevColIdx,jj,i,ii,idx_prev,elem
    real(dp), pointer :: rhs(:)
    integer, allocatable :: list_thp_loc(:),list_thS_loc(:)
    type(FS_port_t), pointer :: FS_ports(:)
    type(SS_src_port_t), pointer :: SS_src_ports(:)

    call krn%update(rhs_or_solution_view=rhs)
    me%big%bvx(1:me%M2D_Prop%nb_elements) => rhs

    if(me%M2D_Prop%Nb_labelCh>0) then
      allocate(me%M2D_Prop%thermP(me%M2D_Prop%Nb_FSports),list_thp_loc(me%M2D_Prop%Nb_FSports))
      idx_prev=0
      do ii=1,me%M2D_Prop%Nb_labelCh
        list_thp_loc(idx_prev+1:idx_prev+me%label%Ch_Arr(ii)%NbElem)= &
                  me%label%Ch_Arr(ii)%array(1:me%label%Ch_Arr(ii)%NbElem)
        idx_prev=idx_prev+me%label%Ch_Arr(ii)%NbElem
      enddo
      call krn%update(FS_p_loc=list_thp_loc,FS_p_view=FS_ports)
      do i=1,me%M2D_Prop%Nb_FSports                          
        me%M2D_Prop%thermP(i)%p => FS_ports(i)
        me%M2D_Prop%thermP(i)%p%typeS = 'mesh2D'
      enddo      
    endif    

    if(me%M2D_Prop%Nb_labelMC>0) then
      allocate(me%M2D_Prop%thermS(me%M2D_Prop%Nb_SSports),list_thS_loc(me%M2D_Prop%Nb_SSports))
      idx_prev=0
      do ii=1,me%M2D_Prop%Nb_labelMC
        list_thS_loc(idx_prev+1:idx_prev+me%label%MC_Arr(ii)%NbElem)= &
                  me%label%MC_Arr(ii)%array(1:me%label%MC_Arr(ii)%NbElem)
        idx_prev=idx_prev+me%label%MC_Arr(ii)%NbElem
      enddo
      call krn%update(SS_src_p_loc=list_thS_loc,SS_src_p_view=SS_src_ports)
      do i=1,me%M2D_Prop%Nb_SSports                          
        me%M2D_Prop%thermS(i)%p => SS_src_ports(i)
        me%M2D_Prop%thermS(i)%p%typeS = 'mesh2D'
      enddo
    endif

    allocate(me%bigLS%ValExp(me%M2D_Prop%nnz),me%bigLS%RowExp(me%M2D_Prop%nnz),me%bigLS%ColExp(me%M2D_Prop%nnz))
    allocate(me%bigLS%ValImp(me%M2D_Prop%nnz),me%bigLS%RowImp(me%M2D_Prop%nnz),me%bigLS%ColImp(me%M2D_Prop%nnz))

    idx=1

    PrevColIdx=0
    PrevRowIdx=0
    do nb_e=1,me%M2D_Prop%nb_elements
        me%bigLS%ValExp(idx)%p => me%big%bmx(nb_e)
        me%bigLS%ColExp(idx)=PrevColIdx+1
        me%bigLS%RowExp(idx)=PrevRowIdx+1

        me%bigLS%ValImp(idx)%p => me%big%bmx(nb_e)
        me%bigLS%ColImp(idx)=PrevColIdx+1
        me%bigLS%RowImp(idx)=PrevRowIdx+1

        idx=idx+1
        PrevColIdx=PrevColIdx+1
        PrevRowIdx=PrevRowIdx+1
    enddo

    PrevRowIdx=0
    do nb_e=1,me%M2D_Prop%nb_elements
        do jj=1,3 ! loop over the 3 potential neighbors of cell nb_e
            if(jj==3 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==2) exit
            if(jj==3 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==1) exit  
            if(jj==2 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==1) exit
  
            me%bigLS%ValExp(idx)%p => me%big%dmx(nb_e,jj)
            me%bigLS%ColExp(idx)=me%M2D_Prop%elem(nb_e)%Neigh(jj)
            me%bigLS%RowExp(idx)=PrevRowIdx+1
    
            me%bigLS%ValImp(idx)%p => me%big%dmx(nb_e,jj)
            me%bigLS%ColImp(idx)=me%M2D_Prop%elem(nb_e)%Neigh(jj)
            me%bigLS%RowImp(idx)=PrevRowIdx+1
  
            idx=idx+1
        enddo
        PrevRowIdx=PrevRowIdx+1        
    enddo

    call krn%coo_add(.true.,me%bigLS%ColExp,me%bigLS%RowExp,me%bigLS%ValExp)
    call krn%coo_add(.false.,me%bigLS%ColImp,me%bigLS%RowImp,me%bigLS%ValImp)

    deallocate(me%bigLS%ValExp,me%bigLS%RowExp,me%bigLS%ColExp)
    deallocate(me%bigLS%ValImp,me%bigLS%RowImp,me%bigLS%ColImp)

end subroutine mesh2D_init_part2

end module cmp_mesh2D_init_m