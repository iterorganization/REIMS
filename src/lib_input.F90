! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module lib_input_m
    !! Library to parse yaml input files
    use hdf5
    use h5lt
    use krn_global_tools_m
    use fortran_yaml_c, only: YamlFile,  type_key_value_pair, &
                              type_node, type_dictionary, type_error, &
                              type_list, type_list_item,  type_scalar
    implicit none
    private
    public :: input_t, get0d_hdf5, get1d_hdf5, get2d_hdf5, continue_if, empty_cfg

    type, extends(type_node) :: input_t
        class(type_node), pointer :: root => null()       !! top level (main) access to the tree
        type(YamlFile)            :: file                 !! used only for top level config
        type(type_dictionary)     :: cmp_by_type          !! dictionary with arrays of components
        type(input_t), pointer    :: cmp_arr(:) => null() !! array of components of the same type
    contains
        procedure :: open    => input_open_yaml  !! Open and parse file and return as input dictionary 
        ! Main get access functions
        procedure :: dict    => input_dict   !! Get input sub-dictionary
        procedure :: dict1d  => input_dict1d !! Get allocatable array of input sub-dictionaries 
        procedure :: dbl     => input_dbl    !! Get 'real' double precision value
        procedure :: dbl1d   => input_dbl1d  !! Get 'real' double precision allocatable one dimension array
        procedure :: dbl2d   => input_dbl2d  !! Get 'real' double precision allocatable two dimensional array
        procedure :: bin     => input_bin    !! Get 'logical' value like: true/false yes/no etc... 
        procedure :: bin1d   => input_bin1d  !! Get 'logical' 1d array of values of true/false yes/no etc... 
        procedure :: int     => input_int    !! Get 'integer' value
        procedure :: int1d   => input_int1d  !! Get 'integer' value
        procedure :: str     => input_str    !! Get allocatable trimmed 'character(:)' value
        procedure :: str1d   => input_str1d  !! Get allocatable trimmed 'character(:)' value
        procedure :: cmp1d   => input_cmp1d  !! Get pointer to array of input components
        procedure :: nb_cmp  => input_nb_cmp !! Get total number of components
        procedure :: has_key => input_has_key!! Return true if dictionary has a key
        procedure :: keys1d  => input_keys1d !! Return 1d array of keys
        ! Private only for debugging
        procedure :: node    => input_node   !! Get "type_node" from yaml lib
        procedure :: dump    => input_dump_dict
    end type input_t

contains
! ---------------------------------------------------------------------------
!                           Accessors Functions
! ---------------------------------------------------------------------------
function input_cmp1d(me,keys) result(return_array)
    class(input_t), intent(in) :: me
    character(*),   intent(in) :: keys(:)
    type(input_t),     pointer :: return_array(:)
    class(type_node),  pointer :: node
    integer :: i,j,nb_cmp

    nb_cmp = 0
    do i = 1, size(keys)
        node => me%cmp_by_type%get(trim(keys(i)))
        if (.not.associated(node)) continue
        select type(node); class is (input_t)
            nb_cmp = nb_cmp + size(node%cmp_arr)
        end select
    enddo
    allocate(return_array(nb_cmp))
    nb_cmp = 0
    do i = 1, size(keys)
        node => me%cmp_by_type%get(trim(keys(i)))
        if (.not.associated(node)) continue
        select type(node); class is (input_t)
            do j = 1, size(node%cmp_arr)
                node%cmp_arr(j)%path = trim(keys(i))
            enddo
            return_array(nb_cmp+1:nb_cmp+size(node%cmp_arr)) = node%cmp_arr
            nb_cmp = nb_cmp + size(node%cmp_arr)
        end select
    enddo
end function input_cmp1d

function input_nb_cmp(me) result(nb_cmp)
    class(input_t), intent(in) :: me
    integer :: nb_cmp
    type(type_key_value_pair), pointer :: cmp

    nb_cmp = 0
    cmp => me%cmp_by_type%first
    do while(associated(cmp))
        select type(cmp_array => cmp%value); class is(input_t)
            nb_cmp = nb_cmp + size(cmp_array%cmp_arr)
        end select
        cmp => cmp%next
    enddo
end function input_nb_cmp

function input_dict(me,key,filter_key,filter_value,empty) result(out_dict)
    class(input_t),         intent(in) :: me
    character(*),           intent(in) :: key
    character(*), optional, intent(in) :: filter_key, filter_value
    logical,      optional, intent(in) :: empty
    class(input_t), allocatable        :: out_dict

    class(type_dictionary), allocatable, target :: empty_dict
    type(type_dictionary)       :: dict
    class(input_t), allocatable :: dict_array(:)

    allocate(out_dict)
    if (present(filter_key) .and. present(filter_value)) then
        dict_array = input_dict1d(me,key,filter_key,filter_value)
        if (size(dict_array) == 0) then
            error stop "Couldn't match key,value: "//filter_key//','//filter_value
        endif
        if (size(dict_array) > 1) then
            error stop "Key,value pair match multiple times: "//filter_key// &
                            ','//filter_value
        endif
        out_dict%root => dict_array(1)%root
    else
        if(present(empty) .and. .not. me%has_key(key)) then
            if(empty) then
                allocate(empty_dict)
                out_dict%root => empty_dict
                return
            endif
        endif
        out_dict%root => me%node(key)
        if (.not.same_type_as(out_dict%root,dict)) &
            error stop 'Value in key: "'//key//'" is not dictionary.'
    endif
end function input_dict

function input_dict1d(me,key,filter_key,filter_value,empty) result(return_array)
    class(input_t),         intent(in) :: me
    character(*),           intent(in) :: key
    character(*), optional, intent(in) :: filter_key, filter_value
    class(input_t), allocatable        :: dict_array(:), return_array(:)
    logical,      optional, intent(in) :: empty

    type(type_error),  allocatable :: err
    class(type_list),      pointer :: list
    type(type_dictionary), pointer :: element
    type(type_list_item),  pointer :: item
    integer :: i
    
    if(present(empty)) then
        if(.not.me%has_key(key) .and. empty) then
            allocate(return_array(0))
            return
        endif
    endif

    list => select_list(me%node(key))
    allocate(dict_array(list%size()))

    i = 0
    item => list%first
    do while(associated(item))
        element => select_dict(item%node);
        if (.not. present(filter_key) .or. &
            element%get_string(filter_key,'',err) == filter_value) then
            i = i + 1
            dict_array(i)%root => element
        endif
        item => item%next
    enddo
    allocate(return_array(i))
    return_array = dict_array(:i)
end function input_dict1d

function input_dbl2d(me,key,default_value) result(val)
    class(input_t), intent(in) :: me
    character(*),    intent(in) :: key
    real(dp), optional, intent(in) :: default_value(:,:) !! default value in case of missing key
    real(dp),  allocatable :: val(:,:)

    type(type_list),       pointer :: element, element1
    type(type_scalar),     pointer :: element2
    type(type_list_item),  pointer :: item1, item2
    integer :: i,j,dim1,dim2
    logical :: ok
    type(input_t) :: tmp

    if(present(default_value)) val = default_value
    if(.not.me%has_key(key) .and. present(default_value)) return
    select type(node=>me%node(key));class is(type_list)
        ! Setting the first element and shaping matrix accordingly 
        dim1 = node%size()
        item1 => node%first
        element => select_list(item1%node)
        dim2 = element%size()
        item2 => element%first
        if (allocated(val)) deallocate(val)
        allocate(val(dim1,dim2))

        i = 1; j = 1
        do while(associated(item1))
            element1 => select_list(item1%node)
            if(element1%size()/=dim2) error stop &
                    "2D array size is inconsistent for: "//element1%path
            j = 1
            item2 => element1%first
            do while(associated(item2))
                element2 => select_scalar(item2%node);
                val(i,j) = element2%to_real(0.0_dp,ok)
                call continue_if(element2,ok)
                item2 => item2%next
                j = j + 1
            enddo
            item1 => item1%next
            i = i + 1
        enddo
        return
        class is(type_dictionary)
        tmp%root => node
        if(get2d_hdf5(tmp,val)) return
    end select
    if(present(default_value)) return
    error stop 'Can not convert to 2d array of "double" for key: '//key
end function input_dbl2d

function input_dbl1d(me,key,default_value,n) result(val)
    class(input_t), intent(in) :: me
    character(*),    intent(in) :: key
    real(dp), optional, intent(in) :: default_value(:) !! default value in case of missing key
    integer,  optional, intent(in) :: n 
    real(dp),  allocatable :: val(:)

    type(type_list_item),  pointer :: item
    type(type_scalar),     pointer :: element
    integer :: i
    logical :: ok
    type(input_t) :: tmp

    if(present(default_value)) val = default_value
    if(.not.me%has_key(key) .and. present(default_value)) return
    select type(node=>me%node(key))
    class is(type_scalar)
        if (.not. present(n)) &
            error stop 'input_dbl1d: scalar node requires "n" for key: '//key
        if (allocated(val)) deallocate(val)
        allocate(val(n))
        val = node%to_real(0.0_dp, ok)
        call continue_if(node, ok)
        return
    class is(type_list)
        if (allocated(val)) deallocate(val)
        allocate(val(node%size()))
        item => node%first
        i = 1
        do while(associated(item))
            element => select_scalar(item%node);
            val(i) = element%to_real(0.0_dp,ok)
            call continue_if(element,ok)
            item => item%next
            i = i + 1
        enddo
        if (present(n)) then
            if (size(val) /= n) then
                print*, 'Expected '//to_str(n)//' elements for "'//key//'" but YAML has '//to_str(size(val))
                error stop 'Array size mismatch in YAML input.'
            endif
        endif
        return
    class is(type_dictionary)
        tmp%root => node
        if(get1d_hdf5(tmp,val)) return
    end select
    if(present(default_value)) return
    error stop 'Can not convert to 1d array of "double" for key: '//key
end function input_dbl1d

function get1d_range(dict,val) result(ok)
    type(input_t), intent(in) :: dict
    integer, intent(inout), allocatable :: val(:)
    logical :: ok

    integer :: i

    ok = .false.
    if (.not.dict%has_key('to')) return
    val = [(i, i = dict%int('from',1),dict%int('to'),dict%int('step',1))]
    ok = .true.
end function get1d_range

function input_int1d(me,key,default_value,n) result(val)
    class(input_t), intent(in) :: me
    character(*),    intent(in) :: key
    integer, optional, intent(in) :: default_value(:) !! default value in case of missing key
    integer, optional, intent(in) :: n
    integer, allocatable :: val(:)

    type(type_list_item),  pointer :: item
    type(type_scalar),     pointer :: element
    integer :: i
    logical :: ok
    type(input_t) :: tmp

    if(present(default_value)) val = default_value
    if(.not.me%has_key(key) .and. present(default_value)) return
    select type(node=>me%node(key))
    class is(type_scalar)
        if (.not. present(n)) &
            error stop 'input_int1d: scalar node requires "n" for key: '//key
        if (allocated(val)) deallocate(val)
        allocate(val(n))
        val = node%to_integer(0, ok)
        call continue_if(node, ok)
        return
    class is(type_list)
        if (allocated(val)) deallocate(val)
        allocate(val(node%size()))
        item => node%first
        i = 1
        do while(associated(item))
            element => select_scalar(item%node);
            val(i) = element%to_integer(0,ok)
            call continue_if(element,ok)
            item => item%next
            i = i + 1
        enddo
        if (present(n)) then
            if (size(val) /= n) then
                print*, 'Expected '//to_str(n)//' elements for "'//key//'" but YAML has '//to_str(size(val))
                error stop 'Array size mismatch in YAML input.'
            endif
        endif
        return
    class is(type_dictionary)
        tmp%root => node
        if(get1d_range(tmp,val)) return
    end select
    if(present(default_value)) return
    error stop 'Can not convert to 1d array of "integer" for key: '//key
end function input_int1d

function input_dbl(me,key,default_value) result(val)
    !! Function returns `double` value given by the key and stops if key doesn't exist.
    !!
    !! Providing one of the optional arguments will not stop program in case of missing
    !! key, default return value is in this case is NaN
    class(input_t),     intent(in)    :: me
    character(*),       intent(in)    :: key !! key to the value in dictionary
    real(dp), optional, intent(in)    :: default_value !! default value in case of missing key
    real(dp)                          :: val !! returned value
    logical :: ok
    type(input_t) :: tmp

    if(present(default_value)) val = default_value
    if(.not.me%has_key(key) .and. present(default_value)) return
    select type(node=>me%node(key));class is(type_scalar)
        val = node%to_real(val, ok)
        call continue_if(node, ok)
        return
    class is(type_dictionary)
        tmp%root => node
        if(get0d_hdf5(tmp,val)) return
    end select
    if(present(default_value)) return
    error stop 'Can not convert to value of "double" for key: '//key
end function input_dbl

function input_str(me,key,default_value) result(val)
    class(input_t), intent(in) :: me
    character(*), intent(in)   :: key
    character(*), optional, intent(in) :: default_value
    character(:), allocatable  :: val
    type(type_scalar), pointer :: scalar

    if(present(default_value)) val = default_value
    if(.not.me%has_key(key) .and. present(default_value)) return
    scalar => select_scalar(me%node(key))
    val = scalar%string
end function input_str

function input_str1d(me,key,default_value,empty) result(val)
    class(input_t), intent(in) :: me
    character(*), intent(in)   :: key
    type(str_ptr), optional, intent(in) :: default_value
    type(str_ptr), allocatable :: val(:)
    logical,      optional, intent(in) :: empty    

    type(type_list_item), pointer :: item
    type(type_scalar),    pointer :: element
    integer :: i

    if(present(empty)) then
        if(.not. me%has_key(key) .and. empty) then
            allocate(val(0))
            return
        endif
    endif

    if(present(default_value)) val = default_value
    if(.not.me%has_key(key) .and. present(default_value)) return
    select type(node=>me%node(key));class is(type_list)
        if (allocated(val)) deallocate(val)
        allocate(val(node%size()))
        item => node%first
        i = 1
        do while(associated(item))
            element => select_scalar(item%node);
            val(i)%p => element%string
            item => item%next
            i = i + 1
        enddo
        return
    end select
    if(present(default_value)) return
    error stop 'Can not convert to 1d array of "string" for key: '//key
end function input_str1d

function input_has_key(me,key)
    class(input_t), intent(in) :: me
    character(*), intent(in) :: key
    logical :: input_has_key

    type(type_dictionary), pointer :: dict
    class(type_node), pointer :: node
    node => null()
    dict => select_dict(me%root)
    node => dict%get(key)
    input_has_key = associated(node)    
end function input_has_key

function input_bin(me,key,default_value) result(val)
    class(input_t), intent(in) :: me
    character(*), intent(in)   :: key
    logical, optional, intent(in) :: default_value
    logical                    :: val
    logical                    :: ok
    type(type_scalar), pointer :: scalar

    if(present(default_value)) val = default_value
    if(.not.me%has_key(key) .and. present(default_value)) return
    scalar => select_scalar(me%node(key))
    val = scalar%to_logical(.true.,ok)
    call continue_if(scalar, ok)
end function input_bin

function input_bin1d(me,key,default_value,n) result(val)
    class(input_t), intent(in) :: me
    character(*), intent(in)   :: key
    logical, optional, intent(in) :: default_value(:) !! default value in case of missing key
    integer, optional, intent(in) :: n
    logical, allocatable :: val(:)

    type(type_list_item), pointer :: item
    type(type_scalar),    pointer :: element
    logical :: ok
    integer :: i

    if(present(default_value)) val = default_value
    if(.not.me%has_key(key) .and. present(default_value)) return
    select type(node=>me%node(key))
    class is(type_scalar)
        if (.not. present(n)) &
            error stop 'input_bin1d: scalar node requires "n" for key: '//key
        if (allocated(val)) deallocate(val)
        allocate(val(n))
        val = node%to_logical(.true., ok)
        call continue_if(node, ok)
        return
    class is(type_list)
        if (allocated(val)) deallocate(val)
        allocate(val(node%size()))
        item => node%first
        i = 1
        do while(associated(item))
            element => select_scalar(item%node);
            val(i) = element%to_logical(.true.,ok)
            call continue_if(element, ok)
            item => item%next
            i = i + 1
        enddo
        if (present(n)) then
            if (size(val) /= n) then
                print*, 'Expected '//to_str(n)//' elements for "'//key//'" but YAML has '//to_str(size(val))
                error stop 'Array size mismatch in YAML input.'
            endif
        endif
        return
    end select
    if(present(default_value)) return
    error stop 'Can not convert to 1d array of "logical" for key: '//key
end function input_bin1d

function input_int(me,key,default_value) result(val)
    class(input_t), intent(in) :: me
    character(*), intent(in)   :: key
    integer, optional, intent(in) :: default_value
    integer                    :: val
    logical                    :: ok
    type(type_scalar), pointer :: scalar

    if(present(default_value)) val = default_value
    if(.not.me%has_key(key) .and. present(default_value)) return
    scalar => select_scalar(me%node(key))
    val = scalar%to_integer(0, ok)
    call continue_if(scalar, ok)
end function input_int

recursive subroutine input_open_yaml(me,file_name,included)
    class(input_t), intent(inout) :: me
    character(*),      intent(in) :: file_name
    logical, optional, intent(in) :: included  !! ignore it - only present for included files

    class(type_node), pointer :: cmp_tree
    character(:), allocatable :: error,folder
    type(type_dictionary)     :: no_cmp
    type(type_dictionary), pointer :: root_dict
    logical :: file_exists

    folder = file_name(:index(file_name,'/',.true.))

    inquire(file=file_name, exist=file_exists)
    if (.not. file_exists) error stop 'Input file not found: '//file_name

    call me%file%parse(file_name, error)
    if (allocated(error)) error stop 'YAML parse error in "'//file_name//'": '//error

    me%root => me%file%root
    call scan_yaml(me%file%root,folder)

    if (.not.present(included)) then
        root_dict => select_dict(me%root);
        ! It is top level input file
        cmp_tree => root_dict%get('components')
        call rename_components(cmp_tree,no_cmp,'','')
        call allocate_components(me,no_cmp)
        call fill_components(me,cmp_tree,no_cmp)
    endif    
end subroutine input_open_yaml

function input_keys1d(me) result(keys)
    class(input_t), intent(in) :: me
    type(str_ptr), allocatable :: keys(:)

    class(type_dictionary),     pointer :: dict
    class(type_key_value_pair), pointer :: kvp
    integer :: i
    
    dict => select_dict(me%root)
    kvp => dict%first
    allocate(keys(dict%size()))
    do i = 1, size(keys)
        keys(i)%p => kvp%key
        kvp => kvp%next
    enddo
end function input_keys1d

! ---------------------------------------------------------------------------
!                           Helper Functions
! ---------------------------------------------------------------------------
recursive subroutine input_dump_dict(self,unit,indent)
    class(input_t), intent(in) :: self
    integer, intent(in) :: unit, indent
    integer :: i
    type(type_key_value_pair),pointer :: pair
    if (associated(self%cmp_by_type%first)) then
        write(unit,*) '----------------- Components by type -------------------'
    endif
    pair=>self%cmp_by_type%first
    do while(associated(pair))
        select type(arr=>pair%value);class is(input_t)
            write(unit,'(a,i4,a)')'>>',size(arr%cmp_arr), &
                ' components of type "'//pair%key//'"'
            do i=1,size(arr%cmp_arr)
                write(unit,'(a)',advance='NO') '  - '
                call arr%cmp_arr(i)%dump(unit,indent+4)
            enddo
        end select
        pair=>pair%next
    enddo
    if (associated(self%cmp_by_type%first)) then
        write(unit,*) '------------------ Main dictionary ---------------------'
    endif
    call self%root%dump(unit,indent)
end subroutine input_dump_dict

subroutine abort_if(err,complementary_message)
    type(type_error), allocatable, intent(in) :: err
    character(*), optional :: complementary_message
    if (allocated(err)) then
        if (.not. present(complementary_message)) then
            complementary_message = ''
        endif
        write(stderr,*) err%message,complementary_message
        stop 1
    endif
end subroutine abort_if

subroutine continue_if(element,ok)
    type(type_scalar), intent(in) :: element
    logical, intent(in) :: ok
    if (.not. ok) then
        error stop "Can't convert value at: '"//element%path//"' value: '"// &
                        element%string//"'"
    endif
end subroutine continue_if

recursive subroutine dictionary_traverse(dict,key)
    !! returns sub-dictionary by traversing from top level dictionary using key 
    !! separation sign: '/'
    type(type_dictionary), pointer, intent(inout) :: dict
    character(*),                   intent(inout) :: key

    integer                        :: key_len
    type(type_error), allocatable  :: err

    key_len = index(key,'/')
    if (key_len /= 0) then
        dict => dict%get_dictionary(key(:key_len-1),.true.,err)
        call abort_if(err,' Key: '//key)
        key = key(key_len+1:)
        call dictionary_traverse(dict,key)
    endif
end subroutine dictionary_traverse

function select_dict(me)
    !! take node and return dictionary, if it's not dictionary it's stops 
    class(type_node), pointer, intent(in) :: me
    type(type_dictionary), pointer :: select_dict
    select type(me)
    class is (type_dictionary)
        select_dict => me
    class default
        error stop 'It is not a dictionary: '//me%path
    end select
end function

function select_list(me)
    !! take node and return list, if it's not list it's stops 
    class(type_node), pointer, intent(in) :: me
    type(type_list), pointer :: select_list
    select type(me)
    class is (type_list)
        select_list => me
    class default
        error stop 'It is not a list: '//me%path
    end select
end function

function select_scalar(me)
    !! take node and return scalar, if it's not scalar it's stops 
    class(type_node), pointer, intent(in) :: me
    type(type_scalar), pointer :: select_scalar
    select type(me)
    class is (type_scalar)
        select_scalar => me
    class default
        error stop 'It is not a scalar: '//me%path
    end select
end function

! ---------------------------------------------------------------------------
!                       Tree processing Functions
! ---------------------------------------------------------------------------

recursive subroutine scan_yaml(me,folder,dict_item,list_item)
    !! Scanning all nodes in yaml tree and attaching to them external files
    !! It also implementing dictionary key marge operator "<<"
    class(type_node),                     intent(inout) :: me
    character(*),                         intent(inout) :: folder
    class(type_key_value_pair), optional, intent(inout) :: dict_item
    class(type_list_item),      optional, intent(inout) :: list_item
    ! TODO: think about having parent (dict_item or list_item) as type_node,
    !       maybe it will make this code nicer.

    character(:),          allocatable :: path
    type(type_error),      allocatable :: err
    type(type_key_value_pair), pointer :: pair
    type(type_list_item),      pointer :: item
    type(type_dictionary),     pointer :: dict
    type(type_list),           pointer :: list
    type(type_list),            target :: tmp_list
    type(input_t) :: file
    
    ! including files external files
    select type(dict_top=>me);class is(type_dictionary)
        path = dict_top%get_string('include','',err)
        if(path/='') then
            call file%open(folder//path,.true.)
            if (present(dict_item)) then
                dict_item%value => file%root
            endif
            if (present(list_item)) then
                list_item%node => file%root
            endif
            nullify(file%file%root) ! to avoid deallocation
        endif
    end select

    ! recursive scan through whole tree
    select type(me)
    class is(type_dictionary)
        pair => me%first
        do while (associated(pair))
            select type (value=>pair%value)
            class is (type_dictionary)
                call scan_yaml(pair%value,folder,dict_item=pair)
            class is (type_list)
                call scan_yaml(pair%value,folder,dict_item=pair)
            class default;end select
            pair => pair%next
        end do
        class is(type_list)
        item => me%first
        do while (associated(item))
            select type (value=>item%node)
            class is (type_dictionary)
                call scan_yaml(item%node,folder,list_item=item)
            class is (type_list)
                call scan_yaml(item%node,folder,list_item=item)
            class default;end select
            item => item%next
        end do
    end select

    ! key merge with operator "<<"
    merge: select type(dict_top=>me);class is(type_dictionary)
        dict => dict_top%get_dictionary('<<',.false.,error=err)
        if(associated(dict)) then
            call tmp_list%append(dict)
            list => tmp_list
        else
            list => dict_top%get_list('<<',.true.,error=err)
        endif
        if(allocated(err)) exit merge
        if(list%size()==0) exit merge 
        select type(first_dict=>list%first%node);class is(type_dictionary)
            item => list%first
            item => item%next
            do while(associated(item))
                select type(sub_dict=>item%node);class is(type_dictionary)
                    pair => sub_dict%first
                    do while (associated(pair))
                        call first_dict%set(pair%key,pair%value)
                        pair => pair%next
                    enddo
                end select
                item => item%next
            enddo

            pair => dict_top%first
            do while (associated(pair))
                if(pair%key/='<<')then
                    call first_dict%set(pair%key,pair%value)
                endif
                pair => pair%next
            enddo

            if (present(dict_item)) then
                dict_item%value => first_dict
            endif
            if (present(list_item)) then
                list_item%node => first_dict
            endif
        end select
        nullify(tmp_list%first)
    end select merge
end subroutine scan_yaml

recursive subroutine rename_components(cmp_tree,no_cmp,prefix,postfix)
    class(type_node), pointer, intent(inout) :: cmp_tree
    type(type_dictionary),     intent(inout) :: no_cmp
    character(*),              intent(in)    :: prefix, postfix

    type(type_error),allocatable :: err
    type(type_key_value_pair), pointer :: pair
    type(type_list_item), pointer :: item
    character(:), allocatable :: prefix1, postfix1, type_name
    type(type_dictionary),pointer :: cmp
    integer :: ln

    select type(cmp_tree)
    class is(type_dictionary)
        pair => cmp_tree%first
        do while(associated(pair))
            ln = len(pair%key)
            prefix1 = prefix
            postfix1 = postfix
            if (ln>2 .and. pair%key(1:1)=='_') then
                postfix1 = postfix//pair%key(2:ln)
            endif
            if (ln>2 .and. pair%key(ln:ln)=='_') then
                prefix1 = pair%key(1:ln-1)//prefix
            endif
            if(pair%key(1:1)/='.') then
                call rename_components(pair%value,no_cmp,prefix1,postfix1)
            endif
            pair=>pair%next
        enddo
    class is(type_list)
        item => cmp_tree%first
        do while(associated(item))
            ! add elements to no_cmp an calculate size
            cmp => select_dict(item%node)
            type_name = cmp%get_string('type',error=err)
            call abort_if(err)
            ln = no_cmp%get_integer(type_name,-1,err)
            if(ln == -1) then
                call no_cmp%set_string(type_name,'0')
            endif
            ln = increment_int_dict(no_cmp,type_name)

            ! renaming main id an connection id
            call rename_id(cmp,prefix,postfix)
    
            item => item%next
        enddo
    end select
end subroutine rename_components

recursive subroutine allocate_components(me,no_cmp)
    type(input_t),            intent(inout) :: me
    type(type_dictionary),     intent(inout) :: no_cmp
    class(input_t), pointer            :: cmp
    class(type_node), pointer            :: cmp_node
    type(type_key_value_pair), pointer :: pair
    integer :: size
    pair => no_cmp%first
    do while(associated(pair))
        select type(size_node => pair%value); class is(type_scalar)
            size = size_node%to_integer(0)
            size_node%string = '0'
            allocate(cmp)
            allocate(cmp%cmp_arr(size))
            cmp_node => cmp
            call me%cmp_by_type%set(pair%key,cmp_node)
        end select
        pair => pair%next
    enddo
end subroutine allocate_components

recursive subroutine fill_components(me,tree,no_cmp)
    type(input_t),            intent(inout) :: me
    class(type_node), pointer, intent(in)    :: tree
    type(type_dictionary),     intent(inout) :: no_cmp

    type(input_t),pointer :: new_dict
    class(type_node),pointer     :: tmp
    type(type_error),allocatable :: err
    type(type_key_value_pair), pointer :: pair
    type(type_list_item), pointer :: item
    character(:), allocatable :: type
    type(type_dictionary),pointer :: cmp
    integer :: nb

    select type(tree)
    class is(type_dictionary)
        pair => tree%first
        do while(associated(pair))
            if(pair%key(1:1)/='!') then
                call fill_components(me,pair%value,no_cmp)
            endif
            pair=>pair%next
        enddo
    class is(type_list)
        item => tree%first
        do while(associated(item))
            ! Adding elements to array in dictionary and increase size
            cmp => select_dict(item%node)
            type = cmp%get_string('type',error=err)
            call abort_if(err)
            allocate(new_dict)
            new_dict%root => cmp
            tmp => me%cmp_by_type%get(type)
            nb = increment_int_dict(no_cmp,type)
            select type(tmp); class is(input_t)
                call new_dict%root%set_path(type)
                tmp%cmp_arr(nb) = new_dict
            end select
            item => item%next
            nullify(new_dict)
        enddo
    end select
end subroutine fill_components

function increment_int_dict(dict,key) result(val)
    type(type_dictionary), intent(inout) :: dict
    character(*),intent(in) :: key
    integer :: val
    type(type_error),allocatable :: err
    val = dict%get_integer(key,error=err) + 1
    call abort_if(err)
    call dict%set_string(key,to_str(val))
end function increment_int_dict

recursive subroutine rename_id(cmp,prefix,postfix)
    type(type_dictionary), intent(inout) :: cmp
    character(*), intent(in) :: prefix,postfix
    type(type_error),allocatable :: err
    character(:), allocatable :: cmp_id
    type(type_dictionary),pointer :: dict
    type(type_list),pointer       :: list
    type(type_list_item),pointer  :: list_item
    list => null()
    dict => null()

    cmp_id = cmp%get_string('id','',err)    
    if (cmp_id/='') then
        cmp_id = prefix//cmp_id//postfix
        call cmp%set_string('id',cmp_id)
    endif

    dict => cmp%get_dictionary('port',.true.,err)
    if (allocated(err)) then; deallocate(err)
    else; call rename_id(dict,prefix,postfix); endif

    dict => cmp%get_dictionary('in',.true.,err)
    if (allocated(err)) then; deallocate(err)
    else; call rename_id(dict,prefix,postfix); endif
    
    dict => cmp%get_dictionary('out',.true.,err)
    if (allocated(err)) then; deallocate(err)
    else; call rename_id(dict,prefix,postfix); endif
    
    list => cmp%get_list('branches',.false.,err)
    if (.not.associated(list)) then
        list => cmp%get_list('link',.false.,err)
        if (.not.associated(list)) then
            return
        endif
    endif
    list_item => list%first
    do while(associated(list_item))
        dict => select_dict(list_item%node)
        call rename_id(dict,prefix,postfix)
        list_item => list_item%next
    enddo
end subroutine rename_id

function input_node(me,key) result(node)
    class(input_t), intent(in) :: me
    character(*),   intent(in) :: key
    class(type_node),  pointer :: node
    character(:),     allocatable  :: new_key
    type(type_dictionary), pointer :: dict

    new_key = key
    dict => select_dict(me%root)
    call dictionary_traverse(dict,new_key)
    node => dict%get(new_key)
    if (.not. associated(node)) error stop 'Missing required key "'//new_key//'" in '//dict%path
end function input_node

function get0d_hdf5(cfg, sig) result(ok)
    class(input_t), intent(in) :: cfg
    real(dp),    intent(inout) :: sig
    logical :: ok

    real(dp) :: sig_tmp(1)
    integer(hId_t) :: file_id, dSet_id
    integer(hSize_t) :: dims(7), type_size
    integer :: error, type_class
    character(:), allocatable :: dSet
    ok = .false.
    if(.not.cfg%has_key('h5')) return
    dSet = cfg%str('data')
    call h5Open_f(error)
    call h5fOpen_f(cfg%str('h5'), H5F_ACC_RDONLY_F, file_id, error)
    call h5ltGet_dataset_info_f(file_id, dSet, dims, type_class, type_size, error)
    call h5dOpen_f(file_id, dSet, dSet_id, error)
    call h5dRead_f(dSet_id, H5T_NATIVE_DOUBLE, sig_tmp, dims, error)
    call h5dClose_f(dSet_id, error)
    call h5fClose_f(file_id, error)
    call h5close_f(error)
    sig = sig_tmp(1)
    ok = .true.
end function get0d_hdf5


function get1d_hdf5(cfg, sig) result(ok)
    class(input_t), intent(in) :: cfg
    real(dp), allocatable, intent(inout) :: sig(:)
    logical :: ok

    integer(hId_t) :: file_id, dSet_id
    integer(hSize_t) :: dims(7), type_size
    integer :: error, type_class
    character(:), allocatable :: dSet

    ok = .false.
    if(.not.cfg%has_key('h5')) return
    dSet = cfg%str('data')
    call h5fOpen_f(cfg%str('h5'), H5F_ACC_RDONLY_F, file_id, error)
    call h5ltGet_dataset_info_f(file_id, dSet, dims, type_class, type_size, error)
    allocate(sig(dims(1)))
    call h5dOpen_f(file_id, dSet, dSet_id, error)
    call h5dRead_f(dSet_id, H5T_NATIVE_DOUBLE, sig, dims, error)
    call h5dClose_f(dSet_id, error)
    call h5fClose_f(file_id, error)
    ok = .true.
end function get1d_hdf5

function get2d_hdf5(cfg, sig) result(ok)
    class(input_t), intent(in) :: cfg
    real(dp), allocatable :: sig(:,:)
    logical :: ok

    integer(hId_t) :: file_id, dSet_id
    integer(hSize_t) :: dims(7), type_size
    integer :: error, type_class
    character(:), allocatable :: dSet

    ok = .false.
    if(.not.cfg%has_key('h5')) return
    dSet = cfg%str('data')
    call h5fOpen_f(cfg%str('h5'), H5F_ACC_RDONLY_F, file_id, error)
    call h5ltGet_dataset_info_f(file_id, dSet, dims, type_class, type_size, error)
    allocate(sig(dims(1),dims(2)))
    call h5dOpen_f(file_id, dSet, dSet_id, error)
    call h5dRead_f(dSet_id, H5T_NATIVE_DOUBLE, sig, dims, error)
    call h5dClose_f(dSet_id, error)
    call h5fClose_f(file_id, error)
    ok = .true.
end function get2d_hdf5

function empty_cfg()
        class(input_t),         allocatable         :: empty_cfg
        class(type_dictionary), allocatable, target :: empty_dict     
        allocate(empty_dict,empty_cfg)
        empty_cfg%root => empty_dict
end function empty_cfg

end module lib_input_m
