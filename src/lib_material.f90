! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module lib_material_m
    use lib_material_nb3sn_m
    use lib_material_nbti_m
    use lib_material_metal_m
    use lib_material_insulation_m
    use lib_input_m
    implicit none

    type material_t
        character(:), allocatable :: material
        real(dp)    :: density = -1.0_dp
        real(dp)    :: E0      = -1.0_dp
        integer     :: nPow    = -1
        type(c_ptr) :: res_cfg      = C_NULL_PTR
        type(c_ptr) :: th_cond_cfg  = C_NULL_PTR
        type(c_ptr) :: stra_cfg     = C_NULL_PTR
        type(c_ptr) :: t_crit_cfg   = C_NULL_PTR
        type(c_ptr) :: bc_cfg       = C_NULL_PTR
        type(c_ptr) :: jc_cfg       = C_NULL_PTR
        type(c_ptr) :: tcs_cfg      = C_NULL_PTR
        type(c_ptr) :: heat_cap_cfg = C_NULL_PTR
        procedure(itf_2var), pointer, nopass :: res      => resistivity_not_implemented
        procedure(itf_2var), pointer, nopass :: th_cond  => thermal_conductivity_not_implemented
        procedure(itf_2var), pointer, nopass :: stra     => strain_not_implemented
        procedure(itf_2var), pointer, nopass :: t_crit   => critical_temperature_not_implemented
        procedure(itf_2var), pointer, nopass :: bc       => critical_field_not_implemented
        procedure(itf_5var), pointer, nopass :: jc       => critical_current_density_not_implemented
        procedure(itf_8var), pointer, nopass :: tcs      => current_sharing_temperature_not_implemented
        procedure(itf_5var), pointer, nopass :: heat_cap => heat_capacity_not_implemented
    contains   
        procedure :: resistivity                  => material_resistivity
        procedure :: resistivity_der              => material_resistivity_der  
        procedure :: thermal_conductivity         => material_thermal_conductivity
        procedure :: thermal_conductivity_der     => material_thermal_conductivity_der  
        procedure :: strain                       => material_strain
        procedure :: critical_temperature         => material_critical_temperature
        procedure :: critical_field               => material_critical_field
        procedure :: critical_current_density     => material_critical_current_density
        procedure :: critical_current_density_der => material_critical_current_density_der  
        procedure :: current_sharing_temperature  => material_current_sharing_temperature    
        procedure :: heat_capacity                => material_heat_capacity
        procedure :: heat_capacity_der            => material_heat_capacity_der          
    end type material_t 
    type(material_t), allocatable, target :: materials(:)
 
    abstract interface
        function init_material_ext_if(cfg_len,cfg_str,c_pt,res_ptr,th_cond_ptr,stra_ptr,t_crit_ptr, &
                                    bc_ptr,jc_ptr,tcs_ptr,heat_cap_ptr,density,E0,nPow) bind(C)
            import :: c_int, c_ptr, c_funptr, dp
            integer(c_int), value :: cfg_len
            character(*)    :: cfg_str 
            type(c_ptr)     :: c_pt      
            type(c_funptr)  :: res_ptr,th_cond_ptr,stra_ptr,t_crit_ptr,bc_ptr,jc_ptr,tcs_ptr,heat_cap_ptr
            real(dp)        :: density,E0
            integer         :: nPow      
            logical         :: init_material_ext_if
        end function
    end interface

    procedure(init_material_ext_if), pointer :: init_material_ext  => null()

contains

    subroutine load_materials(cfgs)
        class(input_t), intent(in) :: cfgs(:)
        
        character(20)  :: names(6) = ['copper','nb3sn','nbti','stainless_steel','glass_kapton_glass','glass_epoxy']
        integer        :: i
        type(c_funptr) :: init_material_ext_ptr

        do i = 1, size(lib_handles)
            init_material_ext_ptr = GetProcAddress(lib_handles(i), 'init_material_ext' // c_null_char)
            if (c_associated(init_material_ext_ptr)) then
                call c_f_procpointer(init_material_ext_ptr, init_material_ext)
                exit
            endif
        enddo

        allocate(materials(size(cfgs) + size(names)))
        do i = 1, size(cfgs)
            call set_material(materials(i), cfgs(i), '')
        enddo 
        do i = 1, size(names)
            call set_material(materials(size(cfgs)+i), empty_cfg(), names(i))
        enddo 
    end subroutine load_materials

    subroutine set_material(mat, cfg, built_in_name)
        type(material_t), intent(inout) :: mat
        class(input_t),   intent(in)    :: cfg
        character(*),     intent(in)    :: built_in_name
        
        type(c_funptr) :: res_ptr, th_cond_ptr, stra_ptr, t_crit_ptr, bc_ptr, jc_ptr, tcs_ptr, heat_cap_ptr
        real(dp) :: dens, E0
        integer  :: nPow,i

        character(:), allocatable :: cfg_str
        type(str_ptr), allocatable :: keys(:)
        type(c_ptr) :: c_cfg
        c_cfg = C_NULL_PTR

        mat%material = cfg%str('material',trim(built_in_name))
        select case (cfg%str('base',mat%material))
        case('copper')
            call material_copper_init(cfg, mat%density, c_cfg)
            mat%res_cfg      =  c_cfg
            mat%res          => resistivity_copper
            mat%th_cond_cfg  =  c_cfg
            mat%th_cond      => thermal_conductivity_copper
            mat%heat_cap_cfg =  c_cfg
            mat%heat_cap     => heat_capacity_copper
        case('nb3sn')
            call material_nb3sn_init(cfg, mat%density, c_cfg)
            call get_power_law_parameters_nb3sn(c_cfg,mat%E0,mat%nPow)
            mat%th_cond_cfg  =  c_cfg
            mat%th_cond      => thermal_conductivity_nb3sn
            mat%stra_cfg     =  c_cfg
            mat%stra         => strain_nb3sn
            mat%t_crit_cfg   =  c_cfg
            mat%t_crit       => critical_temperature_nb3sn
            mat%bc_cfg       =  c_cfg
            mat%bc           => critical_field_nb3sn
            mat%jc_cfg       =  c_cfg
            mat%jc           => critical_current_density_nb3sn
            mat%tcs_cfg      =  c_cfg
            mat%tcs          => current_sharing_temperature_nb3sn
            mat%heat_cap_cfg =  c_cfg
            mat%heat_cap     => heat_capacity_nb3sn
        case('nbti')
            call material_nbti_init(cfg, mat%density, c_cfg)
            call get_power_law_parameters_nbti(c_cfg,mat%E0,mat%nPow)
            mat%th_cond_cfg  =  c_cfg
            mat%th_cond      => thermal_conductivity_nbti
            mat%stra_cfg     =  c_cfg
            mat%stra         => strain_nbti
            mat%t_crit_cfg   =  c_cfg
            mat%t_crit       => critical_temperature_nbti
            mat%bc_cfg       =  c_cfg
            mat%bc           => critical_field_nbti
            mat%jc_cfg       =  c_cfg
            mat%jc           => critical_current_density_nbti
            mat%tcs_cfg      =  c_cfg
            mat%tcs          => current_sharing_temperature_nbti
            mat%heat_cap_cfg =  c_cfg
            mat%heat_cap     => heat_capacity_nbti
        case('stainless_steel')
            call material_stainless_steel_init(cfg, mat%density, c_cfg)
            mat%th_cond_cfg  =  c_cfg
            mat%th_cond      => thermal_conductivity_stainless_steel
            mat%heat_cap_cfg =  c_cfg
            mat%heat_cap     => heat_capacity_stainless_steel
        case('glass_epoxy')
            call material_glass_epoxy_init(cfg, mat%density, c_cfg)
            mat%th_cond_cfg  =  c_cfg
            mat%th_cond      => thermal_conductivity_glass_epoxy
            mat%heat_cap_cfg =  c_cfg
            mat%heat_cap     => heat_capacity_glass_epoxy
        case('glass_kapton_glass')
            call material_glass_kapton_glass_init(cfg, mat%density, c_cfg)
            mat%th_cond_cfg  =  c_cfg
            mat%th_cond      => thermal_conductivity_glass_kapton_glass
            mat%heat_cap_cfg =  c_cfg
            mat%heat_cap     => heat_capacity_glass_kapton_glass               
        end select

        if(.not. associated(init_material_ext)) return

        ! serialization
        keys = cfg%keys1d()
        cfg_str = ''
        do i = 1, size(keys)
            cfg_str =  cfg_str // keys(i)%p // c_null_char // cfg%str(keys(i)%p) // c_null_char
        end do

        res_ptr      = C_NULL_FUNPTR
        th_cond_ptr  = C_NULL_FUNPTR
        stra_ptr     = C_NULL_FUNPTR
        t_crit_ptr   = C_NULL_FUNPTR
        bc_ptr       = C_NULL_FUNPTR
        jc_ptr       = C_NULL_FUNPTR
        tcs_ptr      = C_NULL_FUNPTR
        heat_cap_ptr = C_NULL_FUNPTR
        dens = -1.0_dp
        E0   = -1.0_dp
        nPow = -1
        if(.not. init_material_ext(len(cfg_str), cfg_str, c_cfg, res_ptr, th_cond_ptr, stra_ptr, &
                        t_crit_ptr, bc_ptr, jc_ptr, tcs_ptr, heat_cap_ptr, dens, E0, nPow)) return

        if(dens > 0.0_dp) mat%density = dens
        if(E0   > 0.0_dp) mat%E0      = E0
        if(nPow > 0)      mat%nPow    = nPow          

        if(c_associated(res_ptr)) then
            call c_f_procpointer(res_ptr, mat%res)
            mat%res_cfg = c_cfg
        endif
        if(c_associated(th_cond_ptr)) then 
            call c_f_procpointer(th_cond_ptr, mat%th_cond)
            mat%th_cond_cfg = c_cfg
        endif
        if(c_associated(stra_ptr)) then 
            call c_f_procpointer(stra_ptr, mat%stra)
            mat%stra_cfg = c_cfg
        endif
        if(c_associated(t_crit_ptr)) then 
            call c_f_procpointer(t_crit_ptr, mat%t_crit)
            mat%t_crit_cfg = c_cfg
        endif
        if(c_associated(bc_ptr)) then 
            call c_f_procpointer(bc_ptr, mat%bc)
            mat%bc_cfg = c_cfg
        endif
        if(c_associated(jc_ptr)) then 
            call c_f_procpointer(jc_ptr, mat%jc)
            mat%jc_cfg = c_cfg
        endif
        if(c_associated(tcs_ptr)) then 
            call c_f_procpointer(tcs_ptr, mat%tcs)
            mat%tcs_cfg = c_cfg
        endif
        if(c_associated(heat_cap_ptr)) then 
            call c_f_procpointer(heat_cap_ptr, mat%heat_cap)
            mat%heat_cap_cfg = c_cfg
        endif
        if (mat%density < 0.0_dp) error stop 'Material density not set!'
        
    end subroutine set_material

    subroutine material_init(me,cfg)
        class(material_t), pointer, intent(inout) :: me
        class(input_t), target, intent(in) :: cfg
    
        integer :: i
        character(:), allocatable :: mat   
        
        mat=cfg%str('material')
        do i=1,size(materials)
            if(mat==materials(i)%material) me => materials(i)
        enddo
        if(.NOT.associated(me)) error stop mat//' material does not exist'
    end subroutine material_init

    ! Member functions implementations --- calls to lower-level function pointers
    function material_resistivity(me,T,B)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: T,B
        real(dp) :: material_resistivity
    
        material_resistivity=me%res(me%res_cfg,T,B)
    end function material_resistivity

    function material_resistivity_der(me,T,B)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: T,B
        real(dp) :: material_resistivity_der
        real(dp) :: eps
    
        eps=2.0e-2_dp
        material_resistivity_der=(me%res(me%res_cfg,T+0.5_dp*eps,B) &
                                 -me%res(me%res_cfg,T-0.5_dp*eps,B))/eps
    end function material_resistivity_der

    
    function material_thermal_conductivity(me,T,B)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: T
        real(dp), intent(inout), optional :: B
        real(dp) :: material_thermal_conductivity,B_loc
    
        if(present(B)) then
            B_loc=B
        else
            B_loc=0.0_dp
        endif  
        material_thermal_conductivity=me%th_cond(me%th_cond_cfg,T,B_loc)
    end function material_thermal_conductivity

    function material_thermal_conductivity_der(me,T,B)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: T
        real(dp), intent(inout), optional :: B
        real(dp) :: material_thermal_conductivity_der,B_loc
        real(dp) :: eps
    
        if(present(B)) then
            B_loc=B
        else
            B_loc=0.0_dp
        endif
        eps=2.0e-2_dp
        material_thermal_conductivity_der=(me%th_cond(me%th_cond_cfg,T+0.5_dp*eps,B_loc)- &
                                           me%th_cond(me%th_cond_cfg,T-0.5_dp*eps,B_loc))/eps
    end function material_thermal_conductivity_der   

    function material_strain(me,B,I)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: B,I
        real(dp) :: material_strain
    
        material_strain=me%stra(me%stra_cfg,B,I)
    end function material_strain

    function material_critical_temperature(me,B,St)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: B,St
        real(dp) :: material_critical_temperature
    
        material_critical_temperature=me%t_crit(me%t_crit_cfg,B,St)
    end function material_critical_temperature

    function material_critical_field(me,T,St)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: T,St
        real(dp) :: material_critical_field
    
        material_critical_field=me%bc(me%bc_cfg,T,St)
    end function material_critical_field

    function material_critical_current_density(me,T,B,St,Tc0,Bc)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: T,B,St,Tc0,Bc
        real(dp) :: material_critical_current_density
    
        material_critical_current_density=me%jc(me%jc_cfg,T,B,St,Tc0,Bc)
    end function material_critical_current_density   

    function material_critical_current_density_der(me,T,B,St,Tc0,Bc)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: T,B,St,Tc0,Bc
        real(dp) :: material_critical_current_density_der
        real(dp) :: eps
    
        eps=2.0e-2_dp
        material_critical_current_density_der=(me%jc(me%jc_cfg,T+0.5_dp*eps,B,St,Tc0,Bc) &
                                              -me%jc(me%jc_cfg,T-0.5_dp*eps,B,St,Tc0,Bc))/eps
    end function material_critical_current_density_der   

    function material_current_sharing_temperature(me,B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc 
        real(dp) :: material_current_sharing_temperature
    
        material_current_sharing_temperature=me%tcs(me%tcs_cfg,B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc)
    end function material_current_sharing_temperature   

    function material_heat_capacity(me,T,B,TC,TCS,TC0)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: T
        real(dp), intent(inout), optional :: B,TC,TCS,TC0
        real(dp) :: material_heat_capacity,B_loc,TC_loc,TCS_loc,TC0_loc
    
        if(present(B) .and. present(TC) .and. present(TCS) .and. present(TC0)) then
            B_loc=B; TC_loc=TC; TCS_loc=TCS; TC0_loc=TC0
        else
            B_loc=0.0_dp; TC_loc=0.0_dp; TCS_loc=0.0_dp; TC0_loc=0.0_dp
        endif
        material_heat_capacity=me%heat_cap(me%heat_cap_cfg,T,B_loc,TC_loc,TCS_loc,TC0_loc)
    end function material_heat_capacity

    function material_heat_capacity_der(me,T,B,TC,TCS,TC0)
        class(material_t),  intent(inout) :: me
        real(dp), intent(in) :: T
        real(dp), intent(inout), optional :: B,TC,TCS,TC0
        real(dp) :: material_heat_capacity_der,B_loc,TC_loc,TCS_loc,TC0_loc
        real(dp) :: eps
    
        if(present(B) .and. present(TC) .and. present(TCS) .and. present(TC0)) then
            B_loc=B; TC_loc=TC; TCS_loc=TCS; TC0_loc=TC0
        else
            B_loc=0.0_dp; TC_loc=0.0_dp; TCS_loc=0.0_dp; TC0_loc=0.0_dp
        endif  
        eps=2.0e-2_dp
        material_heat_capacity_der=(me%heat_cap(me%heat_cap_cfg,T+0.5_dp*eps,B_loc,TC_loc,TCS_loc,TC0_loc)-&
                                    me%heat_cap(me%heat_cap_cfg,T-0.5_dp*eps,B_loc,TC_loc,TCS_loc,TC0_loc))/eps
    end function material_heat_capacity_der   

    ! Dummy functions - not implemented, to handle errors when material does not set or override them
    function resistivity_not_implemented(c_pt,T,B) result(ret_val) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: T,B
        real(dp) :: ret_val
        ret_val=0.0_dp
        if(.false.) print*,c_associated(c_pt),T,B
        error stop 'resistivity - not implemented'
    end function resistivity_not_implemented   

    function thermal_conductivity_not_implemented(c_pt,T,B) result(ret_val) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: T,B
        real(dp) :: ret_val
        ret_val=0.0_dp
        if(.false.) print*,c_associated(c_pt),T,B
        error stop 'thermal_conductivity - not implemented'
    end function thermal_conductivity_not_implemented      

    function strain_not_implemented(c_pt,B,I) result(ret_val) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: B,I
        real(dp) :: ret_val
        ret_val=0.0_dp
        if(.false.) print*,c_associated(c_pt),B,I
        error stop 'strain - not implemented'
    end function strain_not_implemented     

    function critical_temperature_not_implemented(c_pt,B,St) result(ret_val) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: B,St
        real(dp) :: ret_val
        ret_val=0.0_dp
        if(.false.) print*,c_associated(c_pt),B,St
        error stop 'critical_temperature - not implemented'
    end function critical_temperature_not_implemented

    function critical_field_not_implemented(c_pt,T,St) result(ret_val) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: T,St
        real(dp) :: ret_val
        ret_val=0.0_dp
        if(.false.) print*,c_associated(c_pt),T,St
        error stop 'critical_field - not implemented'
    end function critical_field_not_implemented

    function critical_current_density_not_implemented(c_pt,T,B,St,Tc0,Bc) result(ret_val) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: T,B,St,Tc0,Bc
        real(dp) :: ret_val
        ret_val=0.0_dp
        if(.false.) print*,c_associated(c_pt),T,B,St,Tc0,Bc
        error stop 'critical_current_density - not implemented'
    end function critical_current_density_not_implemented

    function current_sharing_temperature_not_implemented(c_pt,B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc) result(ret_val) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc
        real(dp) :: ret_val
        ret_val=0.0_dp
        if(.false.) print*,c_associated(c_pt),B,St,Jop,Bc0,Jc0,Tc,Tc0,Bc
        error stop 'current_sharing_temperature - not implemented'
    end function current_sharing_temperature_not_implemented
    
    function heat_capacity_not_implemented(c_pt,T,B,TC,TCS,TC0) result(ret_val) bind(C)
        type(c_ptr), value :: c_pt
        real(dp), value :: T,B,TC,TCS,TC0
        real(dp) :: ret_val
        ret_val=0.0_dp
        if(.false.) print*,c_associated(c_pt),T,B,TC,TCS,TC0
        error stop 'heat_capacity - not implemented'
    end function heat_capacity_not_implemented   

end module lib_material_m