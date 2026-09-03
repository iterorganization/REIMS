! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

!! Standalone dll compilation:
!!     ifx -fpp -fast -dll ext_he_prop.F90 -o he_prop.dll
module krn_global_tools_m
    integer, parameter :: dp = 8
    contains
end module krn_global_tools_m
module krn_simulation_m
    contains
        subroutine set_error(msg)
            character(*), intent(in) :: msg
            print *,"Error: ", msg
            stop 1
        end subroutine set_error
end module krn_simulation_m
#include "../src/lib_ext_math.f90"
#include "../src/lib_He_thermo.f90"
!--------------------------------------------------------------

module ext_he_prop
    use krn_global_tools_m, only: dp
    use krn_simulation_m,   only: set_error
    use lib_He_thermo_m
    implicit none
    
contains

subroutine he_prop_all(type, size, num_cols, table) bind(c, name='he_prop_all')
    !DEC$ ATTRIBUTES DLLEXPORT :: he_prop_all
    use, intrinsic :: iso_c_binding
    integer(c_int), value, intent(in) :: type, size, num_cols
    real(c_double), intent(inout)     :: table(size, num_cols)
    integer :: n

    select case(type)
        case(1) ! r_roT (2 in, 1 out)
            do n = 1, size; table(n,3) = r_roT(table(n,1), table(n,2)); end do     
        case(2) ! ro_pT (2 in, 1 out)
            do n = 1, size; table(n,3) = ro_pT(table(n,1), table(n,2)); end do
        case(3) ! T_roP (2 in, 1 out)
            do n = 1, size; table(n,3) = T_roP(table(n,1), table(n,2)); end do
        case(4) ! T_roE (2 in, 1 out)
            do n = 1, size; table(n,3) = T_roE(table(n,1), table(n,2)); end do
        case(5) ! droeint_droP (2 in, 4 out)
            do n = 1, size
                table(n,3) = droeint_droP(table(n,1), table(n,2), table(n,4), table(n,5), table(n,6))
            end do
        case(6) ! state_roT (2 in, 4 out)
            do n = 1, size
                call state_roT(table(n,1), table(n,2), table(n,3), table(n,4), table(n,5), table(n,6))
            end do
        case(7) ! state_roP (2 in, 4 out)
            do n = 1, size
                call state_roP(table(n,1), table(n,2), table(n,3), table(n,4), table(n,5), table(n,6))
            end do
        case(8) ! state_roP_withR (3 in, 3 out) -> Inputs: ro,p,r (1,2,3) | Outputs: e,T,c (4,5,6)
            do n = 1, size
                call state_roP_withR(table(n,1), table(n,2), table(n,3), table(n,4), table(n,5), table(n,6))
            end do
        case(9) ! state_roE (2 in, 3 out)
            do n = 1, size
                call state_roE(table(n,1), table(n,2), table(n,3), table(n,4), table(n,5))
            end do
        case(10) ! state_roE_withR (3 in, 3 out)
            do n = 1, size
                call state_roE_withR(table(n,1), table(n,2), table(n,3), table(n,4), table(n,5), table(n,6))
            end do
        case(11) ! jacobian_roT (2 in, 10 out)
            do n = 1, size
                call jacobian_roT(table(n,1), table(n,2), table(n,3), table(n,4), table(n,5), table(n,6), &
                                  table(n,7), table(n,8), table(n,9), table(n,10), table(n,11), table(n,12))
            end do
        case(12) ! dc2_roT (2 in, 2 out)
            do n = 1, size
                call dc2_roT(table(n,1), table(n,2), table(n,3), table(n,4))
            end do
        case(13) ! he_prop (2 in, 5 out)
            do n = 1, size
                call he_prop(table(n,1), table(n,2), table(n,3), table(n,4), table(n,5), table(n,6), table(n,7))
            end do
        case default
            call set_error("he_prop_all: unknown type")
    end select
    
end subroutine he_prop_all
    
end module ext_he_prop