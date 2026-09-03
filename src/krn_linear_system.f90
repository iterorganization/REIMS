! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later


module krn_linear_system_m
    use krn_interface_m
    use krn_simulation_m, only: set_error
    implicit none
contains
    
subroutine Pardiso_fact_solve(krn,for_explicit,TmS)
    type(krn_t),intent(inout),target :: krn
    logical, intent(in) :: for_explicit
    type(timer_t(*)), intent(inout) :: TmS

    type(coo_csr_t(:)),pointer :: sch
    integer :: msglvl,error,neq,i

    sch=>krn%imp
    if(for_explicit) sch=>krn%exp          

    neq=size(krn%rhs_or_solution)
    msglvl=0
    sch%perm(:)=0

    sch%csr_val=from_ptr_arr(sch%csr_ptr)

    if (.not. sch%analysed) then
        call pardisoinit(sch%pt, 11, sch%iparm)
        sch%iparm(2)  = 3
        sch%iparm(24) = 1
        call TmS%set(7)
        call pardiso(sch%pt, 1, 1, 11, 11, neq, sch%csr_val, sch%csr_row, &
               sch%csr_col, sch%perm, 1, sch%iparm, msglvl, krn%rhs_or_solution, &
               krn%solution, error) 
        if (error /= 0) then
            print*, 'Pardiso analysis failed, error=', error
            call set_error('Pardiso analysis failed')
            return
        end if
        sch%analysed = .true.
        call TmS%set(8,7,'in Pardiso analysis')
    endif

    call TmS%set(4,11,'Matrix filling')
    call pardiso(sch%pt, 1, 1, 11, 23, neq, sch%csr_val, sch%csr_row, &
            sch%csr_col, sch%perm, 1, sch%iparm, msglvl, krn%rhs_or_solution, &
            krn%solution, error)
    if (error /= 0) then
        print*, 'Pardiso factorization failed, error=', error
        call set_error('Pardiso factorization failed')
        return
    end if
    do i=1,size(krn%rhs_or_solution)
        krn%rhs_or_solution(i)=krn%solution(i)
    enddo
    if (     for_explicit) call TmS%set(5,4,'in explicit Pardiso')
    if (.not.for_explicit) call TmS%set(6,4,'in implicit Pardiso')
    
end subroutine Pardiso_fact_solve         
      
end module krn_linear_system_m