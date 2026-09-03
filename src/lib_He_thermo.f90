! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later
module lib_He_thermo_m
    !! Helium equation of state (Arp reference EOS) and transport properties.
    use krn_global_tools_m, only: dp
    use lib_ext_math_m,     only: brent_t, zero
    use krn_simulation_m,   only: set_error
    use ieee_arithmetic,    only: ieee_value, ieee_quiet_nan
    implicit none
    private

    ! Index definitions (TODO should be in different library)
    integer, parameter, public :: Pri_ro = 1  !! Density
    integer, parameter, public :: Pri_u  = 2  !! Velocity
    integer, parameter, public :: Pri_p  = 3  !! Pressure
    integer, parameter, public :: Pri_e  = 4  !! Internal energy
    integer, parameter, public :: Pri_T  = 5  !! Temperature
    integer, parameter, public :: Pri_c  = 6  !! Sound speed
    integer, parameter, public :: Pri_R  = 7  !! Extra variable for the mechanical equilibrium recovering

    integer, parameter, public :: Con_Mas = 1 !! Mass conservative variable
    integer, parameter, public :: Con_Qdm = 2 !! Momentum conservative variable
    integer, parameter, public :: Con_Ene = 3 !! Energy conservative variable
    integer, parameter, public :: Con_R   = 4 !! Extra variable for conservative variable array

    !  Public interface functions
    public :: state_roT, state_roP, state_roP_withR, state_roE, state_roE_withR
    public :: r_roT, ro_pT, T_roP, T_roE, droeint_droP, jacobian_roT, dc2_roT, he_prop

    real(dp), parameter, public :: T_He_min = 3.0_dp                               !! EOS lower bound (K)
    real(dp), parameter, public :: T_He_max = 1500.0_dp                            !! EOS upper bound (K)
    real(dp), parameter, public :: MolMass = 4.0026e-3_dp                          !! Molar mass of helium (kg/mol)
    real(dp), parameter, public :: R_cte_gaz_Arp = 8.31431e-3_dp                   !! Universal gas constant (kJ/(mol*K))
    real(dp), parameter, public :: cp0_Arp = 5193.16943986d0 * MolMass * 1.0e-3_dp !! Reference specific heat at constant pressure (kJ/(kg*K))
    real(dp), parameter         :: tau = 1.0_dp / (17.399_dp**2)
    real(dp), parameter         :: q0_Arp = 5210.521090595d0 * MolMass * 1.0e-3_dp
    real(dp), parameter         :: nn(32) = [ &
        +0.4558980227431e-04_dp,  0.1260692007853e-02_dp, -0.7139657549318e-02_dp, &
        +0.9728903861441e-02_dp, -0.1589302471562e-01_dp,  0.1454229259623e-05_dp, &
        -0.4708238429298e-04_dp,  0.1132915223587e-02_dp,  0.2410763742104e-02_dp, &
        -0.5093547838381e-08_dp,  0.2699726927900e-05_dp, -0.3954146691114e-04_dp, &
        +0.1551961438127e-08_dp,  0.1050712335785e-07_dp, -0.5501158366750e-07_dp, &
        -0.1037673478521e-09_dp,  0.6446881346448e-12_dp,  0.3298960057071e-10_dp, &
        -0.3555585738784e-12_dp, -0.6885401367690e-02_dp,  0.9166109232806e-02_dp, &
        -0.6544314242937e-05_dp, -0.3315398880031e-04_dp, -0.2067693644676e-07_dp, &
        +0.3850153114958e-07_dp, -0.1399040626999e-10_dp, -0.1888462892389e-11_dp, &
        -0.4595138561035e-14_dp,  0.6872567403738e-14_dp, -0.6097223119177e-18_dp, &
        -0.7636186157005e-17_dp,  0.3848665703556e-17_dp &
    ]

    ! Newton/Brent solver tuning -- shared by every implicit solve in this module.
    real(dp), parameter :: Newton_T_guess0 = 10.0_dp
    real(dp), parameter :: Newton_FunctTol = 1.0e-8_dp
    integer,  parameter :: Newton_MaxIte   = 80
    real(dp), parameter :: Brent_MachEps   = 1.0e-15_dp
    real(dp), parameter :: Brent_Tol       = 1.0e-10_dp
    real(dp), parameter :: Brent_ResidTol  = Newton_FunctTol
    real(dp), parameter :: ro_from_pT_min_Arp = 0.01_dp/(MolMass*1000.0_dp)  ! Search bracket for ro_pT
    real(dp), parameter :: ro_from_pT_max_Arp = 290.0_dp/(MolMass*1000.0_dp) ! mol/L, converted to Arp units.

    type ThArrays_t
        real(dp), dimension(14) :: f,g,Derf,Derg,Dergg,hh,Derhh,Der2f,Der2g,Der2gg,gg,ff
    end type ThArrays_t

    ! Brent-solver objective types
    type, extends(brent_t) :: ro_from_pT_obj_t
        real(dp) :: p, T
      contains
        procedure f => funct_ro_from_pT
    end type
    type, extends(brent_t) :: T_from_roP_obj_t
        real(dp) :: ro, p
      contains
        procedure f => funct_T_from_roP
    end type
    type, extends(brent_t) :: T_from_roE_obj_t
        real(dp) :: ro, e
      contains
        procedure f => funct_T_from_roE
    end type

contains

!  Brent-solver objectives (fallback path for the three implicit solves)
function funct_ro_from_pT(me,x)
    class(ro_from_pT_obj_t), intent(in) :: me
    real(dp),                intent(in) :: x   !! x is ro (Arp units)
    real(dp) :: funct_ro_from_pT
    type(ThArrays_t) :: ThAr

    call eos_pT_terms(x, me%T, ThAr)
    funct_ro_from_pT = me%p - x*R_cte_gaz_Arp*me%T - sum(ThAr%f*ThAr%g)
end function funct_ro_from_pT

function funct_T_from_roP(me,x)
    class(T_from_roP_obj_t), intent(in) :: me
    real(dp),                intent(in) :: x   !! x is T
    real(dp) :: funct_T_from_roP
    type(ThArrays_t) :: ThAr

    call eos_pT_terms(me%ro, x, ThAr)
    funct_T_from_roP = me%p - me%ro*R_cte_gaz_Arp*x - sum(ThAr%f*ThAr%g)
end function funct_T_from_roP

function funct_T_from_roE(me,x)
    class(T_from_roE_obj_t), intent(in) :: me
    real(dp),                intent(in) :: x   !! x is T
    real(dp) :: funct_T_from_roE
    type(ThArrays_t) :: ThAr

    call eos_e_terms(me%ro, x, ThAr)
    funct_T_from_roE = me%e - (cp0_Arp-R_cte_gaz_Arp)*x - q0_Arp - sum(ThAr%gg*ThAr%hh)/me%ro
end function funct_T_from_roE

!  Equation of states (EOS) term evaluators
elemental subroutine fill_f_terms(ro, ThAr)
    real(dp), intent(in)  :: ro
    type(ThArrays_t), intent(inout) :: ThAr

    ThAr%f(1)=ro**2
    ThAr%f(2)=ro**3
    ThAr%f(3)=ro**4
    ThAr%f(4)=ro**5
    ThAr%f(5)=ro**6
    ThAr%f(6)=ro**7
    ThAr%f(7)=ro**8
    ThAr%f(8)=ro**9
    ThAr%f(9) =(ro**3) *exp(-tau*(ro**2))
    ThAr%f(10)=(ro**5) *exp(-tau*(ro**2))
    ThAr%f(11)=(ro**7) *exp(-tau*(ro**2))
    ThAr%f(12)=(ro**9) *exp(-tau*(ro**2))
    ThAr%f(13)=(ro**11)*exp(-tau*(ro**2))
    ThAr%f(14)=(ro**13)*exp(-tau*(ro**2))
end subroutine fill_f_terms

elemental subroutine fill_ff_terms(ro, ThAr)
    real(dp), intent(in)  :: ro
    type(ThArrays_t), intent(inout) :: ThAr

    ThAr%ff(1) =ro
    ThAr%ff(2) =(ro**2)/2.0_dp
    ThAr%ff(3) =(ro**3)/3.0_dp
    ThAr%ff(4) =(ro**4)/4.0_dp
    ThAr%ff(5) =(ro**5)/5.0_dp
    ThAr%ff(6) =(ro**6)/6.0_dp
    ThAr%ff(7) =(ro**7)/7.0_dp
    ThAr%ff(8) =(ro**8)/8.0_dp
    ThAr%ff(9) =(1.0_dp-exp(-tau*(ro**2)))/(2.0_dp*tau)
    ThAr%ff(10)=(1.0_dp-exp(-tau*(ro**2))*(tau*(ro**2)+1.0_dp))/(2.0_dp*(tau**2))
    ThAr%ff(11)=(2.0_dp-exp(-tau*(ro**2))*(tau*(ro**2)*(tau*(ro**2)+2.0_dp)+2.0_dp))/(2.0_dp*(tau**3))
    ThAr%ff(12)=(6.0_dp-exp(-tau*(ro**2))*(tau*(ro**2)*(tau*(ro**2)*(tau*(ro**2)+3.0_dp)+6.0_dp)+6.0_dp))/&
                (2.0_dp*(tau**4))
    ThAr%ff(13)=(exp(-tau*(ro**2))*(-tau*(ro**2)*(tau*(ro**2)*(tau*(ro**2)*(tau*(ro**2)+4.0_dp)+12.0_dp)&
                +24.0_dp)-24.0_dp)+24.0_dp)/(2.0_dp*(tau**5))
    ThAr%ff(14)=(exp(-tau*(ro**2))*(-tau*(ro**2)*(tau*(ro**2)*(tau*(ro**2)*(tau*(ro**2)*(tau*(ro**2)+5.0_dp)+&
                20.0_dp)+60.0_dp)+120.0_dp)-120.0_dp)+120.0_dp)/(2.0_dp*(tau**6))
end subroutine fill_ff_terms

elemental subroutine fill_g_Dreg_terms(T, ThAr)
    real(dp), intent(in)  :: T
    type(ThArrays_t), intent(inout) :: ThAr

    ThAr%g(1) =nn(1) *T+nn(2)*sqrt(T)+nn(3)+nn(4)/T+nn(5)/(T**2)
    ThAr%g(2) =nn(6) *T+nn(7)+nn(8)/T+nn(9)/(T**2)
    ThAr%g(3) =nn(10)*T+nn(11)+nn(12)/T
    ThAr%g(4) =nn(13)
    ThAr%g(5) =nn(14)/T+nn(15)/(T**2)
    ThAr%g(6) =nn(16)/T
    ThAr%g(7) =nn(17)/T+nn(18)/(T**2)
    ThAr%g(8) =nn(19)/(T**2)
    ThAr%g(9) =nn(20)/(T**2)+nn(21)/(T**3)
    ThAr%g(10)=nn(22)/(T**2)+nn(23)/(T**4)
    ThAr%g(11)=nn(24)/(T**2)+nn(25)/(T**3)
    ThAr%g(12)=nn(26)/(T**2)+nn(27)/(T**4)
    ThAr%g(13)=nn(28)/(T**2)+nn(29)/(T**3)
    ThAr%g(14)=nn(30)/(T**2)+nn(31)/(T**3)+nn(32)/(T**4)

    ThAr%Derg(1) =nn(1)+nn(2)/(2.0_dp*sqrt(T))-nn(4)/(T**2)-2.0_dp*nn(5)/(T**3)
    ThAr%Derg(2) =nn(6)-nn(8)/(T**2)-2.0_dp*nn(9)/(T**3)
    ThAr%Derg(3) =nn(10)-nn(12)/(T**2)
    ThAr%Derg(4) =0.0_dp
    ThAr%Derg(5) =-nn(14)/(T**2)-2.0_dp*nn(15)/(T**3)
    ThAr%Derg(6) =-nn(16)/(T**2)
    ThAr%Derg(7) =-nn(17)/(T**2)-2.0_dp*nn(18)/(T**3)
    ThAr%Derg(8) =-2.0_dp*nn(19)/(T**3)
    ThAr%Derg(9) =-2.0_dp*nn(20)/(T**3)-3.0_dp*nn(21)/(T**4)
    ThAr%Derg(10)=-2.0_dp*nn(22)/(T**3)-4.0_dp*nn(23)/(T**5)
    ThAr%Derg(11)=-2.0_dp*nn(24)/(T**3)-3.0_dp*nn(25)/(T**4)
    ThAr%Derg(12)=-2.0_dp*nn(26)/(T**3)-4.0_dp*nn(27)/(T**5)
    ThAr%Derg(13)=-2.0_dp*nn(28)/(T**3)-3.0_dp*nn(29)/(T**4)
    ThAr%Derg(14)=-2.0_dp*nn(30)/(T**3)-3.0_dp*nn(31)/(T**4)-4.0_dp*nn(32)/(T**5)
end subroutine fill_g_Dreg_terms

elemental subroutine eos_pT_terms(ro, T, ThAr) ! ThArr can be smaller: f,g,Derg
    real(dp), intent(in)  :: ro, T
    type(ThArrays_t), intent(inout) :: ThAr
    call fill_f_terms(ro, ThAr)
    call fill_g_Dreg_terms(T, ThAr)
end subroutine eos_pT_terms

elemental subroutine fill_e_Tpart(T, ThAr)
    real(dp), intent(in)  :: T
    type(ThArrays_t), intent(inout) :: ThAr

    ThAr%Der2g(1) =-nn(2)/(4.0_dp*(T**(3.0_dp/2.0_dp)))+2.0_dp*nn(4)/(T**3)+6.0_dp*nn(5)/(T**4)
    ThAr%Der2g(2) =2.0_dp*nn(8)/(T**3)+6.0_dp*nn(9)/(T**4)
    ThAr%Der2g(3) =2.0_dp*nn(12)/(T**3)
    ThAr%Der2g(4) =0.0_dp
    ThAr%Der2g(5) =2.0_dp*nn(14)/(T**3)+6.0_dp*nn(15)/(T**4)
    ThAr%Der2g(6) =2.0_dp*nn(16)/(T**3)
    ThAr%Der2g(7) =2.0_dp*nn(17)/(T**3)+6.0_dp*nn(18)/(T**4)
    ThAr%Der2g(8) =6.0_dp*nn(19)/(T**4)
    ThAr%Der2g(9) =6.0_dp*nn(20)/(T**4)+12.0_dp*nn(21)/(T**5)
    ThAr%Der2g(10)=6.0_dp*nn(22)/(T**4)+20.0_dp*nn(23)/(T**6)
    ThAr%Der2g(11)=6.0_dp*nn(24)/(T**4)+12.0_dp*nn(25)/(T**5)
    ThAr%Der2g(12)=6.0_dp*nn(26)/(T**4)+20.0_dp*nn(27)/(T**6)
    ThAr%Der2g(13)=6.0_dp*nn(28)/(T**4)+12.0_dp*nn(29)/(T**5)
    ThAr%Der2g(14)=6.0_dp*nn(30)/(T**4)+12.0_dp*nn(31)/(T**5)+20.0_dp*nn(32)/(T**6)

    call fill_g_Dreg_terms(T, ThAr)
    ThAr%gg    = ThAr%g - T*ThAr%Derg
    ThAr%Dergg = -T * ThAr%Der2g
end subroutine fill_e_Tpart

elemental subroutine eos_e_terms(ro, T, ThAr)
    real(dp), intent(in)  :: ro, T
    type(ThArrays_t), intent(inout) :: ThAr

    call fill_ff_terms(ro, ThAr)
    ThAr%hh = ro * ThAr%ff
    call fill_e_Tpart(T, ThAr)
end subroutine eos_e_terms

elemental subroutine eos_terms(ro, T, ThAr)
    real(dp), intent(in)  :: ro,T
    type(ThArrays_t), intent(out) :: ThAr
    real(dp) :: Der3g(14)

    ThAr%Derf(1) =2.0_dp*ro
    ThAr%Derf(2) =3.0_dp*(ro**2)
    ThAr%Derf(3) =4.0_dp*(ro**3)
    ThAr%Derf(4) =5.0_dp*(ro**4)
    ThAr%Derf(5) =6.0_dp*(ro**5)
    ThAr%Derf(6) =7.0_dp*(ro**6)
    ThAr%Derf(7) =8.0_dp*(ro**7)
    ThAr%Derf(8) =9.0_dp*(ro**8)
    ThAr%Derf(9) =exp(-tau*(ro**2))*( 3.0_dp*(ro**2 )-2.0_dp*tau*(ro**4 ))
    ThAr%Derf(10)=exp(-tau*(ro**2))*( 5.0_dp*(ro**4 )-2.0_dp*tau*(ro**6 ))
    ThAr%Derf(11)=exp(-tau*(ro**2))*( 7.0_dp*(ro**6 )-2.0_dp*tau*(ro**8 ))
    ThAr%Derf(12)=exp(-tau*(ro**2))*( 9.0_dp*(ro**8 )-2.0_dp*tau*(ro**10))
    ThAr%Derf(13)=exp(-tau*(ro**2))*(11.0_dp*(ro**10)-2.0_dp*tau*(ro**12))
    ThAr%Derf(14)=exp(-tau*(ro**2))*(13.0_dp*(ro**12)-2.0_dp*tau*(ro**14))

    ThAr%Der2f(1) =2.0_dp
    ThAr%Der2f(2) =3.0_dp*(2.0_dp*ro)
    ThAr%Der2f(3) =4.0_dp*(3.0_dp*ro**2)
    ThAr%Der2f(4) =5.0_dp*(4.0_dp*ro**3)
    ThAr%Der2f(5) =6.0_dp*(5.0_dp*ro**4)
    ThAr%Der2f(6) =7.0_dp*(6.0_dp*ro**5)
    ThAr%Der2f(7) =8.0_dp*(7.0_dp*ro**6)
    ThAr%Der2f(8) =9.0_dp*(8.0_dp*ro**7)
    ThAr%Der2f(9) =2.0_dp*ro      *exp(-tau*(ro**2))*(2.0_dp*(tau**2)*(ro**4)- 7.0_dp*tau*(ro**2)+ 3.0_dp)
    ThAr%Der2f(10)=2.0_dp*(ro**3 )*exp(-tau*(ro**2))*(2.0_dp*(tau**2)*(ro**4)-11.0_dp*tau*(ro**2)+10.0_dp)
    ThAr%Der2f(11)=2.0_dp*(ro**5 )*exp(-tau*(ro**2))*(2.0_dp*(tau**2)*(ro**4)-15.0_dp*tau*(ro**2)+21.0_dp)
    ThAr%Der2f(12)=2.0_dp*(ro**7 )*exp(-tau*(ro**2))*(2.0_dp*(tau**2)*(ro**4)-19.0_dp*tau*(ro**2)+36.0_dp)
    ThAr%Der2f(13)=2.0_dp*(ro**9 )*exp(-tau*(ro**2))*(2.0_dp*(tau**2)*(ro**4)-23.0_dp*tau*(ro**2)+55.0_dp)
    ThAr%Der2f(14)=2.0_dp*(ro**11)*exp(-tau*(ro**2))*(2.0_dp*(tau**2)*(ro**4)-27.0_dp*tau*(ro**2)+78.0_dp)

    Der3g(1) = (3.0_dp*(nn(2)*(T**(5.0_dp/2.0_dp))-16.0_dp*nn(4)*T-64.0_dp*nn(5)))/(8.0_dp*(T**5))
    Der3g(2) = -6.0_dp*nn( 8)/(T**4) -24.0_dp*nn(9)/(T**5)
    Der3g(3) = -6.0_dp*nn(12)/(T**4)
    Der3g(4) =  0.0_dp
    Der3g(5) = -6.0_dp*nn(14)/(T**4) -24.0_dp*nn(15)/(T**5)
    Der3g(6) = -6.0_dp*nn(16)/(T**4)
    Der3g(7) = -6.0_dp*nn(17)/(T**4) -24.0_dp*nn(18)/(T**5)
    Der3g(8) =-24.0_dp*nn(19)/(T**5)
    Der3g(9) =-24.0_dp*nn(20)/(T**5) -60.0_dp*nn(21)/(T**6)
    Der3g(10)=-24.0_dp*nn(22)/(T**5)-120.0_dp*nn(23)/(T**7)
    Der3g(11)=-24.0_dp*nn(24)/(T**5) -60.0_dp*nn(25)/(T**6)
    Der3g(12)=-24.0_dp*nn(26)/(T**5)-120.0_dp*nn(27)/(T**7)
    Der3g(13)=-24.0_dp*nn(28)/(T**5) -60.0_dp*nn(29)/(T**6)
    Der3g(14)=-24.0_dp*nn(30)/(T**5) -60.0_dp*nn(31)/(T**6)-120.0_dp*nn(32)/(T**7)

    call fill_f_terms(ro, ThAr)
    call eos_e_terms(ro, T, ThAr)
    ThAr%Derhh = ThAr%ff + ThAr%f/ro
    ThAr%Der2gg = -ThAr%Der2g - T*Der3g
end subroutine eos_terms


function brent(obj, xmin, xmax, label) result(x)
    class(brent_t), intent(in) :: obj
    real(dp),       intent(in) :: xmin, xmax
    character(*),   intent(in) :: label
    real(dp) :: x
    real(dp) :: fmin, fmax

    fmin = obj%f(xmin)
    fmax = obj%f(xmax)
    if (fmin*fmax > 0.0_dp) then
        call set_error(label//': no sign change on bracket, cannot solve')
        x = ieee_value(x, ieee_quiet_nan)
        return
    endif

    x = zero(obj, xmin, xmax, Brent_MachEps, Brent_Tol)
    if (abs(obj%f(x)) > Brent_ResidTol .or. x < xmin .or. x > xmax) then
        call set_error(label//': Brent solve failed to converge')
        x = ieee_value(x, ieee_quiet_nan)
    endif
end function brent


! =====================================================================
!  Public elemental utilities (SI units in, SI units out)
! =====================================================================

elemental function r_roT(ro, T) result(r)
    !! r = rho*h consistent with (ro,T), fast direct evaluation
    real(dp), intent(in) :: ro, T
    real(dp) :: r
    type(ThArrays_t) :: ThAr
    real(dp) :: ro_Arp
    integer  :: n

    ro_Arp = ro/(MolMass*1000.0_dp)
    call fill_f_terms(ro_Arp, ThAr)
    call fill_ff_terms(ro_Arp, ThAr)
    call fill_g_Dreg_terms(T, ThAr)
    ThAr%hh = ro_Arp * ThAr%ff
    ThAr%gg = ThAr%g - T*ThAr%Derg

    r = (cp0_Arp*T + q0_Arp)*ro_Arp
    do n = 1, 14
        r = r + ThAr%gg(n)*ThAr%hh(n) + ThAr%f(n)*ThAr%g(n)
    enddo
    r = r * 1.0e6_dp
end function r_roT


function ro_pT(p, T) result(ro)
    !! Solve p(ro,T) = p for ro.
    real(dp), intent(in) :: p, T
    real(dp) :: ro
    type(ro_from_pT_obj_t) :: obj

    obj%T = T
    obj%p = p*1.0e-6_dp
    ro = brent(obj, ro_from_pT_min_Arp, ro_from_pT_max_Arp, 'ro_pT')
    ro = ro*MolMass*1000.0_dp
end function ro_pT


function T_roP(ro, p) result(T)
    !! Solve p(ro,T) = p for T. 
    real(dp), intent(in) :: ro, p
    real(dp) :: T

    type(T_from_roP_obj_t) :: obj
    type(ThArrays_t) :: ThAr
    real(dp) :: ro_Arp, p_Arp
    real(dp) :: Funct, DerFunct
    logical  :: converged
    integer  :: ite

    ro_Arp = ro/(MolMass*1000.0_dp)
    p_Arp  = p*1.0e-6_dp
    obj%ro = ro_Arp
    obj%p  = p_Arp

    call fill_f_terms(ro_Arp, ThAr)

    T = Newton_T_guess0
    converged = .false.
    do ite = 1, Newton_MaxIte
        call fill_g_Dreg_terms(T, ThAr)
        Funct = p_Arp - ro_Arp*R_cte_gaz_Arp*T - sum(ThAr%f*ThAr%g)
        if (abs(Funct) <= Newton_FunctTol) then
            converged = .true.
            exit
        endif
        DerFunct = -ro_Arp*R_cte_gaz_Arp - sum(ThAr%f*ThAr%Derg)
        T = T - Funct/DerFunct
    enddo
    if (converged .and. T >= T_He_min .and. T <= T_He_max) return
    T = brent(obj, T_He_min, T_He_max, 'T_roP')
end function T_roP


function T_roE(ro, e) result(T)
    !! Solve e(ro,T) = e for T
    real(dp), intent(in) :: ro, e
    real(dp) :: T

    type(T_from_roE_obj_t) :: obj
    type(ThArrays_t) :: ThAr
    real(dp) :: ro_Arp, e_Arp
    real(dp) :: Funct, DerFunct
    logical  :: converged
    integer  :: ite

    ro_Arp = ro/(MolMass*1000.0_dp)
    e_Arp  = e*MolMass*1.0e-3_dp
    obj%ro = ro_Arp
    obj%e  = e_Arp

    call fill_ff_terms(ro_Arp, ThAr)
    ThAr%hh = ro_Arp * ThAr%ff

    T = Newton_T_guess0
    converged = .false.
    do ite = 1, Newton_MaxIte
        call fill_e_Tpart(T, ThAr)
        Funct = e_Arp - (cp0_Arp-R_cte_gaz_Arp)*T - q0_Arp - sum(ThAr%gg*ThAr%hh)/ro_Arp
        if (abs(Funct) <= Newton_FunctTol) then
            converged = .true.
            exit
        endif
        DerFunct = -(cp0_Arp-R_cte_gaz_Arp) - sum(ThAr%Dergg*ThAr%hh)/ro_Arp
        T = T - Funct/DerFunct
    enddo
    if (converged .and. T >= T_He_min .and. T <= T_He_max) return
    T = brent(obj, T_He_min, T_He_max, 'T_roE')
end function T_roE


function droeint_droP(ro, p, e, dedv, dedp) result(droeint_dro)
    !! Energy-related terms from (ro,p):
    !! returns d(rho*e)/d(ro)|p and can also provide e, de/dv|p, and de/dp|v.
    real(dp), intent(in) :: ro, p
    real(dp), intent(out), optional :: e, dedv, dedp
    real(dp) :: droeint_dro
    real(dp) :: ro_Arp, T, cvloc_Arp, dPdRo_Arp, dPdT_Arp, eint_Arp, deintdRo_Arp, dedT_Arp, e_loc
    type(ThArrays_t) :: ThAr
    integer :: n

    T = T_roP(ro, p)
    ro_Arp = ro/(MolMass*1000.0_dp)
    call eos_terms(ro_Arp, T, ThAr)

    cvloc_Arp    = cp0_Arp - R_cte_gaz_Arp
    dPdRo_Arp    = R_cte_gaz_Arp*T
    dPdT_Arp     = ro_Arp*R_cte_gaz_Arp
    eint_Arp     = cp0_Arp*T + q0_Arp - R_cte_gaz_Arp*T
    deintdRo_Arp = 0.0_dp
    dedT_Arp     = cp0_Arp - R_cte_gaz_Arp
    do n = 1, 14
        cvloc_Arp    = cvloc_Arp    + ThAr%Dergg(n)*ThAr%hh(n)/ro_Arp
        dPdRo_Arp    = dPdRo_Arp    + ThAr%Derf(n)*ThAr%g(n)
        dPdT_Arp     = dPdT_Arp     + ThAr%f(n)*ThAr%Derg(n)
        eint_Arp     = eint_Arp     + ThAr%gg(n)*ThAr%ff(n)
        deintdRo_Arp = deintdRo_Arp + ThAr%gg(n)*ThAr%f(n)/(ro_Arp**2)
        dedT_Arp     = dedT_Arp     + ThAr%Dergg(n)*ThAr%ff(n)
    enddo

    e_loc = eint_Arp * 1000.0_dp/MolMass
    droeint_dro = ro_Arp*deintdRo_Arp - ro_Arp*cvloc_Arp*dPdRo_Arp/dPdT_Arp + eint_Arp
    droeint_dro = droeint_dro * 1000.0_dp/MolMass
    if (present(e))    e    = e_loc
    if (present(dedv)) dedv = ro*(e_loc - droeint_dro)
    if (present(dedp)) dedp = (dedT_Arp*1000.0_dp/MolMass)/(dPdT_Arp*1.0e6_dp)
end function droeint_droP


! =====================================================================
!  Public state functions -- one per physically distinct call shape
! =====================================================================

elemental subroutine state_roT(ro, T, r, p, e, c)
    !! Full state from (ro,T): r, p, e, c.
    real(dp), intent(in)  :: ro, T
    real(dp), intent(out) :: r, p, e, c
    real(dp) :: dPdRo, dPdT, dTdRo, cv

    call jacobian_roT(ro, T, cv, dPdT, dTdRo=dTdRo, r=r, p=p)
    dPdRo = -dTdRo*dPdT
    c = sqrt(dPdRo + T*(dPdT**2)/(cv*(ro**2)))
    e = (r - p)/ro
end subroutine state_roT


subroutine state_roP(ro, p, r, e, T, c)
    !! Full state from (ro,p): r, e, T, c. Solves for T internally.
    real(dp), intent(in)  :: ro, p
    real(dp), intent(out) :: r, e, T, c
    real(dp) :: dPdRo, dPdT, dTdRo, cv

    T = T_roP(ro, p)
    call jacobian_roT(ro, T, cv, dPdT, dTdRo=dTdRo, r=r)
    dPdRo = -dTdRo*dPdT
    c = sqrt(dPdRo + T*(dPdT**2)/(cv*(ro**2)))
    e = (r - p)/ro
end subroutine state_roP


subroutine state_roP_withR(ro, p, r, e, T, c)
    !! State from (ro,p) when r is already known (e.g. carried as a
    !! conservative variable): e, T, c. Solves for T internally.
    real(dp), intent(in)  :: ro, p, r
    real(dp), intent(out) :: e, T, c
    real(dp) :: dPdRo, dPdT, dTdRo, cv

    T = T_roP(ro, p)
    call jacobian_roT(ro, T, cv, dPdT, dTdRo=dTdRo)
    dPdRo = -dTdRo*dPdT
    c = sqrt(dPdRo + T*(dPdT**2)/(cv*(ro**2)))
    e = (r - p)/ro
end subroutine state_roP_withR

subroutine state_roE(ro, e, p, T, c)
    !! State from (ro,e): p, T, c for the physical EOS path.
    real(dp), intent(in)  :: ro, e
    real(dp), intent(out) :: p, T, c
    real(dp) :: dPdRo, dPdT, dTdRo, cv

    T = T_roE(ro, e)
    call jacobian_roT(ro, T, cv, dPdT, dTdRo=dTdRo, p=p)
    dPdRo = -dTdRo*dPdT
    c = sqrt(dPdRo + T*(dPdT**2)/(cv*(ro**2)))
end subroutine state_roE

subroutine state_roE_withR(ro, e, r, p, T, c)
    !! State from (ro,e) when r is already known: p, T, c.
    real(dp), intent(in)  :: ro, e, r
    real(dp), intent(out) :: p, T, c
    real(dp) :: dPdRo, dPdT, dTdRo, cv

    p = r - ro*e
    T = T_roP(ro, p)
    call jacobian_roT(ro, T, cv, dPdT, dTdRo=dTdRo)
    dPdRo = -dTdRo*dPdT
    c = sqrt(dPdRo + T*(dPdT**2)/(cv*(ro**2)))
end subroutine state_roE_withR


elemental subroutine jacobian_roT(ro, T, cv, dPdT, dTdP, dTdRo, dRdRo, dRdT, r, p, d2TdP_dT, d2TdP_dRo)
    real(dp), intent(in)  :: ro, T
    real(dp), intent(out) :: cv, dPdT
    real(dp), intent(out), optional :: dTdP, dTdRo, dRdRo, dRdT, r, p
    real(dp), intent(out), optional :: d2TdP_dT, d2TdP_dRo

    type(ThArrays_t) :: ThAr
    real(dp) :: ro_Arp, dPdRo_Arp, dPdT_Arp, cv_Arp, p_Arp, r_Arp, dRdRo_Arp, dRdT_Arp
    real(dp) :: d2PdT2_Arp, d2PdRodT_Arp
    integer  :: n

    ro_Arp = ro/(MolMass*1000.0_dp)
    call eos_terms(ro_Arp, T, ThAr)

    dPdRo_Arp    = R_cte_gaz_Arp*T
    dPdT_Arp     = ro_Arp*R_cte_gaz_Arp
    cv_Arp       = cp0_Arp - R_cte_gaz_Arp
    p_Arp        = ro_Arp*R_cte_gaz_Arp*T
    r_Arp        = (cp0_Arp*T + q0_Arp)*ro_Arp
    dRdRo_Arp    = cp0_Arp*T + q0_Arp
    dRdT_Arp     = cp0_Arp*ro_Arp
    d2PdT2_Arp   = 0.0_dp
    d2PdRodT_Arp = R_cte_gaz_Arp
    do n = 1, 14
        dPdRo_Arp    = dPdRo_Arp    + ThAr%Derf(n) *ThAr%g(n)
        dPdT_Arp     = dPdT_Arp     + ThAr%f(n)    *ThAr%Derg(n)
        cv_Arp       = cv_Arp       + ThAr%Dergg(n)*ThAr%hh(n) / ro_Arp
        p_Arp        = p_Arp        + ThAr%f(n)    *ThAr%g(n)
        r_Arp        = r_Arp        + ThAr%gg(n)   *ThAr%hh(n)    + ThAr%f(n)   *ThAr%g(n)
        dRdRo_Arp    = dRdRo_Arp    + ThAr%gg(n)   *ThAr%Derhh(n) + ThAr%Derf(n)*ThAr%g(n)
        dRdT_Arp     = dRdT_Arp     + ThAr%Dergg(n)*ThAr%hh(n)    + ThAr%f(n)   *ThAr%Derg(n)
        d2PdT2_Arp   = d2PdT2_Arp   + ThAr%f(n)    *ThAr%Der2g(n)
        d2PdRodT_Arp = d2PdRodT_Arp + ThAr%Derf(n) *ThAr%Derg(n)
    enddo

    dPdT  = dPdT_Arp  * 1.0e6_dp
    cv    = cv_Arp    * 1.0e3_dp/MolMass

    if (present(dTdP))      dTdP      = 1.0_dp/dPdT
    if (present(dTdRo))     dTdRo     = -(dPdRo_Arp * 1000.0_dp/MolMass)/dPdT
    if (present(dRdRo))     dRdRo     = dRdRo_Arp * 1.0e6_dp/(MolMass*1000.0_dp)
    if (present(dRdT))      dRdT      = dRdT_Arp  * 1.0e6_dp
    if (present(r))         r         = r_Arp     * 1.0e6_dp
    if (present(p))         p         = p_Arp     * 1.0e6_dp

    if (present(d2TdP_dT))  d2TdP_dT  = -(d2PdT2_Arp   * 1.0e6_dp)                    /(dPdT**2)
    if (present(d2TdP_dRo)) d2TdP_dRo = -(d2PdRodT_Arp * 1.0e6_dp/(MolMass*1000.0_dp))/(dPdT**2)
end subroutine jacobian_roT


elemental subroutine dc2_roT(ro, T, dc2dRo, dc2dP)
    !! d(c^2)/dRo and d(c^2)/dP at constant T / constant Ro respectively.
    real(dp), intent(in)  :: ro, T
    real(dp), intent(out) :: dc2dRo, dc2dP

    real(dp) :: dcSQdRo_T, dcSQdT_Ro, dTdRo, dTdP, AA, BB, DD, cvloc_Arp
    real(dp) :: ro_Arp, DerFunctdRo, DerFunctdP, DerFunctdT
    real(dp) :: DerAAdRo, DerBBdRo, DerDDdRo, DerAAdT, DerBBdT, DerDDdT
    real(dp) :: DercvlocdRo, DercvlocdT
    type(ThArrays_t) :: ThAr
    integer :: n

    ro_Arp = ro/(MolMass*1000.0_dp)
    call eos_terms(ro_Arp, T, ThAr)

    DerFunctdRo = -R_cte_gaz_Arp*T
    DerFunctdP  = 1.0_dp
    DerFunctdT  = -ro_Arp*R_cte_gaz_Arp
    AA = R_cte_gaz_Arp*T
    cvloc_Arp = cp0_Arp - R_cte_gaz_Arp
    DD = ro_Arp*R_cte_gaz_Arp
    DerAAdRo = 0.0_dp
    DerDDdRo = R_cte_gaz_Arp
    DerAAdT  = R_cte_gaz_Arp
    DerDDdT  = 0.0_dp
    DercvlocdRo = 0.0_dp
    DercvlocdT  = 0.0_dp
    do n = 1, 14
        DerFunctdRo = DerFunctdRo - ThAr%Derf(n)*ThAr%g(n)
        DerFunctdT  = DerFunctdT  - ThAr%f(n)*ThAr%Derg(n)
        AA = AA + ThAr%Derf(n)*ThAr%g(n)
        cvloc_Arp = cvloc_Arp + ThAr%Dergg(n)*ThAr%hh(n)/ro_Arp
        DD = DD + ThAr%f(n)*ThAr%Derg(n)
        DerAAdRo = DerAAdRo + ThAr%Der2f(n)*ThAr%g(n)
        DerDDdRo = DerDDdRo + ThAr%Derf(n)*ThAr%Derg(n)
        DerAAdT  = DerAAdT  + ThAr%Derf(n)*ThAr%Derg(n)
        DerDDdT  = DerDDdT  + ThAr%f(n)*ThAr%Der2g(n)
        DercvlocdRo = DercvlocdRo + ThAr%Dergg(n)*(ThAr%Derhh(n)*ro_Arp-ThAr%hh(n))/(ro_Arp**2)
        DercvlocdT  = DercvlocdT  + ThAr%Der2gg(n)*ThAr%hh(n)/ro_Arp
    enddo
    BB = T/(cvloc_Arp*ro_Arp*ro_Arp)
    DerBBdRo = -T*(DercvlocdRo*ro_Arp*ro_Arp+2.0_dp*cvloc_Arp*ro_Arp)/((cvloc_Arp*ro_Arp*ro_Arp)**2)
    DerBBdT  = (1.0_dp/(ro_Arp**2))*(cvloc_Arp-T*DercvlocdT)/(cvloc_Arp**2)
    dTdP  = -DerFunctdP/DerFunctdT
    dTdRo = -DerFunctdRo/DerFunctdT

    dcSQdRo_T = DerAAdRo + DerBBdRo*(DD**2) + 2.0_dp*BB*DD*DerDDdRo
    dcSQdT_Ro = DerAAdT  + DerBBdT*(DD**2)  + 2.0_dp*BB*DD*DerDDdT

    dc2dRo = dcSQdRo_T + dcSQdT_Ro*dTdRo
    dc2dP  = dcSQdT_Ro*dTdP

    dc2dRo = 1000.0_dp*dc2dRo/MolMass/(MolMass*1000.0_dp)
    dc2dP  = 1000.0_dp*dc2dP/MolMass*1.0e-6_dp
end subroutine dc2_roT


subroutine he_prop(ro, T, cv, cp, mu, lambda, dPdT) ! TODO not elemental due to set_error 
    !! Viscosity and thermal conductivity correlations for helium. The critical-region
    !! enhancement term in lambda needs dP/dRo and dP/dT;
    real(dp), intent(in)  :: ro, T
    real(dp), intent(out) :: cv, cp, mu, lambda, dPdT
    real(dp) :: dPdRo, dTdRo

    real(dp), parameter :: CC(4) = [3.739232544_dp, -2.620316969e1_dp, 5.982252246e1_dp, -4.926397634e1_dp]
    real(dp), parameter :: Coef(11) = [1.862970530e-4_dp, -7.275964435e-7_dp, -1.427549651e-4_dp, &
                    3.290833592e-5_dp, -5.213335363e-8_dp, 4.492659933e-8_dp, -5.924416513e-9_dp, &
                    7.087321137e-6_dp, -6.013335678e-6_dp, 8.067145814e-7_dp,  3.995125013e-7_dp]
    real(dp) :: dRodP
    real(dp) :: ro_unit
    real(dp) :: x, EtaPrim0, BBB, CCC, DDD, EtaPrimE, Eta0, EtaE
    real(dp) :: Lambd0, LambdC, LambdE, SumCjSurTj
    real(dp) :: AA, BB, x0, E1, E2, Beta, Gama, Delta, roc, Tc, Pc, RR, m, k
    real(dp) :: DeltaT, DeltaRo, Eta, KT, rostar, LongTerm, KTprim, CoefCrit
    real(dp) :: W, hhh, dhhhdx
    integer  :: j

    call jacobian_roT(ro, T, cv, dPdT, dTdRo=dTdRo)
    dPdRo = -dTdRo*dPdT
    cp = cv + T*(dPdT**2)/(dPdRo*(ro**2))

    if (cp < 0.0_dp .or. cv < 0.0_dp) then
        call set_error('negative cp or cv in he_prop')
        return
    endif

    dRodP = 1.0_dp/dPdRo
    ro_unit = ro/1000.0_dp
    if (T <= 300.0_dp) then
        x = log(T)
    else
        x = 5.7037825_dp  ! log(300.0_dp)
    endif
    EtaPrim0 = -0.135311743_dp/x+1.00347841_dp+1.20654649_dp*x-0.149564551_dp*(x**2)+0.0125208416_dp*(x**3)
    BBB = -47.5295259_dp/x+87.6799309_dp-42.0741589_dp*x+8.33128289_dp*(x**2)-0.589252385_dp*(x**3)
    CCC = 547.309267_dp/x-904.870586_dp+431.404928_dp*x-81.4504854_dp*(x**2)+5.37008433_dp*(x**3)
    DDD = -1684.39324_dp/x+3331.08630_dp-1632.19172_dp*x+308.804413_dp*(x**2)-20.2936367_dp*(x**3)
    EtaPrimE = ro_unit*BBB+(ro_unit**2)*CCC+(ro_unit**3)*DDD
    if (EtaPrimE > 100.0_dp) then
        call set_error('EtaPrimE > 100 during viscosity calculation in he_prop')
        return
    endif
    if (T <= 100.0_dp) then
        mu = exp(EtaPrim0+EtaPrimE)
    else
        Eta0 = 196.0_dp*(T**0.71938_dp)*exp(12.451_dp/T-295.67_dp/(T**2)-4.1249_dp)
        if (T < 110.0_dp) then
            Eta0 = exp(EtaPrim0)+(Eta0-exp(EtaPrim0))*(T-100.0_dp)/10.0_dp
        endif
        EtaE = exp(EtaPrim0+EtaPrimE)-exp(EtaPrim0)
        mu = Eta0+EtaE
    endif
    mu = 1.0e-7_dp*mu

    AA = 2.7870034e-3_dp
    BB = 7.034007057e-1_dp
    SumCjSurTj = 0.0_dp
    do j = 1, 4
        SumCjSurTj = SumCjSurTj + CC(j)/(T**j)
    enddo
    Lambd0 = AA*(T**BB)*exp(SumCjSurTj)

    LambdC = 0.0_dp
    if (T >= 3.5_dp .and. T <= 12.0_dp) then
        x0 = 0.392_dp
        E1 = 2.8461_dp
        E2 = 0.27156_dp
        Beta = 0.3554_dp
        Gama = 1.1743_dp
        Delta = 4.304_dp
        roc = 69.158_dp
        Tc = 5.18992_dp
        Pc = 2.2746e5_dp
        RR = 4.633e-10_dp
        m = 6.6455255e-27_dp
        k = 1.38066e-23_dp
        DeltaT = abs(1.0_dp-T/Tc)
        DeltaRo = abs(1.0_dp-ro/roc)
        Eta = mu
        rostar = ro/roc

        ! sqrt(m*k)/(6*Pi*RR), as in REFPROP:
        CoefCrit = 3.4685233e-17_dp

        KT = dRodP/ro
        W = (DeltaT/0.2_dp)**2+(DeltaRo/0.25_dp)**2

        if (W < 1.0_dp) then
            x = DeltaT/(DeltaRo**(1.0_dp/Beta))
            hhh = E1*(1.0_dp+x/x0)*(1.0_dp+E2*(1.0_dp+x/x0)**(2.0_dp*Beta))**((Gama-1.0_dp)/(2.0_dp*Beta))

            dhhhdx = E1*(Gama-1.0_dp)*E2*(1.0_dp+x/x0)**(2.0_dp*Beta)*(1.0_dp+ &
                E2*(1.0_dp+x/x0)**(2.0_dp*Beta))**((Gama-1.0_dp)/(2.0_dp*Beta)-1.0_dp)+E1*(1.0_dp+ &
                E2*(1.0_dp+x/x0)**(2.0_dp*Beta))**((Gama-1.0_dp)/(2.0_dp*Beta))
            dhhhdx = dhhhdx/x0

            LongTerm = (DeltaRo**(Delta-1.0_dp))*(Delta*hhh-x*dhhhdx/Beta)
            KTprim = 1.0_dp/(rostar**2)/Pc/LongTerm
            KT = W*KT+(1.0_dp-W)*KTprim
        endif

        LambdC = CoefCrit*sqrt(KT*(T**3)/ro)*(dPdT**2)*exp(-18.66_dp*(DeltaT**2)-4.25_dp*(DeltaRo**4))/Eta
    endif

    rostar = ro/68.0_dp
    LambdE=(Coef(1)+Coef(2)*T+Coef(3)*T**(1.0_dp/3.0_dp)+Coef(4)*T**(2.0_dp/3.0_dp))*ro+&
           (Coef(5)+Coef(6)*T**(1.0_dp/3.0_dp)+Coef(7)*T**(2.0_dp/3.0_dp))*(ro**3)+&
           (Coef(8)+Coef(9)*T**(1.0_dp/3.0_dp)+Coef(10)*T**(2.0_dp/3.0_dp)+Coef(11)/T)*(ro**2)*log(rostar)

    lambda = Lambd0 + LambdC + LambdE
end subroutine he_prop

end module lib_He_thermo_m
