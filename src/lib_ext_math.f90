!! GNU LGPL licensed file
module lib_ext_math_m
    !! Mathematical functions downloaded from John Burkardt page fSolve and brent
    !! 
    !!  - https://people.sc.fsu.edu/~jburkardt/f_src/fsolve/fsolve.html
    !!  - https://people.math.sc.edu/Burkardt/f_src/brent/brent.html
    !!
    !! Modified:
    !!
    !!  - 06 April 2010
    !!  - 01 February 2024 - introducing object with virtual fcn
    !!
    !! Author:
    !!
    !!  - Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !!  - FORTRAN90 version by John Burkardt.
    !!  - Fortran03 version by Jacek Kosek and Damien Furfaro 
    !!    
    !! Reference:
    !!
    !!  - Jorge More, Burton Garbow, Kenneth Hillstrom,
    !!  - User Guide for MINPACK-1,
    !!  - Technical Report ANL-80-74,
    !!  - Argonne National Laboratory, 1980.
    use krn_global_tools_m, only: dp
    implicit none
    private
    public zero, brent_t, fSolve, fSolve_t
    type, abstract :: fSolve_t
        !! Procedure 'fcn' of this class is the parameter to fSolve subroutine.
        !!
        !! This abstract object should be extended by user the only requirement
        !! is to overwrite abstract 'fcn' class procedure, which will be used by
        !! fSolve algorithm.
    contains
        procedure (fSolve_function), deferred :: fcn
    end type
    abstract interface
        subroutine fSolve_function(me,n,x,fVec)  !! User-supplied system of equations
            import :: fSolve_t,dp
            class(fSolve_t), intent(in) :: me  !! User-defined object transmitted to the solver 
            integer,  intent(in)    :: n       !! number of functions and variables
            real(dp), intent(inout) :: x(n)    !! result or initial point where the jacobian is evaluated
            real(dp), intent(inout) :: fVec(n) !! the functions evaluated at 'x'
        end subroutine fSolve_function
    end interface

    type, abstract :: brent_t !! class with the function f(x)
    contains
        procedure (brent_function), deferred :: f
    end type
    abstract interface
        function brent_function(me,x) 
            !! user-supplied function, of the form which evaluates the function whose zero is being sought.
            import :: brent_t,dp
            class(brent_t), intent(in) :: me !! User-defined object transmitted to the solver 
            real(dp),       intent(in) :: x  !! 'x' parameter of the function f(x)
            real(dp) :: brent_function       !! returned value
        end function brent_function
    end interface

contains

subroutine dogleg(n, r, lr, diag, qtb, delta, x)
    !! Finds the minimizing combination of Gauss-Newton and gradient steps.
    !!
    !! Given an M by N matrix A, an N by N nonsingular diagonal
    !! matrix D, an M-vector B, and a positive number DELTA, the
    !! problem is to determine the convex combination X of the
    !! Gauss-Newton and scaled gradient directions that minimizes
    !! (A*X - B) in the least squares sense, subject to the
    !! restriction that the euclidean norm of D*X be at most DELTA.
    !! 
    !! This function completes the solution of the problem
    !! if it is provided with the necessary information from the
    !! QR factorization of A.  That is, if A = Q*R, where Q has
    !! orthogonal columns and R is an upper triangular matrix,
    !! then DOGLEG expects the full upper triangle of R and
    !! the first N components of Q'*B.
    integer  :: lr     !! In - size of the R array, which must be no less than (N*(N+1))/2.
    integer  :: n      !! In - order of the matrix R
    real(dp) :: delta  !! In - positive upper bound on the euclidean norm of D*X(1:N).
    real(dp) :: diag(n)!! In - diagonal elements of the matrix D.
    real(dp) :: qtb(n) !! In - first N elements of the vector Q'* B.
    real(dp) :: r(lr)  !! In - upper triangular matrix R stored by rows.
    real(dp) :: x(n)   !! Out - desired convex combination of the Gauss-Newton 
                       !!       direction and the scaled gradient direction.
    real(dp) :: alpha
    real(dp) :: bNorm
    real(dp) :: epsMch
    real(dp) :: gNorm
    integer  :: i
    integer  :: j
    integer  :: jj
    integer  :: k
    integer  :: l
    real(dp) :: qNorm
    real(dp) :: sgNorm
    real(dp) :: sum2
    real(dp) :: temp
    real(dp) :: wa1(n)
    real(dp) :: wa2(n)

    epsMch = epsilon ( epsMch )
    ! Calculate the Gauss-Newton direction.
    jj = ( n * ( n + 1 ) ) / 2 + 1
    do k = 1, n
        j = n - k + 1
        jj = jj - k
        l = jj + 1
        sum2 = 0.0D+00
        do i = j + 1, n
            sum2 = sum2 + r(l) * x(i)
            l = l + 1
        end do
        temp = r(jj)
        if ( temp == 0.0D+00 ) then
            l = j
            do i = 1, j
                temp = max ( temp, abs ( r(l)) )
                l = l + n - i
            end do
            if ( temp == 0.0D+00 ) then
                temp = epsMch
            else
                temp = epsMch * temp
            end if
        end if
        x(j) = ( qtb(j) - sum2 ) / temp
    end do

    ! Test whether the Gauss-Newton direction is acceptable.
    wa1(1:n) = 0.0D+00
    wa2(1:n) = diag(1:n) * x(1:n)
    qNorm = eNorm ( n, wa2 )
    if ( qNorm <= delta ) then
        return
    end if
    
    ! The Gauss-Newton direction is not acceptable.
    ! Calculate the scaled gradient direction.
    l = 1
    do j = 1, n
        temp = qtb(j)
        do i = j, n
            wa1(i) = wa1(i) + r(l) * temp
            l = l + 1
        enddo
        wa1(j) = wa1(j) / diag(j)
    enddo

    ! Calculate the norm of the scaled gradient.
    ! Test for the special case in which the scaled gradient is zero.
    gNorm = eNorm ( n, wa1 )
    sgNorm = 0.0D+00
    alpha = delta / qNorm
    
    if ( gNorm /= 0.0D+00 ) then
        ! Calculate the point along the scaled gradient which minimizes the quadratic.
        wa1(1:n) = ( wa1(1:n) / gNorm ) / diag(1:n)
        l = 1
        do j = 1, n
            sum2 = 0.0D+00
            do i = j, n
                sum2 = sum2 + r(l) * wa1(i)
                l = l + 1
            enddo
            wa2(j) = sum2
        enddo
    
        temp = eNorm ( n, wa2 )
        sgNorm = ( gNorm / temp ) / temp

        ! Test whether the scaled gradient direction is acceptable.
        alpha = 0.0D+00
    
        ! The scaled gradient direction is not acceptable.
        ! Calculate the point along the dogleg at which the quadratic is minimized.
        if ( sgNorm < delta ) then
            bNorm = eNorm ( n, qtb )
            temp = ( bNorm / gNorm ) * ( bNorm / qNorm ) * ( sgNorm / delta )
            temp = temp - ( delta / qNorm ) * ( sgNorm / delta) ** 2 &
                + sqrt ( ( temp - ( delta / qNorm ) ) ** 2 &
                + ( 1.0D+00 - ( delta / qNorm ) ** 2 ) &
                * ( 1.0D+00 - ( sgNorm / delta ) ** 2 ) )
    
            alpha = ( ( delta / qNorm ) * ( 1.0D+00 - ( sgNorm / delta ) ** 2 ) ) &
                / temp
        endif
    endif
    
    ! Form appropriate convex combination of the Gauss-Newton
    ! direction and the scaled gradient direction.
    temp = ( 1.0D+00 - alpha ) * min ( sgNorm, delta )
    x(1:n) = temp * wa1(1:n) + alpha * x(1:n)
end

function eNorm(n, x)
    !! Computes the Euclidean norm of a vector.
    !!
    !! The Euclidean norm is computed by accumulating the sum of
    !! squares in three different sums.  The sums of squares for the
    !! small and large components are scaled so that no overflows
    !! occur.  Non-destructive underflows are permitted.  Underflows
    !! and overflows do not occur in the computation of the unscaled
    !! sum of squares for the intermediate components.
    !! 
    !! The definitions of small, intermediate and large components
    !! depend on two constants, RDWARF and RGIANT.  The main
    !! restrictions on these constants are that RDWARF^2 not
    !! underflow and RGIANT^2 not overflow.
    integer  :: n     !! Length of the vector
    real(dp) :: x(n)  !! Vector whose norm is desired
    real(dp) :: eNorm !! Euclidean norm of the vector

    integer  :: i
    real(dp) :: agiant
    real(dp) :: rdwarf
    real(dp) :: rgiant
    real(dp) :: s1
    real(dp) :: s2
    real(dp) :: s3
    real(dp) :: xAbs
    real(dp) :: x1max
    real(dp) :: x3max
  
    rdwarf = sqrt ( tiny ( rdwarf ) )
    rgiant = sqrt ( huge ( rgiant ) )
  
    s1 = 0.0D+00
    s2 = 0.0D+00
    s3 = 0.0D+00
    x1max = 0.0D+00
    x3max = 0.0D+00
    agiant = rgiant / real ( n, kind = 8 )
  
    do i = 1, n
        xAbs = abs ( x(i) )
        if ( xAbs <= rdwarf ) then
            if ( x3max < xAbs ) then
                s3 = 1.0D+00 + s3 * ( x3max / xAbs ) ** 2
                x3max = xAbs
            else if ( xAbs /= 0.0D+00 ) then
                s3 = s3 + ( xAbs / x3max ) ** 2
            end if
        else if ( agiant <= xAbs ) then
            if ( x1max < xAbs ) then
                s1 = 1.0D+00 + s1 * ( x1max / xAbs ) ** 2
                x1max = xAbs
            else
                s1 = s1 + ( xAbs / x1max ) ** 2
            end if
        else
            s2 = s2 + xAbs ** 2
        end if
    end do
    
    !  Calculation of norm.
    if ( s1 /= 0.0D+00 ) then
        eNorm = x1max * sqrt ( s1 + ( s2 / x1max ) / x1max )
    else if ( s2 /= 0.0D+00 ) then
        if ( x3max <= s2 ) then
            eNorm = sqrt ( s2 * ( 1.0D+00 + ( x3max / s2 ) * ( x3max * s3 ) ) )
        else
            eNorm = sqrt ( x3max * ( ( s2 / x3max ) + ( x3max * s3 ) ) )
        end if
    else
        eNorm = x3max * sqrt ( s3 )
    end if
end

subroutine fdJac1(obj, n, x, fVec, fjac, ldfjac, ml, mu, epsFcn)
    !! Estimates a jacobian matrix using forward differences.
    !!
    !! This function computes a forward-difference approximation
    !! to the N by N jacobian matrix associated with a specified
    !! problem of N functions in N variables. If the jacobian has
    !! a banded form, then function evaluations are saved by only
    !! approximating the nonzero terms.
    class(fSolve_t), intent(in) :: obj !! External object with the user-supplied 
                                       !! subroutine which calculates the functions.
    integer,  intent(in)    :: n       !! Number of functions and variables
    real(dp), intent(inout) :: x(n)    !! Point where the jacobian is evaluated
    real(dp), intent(in)    :: fVec(n) !! Functions evaluated at 'x'
    real(dp), intent(inout) :: fjac(ldfjac,n) !! The N by N approximate
    integer,  intent(in)    :: ldfjac  !! Leading dimension of FJAC, which must not be less than N.
    integer,  intent(in)    ::  ml     !! Number of sub-diagonals within the band of the jacobian matrix
                                       !! If the jacobian is not banded, set ML and MU to N-1.
    integer,  intent(in)    ::  mu     !! Number of super-diagonals within the band of the jacobian matrix
                                       !! If the jacobian is not banded, set ML and MU to N-1.
    real(dp), intent(in)    :: epsFcn
        !! Is used in determining a suitable step length for the forward-difference approximation.
        !! 
        !! This approximation assumes that the relative errors in the functions are of the order of epsFcn.
        !! If epsFcn is less than the machine precision, it is assumed that 
        !! the relative errors in the functions are of the order of the machine precision.

    real(dp) :: eps
    real(dp) :: epsMch
    real(dp) :: h
    integer ::  i
    integer ::  j
    integer ::  k
    integer ::  mSum
    real(dp) :: temp
    real(dp) :: wa1(n)
    real(dp) :: wa2(n)

    epsMch = epsilon ( epsMch )
    eps = sqrt ( max ( epsFcn, epsMch ) )
    mSum = ml + mu + 1

    ! Computation of dense approximate jacobian.
    if ( n <= mSum ) then
        do j = 1, n
            temp = x(j)
            h = eps * abs ( temp )
            if ( h == 0.0D+00 ) then
                h = eps
            end if
            x(j) = temp + h
            call obj%fcn(n, x, wa1)
            x(j) = temp
            fjac(1:n,j) = ( wa1(1:n) - fVec(1:n) ) / h
        end do
    else
        ! Computation of banded approximate jacobian.
        do k = 1, mSum
            do j = k, n, mSum
                wa2(j) = x(j)
                h = eps * abs ( wa2(j) )
                if ( h == 0.0D+00 ) then
                    h = eps
                end if
                x(j) = wa2(j) + h
            end do  
            call obj%fcn(n, x, wa1)  
            do j = k, n, mSum 
                x(j) = wa2(j)  
                h = eps * abs ( wa2(j) )
                if ( h == 0.0D+00 ) then
                    h = eps
                end if  
                fjac(1:n,j) = 0.0D+00 
                do i = 1, n
                    if ( j - mu <= i .and. i <= j + ml ) then
                        fjac(i,j) = ( wa1(i) - fVec(i) ) / h
                    end if
                end do 
            end do 
        end do 
    end if
end

subroutine fsolve(obj, n, x, fVec, tol, info)
    !! Seeks a zero of N non-linear equations in N variables
    !!
    !! Finds a zero of a system of N non-linear functions in N variables
    !! by a modification of the Powell hybrid method. This is done by using the
    !! more general non-linear equation solver HYBRD. The user provides a
    !! subroutine which calculates the functions.  
    !! The jacobian is calculated by a forward-difference approximation.
    class(fSolve_t), intent(in) :: obj !! External object with the user-supplied 
                                       !!   subroutine which calculates the functions.
    integer,  intent(in)    :: n       !! Number of functions and variables
    real(dp), intent(inout) :: x(n)    !! Initial estimate of the solution vector /
                                       !!   the estimate of the solution vector
    real(dp), intent(out)   :: fVec(n) !! The functions evaluated at the output X
    real(dp), intent(in)    :: tol     !! Satisfactory termination occurs when the algorithm
                                       !!   estimates that the relative error between X and
                                       !!   the solution is at most TOL. TOL should be non-negative.
    integer,  intent(out)   :: info    !! Status flag
        !!
        !!  0. Improper input parameters.
        !!  1. Algorithm estimates that the relative error between X and the solution is at most TOL.
        !!  2. Number of calls to FCN has reached or exceeded 200*(N+1).
        !!  3. TOL is too small.  No further improvement in the approximate solution X is possible.
        !!  4. Iteration is not making good progress, as measured by the improvement
        !!     from the last five jacobian evaluations.
        !!  5. Iteration is not making good progress, as measured by the improvement
        !!     from the last ten iterations.
    real(dp) :: diag(n)
    real(dp) :: epsFcn
    real(dp) :: factor
    real(dp) :: fjac(n,n)
    integer ::  ldfjac
    integer ::  lr
    integer ::  maxFev
    integer ::  ml
    integer ::  mode
    integer ::  mu
    integer ::  nFev
    real(dp) :: qtf(n)
    real(dp) :: r((n*(n+1))/2)
    real(dp) :: xTol
  
    if ( n <= 0 ) then
        info = 0
        return
    end if
  
    if ( tol < 0.0D+00 ) then
        info = 0
        return
    end if
  
    xTol = tol
    maxFev = 200 * ( n + 1 )
    ml = n - 1
    mu = n - 1
    epsFcn = 0.0D+00
    diag(1:n) = 1.0D+00
    mode = 2
    factor = 100.0D+00
    info = 0
    nFev = 0
    fjac(1:n,1:n) = 0.0D+00
    ldfjac = n
    r(1:(n*(n+1))/2) = 0.0D+00
    lr = ( n * ( n + 1 ) ) / 2
    qtf(1:n) = 0.0D+00
  
    call hybrd (obj, n, x, fVec, xTol, maxFev, ml, mu, epsFcn, diag, mode, &
        factor, info, nFev, fjac, ldfjac, r, lr, qtf)
end

subroutine hybrd(obj, n, x, fVec, xTol, maxFev, ml, mu, epsFcn, diag, mode, &
    factor, info, nFev, fjac, ldfjac, r, lr, qtf)
    !! Seeks a zero of N non-linear equations in N variables.
    !!
    !! HYBRD finds a zero of a system of N non-linear functions in N variables
    !! by a modification of the Powell hybrid method.  The user must provide a
    !! subroutine which calculates the functions.  
    !!
    !! The jacobian is then calculated by a forward-difference approximation.
    class(fSolve_t), intent(in) :: obj !! External object with the user-supplied 
                                       !!   subroutine which calculates the functions.
    integer,  intent(in)    :: n       !! The number of functions and variables
    real(dp), intent(inout) :: x(n)    !! Input initial estimate of the solution vector /
                                       !!   output X final estimate of the solution vector
    real(dp), intent(out)   :: fVec(n) !! The functions evaluated at the output X
    real(dp), intent(in)    :: xTol    !! Termination occurs when the relative error
            !! between two consecutive iterates is at most XTOL. XTOL should be non-negative.
    integer,  intent(in)    :: maxFev  !! Termination occurs when the number of calls to
                                       !!   FCN is at least MAXFEV by the end of an iteration
    integer,  intent(in)    :: ml      !! Specify the number of sub-diagonals within the band
            !! of the jacobian matrix. If the jacobian is not banded, set ML to at least n - 1.
    integer,  intent(in)    :: mu      !! Specify the number of super-diagonals within the band
            !! of the jacobian matrix. If the jacobian is not banded, set MU to at least n - 1.
    real(dp), intent(in)    :: epsFcn  !! is used in determining a suitable step length for the
            !! forward-difference approximation.  This approximation assumes that the relative
            !! errors in the functions are of the order of EPSFCN.  If EPSFCN is less than the
            !! machine precision, it is assumed that the relative errors in the functions are
            !! of the order of the machine precision.
    real(dp), intent(inout) :: diag(n) !! If MODE = 1, then DIAG is set internally.
            !! If MODE = 2, then DIAG must contain positive entries that
            !! serve as multiplicative scale factors for the variables.
    integer,  intent(in)    :: mode    !! variables will be scaled internally.
                                       !!   2, scaling is specified by the input DIAG vector.
    real(dp), intent(in)    :: factor  !! determines the initial step bound. This bound is set
            !! to the product of FACTOR and the euclidean norm of DIAG*X if non-zero, or else to
            !! FACTOR itself. In most cases, FACTOR should lie in the interval (0.1, 100) with
            !! 100 the recommended value.
    integer,  intent(out)   :: info    !! Error flag
            !!
            !!  0. Improper input parameters.
            !!  1. Relative error between two consecutive iterates is at most XTOL.
            !!  2. Number of calls to FCN has reached or exceeded MAXFEV.
            !!  3. XTOL is too small.  No further improvement in the approximate
            !!     solution X is possible.
            !!  4. Iteration is not making good progress, as measured by the improvement
            !!     from the last five jacobian evaluations.
            !!  5. Iteration is not making good progress, as measured by the improvement
            !!     from the last ten iterations.
    integer,  intent(out) ::  nFev   !! The number of calls to FCN.
                                     !!    produced by the QR factorization of the final
                                     !!    approximate jacobian.
    integer,  intent(in) ::  ldfjac  !! The leading dimension of FJAC. LDFJAC must be at least N.
    integer,  intent(in) ::  lr      !! Size of the R array, which must be no less than (N*(N+1))/2.
    real(dp), intent(out) :: r(lr)   !! Upper triangular matrix produced by the QR factorization
                                     !!    of the final approximate jacobian, stored row-wise.
    real(dp), intent(out) :: qtf(n)  !! Contains the vector Q'*FVEC.
    real(dp), intent(out) :: fjac(ldfjac,n) !! Out - N by N array which contains the orthogonal matrix Q

    real(dp) :: actred
    real(dp) :: delta
    real(dp) :: epsMch
    real(dp) :: fNorm
    real(dp) :: fNorm1
    integer  :: i
    integer  :: iter
    integer  :: iwa(1)
    integer  :: j
    logical  :: jEval
    integer  :: l
    integer  :: mSum
    integer  :: ncFail
    integer  :: nsLow1
    integer  :: nsLow2
    integer  :: ncSuc
    logical  :: pivot
    real(dp) :: pNorm
    real(dp) :: prered
    real(dp) :: ratio
    logical  :: sing
    real(dp) :: sum2
    real(dp) :: temp
    real(dp) :: wa1(n)
    real(dp) :: wa2(n)
    real(dp) :: wa3(n)
    real(dp) :: wa4(n)
    real(dp) :: xNorm
    
    epsMch = epsilon(epsMch)
    info = 0
    nFev = 0

    ! Check the input parameters for errors.
    if ( n <= 0 ) then
        return
    else if ( xTol < 0.0D+00 ) then
        return
    else if ( maxFev <= 0 ) then
        return
    else if ( ml < 0 ) then
        return
    else if ( mu < 0 ) then
        return
    else if ( factor <= 0.0D+00 ) then
        return
    else if ( ldfjac < n ) then
        return
    else if ( lr < ( n * ( n + 1 ) ) / 2 ) then
        return
    end if
  
    if ( mode == 2 ) then
        do j = 1, n
            if ( diag(j) <= 0.0D+00 ) then
                return
            end if
        end do
    end if

    ! Evaluate the function at the starting point
    ! and calculate its norm.
    call obj%fcn(n, x, fVec)
    nFev = 1
    
    fNorm = eNorm(n, fVec)

    ! Determine the number of calls to FCN needed to compute the jacobian matrix.
    mSum = min ( ml + mu + 1, n )

    ! Initialize iteration counter and monitors.
    iter = 1
    ncSuc = 0
    ncFail = 0
    nsLow1 = 0
    nsLow2 = 0
    
    ! Beginning of the outer loop.
    30 continue
        jEval = .true.
    
        ! Calculate the jacobian matrix.
        call fdJac1 (obj, n, x, fVec, fjac, ldfjac, ml, mu, epsFcn )
        nFev = nFev + mSum

        ! Compute the QR factorization of the jacobian.
        pivot = .false.
        call qrFac ( n, n, fjac, ldfjac, pivot, iwa, 1, wa1, wa2 )

        ! On the first iteration, if MODE is 1, scale according
        ! to the norms of the columns of the initial jacobian.
        if ( iter == 1 ) then
            if ( mode /= 2 ) then
                diag(1:n) = wa2(1:n)
                do j = 1, n
                    if ( wa2(j) == 0.0D+00 ) then
                        diag(j) = 1.0D+00
                    end if
                end do
            end if

            ! On the first iteration, calculate the norm of the scaled X
            ! and initialize the step bound DELTA.
            wa3(1:n) = diag(1:n) * x(1:n)
            xNorm = eNorm ( n, wa3 )
            delta = factor * xNorm
            if ( delta == 0.0D+00 ) then
                delta = factor
            end if
        end if

        ! Form Q' * FVEC and store in QTF.
        qtf(1:n) = fVec(1:n)
        do j = 1, n
            if ( fjac(j,j) /= 0.0D+00 ) then
                temp = - dot_product ( qtf(j:n), fjac(j:n,j) ) / fjac(j,j)
                qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
           end if
        end do

        ! Copy the triangular factor of the QR factorization into R.
        sing = .false.
    
        do j = 1, n
            l = j
            do i = 1, j - 1
                r(l) = fjac(i,j)
                l = l + n - i
            end do
            r(l) = wa1(j)
            if ( wa1(j) == 0.0D+00 ) then
                sing = .true.
            end if
        end do

        ! Accumulate the orthogonal factor in FJAC.
        call qform(n, n, fjac, ldfjac)
    
        !  Rescale if necessary.
        if ( mode /= 2 ) then
            do j = 1, n
                diag(j) = max ( diag(j), wa2(j) )
            end do
        end if
    
    ! Beginning of the inner loop.
    180 continue
    
        !  Determine the direction P.
        call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
    
        ! Store the direction P and X + P.
        ! Calculate the norm of P.
        wa1(1:n) = - wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)
    
        pNorm = eNorm ( n, wa3 )
    
        ! On the first iteration, adjust the initial step bound.
        if ( iter == 1 ) then
            delta = min ( delta, pNorm )
        end if
    
        ! Evaluate the function at X + P and calculate its norm.
        call obj%fcn(n, wa2, wa4)
        nFev = nFev + 1
        fNorm1 = eNorm ( n, wa4 )
    
        ! Compute the scaled actual reduction.
        actred = -1.0D+00
        if ( fNorm1 < fNorm ) then
            actred = 1.0D+00 - ( fNorm1 / fNorm ) ** 2
        endif
    
        ! Compute the scaled predicted reduction.
        l = 1
        do i = 1, n
            sum2 = 0.0D+00
            do j = i, n
                sum2 = sum2 + r(l) * wa1(j)
                l = l + 1
            end do
            wa3(i) = qtf(i) + sum2
        end do
    
        temp = eNorm ( n, wa3 )
        prered = 0.0D+00
        if ( temp < fNorm ) then
            prered = 1.0D+00 - ( temp / fNorm ) ** 2
        end if
    
        ! Compute the ratio of the actual to the predicted reduction.
        ratio = 0.0D+00
        if ( 0.0D+00 < prered ) then
            ratio = actred / prered
        end if
    
        ! Update the step bound.
        if ( ratio < 0.1D+00 ) then
    
            ncSuc = 0
            ncFail = ncFail + 1
            delta = 0.5D+00 * delta
        else
            ncFail = 0
            ncSuc = ncSuc + 1
            if ( 0.5D+00 <= ratio .or. 1 < ncSuc ) then
                delta = max ( delta, pNorm / 0.5D+00 )
            end if
            if ( abs ( ratio - 1.0D+00 ) <= 0.1D+00 ) then
                delta = pNorm / 0.5D+00
            end if
        end if
    
        ! Test for successful iteration.
    
        ! Successful iteration.
        ! Update X, FVEC, and their norms.
        if ( 0.0001D+00 <= ratio ) then
            x(1:n) = wa2(1:n)
            wa2(1:n) = diag(1:n) * x(1:n)
            fVec(1:n) = wa4(1:n)
            xNorm = eNorm ( n, wa2 )
            fNorm = fNorm1
            iter = iter + 1
        end if
    
        ! Determine the progress of the iteration.
        nsLow1 = nsLow1 + 1
        if ( 0.001D+00 <= actred ) then
            nsLow1 = 0
        end if
        if ( jEval ) then
            nsLow2 = nsLow2 + 1
        end if
        if ( 0.1D+00 <= actred ) then
            nsLow2 = 0
        end if

        ! Test for convergence.
        if ( delta <= xTol * xNorm .or. fNorm == 0.0D+00 ) then
            info = 1
        end if
        if ( info /= 0 ) then
            return
        end if

        ! Tests for termination and stringent tolerances.
        if ( maxFev <= nFev ) then
            info = 2
        end if
        if ( 0.1D+00 * max ( 0.1D+00 * delta, pNorm ) <= epsMch * xNorm ) then
            info = 3
        end if
        if ( nsLow2 == 5 ) then
            info = 4
        end if
        if ( nsLow1 == 10 ) then
            info = 5
        end if
        if ( info /= 0 ) then
            return
        end if
    
        ! Criterion for recalculating jacobian approximation
        ! by forward differences.
        if ( ncFail == 2 ) then
            go to 290
        end if
    
        ! Calculate the rank one modification to the jacobian
        ! and update QTF if necessary.
        do j = 1, n
            sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
            wa2(j) = ( sum2 - wa3(j) ) / pNorm
            wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pNorm )
            if ( 0.0001D+00 <= ratio ) then
                qtf(j) = sum2
            end if
        end do
    
        ! Compute the QR factorization of the updated jacobian.
        call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
        call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
        call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
    
        ! End of the inner loop.
        jEval = .false.
        go to 180
    290 continue
        ! End of the outer loop.
        go to 30
end

subroutine qform(m, n, q, ldq) 
    !! Produces the explicit QR factorization of a matrix.
    !!
    !! The QR factorization of a matrix is usually accumulated in implicit
    !! form, that is, as a series of orthogonal transformations of the
    !! original matrix.  This routine carries out those transformations,
    !! to explicitly exhibit the factorization constructed by QRFAC.
    integer ::  ldq      !! In - not less than M which specifies the leading dimension of the array Q
    integer ::  m        !! In - number of rows of A and the order of Q
    integer ::  n        !! In - number of columns of A
    real(dp) :: q(ldq,m) !! In/Out - lower trapezoid in the first min(M,N) / accumulated into a square matrix
        !!
        !! On input the full lower trapezoid in the first min(M,N) columns of Q contains the factored form.
        !! On output, Q has been accumulated into a square matrix.
  
    integer ::  j
    integer ::  k
    integer ::  l
    integer ::  minmn
    real(dp) :: temp
    real(dp) :: wa(m)
  
    minmn = min ( m, n )
    do j = 2, minmn
        q(1:j-1,j) = 0.0D+00
    end do
  
    ! Initialize remaining columns to those of the identity matrix.
    q(1:m,n+1:m) = 0.0D+00
    do j = n + 1, m
        q(j,j) = 1.0D+00
    end do
  
    ! Accumulate Q from its factored form.
    do l = 1, minmn
        k = minmn - l + 1
        wa(k:m) = q(k:m,k)
        q(k:m,k) = 0.0D+00
        q(k,k) = 1.0D+00
        if ( wa(k) /= 0.0D+00 ) then
            do j = k, m
                temp = dot_product ( wa(k:m), q(k:m,j) ) / wa(k)
                q(k:m,j) = q(k:m,j) - temp * wa(k:m)
            end do
        end if
    end do
end

subroutine qrFac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm)
    !! computes a QR factorization using Householder transformations.
    !!
    !! This function uses Householder transformations with optional column
    !! pivoting to compute a QR factorization of the
    !! M by N matrix A.  That is, QRFAC determines an orthogonal
    !! matrix Q, a permutation matrix P, and an upper trapezoidal
    !! matrix R with diagonal elements of nonincreasing magnitude,
    !! such that A*P = Q*R.  
    !!
    !! The Householder transformation for column K, K = 1,2,...,min(M,N), 
    !! is of the form:
    !!
    !! I - ( 1 / U(K) ) * U * U'
    !!
    !! where U has zeros in the first K-1 positions.  
    !!
    !! The form of this transformation and the method of pivoting first
    !! appeared in the corresponding LINPACK routine.
    integer  :: lda         !! In - leading dimension of A, which must be no less than M.
    integer  :: lipvt       !! In - dimension of IPVT, which should be N if pivoting is used.
    integer  :: m           !! In - number of rows of A.
    integer  :: n           !! In - number of columns of A.
    logical  :: pivot       !! Out - is TRUE if column pivoting is to be carried out
    real(dp) :: rdiag(n)    !! Out - contains the diagonal elements of R.
    integer ::  ipvt(lipvt) !! Out - defines the permutation matrix P such that A*P = Q*R.
        !!
        !! Column J of P is column IPVT(J) of the identity matrix.  If PIVOT is false,
        !! IPVT is not referenced.
    real(dp) :: ajnorm   !! Out - norms of the corresponding columns of the input matrix A.
        !!
        !! If this information is not needed, then ACNORM can coincide with RDIAG.
    real(dp) :: a(lda,n) !! In/Out - before/after QR factosiation size should be the M by N array.
        !!
        !! On input, A contains the matrix for which the QR factorization is to
        !! be computed.  On output, the strict upper trapezoidal part of A contains
        !! the strict upper trapezoidal part of R, and the lower trapezoidal
        !! part of A contains a factored form of Q, the non-trivial elements of
        !! the U vectors described above.

    real(dp) :: acnorm(n)
    real(dp) :: epsmch
    integer ::  i4_temp
    integer ::  j
    integer ::  k
    integer ::  kmax
    integer ::  minmn
    real(dp) :: r8_temp(m)
    real(dp) :: temp
    real(dp) :: wa(n)
  
    epsmch = epsilon ( epsmch )
    
    ! Compute the initial column norms and initialize several arrays.
    do j = 1, n
        acnorm(j) = eNorm ( m, a(1:m,j) )
    end do
    rdiag(1:n) = acnorm(1:n)
    wa(1:n) = acnorm(1:n)
    if ( pivot ) then
        do j = 1, n
            ipvt(j) = j
        end do
    end if

    ! Reduce A to R with Householder transformations.
    minmn = min ( m, n )
    do j = 1, minmn
        ! Bring the column of largest norm into the pivot position.
        if ( pivot ) then
            kmax = j
            do k = j, n
                if ( rdiag(kmax) < rdiag(k) ) then
                    kmax = k
                end if
            end do
            if ( kmax /= j ) then
                r8_temp(1:m) = a(1:m,j)
                a(1:m,j)     = a(1:m,kmax)
                a(1:m,kmax)  = r8_temp(1:m)      
                rdiag(kmax) = rdiag(j)
                wa(kmax)   = wa(j)
                i4_temp    = ipvt(j)
                ipvt(j)    = ipvt(kmax)
                ipvt(kmax) = i4_temp
            end if
        end if
    
        ! Compute the Householder transformation to reduce the
        ! J-th column of A to a multiple of the J-th unit vector.
        ajnorm = eNorm ( m-j+1, a(j,j) )
        if ( ajnorm /= 0.0D+00 ) then
            if ( a(j,j) < 0.0D+00 ) then
                ajnorm = -ajnorm
            end if
            a(j:m,j) = a(j:m,j) / ajnorm
            a(j,j) = a(j,j) + 1.0D+00
      
            ! Apply the transformation to the remaining columns and update the norms.
            do k = j + 1, n
                temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)
                a(j:m,k) = a(j:m,k) - temp * a(j:m,j)
                if ( pivot .and. rdiag(k) /= 0.0D+00 ) then
                    temp = a(j,k) / rdiag(k)
                    rdiag(k) = rdiag(k) * sqrt ( max ( 0.0D+00, 1.0D+00 - temp ** 2 ) )
                    if ( 0.05D+00 * ( rdiag(k) / wa(k) ) ** 2 <= epsmch ) then
                        rdiag(k) = eNorm ( m-j, a(j+1,k) )
                        wa(k) = rdiag(k)
                    end if
                end if
            end do
        end if
        rdiag(j) = - ajnorm
    end do
end


subroutine r1mpyq(m, n, a, lda, v, w)
    !! computes A*Q, where Q is the product of Householder transformations
    !!
    !! Given an M by N matrix A, this function computes A*Q where
    !! Q is the product of 2*(N - 1) transformations
    !!
    !!      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
    !!
    !! and GV(I), GW(I) are Givens rotations in the (I,N) plane which
    !! eliminate elements in the I-th and N-th planes, respectively.
    !! Q itself is not given, rather the information to recover the
    !! GV, GW rotations is supplied.
    integer ::  lda  !! In - leading dimension of A, which must not be less than M.
    integer ::  m    !! In - number of rows of A
    integer ::  n    !! In - number of columns of A
    real(dp) :: v(n) !! In - contain the information necessary to recover the Givens rotations GV
    real(dp) :: w(n) !! In - contain the information necessary to recover the Givens rotations GW
    real(dp) :: a(lda,n) !! In/Out - the M by N array. On input, the matrix A to be postmultiplied
                     !!              by the orthogonal matrix Q. On output, the value of A*Q.

    real(dp) :: c
    integer ::  i
    integer ::  j
    real(dp) :: s
    real(dp) :: temp
  
    !  Apply the first set of Givens rotations to A.
    do j = n - 1, 1, -1
        if ( 1.0D+00 < abs ( v(j) ) ) then
            c = 1.0D+00 / v(j)
            s = sqrt ( 1.0D+00 - c ** 2 )
        else
            s = v(j)
            c = sqrt ( 1.0D+00 - s ** 2 )
        end if
        do i = 1, m
            temp =   c * a(i,j) - s * a(i,n)
            a(i,n) = s * a(i,j) + c * a(i,n)
            a(i,j) = temp
        end do
    end do
  
    !  Apply the second set of Givens rotations to A.
    do j = 1, n - 1
        if ( 1.0D+00 < abs ( w(j) ) ) then
            c = 1.0D+00 / w(j)
            s = sqrt ( 1.0D+00 - c ** 2 )
        else
            s = w(j)
            c = sqrt ( 1.0D+00 - s ** 2 )
        end if
        do i = 1, m
            temp =     c * a(i,j) + s * a(i,n)
            a(i,n) = - s * a(i,j) + c * a(i,n)
            a(i,j) = temp
        end do
    end do
end


subroutine r1updt(m, n, s, ls, u, v, w, sing)
    !! re-triangularizes a matrix after a rank one update
    !!
    !! Given an M by N lower trapezoidal matrix S, an M-vector U, and an
    !! N-vector V, the problem is to determine an orthogonal matrix Q such that
    !! 
    !!   (S + U * V' ) * Q
    !! 
    !! is again lower trapezoidal.
    !! 
    !! This function determines Q as the product of 2 * (N - 1)
    !! transformations
    !! 
    !!   GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
    !! 
    !! where GV(I), GW(I) are Givens rotations in the (I,N) plane
    !! which eliminate elements in the I-th and N-th planes,
    !! respectively.  Q itself is not accumulated, rather the
    !! information to recover the GV and GW rotations is returned.
    integer ::  ls   !! In - length of the S array.  LS must be at least (N*(2*M-N+1))/2.
    integer ::  m    !! In - number of rows of S.
    integer ::  n    !! In - number of columns of S. N must not exceed M.  
    real(dp) :: s(ls)!! In/Out - On input, the lower trapezoidal matrix S stored by columns.
                     !!      On output S contains the lower trapezoidal
                     !!      matrix produced as described above.
    logical sing     !! Out - is set to TRUE if any of the diagonal elements
                     !!      of the output S are zero.  Otherwise SING is set FALSE.
    real(dp) :: u(m) !! In - the U vector.
    real(dp) :: v(n) !! In/Out - On input, V must contain the vector V.  On output V
                     !!      contains the information necessary to recover the
                     !!      Givens rotations GV described above.
    real(dp) :: w(m) !! Out - contains information necessary to recover the Givens
                     !!      rotations GW described above.
  
    real(dp) :: sin
    real(dp) :: tan
    real(dp) :: tau
    real(dp) :: temp
    real(dp) :: cos
    real(dp) :: cotan
    real(dp) :: giant
    integer ::  i
    integer ::  j
    integer ::  jj
    integer ::  l

    ! GIANT is the largest magnitude.
    giant = huge ( giant )
  
    ! Initialize the diagonal element pointer.
    jj = ( n * ( 2 * m - n + 1 ) ) / 2 - ( m - n )
  
    ! Move the nontrivial part of the last column of S into W.
    l = jj
    do i = n, m
        w(i) = s(l)
        l = l + 1
    end do
  
    ! Rotate the vector V into a multiple of the N-th unit vector
    ! in such a way that a spike is introduced into W.
    do j = n - 1, 1, -1
        jj = jj - ( m - j + 1 )
        w(j) = 0.0D+00
        if ( v(j) /= 0.0D+00 ) then
            ! Determine a Givens rotation which eliminates the J-th element of V.
            if ( abs ( v(n) ) < abs ( v(j) ) ) then
                cotan = v(n) / v(j)
                sin = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan ** 2 )
                cos = sin * cotan
                tau = 1.0D+00
                if ( abs ( cos ) * giant > 1.0D+00 ) then
                    tau = 1.0D+00 / cos
                end if
            else
                tan = v(j) / v(n)
                cos = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan ** 2 )
                sin = cos * tan
                tau = sin
            end if
    
            ! Apply the transformation to V and store the information
            ! necessary to recover the Givens rotation.
            v(n) = sin * v(j) + cos * v(n)
            v(j) = tau
      
            ! Apply the transformation to S and extend the spike in W.
            l = jj
            do i = j, m
                temp = cos * s(l) - sin * w(i)
                w(i) = sin * s(l) + cos * w(i)
                s(l) = temp
                l = l + 1
            enddo
        endif
    enddo
    
    ! Add the spike from the rank 1 update to W.
    w(1:m) = w(1:m) + v(n) * u(1:m)
    
    ! Eliminate the spike.
    sing = .false.
    do j = 1, n-1
        if ( w(j) /= 0.0D+00 ) then
            ! Determine a Givens rotation which eliminates the
            ! J-th element of the spike.
            if ( abs ( s(jj) ) < abs ( w(j) ) ) then
                cotan = s(jj) / w(j)
                sin = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan ** 2 )
                cos = sin * cotan
                if ( 1.0D+00 < abs ( cos ) * giant ) then
                    tau = 1.0D+00 / cos
                else
                    tau = 1.0D+00
                endif
            else
                tan = w(j) / s(jj)
                cos = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan ** 2 )
                sin = cos * tan
                tau = sin
            endif
    
            ! Apply the transformation to S and reduce the spike in W.
            l = jj
            do i = j, m
                temp = cos * s(l) + sin * w(i)
                w(i) = - sin * s(l) + cos * w(i)
                s(l) = temp
                l = l + 1
            end do
        
            ! Store the information necessary to recover the Givens rotation.
            w(j) = tau
        endif
    
        ! Test for zero diagonal elements in the output S.
        if ( s(jj) == 0.0D+00 ) then
            sing = .true.
        end if
        jj = jj + ( m - j + 1 )
    end do

    ! Move W back into the last column of the output S.
    l = jj
    do i = n, m
        s(l) = w(i)
        l = l + 1
    end do
    if ( s(jj) == 0.0D+00 ) then
        sing = .true.
    end if
end

function zero(obj, a, b, machep, t)
    !! Seeks the root of a function F(X) in an interval [A,B].
    !!
    !! The interval [A,B] must be a change of sign interval for F.
    !! That is, F(A) and F(B) must be of opposite signs.  Then
    !! assuming that F is continuous implies the existence of at least
    !! one value C between A and B for which F(C) = 0.
    !!
    !! The location of the zero is determined to within an accuracy
    !! of 6 * MACHEPS * abs ( C ) + 2 * T.
    !!
    !! Thanks to Thomas Secretin for pointing out a transcription error in the
    !! setting of the value of P, 11 February 2013.
    !!
    !! Modified: 11 February 2013
    !!
    !!  Author:
    !!
    !!    - Original FORTRAN77 version by Richard Brent.
    !!    - FORTRAN90 version by John Burkardt.
    !!
    !!  Reference:
    !!
    !!    - Richard Brent,
    !!    - Algorithms for Minimization Without Derivatives,
    !!    - Dover, 2002,
    !!    - ISBN: 0-486-41998-3,
    !!    - LC: QA402.5.B74.
    class(brent_t), intent(in) :: obj !! external object with the user-supplied 
                                   !! subroutine which calculates the functions.
    real(dp), intent(in) :: a,b    !! the endpoints of the change of sign interval
    real(dp), intent(in) :: machep !! an estimate for the relative machine precision.
    real(dp), intent(in) :: t      !! a positive error tolerance
    real(dp) :: zero   !! the estimated value of a zero of the function F

    real(dp) :: c
    real(dp) :: d
    real(dp) :: e
    real(dp) :: fa
    real(dp) :: fb
    real(dp) :: fc
    real(dp) :: m
    real(dp) :: p
    real(dp) :: q
    real(dp) :: r
    real(dp) :: s
    real(dp) :: sa
    real(dp) :: sb
    real(dp) :: tol

    ! Make local copies of A and B.
    sa = a
    sb = b
    fa = obj%f(sa)
    fb = obj%f(sb)
  
    c = sa
    fc = fa
    e = sb - sa
    d = e
    do    
        if (abs(fc) < abs(fb)) then    
            sa = sb
            sb = c
            c = sa
            fa = fb
            fb = fc
            fc = fa    
        endif
  
        tol = 2.0_dp * machep * abs(sb) + t
        m = 0.5_dp * (c - sb)
        if (abs(m) <= tol .or. fb == 0.0_dp) then
            exit
        endif
  
        if (abs(e) < tol .or. abs (fa) <= abs(fb)) then    
            e = m
            d = e    
        else    
            s = fb / fa    
            if (sa == c) then    
                p = 2.0_dp * m * s
                q = 1.0_dp - s    
            else    
                q = fa / fc
                r = fb / fc
                p = s * (2.0_dp * m * q * (q - r) - (sb - sa) * (r - 1.0_dp))
                q = (q - 1.0_dp) * (r - 1.0_dp) * (s - 1.0_dp)    
            endif
  
            if (0.0_dp < p) then
                q = - q
            else
                p = - p
            endif
  
            s = e
            e = d
            if (2.0_dp * p < 3.0_dp * m * q - abs(tol * q) .and. &
                                    p < abs (0.5_dp * s * q)) then
                d = p / q
            else
                e = m
                d = e
            endif
        endif
  
        sa = sb
        fa = fb
        if (tol < abs(d)) then
            sb = sb + d
        else if (0.0_dp < m) then
            sb = sb + tol
        else
            sb = sb - tol
        endif
  
        fb = obj%f(sb)
        if ((0.0_dp < fb .and. 0.0_dp < fc ) .or. &
                (fb <= 0.0_dp .and. fc <= 0.0_dp)) then
            c = sa
            fc = fa
            e = sb - sa
            d = e
        endif
    enddo
    zero = sb
end function zero

end module lib_ext_math_m    
