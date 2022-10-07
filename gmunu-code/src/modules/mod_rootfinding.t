module mod_rootfinding
  implicit none
contains

  recursive subroutine rootfinding_brent(z, zmin, zmax, tolerance, iter_max, return_code, func)
    !> root
    double precision, intent(out):: z
    !> the range of the root
    double precision, intent(in) :: zmin, zmax
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximum iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = the root is not bracketed;
    integer, intent(out) :: return_code
    interface
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in) :: z
       end function func
    end interface

    integer :: it_root
    double precision :: a, b, c, d, e
    double precision :: fa, fb, fc
    double precision :: p, q, r, s
    double precision :: xm, tol
    double precision, parameter :: eps = epsilon(zmin) ! machine precision

    return_code = -1

    a  = zmin
    b  = zmax
    fa = func(a)
    fb = func(b)

    ! check if the bounds are already the root
    if ( dabs(fa) == 0.0d0 ) then
       z = a
       return_code = 0
       return
    else if ( dabs(fb) == 0.0d0 ) then
       z = b
       return_code = 0
       return
    else if ( fa * fb > 0.0d0 ) then
       ! the root is not bracketed
       return_code = 3
       return
       !call mpistop("error in brent: the root is not bracketed!")
    end if

    c  = b
    fc = fb

    do it_root = 1, iter_max

       if ( fb * fc > 0.0d0 ) then
          ! Rename a, b, c and adjust bounding interval d
          c  = a
          fc = fa
          d  = b - a
          e  = d
       end if
 
       if ( dabs(fc) < dabs(fb) ) then
          a  = b
          b  = c
          c  = a
          fa = fb
          fb = fc
          fc = fa
       end if

       ! Convergence check.
       tol = 2.0d0 * eps * dabs(b) + 0.5d0 * tolerance
       xm  = 0.5d0 * ( c - b )
       if ( dabs(xm) <= tol .or. fb == 0.0d0 ) then
          z = b
          return_code = 0
          return
       end if

       if ( dabs(e) >= tol .and. dabs(fa) > dabs(fb) ) then
          ! Attempt inverse quadratic interpolation
          s = fb / fa
          if (a == c) then
             p = 2.0d0 * xm * s 
             q = 1.0d0 - s
          else
             q = fa / fc
             r = fb / fc
             p = s * (2.0d0 * xm * q * (q-r) - (b-a) * (r - 1.0d0))
             q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
          end if
          if (p > 0.0d0) q = -q  ! Check whether in bounds.
          p = dabs(p)
          if ( 2.0d0 * p < min(3.0d0 * xm * q - dabs(tol*q), dabs(e*q)) ) then
             ! Accept interpolation
             e = d
             d = p / q
          else
             ! Interpolation failed, use bisection
             d = xm
             e = d
          end if
       else
          ! Bounds decreasing too slowly, use bisection
          d = xm
          e = d
       end if
       a  = b  ! move last best guess to a
       fa = fb
       b  = b + merge(d, sign(tol,xm), dabs(d) > tol) ! evaluate new trail root
       fb = func(b)
    end do
    ! if the root is not found
    return_code = 1
  end subroutine rootfinding_brent

  ! note: making this recursive since this subroutine will call itself indirectly in some cases.
  recursive subroutine rootfinding_illinois(z, zmin, zmax, tolerance, iter_max, return_code, func)
    implicit none
    !> root
    double precision, intent(out):: z
    !> the range of the root
    double precision, intent(in) :: zmin, zmax
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximun iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = the root is not bracketed;
    integer, intent(out) :: return_code
    interface
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in) :: z
       end function func
    end interface

    integer :: it_root
    integer :: side
    double precision :: zm, zp
    double precision :: f, fm, fp
    double precision :: z_new
    
    side = 0
    return_code = -1

    zm = zmin
    zp = zmax
    fm = func(zm)
    fp = func(zp)

    ! check if the bounds are already the root
    if ( dabs(fm) == 0.0d0 ) then
       z = zm
       return_code = 0
       return
    else if ( dabs(fp) == 0.0d0 ) then
       z = zp
       return_code = 0
       return
    else if ( fm * fp > 0.0d0 ) then
       ! check if the root is bracketed
       !write(*,*) "Warning in illinois: the root is not bracketed!"
       !write(*,*) "root range = ", zm, zp
       !write(*,*) "f(root) range = ", fm, fp
       return_code = 3
       return
       !call mpistop("error in illinois: the root is not bracketed!")
    end if

    do it_root = 1, iter_max

       z = (fm * zp - fp * zm) / (fm - fp)
       f = func(z)
       
       if ( dabs(zp-zm) <= (tolerance * 0.5d0 * dabs(zp+zm)) .or. &
            f == 0.0d0 ) then
          return_code = 0
          return
       end if

       if ( (f * fp) > 0.0d0 ) then
          !f and fp have same sign, copy z to zp
          zp = z
          fp = f
          if (side == 1) fm = fm * 0.5d0
          side = 1
       else if ( (f * fm) > 0.0d0 ) then
          !f and fm have same sign, copy z to zm
          zm = z
          fm = f
          if (side == 2) fp = fp * 0.5d0
          side = 2
       else !it looks like zero
          return_code = 0
          return
       end if

       if ( isnan(z) ) then
          return_code = 2
          return
       end if

    end do
    ! if the root is not found
    return_code = 1
    !call mpistop("fail to find the root")
  end subroutine rootfinding_illinois

  subroutine rootfinding_newton_raphson(z, tolerance, iter_max, return_code, func, dev_func)
    !> input as initial guess, output as root
    double precision, intent(inout):: z
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximun iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    integer, intent(out) :: return_code
    interface
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in) :: z
       end function func
       function dev_func(z)
          implicit none
          double precision :: dev_func
          double precision, intent(in) :: z
       end function dev_func
    end interface

    integer :: it_root
    double precision :: f, df

    return_code = -1

    do it_root = 1, iter_max

       if ( isnan(z) ) then
          return_code = 2
          return
       end if
       
       f  = func(z)
       df = dev_func(z)

       if ( dabs(f)  <  tolerance ) then
          ! the root is found!
          return_code = 0
          return
       end if

       z = z - f / df
    end do
    ! if the root is not found
    return_code = 1
  end subroutine rootfinding_newton_raphson

  subroutine rootfinding_constrained_newton_raphson(z, zm, zp, tolerance, iter_max, return_code, func, dev_func)
    !> input as initial guess, output as root
    double precision, intent(inout):: z
    !> the range of the root
    double precision, intent(inout) :: zm, zp
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximum iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = the root is not bracketed;
    integer, intent(out) :: return_code
    interface
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in) :: z
       end function func
       function dev_func(z)
          implicit none
          double precision :: dev_func
          double precision, intent(in) :: z
       end function dev_func
    end interface

    integer :: it_root
    double precision :: f, df, z_new
    double precision :: fm, fp

    return_code = -1

    fm = func(zm)
    fp = func(zp)
    
    ! check if the bounds are already the root
    if ( dabs(fm) <= tolerance ) then
       z = zm
       return_code = 0
       return
    else if ( dabs(fp) <= tolerance ) then
       z = zp
       return_code = 0
       return
    else if ( fm * fp > 0.0d0 ) then
       ! check if the root is bracketed
       !write(*,*) "root range = ", zm, zp
       !write(*,*) "f(root) range = ", fm, fp
       return_code = 3
       return
       !call mpistop("error in constrained Newton Raphson: the root is not bracketed!")
    end if

    ! initial guess
    if (z > zp .or. z < zm) then
       z_new = ( zp + zm ) * 0.4d0 ! avoid the same initial guess with bisection
    else
       z_new = z ! input as initial guess
    end if

    do it_root = 1, iter_max
       
       z  = z_new
       f  = func(z)
       df = dev_func(z)

       ! Newton-Raphson step
       z_new = z - f / df

       if ( z_new > zp .or. z_new < zm ) then
          ! the root is outside the range, use bisection here
          z_new = ( zp + zm ) * 0.5d0
          f  = func(z_new)
          if ( f * fm > 0.0d0 ) then
             zm = z_new
             fm = f
          else
             zp = z_new
             fp = f
          end if
       end if

       if ( isnan(z_new) ) then
          ! if z_new is NaN, we take the previous z, and quit the iteration
          return_code = 2
          return
       end if

       if ( dabs(f) <= tolerance .or. &
            dabs(z_new-z) <= dabs(0.5d0*(z_new+z)*tolerance) ) then
          z = z_new
          return_code = 0
          return
       end if

    end do
    ! if the root is not found
    return_code = 1
  end subroutine rootfinding_constrained_newton_raphson

  subroutine rootfinding_multid_newton_raphson(z, tolerance, iter_max, &
        return_code, fvec_in)
    use mod_lu
    !> input as initial guess, output as root
    double precision, dimension(:), intent(inout):: z
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximum iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = jacobian is singular;
    integer, intent(out) :: return_code
    interface
       subroutine fvec_in(z_vec, fvec_out, fjac_out)
          implicit none
          double precision, dimension(:), intent(in)               :: z_vec
          double precision, dimension(1:size(z_vec))               :: fvec_out
          double precision, dimension(1:size(z_vec),1:size(z_vec)) :: fjac_out
       end subroutine fvec_in
    end interface

    integer                      :: it_root, n_dim
    integer                      :: err_flag, d_flag
    integer, dimension(size(z))  :: indx
    double precision, dimension(size(z), size(z)) :: jac
    double precision, dimension(size(z))          :: fvec, dz

    return_code = -1
    n_dim = size(z)

    ! set dz as inf at the begining
    dz(:) = huge(0.0d0)

    ! page 280 in F90
    do it_root = 1, iter_max
       
       call fvec_in(z, fvec, jac)

       if ( maxval(dabs(fvec)) <=  tolerance .or. &
            maxval(dabs(dz)) <= tolerance ) then
          ! the root is found!
          return_code = 0
          return
       end if

       dz(:) = -fvec(:)
       call ludcmp(n_dim, jac, indx, d_flag, err_flag)
       if (err_flag==1) then
          return_code = 4
          return
          !call mpistop("Error in Multi-D Newton-Raphson: the jacobian is singular :/")
       end if
       call lubksb(n_dim, jac, indx, dz)
       z = z + dz
    end do
    ! if the root is not found
    return_code = 1
  end subroutine rootfinding_multid_newton_raphson

  ! this combine f and Jac together as a single input
  subroutine rootfinding_global_multid_newton_raphson(z, tolerance, iter_max, &
        return_code, fvec_in, fmin_in)
    use mod_lu
    use mod_lnsrch
    !> input as initial guess, output as root
    double precision, dimension(:), intent(inout):: z
    !> the tolerance
    double precision, intent(in) :: tolerance
    !> maximum iterations
    integer, intent(in)          :: iter_max
    !> return code: 
    !> -1 = nothing happened;
    !> 0 = the root is found;
    !> 1 = the root not found;
    !> 2 = the root is NaN;
    !> 3 = jacobian is singular;
    !> 4 = reached local minimum;
    !> 5 = slope is larger or equals to zero;
    integer, intent(out) :: return_code
    interface
       subroutine fvec_in(z_vec, fvec_out, fjac_out)
          implicit none
          double precision, dimension(:), intent(in)               :: z_vec
          double precision, dimension(1:size(z_vec))               :: fvec_out
          double precision, dimension(1:size(z_vec),1:size(z_vec)) :: fjac_out
       end subroutine fvec_in
       function fmin_in(z_vec)
          implicit none
          double precision :: fmin_in
          double precision, dimension(:), intent(in)  :: z_vec
       end function fmin_in
    end interface

    double precision, parameter                   :: tolz=epsilon(z), STPMX = 1.0d2
    double precision, parameter                   :: tolmin = 1.0d-12
    integer                                       :: it_root, n_dim
    integer                                       :: err_flag, d_flag
    integer, dimension(size(z))                   :: indx
    double precision, dimension(size(z))          :: fvec, dz, zold, g
    double precision, dimension(size(z), size(z)) :: jac
    double precision                              :: f, stpmax, fold

    return_code = -1
    n_dim = size(z)

    ! set dz as inf at the begining
    dz(:)   = huge(0.0d0)
    zold(:) = huge(0.0d0)

    ! calculate the maximum step length allowed in line searches
    !stpmax = huge(0.0d0) ! currently we dont limit the step length
    stpmax = STPMX * max( dsqrt(dot_product(z(:),z(:))), dble(n_dim))

    f = fmin_in(z)

    do it_root = 1, iter_max
       
       call fvec_in(z, fvec, jac)

       if ( (maxval(dabs(fvec)) <= tiny(0.0d0)) .or. &
            (maxval( dabs(z(:)-zold(:))/max(dabs(z(:)), 1.0d0) ) <= tolerance) ) then
          ! the root is found!
          return_code = 0
          return
       end if

       g(:)    = matmul(fvec(:),jac(:,:))  ! grad(f), for the line search
       zold(:) = z(:)
       fold    = f
       dz(:)   = -fvec(:)                  ! put -fvec to the right hand side

       ! solve linear equations by LU decomposition
       call ludcmp(n_dim, jac, indx, d_flag, err_flag)
       if (err_flag==1) then
          !call mpistop("Error in Multi-D Newton-Raphson: the jacobian is singular :/")
          return_code = 4
          return
       end if
       call lubksb(n_dim, jac, indx, dz)

       ! line search, and update new z and fvec
       call lnsrch(zold, fold, g, dz, z, f, stpmax, err_flag, fmin_in)
       ! lnsrch returns new z

       if ( any(isnan(z(:))) ) then
          return_code = 2
          return
       end if

       select case (err_flag)
       case (1)
          !write(*,*) "Warning in Global Multi-D Newton-Raphson: reached local minimum :/"
          if ( (maxval( dabs(g(:))*max(dabs(z(:)), 1.0d0)&
                        / max(f,0.5d0*n_dim) ) <= tolmin) ) then
             return_code = 3 ! suprious convergence, maybe we need another initial guess
          else
             return_code = 0
          end if
          return
       case (2)
          !call mpistop("Error in Multi-D Newton-Raphson: slope is larger then or equals to zero"
          return_code = 2
          return
       end select

    end do
    ! if the root is not found
    return_code = 1
  end subroutine rootfinding_global_multid_newton_raphson

end module mod_rootfinding
