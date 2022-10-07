module mod_grrhd_phys_convert
  use mod_physics
  use mod_grrhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grrhd_phys_convert_init

contains

  !> Initialize the module
  subroutine grrhd_phys_convert_init()
    ! initialize the k_max and v_max
    v_max = dsqrt(1.0d0 - 1.0d0 / lfac_max**2)
    k_max = 2.0d0 * v_max / (1.0d0 + v_max**2) ! k_max < 1
    phys_to_conserved        => grrhd_to_conserved
    phys_to_primitive        => grrhd_to_primitive
  end subroutine grrhd_phys_convert_init

  !> Transform primitive variables into conservative ones
  subroutine grrhd_to_conserved(ixI^L, ixO^L, cons, prim, x)
    use mod_global_parameters
    use mod_eos, only: small_D, small_tau
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: cons(ixI^S, ncons)
    double precision, intent(inout) :: prim(ixI^S, nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision                :: h(ixI^S) ! enthalpy
    double precision                :: lfac(ixI^S) ! Lorentz factor
    double precision                :: psi6(ixI^S) ! conformal factor **6
    double precision                :: gamma(ixI^S,1:3,1:3)
    integer                         :: idir

    if ( fix_small_values ) then
       call grrhd_handle_small_values(prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), ixI^L, ixO^L, .True., 'grrhd_to_conserved')
    end if

    call grrhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                lfac=lfac(ixI^S), &
                h=h(ixI^S) )

    ! Conserved density D
    cons(ixO^S, D_) = lfac(ixO^S) * prim(ixO^S, rho_) 
    where ( cons(ixO^S, D_) < small_D )
       cons(ixO^S, D_) = small_D
    end where

    ! Convert velocity to covariant momentum
    do idir = 1, ndir
       cons(ixO^S, mom(idir)) = cons(ixO^S, D_) * h(ixO^S) * prim(ixO^S, W_vel(idir)) * gamma(ixO^S,idir,idir)
    end do

    ! Conserved energy - D = tau
    cons(ixO^S, tau_) = cons(ixO^S, D_) * lfac(ixO^S) * h(ixO^S) &
                            - prim(ixO^S, press_) - cons(ixO^S, D_) 
    where ( cons(ixO^S, tau_) < small_tau )
       cons(ixO^S, tau_) = small_tau
    end where

    ! conformal transformation
    psi6(ixO^S) = prim(ixO^S, psi_)**6 
    cons(ixO^S, D_) = cons(ixO^S, D_) * psi6(ixO^S)
    cons(ixO^S, tau_) = cons(ixO^S, tau_) * psi6(ixO^S)
    do idir = 1, ndir
       cons(ixO^S, mom(idir)) = cons(ixO^S, mom(idir)) * psi6(ixO^S)
    end do

    ! fixme: conformal factor
    ! Conserved energy density nu_E
    cons(ixO^S, nu_Econs) = prim(ixO^S, nu_E) 
    ! covariant momentum
    ! F_i/E -> F_i
    do idir = 1, ndir
       cons(ixO^S, nu_Fcons(idir)) = prim(ixO^S, nu_F_over_E(idir)) * prim(ixO^S, nu_E) 
    end do
  end subroutine grrhd_to_conserved

  !> Transform conservative variables into primitive ones 
  subroutine grrhd_to_primitive(ixI^L, ixO^L, cons, prim, x)
    use mod_global_parameters
    use mod_geometry
    use mod_small_values
    use mod_eos
    use mod_rootfinding
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
    double precision, intent(inout) :: prim(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer                         :: idir, jdir, flag(ixI^S)
    double precision                :: psi6(ixI^S)
    double precision                :: cons_tmp(ixI^S,ncons)
    integer                         :: ix^D 
    integer                         :: error_code = -1
    logical                         :: adjustment
    ! metric gamma_ij_hat
    double precision                :: gamma(ixI^S,1:3,1:3)

    ! --------------   variables needed in con2prim  ------------- 
    double precision                :: D, S, tau
    double precision                :: r, q, k

    double precision                :: v_bar ! temp velocity
    double precision                :: lfac ! Lorentz factor
    double precision                :: ui_bar(1:ndir)
    double precision                :: rho_bar, eps_bar
    double precision                :: a_bar, h_bar
    double precision                :: eps_min, eps_max
    double precision                :: rescale_factor
    double precision                :: z ! this is the root of the equation
    double precision                :: zp, zm

    ! --------------   variables needed for M1 scheme  ------------- 
    double precision                :: zeta(ixI^S)

    flag(ixO^S) = 0
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
    end do
    
    psi6(ixO^S) = prim(ixO^S, psi_)**6 
    ! conformal transformation back to normal conserved variables
    cons_tmp(ixO^S, D_) = cons(ixO^S, D_) / psi6(ixO^S)
    cons_tmp(ixO^S, tau_) = cons(ixO^S, tau_) / psi6(ixO^S)
    do idir = 1, ndir
       cons_tmp(ixO^S, mom(idir)) = cons(ixO^S, mom(idir)) / psi6(ixO^S)
    end do
   
    ! fixme: conformal factor
    ! Conserved energy density nu_E
    prim(ixO^S, nu_E) = cons(ixO^S, nu_Econs)
    ! covariant momentum
    ! F_i/E -> F_i
    do idir = 1, ndir
       prim(ixO^S, nu_F_over_E(idir)) = cons(ixO^S, nu_Fcons(idir)) / prim(ixO^S, nu_E) 
    end do

    !Note that v is contravariant
    {do ix^D = ixO^LIM^D \}
      adjustment = .False.
      if ( (cons_tmp(ix^D, D_) <= small_D) .or. &
           (cons_tmp(ix^D, tau_) <= small_tau) ) then
         ! atmosphere handling
         ! skip the primitive recovery
         v_bar = 0.0d0
         lfac = 1.0d0
         ui_bar  = 0.0d0
         rho_bar = small_rho
         eps_bar = small_eps
         adjustment = .True. ! adjust all conserved variables
      else
         D = cons_tmp(ix^D, D_)
         S = 0.0d0
         do idir = 1, ndir
            S = S + cons_tmp(ix^D, mom(idir))**2 / gamma(ix^D,idir,idir)
         end do
         S = dsqrt( S )
         tau = cons_tmp(ix^D, tau_)
   
         r = S / D
         q = max( tau / D, 0.0d0 )  ! limit q to be non-negative
         k = S / (tau + D) ! k = r / (1+q)
   
         ! find the root z
         if (dabs(k) <= smalldouble) then
            z = 0.0d0
         else
            if (k > k_max) then 
               ! rescale the momentum such that r=k_max*(1+q)
               rescale_factor = k_max / k
               r = k_max * (1.0d0 + q)
               k = k_max
               cons_tmp(ix^D, mom(:)) = cons_tmp(ix^D, mom(:)) * rescale_factor
               cons(ix^D, mom(:)) = cons(ix^D, mom(:)) * rescale_factor
            endif
          
            ! range of the root
            zm = 0.5d0 * k / dsqrt(1.0d0 - k**2 / 4.0d0)
            zp = k / dsqrt(1.0d0 - k**2)
            if ( dabs(func(zm)) <= tolerance ) then
               z = zm
            else if ( dabs(func(zp)) <= tolerance ) then
               z = zp
            else
               call rootfinding_illinois(z, zm, zp, tolerance, iter_max, error_code, func)
               ! after getting z
               ! check if there are any errors or correction
               select case (error_code)
               !case (0) ! root is found
               !   write(*,*) "z= ", z
               case (-1) ! nothing happened
                  call mpistop("have you ever attemp to find the root in con2prim?")
               case (1) ! failed to find the root
                  call mpistop("Fail to find the root in con2prim")
               case (2) ! z is NaN
                  call mpistop("NaN")
                  flag(ix^D) = -1 ! This special flag is used to indicate NaN error. Normally it wont be negative
               case (3) ! z is not bracketed
                  call mpistop("the root is not bracketed in con2prim")
               end select
            endif 
         endif 
   
         ! calculate all the primitive variables from the solution z
         call get_vars(z, W_out=lfac, rho_out=rho_bar, eps_out=eps_bar, h_out=h_bar)
         do idir=1,ndir
            ui_bar(idir) = cons_tmp(ix^D, mom(idir)) / gamma(ix^D,idir,idir) &
                        / ( h_bar * D ) 
         enddo
   
         ! check if we need any adjustments
         if ( rho_bar <= small_rho_thr ) then
            ! reset the primitive variables
            v_bar = 1.0d0
            lfac = 1.0d0
            ui_bar = 0.0d0
            rho_bar = small_rho
            eps_bar = small_eps
            adjustment = .True. ! adjust conserved variables
         else
            ! check if eps fails into the validity range
            call eos_get_eps_range(rho_bar, eps_min, eps_max)
            if ( eps_bar < eps_min .or. eps_bar > eps_max ) then
               eps_bar = max( min( eps_max, eps_bar ), eps_min )
               adjustment = .True. ! adjust conserved variables
            end if
      
            ! limit the velocities
            v_bar = z / lfac
            if (v_bar > v_max) then 
               ! rescale the velocity such that v = v_max and keeping D constant
               rescale_factor = ( v_max / v_bar ) * ( lfac_max / lfac )
               ui_bar = ui_bar * rescale_factor
               v_bar = v_max
               lfac = lfac_max
      
               ! although D is kept constant, density is changed a bit
               rho_bar = D / lfac
               ! limit eps again based on the updated rho
               call eos_get_eps_range(rho_bar, eps_min, eps_max)
               eps_bar = max( min( eps_max, eps_bar ), eps_min )
               adjustment = .True. ! adjust conserved variables
            end if
         end if

         ! update all the hydro primitive variables here
         prim(ix^D, rho_) = rho_bar
         prim(ix^D, eps_) = eps_bar
         do idir=1,ndir
            prim(ix^D, W_vel(idir)) = ui_bar(idir)
         enddo
         call grrhd_update_eos_one_point(prim(ix^D, 1:nprim))       
       
         ! if any adjustments, recompute conserved varables consistently
         if ( adjustment ) then
            a_bar = prim(ix^D, press_) / (prim(ix^D, rho_)*(1.0d0+prim(ix^D, eps_)))
            a_bar = max( min( 1.0d0, a_bar ) , 0.0d0) ! strictly require a is physical
            h_bar = (1.0d0+prim(ix^D, eps_))*(1.0d0+a_bar)
            ! Conserved density D
            cons(ix^D, D_) = lfac * prim(ix^D, rho_) 
            ! Convert velocity to covariant momentum
            do idir = 1, ndir
               cons(ix^D, mom(idir)) = cons(ix^D, D_) &
                                     * h_bar * prim(ix^D, W_vel(idir)) * gamma(ix^D,idir,idir)
            end do
            ! Conserved energy - D = tau
            cons(ix^D, tau_) = cons(ix^D, D_) * lfac * h_bar &
                                    - prim(ix^D, press_) - cons(ix^D, D_) 
            ! conformal transformation
            cons(ix^D, D_) = cons(ix^D, D_) * psi6(ix^D)
            cons(ix^D, tau_) = cons(ix^D, tau_) * psi6(ix^D)
            do idir = 1, ndir
               cons(ix^D, mom(idir)) = cons(ix^D, mom(idir)) * psi6(ix^D)
            end do
         end if

      end if
    {enddo^D&\}
    
    if ( use_radiation ) then
      call grrhd_get_nu_zeta(ixI^L, ixO^L, prim, x)
    end if 
   
    !check if NaN
    if (any(flag(ixO^S) < 0)) then
       flag(ixO^S) = rho_
       call small_values_error(prim, x, ixI^L, ixO^L, flag, 'grrhd_to_primitive: (NaN z)')
    end if
    
    contains
       subroutine get_vars(z, W_out, rho_out, eps_out, h_out)
          implicit none
          double precision, intent(in)    :: z
          double precision, intent(out), optional   :: W_out
          double precision, intent(out), optional   :: rho_out
          double precision, intent(out), optional   :: eps_out
          double precision, intent(out), optional   :: h_out

          double precision                :: W_bar, eps_bar, a_bar, rho_bar, prs_bar
          W_bar = dsqrt(1.0d0+z**2)
          if (present(W_out)) W_out = W_bar
          rho_bar = D / W_bar
          rho_bar = max( min( eos_rhomax, rho_bar ), eos_rhomin )
          if (present(rho_out)) rho_out = rho_bar
          eps_bar = W_bar * q - z * r + z**2/(1.0d0+W_bar)
          call eos_get_eps_range(rho_bar,eps_min,eps_max)
          eps_bar = max( min( eps_max, eps_bar ), eps_min )
          if (present(eps_out)) eps_out = eps_bar
          if (present(h_out)) then
             call eos_get_pressure_one_grid(prs_bar, rho_bar, eps_bar)
             a_bar = prs_bar / (rho_bar*(1.0d0+eps_bar))
             a_bar = max( min( 1.0d0, a_bar ) , 0.0d0)
             h_out = (1.0d0+eps_bar)*(1.0d0+a_bar); !h_bar = max( hbar, 0.0d0 )
          end if
       end subroutine get_vars          

       !> master function f(z) for finding the root z
       function func(z)
          implicit none
          double precision :: func
          double precision, intent(in)    :: z
          double precision                :: h_bar
          call get_vars(z, h_out=h_bar)
          func = z - r / h_bar
       end function func

  end subroutine grrhd_to_primitive
  
end module mod_grrhd_phys_convert
