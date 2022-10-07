module mod_grhd_phys_convert
  use mod_physics
  use mod_grhd_phys_parameters

  implicit none
  private

  double precision                ::  h0           ! lower bound for the relativistic enthalpy over entire validity region of eos
  double precision                ::  one_over_h0  ! 1 / h0

  ! Public methods
  public :: grhd_phys_convert_init

contains

  !> Initialize the module
  subroutine grhd_phys_convert_init()
    use mod_eos, only: eos_hmin

    ! initialize 1/h0
    h0 = eos_hmin
    one_over_h0 = 1.0d0 / h0

    ! initialize the k_max and v_max
    v_max = dsqrt(1.0d0 - 1.0d0 / lfac_max**2)
    k_max = 2.0d0 * v_max / (1.0d0 + v_max**2) ! k_max < 1

    phys_to_conserved        => grhd_to_conserved
    select case (c2p_method)
    case (c2p_galeazzi)
       phys_to_primitive        => grhd_to_primitive_galeazzi
    case (c2p_kastaun)
       phys_to_primitive        => grhd_to_primitive_kastaun
    end select
  end subroutine grhd_phys_convert_init

  !> Transform primitive variables into conservative ones
  subroutine grhd_to_conserved(ixI^L, ixO^L, cons, prim, x)
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
       call grhd_handle_small_values(prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), ixI^L, ixO^L, .True., 'grhd_to_conserved')
    end if

    call grhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                lfac=lfac(ixI^S), &
                h=h(ixI^S) )

    ! Conserved density D
    cons(ixO^S, D_) = lfac(ixO^S) * prim(ixO^S, rho_) 

    ! Convert velocity to covariant momentum
    do idir = 1, ndir
       cons(ixO^S, mom(idir)) = cons(ixO^S, D_) * h(ixO^S) * prim(ixO^S, W_vel(idir)) * gamma(ixO^S,idir,idir)
    end do

    ! Conserved energy - D = tau
    cons(ixO^S, tau_) = cons(ixO^S, D_) * lfac(ixO^S) * h(ixO^S) &
                            - prim(ixO^S, press_) - cons(ixO^S, D_) 

    ! conformal transformation
    psi6(ixO^S) = prim(ixO^S, psi_)**6 
    cons(ixO^S, D_) = cons(ixO^S, D_) * psi6(ixO^S)
    cons(ixO^S, tau_) = cons(ixO^S, tau_) * psi6(ixO^S)
    do idir = 1, ndir
       cons(ixO^S, mom(idir)) = cons(ixO^S, mom(idir)) * psi6(ixO^S)
    end do
  end subroutine grhd_to_conserved

  !> Transform conservative variables into primitive ones 
  subroutine grhd_to_primitive_galeazzi(ixI^L, ixO^L, cons, prim, x)
    use mod_global_parameters
    use mod_geometry
    use mod_small_values
    use mod_eos
    use mod_rootfinding
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
    double precision, intent(inout) :: prim(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer                         :: itr, idir, flag(ixI^S)
    double precision                :: psi6(ixI^S)
    double precision                :: cons_tmp(ixI^S,ncons)
    integer                         :: ix^D 
    integer                         :: error_code = -1
    logical, save                   :: adjustment
    ! metric gamma_ij_hat
    double precision                :: gamma(ixI^S,1:3,1:3)

    ! temp conserved variables (at one grid)
    double precision                :: S, tau, k
    double precision, save          :: D, q, r

    double precision                :: v_bar ! temp velocity
    double precision                :: lfac ! Lorentz factor
    double precision                :: ui_bar(1:ndir)
    double precision                :: rho_bar, eps_bar
    double precision                :: a_bar, h_bar
    double precision                :: eps_min, eps_max
    double precision                :: rescale_factor
    double precision                :: z ! this is the root of the equation
    double precision                :: zp, zm

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

    !Note that v is contravariant
    {do ix^D = ixO^LIM^D \}
      adjustment = .False.
      if ( (cons_tmp(ix^D, D_) <= small_D) .or. &
           (cons_tmp(ix^D, tau_) <= small_tau) ) then
         ! atmosphere handling
         ! skip the primitive recovery
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
         q = max( tau / D, -1.0d0 )  ! limit q to be non-negative
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
            end if
          
            ! range of the root
            zm = 0.5d0 * k / dsqrt(1.0d0 - k**2 / 4.0d0)
            zp = k / dsqrt(1.0d0 - k**2)
            if ( dabs(func(zm)) <= tolerance ) then
               z = zm
            else if ( dabs(func(zp)) <= tolerance ) then
               z = zp
            else
               !call rootfinding_illinois(z, zm, zp, tolerance, iter_max, error_code, func)
               call rootfinding_brent(z, zm, zp, tolerance, iter_max, error_code, func)
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
         end do
   
         ! check if we need any adjustments
         if ( rho_bar <= small_rho_thr ) then
            ! reset the primitive variables
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
         end do
         call grhd_update_eos_one_point(prim(ix^D, 1:nprim))       
       
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
    {end do^D&\}

    !check if NaN
    if (any(flag(ixO^S) < 0)) then
       flag(ixO^S) = rho_
       call small_values_error(prim, x, ixI^L, ixO^L, flag, 'grhd_to_primitive: (NaN z)')
    end if

    contains
       subroutine get_vars(z, W_out, rho_out, eps_out, h_out)
          implicit none
          double precision, intent(in)              :: z
          double precision, intent(out), optional   :: W_out
          double precision, intent(out), optional   :: rho_out
          double precision, intent(out), optional   :: eps_out
          double precision, intent(out), optional   :: h_out

          double precision                          :: W_bar, eps_bar, a_bar, rho_bar, prs_bar, eps_min, eps_max
          W_bar = dsqrt(1.0d0+z**2)
          if (present(W_out)) W_out = W_bar
          rho_bar = D / W_bar
          if ( rho_bar < eos_rhomin .or. rho_bar > eos_rhomax ) then
             rho_bar = max( min( eos_rhomax, rho_bar ), eos_rhomin )
             adjustment = .True. ! adjust conserved variables
          end if
          if (present(rho_out)) rho_out = rho_bar
          eps_bar = W_bar * q - z * r + z**2/(1.0d0+W_bar)
          call eos_get_eps_range(rho_bar,eps_min,eps_max)
          if ( eps_bar < eps_min .or. eps_bar > eps_max ) then
             eps_bar = max( min( eps_max, eps_bar ), eps_min )
             adjustment = .True. ! adjust conserved variables
          end if
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

  end subroutine grhd_to_primitive_galeazzi

  !> Transform conservative variables into primitive ones 
  !> although this is Kastaun's con2prim, magnetic fields are ignored here
  subroutine grhd_to_primitive_kastaun(ixI^L, ixO^L, cons, prim, x)
    use mod_global_parameters
    use mod_geometry
    use mod_small_values
    use mod_eos
    use mod_rootfinding
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
    double precision, intent(inout) :: prim(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer                         :: idir, flag(ixI^S)
    double precision                :: psi6(ixI^S)
    double precision                :: cons_tmp(ixI^S,1:ncons)
    integer                         :: ix^D 
    double precision                :: gamma(ixI^S,1:3,1:3) ! metric gamma_ij_hat

    ! variables needed during the root founding
    ! rescaled variables
    double precision, save          :: q, D
    double precision                :: r_i(1:ndir) !r_i
    double precision, save          :: r2, v0_sqr
    double precision                :: v_hat_sqr
    double precision                :: vi_hat(1:ndir)

    double precision                :: rescale_factor
    logical, save                   :: adjustment
    integer                         :: error_code
    double precision                :: eps_min, eps_max
    double precision                :: W_hat, eps_hat, rho_hat, prs_hat
    double precision                :: mu_plus ! upper bound of mu
    double precision                :: mu ! this is the root of the equation
    ! extra variables needed after the root founding
    double precision                :: h_hat

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

    !Note that v is contravariant
    {do ix^D = ixO^LIM^D \}
      adjustment = .False.
      if ( (cons_tmp(ix^D, D_) <= small_D) .or. &
           (cons_tmp(ix^D, tau_) <= small_tau) ) then
         ! atmosphere handling 
         ! skip the primitive recovery
         W_hat = 1.0d0
         vi_hat(1:ndir) = 0.0d0
         rho_hat = small_rho
         eps_hat = small_eps
         adjustment = .True. ! adjust all conserved variables
      else 
         D = cons_tmp(ix^D, D_) ! conserved density
         ! rescale the conserved variables
         q = cons_tmp(ix^D, tau_) / cons_tmp(ix^D, D_)
         do idir = 1, ndir
            r_i(idir) = cons_tmp(ix^D, mom(idir)) / cons_tmp(ix^D, D_)
         end do
   
         ! some useful relations
         r2 = 0.0d0
         do idir = 1, ndir
            r2 = r2 + r_i(idir)**2 / gamma(ix^D,idir,idir)
         end do
           
         v0_sqr = 1.0d0 / ( h0**2 + r2 )
         mu_plus = dsqrt(v0_sqr) ! when B fields vanish, mu+ is known directly
         v0_sqr = r2 * v0_sqr
   
         ! solve mu
         if ( dabs(func_for_mu(mu_plus)) <= smalldouble ) then
            mu = mu_plus
         else
            call rootfinding_brent(mu, 0.0d0, mu_plus*(1.0d0+smalldouble), tolerance, iter_max, error_code, func_for_mu)
            ! check mu
            select case (error_code)
            !case (0) ! root is found
            case (-1) ! nothing happened
               call mpistop("have you ever attemp to find mu in con2prim?")
            case (1) ! failed to find the root
               call mpistop("Fail to find the root for mu")
            case (2) ! root is NaN
               ! This special flag is used to indicate NaN error. 
               ! Normally it wont be negative
               flag(ix^D) = -1 
            case (3) ! z is not bracketed
               call mpistop("the root is not bracketed for mu")
            end select
         end if
   
         ! Step 3: calculate all the primitive variables from the solution mu
         call get_vars(mu, W_out=W_hat, v_hat_sqr_out=v_hat_sqr, rho_out=rho_hat, eps_out=eps_hat)
         ! calculate velocities
         do idir = 1, ndir
            vi_hat(idir) = mu * ( r_i(idir) / gamma(ix^D,idir,idir) )
         end do
         
         ! adjust the results if it is invalid
         if ( rho_hat <= small_rho_thr ) then
            ! reset the primitive variables
            W_hat = 1.0d0
            vi_hat(1:ndir) = 0.0d0
            rho_hat = small_rho
            eps_hat = small_eps
            adjustment = .True. ! adjust conserved variables
         else
            ! limit the velocities
            if ( W_hat > lfac_max ) then
               ! rescale the velocity such that v = v_max and keeping D constant
               rescale_factor = v_max / dsqrt( v_hat_sqr )
               vi_hat(1:ndir) = vi_hat(1:ndir) * rescale_factor
               W_hat = lfac_max
               ! although D is kept constant, density is changed a bit
               rho_hat = cons_tmp(ix^D, D_) / W_hat
               adjustment = .True. ! adjust conserved variables
            end if

            ! check if eps fails into the validity range
            call eos_get_eps_range(rho_hat, eps_min, eps_max)
            if ( eps_hat < eps_min ) then
               eps_hat = max( eps_hat, eps_min )
               ! This special flag is used to indicate that 
               ! we will only update tau without changing D and S
               flag(ix^D) = 1
               adjustment = .True. ! adjust conserved variables
            else if ( eps_hat > eps_max ) then
               eps_hat = min( eps_max, eps_hat )
               adjustment = .True. ! adjust all conserved variables
            end if

         end if
      end if ! end of con2prim
   
      ! update the primitive variables,
      prim(ix^D, rho_)          = rho_hat
      prim(ix^D, eps_)          = eps_hat
      prim(ix^D, W_vel(1:ndir)) = vi_hat(1:ndir) * W_hat
      call grhd_update_eos_one_point(prim(ix^D, 1:nprim))       

      ! if any adjustments, recompute pressure and conserved varables consistently
      if ( adjustment ) then
         h_hat = 1.0d0 + eps_hat + prim(ix^D, press_) / rho_hat ! enthalpy
         prs_hat = prim(ix^D, press_) ! pressure

         ! update conserved variables
         cons_tmp(ix^D, D_) = W_hat * rho_hat
         cons_tmp(ix^D, tau_) = cons_tmp(ix^D, D_) * ( W_hat * h_hat - 1.0d0 ) &
                               - prs_hat

         ! conformal transformation
         cons(ix^D, tau_) = psi6(ix^D) * cons_tmp(ix^D, tau_) 

         ! we need to update all conserved variables
         if ( flag(ix^D) /= 1 ) then
            ! conformal transformation
            cons(ix^D, D_) = psi6(ix^D) * cons_tmp(ix^D, D_) 
            do idir = 1, ndir
               cons(ix^D, mom(idir)) = cons(ix^D, D_) * h_hat * prim(ix^D, W_vel(idir)) &
                                    * gamma(ix^D,idir,idir)
            end do
         end if
      end if ! end adjustments for conserved

    {end do^D&\}

    !check if NaN
    if (any(flag(ixO^S) < 0)) then
       flag(ixO^S) = rho_
       call small_values_error(prim, x, ixI^L, ixO^L, flag, 'grhd_to_primitive: (NaN mu)')
    end if

    contains
       subroutine get_vars(mu, f_mu, W_out, v_hat_sqr_out, rho_out, eps_out, h_out)
          implicit none
          double precision, intent(in)              :: mu
          double precision, intent(out), optional   :: f_mu
          double precision, intent(out), optional   :: W_out, v_hat_sqr_out
          double precision, intent(out), optional   :: rho_out
          double precision, intent(out), optional   :: eps_out
          double precision, intent(out), optional   :: h_out

          double precision                          :: eps_min, eps_max
          double precision                          :: W_hat, eps_hat, a_hat, rho_hat, prs_hat
          double precision                          :: v_hat_sqr
          double precision                          :: mu_hat
          double precision                          :: nu_A, nu_B, nu_hat

          v_hat_sqr = min( mu**2 * r2 , v0_sqr)
          if (present(v_hat_sqr_out)) v_hat_sqr_out = v_hat_sqr
          W_hat = 1.0d0 / dsqrt( 1.0d0 - v_hat_sqr )
          if (present(W_out)) W_out = W_hat
          rho_hat = D / W_hat
          if ( rho_hat < eos_rhomin .or. rho_hat > eos_rhomax ) then
             rho_hat = max( min( eos_rhomax, rho_hat ), eos_rhomin )
             adjustment = .True. ! adjust conserved variables
          end if
          if (present(rho_out)) rho_out = rho_hat
          eps_hat = W_hat * ( q - mu * r2 ) &
                    + v_hat_sqr * W_hat**2 / ( 1.0d0 + W_hat )
          call eos_get_eps_range(rho_hat,eps_min,eps_max)
          if ( eps_hat < eps_min .or. eps_hat > eps_max ) then
             eps_hat = max( min( eps_max, eps_hat ), eps_min )
             adjustment = .True. ! adjust conserved variables
          end if
          if (present(eps_out)) eps_out = eps_hat
          call eos_get_pressure_one_grid(prs_hat, rho_hat, eps_hat)
 
          if (present(f_mu)) then
             a_hat = prs_hat / (rho_hat*(1.0d0+eps_hat))
             nu_A = (1.0d0 + a_hat) * (1.0d0+eps_hat) / W_hat ! h/W
             nu_B = (1.0d0 + a_hat) * ( 1.0d0 + q - mu * r2 ) 
             nu_hat = max(nu_A,nu_B)
             mu_hat = 1.0d0 / ( nu_hat + r2 * mu )
             f_mu = mu - mu_hat
          end if
       end subroutine get_vars          

       !> master function f(mu) for finding the root mu
       double precision function func_for_mu(mu)
          double precision, intent(in)    :: mu
          call get_vars(mu, f_mu=func_for_mu)
       end function func_for_mu
  end subroutine grhd_to_primitive_kastaun
  
end module mod_grhd_phys_convert
