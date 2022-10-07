module mod_gremhd_phys_convert
  use mod_physics
  use mod_gremhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: gremhd_phys_convert_init

contains

  !> Initialize the module
  subroutine gremhd_phys_convert_init()
    ! initialize the k_max and v_max
    v_max = dsqrt(1.0d0 - 1.0d0 / lfac_max**2)
    k_max = 2.0d0 * v_max / (1.0d0 + v_max**2) ! k_max < 1
    phys_to_conserved        => gremhd_to_conserved
    phys_to_primitive        => gremhd_to_primitive
  end subroutine gremhd_phys_convert_init

  !> Transform primitive variables into conservative ones
  subroutine gremhd_to_conserved(ixI^L, ixO^L, cons, prim, x)
    use mod_global_parameters
    use mod_eos, only: small_tau
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
    double precision, intent(inout) :: prim(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision                :: h(ixI^S)    ! enthalpy
    double precision                :: lfac(ixI^S) ! Lorentz factor
    double precision                :: E2(ixI^S), B2(ixI^S) 
    double precision                :: psi6(ixI^S) ! conformal factor **6
    double precision                :: gamma(ixI^S,1:3,1:3)
    double precision                :: sqrt_gamma(ixI^S)
    integer                         :: idir, jdir, kdir

    if ( fix_small_values ) then
       call gremhd_handle_small_values(prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), ixI^L, ixO^L, .True., 'gremhd_to_conserved')
    end if
    psi6(ixO^S) = prim(ixO^S, psi_)**6 

    call gremhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), sqrt_gamma=sqrt_gamma(ixI^S), &
                E2=E2(ixI^S), B2=B2(ixI^S), &
                lfac=lfac(ixI^S), h=h(ixI^S) )

    ! Conserved density D
    cons(ixO^S, D_) = lfac(ixO^S) * prim(ixO^S, rho_) 

    ! Convert velocity to covariant momentum
    do idir = 1, ndir
       cons(ixO^S, mom(idir)) = ( cons(ixO^S, D_) * h(ixO^S) * prim(ixO^S, W_vel(idir)) ) * gamma(ixO^S,idir,idir)
       do jdir = 1, ndir
          do kdir = 1, ndir
             cons(ixO^S, mom(idir)) = cons(ixO^S, mom(idir)) &
                         + sqrt_gamma(ixO^S) * lvc(idir,jdir,kdir) * prim(ixO^S, Evec(jdir)) * prim(ixO^S, Bvec(kdir))
          end do
       end do
    end do

    ! Conserved energy - D = tau
    cons(ixO^S, tau_) = cons(ixO^S, D_) * ( lfac(ixO^S) * h(ixO^S) - 1.0d0 ) - prim(ixO^S, press_) &
                        + 0.5d0 * ( E2(ixO^S) + B2(ixO^S) )
    !where ( cons(ixO^S, tau_) < small_tau )
    !   cons(ixO^S, tau_) = small_tau
    !end where

    ! conformal transformation
    cons(ixO^S, D_) = cons(ixO^S, D_) * psi6(ixO^S)
    cons(ixO^S, tau_) = cons(ixO^S, tau_) * psi6(ixO^S)
    do idir = 1, ndir
       cons(ixO^S, mom(idir))   = cons(ixO^S, mom(idir)) * psi6(ixO^S)
       cons(ixO^S, Bcons(idir)) = prim(ixO^S, Bvec(idir)) * psi6(ixO^S)
       cons(ixO^S, Econs(idir)) = prim(ixO^S, Evec(idir)) * psi6(ixO^S)
    end do
    if (type_divb == glm)  cons(ixO^S, Bphi_cons_) = prim(ixO^S, Bphi_) * psi6(ixO^S)
    if (type_dive == glm)  cons(ixO^S, Ephi_cons_) = prim(ixO^S, Ephi_) * psi6(ixO^S)
  end subroutine gremhd_to_conserved

  !> Transform conservative variables into primitive ones 
  subroutine gremhd_to_primitive(ixI^L, ixO^L, cons, prim, x)
    use mod_global_parameters
    use mod_geometry
    use mod_small_values
    use mod_eos
    use mod_rootfinding
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
    double precision, intent(inout) :: prim(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer                         :: idir, jdir, kdir
    integer                         :: flag(ixI^S)
    double precision                :: psi6(ixI^S)
    double precision                :: cons_tmp(ixI^S,ncons)
    integer                         :: ix^D 
    double precision                :: gamma(ixI^S,1:3,1:3) ! metric gamma_ij
    double precision                :: sqrt_gamma(ixI^S)
    double precision                :: sqrt_gamma_hat(ixI^S)
    ! variables needed for the root founding
    logical                         :: adjustment
    integer                         :: error_code
    ! temp conserved variables (at one grid)
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

    flag(ixO^S) = 0
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
    end do
    psi6(ixO^S) = prim(ixO^S, psi_)**6 
    call get_sqrt_gamma_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, sqrt_gamma_hat(ixI^S))
    sqrt_gamma(ixO^S) = psi6(ixO^S) * sqrt_gamma_hat(ixO^S)

    ! Here, we assume that both E and B are known 
    do idir = 1, ndir
       prim(ixO^S, Bvec(idir)) = cons(ixO^S, Bcons(idir)) / psi6(ixO^S)
       prim(ixO^S, Evec(idir)) = cons(ixO^S, Econs(idir)) / psi6(ixO^S)
    end do
    if (type_divb == glm) &
        prim(ixO^S, Bphi_) = cons(ixO^S, Bphi_cons_) / psi6(ixO^S)
    if (type_dive == glm) &
        prim(ixO^S, Ephi_) = cons(ixO^S, Ephi_cons_) / psi6(ixO^S)

    ! conformal transformation back to normal conserved variables
    cons_tmp(ixO^S, D_) = cons(ixO^S, D_) / psi6(ixO^S)
    cons_tmp(ixO^S, tau_) = cons(ixO^S, tau_) / psi6(ixO^S)
    do idir = 1, ndir
       cons_tmp(ixO^S, mom(idir)) = cons(ixO^S, mom(idir)) / psi6(ixO^S)
    end do

    ! remove EM part in the conserved variables
    do idir = 1, ndir
       cons_tmp(ixO^S, tau_) = cons_tmp(ixO^S, tau_) &
         - 0.5d0 * gamma(ixO^S, idir, idir) * ( prim(ixO^S, Evec(idir))**2 + prim(ixO^S, Bvec(idir))**2 )
    end do
    do idir = 1, ndir
       jdir=mod(idir,ndir)+1   ! 'Next' direction
       kdir=mod(idir+1,ndir)+1 ! 'Next Next' direction

       cons_tmp(ixO^S, mom(idir)) = cons_tmp(ixO^S, mom(idir)) &
            - sqrt_gamma(ixO^S) * ( prim(ixO^S, Evec(jdir)) * prim(ixO^S, Bvec(kdir)) &
             - prim(ixO^S, Evec(kdir)) * prim(ixO^S, Bvec(jdir)) )
!       do jdir = 1, ndir
!          do kdir = 1, ndir
!             cons_tmp(ixO^S, mom(idir)) = cons_tmp(ixO^S, mom(idir)) &
!                 - sqrt_gamma(ixO^S) * lvc(idir, jdir, kdir) * prim(ixO^S, Evec(jdir)) * prim(ixO^S, Bvec(kdir))
!          end do
!       end do
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
               adjustment = .True. ! adjust all conserved variables
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
               ! rescale the velocity such that v_bar = v_max and keeping D constant
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
      end if
   
      ! update all the hydro primitive variables here
      prim(ix^D, rho_) = rho_bar
      prim(ix^D, eps_) = eps_bar
      do idir=1,ndir
         prim(ix^D, W_vel(idir)) = ui_bar(idir)
      enddo
      call gremhd_update_eos_one_point(prim(ix^D, 1:nprim))       

      ! if any adjustments, recompute conserved varables consistently
      if ( adjustment ) then
         h_bar = 1.0d0 + eps_bar + prim(ix^D, press_) / rho_bar
         cons(ix^D, D_) = lfac * rho_bar
         ! Convert velocity to covariant momentum
         do idir = 1, ndir
            cons(ix^D, mom(idir)) = ( cons(ix^D, D_) * h_bar * prim(ix^D, W_vel(idir)) ) * gamma(ix^D,idir,idir)
            do jdir = 1, ndir
               do kdir = 1, ndir
                  cons(ix^D, mom(idir)) = cons(ix^D, mom(idir)) &
                              + sqrt_gamma(ix^D) * lvc(idir,jdir,kdir) * prim(ix^D, Evec(jdir)) * prim(ix^D, Bvec(kdir))
               end do
            end do
         end do
     
         ! Conserved energy - D = tau
         cons(ix^D, tau_) = cons(ix^D, D_) * ( lfac * h_bar - 1.0d0 ) - prim(ix^D, press_)
         do idir = 1, ndir
            cons(ix^D, tau_) = cons(ix^D, tau_) &
                         + 0.5d0 * ( prim(ix^D, Evec(idir))**2 + prim(ix^D, Bvec(idir))**2 ) * gamma(ix^D, idir, idir)
         end do
     
         ! conformal transformation
         cons(ix^D, D_) = cons(ix^D, D_) * psi6(ix^D)
         cons(ix^D, tau_) = cons(ix^D, tau_) * psi6(ix^D)
         do idir = 1, ndir
            cons(ix^D, mom(idir))   = cons(ix^D, mom(idir)) * psi6(ix^D)
         end do
      end if ! end adjustments for conserved

    {end do^D&\}

    !check if NaN
    if (any(flag(ixO^S) < 0)) then
       flag(ixO^S) = rho_
       call small_values_error(prim, x, ixI^L, ixO^L, flag, 'gremhd_to_primitive: (NaN z)')
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
       double precision function func(z)
          double precision, intent(in)    :: z
          double precision                :: h_bar
          call get_vars(z, h_out=h_bar)
          func = z - r / h_bar
       end function func
  end subroutine gremhd_to_primitive

end module mod_gremhd_phys_convert
