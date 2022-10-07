module mod_grrhd_phys_parameters
  use mod_physics
  implicit none
  public

  !-------------------------------------------------------------------!
  ! Parameters for global settings
  !-------------------------------------------------------------------!
  logical                      :: evolve_hydro = .True.
  logical                      :: use_process_adv_global = .false.
  logical                      :: use_GR = .false.
  logical                      :: use_radiation = .false.

  !-------------------------------------------------------------------!
  ! Parameters for con2prim
  !-------------------------------------------------------------------!
  !> tolerance for the root finding
  double precision             :: tolerance = 1.0d-15
  !> maximum iteration for the root finding
  integer                      :: iter_max = 50000
  !> maximum Lorentz factor
  double precision             :: lfac_max = 11.0d0
  double precision             :: v_max, k_max

  !-------------------------------------------------------------------!
  ! Parameters for radiation transport
  !-------------------------------------------------------------------!
  !> tolerance for the root finding
  double precision             :: tol_rad = 1.0d-15
  !> maximum iteration for the root finding
  integer                      :: iter_rad_max = 50000

  !-------------------------------------------------------------------!

  contains

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grrhd_update_eos_one_point(prim)
    use mod_global_parameters
    use mod_eos
    double precision, intent(inout)           :: prim(:)

    ! rho and eps are given, update the rest of the primitive variables
    call eos_get_pressure_one_grid(prim(press_),prim( rho_),prim( eps_))
    call eos_get_cs2_one_grid(prim(cs2_),prim( rho_),prim( eps_))
    ! strictly require cs2 is physical
    prim(cs2_) = max( min( 1.0d0, prim(cs2_) ) , 0.0d0)
  end subroutine grrhd_update_eos_one_point

  subroutine grrhd_get_lambda(ixI^L, ixO^L, idim, prim, x, lambda)
     use mod_global_parameters
     integer, intent(in)             :: ixI^L, ixO^L, idim
     double precision, intent(in)    :: prim(ixI^S, 1:nprim)
     double precision, intent(in)    :: x(ixI^S, 1:ndim)
     double precision, intent(out)   :: lambda(ixI^S, 1:2)

     integer                         :: idir
     integer                         :: ix^D
     double precision                :: tmp1(ixI^S), tmp2(ixI^S)
     double precision                :: gamma(ixI^S,1:3,1:3)
     double precision                :: cs2(ixI^S), v2(ixI^S)
     double precision                :: lfac(ixI^S), lfac2(ixI^S)
     double precision                :: pi(ixI^S)
     double precision                :: chi(ixI^S) ! closure function
     double precision                :: coeff_thin(ixI^S), coeff_thick(ixI^S)
     double precision                :: lambda_thin(ixI^S, 1:2), lambda_thick(ixI^S, 1:2)
     double precision                :: vel(ixI^S, 1:ndir)

     call grrhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                   gamma=gamma(ixI^S,1:3,1:3), v2=v2(ixI^S), lfac2=lfac2(ixI^S), lfac=lfac(ixI^S))
     if (.not.use_radiation) then
        ! normal GRHD lambda
        ! sound speed square
        cs2(ixO^S) = prim(ixO^S, cs2_)
   
        do idir = 1, ndir
           vel(ixO^S,idir) = prim(ixO^S, W_vel(idir)) / lfac(ixO^S)
        end do
   
        tmp1(ixO^S) = vel(ixO^S, idim) * ( 1.0d0 - cs2(ixO^S) )
        tmp2(ixO^S) = dsqrt( cs2(ixO^S) * ( 1.0d0 - v2(ixO^S) ) * &
                       ( ( 1.0d0 - v2(ixO^S) * cs2(ixO^S) ) / gamma(ixO^S, idim, idim) &
                       - vel(ixO^S, idim)**2 * ( 1.0d0 - cs2(ixO^S) ) ) )
   
        lambda(ixO^S,1) = ( tmp1(ixO^S) + tmp2(ixO^S) ) / ( 1.0d0 - v2(ixO^S) * cs2(ixO^S) )
        lambda(ixO^S,2) = ( tmp1(ixO^S) - tmp2(ixO^S) ) / ( 1.0d0 - v2(ixO^S) * cs2(ixO^S) )
   
        ! limit with speed of light, here we use tmp1 to store it
        tmp1(ixO^S) = dsqrt( 1.0d0 / gamma(ixO^S, idim, idim) )
        lambda(ixO^S,1) = max( min( lambda(ixO^S,1), tmp1(ixO^S) ), -tmp1(ixO^S) )
        lambda(ixO^S,2) = max( min( lambda(ixO^S,2), tmp1(ixO^S) ), -tmp1(ixO^S) )
   
        !lambda(ixO^S,0) = prim(ixO^S, alp_) * prim(ixO^S, vel(idim)) - prim(ixO^S, beta(idim))  
        lambda(ixO^S,1) = prim(ixO^S, alp_) * lambda(ixO^S,1) - prim(ixO^S, beta(idim)) 
        lambda(ixO^S,2) = prim(ixO^S, alp_) * lambda(ixO^S,2) - prim(ixO^S, beta(idim)) 
     else
        ! get F -> F * E
        do idir = 1, ndir
           vel(ixO^S, idir) = prim(ixO^S, nu_F_over_E(idir)) * prim(ixO^S, nu_E) 
        end do
        ! get sqrt F
        tmp2 = 0.0d0
        do idir = 1, ndir
           tmp2(ixO^S) = tmp2(ixO^S) + gamma(ixO^S,idir,idir) * vel(ixO^S, idir)**2
        end do
        tmp2(ixO^S) = dsqrt( tmp2(ixO^S) )
    
        ! fixme: need to up index
        tmp1(ixO^S) = dabs( vel(ixO^S, idim) ) / tmp2(ixO^S)
    
        lambda_thin(ixO^S,1) = - prim(ixO^S, beta(idim)) + prim(ixO^S, alp_) * tmp1(ixO^S) 
        lambda_thin(ixO^S,2) = - prim(ixO^S, beta(idim)) - prim(ixO^S, alp_) * tmp1(ixO^S) 

        ! now work out lambda_thick
        pi(ixO^S) = prim(ixO^S, alp_) * prim(ixO^S, W_vel(idim)) / lfac2(ixO^S)
        tmp1(ixO^S) = 2.0d0 * prim(ixO^S, alp_) * prim(ixO^S, W_vel(idim))
        tmp2(ixO^S) = dsqrt( prim(ixO^S, alp_) * ( 2.0d0 * lfac2(ixO^S) + 1.0d0 ) / gamma(ixO^S, idim, idim) &
                       - 2.0d0 * lfac2(ixO^S) * pi(ixO^S)**2 ) 

        lambda_thick(ixO^S,1) = - prim(ixO^S, beta(idim)) &
                                 + ( tmp1(ixO^S) + tmp2(ixO^S) ) / ( 2.0d0 * lfac2(ixO^S) + 1.0d0 )
        lambda_thick(ixO^S,2) = - prim(ixO^S, beta(idim)) &
                                 + ( tmp1(ixO^S) - tmp2(ixO^S) ) / ( 2.0d0 * lfac2(ixO^S) + 1.0d0 )  

        lambda_thick(ixO^S,1) = min( - prim(ixO^S, beta(idim)) + pi(ixO^S), lambda_thick(ixO^S,1) )
        lambda_thick(ixO^S,2) = min( - prim(ixO^S, beta(idim)) + pi(ixO^S), lambda_thick(ixO^S,2) )
    
        ! now combine two lambdas
        ! loop over points
        {do ix^D = ixO^LIM^D \}
           chi(ix^D) = minerbo_closure_func( prim(ix^D, nu_zeta) )
           coeff_thin(ix^D) = c_thin(chi(ix^D))
           coeff_thick(ix^D) = 1.0d0 - coeff_thin(ix^D)
        {end do^D&\} ! end loop for every points
        lambda(ixO^S,1) = coeff_thin(ixO^S) * lambda_thin(ixO^S,1) &
                         + coeff_thick(ixO^S) * lambda_thick(ixO^S,1)
        lambda(ixO^S,2) = coeff_thin(ixO^S) * lambda_thin(ixO^S,2) &
                         + coeff_thick(ixO^S) * lambda_thick(ixO^S,2)
     end if
  end subroutine grrhd_get_lambda

  !> This subroutine fix the abnormal values in conserved variables !
  subroutine grrhd_handle_small_values(prim, x, ixI^L, ixO^L, update_eos, subname)
    use mod_global_parameters
    use mod_small_values
    use mod_eos
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: prim(ixI^S,1:nprim)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    logical, intent(in)             :: update_eos

    integer                      :: idir
    integer                      :: ix^D
    double precision             :: eps_min, eps_max
    logical                      :: need_to_fix_eos

    select case (small_values_method)
    case ("replace")
       ! check the prim variables one by one
       {do ix^D = ixO^LIM^D \}
          need_to_fix_eos = .False.
          if ( prim(ix^D, rho_) < small_rho_thr ) then
             ! atmosphere handling
             prim(ix^D, rho_) = small_rho
             prim(ix^D, eps_) = small_eps
             prim(ix^D, W_vel(:)) = 0.0d0
             need_to_fix_eos = .True.
          else
             ! this is not atmosphere
             call eos_get_eps_range( prim(ix^D, rho_), eps_min, eps_max)
             if ( ( prim(ix^D, eps_) < eps_min ) .or. ( prim(ix^D, eps_) > eps_max ) ) then
                prim(ix^D, eps_) = max( min( eps_max, prim(ix^D, eps_) ), eps_min )
                need_to_fix_eos = .True.
             end if
          end if
          if ( need_to_fix_eos .and. update_eos ) then
             call grrhd_update_eos_one_point(prim(ix^D,1:nprim))
          end if
       {enddo^D&\}
       !where ( dabs(prim(ixO^S, W_vel(ndir))) <= smalldouble )
       !   prim(ixO^S, W_vel(ndir)) = 0.0d0
       !end where
    case default
       ! nothing to do here
       !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
       return
    end select
  end subroutine grrhd_handle_small_values

  !> get some useful variables from primitive
  subroutine grrhd_get_intermediate_variables(ixI^L, ixO^L, prim, x, &
             Pij, &
             gamma, lfac2, lfac, v2, h )
    use mod_global_parameters
    use mod_geometry
    use mod_rootfinding
    integer, intent(in)                     :: ixI^L, ixO^L
    double precision, intent(in)            :: prim(ixI^S, 1:nprim)
    double precision, intent(in)            :: x(ixI^S, 1:ndim)

    double precision, intent(out), optional :: gamma(ixI^S,1:3,1:3)
    double precision, intent(out), optional :: v2(ixI^S)        ! 
    double precision, intent(out), optional :: lfac2(ixI^S) ! Lorentz factor square
    double precision, intent(out), optional :: lfac(ixI^S) ! Lorentz factor
    double precision, intent(out), optional :: h(ixI^S) ! enthalpy: h
    double precision, intent(out), optional :: Pij(ixI^S,1:3,1:3)

    integer                                 :: idir, jdir
    double precision                        :: W2v2(ixI^S), v2_tmp(ixI^S)        ! 
    !double precision                        :: v_bar(ixI^S,1:3) 
    double precision                        :: lfac2_tmp(ixI^S) ! Lorentz factor square
    double precision                        :: lfac_tmp(ixI^S) ! Lorentz factor
    double precision                        :: gamma_tmp(ixI^S,1:3,1:3)
    double precision                        :: E(ixI^S), F2(ixI^S), Fi(ixI^S,1:ndir) ! F_i
    double precision                        :: v_dot_F(ixI^S)  !  v^i F_i
    double precision                        :: gamma_H(ixI^S,1:ndir)
    integer                                 :: ix^D 
    ! --------------   variables needed to get Pij  ------------- 
    double precision                        :: chi(ixI^S) ! closure function
    double precision                        :: coeff_thin(ixI^S), coeff_thick(ixI^S)
    double precision                        :: J_over_3(ixI^S)

    if ( present(h) ) then
       ! Calculate the magnetically modified specific enthalpy
       h(ixO^S) = 1.0d0 + prim(ixO^S, eps_) & 
                + prim(ixO^S, press_) / prim(ixO^S, rho_) 
    end if

    if ( .not. ( present(gamma) .or. &
                 present(lfac2) .or. present(lfac) .or. present(v2) .or. &
                 present(Pij) ) ) then
       return
    end if

    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_tmp(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma_tmp(ixO^S,idir,idir) = gamma_tmp(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
    end do
    if ( present(gamma) ) then
       gamma(ixO^S,1:3,1:3) = gamma_tmp(ixO^S,1:3,1:3) 
    end if

    W2v2 = 0.0d0
    ! calculate W^2 * v^2 first
    do idir = 1, ndir
       W2v2(ixO^S) = W2v2(ixO^S) + gamma_tmp(ixO^S,idir,idir)*prim(ixO^S, W_vel(idir))**2
    end do
    ! Calculate the Lorentz factor from velocity
    lfac2_tmp(ixO^S) = W2v2(ixO^S) + 1.0d0 
    if ( present(lfac2) ) lfac2(ixO^S) = lfac2_tmp(ixO^S)
    lfac_tmp(ixO^S) = dsqrt( lfac2_tmp(ixO^S) )
    if ( present(lfac) ) lfac(ixO^S) = lfac_tmp(ixO^S)
    v2_tmp(ixO^S) = W2v2(ixO^S) / lfac2_tmp(ixO^S)
    if ( present(v2) ) v2(ixO^S) = v2_tmp(ixO^S)

    if ( .not. present(Pij)  ) then
       return
    end if

    ! some useful relations
    E(ixO^S) = prim(ixO^S, nu_E)
    ! Note that the momentum density in prim is F_i / E
    do idir = 1, ndir
       Fi(ixO^S, idir) = prim(ixO^S, nu_F_over_E(idir)) * E(ixO^S)
    end do
    ! calculate the F^i F_i
    F2 = 0.0d0
    do idir = 1, ndir
       F2(ixO^S) = F2(ixO^S) + Fi(ixO^S, idir)**2 / gamma(ixO^S, idir, idir)
    end do
    ! v^i * F_i
    v_dot_F = 0.0d0
    do idir = 1, ndir
       v_dot_F(ixO^S) = v_dot_F(ixO^S) + prim(ixO^S, W_vel(idir)) / lfac_tmp(ixO^S) * Fi(ixO^S, idir)
    end do

    if ( present(Pij) ) then
       ! loop over points
       {do ix^D = ixO^LIM^D \}
          chi(ix^D) = minerbo_closure_func( prim(ix^D, nu_zeta) )
          coeff_thin(ix^D) = c_thin(chi(ix^D))
          coeff_thick(ix^D) = 1.0d0 - coeff_thin(ix^D)
       {end do^D&\} ! end loop for every points

       ! J and H can be computed here if needed

       ! optically thin part of pressure tensor
       do idir = 1, ndir; do jdir = 1, ndir
          Pij(ixO^S,idir,jdir) = coeff_thin(ixO^S) * E(ixO^S)/F2(ixO^S) * Fi(ixO^S,idir) * Fi(ixO^S,jdir)
       end do; end do

       ! optically thick part of pressure tensor
       J_over_3(ixO^S) = 1.0d0 / (2.0d0 * lfac2_tmp(ixO^S) + 1.0d0) &
                  * ( (2.0d0 * lfac2_tmp(ixO^S) - 1.0d0) * E(ixO^S) - 2.0d0 * lfac2_tmp(ixO^S) * v_dot_F(ixO^S) )
       do idir = 1, ndir
          gamma_H(ixO^S,idir) = Fi(ixO^S,idir)/gamma(ixO^S,idir,idir)/lfac_tmp(ixO^S) &
               + prim(ixO^S,W_vel(idir)) / (2.0d0 * lfac2_tmp(ixO^S) + 1.0d0) &
               * ( (4.0d0 * lfac2_tmp(ixO^S) + 1.0d0) * v_dot_F(ixO^S) - 4.0d0 * lfac2_tmp(ixO^S) * E(ixO^S) )
       end do

       do idir = 1, ndir
          Pij(ixO^S,idir,idir) = Pij(ixO^S,idir,idir) + &
           coeff_thick(ixO^S) * ( J_over_3(ixO^S) * &
             ( 4.0d0 * prim(ixO^S,W_vel(idir)) * prim(ixO^S,W_vel(idir)) + 1.0d0 / gamma(ixO^S,idir,idir)) )
       end do

       do idir = 1, ndir; do jdir = 1, ndir
          Pij(ixO^S,idir,jdir) = Pij(ixO^S,idir,jdir) + &
              coeff_thick(ixO^S) * ( lfac_tmp(ixO^S) &
                * (gamma_H(ixO^S,idir) * prim(ixO^S, W_vel(jdir)) / lfac_tmp(ixO^S) &
                  +gamma_H(ixO^S,jdir) * prim(ixO^S, W_vel(idir)) / lfac_tmp(ixO^S) &
                  ) ) 
       end do; end do
    end if
  end subroutine grrhd_get_intermediate_variables

  !> the implementation here follows SpECTRE https://spectre-code.org/
  subroutine grrhd_get_nu_zeta(ixI^L, ixO^L, prim, x, &
             v2_in, lfac_in, lfac2_in, gamma_in )
    use mod_global_parameters
    use mod_geometry
    use mod_rootfinding
    implicit none
    integer, intent(in)                     :: ixI^L, ixO^L
    double precision, intent(inout)         :: prim(ixI^S, 1:nprim)
    double precision, intent(in)            :: x(ixI^S, 1:ndim)

    double precision, intent(in), optional  :: gamma_in(ixI^S,1:3,1:3)
    double precision, intent(in), optional  :: v2_in(ixI^S)        ! 
    double precision, intent(in), optional  :: lfac2_in(ixI^S) ! Lorentz factor square
    double precision, intent(in), optional  :: lfac_in(ixI^S) ! Lorentz factor

    integer                                 :: ix^D 
    integer                                 :: idir, jdir
    double precision                        :: gamma(ixI^S,1:3,1:3)
    double precision                        :: v2(ixI^S), lfac2(ixI^S), lfac(ixI^S) 
    integer                                 :: error_code = -1
    ! --------------   variables needed to solve zeta  ------------- 
    double precision                        :: E, F2, Fi(1:ndir) ! F_i
    double precision                        :: v_dot_F  !  v^i F_i
    double precision                        :: zeta ! variable Eddington factor
    double precision                        :: zeta_m, zeta_p
    double precision                        :: chi ! closure function    
    double precision                        :: coeff_thin, coeff_thick

    ! fluid-frame energy density:
    double precision                        :: J_0, J_thin, J_thick
    ! fluid-frame momentum density:
    double precision                        :: H_0_t, H_thin_t, H_thick_t
    double precision                        :: H_0_v, H_thin_v, H_thick_v
    double precision                        :: H_0_F, H_thin_F, H_thick_F
    ! fluid-frame momentum density square and the relevent variables:
    double precision                        :: H2_0
    double precision                        :: H2_thin, H2_thick
    double precision                        :: H2_thin_thin, H2_thick_thick
    double precision                        :: H2_thin_thick

    if ( .not. ( present(gamma_in) &
          .or. present(lfac_in) .or. present(lfac2_in) .or. present(v2_in) ) ) then
     call grrhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                   gamma=gamma(ixI^S,1:3,1:3), v2=v2(ixI^S), lfac2=lfac2(ixI^S), lfac=lfac(ixI^S))
    else
       gamma(ixI^S,1:3,1:3) = gamma_in(ixI^S,1:3,1:3) 
       lfac(ixI^S) = lfac_in(ixI^S) 
       lfac2(ixI^S) = lfac2_in(ixI^S) 
       v2(ixI^S) = v2_in(ixI^S) 
    end if

    ! loop over points
    {do ix^D = ixO^LIM^D \}
       E = prim(ix^D, nu_E)
       ! Note that the momentum density in prim is F_i / E
       Fi(1:ndir) = prim(ix^D, nu_F_over_E(1:ndir)) * E
   
       ! calculate the F^i F_i
       F2 = 0.0d0
       do idir = 1, ndir
          F2 = F2 + Fi(idir)**2 / gamma(ix^D,idir,idir)
       end do
       ! if the velocity can be ignored, assume that
       ! fluid frame = inertial frame
       if ( v2(ix^D) <= smalldouble ) then
          zeta = min( max( dsqrt( F2 )/ E, 0.0d0), 1.0d0)
      
       ! if the velocity can NOT be ignored, we need to go through a more
       ! expensive closure calculation
       else
          ! initial guess of the closure factor
          zeta = min( max( prim(ix^D, nu_zeta), 0.0d0), 1.0d0)
      
          ! v^i * F_i
          v_dot_F = 0.0d0
          do idir = 1, ndir
             v_dot_F = v_dot_F + prim(ix^D, W_vel(idir)) / lfac(ix^D) * Fi(idir)
          end do
          ! Decomposition of the fluid-frame energy density:
          ! J = J_0 + coeff_thin * J_Thin + coeff_thick * J_Thick
          J_0 = lfac2(ix^D) * ( E - 2.0d0 * v_dot_F )
          J_thin = lfac2(ix^D) * E * v_dot_F**2 / F2 
          J_thick = ( lfac2(ix^D) - 1.0d0 ) / ( 1.0d0 + 2.0d0 * lfac2(ix^D) ) &
               * ( 4.0d0 * lfac2(ix^D) * v_dot_F + E * ( 3.0d0 - 2.0d0 * lfac2(ix^D) ) )
      
          ! Decomposition of the fluid-frame momentum density:
          ! H_a = -( H_0_t + coeff_thick H_thick_t + coeff_thin H_thin_t) t_a
          !  - ( H_0_v + coeff_thick H_thick_v + coeff_thin H_thin_v) v_a
          !  - ( H_0_F + coeff_thick H_thick_F + coeff_thin H_thin_F) F_a
          ! with t_a the unit normal, v_a the 3-velocity, and F_a the
          ! inertial frame momentum density. 
          !
          ! This is a decomposition of convenience, 
          ! which is not unique: F_a and v_a are not
          ! orthogonal vectors, but both are normal to t_a.
          H_0_t = lfac(ix^D) * ( J_0 + v_dot_F - E )
          H_0_v = lfac(ix^D) * J_0 
          H_0_F = - lfac(ix^D)
          H_thin_t = lfac(ix^D) * J_thin
          H_thin_v = H_thin_t
          H_thin_F = lfac(ix^D) * E * v_dot_F / F2
          H_thick_t = lfac(ix^D) * J_thick
          H_thick_v = H_thick_t &
             + lfac(ix^D) / ( 2.0d0 * lfac2(ix^D) + 1.0d0 ) * &
                ( (3.0d0 - 2.0d0 * lfac2(ix^D)) * E - (1.0d0 - 2.0d0 * lfac2(ix^D)) * v_dot_F )
          H_thick_F = lfac2(ix^D) * v2(ix^D)
          ! Quantities needed for the computation of H^2 = H^a H_a,
          ! independent of zeta. We write:
          ! H^2 = H2_0 + H2_thin * coeff_thin + H2_thick*coeff_thick
          ! + H2_thin_thin * coeff_thin^2 + H2_thick_thick * coeff_thick^2
          ! + H2_thin_thick * coeff_thin * coeff_thick;
          H2_0 = - H_0_t**2 + H_0_v**2 * v2(ix^D) + H_0_F**2 * F2 &
                 + 2.0d0 * H_0_v * H_0_F * v_dot_F
          H2_thin = 2.0d0 * &
               (H_0_v * H_thin_v * v2(ix^D) + H_0_F * H_thin_F * F2 + &
                H_0_v * H_thin_F * v_dot_F + H_0_F * H_thin_v * v_dot_F - &
                H_0_t * H_thin_t)
          H2_thick = 2.0d0 * &
               (H_0_v * H_thick_v * v2(ix^D) + H_0_F * H_thick_F * F2 + &
                H_0_v * H_thick_F * v_dot_F + &
                H_0_F * H_thick_v * v_dot_F - H_0_t * H_thick_t)
          H2_thin_thick = 2.0d0 * &
               (H_thin_v * H_thick_v * v2(ix^D) + H_thin_F * H_thick_F * F2 + &
                H_thin_v * H_thick_F * v_dot_F + &
                H_thin_F * H_thick_v * v_dot_F - H_thin_t * H_thick_t)
          H2_thick_thick = (H_thick_v**2) * v2(ix^D) + (H_thick_F**2) * F2 + &
                2.0d0 * H_thick_v * H_thick_F * v_dot_F - (H_thick_t**2) 
          H2_thin_thin = (H_thin_v**2) * v2(ix^D) + (H_thin_F**2) * F2 + &
               2.0d0 * H_thin_v * H_thin_F * v_dot_F - (H_thin_t**2)
      
          ! find zeta though root finding
      
          ! atempt Newton-Raphson method to find the root
          ! To avoid failures in Newton-Raphson solver at the boundary of
          ! the allowed domain for zeta, test the edge values first.
          zeta_m = 0.0d0
          zeta_p = 1.0d0
          if ( dabs(zeta_J2_minus_H2(zeta_m)) <= tol_rad ) then
             zeta = zeta_m
          else if ( dabs(zeta_J2_minus_H2(zeta_p)) <= tol_rad ) then
             zeta = zeta_p
          else
             call rootfinding_constrained_newton_raphson(zeta, zeta_m, zeta_p, &
                   tol_rad, 50, error_code, zeta_J2_minus_H2, dev_zeta_J2_minus_H2)
             ! just in case if newton raphson method fails, use illinois method
             if( error_code/=0 ) then
                write(*,*) "Warning: Newton-Raphson method failed to find zeta, now try illinois method"
                call rootfinding_illinois(zeta, zeta_m, zeta_p, tol_rad, iter_rad_max, error_code, zeta_J2_minus_H2)
             end if
             ! check the root
             select case (error_code)
             case (-1) ! nothing happened
                call mpistop("have you ever attemp to find the root in con2prim?")
             case (1) ! failed to find the root
                call mpistop("cannot find the root in zeta_J2_minus_H2")
             case (2) ! z is NaN
                call mpistop("NaN in zeta_J2_minus_H2")
             case (3) ! z is not bracketed
                call mpistop("the root is not bracketed in con2prim")
             end select
          end if
       end if
       prim(ix^D, nu_zeta) = zeta
    {end do^D&\} ! end loop for every points

    contains
       double precision function zeta_J2_minus_H2(zeta) result(f)
          double precision, intent(in)     :: zeta
          double precision                 :: chi
          double precision                 :: coeff_thin
          double precision                 :: coeff_thick
          double precision                 :: J
          double precision                 :: H2
   
          chi = minerbo_closure_func(zeta)
          coeff_thin = c_thin(chi)
          coeff_thick = 1.0d0 - coeff_thin
   
          J = J_0 + coeff_thin * J_Thin + coeff_thick * J_Thick
   
          H2 = H2_0 + H2_thick * coeff_thick + &
                      H2_thin * coeff_thin + &
                      H2_thin_thin * (coeff_thin**2) + &
                      H2_thick_thick * (coeff_thick**2) + &
                      H2_thin_thick * coeff_thin * coeff_thick
   
          f = ( ( zeta * J )**2 - H2 ) / E**2
       end function zeta_J2_minus_H2
   
       double precision function dev_zeta_J2_minus_H2(zeta) result(f)
          double precision, intent(in)     :: zeta
          double precision                 :: chi, dchi
          double precision                 :: coeff_thin, coeff_thick
          double precision                 :: dcoeff_thin, dcoeff_thick
          double precision                 :: J, H2
          double precision                 :: dJ, dH2, dH2_dcoeff_thin, dH2_dcoeff_thick
   
          chi = minerbo_closure_func(zeta)
          dchi = minerbo_closure_deriv(zeta)
          coeff_thin = c_thin(chi)
          coeff_thick = 1.0d0 - coeff_thin
          dcoeff_thin = dc_thin_dchi(chi) * dchi
          dcoeff_thick = dc_thick_dchi(chi) * dchi
   
          J = J_0 + coeff_thin * J_Thin + coeff_thick * J_Thick
   
          H2 = H2_0 + H2_thick * coeff_thick + &
                      H2_thin * coeff_thin + &
                      H2_thin_thin * (coeff_thin**2) + &
                      H2_thick_thick * (coeff_thick**2) + &
                      H2_thin_thick * coeff_thin * coeff_thick
   
          dJ = dcoeff_thin * J_Thin + dcoeff_thick * J_Thick
   
          dH2_dcoeff_thin = &
                      H2_thin + &
                      H2_thin_thin * (2.0d0 * coeff_thin) + &
                      H2_thin_thick * coeff_thick
   
          dH2_dcoeff_thick = &
                      H2_thick + &
                      H2_thick_thick * (2.0d0 * coeff_thick) + &
                      H2_thin_thick * coeff_thin
   
          dH2 = dH2_dcoeff_thin * dcoeff_thin + dH2_dcoeff_thick * dcoeff_thick
   
          f = ( 2.0d0 * ( zeta * J**2 + zeta**2 * J * dJ) &
                - dH2 ) / E**2
       end function dev_zeta_J2_minus_H2

 end subroutine grrhd_get_nu_zeta

  ! --------------------------
  !  some useful functions
  ! --------------------------
  double precision function minerbo_closure_func(zeta) result(chi)
    double precision, intent(in)     :: zeta
    chi = 1.0d0/3.0d0 &
          + zeta**2 * (0.4d0 - 2.0d0/15.0d0 * zeta + 0.4d0 * zeta**2)
    !if ( (chi > 1.0d0) .or. (chi < 1.0d0/3.0d0) ) stop "wrong chi"
  end function minerbo_closure_func

  double precision function minerbo_closure_deriv(zeta) result(dchi)
    double precision, intent(in)     :: zeta
    dchi = 0.4d0 * zeta * ( 2.0d0 - zeta + 4.0d0 * zeta**2)
  end function minerbo_closure_deriv

  double precision function c_thin(chi) result(coeff)
    double precision, intent(in)     :: chi 
    coeff = 1.5d0 * chi - 0.5d0
  end function c_thin

  double precision function c_thick(chi) result(coeff)
    double precision, intent(in)     :: chi 
    coeff = 1.5d0 * ( 1.0d0 - chi )
  end function c_thick

  double precision function dc_thin_dchi(chi) result(coeff)
    double precision, intent(in)     :: chi 
    coeff = 1.5d0 
  end function dc_thin_dchi

  double precision function dc_thick_dchi(chi) result(coeff)
    double precision, intent(in)     :: chi 
    coeff = - 1.5d0
  end function dc_thick_dchi

end module mod_grrhd_phys_parameters
