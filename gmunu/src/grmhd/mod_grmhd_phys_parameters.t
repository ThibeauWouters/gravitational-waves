module mod_grmhd_phys_parameters
  use mod_global_parameters, only: std_len
  use mod_physics
  implicit none
  public

  !-------------------------------------------------------------------!
  ! Parameters for global settings
  !-------------------------------------------------------------------!
  logical                      :: evolve_hydro = .True.
  logical                      :: evolve_EM = .True.
  logical                      :: use_GR = .false.

  !-------------------------------------------------------------------!
  ! Parameters for con2prim
  !-------------------------------------------------------------------!
  !> tolerance for the root finding
  double precision             :: tolerance = 1.0d-15
  !> maximum iteration for the root finding
  integer                      :: iter_max = 50000
  !> maximum Lorentz factor
  double precision             :: lfac_max = 1.1d1
  double precision             :: v_max

  !-------------------------------------------------------------------!
  ! Parameters for magnetic field
  !-------------------------------------------------------------------!
  !> minimum allowed square amplitude of the magnetic field
  double precision             :: B2_min = 0.0d0

  !-------------------------------------------------------------------!
  ! Parameters for divergence handle
  !-------------------------------------------------------------------!
  logical                      :: divb_4thorder = .False.
  ! DivB cleaning methods
  integer, parameter           :: none           = 0
  integer, parameter           :: divb_multigrid = 1
  integer, parameter           :: divb_glm       = 2
  integer, parameter           :: divb_ct        = 3

  integer                      :: type_divb = 0

  !-------------------------------------------------------------------!
  ! Parameters for divergence with Multigrid
  !-------------------------------------------------------------------!
  logical                      :: elliptic_cleaning = .False.
  ! N smooths will be applied for 1: cycling up, 2: cycling down
  integer, dimension(1:2)      :: divB_mg_n_cycle = (/5,5/)
  logical                      :: divB_mg_redblack = .True.
  double precision             :: divB_mg_tol = 1.0d-8
  integer                      :: divB_mg_it_max = 50

  !-------------------------------------------------------------------!
  ! Parameters for divergence with GLM
  !-------------------------------------------------------------------!
  double precision             :: divB_glm_kappa = 5.0d0

  !-------------------------------------------------------------------!
  ! Parameters for divergence with constrant transport
  !-------------------------------------------------------------------!
  integer                      :: ct_type     = 2
  integer, parameter           :: ct_ave      = 0
  integer, parameter           :: uct_hll     = 1
  integer, parameter           :: uct_contact = 2

  integer                      :: face_to_center_order = 2


  contains

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grmhd_update_eos_one_point(prim)
    use mod_global_parameters
    use mod_eos
    double precision, intent(inout)           :: prim(:)

    ! rho and eps are given, update the rest of the primitive variables
    call eos_get_pressure_one_grid(prim(press_),prim( rho_),prim( eps_))
    call eos_get_cs2_one_grid(prim(cs2_),prim( rho_),prim( eps_))
    ! strictly require cs2 is physical
    if ( (prim(cs2_) > 1.0d0).or.(prim(cs2_) < 0.0d0) ) then
       prim(cs2_) = 0.0d0
    end if
  end subroutine grmhd_update_eos_one_point

  subroutine grmhd_get_lfac2(ixI^L,ixO^L,ps_in,lfac2)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixI^L, ixO^L
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: lfac2(ixI^S)
    integer                         :: idir, ix^D
    double precision                :: gamma(ixI^S,1:3,1:3)

    associate(prim=>ps_in%prim,x=>ps_in%x)
    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
    end do
    
    lfac2 = 1.0d0
    ! calculate W^2 = 1 + W^2 * v^2
    do idir = 1, ndir
       lfac2(ixO^S) = lfac2(ixO^S) + gamma(ixO^S,idir,idir)*prim(ixO^S, W_vel(idir))**2
    end do
    end associate
  end subroutine grmhd_get_lfac2

  !> get some useful variables from primitive
  subroutine grmhd_get_intermediate_variables(ixI^L, ixO^L, prim, x, &
             gamma, lfac2, lfac, v2, v_hat, b2, B_dot_v, Bvec2, bmu, Ptot, htot )
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)                     :: ixI^L, ixO^L
    double precision, intent(in)            :: prim(ixI^S, 1:nprim)
    double precision, intent(in)            :: x(ixI^S, 1:ndim)

    double precision, intent(out), optional :: B_dot_v(ixI^S)   
    double precision, intent(out), optional :: v_hat(ixI^S,1:ndir)   
    double precision, intent(out), optional :: bmu(ixI^S,0:3)   ! projection of B^mu along fluid four-velocity u^nu
    double precision, intent(out), optional :: b2(ixI^S)        ! 
    double precision, intent(out), optional :: Bvec2(ixI^S)        ! 
    double precision, intent(out), optional :: gamma(ixI^S,1:3,1:3)
    double precision, intent(out), optional :: v2(ixI^S)        ! 
    double precision, intent(out), optional :: lfac2(ixI^S) ! Lorentz factor square
    double precision, intent(out), optional :: lfac(ixI^S) ! Lorentz factor
    double precision, intent(out), optional :: htot(ixI^S) ! modified enthalpy:(h + b2/rho)
    double precision, intent(out), optional :: Ptot(ixI^S) ! total pressure

    integer                                 :: idir
    double precision                        :: B_dot_v_tmp(ixI^S)   
    double precision                        :: b2_tmp(ixI^S)        ! 
    double precision                        :: W2v2(ixI^S)        ! 
    double precision                        :: v_hat_tmp(ixI^S,1:ndir) 
    double precision                        :: lfac2_tmp(ixI^S) ! Lorentz factor square
    double precision                        :: lfac_tmp(ixI^S) ! Lorentz factor
    double precision                        :: gamma_tmp(ixI^S,1:3,1:3)

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
    lfac_tmp(ixO^S) = dsqrt( lfac2_tmp(ixO^S) )
    if ( present(lfac2) ) lfac2(ixO^S) = lfac2_tmp(ixO^S)
    if ( present(lfac) ) lfac(ixO^S) = lfac_tmp(ixO^S)
    if ( present(v2) ) v2(ixO^S) = W2v2(ixO^S) / lfac2_tmp(ixO^S)

    ! v_hat^i = v^i * alp - beta
    do idir = 1, ndir
       v_hat_tmp(ixO^S, idir) = prim(ixO^S, alp_) * prim(ixO^S, W_vel(idir)) / lfac_tmp(ixO^S) - prim(ixO^S, beta(idir))
    end do
    if ( present(v_hat) ) v_hat(ixO^S,1:ndir) = v_hat_tmp(ixO^S,1:ndir)

    B_dot_v_tmp(ixO^S) = 0.0d0
    do idir = 1, ndir
       B_dot_v_tmp(ixO^S) = B_dot_v_tmp(ixO^S) & 
              + prim(ixO^S, Bvec(idir)) * gamma_tmp(ixO^S,idir,idir) & 
               * prim(ixO^S, W_vel(idir)) / lfac_tmp(ixO^S) 
    end do
    if ( present(B_dot_v) ) then
       B_dot_v(ixO^S) = B_dot_v_tmp(ixO^S)
    end if

    b2_tmp(ixO^S) = 0.0d0
    do idir = 1, ndir
       b2_tmp(ixO^S) = b2_tmp(ixO^S) & 
              + prim(ixO^S, Bvec(idir))**2 * gamma_tmp(ixO^S,idir,idir)
    end do
    if ( present(Bvec2) ) then
       Bvec2(ixO^S) = b2_tmp(ixO^S)
    end if
    b2_tmp(ixO^S) = b2_tmp(ixO^S) / lfac2_tmp(ixO^S) + B_dot_v_tmp(ixO^S)**2
    if ( present(b2) ) then
       b2(ixO^S) = b2_tmp(ixO^S)
    end if

    if ( present(bmu) ) then
       ! Calculate the projection of B^mu along fluid four-velocity u^nu
       bmu(ixO^S,0) = B_dot_v_tmp(ixO^S) * lfac_tmp(ixO^S) / prim(ixO^S, alp_)
       do idir = 1, ndir
          bmu(ixO^S,idir) = prim(ixO^S, Bvec(idir)) / lfac_tmp(ixO^S) &
                     + bmu(ixO^S,0) * v_hat_tmp(ixO^S, idir)
       end do
       if (ndir>ndim) bmu(ixO^S,ndir+1:3) = 0.0d0
    end if

    if ( present(Ptot) ) then
       ! Calculate the total pressure
       Ptot(ixO^S) = prim(ixO^S, press_) + 0.5d0 * b2_tmp(ixO^S)
    end if

    if ( present(htot) ) then
       ! Calculate the magnetically modified specific enthalpy
       htot(ixO^S) = 1.0d0 + prim(ixO^S, eps_) & 
                + ( prim(ixO^S, press_) + b2_tmp(ixO^S) ) / prim(ixO^S, rho_) 
    end if
  end subroutine grmhd_get_intermediate_variables

  subroutine grmhd_get_lambda(ixI^L, ixO^L, idim, prim, x, lambda, v_hat)
     use mod_global_parameters
     integer, intent(in)             :: ixI^L, ixO^L, idim
     double precision, intent(in)    :: prim(ixI^S, 1:nprim)
     double precision, intent(in)    :: x(ixI^S, 1:ndim)
     double precision, intent(out)   :: lambda(ixI^S, 1:2)
     double precision, intent(out), optional   :: v_hat(ixI^S, 1:ndir)

     double precision                :: gamma(ixI^S,1:3,1:3)
     double precision                :: htot(ixI^S)    ! modified enthalpy:(h + b2/rho)
     double precision                :: cs2(ixI^S), ca2(ixI^S)
     double precision                :: v2(ixI^S)
     double precision                :: lfac(ixI^S)
     double precision                :: b2(ixI^S)
     double precision                :: a2(ixI^S) ! upper bound for the fast wave speed
     double precision                :: tmp1(ixI^S), tmp2(ixI^S), tmp3(ixI^S)
     double precision                :: vel(ixI^S)
     double precision                :: clight(ixI^S)
     double precision                :: inv_gamma_ii(ixI^S)
     integer                         :: idir

     if (present (v_hat)) then
        ! v_hat is needed when in CT schemes
        call grmhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                   gamma=gamma(ixI^S,1:3,1:3), &
                   v2=v2(ixI^S), lfac=lfac(ixI^S), v_hat=v_hat(ixI^S, 1:ndir), &
                   b2=b2(ixI^S), &
                   htot=htot(ixI^S) )
     else
        call grmhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                   gamma=gamma(ixI^S,1:3,1:3), &
                   v2=v2(ixI^S), lfac=lfac(ixI^S), &
                   b2=b2(ixI^S), &
                   htot=htot(ixI^S) )
     end if

     ! beware of coordinate singularities
     where ( dabs(gamma(ixO^S, idim, idim)) > smalldouble )
        ! sound speed square
        cs2(ixO^S) = prim(ixO^S, cs2_)
        ! Calculate Alfven speed
        ca2(ixO^S) =  b2(ixO^S) / ( prim(ixO^S, rho_) * htot(ixO^S) )
        vel(ixO^S) = prim(ixO^S, W_vel(idim)) / lfac(ixO^S)
        inv_gamma_ii(ixO^S) = 1.0d0 / gamma(ixO^S, idim, idim)
     else where
        cs2(ixO^S) = 0.0d0
        ca2(ixO^S) = 0.0d0
        vel(ixO^S) = 0.0d0
        inv_gamma_ii(ixO^S) = 0.0d0
        v2(ixO^S) = 0.0d0
        b2(ixO^S) = 0.0d0
     end where

     ! upper bound for the fast wave speed
     a2(ixO^S) = cs2(ixO^S) + ca2(ixO^S) - ca2(ixO^S) * cs2(ixO^S)
 
     tmp1(ixO^S) = ( 1.0d0 - a2(ixO^S) ) * vel(ixO^S) 
     tmp2(ixO^S) = ( 1.0d0 - v2(ixO^S) * a2(ixO^S) ) 
     tmp3(ixO^S) = ( a2(ixO^S) * ( 1.0d0 - v2(ixO^S) ) * &
                     ( tmp2(ixO^S) * inv_gamma_ii(ixO^S) &
                       - ( 1.0d0 - a2(ixO^S) ) * vel(ixO^S)**2 ) )

     where ( tmp2(ixO^S)==0.0d0 .or. tmp3(ixO^S)<0.0d0 )
        lambda(ixO^S,1) = 0.0d0
        lambda(ixO^S,2) = 0.0d0
     else where
        tmp3(ixO^S) = dsqrt( tmp3(ixO^S) )
        lambda(ixO^S,1) = ( tmp1(ixO^S) + tmp3(ixO^S) ) / tmp2(ixO^S)
        lambda(ixO^S,2) = ( tmp1(ixO^S) - tmp3(ixO^S) ) / tmp2(ixO^S)
     end where

     ! limit with speed of light
     clight(ixO^S) = dsqrt( 1.0d0 / gamma(ixO^S, idim, idim) )
     lambda(ixO^S,1) = max( min( lambda(ixO^S,1), clight(ixO^S) ), -clight(ixO^S) )
     lambda(ixO^S,2) = max( min( lambda(ixO^S,2), clight(ixO^S) ), -clight(ixO^S) )

     ! when using GLM, we need to include two additional modes with speed of
     ! light, to without violating causality
     if (type_divb == divb_glm) then
        lambda(ixO^S,1) = clight(ixO^S)
        lambda(ixO^S,2) = -clight(ixO^S)
     end if

     lambda(ixO^S,1) = prim(ixO^S, alp_) * lambda(ixO^S,1) - prim(ixO^S, beta(idim))
     lambda(ixO^S,2) = prim(ixO^S, alp_) * lambda(ixO^S,2) - prim(ixO^S, beta(idim))
  end subroutine grmhd_get_lambda

  !> Calculate div B within ixO
  subroutine grmhd_get_divb(ixI^L,ixO^L,s,divb)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    type(state), intent(in)         :: s
    double precision, intent(out)   :: divb(ixI^S)

    double precision                :: bvector(ixI^S,1:ndir)
    integer                         :: ixC^L, idir

    associate(prim=>s%prim)

    if ( stagger_grid ) then
      ! note: when use stagger_grid, we use block pointer directly
      divb=0.d0
      do idir=1,ndim
        ixC^L=ixO^L-kr(idir,^D);
        divb(ixO^S)=divb(ixO^S)+s%conss(ixO^S,idir)*s%surfaceC(ixO^S,idir)-&
                                s%conss(ixC^S,idir)*s%surfaceC(ixC^S,idir)
      end do
      divb(ixO^S)=divb(ixO^S)/s%dvolume(ixO^S)
    else
      do idir = 1, ndir
         bvector(ixI^S,idir) = prim(ixI^S, psi_)**6 * prim(ixI^S, Bvec(idir))
      end do
      select case(typediv)
      case("central")
        call divvector(bvector(ixI^S,1:ndir),ixI^L,ixO^L,divb(ixI^S),divb_4thorder)
      case("limited") ! fixme: need to test
        call divvectorS(bvector(ixI^S,1:ndir),ixI^L,ixO^L,divb)
      end select
    end if
    end associate
  end subroutine grmhd_get_divb

  !> This subroutine fix the abnormal values in primitive variables !
  subroutine grmhd_handle_small_values(prim, x, ixI^L, ixO^L, update_eos, subname)
    use mod_global_parameters
    use mod_small_values
    use mod_eos
    use mod_geometry
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: prim(ixI^S,1:nprim)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    logical, intent(in)             :: update_eos

    integer                         :: idir
    integer                         :: ix^D
    double precision                :: eps_min, eps_max
    logical                         :: need_to_fix_eos

    ! avoid coordinate singularities
    ! fixme: this might depends on different bc, but in general this should work
    if ( coordinate /= cartesian ) then
       where ( dabs(x(ixO^S,1)) < smalldouble ) 
          prim(ixO^S, W_vel(1)) = 0.0d0
          prim(ixO^S, beta(1)) = 0.0d0
       end where
    end if

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
             call grmhd_update_eos_one_point(prim(ix^D,1:nprim))
          end if
       {end do^D&\}
    case default
       ! nothing to do here
       !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
       return
    end select
  end subroutine grmhd_handle_small_values

end module mod_grmhd_phys_parameters
