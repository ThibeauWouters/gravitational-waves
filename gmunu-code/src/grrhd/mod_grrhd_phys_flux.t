module mod_grrhd_phys_flux
  use mod_physics
  use mod_grrhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grrhd_phys_flux_init

contains
  !> Initialize the module
  subroutine grrhd_phys_flux_init()
    phys_get_cbounds         => grrhd_get_cbounds
    phys_get_flux            => grrhd_get_flux
  end subroutine grrhd_phys_flux_init

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine grrhd_get_cbounds(consL, consR, primL, primR, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative left and right status
    double precision, intent(in)    :: consL(ixI^S, 1:ncons), consR(ixI^S, 1:ncons)
    ! primitive left and right status
    double precision, intent(in)    :: primL(ixI^S, 1:nprim), primR(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision, dimension(ixI^S) :: tmp1,tmp2
    double precision, dimension(ixI^S,1:2) :: lambdaL, lambdaR

    call grrhd_get_lambda(ixI^L, ixO^L, idim, primL(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambdaL(ixI^S,1:2))
    call grrhd_get_lambda(ixI^L, ixO^L, idim, primR(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambdaR(ixI^S,1:2))
  
    tmp1(ixO^S)=max(0.0d0, lambdaL(ixO^S,1), lambdaR(ixO^S,1) )
    tmp2(ixO^S)=min(0.0d0, lambdaL(ixO^S,2), lambdaR(ixO^S,2) )

    if(present(cmin)) then
      cmax(ixO^S) = tmp1(ixO^S)
      cmin(ixO^S) = tmp2(ixO^S)
    else
      cmax(ixO^S) = max(abs(tmp1(ixO^S)), abs(tmp2(ixO^S)))
    end if

  end subroutine grrhd_get_cbounds

  ! Calculate flux f_idim[iw]
  subroutine grrhd_get_flux(cons, prim, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in)    :: cons(ixI^S, 1:ncons)
    ! primitive w
    double precision, intent(in)    :: prim(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(out)   :: f(ixI^S, 1:ncons)

    double precision                :: gamma(ixI^S,1:3,1:3)
    double precision                :: v_hat(ixI^S) 
    double precision                :: lfac(ixI^S) 
    double precision                :: alppsi6(ixI^S) 
    integer                         :: idir, jdir, kdir
    double precision                :: Pij(ixI^S,1:3,1:3)
    double precision                :: P_ij(ixI^S,1:3,1:3)

    call grrhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                   gamma=gamma(ixI^S,1:3,1:3),  lfac=lfac(ixI^S), &
                   Pij=P_ij(ixI^S,1:3,1:3)  )

    ! v_hat = v * alp - beta
    v_hat(ixO^S) = prim(ixO^S,W_vel(idim)) / lfac(ixO^S) * prim(ixO^S,alp_) - prim(ixO^S,beta(idim))

    ! alp * psi**6
    alppsi6(ixO^S) = prim(ixO^S,alp_) * prim(ixO^S,psi_)**6

    if ( evolve_hydro ) then 
       ! Density flux
       f(ixO^S, D_) = v_hat(ixO^S) * cons(ixO^S, D_)
   
       ! Momentum flux
       do idir = 1, ndir
         f(ixO^S, mom(idir)) = v_hat(ixO^S) * cons(ixO^S, mom(idir))
       end do
       f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) + alppsi6(ixO^S) * prim(ixO^S, press_) 
   
       ! Energy flux 
       f(ixO^S, tau_) = v_hat(ixO^S) * cons(ixO^S, tau_) &
                         + alppsi6(ixO^S) * prim(ixO^S, press_) * prim(ixO^S,W_vel(idim)) / lfac(ixO^S)
    else
       f(ixO^S, D_) = 0.0d0
       f(ixO^S, mom(1:ndir)) = 0.0d0
       f(ixO^S, tau_) = 0.0d0
    end if

    if ( use_radiation ) then 
       ! first we get P^{ij}, here we denote it as P_ij
       !call grrhd_get_Pij(ixI^L, ixO^L, cons, prim, x, P_ij)
       ! here we need P^j_i, here we denote it as Pij
       Pij = 0.0d0
       do idir = 1, ndir; do jdir = 1, ndir
          do kdir = 1, ndir
             Pij(ixO^S,jdir,idir) = Pij(ixO^S,jdir,idir) &
                        + gamma(ixO^S,kdir,idir) * P_ij(ixO^S, jdir, kdir)
          end do
       end do; end do
   
       ! Momentum flux
       do idir = 1, ndir
         f(ixO^S, nu_Fcons(idir)) = &
                     prim(ixO^S, alp_) * Pij(ixO^S, idim, idir) &
                   - prim(ixO^S, beta(idim)) * cons(ixO^S, nu_Fcons(idir))
       end do
   
       ! Energy flux 
       f(ixO^S, nu_Econs) = prim(ixO^S,alp_) * cons(ixO^S, nu_Fcons(idim)) &
                                      / gamma(ixO^S, idim, idim)  &
                           - prim(ixO^S, beta(idim)) * cons(ixO^S,nu_Econs) 
    end if

  end subroutine grrhd_get_flux

  subroutine grrhd_modify_flux(ixI^L, ixO^L, consL, consR, primL, primR, xi, s, idir, f)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(inout) :: consL(ixI^S,1:ncons), consR(ixI^S,1:ncons)
    double precision, intent(inout) :: primL(ixI^S,1:nprim), primR(ixI^S,1:nprim)
    double precision, intent(in)    :: xi(ixI^S, 1:^ND)
    type(state)                     :: s
    double precision, intent(inout) :: f(ixI^S, ncons)

    double precision                :: weight_factor(ixI^S, 1:2)
    double precision                :: opacities(ixI^S, 1:3) ! they are eta, and two kappas
    double precision                :: kappa_total(ixI^S) ! total opacity

    !call usr_get_resistivity(ixI^L, ixO^L, cons(ixI^S, 1:ncons), x(ixI^S, 1:ndim), eta(ixI^S))

  end subroutine grrhd_modify_flux

end module mod_grrhd_phys_flux
