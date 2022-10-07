module mod_grhd_phys_flux
  use mod_physics
  use mod_grhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grhd_phys_flux_init

contains
  !> Initialize the module
  subroutine grhd_phys_flux_init()
    phys_get_cbounds         => grhd_get_cbounds
    phys_get_flux            => grhd_get_flux
  end subroutine grhd_phys_flux_init

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine grhd_get_cbounds(consL, consR, primL, primR, x, ixI^L, ixO^L, idim, cmax, cmin)
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

    call grhd_get_lambda(ixI^L, ixO^L, idim, primL(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambdaL(ixI^S,1:2))
    call grhd_get_lambda(ixI^L, ixO^L, idim, primR(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambdaR(ixI^S,1:2))
  
    tmp1(ixO^S)=max(0.0d0, lambdaL(ixO^S,1), lambdaR(ixO^S,1) )
    tmp2(ixO^S)=min(0.0d0, lambdaL(ixO^S,2), lambdaR(ixO^S,2) )

    if(present(cmin)) then
      cmax(ixO^S) = tmp1(ixO^S)
      cmin(ixO^S) = tmp2(ixO^S)
    else
      cmax(ixO^S) = max(abs(tmp1(ixO^S)), abs(tmp2(ixO^S)))
    end if

  end subroutine grhd_get_cbounds

  ! Calculate flux f_idim[iw]
  subroutine grhd_get_flux(cons, prim, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in)    :: cons(ixI^S, 1:ncons)
    ! primitive w
    double precision, intent(in)    :: prim(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(out)   :: f(ixI^S, 1:ncons)

    double precision                :: v_hat(ixI^S) 
    double precision                :: lfac(ixI^S) 
    double precision                :: alppsi6(ixI^S) 
    integer                         :: idir

    if ( .not. evolve_hydro ) then 
       f(ixO^S, D_) = 0.0d0
       f(ixO^S, mom(1:ndir)) = 0.0d0
       f(ixO^S, tau_) = 0.0d0
       return
    end if

    call grhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                lfac=lfac(ixI^S) )

    ! v_hat = alp * v  - beta
    v_hat(ixO^S) = prim(ixO^S,W_vel(idim)) / lfac(ixO^S) * prim(ixO^S,alp_) - prim(ixO^S,beta(idim))

    ! alp * psi**6
    alppsi6(ixO^S) = prim(ixO^S,alp_) * prim(ixO^S,psi_)**6

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

  end subroutine grhd_get_flux

end module mod_grhd_phys_flux
