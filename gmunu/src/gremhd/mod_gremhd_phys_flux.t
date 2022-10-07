module mod_gremhd_phys_flux
  use mod_physics
  use mod_gremhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: gremhd_phys_flux_init

contains
  !> Initialize the module
  subroutine gremhd_phys_flux_init()
    phys_get_cbounds         => gremhd_get_cbounds
    phys_get_flux            => gremhd_get_flux
  end subroutine gremhd_phys_flux_init

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine gremhd_get_cbounds(consL, consR, primL, primR, x, ixI^L, ixO^L, idim, cmax, cmin)
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

    double precision, dimension(ixI^S,1:3,1:3) :: gamma_hat
    double precision, dimension(ixI^S,1:2) :: lambdaL, lambdaR
    double precision, dimension(ixI^S,1:2) :: tmp_c

    call gremhd_get_lambda(ixI^L, ixO^L, idim, primL(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambdaL(ixI^S,1:2))
    call gremhd_get_lambda(ixI^L, ixO^L, idim, primR(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambdaR(ixI^S,1:2))
    
    tmp_c(ixO^S,1)=max(0.0d0, lambdaL(ixO^S,1), lambdaR(ixO^S,1) ) 
    tmp_c(ixO^S,2)=min(0.0d0, lambdaL(ixO^S,2), lambdaR(ixO^S,2) )

    if(present(cmin)) then
      cmax(ixO^S) = tmp_c(ixO^S,1)
      cmin(ixO^S) = tmp_c(ixO^S,2)
    else
      cmax(ixO^S) = max(abs(tmp_c(ixO^S,1)), abs(tmp_c(ixO^S,2)))
    end if
  end subroutine gremhd_get_cbounds

  ! Calculate flux f_idim[iw]
  subroutine gremhd_get_flux(cons, prim, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in)    :: cons(ixI^S, 1:ncons)
    ! primitive w
    double precision, intent(in)    :: prim(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(out)   :: f(ixI^S, ncons)

    double precision                :: gamma(ixI^S,1:3,1:3)
    double precision                :: v_hat(ixI^S, 1:ndir) 
    double precision                :: alppsi6(ixI^S) 
    double precision                :: h(ixI^S) ! 
    double precision                :: lfac(ixI^S) ! Lorentz factor
    double precision                :: sqrt_gamma(ixI^S)
    double precision                :: E2(ixI^S), B2(ixI^S) 
    double precision                :: inv_sqrt_gamma(ixI^S), inv_gamma_ii(ixI^S)
    integer                         :: idir, jdir, kdir

    call gremhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), sqrt_gamma=sqrt_gamma(ixI^S), &
                lfac=lfac(ixI^S), v_hat=v_hat(ixI^S,1:ndir), &
                h=h(ixI^S), E2=E2(ixI^S), B2=B2(ixI^S))

    where ( sqrt_gamma(ixO^S) > smalldouble )
       inv_sqrt_gamma(ixO^S) = 1.0d0 / sqrt_gamma(ixO^S)
       inv_gamma_ii(ixO^S) = 1.0d0 / gamma(ixO^S, idim, idim)
    else where
       inv_sqrt_gamma(ixO^S) = 0.0d0
       inv_gamma_ii(ixO^S) = 0.0d0
    end where

    ! alp * psi**6
    alppsi6(ixO^S) = prim(ixO^S,alp_) * prim(ixO^S,psi_)**6

    if ( evolve_hydro ) then
       ! Density flux
       f(ixO^S, D_) = v_hat(ixO^S, idim) * cons(ixO^S, D_)
   
       ! Momentum flux f^idim_idir
       do idir = 1, ndir
         f(ixO^S, mom(idir)) =  - prim(ixO^S, beta(idim)) * cons(ixO^S, mom(idir)) &
                    + alppsi6(ixO^S) * gamma(ixO^S, idir, idir) * ( prim(ixO^S, rho_) * h(ixO^S) * &
                                        prim(ixO^S, W_vel(idim)) * prim(ixO^S, W_vel(idir))&
                                      - prim(ixO^S, Evec(idim)) * prim(ixO^S, Evec(idir)) &
                                      - prim(ixO^S, Bvec(idim)) * prim(ixO^S, Bvec(idir)) ) 
       end do
       f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) &
              + alppsi6(ixO^S) * ( prim(ixO^S, press_) + 0.5d0 * ( E2(ixO^S) + B2(ixO^S)) )
   
       ! Energy flux 
       f(ixO^S, tau_) = - prim(ixO^S, beta(idim)) * cons(ixO^S, tau_) &
                    + prim(ixO^S, alp_) * ( cons(ixO^S, mom(idim)) * inv_gamma_ii(ixO^S) &
                                           - cons(ixO^S, D_) * prim(ixO^S, W_vel(idim)) / lfac(ixO^S) ) 
    else
       f(ixO^S, D_) = 0.0d0
       f(ixO^S, tau_) = 0.0d0
       f(ixO^S, mom(:)) = 0.0d0
    end if

    if ( evolve_EM ) then
       ! flux of E and B field
       do idir = 1, ndir
         if (idir /= idim) then
            f(ixO^S, Econs(idir)) = - prim(ixO^S, beta(idim)) * cons(ixO^S, Econs(idir)) &
                                    + prim(ixO^S, beta(idir)) * cons(ixO^S, Econs(idim))
            f(ixO^S, Bcons(idir)) = - prim(ixO^S, beta(idim)) * cons(ixO^S, Bcons(idir)) &
                                    + prim(ixO^S, beta(idir)) * cons(ixO^S, Bcons(idim))
         else
            f(ixO^S, Econs(idir)) = 0.0d0
            f(ixO^S, Bcons(idir)) = 0.0d0
         end if
       end do

       do idir = 1, ndir
          do jdir = 1, ndir
             f(ixO^S, Econs(idir)) = f(ixO^S, Econs(idir)) &
                        + prim(ixO^S, alp_) * inv_sqrt_gamma(ixO^S) * lvc(idim, idir, jdir) * gamma(ixO^S, jdir, jdir) * cons(ixO^S, Bcons(jdir))
             f(ixO^S, Bcons(idir)) = f(ixO^S, Bcons(idir)) &
                        - prim(ixO^S, alp_) * inv_sqrt_gamma(ixO^S) * lvc(idim, idir, jdir) * gamma(ixO^S, jdir, jdir) * cons(ixO^S, Econs(jdir))
          end do
       end do

       if (type_divb == glm) then
         f(ixO^S, Bphi_cons_) = prim(ixO^S,alp_) * cons(ixO^S, Bcons(idim)) &
                           - prim(ixO^S,beta(idim)) * cons(ixO^S, Bphi_cons_)
       end if
       if (type_dive == glm) then
         f(ixO^S, Ephi_cons_) = prim(ixO^S,alp_) * cons(ixO^S, Econs(idim)) &
                           - prim(ixO^S,beta(idim)) * cons(ixO^S, Ephi_cons_)
       end if
    else
       f(ixO^S, Bcons(:)) = 0.0d0
       f(ixO^S, Econs(:)) = 0.0d0
       if (type_divb == glm) &
          f(ixO^S, Bphi_cons_) = 0.0d0
       if (type_dive == glm) &
          f(ixO^S, Ephi_cons_) = 0.0d0
    end if

    ! beware of sqrt_gamma = 0 cases, which means it hits the coordinate singularities!
    ! fixme: maybe no needed
    where( sqrt_gamma(ixO^S) < smalldouble )
       f(ixO^S, D_) = 0.0d0
       f(ixO^S, tau_) = 0.0d0
       f(ixO^S, mom(idim)) = 0.0d0
       f(ixO^S, Econs(idim)) = 0.0d0
       f(ixO^S, Bcons(idim)) = 0.0d0
    end where
  end subroutine gremhd_get_flux

end module mod_gremhd_phys_flux
