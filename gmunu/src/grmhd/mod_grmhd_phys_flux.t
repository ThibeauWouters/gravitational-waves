module mod_grmhd_phys_flux
  use mod_physics
  use mod_grmhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grmhd_phys_flux_init

contains
  !> Initialize the module
  subroutine grmhd_phys_flux_init()
    phys_get_cbounds         => grmhd_get_cbounds
    phys_get_flux            => grmhd_get_flux
  end subroutine grmhd_phys_flux_init

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine grmhd_get_cbounds(consL, consR, primL, primR, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative left and right status
    double precision, intent(in)    :: consL(ixI^S, 1:ncons), consR(ixI^S, 1:ncons)
    ! primitive left and right status
    double precision, intent(in)    :: primL(ixI^S, 1:nprim), primR(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision, dimension(ixI^S,1:2) :: lambdaL, lambdaR
    double precision, dimension(ixI^S,1:2) :: tmp_c

     call grmhd_get_lambda(ixI^L, ixO^L, idim, primL(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambdaL(ixI^S,1:2))
     call grmhd_get_lambda(ixI^L, ixO^L, idim, primR(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambdaR(ixI^S,1:2))

    tmp_c(ixO^S,1)=max(0.0d0, lambdaL(ixO^S,1), lambdaR(ixO^S,1) )
    tmp_c(ixO^S,2)=min(0.0d0, lambdaL(ixO^S,2), lambdaR(ixO^S,2) ) 

    if(present(cmin)) then
      cmax(ixO^S) = tmp_c(ixO^S,1)
      cmin(ixO^S) = tmp_c(ixO^S,2)
    else
      cmax(ixO^S) = max(abs(tmp_c(ixO^S,1)), abs(tmp_c(ixO^S,2)))
    end if

  end subroutine grmhd_get_cbounds

  ! Calculate flux f_idim[iw]
  subroutine grmhd_get_flux(cons, prim, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

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
    double precision                :: B_dot_v(ixI^S)
    double precision                :: Ptot(ixI^S) ! total pressure
    double precision                :: lfac(ixI^S) ! Lorentz factor
    double precision                :: B_over_W(ixI^S, 1:ndir) 
    integer                         :: idir

    call grmhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                lfac=lfac(ixI^S), &
                v_hat=v_hat(ixI^S,1:ndir), &
                B_dot_v=B_dot_v(ixI^S), &
                Ptot=Ptot(ixI^S) )

    ! alp * psi**6
    alppsi6(ixO^S) = prim(ixO^S,alp_) * prim(ixO^S,psi_)**6

    ! B^i / W
    do idir = 1, ndir
       B_over_W(ixO^S, idir) = prim(ixO^S, Bvec(idir)) / lfac(ixO^S)
    end do

    if ( evolve_hydro ) then 
       ! Density flux
       f(ixO^S, D_) = v_hat(ixO^S, idim) * cons(ixO^S, D_)
   
       ! Momentum flux f^idim_idir
       do idir = 1, ndir
         f(ixO^S, mom(idir)) = v_hat(ixO^S, idim) * cons(ixO^S, mom(idir)) &
                        - alppsi6(ixO^S) * ( B_over_W(ixO^S, idir) + B_dot_v(ixO^S) * prim(ixO^S, W_vel(idir)) ) &
                          * gamma(ixO^S, idir, idir) * B_over_W(ixO^S, idim)
       end do
       f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) + alppsi6(ixO^S) * Ptot(ixO^S) 
   
       ! Energy flux 
       f(ixO^S, tau_) = v_hat(ixO^S, idim) * cons(ixO^S, tau_) &
                         + alppsi6(ixO^S) * ( Ptot(ixO^S) * prim(ixO^S, W_vel(idim)) / lfac(ixO^S) &
                                - B_dot_v(ixO^S) * prim(ixO^S, Bvec(idim)) )
    else
       f(ixO^S, D_) = 0.0d0
       f(ixO^S, mom(1:ndir)) = 0.0d0
       f(ixO^S, tau_) = 0.0d0
    end if

    if ( evolve_EM ) then 
       ! flux of B field
       do idir = 1, ndir
         f(ixO^S, Bcons(idir)) = v_hat(ixO^S, idim) * cons(ixO^S, Bcons(idir)) &
                               - v_hat(ixO^S, idir) * cons(ixO^S, Bcons(idim))
       end do
   
       ! flux terms when using GLM 
       if (type_divb == divb_glm) then
         f(ixO^S, Bphi_cons_) = prim(ixO^S,alp_) * cons(ixO^S, Bcons(idim)) &
                           - prim(ixO^S,beta(idim)) * cons(ixO^S, Bphi_cons_)
       end if
    else
       f(ixO^S, Bcons(1:ndir)) = 0.0d0
       if (type_divb == divb_glm) &
          f(ixO^S, Bphi_cons_) = 0.0d0
    end if

  end subroutine grmhd_get_flux

end module mod_grmhd_phys_flux
