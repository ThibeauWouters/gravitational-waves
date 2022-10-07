!> General relativistic hydrodynamics PPM module
module mod_grhd_ppm
  use mod_physics

  implicit none
  private

  ! Public methods
  public :: grhd_ppm_init

contains

  !> Initialize the module
  subroutine grhd_ppm_init()
    use mod_physics_ppm

    phys_ppm_flatcd => grhd_ppm_flatcd
    phys_ppm_flatsh => grhd_ppm_flatsh
  end subroutine grhd_ppm_init

  subroutine grhd_ppm_flatcd(rec_from, rec_to, ixI^L, ixO^L, ixL^L, ixR^L, w, d2w, drho, dp)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)             :: rec_from, rec_to, ixI^L, ixO^L, ixL^L, ixR^L
    double precision, intent(in)    :: w(ixI^S, 1:nprim),d2w(ixG^T,rec_from:rec_to)
    double precision, intent(inout) :: drho(ixG^T), dp(ixG^T)

    double precision                :: gamma_eos_tmp(ixI^S), d2w_press(ixI^S)
    integer                         :: ix^D 

    ! cs2 needed to be updated based on the updated prim vars
    call phys_get_csound2(ixI^L, ixO^L,w(ixI^S, 1:nprim), gamma_eos_tmp(ixO^S))
    gamma_eos_tmp(ixO^S) = gamma_eos_tmp(ixO^S) * w(ixO^S,rho_)/w(ixO^S,press_)& 
                           * (1.0d0 + w(ixO^S,press_)/w(ixO^S,rho_) + w(ixO^S,eps_))
    {do ix^D = ixO^LIM^D \}
    ! pressure is needed to be updated based on the d2w
       call eos_get_pressure_one_grid(d2w_press(ix^D), d2w(ix^D, rho_), d2w(ix^D, eps_))
    {enddo^D&\}

    drho(ixO^S) = gamma_eos_tmp(ixO^S) * abs(d2w(ixO^S,rho_)) &
                  / min( w(ixL^S, rho_), w(ixR^S, rho_) )

    dp(ixO^S) = gamma_eos_tmp(ixO^S) * abs(d2w(ixO^S,press_)) &
                  / min( w(ixL^S, press_), w(ixR^S, press_) )
  end subroutine grhd_ppm_flatcd

  subroutine grhd_ppm_flatsh(ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L, idims, w, beta, z, dv)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S, 1:nprim)
    double precision, intent(inout) :: beta(ixG^T), z(ixG^T), dv(ixG^T)

    double precision                :: cs2(ixI^S)
    integer                         :: jxO^L, hxO^L, ix^L

    ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
    where ( abs( w(ixRR^S, press_) - w(ixLL^S, press_) )> smalldouble )
       beta(ixO^S) = abs( w(ixR^S, press_) - w(ixL^S, press_) ) &
                     /( w(ixRR^S, press_) - w(ixLL^S, press_) )
    else where
       beta(ixO^S) = 0.0d0
    end where

    ! eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
    call phys_get_csound2(ixI^L, ixO^L, w(ixI^S, 1:nprim), cs2(ixO^S))
    z(ixO^S) = abs( w(ixR^S, press_) - w(ixL^S, press_) ) &
                     /( w(ixO^S, rho_)*cs2(ixO^S) )

    hxO^L=ixO^L-kr(idims,^D);
    jxO^L=ixO^L+kr(idims,^D);
    ix^L=ixO^L^LADD1;

    !dv(ixO^S)=(w(jxO^S,veloc(idims))-w(hxO^S,veloc(idims)) )
    dv(ixO^S)=(w(jxO^S,W_vel(idims))-w(hxO^S,W_vel(idims)) )

  end subroutine grhd_ppm_flatsh

end module mod_grhd_ppm
