module mod_usr
  use mod_physics
  use mod_grhd
  use mod_eos_idealgas

  implicit none
  private

  double precision, parameter :: A = 2.0d-1
  double precision, parameter :: v0 = 2.0d-1

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid

    call set_coordinate_system("Cartesian_1D")
    call grhd_activate()
    call eos_idealgas_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm_init_one_grid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics
    use mod_grhd
    use mod_eos

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim) ! primitive variables
    integer  ::  ix^D

    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0

    w(ixO^S,rho_) = 1.0d0 + A * dsin ( 2.0d0 * dpi &
               * ( x(ixO^S,1) ) )
    w(ixO^S,press_) = 1.0d0
    ! get eps
    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

    ! get Wv
    w(ixO^S,W_vel(1))   = v0 / dsqrt(1.0d0 - v0**2)

    ! update rest of the prim
    call phys_update_eos(ixI^L, ixO^L, w)

  end subroutine rm_init_one_grid

end module mod_usr
