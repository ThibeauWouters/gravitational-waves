module mod_usr
  use mod_physics
  use mod_grhd
  use mod_eos_idealgas

  implicit none
  private

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid

    call set_coordinate_system("Cartesian_2D")

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
    double precision                :: v(ixI^S,1:ndir)
    double precision                :: lfac(ixI^S)

    integer  ::  ix^D, idir, max_loc(1:ndir)

    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0

    {^IFONED call mpistop("This is a multi-D HD problem") }

    ! the left top  quadrant
    where (x(ixO^S,1)<0.0d0 .and. x(ixO^S,2)>0.0d0)
       w(ixO^S,press_)     = 1.0d0
       w(ixO^S,rho_)   = 0.1d0
       v(ixO^S,1) = 0.99d0
       v(ixO^S,2) = 0.0d0
    end where

    ! the right top  quadrant
    where (x(ixO^S,1)>0.0d0 .and. x(ixO^S,2)>0.d0)
       w(ixO^S,press_)     = 2.762987d-3
       w(ixO^S,rho_)   = 5.477875d-3
       v(ixO^S,1) = 0.0d0
       v(ixO^S,2) = 0.0d0
    end where

    ! the left bottom quadrant
    where (x(ixO^S,1)<0.0d0 .and. x(ixO^S,2)<0.0d0)
       w(ixO^S,press_)     = 1.0d0
       w(ixO^S,rho_)   = 0.5d0
       v(ixO^S,1) = 0.0d0
       v(ixO^S,2) = 0.0d0
    end where

    ! the right bottom quadrant
    where (x(ixO^S,1)>0.0d0 .and. x(ixO^S,2)<0.0d0)
       w(ixO^S,press_)     = 1.0d0
       w(ixO^S,rho_)   = 0.1d0
       v(ixO^S,1) = 0.0d0
       v(ixO^S,2) = 0.99d0
    end where

    lfac(ixO^S) = v(ixO^S,1)**2 + v(ixO^S,2)**2
    lfac(ixO^S) = 1.0d0 / dsqrt( 1.0d0 - lfac(ixO^S) )

    w(ixO^S,W_vel(1))   = v(ixO^S, 1) * lfac(ixO^S)
    w(ixO^S,W_vel(2))   = v(ixO^S, 2) * lfac(ixO^S)

    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_pressure_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_cs2_one_grid(w(ix^D, cs2_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

  end subroutine rm_init_one_grid

end module mod_usr
