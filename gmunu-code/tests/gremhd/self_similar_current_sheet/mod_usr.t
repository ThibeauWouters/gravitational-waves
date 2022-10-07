module mod_usr
  use mod_physics
  use mod_gremhd
  use mod_eos_idealgas

  implicit none
  private

  double precision, parameter :: resistivity = 1.0d-2

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid
    usr_get_resistivity => my_get_resistivity

    call set_coordinate_system("Cartesian_1.75D")

    call gremhd_activate()
    call eos_idealgas_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm_init_one_grid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics
    use mod_gremhd
    use mod_eos

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim) ! primitive variables

    integer  ::  ix^D

    w = 0.0d0

    ! no gauge effects
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(:)) = 0.0d0

    w(ixO^S,rho_)       = 1.0d0
    w(ixO^S,press_)     = 5.0d3
    w(ixO^S,W_vel(:))   = 0.0d0 

    w(ixO^S,Bvec(2))   = erf( 0.5d0 * x(ixO^S,1) / dsqrt(resistivity) )

    w(ixO^S,Evec(3))   = dsqrt(resistivity/dpi) * exp( - 2.5d-1 * x(ixO^S,1)**2 / resistivity )

    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       !call eos_get_pressure_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_cs2_one_grid(w(ix^D, cs2_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

  end subroutine rm_init_one_grid

  ! get resistivity
  subroutine my_get_resistivity(ixI^L,ixO^L,cons,x,eta)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(in) :: cons(ixI^S,1:ncons) 
    double precision, intent(out):: eta(ixI^S) 
    eta(ixO^S) = resistivity
  end subroutine my_get_resistivity


end module mod_usr
