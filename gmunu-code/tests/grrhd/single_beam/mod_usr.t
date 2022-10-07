module mod_usr
  use mod_physics
  use mod_grrhd
  use mod_eos_polytrope

  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid => my_init_one_grid
    usr_internal_bc    => internalbound_usr

    call set_coordinate_system("Cartesian_2.5D")
    call grrhd_activate()
    ! use atmosphere
    call eos_atmo_activate()
    call eos_polytrope_activate()

    !call initialize_gmunu()
    ! setup atmosphere density from the profile
      call eos_initialize_atmo(1.0d0)
  end subroutine usr_init

  ! Initialize one grid
  subroutine my_init_one_grid(ixG^L,ix^L,w,x)
    use mod_global_parameters
    use mod_eos

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nprim)

    integer  ::  idir
    integer  ::  n_theta_offset, n2
    double precision :: delta_theta

    {^IFONED call mpistop("This is a multi-D Rad problem") }

    !w = 1.0d-15
    w = 0.0d0
    
    ! density cannot be zero
    w(ix^S, rho_) = 1.0d-15

    ! flat space initial condition
    w(ix^S, alp_) = 1.0d0
    w(ix^S, psi_) = 1.0d0
    w(ix^S, beta(1:3)) = 0.0d0

    where ( (x(ix^S, 1) <= - 0.4d0) .and. &
            (abs(x(ix^S, 2)) <= 0.12d0) )
       w(ix^S, nu_E)      = 1.0d0
       w(ix^S, nu_F_over_E(1)) = 0.999999d0
    else where
       w(ix^S, nu_E)      = 1.0d-10
       w(ix^S, nu_F_over_E(1)) = 0.999999d-10
    end where
    w(ix^S, nu_F_over_E(2))      = 0.0d0
    w(ix^S, nu_F_over_E(3))      = 0.0d0

    ! in our notation, prim Fi = Fi/Ei
    do idir=1, ndir
       w(ix^S, nu_F_over_E(idir)) = w(ix^S, nu_F_over_E(idir)) / w(ix^S, nu_E)
    end do

  end subroutine my_init_one_grid

  subroutine internalbound_usr(level,qt,ixI^L,ixO^L,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixI^L, ixO^L, level
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim)

    where ( (x(ixI^S, 1) <= - 0.4d0) .and. &
            (abs(x(ixI^S, 2)) <= 0.12d0) )
       w(ixI^S, nu_E)           = 1.0d0
       w(ixI^S, nu_F_over_E(1))      = 0.9999999d0
       w(ixI^S, nu_F_over_E(2))      = 0.0d0
       w(ixI^S, nu_F_over_E(3))      = 0.0d0
    else where ( (x(ixI^S, 1) <= - 0.4d0) .and. &
            (abs(x(ixI^S, 2)) > 0.12d0) )
       w(ixI^S, nu_E)           = 1.0d-10
       w(ixI^S, nu_F_over_E(1))      = 0.9999999d0
       w(ixI^S, nu_F_over_E(2))      = 0.0d0
       w(ixI^S, nu_F_over_E(3))      = 0.0d0
    end where

  end subroutine internalbound_usr

end module mod_usr
