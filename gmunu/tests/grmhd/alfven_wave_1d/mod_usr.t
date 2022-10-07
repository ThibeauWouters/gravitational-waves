module mod_usr
  use mod_physics
  use mod_grmhd
  use mod_eos_idealgas

  implicit none
  private

  ! aux output variables
  integer, protected  :: divB_

  double precision, parameter :: rho_0 = 1.0d0
  double precision, parameter :: P_0 = 0.5d0 ! this can be 0.1 to 1, we pick 0.5 so that va = 0.5
  double precision, parameter :: A_0 = 1.0d0
  double precision, parameter :: B_0 = 1.0d0

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid

    call set_coordinate_system("Cartesian_1.75D")
    call grmhd_activate()
    call eos_idealgas_activate()

    divB_ = var_set_auxvar('divB')
  end subroutine usr_init

  ! Initialize one grid
  subroutine rm_init_one_grid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics
    use mod_grmhd
    use mod_eos

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim) ! primitive variables

    double precision :: rho_h(ixI^S) ! rho(1+eps) + P
    double precision :: va2(ixI^S) ! Alfven speed
    double precision :: tmp(ixI^S)
    double precision :: k ! wave vector 2pi / L_x
    integer  ::  ix^D

    k = 2.0d0 * dpi / dabs( xprobmax1 - xprobmin1 )
    !k = 1.0d0
    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0
    w(ixO^S,beta(3)) = 0.0d0

    w(ixO^S,rho_) = rho_0
    w(ixO^S,press_) = P_0
    ! get eps
    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

    rho_h(ixO^S) = w(ixO^S,rho_) * ( 1.0d0 + w(ixO^S,eps_) ) + w(ixO^S,press_)

    select case (iprob)

    case (0)
       ! small amplitude cp alfven wave test
       w(ixO^S,Bvec(1))    = B_0
       w(ixO^S,Bvec(2))    = A_0 * w(ixO^S,Bvec(1)) * dcos( k * x(ixO^S,1) )
       w(ixO^S,Bvec(3))    = A_0 * w(ixO^S,Bvec(1)) * dsin( k * x(ixO^S,1) )

    case default
       error stop "Unknown iprob"
    end select

    ! get Alfven speed
    tmp(ixO^S) = 2.0d0 * w(ixO^S,Bvec(1))**2 &
               / ( rho_h(ixO^S) + w(ixO^S,Bvec(1))**2 * (1.0d0 + A_0**2))
    va2(ixO^S) = tmp(ixO^S) / ( 1.0d0 + dsqrt( 1.0d0 - A_0**2 * tmp(ixO^S)**2 ) )

    ! get W
    tmp(ixO^S) = 1.0d0 / dsqrt(1.0d0 - va2(ixO^S) * A_0**2)
    ! get Wv
    w(ixO^S,W_vel(1))   = 0.0d0
    w(ixO^S,W_vel(2))   = - tmp(ixO^S) * dsqrt(va2(ixO^S)) * A_0 * dcos( k * x(ixO^S,1) )  
    w(ixO^S,W_vel(3))   = - tmp(ixO^S) * dsqrt(va2(ixO^S)) * A_0 * dsin( k * x(ixO^S,1) )  

    ! update rest of the prim
    call phys_update_eos(ixI^L, ixO^L, w)

  end subroutine rm_init_one_grid

end module mod_usr
