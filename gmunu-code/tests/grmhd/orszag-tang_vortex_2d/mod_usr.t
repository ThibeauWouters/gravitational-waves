module mod_usr
  use mod_physics
  use mod_grmhd
  use mod_eos_idealgas

  implicit none
  private

  ! aux output variables
  integer, protected  :: divB_ = -1

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid
    usr_before_main_loop => rm_before_main_loop
    usr_process_adv_grid => cal_div_B

    call set_coordinate_system("Cartesian_2D")

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
    double precision                :: v(ixI^S,1:ndir)
    double precision                :: lfac(ixI^S)

    integer  ::  ix^D, idir, max_loc(1:ndir)

    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0

    v(ixO^S, 1)   = -3.0d0 / (4.0d0 * dsqrt(2.0d0)) * dsin( x(ixO^S,2) )
    v(ixO^S, 2)   =  3.0d0 / (4.0d0 * dsqrt(2.0d0)) * dsin( x(ixO^S,1) )

    lfac(ixO^S) = v(ixO^S,1)**2 + v(ixO^S,2)**2
    lfac(ixO^S) = 1.0d0 / dsqrt( 1.0d0 - lfac(ixO^S) )
    ! so that it wont less than 1
    where ( lfac(ixO^S) < 1.0d0 )
    !   lfac(ixO^S) = 1.0d0
    end where

    w(ixO^S,W_vel(1))   = v(ixO^S, 1) * lfac(ixO^S)
    w(ixO^S,W_vel(2))   = v(ixO^S, 2) * lfac(ixO^S)

    w(ixO^S,rho_)       = 1.0d0
    w(ixO^S,press_)     = 1.0d0
    w(ixO^S,Bvec(1))    = - sin( x(ixO^S,2) )
    w(ixO^S,Bvec(2))    =   sin( 2.0d0 * x(ixO^S,1) )

    do idir = 1, ndir
       where ( dabs( w(ixO^S,Bvec(idir)) ) <= smalldouble ) 
       !   w(ixO^S,Bvec(idir))    = smalldouble
       end where
    end do

    w(ixO^S, divB_) = 0.0d0

    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_pressure_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_cs2_one_grid(w(ix^D, cs2_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

  end subroutine rm_init_one_grid

  subroutine rm_before_main_loop()
    use mod_grmhd_phys_parameters
    use mod_grmhd_phys_divb_mg

    if (type_divb == divb_multigrid) call grmhd_clean_divb_multigrid()
  end subroutine rm_before_main_loop

  !> calculate div B after the advance 
  subroutine cal_div_B(igrid,level,ixI^L,ixO^L,qt,prim,x)
    use mod_global_parameters
    use mod_grmhd_phys_parameters

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: prim(ixI^S,1:nprim)

    double precision                :: divb(ixI^S)

    call grmhd_get_divb(ixI^L, ixO^L, prim(ixI^S,1:nprim), divb(ixI^S), divb_4thorder )
    prim(ixO^S, divB_) = dabs( divb(ixO^S) )

  end subroutine cal_div_B

end module mod_usr
