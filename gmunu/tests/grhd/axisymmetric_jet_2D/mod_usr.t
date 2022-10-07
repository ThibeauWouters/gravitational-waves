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
    usr_special_bc    => specialbound_usr

    call set_coordinate_system("cylindrical_2D")

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

    integer  ::  ix^D, idir

    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0

    {^IFONED call mpistop("This is a multi-D HD problem") }

    w(ixO^S,press_)     = 5.0d-3
    w(ixO^S,rho_)   = 1.0d2
    v(ixO^S,1) = 0.0d0
    v(ixO^S,2) = 0.0d0

    where (x(ixO^S,1)<1.0d0 .and. x(ixO^S,2)<1.0d0)
       w(ixO^S,rho_)   = 1.0d0
       v(ixO^S,2) = 0.995d0
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

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)
    use mod_eos

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nprim)
    double precision                :: v(ixG^S,1:ndir)
    double precision                :: lfac(ixG^S)
    integer :: ixI^L, ix^D

    ixImin^DD=ixBmin^DD;
    ixImax^DD=ixBmin^D-1+nghostcells^D%ixImax^DD=ixBmax^DD;
    ! Outflow:
    do ix2=ixImin2,ixImax2
       w(ixImin1:ixImax1,ix2,rho_) = w(ixImin1:ixImax1,ixImax2+1,rho_) 
       w(ixImin1:ixImax1,ix2,press_)   = w(ixImin1:ixImax1,ixImax2+1,press_) 
       w(ixImin1:ixImax1,ix2,W_vel(1))  = w(ixImin1:ixImax1,ixImax2+1,W_vel(1))
       w(ixImin1:ixImax1,ix2,W_vel(2))  = w(ixImin1:ixImax1,ixImax2+1,W_vel(2))
    end do

    where(dabs(x(ixG^S,1))<1.0d0)
       w(ixG^S,rho_)       = 1.0d0
       w(ixG^S,press_)     = 5.0d-3

       v(ixG^S,1) = 0.0d0
       v(ixG^S,2) = 0.995d0
       lfac(ixG^S) = v(ixG^S,1)**2 + v(ixG^S,2)**2
       lfac(ixG^S) = 1.0d0 / dsqrt( 1.0d0 - lfac(ixG^S) )
   
       w(ixG^S,W_vel(1))   = v(ixG^S, 1) * lfac(ixG^S)
       w(ixG^S,W_vel(2))   = v(ixG^S, 2) * lfac(ixG^S)
    else where
       ! Reflective:
       !   w(ixI^S,rho_) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,rho_) 
       !   w(ixI^S,e_) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,e_) 
       !   w(ixI^S,W_vel(1)) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,W_vel(1))
       !   w(ixI^S,W_vel(2)) =-w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,W_vel(2))
    end where

    ! dont forget to update these primitive variables
    {do ix^D = ixG^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_pressure_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_cs2_one_grid(w(ix^D, cs2_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

  end subroutine specialbound_usr

end module mod_usr
