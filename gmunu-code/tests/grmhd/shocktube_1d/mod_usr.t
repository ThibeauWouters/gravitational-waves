module mod_usr
  use mod_physics
  use mod_grmhd
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

    call set_coordinate_system("Cartesian_1.5D")

    call grmhd_activate()
    call eos_idealgas_activate()

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

    integer  ::  ix^D

    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0

    w(ixO^S,W_vel(1))   = 0.0d0 
    w(ixO^S,W_vel(2))   = 0.0d0 

    select case (iprob)

    case (0)
       ! normal SRHD shocktube problem
       where ( x(ixO^S,1) < 0.5d0 )
          w(ixO^S,rho_)       = 1.0d1
          w(ixO^S,press_)     = 4.0d1 / 3.0d0
       else where
          w(ixO^S,rho_)       = 1.0d0
          w(ixO^S,press_)     = 1.0d-6
       end where

    case (1)
       ! SRMHD shocktube problem
       where ( x(ixO^S,1) < 0.0d0 )
          w(ixO^S,rho_)       = 1.0d0
          w(ixO^S,eps_)       = 1.0d0
          w(ixO^S,press_)     = 1.0d0
          w(ixO^S,Bvec(1))    = 0.5d0
          w(ixO^S,Bvec(2))    = 1.0d0
       else where
          w(ixO^S,rho_)       = 1.25d-1
          w(ixO^S,eps_)       = 0.8d0
          w(ixO^S,press_)     = 1.0d-1
          w(ixO^S,Bvec(1))    = 0.5d0
          w(ixO^S,Bvec(2))    = -1.0d0
       end where

    case default
       error stop "Unknown iprob"
    end select


    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_pressure_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_cs2_one_grid(w(ix^D, cs2_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

  end subroutine rm_init_one_grid

end module mod_usr

!@ARTICLE{2001ApJS..132...83B,
!       author = {{Balsara}, Dinshaw},
!        title = "{Total Variation Diminishing Scheme for Relativistic
!Magnetohydrodynamics}",
!      journal = {\apjs},
!     keywords = {Hydrodynamics, Magnetohydrodynamics: MHD, Relativity},
!         year = 2001,
!        month = jan,
!       volume = {132},
!       number = {1},
!        pages = {83-101},
!          doi = {10.1086/318941},
!       adsurl = {https://ui.adsabs.harvard.edu/abs/2001ApJS..132...83B},
!      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
!}
