!> Module for eos
module mod_eos
  use mod_global_parameters, only: name_len, tinydouble, smalldouble, bigdouble

  implicit none
  public

  !> describing the eos type of the simulation
  integer                              :: eos_type = -1

  !> EOS types
  integer, parameter                   :: vacuum    = 0
  integer, parameter                   :: polytrope = 1
  integer, parameter                   :: idealgas  = 2
  integer, parameter                   :: hybrid    = 3
  integer, parameter                   :: tabulated = 4
  integer, parameter                   :: piecewise = 5
  integer, parameter                   :: hqs       = 6

  !> The smallest allowed density
  double precision, public             :: small_rho_fac = 0.1d0
  double precision, public             :: small_rho_thr = 0.0d0
  double precision, public             :: small_rho     = smalldouble
  !> The smallest allowed eps
  double precision, public             :: small_eps     = smalldouble 
  !> The smallest allowed press
  double precision, public             :: small_press   = smalldouble 

  !> The smallest allowed conserved density D
  double precision, public             :: small_D   = smalldouble
  !> The smallest allowed conserved variables tau ~ eps * rho
  double precision, public             :: small_tau = smalldouble

  !-------------------------------------------------------------------!
  ! About small values & atmosphere treatments
  !-------------------------------------------------------------------!
  ! When atmosphere treactment is used, small_rho, small_eps and 
  ! small_press are the density, eps and pressure of the atmosphere
  ! here.
  !
  ! IF YOU WANT TO USE ANOTHER SMALL VALUE TREATMENTS HERE, PLS WRITE DOWN WT A SUMMARY OF TREATMENT YOU USED HERE!!!!
  !
  ! rho < small_rho_thr 	      => all atmosphere outputs
  ! eps < eos_epsmin or < small_eps   => set it to small_eps value to avoid problems (generally set smapp_eps >0.0d0)
  ! ye  > eos_yemax      	      => set it to eos_yemax
  ! temp< eos_tempmin  		      => set it to eos_tempmin 

  !> parameters for atmosphere
  integer                              :: atmo_eos_type = -1
  double precision                     :: atmo_gamma    = 2.0d0
  double precision                     :: atmo_adiab    = 1.0d2
  double precision                     :: atmo_cs2      = 0.0d0
  double precision                     :: atmo_ye       = 0.5d0
  double precision                     :: atmo_temp     = smalldouble
  logical                              :: atmo_beta_eqm = .False.

  !> min-max values for tabulated eos:
  double precision, public             :: eos_rhomin  = smalldouble
  double precision, public             :: eos_rhomax  = bigdouble
  double precision, public             :: eos_yemin   = 0.0d0
  double precision, public             :: eos_yemax   = 1.0d0
  double precision, public             :: eos_tempmin = smalldouble
  double precision, public             :: eos_tempmax = bigdouble
  double precision, public             :: eos_epsmin  = 0.0d0
  double precision, public             :: eos_epsmax  = bigdouble
  double precision, public             :: eos_hmin    = 1.0d0   ! used in Kcon2prim

  ! constants
  double precision, parameter :: ggrav              = 6.673d-8            ! gravitational constant, cm^3 g^(-1) s^(-2)
  double precision, parameter :: clight             = 2.99792458d10       ! speed of light, cm s^(-1)

  double precision, parameter :: rho_gf             = 1.61930347d-18
  double precision, parameter :: press_gf           = 1.80171810d-39
  double precision, parameter :: eps_gf             = 1.11265006d-21
  double precision, parameter :: time_gf            = 2.03001708d+05
  double precision, parameter :: mass_gf            = 5.02765209d-34
  double precision, parameter :: length_gf          = 6.77140812d-06
  double precision, parameter :: energy_gf          = 5.59424238d-55
  double precision, parameter :: lum_gf             = 2.7556091d-60

  double precision, parameter :: mev_to_erg         = 1.60217733d-6
  double precision, parameter :: erg_to_mev         = 6.24150636d5
  double precision, parameter :: amu_cgs            = 1.66053873d-24
  double precision, parameter :: massn_cgs          = 1.674927211d-24
  double precision, parameter :: amu_mev            = 931.49432d0
  double precision, parameter :: kb_erg             = 1.380658d-16
  double precision, parameter :: kb_mev             = 8.61738568d-11
  double precision, parameter :: temp_mev_to_kelvin = 1.1604447522806d10
  double precision, parameter :: planck             = 6.626176d-27
  double precision, parameter :: avo                = 6.0221367d23
  double precision, parameter :: hbarc_mevcm        = 1.97326966d-11


  !public

  procedure(sub_get_pressure_one_grid), pointer        :: eos_get_pressure_one_grid        => null()
  procedure(sub_get_eps_one_grid), pointer             :: eos_get_eps_one_grid             => null()
  procedure(sub_get_cs2_one_grid), pointer             :: eos_get_cs2_one_grid             => null()

  procedure(sub_get_eps_range), pointer                :: eos_get_eps_range                => null()
  procedure(sub_get_temp_one_grid), pointer            :: eos_get_temp_one_grid            => null()
  procedure(sub_get_beta_eqm_ye_one_grid), pointer     :: eos_get_beta_eqm_ye_one_grid     => null()
  procedure(sub_eps_get_all_one_grid), pointer         :: eos_eps_get_all_one_grid         => null()
  procedure(sub_temp_get_all_one_grid), pointer        :: eos_temp_get_all_one_grid        => null()

  abstract interface

     subroutine sub_get_eps_range(rho,epsmin,epsmax,ye)
       double precision, intent(in) :: rho
       double precision, intent(out):: epsmin, epsmax
       double precision, intent(in), optional :: ye
     end subroutine sub_get_eps_range

     subroutine sub_get_pressure_one_grid(prs,rho,eps,temp,ye)
       double precision, intent(inout) :: prs
       double precision, intent(in) :: rho
       double precision, intent(in) :: eps
       double precision, intent(inout), optional :: temp ! temp may be updated when find prs
       double precision, intent(in), optional :: ye
     end subroutine sub_get_pressure_one_grid
   
     subroutine sub_get_eps_one_grid(prs,rho,eps,temp,ye)
       double precision, intent(in) :: prs
       double precision, intent(in) :: rho
       double precision, intent(inout) :: eps
       double precision, intent(in), optional :: temp, ye
     end subroutine sub_get_eps_one_grid
   
     subroutine sub_get_cs2_one_grid(cs2,rho,eps,temp,ye)
       double precision, intent(inout) :: cs2
       double precision, intent(in) :: rho
       double precision, intent(in) :: eps
       double precision, intent(in), optional :: temp, ye
     end subroutine sub_get_cs2_one_grid

     subroutine sub_get_temp_one_grid(rho,eps,temp,ye)
       double precision, intent(in) :: rho
       double precision, intent(in) :: eps
       double precision, intent(in) :: ye
       double precision, intent(inout) :: temp
     end subroutine sub_get_temp_one_grid

     subroutine sub_get_beta_eqm_ye_one_grid(rho,temp,ye)
       double precision, intent(in)    :: rho
       double precision, intent(in)    :: temp
       double precision, intent(inout) :: ye
     end subroutine sub_get_beta_eqm_ye_one_grid

     subroutine sub_eps_get_all_one_grid(rho,eps,ye,temp,prs,ent,cs2,dedt,&
                 dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu)
       double precision, intent(inout) :: rho, eps, ye
       double precision, intent(inout), optional :: ent, prs, temp, cs2, dedt
       double precision, intent(inout), optional :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar
       double precision, intent(inout), optional :: mu_e,mu_n,mu_p,muhat,munu
     end subroutine sub_eps_get_all_one_grid

     subroutine sub_temp_get_all_one_grid(rho,temp,ye,eps,prs,ent,cs2,dedt,&
                 dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu)
       double precision, intent(inout) :: rho
       double precision, intent(in) :: ye
       double precision, intent(inout) :: eps, temp

       double precision, intent(inout), optional :: prs,ent,cs2,dedt,&
                   dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu

     end subroutine sub_temp_get_all_one_grid

  end interface

contains

  subroutine eos_check

    if (eos_type == -1) call mpistop("Error: no eos module is loaded")

    ! Checks whether the required physics methods have been defined
    if (.not. associated(eos_get_pressure_one_grid)) &
         call mpistop("Error: eos_get_pressure_one_grid not defined")

    if (.not. associated(eos_get_eps_one_grid)) &
         call mpistop("Error: eos_get_eps_one_grid not defined")

    if (.not. associated(eos_get_cs2_one_grid)) &
         call mpistop("Error: eos_get_cs2_one_grid not defined")

    if (.not. associated(eos_get_eps_range)) then
       if (eos_type == tabulated) &
         call mpistop("Error: eos_get_eps_range is not defined but using tabulated eos")
       eos_get_eps_range => unlimited_eps_range
    end if
  end subroutine eos_check

  subroutine eos_atmo_activate()
    use mod_global_parameters
    use mod_gmunu
    integer                      :: n
    character(len=name_len)      :: atmo_type = ""

    namelist /atmo_list/ small_rho, small_rho_fac, atmo_type, &
                         atmo_gamma, atmo_adiab, &
                         atmo_temp, atmo_ye, atmo_beta_eqm

    do n = 1, size(par_files)
       open(unitpar, file=trim(par_files(n)), status="old")
       read(unitpar, atmo_list, end=111)
111    close(unitpar)
    end do

    if (small_rho_fac == 0.0d0 .or. small_rho_fac >= 1.0d0 ) call mpistop ("Error: small_rho_fac >= 1 or == 0")
    if (small_rho <= 0.0d0) call mpistop  ("Error: small_rho <= 0")

    select case (atmo_type)
    case ('vacuum') 
      atmo_eos_type = vacuum
    case ('polytrope','idealgas') 
      atmo_eos_type = polytrope
      ! check if the parameters are make sences
      if (atmo_gamma <= 0.0d0) call mpistop ("Error: atmo_gamma <= 0")
      if (atmo_adiab < 0.0d0) call mpistop  ("Error: atmo_adiab < 0")
    case ('tabulated') 
      atmo_eos_type = tabulated
      if (atmo_ye < 0.0d0 .or. atmo_ye > 1.0d0 ) call mpistop ("Error: atmo_ye doesnt make sense")
      if (atmo_temp < 0.0d0) call mpistop  ("Error: atmo_temp < 0")
    case default
      error stop "this atmo_eos_type is not supported"
    end select
  end subroutine eos_atmo_activate

  subroutine eos_initialize_atmo(rho_ref)
    implicit none
    double precision, intent(in) :: rho_ref

    if (eos_type == -1) &
      call mpistop("Error: atmosphere can only be initialized after activitating eos")

    !Note: the input is the rho_ref
    ! limit atmosphere value
    small_rho_thr = max(rho_ref * small_rho, eos_rhomin)

    select case (atmo_eos_type)
    case (vacuum)
      small_rho_thr = epsilon(0.0d0)
      small_rho     = smalldouble
      small_press   = 0.0d0
      small_eps     = eos_epsmin
      atmo_cs2      = 0.0d0

    case (polytrope)
      ! set small_rho_thr
      small_rho     = small_rho_thr
      small_rho_thr = small_rho / small_rho_fac

      ! based on polytrope
      small_press = atmo_adiab * small_rho**atmo_gamma
      small_eps   = atmo_adiab * small_rho**( atmo_gamma - 1.0d0 ) &
                    / ( atmo_gamma - 1.0d0 )
      atmo_cs2 = atmo_gamma * ( atmo_gamma - 1.0d0 ) * small_press &
           / ( small_rho * ( atmo_gamma - 1.0d0 ) + atmo_gamma * small_press )

      ! make sure all the allowed min value are equal to atmosphere
      eos_rhomin = small_rho
      eos_epsmin = small_eps

    case (tabulated)
      ! limit atmo_temp
      atmo_temp     = max(atmo_temp, eos_tempmin)

      ! get beta equilibrium ye
      if (atmo_beta_eqm) call eos_get_beta_eqm_ye_one_grid(small_rho_thr, atmo_temp, atmo_ye)
      call eos_temp_get_all_one_grid(small_rho_thr, atmo_temp, atmo_ye,&
         small_eps, prs=small_press)

      ! set small_rho_thr
      small_rho     = small_rho_thr
      small_rho_thr = small_rho / small_rho_fac

      ! fixme: this may cause problem
      atmo_cs2      = 0.0d0
    case default
      call mpistop("this atmo_eos_type is not supported")
    end select

    small_D   = small_rho_thr
    small_tau = small_rho_thr * small_eps ! this will no longer use
  
  end subroutine eos_initialize_atmo

  subroutine atmo_get_pressure_one_grid(prs,rho,eps)
    implicit none
    double precision, intent(inout) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    prs = small_press
    !select case (atmo_eos_type)
    !case (polytrope) 
    !  prs = atmo_adiab * rho**atmo_gamma
    !case (idealgas)
    !  prs = ( atmo_gamma - 1.0d0 ) * rho * eps
    !case default
    !  error stop "this atmo_eos is not supported"
    !end select
  end subroutine atmo_get_pressure_one_grid

  subroutine atmo_get_eps_one_grid(prs,rho,eps)
    implicit none
    double precision, intent(in) :: prs
    double precision, intent(in) :: rho
    double precision, intent(inout) :: eps
    ! this is always polytrope
    eps = small_eps
    !eps = prs / rho / ( atmo_gamma - 1.0d0 )
  end subroutine atmo_get_eps_one_grid

  subroutine atmo_get_cs2_one_grid(cs2,rho,eps)
    implicit none
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision             :: prs
    double precision             :: dpde, h
    double precision             :: dpdrho

    cs2 = atmo_cs2
    !select case (atmo_eos_type)
    !case (polytrope) 
    !  prs = atmo_adiab * rho**atmo_gamma
    !  h = 1.0d0 + eps + prs/rho
    !  ! use prs and h
    !  cs2= atmo_gamma * prs / ( rho * h )
    !  ! use prs only
    !  !cs2= atmo_gamma*( atmo_gamma - 1.0d0 ) * prs &
    !  !     /( rho*( atmo_gamma - 1.0d0 ) + atmo_gamma*prs )
    !case (idealgas)
    !  prs = ( atmo_gamma - 1.0d0 ) * rho * eps
    !  dpde = (atmo_gamma - 1.0d0 ) * rho
    !  dpdrho = (atmo_gamma - 1.0d0 ) * eps
  
    !  cs2= dpdrho+dpde*prs/rho**2
  
    !  cs2 = cs2 / ( (1.0d0 + prs/rho) + eps )
    !case default
    !  error stop "this atmo_eos is not supported"
    !end select
  end subroutine atmo_get_cs2_one_grid

  subroutine unlimited_eps_range(rho, epsmin, epsmax, ye)
     double precision, intent(in)    :: rho
     double precision, intent(out)   :: epsmin, epsmax
     double precision, intent(in), optional :: ye
     ! if it is not using tabulated eos
     epsmin = small_eps
     epsmax = bigdouble
  end subroutine unlimited_eps_range

end module mod_eos
