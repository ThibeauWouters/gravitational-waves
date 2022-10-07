!> Module for eos
module mod_eos_idealgas

  use mod_gmunu

  implicit none
  public

  double precision, public                :: eos_gamma = 5.d0/3.0d0
  double precision, public                :: eos_adiab = 1.0d0

contains

  !> Read this module's parameters from a file
  subroutine eos_idealgas_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /eos_idealgas_list/ eos_gamma, eos_adiab

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, eos_idealgas_list, end=111)
111    close(unitpar)
    end do

  end subroutine eos_idealgas_read_params

  subroutine eos_idealgas_activate()
    use mod_eos

    call eos_idealgas_read_params(par_files)
    eos_type = idealgas
    ! check if the parameters are make sences
    if (eos_gamma <= 0.0d0) call mpistop ("Error: eos_gamma <= 0")
    if (eos_adiab < 0.0d0) call mpistop  ("Error: eos_adiab < 0")

    eos_get_pressure_one_grid         => idealgas_get_pressure_one_grid
    eos_get_eps_one_grid              => idealgas_get_eps_one_grid
    eos_get_cs2_one_grid              => idealgas_get_cs2_one_grid


  end subroutine eos_idealgas_activate

  subroutine idealgas_get_pressure_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(inout) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(inout), optional :: temp
    double precision, intent(in), optional    :: ye

    if (rho<small_rho_thr) then
       call atmo_get_pressure_one_grid(prs,rho,eps)
       return
    end if
    prs = ( eos_gamma - 1.0d0 ) * rho * eps

  end subroutine idealgas_get_pressure_one_grid

  subroutine idealgas_get_eps_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(in) :: prs
    double precision, intent(in) :: rho
    double precision, intent(inout) :: eps
    double precision, intent(in), optional :: temp, ye

    if (rho<small_rho_thr) then
       call atmo_get_eps_one_grid(prs,rho,eps)
       return
    end if
    eps = prs / rho / ( eos_gamma - 1.0d0 )

  end subroutine idealgas_get_eps_one_grid

  subroutine idealgas_get_cs2_one_grid(cs2,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: temp, ye

    double precision             :: prs
    double precision             :: dpde
    double precision             :: dpdrho

    if (rho<small_rho_thr) then
       call atmo_get_cs2_one_grid(cs2,rho,eps)
       return
    end if

    prs = ( eos_gamma - 1.0d0 ) * rho * eps
    dpde = (eos_gamma - 1.0d0 ) * rho
    dpdrho = (eos_gamma - 1.0d0 ) * eps

    cs2= dpdrho+dpde*prs/rho**2

    cs2 = cs2 / ( (1.0d0 + prs/rho) + eps )

  end subroutine idealgas_get_cs2_one_grid


end module mod_eos_idealgas
