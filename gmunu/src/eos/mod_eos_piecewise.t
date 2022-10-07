!> Module for eos
module mod_eos_piecewise

  use mod_gmunu

  implicit none
  public

  double precision, public                :: eos_gamma_1 = 1.3569d0 !1.3569
  double precision, public                :: eos_gamma_2 = 2.664d0  !2.664
  double precision, public                :: eos_gamma_3 = 2.194d0  !2.194
  double precision, public                :: eos_gamma_4 = 2.304d0  !2.304
  double precision, public                :: eos_adiab_1 = 8.951d-2 !0.08951
  double precision, public                :: eos_adiab_2 = 1.371d4  !13710
  double precision, public                :: eos_adiab_3 = 4.836d2  !483.6
  double precision, public                :: eos_adiab_4 = 9.805d2  !980.5 
  double precision, public                :: eos_rho_1 = 1.07896d-4 
  double precision, public                :: eos_rho_2 = 8.11905d-4 
  double precision, public                :: eos_rho_3 = 1.61996d-3 
  double precision, public                :: eos_a_1 = 0.0d0      !0.0000
  double precision, public                :: eos_a_2 = 7.5597d-3  
  double precision, public                :: eos_a_3 = -1.5795d-2 !-0.015795
  double precision, public                :: eos_a_4 = 1.17801d-4  

contains

  !> Read this module's parameters from a file
  subroutine eos_piecewise_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /eos_piecewise_list/ eos_gamma_1, eos_gamma_2, eos_gamma_3, &
                                  eos_gamma_4, eos_adiab_1, eos_adiab_2, &
                                  eos_adiab_3, eos_adiab_4, eos_rho_1,   &
                                  eos_rho_2, eos_rho_3, eos_a_1, &
                                  eos_a_2, eos_a_3, eos_a_4

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, eos_piecewise_list, end=111)
111    close(unitpar)
    end do

  end subroutine eos_piecewise_read_params

  subroutine eos_piecewise_activate()
    use mod_eos

    call eos_piecewise_read_params(par_files)
    eos_type = piecewise
    ! check if the parameters are make sences
    if (eos_gamma_1 <= 0.0d0) call mpistop ("Error: eos_gamma_1 <= 0")
    if (eos_adiab_1 < 0.0d0) call mpistop  ("Error: eos_adiab_1 < 0")
    if (eos_gamma_2 <= 0.0d0) call mpistop ("Error: eos_gamma_2 <= 0")
    if (eos_adiab_2 < 0.0d0) call mpistop  ("Error: eos_adiab_2 < 0")
    if (eos_gamma_3 <= 0.0d0) call mpistop ("Error: eos_gamma_3 <= 0")
    if (eos_adiab_3 < 0.0d0) call mpistop  ("Error: eos_adiab_3 < 0")
    if (eos_gamma_4 <= 0.0d0) call mpistop ("Error: eos_gamma_4 <= 0")
    if (eos_adiab_4 < 0.0d0) call mpistop  ("Error: eos_adiab_4 < 0")
    eos_get_pressure_one_grid         => piecewise_get_pressure_one_grid
    eos_get_eps_one_grid              => piecewise_get_eps_one_grid
    eos_get_cs2_one_grid              => piecewise_get_cs2_one_grid


  end subroutine eos_piecewise_activate

  subroutine piecewise_get_pressure_one_grid(prs,rho,eps,temp,ye)

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
    
    if (rho<=eos_rho_1) then
        prs = eos_adiab_1 * rho**eos_gamma_1
    
    else if (rho<=eos_rho_2) then
        prs = eos_adiab_2 * rho**eos_gamma_2
    
    else if (rho<=eos_rho_3) then
        prs = eos_adiab_3 * rho**eos_gamma_3

    else
        prs = eos_adiab_4 * rho**eos_gamma_4
    
    end if
  end subroutine piecewise_get_pressure_one_grid

  subroutine piecewise_get_eps_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(in) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in), optional :: temp, ye
    double precision, intent(inout) :: eps

    if (rho<small_rho_thr) then
       call atmo_get_eps_one_grid(prs,rho,eps)
       return
    end if
    !eps = prs / rho / ( eos_gamma - 1.0d0 )
    
    if (rho<=eos_rho_1) then 
        eps = eos_a_1 + (prs / rho / (eos_gamma_1 - 1.0d0)) 

    else if (rho<=eos_rho_2) then 
        eps = eos_a_2 + (prs / rho / (eos_gamma_2 - 1.0d0)) 
    
    else if (rho<=eos_rho_3) then 
        eps = eos_a_3 + (prs / rho / (eos_gamma_3 - 1.0d0)) 

    else  
        eps = eos_a_4 + (prs / rho / (eos_gamma_4 - 1.0d0)) 
    
    end if
    !eps = eos_adiab * rho**( eos_gamma - 1.0d0 ) &
    !      / ( eos_gamma - 1.0d0 )

  end subroutine piecewise_get_eps_one_grid

  subroutine piecewise_get_cs2_one_grid(cs2,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: temp, ye

    double precision             :: prs, h

    if (rho<small_rho_thr) then
       call atmo_get_cs2_one_grid(cs2,rho,eps)
       return
    end if

    if (rho<=eos_rho_1) then
        prs = eos_adiab_1 * rho**eos_gamma_1
    
    else if (rho<=eos_rho_2) then
        prs = eos_adiab_2 * rho**eos_gamma_2
    
    else if (rho<=eos_rho_3) then
        prs = eos_adiab_3 * rho**eos_gamma_3

    else
        prs = eos_adiab_4 * rho**eos_gamma_4
    
    end if
   
    h = 1.0d0 + eps + prs/rho

    ! use prs and h
    !cs2= eos_gamma * prs / ( rho * h )

    if (rho<=eos_rho_1) then        
        cs2= eos_gamma_1 * prs / ( rho * h )
    
    else if (rho<=eos_rho_2) then        
        cs2= eos_gamma_2 * prs / ( rho * h )
    
    else if (rho<=eos_rho_3) then        
        cs2= eos_gamma_3 * prs / ( rho * h )

    else        
        cs2= eos_gamma_4 * prs / ( rho * h )
    
    end if

    ! use prs only
    !cs2= eos_gamma*( eos_gamma - 1.0d0 ) * prs &
    !     /( rho*( eos_gamma - 1.0d0 ) + eos_gamma*prs )

  end subroutine piecewise_get_cs2_one_grid

end module mod_eos_piecewise
