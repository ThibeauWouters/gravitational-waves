!> Module for eos
! Note: this module is used only for code testing, so it might not efficient
module mod_eos_hybrid

  use mod_gmunu

  implicit none
  public

  double precision, public                :: eos_gamma_1 = 4.d0/3.0d0
  double precision, public                :: eos_gamma_2 = 2.5d0
  double precision, public                :: eos_gamma_th = 1.5d0

  double precision, public                :: coeff_E(1:3)
  double precision, public                :: coeff_K(1:2)
  
  ! parameters in cgs unit
  double precision, public                :: eos_adiab = 1.2435d15 * (0.5d0)**(4.0d0/3.0d0)
  double precision, public                :: rho_nuc = 2.0d14
  !double precision, public                :: eos_adiab = 0.46757897525669273D0
  !double precision, public                :: rho_nuc = 3.2386069399999997D-004

contains

  !> Read this module's parameters from a file
  subroutine eos_hybrid_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /eos_hybrid_list/ eos_gamma_1, eos_gamma_2, eos_gamma_th, &
                   eos_adiab, rho_nuc

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, eos_hybrid_list, end=111)
111    close(unitpar)
    end do

  end subroutine eos_hybrid_read_params

  subroutine eos_hybrid_activate()
    use mod_eos

    call eos_hybrid_read_params(par_files)
    eos_type = hybrid

    ! cgs to code unit
    !eos_adiab = 1.2435d15 * (0.5d0)**(4.0d0/3.0d0)
    !eos_adiab = eos_adiab * 1.80171810d-39 / ((1.61930347d-18)**(4.0d0/3.0d0))
    !rho_nuc = 2.0d14 * rho_gf
    !write(*,*) eos_adiab, rho_nuc
    !write(*,*) eos_gamma_1, eos_gamma_2, eos_gamma_th
    !stop

    if (eos_gamma_1 <= 0.0d0) call mpistop ("Error: eos_gamma_1 <= 0")
    if (eos_gamma_2 <= 0.0d0) call mpistop ("Error: eos_gamma_2 <= 0")
    if (eos_gamma_th <= 0.0d0) call mpistop ("Error: eos_gamma_th <= 0")
    if (eos_adiab < 0.0d0) call mpistop  ("Error: eos_adiab < 0")

    coeff_K(1) = eos_adiab
    coeff_E(1) = coeff_K(1) / ( eos_gamma_1 - 1.0d0 )
    coeff_K(2) = ( eos_gamma_1 - 1.0d0 ) * coeff_E(1) * rho_nuc**( eos_gamma_1 - eos_gamma_2 )
    coeff_E(2) = coeff_K(2) / ( eos_gamma_2 - 1.0d0 )
    coeff_E(3) = ( eos_gamma_2 - eos_gamma_1 )  / ( eos_gamma_2 - 1.0d0 ) * coeff_E(1) * rho_nuc**( eos_gamma_1 - 1.0d0 )

    eos_get_pressure_one_grid         => hybrid_get_pressure_one_grid
    eos_get_eps_one_grid              => hybrid_get_eps_one_grid
    eos_get_cs2_one_grid              => hybrid_get_cs2_one_grid

  end subroutine eos_hybrid_activate

  subroutine hybrid_get_pressure_one_grid(prs_out,rho_in,eps_in,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(inout) :: prs_out
    double precision, intent(in) :: rho_in
    double precision, intent(in) :: eps_in
    double precision, intent(inout), optional :: temp
    double precision, intent(in), optional    :: ye

    double precision             :: prs
    double precision             :: rho
    double precision             :: eps

    double precision             :: Kx, Ex, Gx, Ex3
    double precision             :: up
    double precision             :: prs_cold, prs_th
    double precision             :: eps_th, eps_cold
    double precision             :: dpde
    double precision             :: dpdrho

    if (rho_in<small_rho_thr) then
       call atmo_get_pressure_one_grid(prs_out,rho_in,eps_in)
       return
    end if

    ! code unit to cgs unit
    rho = rho_in / rho_gf
    eps = eps_in / eps_gf

!    prs = ( eos_gamma(rho) - eos_gamma_th )/( eos_gamma(rho) - 1.0d0 ) &
!              * eos_adiab * rho_nuc**(eos_gamma_1-eos_gamma(rho)) * rho**eos_gamma(rho) &
!          - ( eos_gamma(rho) - eos_gamma_1 ) * ( eos_gamma_th - 1.0d0 ) &
!                  /( eos_gamma_1 - 1.0d0 ) / ( eos_gamma_2 - 1.0d0 ) &
!              * eos_adiab * rho_nuc**(eos_gamma_1-1.0d0) * rho & 
!          + ( eos_gamma_th - 1.0d0 ) * rho * eps
    if ( rho < rho_nuc ) then
       Kx = coeff_K(1)
       Ex = coeff_E(1)
       Gx = eos_gamma_1
       Ex3 = 0.0d0
    else
       Kx = coeff_K(2)
       Ex = coeff_E(2)
       Gx = eos_gamma_2
       Ex3 = coeff_E(3)
    end if

    ! thermal part
    up = Ex * rho**Gx + Ex3 * rho
    eps_th = eps * rho - up
    prs_th = ( eos_gamma_th - 1.0d0 ) * eps_th
    if ( prs_th < 0.0d0 ) then
       prs_th = 0.0d0
    end if

    ! cold part
    prs_cold = Kx * rho**Gx
    
    prs = prs_cold + prs_th
    
    prs_out = prs * press_gf

  end subroutine hybrid_get_pressure_one_grid

  subroutine hybrid_get_eps_one_grid(prs_in,rho_in,eps_out,temp,ye)
    use mod_eos
    implicit none
    double precision, intent(in) :: prs_in
    double precision, intent(in) :: rho_in
    double precision, intent(inout) :: eps_out
    double precision, intent(in), optional :: temp, ye

    double precision             :: prs
    double precision             :: rho
    double precision             :: eps

    if (rho_in<small_rho_thr) then
       call atmo_get_eps_one_grid(prs_in,rho_in,eps_out)
       return
    end if

    ! code unit to cgs unit
    rho = rho_in / rho_gf
    prs = prs_in / press_gf

!    ! Note: this contains only the eps_cold part!
!    if (rho <= rho_nuc) then
!       eps = eos_adiab * rho**( eos_gamma_1 - 1.0d0 ) &
!             / ( eos_gamma_1 - 1.0d0 )
!    else
!       eps = eos_adiab * rho_nuc**( eos_gamma_1 - eos_gamma_2 ) &
!             / ( eos_gamma_2 - 1.0d0 ) &
!             * ( rho**(eos_gamma_2-1.0d0) &
!                - ( eos_gamma_1 - eos_gamma_2 ) / ( eos_gamma_1 - 1.0d0 )*rho_nuc**(eos_gamma_2-1.0d0))
!    end if
    eps = eos_adiab * rho**( eos_gamma_1 - 1.0d0 ) &
             / ( eos_gamma_1 - 1.0d0 )

    eps_out = eps * eps_gf
  end subroutine hybrid_get_eps_one_grid

  subroutine hybrid_get_cs2_one_grid(cs2,rho_in,eps_in,temp,ye)
    use mod_eos
    implicit none
    
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho_in
    double precision, intent(in) :: eps_in
    double precision, intent(in), optional :: temp, ye

    double precision             :: prs, eps, rho
    double precision             :: Kx, Ex, Gx, Ex3
    double precision             :: up
    double precision             :: prs_cold, prs_th
    double precision             :: eps_th, eps_cold
    double precision             :: dpde
    double precision             :: dpdrho

    if (rho_in<small_rho_thr) then
       call atmo_get_cs2_one_grid(cs2,rho_in,eps_in)
       return
    end if

    ! code unit to cgs unit
    rho = rho_in / rho_gf
    eps = eps_in / eps_gf

!    if (rho <= rho_nuc) then
!       prs_cold = eos_adiab * rho**eos_gamma_1
!    else
!       prs_cold = eos_adiab * rho**eos_gamma_2 &
!                  * rho_nuc*(eos_gamma_1 - eos_gamma_2)
!    end if
!
!    if (rho <= rho_nuc) then
!       eps_cold = eos_adiab * rho**( eos_gamma_1 - 1.0d0 ) &
!             / ( eos_gamma_1 - 1.0d0 )
!    else
!       eps_cold = eos_adiab * rho_nuc**( eos_gamma_1 - eos_gamma_2 ) &
!             / ( eos_gamma_2 - 1.0d0 ) &
!             * ( rho**(eos_gamma_2-1.0d0) &
!                - ( eos_gamma_1 - eos_gamma_2 ) / ( eos_gamma_1 - 1.0d0 )*rho_nuc**(eos_gamma_2-1.0d0))
!    end if
!
!    eps_th = eps - eps_cold
!    eps_th = max(eps_th, 0.0d0)
!    prs_th = ( eos_gamma_th - 1.0d0 ) * rho * eps_th
!    if (prs_th < 0.0d0) stop "WTF"
!
!    prs_cold = prs_cold * press_gf
!    prs_th = prs_th * press_gf
!
!    prs = prs_cold + prs_th
!
!    cs2 = eos_gamma(rho) * prs_cold + eos_gamma_th * prs_th
!    cs2 = cs2 / ( 1.0d0 + prs/rho_in + eps_in ) / rho_in
!
!    cs2 = cs2 / clight**2

    if ( rho < rho_nuc ) then
       Kx = coeff_K(1)
       Ex = coeff_E(1)
       Gx = eos_gamma_1
       Ex3 = 0.0d0
    else
       Kx = coeff_K(2)
       Ex = coeff_E(2)
       Gx = eos_gamma_2
       Ex3 = coeff_E(3)
    end if

    ! thermal part
    up = Ex * rho**Gx + Ex3 * rho
    eps_th = eps * rho - up
    prs_th = ( eos_gamma_th - 1.0d0 ) * eps_th
    dpdrho = ( eos_gamma_th - 1.0d0 ) &
            * ( eps - Ex * Gx * rho**(Gx - 1.0d0) - Ex3 )
    dpde = ( eos_gamma_th - 1.0d0 ) * rho
    if ( prs_th < 0.0d0 ) then
       prs_th = 0.0d0
       dpdrho = 0.0d0
       dpde = 0.0d0
    end if

    ! cold part
    prs_cold = Kx * rho**Gx
    dpdrho = dpdrho + Gx * prs_cold / rho
    
    prs = prs_cold + prs_th
    cs2 = dpdrho + dpde * prs / rho**2

    prs = prs * press_gf
    cs2 = cs2 / clight**2 / ( 1.0d0 + prs/rho_in + eps_in ) 

  end subroutine hybrid_get_cs2_one_grid

  double precision function eos_gamma(rho)
    use mod_eos, only: rho_gf
    implicit none
    double precision,intent(in)         :: rho
    ! note that the input rho is in cgs
    if (rho <= rho_nuc) then
       eos_gamma = eos_gamma_1
    else
       eos_gamma = eos_gamma_2
    end if
  end function eos_gamma

end module mod_eos_hybrid
