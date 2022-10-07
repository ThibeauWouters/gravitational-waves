!> Module for eos
module mod_eos_hqs

  use mod_gmunu

  implicit none
  public

  double precision, public                :: eos_gamma = 5.d0/3.0d0
  double precision, public                :: eos_adiab = 1.0d2
  double precision, public                :: eos_bag_oneforth = 1.7d2 ! in MeV
  double precision, public                :: eos_rho_hm = 1.0d0
  double precision, public                :: eos_rho_qm = 1.0d0
  double precision, public                :: eos_delta  = 2.0d0

  double precision, parameter             :: press_natural_to_code_unit = 3.75697085d-1  ! for 1 GeV^4
  double precision                        :: eos_bag_constant

contains

  !> Read this module's parameters from a file
  subroutine eos_hqs_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /eos_hqs_list/ eos_gamma, eos_adiab, eos_bag_oneforth, &
                            eos_rho_hm, eos_rho_qm, eos_delta

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, eos_hqs_list, end=111)
111    close(unitpar)
    end do

  end subroutine eos_hqs_read_params

  subroutine eos_hqs_activate()
    use mod_eos

    call eos_hqs_read_params(par_files)
    eos_type = hqs
    ! check if the parameters are make sences
    if (eos_gamma <= 0.0d0) call mpistop ("Error: eos_gamma <= 0")
    if (eos_adiab < 0.0d0) call mpistop  ("Error: eos_adiab < 0")
    ! fixme: not detaily check the parameter

    eos_get_pressure_one_grid         => hqs_get_pressure_one_grid
    eos_get_eps_one_grid              => hqs_get_eps_one_grid
    eos_get_cs2_one_grid              => hqs_get_cs2_one_grid

    eos_bag_constant = (eos_bag_oneforth * 1.0d-3) ** 4.0d0 ! Bag constant in GeV^4
    eos_bag_constant = eos_bag_constant * press_natural_to_code_unit ! transform it into code unit

  end subroutine eos_hqs_activate

  subroutine hqs_get_pressure_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(inout) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(inout), optional :: temp
    double precision, intent(in), optional :: ye
    double precision :: prs_hm, prs_qm, alp_hqs
     
    if (rho<small_rho_thr) then
       call atmo_get_pressure_one_grid(prs,rho,eps)
       return
    end if

    prs_hm = ( eos_gamma - 1.0d0 ) * rho * eps
    prs_qm = 1.0d0/3.0d0 * (rho + rho * eps - 4.0d0 * eos_bag_constant) 
    
    if (rho<eos_rho_hm) then            
        prs = prs_hm

    else if (rho<=eos_rho_qm) then 
        alp_hqs = alp_of_rho(rho)
        prs = alp_hqs * prs_qm + (1.0d0 - alp_hqs) * prs_hm
 
    else
        prs = prs_qm
         
    end if

  end subroutine hqs_get_pressure_one_grid

  subroutine hqs_get_eps_one_grid(prs,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(in) :: prs
    double precision, intent(in) :: rho
    double precision, intent(inout) :: eps
    double precision, intent(in), optional :: temp, ye
    double precision :: eps_nom, eps_denom
    double precision :: alp_hqs
        
    if (rho<small_rho_thr) then
       call atmo_get_eps_one_grid(prs,rho,eps)
       return
    end if

    ! The initial star is assumed to be polytropic, 
    ! so only polytrope eos is used here
    eps = prs / rho / ( eos_gamma - 1.0d0 )
    return

!    if (rho<eos_rho_hm) then            
!        eps = prs / rho / ( eos_gamma - 1.0d0 )
!
!    else if (rho<=eos_rho_qm) then
!        alp_hqs = alp_of_rho(rho)
!        eps_nom = prs - alp_hqs/3.0d0 * (rho - 4.0d0 * eos_bag_constant)
!        eps_denom = (alp_hqs/3.0d0 + (1.0d0 - alp_hqs) * & 
!                    (eos_gamma - 1.0d0)) * rho 
!        eps = eps_nom / eps_denom
!
!    else
!        eps = 3.0d0 * prs / rho + 4.0d0 * eos_bag_constant / rho - 1.0d0
!
!    end if

  end subroutine hqs_get_eps_one_grid

  subroutine hqs_get_cs2_one_grid(cs2,rho,eps,temp,ye)

    use mod_eos
    implicit none
    
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: temp, ye

    double precision             :: prs
    !double precision             :: h

    double precision :: prs_hm, prs_qm, alp_hqs
    double precision :: dpde, dpde_hm, dpde_qm
    double precision :: dpdrho, dpdrho_hm, dpdrho_qm
    double precision :: dalpdrho

    if (rho<small_rho_thr) then
       call atmo_get_cs2_one_grid(cs2,rho,eps)
       return
    end if

    prs_hm = ( eos_gamma - 1.0d0 ) * rho * eps
    prs_qm = 1.0d0/3.0d0 * (rho + rho * eps - 4.0d0 * eos_bag_constant) 

    dpde_hm = (eos_gamma - 1.0d0 ) * rho
    dpdrho_hm = (eos_gamma - 1.0d0 ) * eps

    dpde_qm = 1.0d0/3.0d0 * rho
    dpdrho_qm = 1.0d0/3.0d0 * (1.0d0 + eps)

    if (rho<eos_rho_hm) then            
        prs = prs_hm 
        
        dpde = dpde_hm
        dpdrho = dpdrho_hm
    
    else if (rho<=eos_rho_qm) then
        alp_hqs = alp_of_rho(rho)
        prs = alp_hqs * prs_qm + (1.0d0 - alp_hqs) * prs_hm 

        dalpdrho = (eos_rho_qm - rho) / (eos_rho_qm - eos_rho_hm)
        dalpdrho = eos_delta * (dalpdrho ** eos_delta)
        dalpdrho = dalpdrho / (eos_rho_qm - rho)
    
        dpde = alp_hqs * dpde_qm + (1.0d0 - alp_hqs) * dpde_hm
        dpdrho = dalpdrho * prs_qm + alp_hqs * dpdrho_qm - &
                 dalpdrho * prs_hm + (1.0d0 - alp_hqs) * dpdrho_hm 
        
    else
        prs = prs_qm      
      
        dpde = dpde_qm
        dpdrho = dpdrho_qm
 
    end if

    cs2= dpdrho+dpde*prs/rho**2
    cs2 = cs2 / ( (1.0d0 + prs/rho) + eps )

  end subroutine hqs_get_cs2_one_grid

  double precision function alp_of_rho(rho)
    double precision, intent(in) :: rho
    alp_of_rho = (eos_rho_qm - rho)/(eos_rho_qm - eos_rho_hm)  
    alp_of_rho = 1.0d0 - alp_of_rho ** eos_delta   
  end function alp_of_rho


end module mod_eos_hqs
