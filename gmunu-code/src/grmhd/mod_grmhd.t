!> Module containing all hydrodynamics
module mod_grmhd
  use mod_grmhd_phys
  use mod_grmhd_ppm

  use mod_gmunu

  implicit none
  public

contains

  subroutine grmhd_activate()
    call grmhd_phys_init()
    call grmhd_ppm_init()
  end subroutine grmhd_activate

end module mod_grmhd
