!> Module containing all hydrodynamics
module mod_grhd
  use mod_grhd_phys
  use mod_grhd_ppm

  use mod_gmunu

  implicit none
  public

contains

  subroutine grhd_activate()
    call grhd_phys_init()
    call grhd_ppm_init()
  end subroutine grhd_activate

end module mod_grhd
