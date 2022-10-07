!> Module containing all hydrodynamics
module mod_grrhd
  use mod_grrhd_phys
  use mod_grrhd_ppm

  use mod_gmunu

  implicit none
  public

contains

  subroutine grrhd_activate()
    call grrhd_phys_init()
    call grrhd_ppm_init()
  end subroutine grrhd_activate

end module mod_grrhd
