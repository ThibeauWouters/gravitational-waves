!> Module containing all hydrodynamics
module mod_gremhd
  use mod_gremhd_phys
  use mod_gremhd_ppm

  use mod_gmunu

  implicit none
  public

contains

  subroutine gremhd_activate()
    call gremhd_phys_init()
    call gremhd_ppm_init()
  end subroutine gremhd_activate

end module mod_gremhd
