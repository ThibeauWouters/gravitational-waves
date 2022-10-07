module mod_grhd_phys_one_d_core
  use mod_physics
  use mod_grhd_phys_parameters, only: r_core, oneDcore

  implicit none
  private

  ! Public methods
  public :: grhd_phys_one_d_core_init

contains

  !> Initialize the module
  subroutine grhd_phys_one_d_core_init()
    use mod_usr_methods

    integer :: itr, idir

    {^NOONED
    if ( oneDcore ) then
      usr_set_flux => grhd_set_flux
      usr_internal_bc => grhd_internal_bc
    end if
    }

  end subroutine grhd_phys_one_d_core_init

  {^NOONED
  !> 1D core treatment: internal_bc
  subroutine grhd_internal_bc(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nprim)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    ! fixme: should be averaged radial direction

    ! set vtheta and vphi = zero
    where ( x(ixO^S,r_) <= r_core )
      w(ixO^S,W_vel(2)) = 0.0d0
      ! Note: v3 is large for rapity rotating NS
      !w(ixO^S,W_vel(3)) = 0.0d0
    end where
  end subroutine grhd_internal_bc

  !> 1D core treatment: flux
  subroutine grhd_set_flux(ixI^L,ixC^L,idim,fC,xi)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixC^L, idim
    ! face-center xi
    double precision,intent(in)  :: xi(ixI^S,1:ndim)
    ! face-center flux
    double precision,intent(inout) :: fC(ixI^S,1:ncons,1:ndim)

    integer                      :: i_cons
    
    ! nothing to do for radial direction flux
    ! fixme: should be averaged
    if (idim == 1) return

    ! non-radial direction flux are set to be zero
    do i_cons = 1, ncons
       where (xi(ixC^S,r_) <= r_core )
         fC(ixC^S,i_cons,idim) = 0.0d0
       end where
    end do
  end subroutine grhd_set_flux
  }

end module mod_grhd_phys_one_d_core
