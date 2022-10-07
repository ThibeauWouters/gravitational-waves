#include "cpp_macros.h"
module m_finite_difference

  implicit none
  public

  integer, parameter  :: dxdy = 1
  integer, parameter  :: dydx = -1

contains

  subroutine Laplacian(lop, u_in, idr2, dLdu, face_coeff_in, pre_fac_in)
    real(dp), intent(out)               :: lop
    real(dp), intent(out), optional     :: dLdu
    real(dp), intent(in)                :: u_in(-1:1, NDIM)
    real(dp), intent(in)                :: idr2(NDIM)
    real(dp), intent(in), optional      :: face_coeff_in(0:1, NDIM)
    real(dp), intent(in), optional      :: pre_fac_in(NDIM)

    integer                             :: idir
    real(dp)                            :: u(-1:1, NDIM) = 0.0d0
    real(dp)                            :: face_coeff(0:1, NDIM) = 1.0d0
    real(dp)                            :: pre_fac(NDIM) = 1.0d0

    if (present(face_coeff_in)) then
       face_coeff = face_coeff_in
    end if

    if (present(pre_fac_in)) then
       pre_fac = pre_fac_in
    end if

    if (present(dLdu)) then
       dLdu = 0.0d0
       do idir = 1, NDIM
          dLdu = dLdu + pre_fac(idir) * idr2(idir) * ( &
                       - face_coeff(1, idir) - face_coeff(0, idir) )
       end do
    end if

    lop = 0.0d0
    u = u_in
    do idir = 1, NDIM
       lop = lop + pre_fac(idir) * idr2(idir) * ( &
                    face_coeff(1, idir) * (u(1, idir) - u(0, idir))   - &
                    face_coeff(0, idir) * (u(0, idir) - u(-1, idir)) )
    end do

  end subroutine Laplacian

  real(dp) function dudr(u, dr)
    real(dp), intent(in)                :: u(-1:1)
    real(dp), intent(in)                :: dr
    dudr = ( u(1) - u(-1) ) / (2.0d0 * dr)
  end function dudr

  subroutine derivatives_2nd(d2udr2, u, idr2, dLdu, face_coeff_in, pre_fac_in)
    real(dp), intent(out)               :: d2udr2 
    real(dp), intent(out), optional     :: dLdu
    real(dp), intent(in)                :: u(-1:1)
    real(dp), intent(in)                :: idr2
    real(dp), intent(in), optional      :: face_coeff_in(0:1)
    real(dp), intent(in), optional      :: pre_fac_in

    integer                             :: idir
    real(dp)                            :: face_coeff(0:1) = 1.0d0
    real(dp)                            :: pre_fac = 1.0d0

    if (present(face_coeff_in)) then
       face_coeff = face_coeff_in
    end if

    if (present(pre_fac_in)) then
       pre_fac = pre_fac_in
    end if

    if (present(dLdu)) then
       dLdu = pre_fac * idr2 * ( &
                       - face_coeff(1) - face_coeff(0) )
    end if

    d2udr2 = pre_fac * idr2 * ( &
                    face_coeff(1) * (u(1) - u(0))   - &
                    face_coeff(0) * (u(0) - u(-1)) )

  end subroutine derivatives_2nd

  subroutine mixed_derivatives_2nd(d2udxdy, u, dx, dy, i, j, nc, dLdu)
    real(dp), intent(out)               :: d2udxdy
    real(dp), intent(out), optional     :: dLdu
    real(dp), intent(in)                :: u(-1:1,-1:1)
    real(dp), intent(in)                :: dx, dy
    integer, intent(in)                 :: i,j,nc

    real(dp)                            :: idxdy

    idxdy = 1.0_dp/ ( dx * dy )
    if ( (i==1 .and. j==1).or.(i==nc .and. j==nc) ) then
       ! at lower left or upper right corner
       d2udxdy = - u(1,-1) - u(-1,1) + u(1,0) + u(0,1) + u(-1,0) + u(0,-1) - 2.0_dp * u(0,0)
       d2udxdy = d2udxdy * 0.5_dp * idxdy
       if (present(dLdu)) then
          dLdu =  - idxdy
       end if
    else if ( (i==1 .and. j==nc).or.(i==nc .and. j==1) ) then
       ! at lower right or upper left corner
       d2udxdy = - u(1,1) - u(-1,-1) + u(1,0) + u(0,1) + u(-1,0) + u(0,-1) - 2.0_dp * u(0,0)
       d2udxdy = d2udxdy * 0.5_dp * idxdy
       if (present(dLdu)) then
          dLdu =  - idxdy
       end if
    else
       ! normal cases
       d2udxdy = u(1,1) - u(1,-1) - u(-1,1) + u(-1,-1)
       d2udxdy = d2udxdy * 0.25_dp * idxdy
       if (present(dLdu)) then
          dLdu = 0.0d0
       end if
    end if
  end subroutine mixed_derivatives_2nd

end module m_finite_difference
