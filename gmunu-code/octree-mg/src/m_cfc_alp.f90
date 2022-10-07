#include "cpp_macros.h"
!> Module which contains multigrid procedures for the L operator for (alp*psi - 1)
module m_cfc_alp
  use m_data_structures
  use m_finite_difference

  implicit none
  private

  public :: cfc_alp_set_methods

contains

  subroutine cfc_alp_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    mg%vector_equation = .False.

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_lalp

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_lalp
       case default
          error stop "cfc_alp_set_methods: unsupported smoother type"
       end select
#if NDIM != 1
    case (mg_cylindrical)
       mg%box_op => box_clalp

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_clalp
       case default
          error stop "cfc_alp_set_methods: unsupported smoother type"
       end select
#endif
    case (mg_spherical)
       mg%box_op => box_slalp

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_slalp
       case default
          error stop "cfc_alp_set_methods: unsupported smoother type"
       end select
    case default
       error stop "cfc_alp_set_methods: unsupported geometry"
    end select

  end subroutine cfc_alp_set_methods

  !> Perform L operator on a box cartesian geometry
  subroutine box_lalp(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK

    real(dp)                  :: u(-1:1, NDIM)
    real(dp)                  :: dr(NDIM), idr2(NDIM)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 1
      do i = 1, nc
         ! Laplacian
         u(-1:1, 1) = cc(i-1:i+1, n)
         call Laplacian(cc(i, i_out), u, idr2)

         ! nonlinear source terms
         cc(i, i_out) = cc(i, i_out) &
              + cc(i, f1) * ( 1.0d0 + cc(i, n) )
      end do
#elif NDIM == 2
      do j = 1, nc
         do i = 1, nc
            ! Laplacian
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            call Laplacian(cc(i, j, i_out), u, idr2)

            ! nonlinear source terms
            cc(i, j, i_out) = cc(i, j, i_out) &
                       + cc(i, j, f1) * ( 1.0d0 + cc(i, j, n) )
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               ! Laplacian
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               call Laplacian(cc(i, j, k, i_out), u, idr2)
   
               ! nonlinear source terms
               cc(i, j, k, i_out) = cc(i, j, k, i_out) &
                       + cc(i, j, k, f1) * ( 1.0d0 + cc(i, j, k, n) )
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_lalp

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> cartesian geometry.
  subroutine box_gs_lalp(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di
    integer                   :: lo(NDIM), hi(NDIM), pm(NDIM)

    real(dp)                  :: u(-1:1, NDIM)

    logical                   :: redblack
    real(dp)                  :: idr2(NDIM), dr(NDIM), dr2(NDIM)
    real(dp)                  :: Lop, dLdu

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2

    where ( mg%boxes(id)%r_min + 0.5_dp * nc * mg%boxes(id)%dr.ge.0.0_dp .and. mg%symm)
      pm = 1
      lo = 1
      hi = nc
    else where
      pm = -1
      lo = nc
      hi = 1
    end where

    i0 = lo(1)
    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2 * pm(1)
    else
       di = pm(1)
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 1
      if (redblack) i0 = lo(1) + iand(lo(1),1) - iand(redblack_cntr, 1)

      do i = i0, hi(1), di
         ! Laplacian
         u(-1:1, 1) = cc(i-1:i+1, n)
         call Laplacian(Lop, u, idr2, dLdu = dLdu)

         ! nonlinear source terms
         Lop = Lop + cc(i, f1) * ( 1.0d0 + cc(i, n) )
         dLdu = dLdu + cc(i, f1)
 
         cc(i, n) = cc(i, n) - ( Lop - cc(i, mg_irhs) ) / dLdu
      end do

#elif NDIM == 2
      do j = lo(2), hi(2), pm(2)
         if (redblack) &
              i0 = lo(1) + iand(lo(1),1) - iand(ieor(redblack_cntr, j), 1)

         do i = i0, hi(1), di

            ! Laplacian 
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            call Laplacian(Lop, u, idr2, dLdu = dLdu)

            ! nonlinear source terms
            Lop = Lop + cc(i, j, f1) * ( 1.0d0 + cc(i, j, n) )
            dLdu = dLdu + cc(i, j, f1) 
    
            cc(i, j, n) = cc(i, j, n) - ( Lop - cc(i, j, mg_irhs) ) / dLdu
         end do
      end do

#elif NDIM == 3
      do k = lo(3), hi(3), pm(3)
         do j = lo(2), hi(2), pm(2)
            if (redblack) &
                 i0 = lo(1) + iand(lo(1),1) - iand(ieor(redblack_cntr, j+k), 1)
   
            do i = i0, hi(1), di
   
               ! Laplacian 
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               call Laplacian(Lop, u, idr2, dLdu = dLdu)
   
               ! nonlinear source terms
               Lop = Lop + cc(i, j, k, f1) * ( 1.0d0 + cc(i, j, k, n) )
               dLdu = dLdu + cc(i, j, k, f1) 
       
               cc(i, j, k, n) = cc(i, j, k, n) - ( Lop - cc(i, j, k, mg_irhs) ) / dLdu
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_lalp

  !> Perform L operator on a box in cylindrical geometry, using (r,z,phi)
  subroutine box_clalp(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0d0
    real(dp)                  :: pre_fac(NDIM) = 1.0d0
    real(dp)                  :: u(-1:1, NDIM)

    real(dp)                  :: dr(NDIM), idr2(NDIM)
    real(dp)                  :: r_face(1:nc+1)
    real(dp)                  :: r_cc(nc)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])

    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 2
      do j = 1, nc
         do i = 1, nc
            ! Laplacian
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)
            pre_fac(1) = r_cc(i)
            pre_fac(2) = r_cc(i)**2
            call Laplacian(cc(i, j, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
            ! nonlinear source terms
            cc(i, j, i_out) = cc(i, j, i_out) &
                       + r_cc(i)**2 * &
                         cc(i, j, f1) * ( 1.0d0 + cc(i, j, n) )
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               ! Laplacian
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)
               pre_fac(1) = r_cc(i)
               pre_fac(2) = r_cc(i)**2
               call Laplacian(cc(i, j, k, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
   
               ! nonlinear source terms
               cc(i, j, k, i_out) = cc(i, j, k, i_out) &
                       + r_cc(i)**2 * cc(i, j, k, f1) * ( 1.0d0 + cc(i, j, k, n) )
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_clalp

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> cylindrical geometry.
  subroutine box_gs_clalp(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di
    integer                   :: lo(NDIM), hi(NDIM), pm(NDIM)

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0d0
    real(dp)                  :: pre_fac(NDIM) = 1.0d0
    real(dp)                  :: u(-1:1, NDIM)

    logical                   :: redblack
    real(dp)                  :: idr2(NDIM), dr(NDIM), dr2(NDIM)
    real(dp)                  :: Lop, dLdu
    real(dp)                  :: r_face(1:nc+1)
    real(dp)                  :: r_cc(nc)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])

    where ( mg%boxes(id)%r_min + 0.5_dp * nc * mg%boxes(id)%dr.ge.0.0_dp .and. mg%symm)
      pm = 1
      lo = 1
      hi = nc
    else where
      pm = -1
      lo = nc
      hi = 1
    end where

    i0 = lo(1)
    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2 * pm(1)
    else
       di = pm(1)
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 2
      do j = lo(2), hi(2), pm(2)
         if (redblack) &
              i0 = lo(1) + iand(lo(1),1) - iand(ieor(redblack_cntr, j), 1)

         do i = i0, hi(1), di

            ! Laplacian 
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)
            pre_fac(1) = r_cc(i)
            pre_fac(2) = r_cc(i)**2
            call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
            ! nonlinear source terms
            Lop = Lop + r_cc(i)**2 * &
                       cc(i, j, f1) * ( 1.0d0 + cc(i, j, n) )
            dLdu = dLdu + r_cc(i)**2 * cc(i, j, f1) 
    
            cc(i, j, n) = cc(i, j, n) - ( Lop - cc(i, j, mg_irhs) ) / dLdu
         end do
      end do

#elif NDIM == 3
      do k = lo(3), hi(3), pm(3)
         do j = lo(2), hi(2), pm(2)
            if (redblack) &
                 i0 = lo(1) + iand(lo(1),1) - iand(ieor(redblack_cntr, j+k), 1)
   
            do i = i0, hi(1), di
   
               ! Laplacian 
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)
               pre_fac(1) = r_cc(i)
               pre_fac(2) = r_cc(i)**2
               call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
   
               ! nonlinear source terms
               Lop = Lop + r_cc(i)**2 * &
                           cc(i, j, k, f1) * ( 1.0d0 + cc(i, j, k, n) )
               dLdu = dLdu + r_cc(i)**2 * &
                           cc(i, j, k, f1) 
       
               cc(i, j, k, n) = cc(i, j, k, n) - ( Lop - cc(i, j, k, mg_irhs) ) / dLdu
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_clalp

  !> Perform L operator on a box in spherical geometry, using (r,theta,phi)
  subroutine box_slalp(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0d0
    real(dp)                  :: pre_fac(NDIM) = 1.0d0
    real(dp)                  :: u(-1:1, NDIM)

    real(dp)                  :: dr(NDIM), idr2(NDIM)
    real(dp)                  :: r_face(1:nc+1)
    real(dp)                  :: r_cc(nc)
#if NDIM != 1
    real(dp)                  :: sin_theta_face(1:nc+1)
    real(dp)                  :: sin_theta_cc(nc)
#endif

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])
#if NDIM != 1
    sin_theta_face  = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-1.0_dp, i=1,nc+1)])
    sin_theta_cc    = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-0.5_dp, i=1,nc)])
#endif

    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 1
      do i = 1, nc
         ! Laplacian
         u(-1:1, 1) = cc(i-1:i+1, n)
         face_coeff(0:1, 1) = r_face(i:i+1)**2
         call Laplacian(cc(i, i_out), u, idr2, face_coeff_in = face_coeff)
         ! nonlinear source terms
         cc(i, i_out) = cc(i, i_out) &
              + r_cc(i)**2 * cc(i, f1) * ( 1.0d0 + cc(i, n) )
      end do
#elif NDIM == 2
      do j = 1, nc
         do i = 1, nc
            ! Laplacian
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)**2
            face_coeff(0:1, 2) = sin_theta_face(j:j+1)
            pre_fac(1) = sin_theta_cc(j)
            call Laplacian(cc(i, j, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
            ! nonlinear source terms
            cc(i, j, i_out) = cc(i, j, i_out) &
                       + r_cc(i)**2 * sin_theta_cc(j) *  &
                         cc(i, j, f1) * ( 1.0d0 + cc(i, j, n) )
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               ! Laplacian
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)**2
               face_coeff(0:1, 2) = sin_theta_face(j:j+1)
               pre_fac(1) = sin_theta_cc(j)**2
               pre_fac(2) = sin_theta_cc(j)
               call Laplacian(cc(i, j, k, i_out), u, idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
   
               ! nonlinear source terms
               cc(i, j, k, i_out) = cc(i, j, k, i_out) &
                       + r_cc(i)**2 * sin_theta_cc(j)**2 *  &
                          cc(i, j, k, f1) * ( 1.0d0 + cc(i, j, k, n) )
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_slalp

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> spherical geometry.
  subroutine box_gs_slalp(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0d0
    real(dp)                  :: pre_fac(NDIM) = 1.0d0
    real(dp)                  :: u(-1:1, NDIM)

    logical                   :: redblack
    real(dp)                  :: idr2(NDIM), dr(NDIM), dr2(NDIM)
    real(dp)                  :: Lop, dLdu
    real(dp)                  :: r_face(1:nc+1)
    real(dp)                  :: r_cc(nc)
#if NDIM != 1
    real(dp)                  :: sin_theta_face(1:nc+1)
    real(dp)                  :: sin_theta_cc(nc)
#endif

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])
#if NDIM != 1
    sin_theta_face  = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-1.0_dp, i=1,nc+1)])
    sin_theta_cc  = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-0.5_dp, i=1,nc)])
#endif

    i0  = 1
    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
               f1 => mg_itmp1, f2 => mg_itmp2)
#if NDIM == 1
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)

      do i = i0, nc, di
         ! Laplacian
         u(-1:1, 1) = cc(i-1:i+1, n)
         face_coeff(0:1, 1) = r_face(i:i+1)**2
         call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff)

         ! nonlinear source terms
         Lop = Lop + r_cc(i)**2 * cc(i, f1) * ( 1.0d0 + cc(i, n) )
         dLdu = dLdu + r_cc(i)**2 * cc(i, f1)
 
         cc(i, n) = cc(i, n) - ( Lop - cc(i, mg_irhs) ) / dLdu
      end do

#elif NDIM == 2
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)

         do i = i0, nc, di

            ! Laplacian 
            u(-1:1, 1) = cc(i-1:i+1, j, n)
            u(-1:1, 2) = cc(i, j-1:j+1, n)
            face_coeff(0:1, 1) = r_face(i:i+1)**2
            face_coeff(0:1, 2) = sin_theta_face(j:j+1)
            pre_fac(1) = sin_theta_cc(j)
            call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
            ! nonlinear source terms
            Lop = Lop + r_cc(i)**2 * sin_theta_cc(j) * &
                       cc(i, j, f1) * ( 1.0d0 + cc(i, j, n) )
            dLdu = dLdu + r_cc(i)**2 * sin_theta_cc(j) * cc(i, j, f1) 
    
            cc(i, j, n) = cc(i, j, n) - ( Lop - cc(i, j, mg_irhs) ) / dLdu
         end do
      end do

#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, j+k), 1)
   
            do i = i0, nc, di
   
               ! Laplacian 
               u(-1:1, 1) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)**2
               face_coeff(0:1, 2) = sin_theta_face(j:j+1)
               pre_fac(1) = sin_theta_cc(j)**2
               pre_fac(2) = sin_theta_cc(j)
               call Laplacian(Lop, u, idr2, dLdu = dLdu, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
   
               ! nonlinear source terms
               Lop = Lop + r_cc(i)**2 * sin_theta_cc(j)**2 * &
                           cc(i, j, k, f1) * ( 1.0d0 + cc(i, j, k, n) )
               dLdu = dLdu + r_cc(i)**2 * sin_theta_cc(j)**2 * &
                           cc(i, j, k, f1) 
       
               cc(i, j, k, n) = cc(i, j, k, n) - ( Lop - cc(i, j, k, mg_irhs) ) / dLdu
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_slalp

end module m_cfc_alp
