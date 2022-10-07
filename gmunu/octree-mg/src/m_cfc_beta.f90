#include "cpp_macros.h"
!> Module which contains multigrid procedures for the L operator for (beta)
module m_cfc_beta
  use m_data_structures

  implicit none
  private

  public :: cfc_beta_set_methods

contains

  subroutine cfc_beta_set_methods(mg)
    type(mg_t), intent(inout) :: mg
    
    mg%vector_equation = .true.

    select case (mg%geometry_type)
#if NDIM == 3
    case (mg_cartesian)
       mg%box_vec_op => box_lbeta

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_lbeta
       case default
          error stop "cfc_beta_set_methods: unsupported smoother type"
       end select
#endif
#if NDIM != 1
    case (mg_cylindrical)
       mg%box_vec_op => box_clbeta

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_clbeta
       case default
          error stop "cfc_beta_set_methods: unsupported smoother type"
       end select
#endif
    case (mg_spherical)
       mg%box_vec_op => box_slbeta

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_slbeta
       case default
          error stop "cfc_beta_set_methods: unsupported smoother type"
       end select
    case default
       error stop "cfc_beta_set_methods: unsupported geometry"
    end select

  end subroutine cfc_beta_set_methods

#if NDIM == 3
  !> Perform L operator on a box in cartesian geometry
  subroutine box_lbeta(mg, id, nc, i_dir, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_dir
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK

    real(dp)                  :: u(-1:1, NDIM, 3)
    real(dp)                  :: du(3, NDIM)
    real(dp)                  :: d2u(3, NDIM, NDIM)

    real(dp)                  :: idr2(NDIM), dr(NDIM)
    real(dp)                  :: Lop(3), dLdu(3)
    real(dp)                  :: u_tmp(-1:1)
    real(dp)                  :: u2_tmp(-1:1,-1:1)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2

    associate (cc => mg%boxes(id)%cc, n => mg_vec_iphi(i_dir))
    ! Laplacian term 
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             u(-1:1, 1, i_dir) = cc(i-1:i+1, j, k, n)
             u(-1:1, 2, i_dir) = cc(i, j-1:j+1, k, n)
             u(-1:1, 3, i_dir) = cc(i, j, k-1:k+1, n)
             call Laplacian(cc(i, j, k, i_out), u(:,:,i_dir), idr2)
          end do
       end do
    end do

    select case (i_dir)
    case (1)
       do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u_tmp(-1:1) = cc(i-1:i+1, j, k, mg_vec_iphi(i_dir))
               call derivatives_2nd( d2u(i_dir,i_dir,i_dir), u_tmp, idr2(i_dir) )
               u2_tmp(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(2))
               call mixed_derivatives_2nd(d2u(2,i_dir,2), u2_tmp, dr(i_dir), dr(2), i, j, nc)
               u2_tmp(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(3))
               call mixed_derivatives_2nd(d2u(3,i_dir,3), u2_tmp, dr(i_dir), dr(3), i, k, nc)

               ! div_X term
               cc(i, j, k, i_out) = cc(i, j, k, i_out) + 1.0_dp/3.0_dp * (  &
                  d2u(1,i_dir,1) + d2u(2,i_dir,2) + d2u(3,i_dir,3)      )
            end do
         end do
       end do

    case (2)
       do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u_tmp(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(i_dir))
               call derivatives_2nd( d2u(i_dir,i_dir,i_dir), u_tmp, idr2(i_dir) )
               u2_tmp(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(1))
               call mixed_derivatives_2nd(d2u(1,i_dir,1), u2_tmp, dr(i_dir), dr(1), j, i, nc)
               u2_tmp(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(3))
               call mixed_derivatives_2nd(d2u(3,i_dir,3), u2_tmp, dr(i_dir), dr(3), j, k, nc)

               ! div_X term
               cc(i, j, k, i_out) = cc(i, j, k, i_out) + 1.0_dp/3.0_dp * (  &
                  d2u(1,i_dir,1) + d2u(2,i_dir,2) + d2u(3,i_dir,3)     )
            end do
         end do
       end do

    case (3)
       do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u_tmp(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(i_dir))
               call derivatives_2nd( d2u(i_dir,i_dir,i_dir), u_tmp, idr2(i_dir) )
               u2_tmp(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(1))
               call mixed_derivatives_2nd(d2u(1,i_dir,1), u2_tmp, dr(i_dir), dr(1), k, i, nc)
               u2_tmp(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(2))
               call mixed_derivatives_2nd(d2u(2,i_dir,2), u2_tmp, dr(i_dir), dr(2), k, j, nc)

               ! div_X term
               cc(i, j, k, i_out) = cc(i, j, k, i_out) + 1.0_dp/3.0_dp * (  &
                  d2u(1,i_dir,1) + d2u(2,i_dir,2) + d2u(3,i_dir,3)     )
            end do
         end do
       end do
    case default
     cc(1:nc, 1:nc, 1:nc, i_out) = 0.0_dp
    end select
    end associate
  end subroutine box_lbeta

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> spherical geometry.
  subroutine box_gs_lbeta(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di, n
    integer                   :: lo(NDIM), hi(NDIM), pm(NDIM)

    real(dp)                  :: u(-1:1, NDIM, 3)
    real(dp)                  :: u_tmp_1(-1:1)
    real(dp)                  :: du(3, NDIM)
    real(dp)                  :: d2u(3, NDIM, NDIM)
    real(dp)                  :: dLdu_d2u(3, NDIM, NDIM)

    logical                   :: redblack
    real(dp)                  :: idr2(NDIM), dr(NDIM)
    real(dp)                  :: Lop(3), dLdu(3)
    integer                   :: i_dir
    real(dp)                  :: u_tmp(-1:1)
    real(dp)                  :: u2_tmp(-1:1,-1:1)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2

    where ( mg%boxes(id)%r_min + 0.5_dp * nc * mg%boxes(id)%dr.ge.0.0_dp .and. mg%symm )
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
    associate (cc => mg%boxes(id)%cc)
    do k = lo(3), hi(3), pm(3)
      do j = lo(2), hi(2), pm(2)
         if (redblack) &
              i0 = lo(1) + iand(lo(1),1) - iand(ieor(redblack_cntr, j+k), 1)
         do i = i0, hi(1), di

            do i_dir = 1, 3
               n = mg_vec_iphi(i_dir)
               ! Laplacian term
               u(-1:1, 1, i_dir) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2, i_dir) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3, i_dir) = cc(i, j, k-1:k+1, n)
               call Laplacian( Lop(i_dir), u(:,:,i_dir), idr2, dLdu = dLdu(i_dir) )

               select case (i_dir)
               case (1)
                  u_tmp(-1:1) = cc(i-1:i+1, j, k, mg_vec_iphi(i_dir))
                  call derivatives_2nd( d2u(i_dir,i_dir,i_dir), u_tmp, idr2(i_dir), dLdu=dLdu_d2u(i_dir,i_dir,i_dir))
                  u2_tmp(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(2))
                  call mixed_derivatives_2nd(d2u(2,i_dir,2), u2_tmp, dr(i_dir), dr(2), i, j, nc, dLdu=dLdu_d2u(2,i_dir,2))
                  u2_tmp(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(3))
                  call mixed_derivatives_2nd(d2u(3,i_dir,3), u2_tmp, dr(i_dir), dr(3), i, k, nc, dLdu=dLdu_d2u(3,i_dir,3))
               case (2)
                  u_tmp(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(i_dir))
                  call derivatives_2nd( d2u(i_dir,i_dir,i_dir), u_tmp, idr2(i_dir), dLdu=dLdu_d2u(i_dir,i_dir,i_dir))
                  u2_tmp(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(1))
                  call mixed_derivatives_2nd(d2u(1,i_dir,1), u2_tmp, dr(i_dir), dr(1), j, i, nc, dLdu=dLdu_d2u(1,i_dir,1))
                  u2_tmp(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(3))
                  call mixed_derivatives_2nd(d2u(3,i_dir,3), u2_tmp, dr(i_dir), dr(3), j, k, nc, dLdu=dLdu_d2u(3,i_dir,3))
               case (3)
                  u_tmp(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(i_dir))
                  call derivatives_2nd( d2u(i_dir,i_dir,i_dir), u_tmp, idr2(i_dir), dLdu=dLdu_d2u(i_dir,i_dir,i_dir))
                  u2_tmp(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(1))
                  call mixed_derivatives_2nd(d2u(1,i_dir,1), u2_tmp, dr(i_dir), dr(1), k, i, nc, dLdu=dLdu_d2u(1,i_dir,1))
                  u2_tmp(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(2))
                  call mixed_derivatives_2nd(d2u(2,i_dir,2), u2_tmp, dr(i_dir), dr(2), k, j, nc, dLdu=dLdu_d2u(2,i_dir,2))
               case default
                 Lop(i_dir) = 0.0_dp
                 dLdu(i_dir) = 1.0_dp
               end select
   
               ! div_X term
               Lop(i_dir) = Lop(i_dir) &
                  + 1.0_dp/3.0_dp * (   d2u(1,i_dir,1) + d2u(2,i_dir,2) + d2u(3,i_dir,3)    )
               ! div_X for dLdu
               dLdu(i_dir) = dLdu(i_dir) &
                  + 1.0_dp/3.0_dp * (  dLdu_d2u(1,i_dir,1) + dLdu_d2u(2,i_dir,2) + dLdu_d2u(3,i_dir,3)   )

            end do ! end i_dir

            do i_dir = 1, 3
               cc(i, j, k, mg_vec_iphi(i_dir)) = cc(i, j, k, mg_vec_iphi(i_dir)) &
                                       - (Lop(i_dir) - cc(i, j, k, mg_vec_irhs(i_dir))) / dLdu(i_dir)
            end do ! end i_dir

         end do ! end i
      end do ! end j
    end do ! end k
    end associate
  end subroutine box_gs_lbeta
#endif

  !> Perform L operator on a box in cylindrical geometry, using (r,z,phi)
  subroutine box_clbeta(mg, id, nc, i_dir, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_dir
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0_dp
    real(dp)                  :: pre_fac(NDIM) = 1.0_dp
    real(dp)                  :: u(-1:1, NDIM, 3)
    real(dp)                  :: u_tmp_1(-1:1)
    real(dp)                  :: du(3, NDIM)
    real(dp)                  :: d2u(3, NDIM, NDIM)

    real(dp)                  :: idr2(NDIM), dr(NDIM)
    real(dp)                  :: r_face(1:nc+1), r_cc(nc)
    real(dp)                  :: Lop(3), dLdu(3)
    real(dp)                  :: u_tmp_2(-1:1,-1:1)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])

    associate (cc => mg%boxes(id)%cc, n => mg_vec_iphi(i_dir))
#if NDIM == 2
    ! Laplacian term * r**2 
    do j = 1, nc
       do i = 1, nc
          u(-1:1, 1, i_dir) = cc(i-1:i+1, j, n)
          u(-1:1, 2, i_dir) = cc(i, j-1:j+1, n)
          face_coeff(0:1, 1) = r_face(i:i+1)
          pre_fac(1) = r_cc(i)
          pre_fac(2) = r_cc(i)**2
          call Laplacian(cc(i, j, i_out), u(:,:,i_dir), idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
       end do
    end do

    select case (i_dir)
    case (1)
      do j = 1, nc
         do i = 1, nc
            u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, mg_vec_iphi(2))
            call mixed_derivatives_2nd(d2u(2,1,2), u_tmp_2, dr(1), dr(2), i, j, nc)
            u_tmp_1(-1:1) = cc(i-1:i+1, j, mg_vec_iphi(1))
            du(1,1) = dudr(u_tmp_1, dr(1))
            call derivatives_2nd(d2u(1,1,1), u_tmp_1, idr2(1))

            ! extra Laplacian term for vectors
            cc(i, j, i_out) = cc(i, j, i_out) &
                - cc(i, j, mg_vec_iphi(1))  

            ! div_X term
            cc(i, j, i_out) = cc(i, j, i_out) + &
                1.0_dp/3.0_dp * ( du(1,1) * r_cc(i) - cc(i, j, mg_vec_iphi(1)) &
                 + r_cc(i)**2 * ( d2u(1,1,1) + d2u(2,1,2) )  &
               )
         end do
      end do
    case (2)
      do j = 1, nc
         do i = 1, nc
            u_tmp_1(-1:1) = cc(i, j-1:j+1, mg_vec_iphi(1))
            du(1,2) = dudr(u_tmp_1, dr(2))
            u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, mg_vec_iphi(1))
            call mixed_derivatives_2nd(d2u(1,1,2), u_tmp_2, dr(1), dr(2), i, j, nc)
            u_tmp_1(-1:1) = cc(i, j-1:j+1, mg_vec_iphi(2))
            call derivatives_2nd(d2u(2,2,2), u_tmp_1, idr2(2))

            ! div_X term
            cc(i, j, i_out) = cc(i, j, i_out) + &
                1.0_dp/3.0_dp * ( du(1,2) * r_cc(i) +  &
                 r_cc(i)**2 * ( d2u(1,1,2) + d2u(2,2,2) ) &
                ) 
         end do
      end do
    case (3)
      do j = 1, nc
         do i = 1, nc
            ! extra Laplacian term for vectors
            cc(i, j, i_out) = cc(i, j, i_out) &
                - cc(i, j, mg_vec_iphi(3)) 
         end do
      end do
     case default
      cc(DTIMES(1:nc), i_out) = 0.0_dp
     end select
#elif NDIM == 3
    ! Laplacian term 
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             u(-1:1, 1, i_dir) = cc(i-1:i+1, j, k, n)
             u(-1:1, 2, i_dir) = cc(i, j-1:j+1, k, n)
             u(-1:1, 3, i_dir) = cc(i, j, k-1:k+1, n)
             face_coeff(0:1, 1) = r_face(i:i+1)
             pre_fac(1) = r_cc(i)
             pre_fac(2) = r_cc(i)**2
             call Laplacian(cc(i, j, k, i_out), u(:,:,i_dir), idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
          end do
       end do
    end do

    select case (i_dir)
    case (1)
       do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(2))
               call mixed_derivatives_2nd(d2u(2,1,2), u_tmp_2, dr(1), dr(2), i, j, nc)
               u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(3))
               call mixed_derivatives_2nd(d2u(3,1,3), u_tmp_2, dr(1), dr(3), i, k, nc)
               u_tmp_1(-1:1) = cc(i-1:i+1, j, k, mg_vec_iphi(1))
               du(1,1) = dudr(u_tmp_1, dr(1))
               u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(3))
               du(3,3) = dudr(u_tmp_1, dr(3))
               call derivatives_2nd(d2u(1,1,1), u_tmp_1, idr2(1))
   
               ! extra Laplacian term for vectors
               cc(i, j, k, i_out) = cc(i, j, k, i_out) &
                   - cc(i, j, k, mg_vec_iphi(1)) - 2.0_dp * du(3,3)  
   
               ! div_X term
               cc(i, j, k, i_out) = cc(i, j, k, i_out) + &
                   1.0_dp/3.0_dp * ( du(1,1) * r_cc(i) - cc(i, j, k, mg_vec_iphi(1)) - du(3,3) &
                    + r_cc(i)**2 * ( d2u(1,1,1) + d2u(2,1,2)  + d2u(3,1,3) )  &
                  )
            end do
         end do
       end do

    case (2)
       do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u_tmp_1(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(1))
               du(1,2) = dudr(u_tmp_1, dr(2))
               u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(1))
               call mixed_derivatives_2nd(d2u(1,1,2), u_tmp_2, dr(1), dr(2), i, j, nc)
               u_tmp_2(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(3))
               call mixed_derivatives_2nd(d2u(3,2,3), u_tmp_2, dr(2), dr(3), j, k, nc)
               u_tmp_1(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(2))
               call derivatives_2nd(d2u(2,2,2), u_tmp_1, idr2(2))
   
               ! div_X term
               cc(i, j, k, i_out) = cc(i, j, k, i_out) + &
                   1.0_dp/3.0_dp * ( du(1,2) * r_cc(i) +  &
                    r_cc(i)**2 * ( d2u(1,1,2) + d2u(2,2,2)  + d2u(3,2,3) ) &
                   ) 
            end do
         end do
       end do

    case (3)
       do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(1))
               du(1,3) = dudr(u_tmp_1, dr(3))
               u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(1))
               call mixed_derivatives_2nd(d2u(1,1,3), u_tmp_2, dr(1), dr(3), i, k, nc)
               u_tmp_2(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(2))
               call mixed_derivatives_2nd(d2u(2,2,3), u_tmp_2, dr(2), dr(3), j, k, nc)
               u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(3))
               call derivatives_2nd(d2u(3,3,3), u_tmp_1, idr2(3))

               ! extra Laplacian term for vectors
               cc(i, j, k, i_out) = cc(i, j, k, i_out) &
                   - cc(i, j, k, mg_vec_iphi(3)) + 2.0_dp * du(1,3)  
   
               ! div_X term
               cc(i, j, k, i_out) = cc(i, j, k, i_out) + &
                   1.0_dp/3.0_dp * ( du(1,3) + &
                    r_cc(i) * ( d2u(1,1,3) + d2u(3,3,3)  + d2u(2,2,3) ) &
                   ) 
            end do
         end do
       end do
    case default
      cc(DTIMES(1:nc), i_out) = 0.0_dp
    end select
#endif
    end associate
  end subroutine box_clbeta

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> cylindrical geometry.
  subroutine box_gs_clbeta(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di, n
    integer                   :: lo(NDIM), hi(NDIM), pm(NDIM)

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0_dp
    real(dp)                  :: pre_fac(NDIM) = 1.0_dp
    real(dp)                  :: u(-1:1, NDIM, 3)
    real(dp)                  :: u_tmp_1(-1:1)
    real(dp)                  :: du(3, NDIM)
    real(dp)                  :: d2u(3, NDIM, NDIM)
    real(dp)                  :: dLdu_d2u(3, NDIM, NDIM)

    logical                   :: redblack
    real(dp)                  :: idr2(NDIM), dr(NDIM)
    real(dp)                  :: r_face(1:nc+1), r_cc(nc)
    real(dp)                  :: Lop(3), dLdu(3)
    integer                   :: i_dir
    real(dp)                  :: u_tmp_2(-1:1,-1:1)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])

    where ( mg%boxes(id)%r_min + 0.5_dp * nc * mg%boxes(id)%dr.ge.0.0_dp .and. mg%symm )
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
    associate (cc => mg%boxes(id)%cc)
#if NDIM == 2
      do j = lo(2), hi(2), pm(2)
         if (redblack) &
              i0 = lo(1) + iand(lo(1),1) - iand(ieor(redblack_cntr, j), 1)
         do i = i0, hi(1), di
            do i_dir = 1, 3
               n = mg_vec_iphi(i_dir)
               ! Laplacian term
               u(-1:1, 1, i_dir) = cc(i-1:i+1, j, n)
               u(-1:1, 2, i_dir) = cc(i, j-1:j+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)
               pre_fac(1) = r_cc(i)
               pre_fac(2) = r_cc(i)**2
               call Laplacian(Lop(i_dir), u(:,:,i_dir), idr2, dLdu = dLdu(i_dir), face_coeff_in = face_coeff, pre_fac_in = pre_fac)

               select case (i_dir)
               case (1)
                 u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, mg_vec_iphi(2))
                 call mixed_derivatives_2nd(d2u(2,1,2), u_tmp_2, dr(1), dr(2), i, j, nc, dLdu=dLdu_d2u(2,1,2))
                 u_tmp_1(-1:1) = cc(i-1:i+1, j, mg_vec_iphi(1))
                 du(1,1) = dudr(u_tmp_1, dr(1))
                 call derivatives_2nd(d2u(1,1,1), u_tmp_1, idr2(1), &
                       dLdu=dLdu_d2u(1,1,1))

                 ! extra Laplacian term for vectors
                 Lop(i_dir) = Lop(i_dir) - cc(i, j, mg_vec_iphi(1)) 
                 dLdu(i_dir) = dLdu(i_dir) - 1.0_dp 

                 ! div_X term
                 Lop(i_dir) = Lop(i_dir) + &
                     1.0_dp/3.0_dp * ( du(1,1) * r_cc(i) - cc(i, j, mg_vec_iphi(1)) &
                     + r_cc(i)**2 * ( d2u(1,1,1) + d2u(2,1,2) )  &
                    )
                 ! div_X for dLdu
                 dLdu(i_dir) = dLdu(i_dir) + &
                  1.0_dp/3.0_dp * ( r_cc(i)**2 * ( dLdu_d2u(1,1,1) + dLdu_d2u(2,1,2) ) - 1.0_dp )  
               case (2)
                 u_tmp_1(-1:1) = cc(i, j-1:j+1, mg_vec_iphi(1))
                 du(1,2) = dudr(u_tmp_1, dr(2))
                 u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, mg_vec_iphi(1))
                 call mixed_derivatives_2nd(d2u(1,1,2), u_tmp_2, dr(1), dr(2), i, j, nc, dLdu=dLdu_d2u(1,1,2))
                 u_tmp_1(-1:1) = cc(i, j-1:j+1, mg_vec_iphi(2))
                 call derivatives_2nd(d2u(2,2,2), u_tmp_1, idr2(2), &
                      dLdu=dLdu_d2u(2,2,2))

                 ! div_X term
                 Lop(i_dir) = Lop(i_dir) + &
                      1.0_dp/3.0_dp * ( du(1,2) * r_cc(i) + &
                       r_cc(i)**2 * ( d2u(1,1,2) + d2u(2,2,2) ) &
                      ) 

                 ! div_X for dLdu
                 dLdu(i_dir) = dLdu(i_dir)&
                     + 1.0_dp/3.0_dp *  &
                         r_cc(i)**2 * ( dLdu_d2u(1,1,2) + dLdu_d2u(2,2,2) )
               case (3)
                 ! extra Laplacian term for vectors
                 Lop(i_dir) = Lop(i_dir) - cc(i, j, mg_vec_iphi(3)) 
                 dLdu(i_dir) = dLdu(i_dir) - 1.0_dp
               case default
                 Lop(i_dir) = 0.0_dp
                 dLdu(i_dir) = 1.0_dp
               end select
            end do ! end i_dir

            do i_dir = 1, 3
               cc(i,j, mg_vec_iphi(i_dir)) = cc(i,j, mg_vec_iphi(i_dir)) &
                                       - (Lop(i_dir) - cc(i,j, mg_vec_irhs(i_dir))) / dLdu(i_dir)
            end do ! end i_dir

         end do ! end i
      end do ! end j

#elif NDIM == 3
    do k = lo(3), hi(3), pm(3)
      do j = lo(2), hi(2), pm(2)
         if (redblack) &
              i0 = lo(1) + iand(lo(1),1) - iand(ieor(redblack_cntr, j+k), 1)
         do i = i0, hi(1), di

            do i_dir = 1, 3
               n = mg_vec_iphi(i_dir)
               ! Laplacian term
               u(-1:1, 1, i_dir) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2, i_dir) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3, i_dir) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)
               pre_fac(1) = r_cc(i)
               pre_fac(2) = r_cc(i)**2
               call Laplacian( Lop(i_dir), u(:,:,i_dir), idr2, dLdu = dLdu(i_dir) , face_coeff_in = face_coeff, pre_fac_in = pre_fac)

               select case (i_dir)
               case (1)
                  u_tmp_1(-1:1) = cc(i-1:i+1, j, k, mg_vec_iphi(1))
                  du(1,1) = dudr(u_tmp_1, dr(1))
                  u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(3))
                  du(3,3) = dudr(u_tmp_1, dr(3))
                  u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(2))
                  call mixed_derivatives_2nd(d2u(2,1,2), u_tmp_2, dr(1), dr(2), i, j, nc, dLdu=dLdu_d2u(2,1,2))
                  u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(3))
                  call mixed_derivatives_2nd(d2u(3,1,3), u_tmp_2, dr(1), dr(3), i, k, nc, dLdu=dLdu_d2u(3,1,3))
                  call derivatives_2nd(d2u(1,1,1), u_tmp_1, idr2(1), dLdu=dLdu_d2u(1,1,1))
      
                  ! extra Laplacian term for vectors
                  Lop(i_dir) = Lop(i_dir) &
                      - cc(i, j, k, mg_vec_iphi(1)) - 2.0_dp * du(3,3)  
                  dLdu(i_dir) = dLdu(i_dir) - 1.0_dp 
      
                  ! div_X term
                  Lop(i_dir) = Lop(i_dir) + &
                      1.0_dp/3.0_dp * ( du(1,1) * r_cc(i) - cc(i, j, k, mg_vec_iphi(1)) - du(3,3) &
                       + r_cc(i)**2 * ( d2u(1,1,1) + d2u(2,1,2)  + d2u(3,1,3) )  &
                     )
                  ! div_X for dLdu
                  dLdu(i_dir) = dLdu(i_dir) + &
                      1.0_dp/3.0_dp * ( r_cc(i)**2 * ( dLdu_d2u(1,1,1) + dLdu_d2u(2,1,2) + dLdu_d2u(3,1,3) ) - 1.0_dp )  

               case (2)
                 u_tmp_1(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(1))
                 du(1,2) = dudr(u_tmp_1, dr(2))
                 u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(1))
                 call mixed_derivatives_2nd(d2u(1,1,2), u_tmp_2, dr(1), dr(2), i, j, nc, dLdu=dLdu_d2u(1,1,2))
                 u_tmp_2(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(3))
                 call mixed_derivatives_2nd(d2u(3,2,3), u_tmp_2, dr(2), dr(3), j, k, nc, dLdu=dLdu_d2u(3,2,3))
                 u_tmp_1(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(2))
                 call derivatives_2nd(d2u(2,2,2), u_tmp_1, idr2(2), dLdu=dLdu_d2u(2,2,2))

                 ! div_X term
                 Lop(i_dir) = Lop(i_dir) + &
                   1.0_dp/3.0_dp * ( du(1,2) * r_cc(i) +  &
                    r_cc(i)**2 * ( d2u(1,1,2) + d2u(2,2,2)  + d2u(3,2,3) ) &
                   ) 

                 ! div_X for dLdu
                 dLdu(i_dir) = dLdu(i_dir)&
                     + 1.0_dp/3.0_dp *  &
                         r_cc(i)**2 * ( dLdu_d2u(1,1,2) + dLdu_d2u(2,2,2) + dLdu_d2u(3,2,3) )

               case (3)
                 u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(1))
                 du(1,3) = dudr(u_tmp_1, dr(3))
                 u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(1))
                 call mixed_derivatives_2nd(d2u(1,1,3), u_tmp_2, dr(1), dr(3), i, k, nc, dLdu=dLdu_d2u(1,1,3))
                 u_tmp_2(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(2))
                 call mixed_derivatives_2nd(d2u(2,2,3), u_tmp_2, dr(2), dr(3), j, k, nc, dLdu=dLdu_d2u(2,2,3))
                 u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(3))
                 call derivatives_2nd(d2u(3,3,3), u_tmp_1, idr2(3), dLdu=dLdu_d2u(3,3,3))

                 ! extra Laplacian term for vectors
                 Lop(i_dir) = Lop(i_dir) &
                      - cc(i, j, k, mg_vec_iphi(3)) - 2.0_dp * du(1,3)  
                 dLdu(i_dir) = dLdu(i_dir) - 1.0_dp 

                 ! div_X term
                 Lop(i_dir) = Lop(i_dir) + &
                   1.0_dp/3.0_dp * ( du(1,3) +  &
                    r_cc(i) * ( d2u(1,1,3) + d2u(3,3,3)  + d2u(2,2,3) ) &
                   ) 

                 ! div_X for dLdu
                 dLdu(i_dir) = dLdu(i_dir)&
                     + 1.0_dp/3.0_dp *  &
                         r_cc(i) * ( dLdu_d2u(1,1,3) + dLdu_d2u(2,2,3) + dLdu_d2u(3,3,3) )

               case default
                 Lop(i_dir) = 0.0_dp
                 dLdu(i_dir) = 1.0_dp
               end select
            end do ! end i_dir

            do i_dir = 1, 3
               cc(i, j, k, mg_vec_iphi(i_dir)) = cc(i, j, k, mg_vec_iphi(i_dir)) &
                                       - (Lop(i_dir) - cc(i, j, k, mg_vec_irhs(i_dir))) / dLdu(i_dir)
            end do ! end i_dir

         end do ! end i
      end do ! end j
    end do ! end k
#endif
    end associate
  end subroutine box_gs_clbeta

  !> Perform L operator on a box in spherical geometry, using (r,theta,phi)
  subroutine box_slbeta(mg, id, nc, i_dir, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_dir
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0_dp
    real(dp)                  :: pre_fac(NDIM) = 1.0_dp
    real(dp)                  :: u(-1:1, NDIM, 3)
    real(dp)                  :: u_tmp_1(-1:1)
    real(dp)                  :: du(3, NDIM)
    real(dp)                  :: d2u(3, NDIM, NDIM)

    real(dp)                  :: idr2(NDIM), dr(NDIM)
    real(dp)                  :: r_face(1:nc+1), r_cc(0:nc+1)
    real(dp)                  :: Lop(3), dLdu(3)
#if NDIM != 1
    real(dp)                  :: u_tmp_2(-1:1,-1:1)
    real(dp)                  :: sin_theta_cc(0:nc+1)
    real(dp)                  :: sin_theta_face(1:nc+1)
    real(dp)                  :: cos_theta_cc(0:nc+1)
#endif

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=0,nc+1)])
#if NDIM != 1
    sin_theta_cc    = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-0.5_dp, i=0,nc+1)])
    sin_theta_face  = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-1.0_dp, i=1,nc+1)])
    cos_theta_cc    = dcos(mg%boxes(id)%r_min(2) + dr(2) * [(i-0.5_dp, i=0,nc+1)]) 
#endif

    associate (cc => mg%boxes(id)%cc, n => mg_vec_iphi(i_dir))
#if NDIM == 1
    select case (i_dir)
    case (1)
      do i = 1, nc
         ! Laplacian
         u(-1:1, 1, 1) = cc(i-1:i+1, n)
         face_coeff(0:1, 1) = r_face(i:i+1)**2
         call Laplacian(cc(i, i_out), u(:,:,1), idr2, face_coeff_in = face_coeff)

         ! source terms
         cc(i, i_out) = 4.0_dp/3.0_dp * ( cc(i, i_out) - 2.0_dp * cc(i, n) )
      end do
     case default
      cc(1:nc, i_out) = 0.0_dp
     end select
#elif NDIM == 2
    ! Laplacian term * r**2 sin**2(theta)
    do j = 1, nc
       do i = 1, nc
          u(-1:1, 1, i_dir) = cc(i-1:i+1, j, n)
          u(-1:1, 2, i_dir) = cc(i, j-1:j+1, n)
          face_coeff(0:1, 1) = r_face(i:i+1)**2
          face_coeff(0:1, 2) = sin_theta_face(j:j+1)
          ! note: here is different from the scalar cases
          pre_fac(1) = sin_theta_cc(j)**2
          pre_fac(2) = sin_theta_cc(j)
          call Laplacian(cc(i, j, i_out), u(:,:,i_dir), idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
       end do
    end do

    select case (i_dir)
    case (1)
      do j = 1, nc
         do i = 1, nc
            u_tmp_1(-1:1) = cc(i-1:i+1, j, mg_vec_iphi(2))
            du(2,1) = dudr(u_tmp_1, dr(1))
            u_tmp_1(-1:1) = cc(i, j-1:j+1, mg_vec_iphi(2))
            du(2,2) = dudr(u_tmp_1, dr(2))
            u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, mg_vec_iphi(2))
            call mixed_derivatives_2nd(d2u(2,1,2), u_tmp_2, dr(1), dr(2), i, j, nc)
            u_tmp_1(-1:1) = cc(i-1:i+1, j, mg_vec_iphi(1))
            face_coeff(0:1, 1) = r_face(i:i+1)**2
            call derivatives_2nd(d2u(1,1,1), u_tmp_1, idr2(1), face_coeff_in=face_coeff(0:1, 1))

            ! extra Laplacian term for vectors
            u_tmp_1(-1:1) = sin_theta_cc(j-1:j+1) * cc(i, j-1:j+1, mg_vec_iphi(2))
            cc(i, j, i_out) = cc(i, j, i_out) - 2.0_dp * ( &
                + sin_theta_cc(j)**2 * cc(i, j, mg_vec_iphi(1)) &
                + sin_theta_cc(j) * dudr(u_tmp_1, dr(2)) ) 

            ! div_X term
            cc(i, j, i_out) = cc(i, j, i_out) + &
                1.0_dp/3.0_dp * ( sin_theta_cc(j)**2 * ( d2u(1,1,1) &
               - 2.0_dp * cc(i, j, mg_vec_iphi(1)) - du(2,2) )   &
               - sin_theta_cc(j) * cos_theta_cc(j) * cc(i, j, mg_vec_iphi(2))  &
               + r_cc(i) * ( sin_theta_cc(j)**2 * d2u(2,1,2) &
                          + sin_theta_cc(j) * cos_theta_cc(j) * du(2,1) ) &
               )
         end do
      end do
    case (2)
      do j = 1, nc
         do i = 1, nc
            u_tmp_1(-1:1) = cc(i, j-1:j+1, mg_vec_iphi(1))
            du(1,2) = dudr(u_tmp_1, dr(2))
            u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, mg_vec_iphi(1))
            call mixed_derivatives_2nd(d2u(1,1,2), u_tmp_2, dr(1), dr(2), i, j, nc)
            u_tmp_1(-1:1) = cc(i, j-1:j+1, mg_vec_iphi(2))
            face_coeff(0:1, 1) = sin_theta_face(j:j+1)
            call derivatives_2nd(d2u(2,2,2), u_tmp_1, idr2(2), face_coeff_in=face_coeff(0:1, 1))

            ! extra Laplacian term for vectors
            cc(i, j, i_out) = cc(i, j, i_out) + ( &
                + 2.0_dp * sin_theta_cc(j)**2 * du(1,2)  &
                - cc(i, j, mg_vec_iphi(2)) &
                 ) 

            ! div_X term
            cc(i, j, i_out) = cc(i, j, i_out) &
                + 1.0_dp/3.0_dp * ( &
                   sin_theta_cc(j)**2 * ( r_cc(i) * d2u(1,1,2) + 2.0_dp * du(1,2) ) &
                 + sin_theta_cc(j) * d2u(2,2,2) - cc(i, j, mg_vec_iphi(2))  &
                   ) 
         end do
      end do
    case (3)
      do j = 1, nc
         do i = 1, nc
            ! extra Laplacian term for vectors
            cc(i, j, i_out) = cc(i, j, i_out) &
                - cc(i, j, mg_vec_iphi(3)) 
         end do
      end do
     case default
      cc(1:nc, 1:nc, i_out) = 0.0_dp
     end select
#elif NDIM == 3
    ! Laplacian term * r**2 sin**2(theta)
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             u(-1:1, 1, i_dir) = cc(i-1:i+1, j, k, n)
             u(-1:1, 2, i_dir) = cc(i, j-1:j+1, k, n)
             u(-1:1, 3, i_dir) = cc(i, j, k-1:k+1, n)
             face_coeff(0:1, 1) = r_face(i:i+1)**2
             face_coeff(0:1, 2) = sin_theta_face(j:j+1)
             ! note: here is different from the scalar cases
             pre_fac(1) = sin_theta_cc(j)**2
             pre_fac(2) = sin_theta_cc(j)
             call Laplacian(cc(i, j, k, i_out), u(:,:,i_dir), idr2, face_coeff_in = face_coeff, pre_fac_in = pre_fac)
          end do
       end do
    end do

    select case (i_dir)
    case (1)
       do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u_tmp_1(-1:1) = cc(i-1:i+1, j, k, mg_vec_iphi(2))
               du(2,1) = dudr(u_tmp_1, dr(1))
               u_tmp_1(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(2))
               du(2,2) = dudr(u_tmp_1, dr(2))
               u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(2))
               du(2,3) = dudr(u_tmp_1, dr(3))
               u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(3))
               du(3,3) = dudr(u_tmp_1, dr(3))
               u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(2))
               call mixed_derivatives_2nd(d2u(2,1,2), u_tmp_2, dr(1), dr(2), i, j, nc)
               u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(3))
               call mixed_derivatives_2nd(d2u(3,1,3), u_tmp_2, dr(1), dr(3), i, k, nc)

               u_tmp_1(-1:1) = cc(i-1:i+1, j, k, mg_vec_iphi(1))
               face_coeff(0:1, 1) = r_face(i:i+1)**2
               call derivatives_2nd(d2u(1,1,1), u_tmp_1, idr2(1), face_coeff_in=face_coeff(0:1, 1))
   
               u_tmp_1(-1:1) = sin_theta_cc(j-1:j+1) * cc(i, j-1:j+1, k, mg_vec_iphi(2))
               ! extra Laplacian term for vectors
               cc(i, j, k, i_out) = cc(i, j, k, i_out) &
                    - 2.0_dp * ( &
                   + sin_theta_cc(j)**2 * cc(i, j, k, mg_vec_iphi(1)) &
                   + sin_theta_cc(j) * ( dudr(u_tmp_1, dr(2)) + du(3,3) ) &
                              ) 
   
               ! div_X term
               cc(i, j, k, i_out) = cc(i, j, k, i_out) + &
                   1.0_dp/3.0_dp * ( &
                     sin_theta_cc(j)**2 * ( d2u(1,1,1) - 2.0_dp * cc(i, j, k, mg_vec_iphi(1)) - du(2,2) + r_cc(i) * d2u(2,1,2) )   &
                       + sin_theta_cc(j) * ( cos_theta_cc(j) * ( r_cc(i) * du(2,1) - cc(i, j, k, mg_vec_iphi(2)) ) - du(2,3) + r_cc(i) * d2u(3,1,3)  )  &
                  )
            end do
         end do
       end do

    case (2)
       do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u_tmp_1(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(1))
               du(1,2) = dudr(u_tmp_1, dr(2))
               u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(3))
               du(3,3) = dudr(u_tmp_1, dr(3))
               u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(1))
               call mixed_derivatives_2nd(d2u(1,1,2), u_tmp_2, dr(1), dr(2), i, j, nc)
               u_tmp_2(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(3))
               call mixed_derivatives_2nd(d2u(3,2,3), u_tmp_2, dr(2), dr(3), j, k, nc)
   
               u_tmp_1(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(2))
               face_coeff(0:1, 1) = sin_theta_face(j:j+1)
               call derivatives_2nd(d2u(2,2,2), u_tmp_1, idr2(2), face_coeff_in=face_coeff(0:1, 1))

               ! extra Laplacian term for vectors
               cc(i, j, k, i_out) = cc(i, j, k, i_out) & 
                 + ( 2.0_dp * sin_theta_cc(j)**2 * du(1,2)  &
                   - cc(i, j, k, mg_vec_iphi(2)) &
                   - 2.0_dp * cos_theta_cc(j) * du(3,3) &
                    ) 
   
               ! div_X term
               cc(i, j, k, i_out) = cc(i, j, k, i_out) &
                   + 1.0_dp/3.0_dp * ( &
                      sin_theta_cc(j)**2 * ( r_cc(i) * d2u(1,1,2) + 2.0_dp * du(1,2) ) &
                    + sin_theta_cc(j) * ( d2u(2,2,2) + d2u(3,2,3) ) & 
                    - cc(i, j, k, mg_vec_iphi(2)) - cos_theta_cc(j) * du(3,3) &
                      ) 
            end do
         end do
       end do

    case (3)
       do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(1))
               du(1,3) = dudr(u_tmp_1, dr(3))
               u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(2))
               du(2,3) = dudr(u_tmp_1, dr(3))
               u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(1))
               call mixed_derivatives_2nd(d2u(1,1,3), u_tmp_2, dr(1), dr(3), i, k, nc)
               u_tmp_2(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(2))
               call mixed_derivatives_2nd(d2u(2,2,3), u_tmp_2, dr(2), dr(3), j, k, nc)
               u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(3))
               call derivatives_2nd(d2u(3,3,3), u_tmp_1, idr2(3))

               ! extra Laplacian term for vectors
               cc(i, j, k, i_out) = cc(i, j, k, i_out) &
                   - cc(i, j, k, mg_vec_iphi(3)) &
                   + 2.0_dp * ( sin_theta_cc(j) * du(1,3) + cos_theta_cc(j) * du(2,3) )
   
               ! div_X term
               cc(i, j, k, i_out) = cc(i, j, k, i_out) + &
                   1.0_dp/3.0_dp * ( d2u(3,3,3) + cos_theta_cc(j) * du(2,3) &
                    + sin_theta_cc(j) * ( r_cc(i) * d2u(1,1,3) + 2.0_dp * du(1,3) + d2u(2,2,3) ) &
                   ) 
            end do
         end do
       end do
    case default
      cc(DTIMES(1:nc), i_out) = 0.0_dp
    end select
#endif
    end associate
  end subroutine box_slbeta

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> spherical geometry.
  subroutine box_gs_slbeta(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di, n

    real(dp)                  :: face_coeff(0:1, NDIM) = 1.0_dp
    real(dp)                  :: pre_fac(NDIM) = 1.0_dp
    real(dp)                  :: u(-1:1, NDIM, 3)
    real(dp)                  :: u_tmp_1(-1:1)
    real(dp)                  :: du(3, NDIM)
    real(dp)                  :: d2u(3, NDIM, NDIM)
    real(dp)                  :: dLdu_d2u(3, NDIM, NDIM)

    logical                   :: redblack
    real(dp)                  :: idr2(NDIM), dr(NDIM)
    real(dp)                  :: r_face(1:nc+1), r_cc(0:nc+1)
    real(dp)                  :: Lop(3), dLdu(3)
#if NDIM != 1
    integer                   :: i_dir
    real(dp)                  :: u_tmp_2(-1:1,-1:1)
    real(dp)                  :: sin_theta_cc(0:nc+1)
    real(dp)                  :: sin_theta_face(1:nc+1)
    real(dp)                  :: cos_theta_cc(0:nc+1)
#endif

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i-1.0_dp, i=1,nc+1)]
    r_cc  = (mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=0,nc+1)])
#if NDIM != 1
    sin_theta_cc    = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-0.5_dp, i=0,nc+1)])
    sin_theta_face  = dsin(mg%boxes(id)%r_min(2) + dr(2) * [(i-1.0_dp, i=1,nc+1)])
    cos_theta_cc    = dcos(mg%boxes(id)%r_min(2) + dr(2) * [(i-0.5_dp, i=0,nc+1)]) 
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
    associate (cc => mg%boxes(id)%cc)
#if NDIM == 1
      n = mg_vec_iphi(1)
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)
      do i = i0, nc, di
         ! Laplacian
         u(-1:1, 1, 1) = cc(i-1:i+1, n)
         face_coeff(0:1, 1) = r_face(i:i+1)**2
         call Laplacian(Lop(1), u(:,:,1), idr2, dLdu = dLdu(1), face_coeff_in = face_coeff)

         ! source terms
         Lop(1)  = 4.0_dp/3.0_dp * ( Lop(1) - 2.0_dp * cc(i, n) )
         dLdu(1) = 4.0_dp/3.0_dp * ( dLdu(1) - 2.0_dp )
 
         cc(i, n) = cc(i, n) - ( Lop(1) - cc(i, mg_vec_irhs(1)) ) / dLdu(1)
      end do
      cc(:, mg_vec_iphi(2)) = 0.0_dp
      cc(:, mg_vec_iphi(3)) = 0.0_dp

#elif NDIM == 2
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)
         do i = i0, nc, di

            do i_dir = 1, 3
               n = mg_vec_iphi(i_dir)
               ! Laplacian term
               u(-1:1, 1, i_dir) = cc(i-1:i+1, j, n)
               u(-1:1, 2, i_dir) = cc(i, j-1:j+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)**2
               face_coeff(0:1, 2) = sin_theta_face(j:j+1)
               ! note: here is different from the scalar cases
               pre_fac(1) = sin_theta_cc(j)**2
               pre_fac(2) = sin_theta_cc(j)
               call Laplacian(Lop(i_dir), u(:,:,i_dir), idr2, dLdu = dLdu(i_dir), face_coeff_in = face_coeff, pre_fac_in = pre_fac)

               select case (i_dir)
               case (1)
                 u_tmp_1(-1:1) = cc(i-1:i+1, j, mg_vec_iphi(2))
                 du(2,1) = dudr(u_tmp_1, dr(1))
                 u_tmp_1(-1:1) = cc(i, j-1:j+1, mg_vec_iphi(2))
                 du(2,2) = dudr(u_tmp_1, dr(2))
                 u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, mg_vec_iphi(2))
                 call mixed_derivatives_2nd(d2u(2,1,2), u_tmp_2, dr(1), dr(2), i, j, nc, dLdu=dLdu_d2u(2,1,2))
                 u_tmp_1(-1:1) = cc(i-1:i+1, j, mg_vec_iphi(1))
                 face_coeff(0:1, 1) = r_face(i:i+1)**2
                 call derivatives_2nd(d2u(1,1,1), u_tmp_1, idr2(1), &
                       dLdu=dLdu_d2u(1,1,1), face_coeff_in=face_coeff(0:1, 1))

                 ! extra Laplacian term for vectors
                 u_tmp_1(-1:1) = sin_theta_cc(j-1:j+1) * cc(i, j-1:j+1, mg_vec_iphi(2))
                 Lop(i_dir) = Lop(i_dir) - 2.0_dp * ( &
                     + sin_theta_cc(j)**2 * cc(i, j, mg_vec_iphi(1)) &
                     + sin_theta_cc(j) * dudr(u_tmp_1, dr(2)) ) 
                 dLdu(i_dir) = dLdu(i_dir) - 2.0_dp * sin_theta_cc(j)**2

                 ! div_X term
                 Lop(i_dir) = Lop(i_dir) + &
                     1.0_dp/3.0_dp * ( sin_theta_cc(j)**2 * ( d2u(1,1,1) &
                    - 2.0_dp * cc(i, j, mg_vec_iphi(1)) - du(2,2) )   &
                    - sin_theta_cc(j) * cos_theta_cc(j) * cc(i, j, mg_vec_iphi(2))  &
                    + r_cc(i) * ( sin_theta_cc(j)**2 * d2u(2,1,2) &
                               + sin_theta_cc(j) * cos_theta_cc(j) * du(2,1) ) &
                    )
                 ! div_X for dLdu
                 dLdu(i_dir) = dLdu(i_dir) + &
                  1.0_dp/3.0_dp * ( sin_theta_cc(j)**2 * ( dLdu_d2u(1,1,1) - 2.0_dp &
                                      + r_cc(i) * dLdu_d2u(2,1,2) ) ) 
               case (2)
                 u_tmp_1(-1:1) = cc(i, j-1:j+1, mg_vec_iphi(1))
                 du(1,2) = dudr(u_tmp_1, dr(2))
                 u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, mg_vec_iphi(1))
                 call mixed_derivatives_2nd(d2u(1,1,2), u_tmp_2, dr(1), dr(2), i, j, nc, dLdu=dLdu_d2u(1,1,2))
                 u_tmp_1(-1:1) = cc(i, j-1:j+1, mg_vec_iphi(2))
                 face_coeff(0:1, 1) = sin_theta_face(j:j+1)
                 call derivatives_2nd(d2u(2,2,2), u_tmp_1, idr2(2), &
                      dLdu=dLdu_d2u(2,2,2), face_coeff_in=face_coeff(0:1, 1))

                 ! extra Laplacian term for vectors
                 Lop(i_dir) = Lop(i_dir) + ( &
                     + 2.0_dp * sin_theta_cc(j)**2 * du(1,2)  &
                     - cc(i, j, mg_vec_iphi(2)) &
                      ) 
                 dLdu(i_dir) = dLdu(i_dir) - 1.0_dp

                 ! div_X term
                 Lop(i_dir) = Lop(i_dir) &
                     + 1.0_dp/3.0_dp * ( &
                        sin_theta_cc(j)**2 * ( r_cc(i) * d2u(1,1,2) + 2.0_dp * du(1,2) ) &
                      + sin_theta_cc(j) * d2u(2,2,2) - cc(i, j, mg_vec_iphi(2))  &
                        ) 

                 ! div_X for dLdu
                 dLdu(i_dir) = dLdu(i_dir)&
                     + 1.0_dp/3.0_dp * ( &
                        sin_theta_cc(j)**2 * ( r_cc(i) * dLdu_d2u(1,1,2) ) &
                      + sin_theta_cc(j) * dLdu_d2u(2,2,2) - 1.0_dp  &
                        ) 
               case (3)
                 ! extra Laplacian term for vectors
                 Lop(i_dir) = Lop(i_dir) - cc(i, j, mg_vec_iphi(3)) 
                 dLdu(i_dir) = dLdu(i_dir) - 1.0_dp
               case default
                 Lop(i_dir) = 0.0_dp
                 dLdu(i_dir) = 1.0_dp
               end select
            end do ! end i_dir

            do i_dir = 1, 3
               cc(i,j, mg_vec_iphi(i_dir)) = cc(i,j, mg_vec_iphi(i_dir)) &
                                       - (Lop(i_dir) - cc(i,j, mg_vec_irhs(i_dir))) / dLdu(i_dir)
            end do ! end i_dir

         end do ! end i
      end do ! end j

#elif NDIM == 3
    do k = 1, nc
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)
         do i = i0, nc, di

            do i_dir = 1, 3
               n = mg_vec_iphi(i_dir)
               ! Laplacian term * r**2 sin**2(theta)
               u(-1:1, 1, i_dir) = cc(i-1:i+1, j, k, n)
               u(-1:1, 2, i_dir) = cc(i, j-1:j+1, k, n)
               u(-1:1, 3, i_dir) = cc(i, j, k-1:k+1, n)
               face_coeff(0:1, 1) = r_face(i:i+1)**2
               face_coeff(0:1, 2) = sin_theta_face(j:j+1)
               pre_fac(1) = sin_theta_cc(j)**2
               pre_fac(2) = sin_theta_cc(j)
               call Laplacian( Lop(i_dir), u(:,:,i_dir), idr2, dLdu = dLdu(i_dir) , face_coeff_in = face_coeff, pre_fac_in = pre_fac)

               select case (i_dir)
               case (1)
                  u_tmp_1(-1:1) = cc(i-1:i+1, j, k, mg_vec_iphi(2))
                  du(2,1) = dudr(u_tmp_1, dr(1))
                  u_tmp_1(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(2))
                  du(2,2) = dudr(u_tmp_1, dr(2))
                  u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(2))
                  du(2,3) = dudr(u_tmp_1, dr(3))
                  u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(3))
                  du(3,3) = dudr(u_tmp_1, dr(3))
                  u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(2))
                  call mixed_derivatives_2nd(d2u(2,1,2), u_tmp_2, dr(1), dr(2), i, j, nc, dLdu=dLdu_d2u(2,1,2))
                  u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(3))
                  call mixed_derivatives_2nd(d2u(3,1,3), u_tmp_2, dr(1), dr(3), i, k, nc, dLdu=dLdu_d2u(3,1,3))
   
                  u_tmp_1(-1:1) = cc(i-1:i+1, j, k, mg_vec_iphi(1))
                  face_coeff(0:1, 1) = r_face(i:i+1)**2
                  call derivatives_2nd(d2u(1,1,1), u_tmp_1, idr2(1), face_coeff_in=face_coeff(0:1, 1), dLdu=dLdu_d2u(1,1,1))
      
                  u_tmp_1(-1:1) = sin_theta_cc(j-1:j+1) * cc(i, j-1:j+1, k, mg_vec_iphi(2))

                  ! extra Laplacian term for vectors
                  Lop(i_dir) = Lop(i_dir) &
                    - 2.0_dp * ( &
                   + sin_theta_cc(j)**2 * cc(i, j, k, mg_vec_iphi(1)) &
                   + sin_theta_cc(j) * ( dudr(u_tmp_1, dr(2)) + du(3,3) ) &
                              ) 
                  dLdu(i_dir) = dLdu(i_dir) &
                    - 2.0_dp * ( sin_theta_cc(j)**2 ) 
      
                  ! div_X term
                  Lop(i_dir) = Lop(i_dir) + &
                   1.0_dp/3.0_dp * ( &
                     sin_theta_cc(j)**2 * ( d2u(1,1,1) - 2.0_dp * cc(i, j, k, mg_vec_iphi(1)) - du(2,2) + r_cc(i) * d2u(2,1,2) )   &
                       + sin_theta_cc(j) * ( cos_theta_cc(j) * ( r_cc(i) * du(2,1) - cc(i, j, k, mg_vec_iphi(2)) ) - du(2,3) + r_cc(i) * d2u(3,1,3)  ) & 
                  )
                  ! div_X for dLdu
                  dLdu(i_dir) = dLdu(i_dir) + &
                   1.0_dp/3.0_dp * ( &
                     sin_theta_cc(j)**2 * ( dLdu_d2u(1,1,1) - 2.0_dp + r_cc(i) * dLdu_d2u(2,1,2) )   &
                       + sin_theta_cc(j) * ( r_cc(i) * dLdu_d2u(3,1,3) ) & 
                  )

               case (2)
                  u_tmp_1(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(1))
                  du(1,2) = dudr(u_tmp_1, dr(2))
                  u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(3))
                  du(3,3) = dudr(u_tmp_1, dr(3))
                  u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j-1:j+1, k, mg_vec_iphi(1))
                  call mixed_derivatives_2nd(d2u(1,1,2), u_tmp_2, dr(1), dr(2), i, j, nc, dLdu=dLdu_d2u(1,1,2))
                  u_tmp_2(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(3))
                  call mixed_derivatives_2nd(d2u(3,2,3), u_tmp_2, dr(2), dr(3), j, k, nc, dLdu=dLdu_d2u(3,2,3))
      
                  u_tmp_1(-1:1) = cc(i, j-1:j+1, k, mg_vec_iphi(2))
                  face_coeff(0:1, 1) = sin_theta_face(j:j+1)
                  call derivatives_2nd(d2u(2,2,2), u_tmp_1, idr2(2), face_coeff_in=face_coeff(0:1, 1), dLdu=dLdu_d2u(2,2,2))


                  ! extra Laplacian term for vectors
                  Lop(i_dir) = Lop(i_dir) &
                   + ( 2.0_dp * sin_theta_cc(j)**2 * du(1,2)  &
                   - cc(i, j, k, mg_vec_iphi(2)) &
                   - 2.0_dp * cos_theta_cc(j) * du(3,3) &
                    ) 
                  dLdu(i_dir) = dLdu(i_dir) - 1.0_dp 

                  ! div_X term
                  Lop(i_dir) = Lop(i_dir) &
                    + 1.0_dp/3.0_dp * ( &
                       sin_theta_cc(j)**2 * ( r_cc(i) * d2u(1,1,2) + 2.0_dp * du(1,2) ) &
                     + sin_theta_cc(j) * ( d2u(2,2,2) + d2u(3,2,3) ) & 
                     - cc(i, j, k, mg_vec_iphi(2)) - cos_theta_cc(j) * du(3,3) &
                       ) 
 
                  ! div_X for dLdu
                  dLdu(i_dir) = dLdu(i_dir)&
                    + 1.0_dp/3.0_dp * ( &
                       sin_theta_cc(j)**2 * ( r_cc(i) * dLdu_d2u(1,1,2) ) &
                     + sin_theta_cc(j) * ( dLdu_d2u(2,2,2) + dLdu_d2u(3,2,3) ) & 
                     - 1.0_dp  &
                       ) 

               case (3)
                  u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(1))
                  du(1,3) = dudr(u_tmp_1, dr(3))
                  u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(2))
                  du(2,3) = dudr(u_tmp_1, dr(3))
                  u_tmp_2(-1:1,-1:1) = cc(i-1:i+1, j, k-1:k+1, mg_vec_iphi(1))
                  call mixed_derivatives_2nd(d2u(1,1,3), u_tmp_2, dr(1), dr(3), i, k, nc, dLdu=dLdu_d2u(1,1,3))
                  u_tmp_2(-1:1,-1:1) = cc(i, j-1:j+1, k-1:k+1, mg_vec_iphi(2))
                  call mixed_derivatives_2nd(d2u(2,2,3), u_tmp_2, dr(2), dr(3), j, k, nc, dLdu=dLdu_d2u(2,2,3))
                  u_tmp_1(-1:1) = cc(i, j, k-1:k+1, mg_vec_iphi(3))
                  call derivatives_2nd(d2u(3,3,3), u_tmp_1, idr2(3), dLdu=dLdu_d2u(3,3,2))

                  ! extra Laplacian term for vectors
                  Lop(i_dir) = Lop(i_dir) &
                   - cc(i, j, k, mg_vec_iphi(3)) &
                   + 2.0_dp * ( sin_theta_cc(j) * du(1,3) + cos_theta_cc(j) * du(2,3) )
                  dLdu(i_dir) = dLdu(i_dir) - 1.0_dp 
 
                  ! div_X term
                  Lop(i_dir) = Lop(i_dir) + &
                   1.0_dp/3.0_dp * ( d2u(3,3,3) + cos_theta_cc(j) * du(2,3) &
                    + sin_theta_cc(j) * ( r_cc(i) * d2u(1,1,3) + 2.0_dp * du(1,3) + d2u(2,2,3) ) &
                   ) 
                  ! div_X for dLdu
                  dLdu(i_dir) = dLdu(i_dir) &
                   + 1.0_dp/3.0_dp * ( dLdu_d2u(3,3,3) &
                    + sin_theta_cc(j) * ( r_cc(i) * dLdu_d2u(1,1,3) + dLdu_d2u(2,2,3) ) &
                   ) 
               case default
                 Lop(i_dir) = 0.0_dp
                 dLdu(i_dir) = 1.0_dp
               end select
            end do ! end i_dir

            do i_dir = 1, 3
               cc(i, j, k, mg_vec_iphi(i_dir)) = cc(i, j, k, mg_vec_iphi(i_dir)) &
                                       - (Lop(i_dir) - cc(i, j, k, mg_vec_irhs(i_dir))) / dLdu(i_dir)
            end do ! end i_dir

         end do ! end i
      end do ! end j
    end do ! end k
#endif
    end associate
  end subroutine box_gs_slbeta
end module m_cfc_beta
