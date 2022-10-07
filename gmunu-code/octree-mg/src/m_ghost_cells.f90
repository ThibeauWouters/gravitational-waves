#include "cpp_macros.h"
module m_ghost_cells
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: mg_ghost_cell_buffer_size
  public :: mg_fill_ghost_cells
  public :: mg_fill_ghost_cells_lvl
  public :: mg_phi_bc_store

contains

  !> Specify minimum buffer size (per process) for communication
  subroutine mg_ghost_cell_buffer_size(mg, n_send, n_recv, dsize)
    type(mg_t), intent(inout) :: mg
    integer, intent(out)      :: n_send(0:mg%n_cpu-1)
    integer, intent(out)      :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)      :: dsize
    integer                   :: i, id, lvl, nc

    allocate(mg%comm_ghostcell%n_send(0:mg%n_cpu-1, &
         mg%first_normal_lvl:mg%highest_lvl))
    allocate(mg%comm_ghostcell%n_recv(0:mg%n_cpu-1, &
         mg%first_normal_lvl:mg%highest_lvl))

    dsize = mg%box_size**(NDIM-1)

    do lvl = mg%first_normal_lvl, mg%highest_lvl
       nc               = mg%box_size_lvl(lvl)
       mg%buf(:)%i_send = 0
       mg%buf(:)%i_recv = 0
       mg%buf(:)%i_ix   = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call buffer_ghost_cells(mg, id, nc, 1, dry_run=.true.)
       end do

       if (lvl > 1) then
          do i = 1, size(mg%lvls(lvl-1)%my_ref_bnds)
             id = mg%lvls(lvl-1)%my_ref_bnds(i)
             call buffer_refinement_boundaries(mg, id, nc, 1, dry_run=.true.)
          end do
       end if

       ! Set ghost cells to received data
       mg%buf(:)%i_recv = 0
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call set_ghost_cells(mg, id, nc, 1, dry_run=.true.)
       end do

       mg%comm_ghostcell%n_send(:, lvl) = mg%buf(:)%i_send/dsize
       mg%comm_ghostcell%n_recv(:, lvl) = mg%buf(:)%i_recv/dsize
    end do

    n_send = maxval(mg%comm_ghostcell%n_send, dim=2)
    n_recv = maxval(mg%comm_ghostcell%n_recv, dim=2)
  end subroutine mg_ghost_cell_buffer_size

  !> Store boundary conditions for the solution variable, this can speed up
  !> calculations if the same boundary conditions are re-used.
  subroutine mg_phi_bc_store(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: lvl, nc

    nc = mg%box_size

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       call mg_phi_bc_store_lvl(mg, lvl, nc)
    end do

    mg%phi_bc_data_stored = .true.
  end subroutine mg_phi_bc_store

  subroutine mg_phi_bc_store_lvl(mg, lvl, nc)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer, intent(in)       :: nc
#if NDIM == 1
    real(dp)                  :: bc(1)
#elif NDIM == 2
    real(dp)                  :: bc(nc)
#elif NDIM == 3
    real(dp)                  :: bc(nc, nc)
#endif
    integer                   :: i, id, nb, nb_id, bc_type

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       do nb = 1, mg_num_neighbors
          nb_id = mg%boxes(id)%neighbors(nb)
          if (nb_id < mg_no_box) then
             ! Physical boundary
             if (associated(mg%bc(nb, mg_iphi)%boundary_cond)) then
                call mg%bc(nb, mg_iphi)%boundary_cond(mg%boxes(id), nc, &
                     mg_iphi, nb, bc_type, bc)
             else
                bc_type = mg%bc(nb, mg_iphi)%bc_type
                bc      = mg%bc(nb, mg_iphi)%bc_value
             end if

             ! Store the boundary condition type. This is not globally set in
             ! the tree, but all negative values are treated the same in
             ! other parts of the code
             mg%boxes(id)%neighbors(nb) = bc_type

             ! Store ghost cell data in the right-hand side
             call box_set_gc(mg%boxes(id), nb, nc, mg_irhs, bc)
          end if
       end do
    end do
  end subroutine mg_phi_bc_store_lvl

  !> Fill ghost cells at all grid levels
  subroutine mg_fill_ghost_cells(mg, iv)
    type(mg_t)          :: mg
    integer, intent(in) :: iv !< Index of variable
    integer             :: lvl

    do lvl = mg%lowest_lvl, mg%highest_lvl
       call mg_fill_ghost_cells_lvl(mg, lvl, iv)
    end do
  end subroutine mg_fill_ghost_cells

  !> Fill ghost cells at a grid level
  subroutine mg_fill_ghost_cells_lvl(mg, lvl, iv)
    use m_communication
    type(mg_t)                   :: mg
    integer, intent(in)          :: lvl
    integer, intent(in)          :: iv !< Index of variable
    integer                      :: i, id, dsize, nc

    if (lvl < mg%lowest_lvl) &
         error stop "fill_ghost_cells_lvl: lvl < lowest_lvl"
    if (lvl > mg%highest_lvl) &
         error stop "fill_ghost_cells_lvl: lvl > highest_lvl"

    nc               = mg%box_size_lvl(lvl)

    if (lvl >= mg%first_normal_lvl) then
       dsize            = nc**(NDIM-1)
       mg%buf(:)%i_send = 0
       mg%buf(:)%i_recv = 0
       mg%buf(:)%i_ix   = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call buffer_ghost_cells(mg, id, nc, iv, .false.)
       end do

       if (lvl > 1) then
          do i = 1, size(mg%lvls(lvl-1)%my_ref_bnds)
             id = mg%lvls(lvl-1)%my_ref_bnds(i)
             call buffer_refinement_boundaries(mg, id, nc, iv, .false.)
          end do
       end if

       ! Transfer data between processes
       mg%buf(:)%i_recv = mg%comm_ghostcell%n_recv(:, lvl) * dsize
       call sort_and_transfer_buffers(mg, dsize)

       ! Set ghost cells to received data
       mg%buf(:)%i_recv = 0
    end if

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call set_ghost_cells(mg, id, nc, iv, .false.)
    end do
  end subroutine mg_fill_ghost_cells_lvl

  subroutine buffer_ghost_cells(mg, id, nc, iv, dry_run)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    logical, intent(in)       :: dry_run
    integer                   :: nb, nb_id, nb_rank

    do nb = 1, mg_num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)

       if (nb_id > mg_no_box) then
          ! There is a neighbor
          nb_rank    = mg%boxes(nb_id)%rank

          if (nb_rank /= mg%my_rank) then
             call buffer_for_nb(mg, mg%boxes(id), nc, iv, nb_id, nb_rank, &
                  nb, dry_run)
          end if
       end if
    end do
  end subroutine buffer_ghost_cells

  subroutine buffer_refinement_boundaries(mg, id, nc, iv, dry_run)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    logical, intent(in)       :: dry_run
    integer                   :: nb, nb_id, c_ids(2**(NDIM-1))
    integer                   :: n, c_id, c_rank

    do nb = 1, mg_num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)
       if (nb_id > mg_no_box) then
          if (mg_has_children(mg%boxes(nb_id))) then
             c_ids = mg%boxes(nb_id)%children(&
                  mg_child_adj_nb(:, mg_neighb_rev(nb)))

             do n = 1, mg_num_children/2
                c_id = c_ids(n)
                c_rank = mg%boxes(c_id)%rank

                if (c_rank /= mg%my_rank) then
                   ! Send all coarse ghost cells
                   call buffer_for_fine_nb(mg, mg%boxes(id), nc, iv, c_id, &
                        c_rank, nb, dry_run)
                end if
             end do
          end if
       end if
    end do
  end subroutine buffer_refinement_boundaries

  !> The routine that actually fills the ghost cells
  subroutine set_ghost_cells(mg, id, nc, iv, dry_run)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    logical, intent(in)       :: dry_run
#if NDIM == 1
    real(dp)                  :: bc(1)
#elif NDIM == 2
    real(dp)                  :: bc(nc)
#elif NDIM == 3
    real(dp)                  :: bc(nc, nc)
#endif
    integer                   :: nb, nb_id, nb_rank, bc_type

    do nb = 1, mg_num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)

       if (nb_id > mg_no_box) then
          ! There is a neighbor
          nb_rank = mg%boxes(nb_id)%rank

          if (nb_rank /= mg%my_rank) then
             call fill_buffered_nb(mg, mg%boxes(id), nb_rank, &
                  nb, nc, iv, dry_run)
          else if (.not. dry_run) then
             call copy_from_nb(mg%boxes(id), mg%boxes(nb_id), &
                  nb, nc, iv)
          end if
       else if (nb_id == mg_no_box) then
          ! Refinement boundary
          call fill_refinement_bnd(mg, id, nb, nc, iv, dry_run)
       else if (.not. dry_run) then
          ! Physical boundary
          if (mg%phi_bc_data_stored .and. iv == mg_iphi) then
             ! Copy the boundary conditions stored in the ghost cells of the
             ! right-hand side
             call box_get_gc(mg%boxes(id), nb, nc, mg_irhs, bc)
             bc_type = nb_id
          else
             if (associated(mg%bc(nb, iv)%boundary_cond)) then
                call mg%bc(nb, iv)%boundary_cond(mg%boxes(id), nc, iv, &
                     nb, bc_type, bc)
             else
                bc_type = mg%bc(nb, iv)%bc_type
                bc = mg%bc(nb, iv)%bc_value
             end if
          end if

          call box_set_gc(mg%boxes(id), nb, nc, iv, bc)
          call bc_to_gc(mg, id, nc, iv, nb, bc_type)
       end if
    end do
  end subroutine set_ghost_cells

  subroutine fill_refinement_bnd(mg, id, nb, nc, iv, dry_run)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer, intent(in)       :: nb
    logical, intent(in)       :: dry_run
#if NDIM == 1
    real(dp)                  :: gc(1)
#elif NDIM == 2
    real(dp)                  :: gc(nc)
#elif NDIM == 3
    real(dp)                  :: gc(nc, nc)
#endif
    integer                   :: p_id, p_nb_id, ix_offset(NDIM)
    integer                   :: i, dsize, p_nb_rank

    dsize     = nc**(NDIM-1)
    p_id      = mg%boxes(id)%parent
    p_nb_id   = mg%boxes(p_id)%neighbors(nb)
    p_nb_rank = mg%boxes(p_nb_id)%rank

    if (p_nb_rank /= mg%my_rank) then
       i = mg%buf(p_nb_rank)%i_recv
       if (.not. dry_run) then
          gc = reshape(mg%buf(p_nb_rank)%recv(i+1:i+dsize), shape(gc))
       end if
       mg%buf(p_nb_rank)%i_recv = mg%buf(p_nb_rank)%i_recv + dsize
    else if (.not. dry_run) then
       ix_offset = mg_get_child_offset(mg, id)
       call box_gc_for_fine_neighbor(mg%boxes(p_nb_id), mg_neighb_rev(nb), &
            ix_offset, nc, iv, gc)
    end if

    if (.not. dry_run) then
       if (associated(mg%bc(nb, iv)%refinement_bnd)) then
          call mg%bc(nb, iv)%refinement_bnd(mg%boxes(id), nc, iv, nb, gc)
       else
          call sides_rb(mg%boxes(id), nc, iv, nb, gc)
       end if
    end if
  end subroutine fill_refinement_bnd

  subroutine copy_from_nb(box, box_nb, nb, nc, iv)
    type(mg_box_t), intent(inout) :: box
    type(mg_box_t), intent(in)    :: box_nb
    integer, intent(in)           :: nb
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv
#if NDIM == 1
    real(dp)                      :: gc(1)
#elif NDIM == 2
    real(dp)                      :: gc(nc)
#elif NDIM == 3
    real(dp)                      :: gc(nc, nc)
#endif

    call box_gc_for_neighbor(box_nb, mg_neighb_rev(nb), nc, iv, gc)
    call box_set_gc(box, nb, nc, iv, gc)
  end subroutine copy_from_nb

  subroutine buffer_for_nb(mg, box, nc, iv, nb_id, nb_rank, nb, dry_run)
    type(mg_t), intent(inout)  :: mg
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)        :: nc
    integer, intent(in)        :: iv
    integer, intent(in)        :: nb_id
    integer, intent(in)        :: nb_rank
    integer, intent(in)        :: nb
    logical, intent(in)        :: dry_run
    integer                    :: i, dsize
#if NDIM == 1
    real(dp)                   :: gc(1)
#elif NDIM == 2
    real(dp)                   :: gc(nc)
#elif NDIM == 3
    real(dp)                   :: gc(nc, nc)
#endif

    i     = mg%buf(nb_rank)%i_send
    dsize = nc**(NDIM-1)

    if (.not. dry_run) then
       call box_gc_for_neighbor(box, nb, nc, iv, gc)
       mg%buf(nb_rank)%send(i+1:i+dsize) = pack(gc, .true.)
    end if

    ! Later the buffer is sorted, using the fact that loops go from low to high
    ! box id, and we fill ghost cells according to the neighbor order
    i = mg%buf(nb_rank)%i_ix
    if (.not. dry_run) then
       mg%buf(nb_rank)%ix(i+1) = mg_num_neighbors * nb_id + mg_neighb_rev(nb)
    end if

    mg%buf(nb_rank)%i_send = mg%buf(nb_rank)%i_send + dsize
    mg%buf(nb_rank)%i_ix   = mg%buf(nb_rank)%i_ix + 1
  end subroutine buffer_for_nb

  subroutine buffer_for_fine_nb(mg, box, nc, iv, fine_id, fine_rank, nb, dry_run)
    type(mg_t), intent(inout)  :: mg
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)        :: nc
    integer, intent(in)        :: iv
    integer, intent(in)        :: fine_id
    integer, intent(in)        :: fine_rank
    integer, intent(in)        :: nb
    logical, intent(in)        :: dry_run
    integer                    :: i, dsize, ix_offset(NDIM)
#if NDIM == 1
    real(dp)                   :: gc(1)
#elif NDIM == 2
    real(dp)                   :: gc(nc)
#elif NDIM == 3
    real(dp)                   :: gc(nc, nc)
#endif

    i     = mg%buf(fine_rank)%i_send
    dsize = nc**(NDIM-1)

    if (.not. dry_run) then
       ix_offset = mg_get_child_offset(mg, fine_id)
       call box_gc_for_fine_neighbor(box, nb, ix_offset, nc, iv, gc)
       mg%buf(fine_rank)%send(i+1:i+dsize) = pack(gc, .true.)
    end if

    ! Later the buffer is sorted, using the fact that loops go from low to high
    ! box id, and we fill ghost cells according to the neighbor order
    i = mg%buf(fine_rank)%i_ix
    if (.not. dry_run) then
       mg%buf(fine_rank)%ix(i+1) = mg_num_neighbors * fine_id + &
            mg_neighb_rev(nb)
    end if

    mg%buf(fine_rank)%i_send = mg%buf(fine_rank)%i_send + dsize
    mg%buf(fine_rank)%i_ix   = mg%buf(fine_rank)%i_ix + 1
  end subroutine buffer_for_fine_nb

  subroutine fill_buffered_nb(mg, box, nb_rank, nb, nc, iv, dry_run)
    type(mg_t), intent(inout)  :: mg
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)        :: nb_rank
    integer, intent(in)        :: nb
    integer, intent(in)        :: nc
    integer, intent(in)        :: iv
    logical, intent(in)        :: dry_run
    integer                    :: i, dsize
#if NDIM == 1
    real(dp)                   :: gc(1)
#elif NDIM == 2
    real(dp)                   :: gc(nc)
#elif NDIM == 3
    real(dp)                   :: gc(nc, nc)
#endif

    i     = mg%buf(nb_rank)%i_recv
    dsize = nc**(NDIM-1)

    if (.not. dry_run) then
#if NDIM > 1
       gc = reshape(mg%buf(nb_rank)%recv(i+1:i+dsize), shape(gc))
#else
       gc = mg%buf(nb_rank)%recv(i+1)
#endif
       call box_set_gc(box, nb, nc, iv, gc)
    end if
    mg%buf(nb_rank)%i_recv = mg%buf(nb_rank)%i_recv + dsize

  end subroutine fill_buffered_nb

  subroutine box_gc_for_neighbor(box, nb, nc, iv, gc)
    type(mg_box_t), intent(in) :: box
    integer, intent(in)     :: nb, nc, iv
#if NDIM == 1
    real(dp), intent(out)   :: gc(1)
#elif NDIM == 2
    real(dp), intent(out)   :: gc(nc)
#elif NDIM == 3
    real(dp), intent(out)   :: gc(nc, nc)
#endif

    select case (nb)
#if NDIM == 1
    case (mg_neighb_lowx)
       gc = box%cc(1, iv)
    case (mg_neighb_highx)
       gc = box%cc(nc, iv)
#elif NDIM == 2
    case (mg_neighb_lowx)
       gc = box%cc(1, 1:nc, iv)
    case (mg_neighb_highx)
       gc = box%cc(nc, 1:nc, iv)
    case (mg_neighb_lowy)
       gc = box%cc(1:nc, 1, iv)
    case (mg_neighb_highy)
       gc = box%cc(1:nc, nc, iv)
#elif NDIM == 3
    case (mg_neighb_lowx)
       gc = box%cc(1, 1:nc, 1:nc, iv)
    case (mg_neighb_highx)
       gc = box%cc(nc, 1:nc, 1:nc, iv)
    case (mg_neighb_lowy)
       gc = box%cc(1:nc, 1, 1:nc, iv)
    case (mg_neighb_highy)
       gc = box%cc(1:nc, nc, 1:nc, iv)
    case (mg_neighb_lowz)
       gc = box%cc(1:nc, 1:nc, 1, iv)
    case (mg_neighb_highz)
       gc = box%cc(1:nc, 1:nc, nc, iv)
#endif
    end select
  end subroutine box_gc_for_neighbor

  !> Get ghost cells for a fine neighbor
  subroutine box_gc_for_fine_neighbor(box, nb, di, nc, iv, gc)
    type(mg_box_t), intent(in) :: box
    integer, intent(in)     :: nb       !< Direction of fine neighbor
    integer, intent(in)     :: di(NDIM) !< Index offset of fine neighbor
    integer, intent(in)     :: nc, iv
#if NDIM == 1
    real(dp), intent(out)   :: gc(1)
    real(dp)                :: tmp
    integer                 :: hnc
#elif NDIM == 2
    real(dp), intent(out)   :: gc(nc)
    real(dp)                :: tmp(0:nc/2+1), grad(NDIM-1)
    integer                 :: i, hnc
#elif NDIM == 3
    real(dp), intent(out)   :: gc(nc, nc)
    real(dp)                :: tmp(0:nc/2+1, 0:nc/2+1), grad(NDIM-1)
    integer                 :: i, j, hnc
#endif

    hnc = nc/2

    ! First fill a temporary array with data next to the fine grid
    select case (nb)
#if NDIM == 1
    case (mg_neighb_lowx)
       tmp = box%cc(1, iv)
    case (mg_neighb_highx)
       tmp = box%cc(nc, iv)
#elif NDIM == 2
    case (mg_neighb_lowx)
       tmp = box%cc(1, di(2):di(2)+hnc+1, iv)
    case (mg_neighb_highx)
       tmp = box%cc(nc, di(2):di(2)+hnc+1, iv)
    case (mg_neighb_lowy)
       tmp = box%cc(di(1):di(1)+hnc+1, 1, iv)
    case (mg_neighb_highy)
       tmp = box%cc(di(1):di(1)+hnc+1, nc, iv)
#elif NDIM == 3
    case (mg_neighb_lowx)
       tmp = box%cc(1, di(2):di(2)+hnc+1, di(3):di(3)+hnc+1, iv)
    case (mg_neighb_highx)
       tmp = box%cc(nc, di(2):di(2)+hnc+1, di(3):di(3)+hnc+1, iv)
    case (mg_neighb_lowy)
       tmp = box%cc(di(1):di(1)+hnc+1, 1, di(3):di(3)+hnc+1, iv)
    case (mg_neighb_highy)
       tmp = box%cc(di(1):di(1)+hnc+1, nc, di(3):di(3)+hnc+1, iv)
    case (mg_neighb_lowz)
       tmp = box%cc(di(1):di(1)+hnc+1, di(2):di(2)+hnc+1, 1, iv)
    case (mg_neighb_highz)
       tmp = box%cc(di(1):di(1)+hnc+1, di(2):di(2)+hnc+1, nc, iv)
#endif
    case default
       error stop
    end select

    ! Now interpolate the coarse grid data to obtain values 'straight' next to
    ! the fine grid points
#if NDIM == 1
    gc = tmp
#elif NDIM == 2
    do i = 1, hnc
       grad(1) = 0.125_dp * (tmp(i+1) - tmp(i-1))
       gc(2*i-1) = tmp(i) - grad(1)
       gc(2*i) = tmp(i) + grad(1)
    end do
#elif NDIM == 3
    do j = 1, hnc
       do i = 1, hnc
          grad(1)          = 0.125_dp * (tmp(i+1, j) - tmp(i-1, j))
          grad(2)          = 0.125_dp * (tmp(i, j+1) - tmp(i, j-1))
          gc(2*i-1, 2*j-1) = tmp(i, j) - grad(1) - grad(2)
          gc(2*i, 2*j-1)   = tmp(i, j) + grad(1) - grad(2)
          gc(2*i-1, 2*j)   = tmp(i, j) - grad(1) + grad(2)
          gc(2*i, 2*j)     = tmp(i, j) + grad(1) + grad(2)
       end do
    end do
#endif
  end subroutine box_gc_for_fine_neighbor

  subroutine box_get_gc(box, nb, nc, iv, gc)
    type(mg_box_t), intent(in) :: box
    integer, intent(in)        :: nb, nc, iv
#if NDIM == 1
    real(dp), intent(out)       :: gc(1)
#elif NDIM == 2
    real(dp), intent(out)       :: gc(nc)
#elif NDIM == 3
    real(dp), intent(out)       :: gc(nc, nc)
#endif

    select case (nb)
#if NDIM == 1
    case (mg_neighb_lowx)
       gc = box%cc(0, iv)
    case (mg_neighb_highx)
       gc = box%cc(nc+1, iv)
#elif NDIM == 2
    case (mg_neighb_lowx)
       gc = box%cc(0, 1:nc, iv)
    case (mg_neighb_highx)
       gc = box%cc(nc+1, 1:nc, iv)
    case (mg_neighb_lowy)
       gc = box%cc(1:nc, 0, iv)
    case (mg_neighb_highy)
       gc = box%cc(1:nc, nc+1, iv)
#elif NDIM == 3
    case (mg_neighb_lowx)
       gc = box%cc(0, 1:nc, 1:nc, iv)
    case (mg_neighb_highx)
       gc = box%cc(nc+1, 1:nc, 1:nc, iv)
    case (mg_neighb_lowy)
       gc = box%cc(1:nc, 0, 1:nc, iv)
    case (mg_neighb_highy)
       gc = box%cc(1:nc, nc+1, 1:nc, iv)
    case (mg_neighb_lowz)
       gc = box%cc(1:nc, 1:nc, 0, iv)
    case (mg_neighb_highz)
       gc = box%cc(1:nc, 1:nc, nc+1, iv)
#endif
    end select
  end subroutine box_get_gc

  subroutine box_set_gc(box, nb, nc, iv, gc)
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)        :: nb, nc, iv
#if NDIM == 1
    real(dp), intent(in)       :: gc(1)
#elif NDIM == 2
    real(dp), intent(in)       :: gc(nc)
#elif NDIM == 3
    real(dp), intent(in)       :: gc(nc, nc)
#endif

    select case (nb)
#if NDIM == 1
    case (mg_neighb_lowx)
       box%cc(0, iv)    = gc(1)
    case (mg_neighb_highx)
       box%cc(nc+1, iv) = gc(1)
#elif NDIM == 2
    case (mg_neighb_lowx)
       box%cc(0, 1:nc, iv)    = gc
    case (mg_neighb_highx)
       box%cc(nc+1, 1:nc, iv) = gc
    case (mg_neighb_lowy)
       box%cc(1:nc, 0, iv)    = gc
    case (mg_neighb_highy)
       box%cc(1:nc, nc+1, iv) = gc
#elif NDIM == 3
    case (mg_neighb_lowx)
       box%cc(0, 1:nc, 1:nc, iv)    = gc
    case (mg_neighb_highx)
       box%cc(nc+1, 1:nc, 1:nc, iv) = gc
    case (mg_neighb_lowy)
       box%cc(1:nc, 0, 1:nc, iv)    = gc
    case (mg_neighb_highy)
       box%cc(1:nc, nc+1, 1:nc, iv) = gc
    case (mg_neighb_lowz)
       box%cc(1:nc, 1:nc, 0, iv)    = gc
    case (mg_neighb_highz)
       box%cc(1:nc, 1:nc, nc+1, iv) = gc
#endif
    end select
  end subroutine box_set_gc

  subroutine bc_to_gc(mg, id, nc, iv, nb, bc_type)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer, intent(in)       :: nb      !< Neighbor direction
    integer, intent(in)       :: bc_type !< Type of b.c.
    real(dp)                  :: c0, c1, c2

    !< variables needed for Neumann and Robin b.c.
    real(dp)                  :: dr, r0, r1, pow
#if NDIM == 1
    real(dp)                  :: gc(1)
#elif NDIM == 2
    real(dp)                  :: gc(nc)
#elif NDIM == 3             
    real(dp)                  :: gc(nc,nc)
#endif
    integer                   :: idir 

    ! The Robin B.C. in cart and cylindrical coordinate is quite different, 
    ! we here sperate the implementation
    if ( bc_type == mg_bc_robin .and. &
         mg%geometry_type /= mg_spherical ) then

       call robin_bc_in_non_spherical(mg, id, nc, iv, nb, bc_type)

       return
    end if

    ! If we call the interior point phi1, phi2 and the ghost point phi0, then a
    ! Dirichlet boundary value b can be imposed as:
    ! phi0 = -phi1 + 2*b
    ! A Neumann b.c. can be imposed as:
    ! phi0 = phi1 -/+ dx * b
    ! A Robin b.c. [ d/dx (x*phi) = b ] can be imposed as:
    ! phi0 = ( phi1 * x1 -/+ dx * b ) / x0
    !      =  phi1 * (x1 / x0) -/+ dx / x0 * b 
    ! A continuous boundary (same slope) as:
    ! phi0 = 2 * phi1 - phi2
    !
    ! Below, we set coefficients to handle these cases
    ! where
    ! c0 is the coefficient for boundary value b,
    ! c1 is the coefficient for phi1,
    ! c2 is the coefficient for phi2.
    ! Note that boundary value b is stored in the ghost cell phi0 originally.

    select case (bc_type)
    case (mg_bc_dirichlet)
       c0 = 2
       c1 = -1
       c2 = 0
    case (mg_bc_neumann)
       dr = mg%dr(mg_neighb_dim(nb), mg%boxes(id)%lvl)
       c0 = dr * mg_neighb_high_pm(nb) ! This gives a - or + sign
       c1 = 1
       c2 = 0
    case (mg_bc_continuous)
       c0 = 0
       c1 = 2
       c2 = -1
    case (mg_bc_robin)
       call box_get_gc(mg%boxes(id), nb, nc, iv, gc)
#if NDIM <= 2
       pow = gc(1)
#elif NDIM == 3
       pow = gc(1,1)
#endif
       idir = mg_neighb_dim(nb)
       dr = mg%dr(idir, mg%boxes(id)%lvl) 
       if ( mod(nb,2) == 1 ) then ! it is lowx/y/z
         r0 = mg%boxes(id)%r_min(idir) + dr * (-0.5_dp)
         r1 = r0 + dr
       else ! it is highx/y/z
         r1 = mg%boxes(id)%r_min(idir) + dr * (nc-0.5_dp)
         r0 = r1 + dr
       end if
       c1 = ( r1 / r0 )**pow
       c2 = 0
       !--------------------------------------------------------------------
       ! special Robin B.C. treatment only for alp*psi equation
       ! such that the Robin B.C. is applied on (alp-1) only.
       !--------------------------------------------------------------------
       if ( mg%operator_type==mg_cfc_alp ) then
         ! phi0 =  phi1 * (x1 / x0) + (1-psi1)*(r1/r0-1) 
         ! so we need to modify the ghost cell phi0 = (1-psi1)*(r1/r0-1) 
         c0 = 1.0_dp
         call robin_bc_for_alppsi(mg, id, nc, iv, nb, c1)
       else
         !c0 = dr / r0 * mg_neighb_high_pm(nb) ! This gives a - or + sign
         c0 = 0.0_dp
       end if
    case default
       error stop "bc_to_gc: unknown boundary condition"
    end select

    select case (nb)
#if NDIM == 1
    case (mg_neighb_lowx)
       mg%boxes(id)%cc(0, iv) = &
            c0 * mg%boxes(id)%cc(0, iv) + &
            c1 * mg%boxes(id)%cc(1, iv) + &
            c2 * mg%boxes(id)%cc(2, iv)
    case (mg_neighb_highx)
       mg%boxes(id)%cc(nc+1, iv) = &
            c0 * mg%boxes(id)%cc(nc+1, iv) + &
            c1 * mg%boxes(id)%cc(nc, iv) + &
            c2 * mg%boxes(id)%cc(nc-1, iv)
#elif NDIM == 2
    case (mg_neighb_lowx)
       mg%boxes(id)%cc(0, 1:nc, iv) = &
            c0 * mg%boxes(id)%cc(0, 1:nc, iv) + &
            c1 * mg%boxes(id)%cc(1, 1:nc, iv) + &
            c2 * mg%boxes(id)%cc(2, 1:nc, iv)
    case (mg_neighb_highx)
       mg%boxes(id)%cc(nc+1, 1:nc, iv) = &
            c0 * mg%boxes(id)%cc(nc+1, 1:nc, iv) + &
            c1 * mg%boxes(id)%cc(nc, 1:nc, iv) + &
            c2 * mg%boxes(id)%cc(nc-1, 1:nc, iv)
    case (mg_neighb_lowy)
       mg%boxes(id)%cc(1:nc, 0, iv) = &
            c0 * mg%boxes(id)%cc(1:nc, 0, iv) + &
            c1 * mg%boxes(id)%cc(1:nc, 1, iv) + &
            c2 * mg%boxes(id)%cc(1:nc, 2, iv)
    case (mg_neighb_highy)
       mg%boxes(id)%cc(1:nc, nc+1, iv) = &
            c0 * mg%boxes(id)%cc(1:nc, nc+1, iv) + &
            c1 * mg%boxes(id)%cc(1:nc, nc, iv) + &
            c2 * mg%boxes(id)%cc(1:nc, nc-1, iv)
#elif NDIM == 3
    case (mg_neighb_lowx)
       mg%boxes(id)%cc(0, 1:nc, 1:nc, iv) = &
            c0 * mg%boxes(id)%cc(0, 1:nc, 1:nc, iv) + &
            c1 * mg%boxes(id)%cc(1, 1:nc, 1:nc, iv) + &
            c2 * mg%boxes(id)%cc(2, 1:nc, 1:nc, iv)
    case (mg_neighb_highx)
       mg%boxes(id)%cc(nc+1, 1:nc, 1:nc, iv) = &
            c0 * mg%boxes(id)%cc(nc+1, 1:nc, 1:nc, iv) + &
            c1 * mg%boxes(id)%cc(nc, 1:nc, 1:nc, iv) + &
            c2 * mg%boxes(id)%cc(nc-1, 1:nc, 1:nc, iv)
    case (mg_neighb_lowy)
       mg%boxes(id)%cc(1:nc, 0, 1:nc, iv) = &
            c0 * mg%boxes(id)%cc(1:nc, 0, 1:nc, iv) + &
            c1 * mg%boxes(id)%cc(1:nc, 1, 1:nc, iv) + &
            c2 * mg%boxes(id)%cc(1:nc, 2, 1:nc, iv)
    case (mg_neighb_highy)
       mg%boxes(id)%cc(1:nc, nc+1, 1:nc, iv) = &
            c0 * mg%boxes(id)%cc(1:nc, nc+1, 1:nc, iv) + &
            c1 * mg%boxes(id)%cc(1:nc, nc, 1:nc, iv) + &
            c2 * mg%boxes(id)%cc(1:nc, nc-1, 1:nc, iv)
    case (mg_neighb_lowz)
       mg%boxes(id)%cc(1:nc, 1:nc, 0, iv) = &
            c0 * mg%boxes(id)%cc(1:nc, 1:nc, 0, iv) + &
            c1 * mg%boxes(id)%cc(1:nc, 1:nc, 1, iv) + &
            c2 * mg%boxes(id)%cc(1:nc, 1:nc, 2, iv)
    case (mg_neighb_highz)
       mg%boxes(id)%cc(1:nc, 1:nc, nc+1, iv) = &
            c0 * mg%boxes(id)%cc(1:nc, 1:nc, nc+1, iv) + &
            c1 * mg%boxes(id)%cc(1:nc, 1:nc, nc, iv) + &
            c2 * mg%boxes(id)%cc(1:nc, 1:nc, nc-1, iv)
#endif
    end select
  end subroutine bc_to_gc

  ! fixme: transform alp*psi-1 to alp-1 seems better
  subroutine robin_bc_for_alppsi(mg, id, nc, iv, nb, c1)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer, intent(in)       :: nb      !< Neighbor direction
    real(dp), intent(in)      :: c1

    ! special Robin B.C. treatment only for alp*psi equation
    ! such that the Robin B.C. is applied on (alp-1) only.
    ! phi0 =  phi1 * (x1 / x0) + (1-psi1)*(1-alp1)*(r1/r0-1) 
    ! so we need to modify the ghost cell phi0 = (1-psi1)*(1-alp1)*(r1/r0-1) 
    ! Note: only highx is implemented
    select case (nb)
#if NDIM == 1
    case (mg_neighb_lowx)
       mg%boxes(id)%cc(0, iv) = ( mg%boxes(id)%cc(1, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1) * &
                           ( 1.0d0 - ( 1.0d0 + mg%boxes(id)%cc(1, iv ) ) / mg%boxes(id)%cc(1, mg_itmp2) )
    case (mg_neighb_highx)
       mg%boxes(id)%cc(nc+1, iv) = ( mg%boxes(id)%cc(nc, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1) * &
                           ( 1.0d0 - ( 1.0d0 + mg%boxes(id)%cc(nc, iv ) ) / mg%boxes(id)%cc(nc, mg_itmp2) )
#elif NDIM == 2
    case (mg_neighb_lowx)
       mg%boxes(id)%cc(0, 1:nc, iv) = ( mg%boxes(id)%cc(1, 1:nc, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1)
    case (mg_neighb_highx)
       mg%boxes(id)%cc(nc+1, 1:nc, iv) = ( mg%boxes(id)%cc(nc, 1:nc, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1) * &
                      ( 1.0d0 - ( 1.0d0 + mg%boxes(id)%cc(nc, 1:nc, iv ) ) / mg%boxes(id)%cc(nc, 1:nc, mg_itmp2) )
    case (mg_neighb_lowy)
       mg%boxes(id)%cc(1:nc, 0, iv) = ( mg%boxes(id)%cc(1:nc, 1, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1)
    case (mg_neighb_highy)
       mg%boxes(id)%cc(1:nc, nc+1, iv) = ( mg%boxes(id)%cc(1:nc, nc, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1)
#elif NDIM == 3
    case (mg_neighb_lowx)
       mg%boxes(id)%cc(0, 1:nc, 1:nc, iv) = ( mg%boxes(id)%cc(1, 1:nc, 1:nc, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1)
    case (mg_neighb_highx)
       mg%boxes(id)%cc(nc+1, 1:nc, 1:nc, iv) = ( mg%boxes(id)%cc(nc, 1:nc, 1:nc, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1) * &
                      ( 1.0d0 - ( 1.0d0 + mg%boxes(id)%cc(nc, 1:nc, 1:nc, iv ) ) / mg%boxes(id)%cc(nc, 1:nc, 1:nc, mg_itmp2) )
    case (mg_neighb_lowy)
       mg%boxes(id)%cc(1:nc, 0, 1:nc, iv) = ( mg%boxes(id)%cc(1:nc, 1, 1:nc, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1)
    case (mg_neighb_highy)
       mg%boxes(id)%cc(1:nc, nc+1, 1:nc, iv) = ( mg%boxes(id)%cc(1:nc, nc, 1:nc, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1)
    case (mg_neighb_lowz)
       mg%boxes(id)%cc(1:nc, 1:nc, 0, iv) = ( mg%boxes(id)%cc(1:nc, 1:nc, 1, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1)
    case (mg_neighb_highz)
       mg%boxes(id)%cc(1:nc, 1:nc, nc+1, iv) = ( mg%boxes(id)%cc(1:nc, 1:nc, nc, mg_itmp2) - 1.0d0 ) * (1.0d0 - c1)
#endif
    end select
  end subroutine robin_bc_for_alppsi

  subroutine robin_bc_in_non_spherical(mg, id, nc, iv, nb, bc_type)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer, intent(in)       :: nb      !< Neighbor direction
    integer, intent(in)       :: bc_type !< Type of b.c.

    real(dp)                  :: dx(NDIM)
    real(dp)                  :: x(NDIM)
    integer                   :: pm_one
    integer                   :: IJK
    real(dp)                  :: tmp

    ! the reason why we need to do such complicated LHS in alppsi-1 case is that,
    ! the stored psi (mg_itmp2) will not be stored at the ghost cell at each level,
    ! so we need to avoid using psi at the ghost cell.
    do i = 1, NDIM
       dx(i) = mg%dr(i, mg%boxes(id)%lvl) 
    end do

    select case (mg%geometry_type)
#if NDIM == 3
    case (mg_cartesian)
       if ( nb <= 2 ) then
          ! lowx or highx
          pm_one = mg_neighb_high_pm(nb) ! This gives a -1 or +1 
          if (pm_one < 0) then
             i = 1
          else
             i = nc
          end if
          x(1) = mg%boxes(id)%r_min(1) + dx(1) * (i-0.5_dp)
          do k = 1, nc
             x(3) = mg%boxes(id)%r_min(3) + dx(3) * (k-0.5_dp)
             do j = 1, nc
                x(2) = mg%boxes(id)%r_min(2) + dx(2) * (j-0.5_dp)
                tmp = get_tmp( mg%operator_type==mg_cfc_alp ) * mg%boxes(id)%cc(i+pm_one, j, k, iv) &
                      + x(2) * ( mg%boxes(id)%cc(i, j+1, k, iv) - mg%boxes(id)%cc(i, j-1, k, iv) ) / (2.0_dp * dx(2)) &
                      + x(3) * ( mg%boxes(id)%cc(i, j, k+1, iv) - mg%boxes(id)%cc(i, j, k-1, iv) ) / (2.0_dp * dx(3)) 
                tmp = - 2.0_dp * dx(1) / x(1) * ( tmp ) * pm_one
                mg%boxes(id)%cc(i+pm_one, j, k, iv) = tmp + mg%boxes(id)%cc(i-pm_one, j, k, iv)
             end do
          end do
       else if ( nb <= 4 ) then
          ! lowy or highy
          pm_one = mg_neighb_high_pm(nb) ! This gives a -1 or +1 
          if (pm_one < 0) then
             j = 1
          else
             j = nc
          end if
          x(2) = mg%boxes(id)%r_min(2) + dx(2) * (j-0.5_dp)
          do k = 1, nc
             x(3) = mg%boxes(id)%r_min(3) + dx(3) * (k-0.5_dp)
             do i = 1, nc
                x(1) = mg%boxes(id)%r_min(1) + dx(1) * (i-0.5_dp)
                tmp = get_tmp( mg%operator_type==mg_cfc_alp ) * mg%boxes(id)%cc(i, j+pm_one, k, iv) &
                      + x(1) * ( mg%boxes(id)%cc(i+1, j, k, iv) - mg%boxes(id)%cc(i-1, j, k, iv) ) / (2.0_dp * dx(1)) &
                      + x(3) * ( mg%boxes(id)%cc(i, j, k+1, iv) - mg%boxes(id)%cc(i, j, k-1, iv) ) / (2.0_dp * dx(3)) 
                tmp = - 2.0_dp * dx(2) / x(2) * ( tmp ) * pm_one
                mg%boxes(id)%cc(i, j+pm_one, k, iv) = tmp + mg%boxes(id)%cc(i, j-pm_one, k, iv)
             end do
          end do
       else if ( nb <= 6 ) then
          ! lowz or highz
          pm_one = mg_neighb_high_pm(nb) ! This gives a -1 or +1 
          if (pm_one < 0) then
             k = 1
          else
             k = nc
          end if
          x(3) = mg%boxes(id)%r_min(3) + dx(3) * (k-0.5_dp)
          do j = 1, nc
             x(2) = mg%boxes(id)%r_min(2) + dx(2) * (j-0.5_dp)
             do i = 1, nc
                x(1) = mg%boxes(id)%r_min(1) + dx(1) * (i-0.5_dp)
                tmp = get_tmp( mg%operator_type==mg_cfc_alp ) * mg%boxes(id)%cc(i, j, k+pm_one, iv) &
                      + x(1) * ( mg%boxes(id)%cc(i+1, j, k, iv) - mg%boxes(id)%cc(i-1, j, k, iv) ) / (2.0_dp * dx(1)) &
                      + x(2) * ( mg%boxes(id)%cc(i, j+1, k, iv) - mg%boxes(id)%cc(i, j-1, k, iv) ) / (2.0_dp * dx(2)) 
                tmp = - 2.0_dp * dx(3) / x(3) * ( tmp ) * pm_one
                mg%boxes(id)%cc(i, j, k+pm_one, iv) = tmp + mg%boxes(id)%cc(i, j, k-pm_one, iv)
             end do
          end do
       else
          error stop "robin_bc_in_non_spherical: shouldnt be here"
       end if
#endif

    case (mg_cylindrical)
#if NDIM == 2
       if ( nb <= 2 ) then
          ! lowx or highx
          pm_one = mg_neighb_high_pm(nb) ! This gives a -1 or +1 
          if (pm_one < 0) then
             i = 1
          else
             i = nc
          end if
          x(1) = mg%boxes(id)%r_min(1) + dx(1) * (i-0.5_dp)
          do j = 1, nc
             x(2) = mg%boxes(id)%r_min(2) + dx(2) * (j-0.5_dp)
             tmp = get_tmp( mg%operator_type==mg_cfc_alp ) * mg%boxes(id)%cc(i+pm_one, j, iv) &
                   + x(2) * ( mg%boxes(id)%cc(i, j+1, iv) - mg%boxes(id)%cc(i, j-1, iv) ) / (2.0_dp * dx(2))
             tmp = - 2.0_dp * dx(1) / x(1) * ( tmp ) * pm_one
             mg%boxes(id)%cc(i+pm_one, j, iv) = tmp + mg%boxes(id)%cc(i-pm_one, j, iv)
          end do
       else if ( nb <= 4 ) then
          ! lowy or highy
          pm_one = mg_neighb_high_pm(nb) ! This gives a -1 or +1 
          if (pm_one < 0) then
             j = 1
          else
             j = nc
          end if
          x(2) = mg%boxes(id)%r_min(2) + dx(2) * (j-0.5_dp)
          do i = 1, nc
             x(1) = mg%boxes(id)%r_min(1) + dx(1) * (i-0.5_dp)
             tmp = get_tmp( mg%operator_type==mg_cfc_alp ) * mg%boxes(id)%cc(i, j+pm_one, iv) &
                     + x(1) * ( mg%boxes(id)%cc(i+1, j, iv) - mg%boxes(id)%cc(i-1, j, iv) ) / (2.0_dp * dx(1))
             tmp = - 2.0_dp * dx(2) / x(2) * ( tmp ) * pm_one
             mg%boxes(id)%cc(i, j+pm_one, iv) = tmp + mg%boxes(id)%cc(i, j-pm_one, iv)
          end do
       else
          error stop "robin_bc_in_non_spherical: shouldnt be here"
       end if

#elif NDIM == 3
       if ( nb <= 2 ) then
          ! lowx or highx
          pm_one = mg_neighb_high_pm(nb) ! This gives a -1 or +1 
          if (pm_one < 0) then
             i = 1
          else
             i = nc
          end if
          x(1) = mg%boxes(id)%r_min(1) + dx(1) * (i-0.5_dp)
          do k = 1, nc
             !x(3) = mg%boxes(id)%r_min(3) + dx(3) * (k-0.5_dp)
             do j = 1, nc
                x(2) = mg%boxes(id)%r_min(2) + dx(2) * (j-0.5_dp)
                tmp = get_tmp( mg%operator_type==mg_cfc_alp ) * mg%boxes(id)%cc(i+pm_one, j, k, iv) &
                      + x(2) * ( mg%boxes(id)%cc(i, j+1, k, iv) - mg%boxes(id)%cc(i, j-1, k, iv) ) / (2.0_dp * dx(2))
                tmp = - 2.0_dp * dx(1) / x(1) * ( tmp ) * pm_one
                mg%boxes(id)%cc(i+pm_one, j, k, iv) = tmp + mg%boxes(id)%cc(i-pm_one, j, k, iv)
             end do
          end do
       else if ( nb <= 4 ) then
          ! lowy or highy
          pm_one = mg_neighb_high_pm(nb) ! This gives a -1 or +1 
          if (pm_one < 0) then
             j = 1
          else
             j = nc
          end if
          x(2) = mg%boxes(id)%r_min(2) + dx(2) * (j-0.5_dp)
          do k = 1, nc
             !x(3) = mg%boxes(id)%r_min(3) + dx(3) * (k-0.5_dp)
             do i = 1, nc
                x(1) = mg%boxes(id)%r_min(1) + dx(1) * (i-0.5_dp)
                tmp = get_tmp( mg%operator_type==mg_cfc_alp ) * mg%boxes(id)%cc(i, j+pm_one, k, iv) &
                      + x(1) * ( mg%boxes(id)%cc(i+1, j, k, iv) - mg%boxes(id)%cc(i-1, j, k, iv) ) / (2.0_dp * dx(1))
                tmp = - 2.0_dp * dx(2) / x(2) * ( tmp ) * pm_one
                mg%boxes(id)%cc(i, j+pm_one, k, iv) = tmp + mg%boxes(id)%cc(i, j-pm_one, k, iv)
             end do
          end do
       else
          error stop "robin_bc_in_non_spherical: shouldnt be here"
       end if
#endif
    case default
       error stop "robin_bc_in_non_spherical: unknown boundary condition"
    end select

    CONTAINS
       real(dp) function get_tmp(flag) 
          implicit none
          logical, intent(in)       :: flag
          if (flag) then
             ! flag is: mg%operator_type==mg_cfc_alp 
             ! LHS = (u+1)(2 - 1/psi) - psi
             get_tmp = ( mg%boxes(id)%cc(IJK, iv) + 1.0_dp ) &
                  * ( 2.0_dp - 1.0_dp / mg%boxes(id)%cc(IJK, mg_itmp2) ) &
                  - mg%boxes(id)%cc(IJK, mg_itmp2)
          else
             get_tmp = mg%boxes(id)%cc(IJK, iv) 
          end if
       end function get_tmp
  end subroutine robin_bc_in_non_spherical

  !> Fill ghost cells near refinement boundaries which preserves diffusive fluxes.
  subroutine sides_rb(box, nc, iv, nb, gc)
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer, intent(in)       :: nb !< Ghost cell direction
    !> Interpolated coarse grid ghost cell data (but not yet in the nb direction)
#if NDIM == 1
    real(dp), intent(in)      :: gc(1)
    integer                   :: di
#elif NDIM == 2
    real(dp), intent(in)      :: gc(nc)
    integer                   :: di, dj
#elif NDIM == 3
    real(dp), intent(in)      :: gc(nc, nc)
    integer                   :: di, dj, dk
#endif
    integer                   :: IJK, ix, dix

    if (mg_neighb_low(nb)) then
       ix = 1
       dix = 1
    else
       ix = nc
       dix = -1
    end if

    select case (mg_neighb_dim(nb))
#if NDIM == 1
    case (1)
       i = ix
       di = dix
       box%cc(i-di, iv) = (2 * gc(1) + box%cc(i, iv))/3.0_dp
#elif NDIM == 2
    case (1)
       i = ix
       di = dix

       do j = 1, nc
          dj = -1 + 2 * iand(j, 1)
          box%cc(i-di, j, iv) = 0.5_dp * gc(j) &
               + 0.75_dp * box%cc(i, j, iv) &
               - 0.25_dp * box%cc(i+di, j, iv)
       end do
    case (2)
       j = ix
       dj = dix
       do i = 1, nc
          di = -1 + 2 * iand(i, 1)
          box%cc(i, j-dj, iv) = 0.5_dp * gc(i) &
               + 0.75_dp * box%cc(i, j, iv) &
               - 0.25_dp * box%cc(i, j+dj, iv)
       end do
#elif NDIM == 3
    case (1)
       i = ix
       di = dix
       do k = 1, nc
          dk = -1 + 2 * iand(k, 1)
          do j = 1, nc
             dj = -1 + 2 * iand(j, 1)
             box%cc(i-di, j, k, iv) = 0.5_dp * gc(j, k) &
                  + 0.75_dp * box%cc(i, j, k, iv) &
                  - 0.25_dp * box%cc(i+di, j, k, iv)
          end do
       end do
    case (2)
       j = ix
       dj = dix
       do k = 1, nc
          dk = -1 + 2 * iand(k, 1)
          do i = 1, nc
             di = -1 + 2 * iand(i, 1)
             box%cc(i, j-dj, k, iv) = 0.5_dp * gc(i, k) &
                  + 0.75_dp * box%cc(i, j, k, iv) &
                  - 0.25_dp * box%cc(i, j+dj, k, iv)
          end do
       end do
    case (3)
       k = ix
       dk = dix
       do j = 1, nc
          dj = -1 + 2 * iand(j, 1)
          do i = 1, nc
             di = -1 + 2 * iand(i, 1)
             box%cc(i, j, k-dk, iv) = 0.5_dp * gc(i, j) &
                  + 0.75_dp * box%cc(i, j, k, iv) &
                  - 0.25_dp * box%cc(i, j, k+dk, iv)
          end do
       end do
#endif
    end select

  end subroutine sides_rb

end module m_ghost_cells
