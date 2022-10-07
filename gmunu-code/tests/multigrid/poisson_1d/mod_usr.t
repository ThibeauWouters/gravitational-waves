module mod_usr
  use mod_multigrid_coupling
  use mod_grhd
  use mod_eos_polytrope

  implicit none

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    use_multigrid = .true.
    usr_init_one_grid => mg_init_one_grid
    usr_process_global => compute_phi
    !usr_before_main_loop => compute_phi

    call set_coordinate_system("Cartesian")
    call grhd_activate()
    call eos_polytrope_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine mg_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    use mod_global_parameters
    use mod_eos

    integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nprim)

    integer  ::  ix1, n
    double precision, parameter :: rho0 = 1.0d0

    w = 0.0d0
    
    where (x(ix^S,1)<1.0d0)
      w(ix^S,rho_) = rho0*(1.0d0-x(ix^S,1)**2)
    end where

    do n = 1, mg_num_neighbors
       !mg%bc(n, mg_iphi)%boundary_cond => my_bc
       !mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
    end do
    !mg%bc(1, mg_iphi)%bc_type = mg_bc_neumann
    !mg%bc(1, mg_iphi)%bc_value = dpi
    !mg%bc(2, mg_iphi)%bc_value = 8.0d0 *dpi / 150.0d0

    !mg%smoother_type = mg_smoother_gsrb
    mg%smoother_type = mg_smoother_gs
    mg%n_cycle_down = 2
    mg%n_cycle_up = 2

    do ix1 =ixmin1, ixmax1 
       ! get_eps
       !call eos_get_pressure_one_grid(w(ix1, press_),w(ix1, rho_),w(ix1, eps_))
       !call eos_get_eps_one_grid(w(ix1, press_),w(ix1, rho_),w(ix1, eps_))
    end do
    w(ix^S,psi_) = 1.0d0
    w(ix^S,alp_) = 1.0d0
    !w(ix^S,alp_) = 4.0d0 * dpi *w(ix^S,rho_)
  end subroutine mg_init_one_grid

  !> To fill ghost cells near physical boundaries
  subroutine my_bc(box, nc, iv, nb, bc_type, bc)
    type(mg_box_t), intent(in)    :: box
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv      !< Index of variable
    integer, intent(in)           :: nb      !< Direction
    integer, intent(out)          :: bc_type !< Type of b.c.
    double precision, intent(out) :: bc(1)
    double precision              :: x(nc, 1)

    !call mg_get_face_coords(box, nb, nc, x)
    bc_type = mg_bc_neumann
    bc      = 0.0_dp
    !bc      = solution(x(:, 1))
  end subroutine my_bc

  subroutine compute_phi(iit,qt)
    use mod_input_output
    use mod_forest
    use mod_multigrid_coupling
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer                      :: id, iigrid, igrid
    double precision             :: max_res, wmax(nprim), modes(nprim), vol
    integer                        :: nc, lvl, it
    type(tree_node), pointer       :: pnode
    real(dp)                       :: fac
    logical                        :: do_restrict, do_gc

    !if (iit == 0) then
       !call mg_copy_to_tree(alp_, mg_irhs, .false., .false.)
       !call mg_phi_bc_store(mg)
    !end if

    !if (multigrid_cycle_name == 'fmg') then
    !   call mg_fas_fmg(mg, .true., max_res)
    !else
    !end if

    ! disable time update
    time_advance = .false.

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Include one layer of ghost cells on grid leaves
       !mg%boxes(id)%cc(0:nc+1, mg_iphi) = &
       !     ps(igrid)%prim(ixMlo1-1:ixMhi1+1, psi_)
       mg%boxes(id)%cc(0:nc+1, mg_iphi) = 0.0d0
       mg%boxes(id)%cc(0:nc+1, mg_irhs) = 0.5d0*(dsin(dpi*ps(igrid)%x(ixMlo1-1:ixMhi1+1, 1)) &
            +dsin(16.0d0*dpi*ps(igrid)%x(ixMlo1-1:ixMhi1+1, 1)))
    end do

    do it = 1, 50
    call mg_fas_vcycle(mg, max_res=max_res)
    write(*,*) max_res
    end do
    call mg_copy_from_tree(mg_iphi, psi_)

    !call get_global_maxima(wmax)
    !call get_volume_average(2, modes, vol)
    !if (mype == 0) print *, "CONV", it, max_res, wmax(i_err), sqrt(modes(i_err)), sqrt(modes(i_res))

    !stop "solved"
  end subroutine compute_phi

end module mod_usr
