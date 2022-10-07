module mod_usr
  use mod_physics
  use mod_grhd
  use mod_eos_polytrope
  use mod_XNS

  implicit none
   
  logical, parameter  ::  initialize_metric = .True.
  !logical, parameter  ::  initialize_metric = .False.

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos
    use mod_cfc

    character*128  ::  profile_path = './BU0_1D_640'
    integer  ::  Nth, Nr   ! number of grid of the profile


    usr_init_one_grid => ns_init_one_grid
    usr_improve_initial_condition => ns_improve_initial_condition
    usr_refine_grid => my_refine
    usr_print_log => printlog

    call set_coordinate_system("spherical")
    call grhd_activate()
    call eos_polytrope_activate()
    ! use atmosphere
    call eos_atmo_activate()

    !call initialize_gmunu()

    call mod_XNS_read_profile(profile_path)
    if ( prof%Nth == 1) then  
    ! This is a 1D profile.
       if ( ndir .ne. 1 ) stop "! Dimension = 2 but the profile is 1D."
    else
    ! This is a 2D profile.
       if ( ndir < 2 ) stop "! Dimension = 1 but the profile is 2D."
    endif

    ! setup atmosphere density from the profile
    call eos_initialize_atmo(maxval(prof%rho))

    ! to use metric solver, we need to activate multigrid
    call cfc_solver_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine ns_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    use mod_global_parameters
    use mod_eos
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nprim)

    double precision                :: v(ixG^S,1:ndir)
    double precision                :: lfac(ixG^S)
    double precision                :: gamma(ixG^S,1:3,1:3)
    double precision                :: r_star, delta(0:2) ! for perturbation
    integer  ::  ix^D, idir

    w(ixG^S,:) = 0.0d0
    v(ixG^S,:) = 0.0d0

    do ix1 =ixmin1, ixmax1 
       call mod_XNS_map_1D(w(ix1,rho_),x(ix1,r_),prof%rho(:,1),prof%radius,prof%Nr)
       call mod_XNS_map_1D(w(ix1,press_),x(ix1,r_),prof%press(:,1),prof%radius,prof%Nr)
       call mod_XNS_map_1D(w(ix1,alp_),x(ix1,r_),prof%alp(:,1),prof%radius,prof%Nr)
       call mod_XNS_map_1D(w(ix1,psi_),x(ix1,r_),prof%psi(:,1),prof%radius,prof%Nr)
    end do
    w(ixmin1:ixmax1,beta(1))   = 0.0d0

    ! get the metric
    call get_gamma_ij_hat(x(ixG^S, 1:ndim), ixG^L, ix^L, gamma(ixG^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ix^S,idir,idir) = gamma(ix^S,idir,idir) * w(ix^S, psi_)**4 
    end do

    !do ix1 =ixmin1, ixmax1 
       ! for l=0 mode, here we perturbe v^r
       !v(ix^D,1) = v(ix^D,1) &
       !         - 1.0d-2 / gamma(ix^D, 1, 1) * dsin( dpi * x(ix^D,r_) / prof%r_star(1) )
    !end do
    !w(ixG^S,rho_) = w(ixG^S,rho_) * 1.001d0

    ! get W
    lfac(ix^S) = 1.0d0
    ! calculate v^2 first
    do idir = 1, ndir
       lfac(ix^S) = lfac(ix^S) - gamma(ix^S,idir,idir) * v(ix^S, idir)**2
    end do
    ! Calculate the Lorentz factor from velocity
    lfac(ix^S) = dsqrt( 1.0d0 / lfac(ix^S) )

    ! set veloc -> W_vel
    w(ix^S,W_vel(1))   = v(ix^S, 1) * lfac(ix^S)

    ! reset eps
    {do ix^D = ix^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
    ! update rest of the prim
    call phys_update_eos(ixG^L, ix^L, w)

  end subroutine ns_init_one_grid
  
  !> before the main loop, we improve the initial data here
  subroutine ns_improve_initial_condition()
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_physics
    use mod_eos
    use mod_cfc

    integer              :: ix^D
    integer              :: igrid, iigrid

    ! deallocate profile varibles to save memories
    call mod_XNS_deallocate_var()

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       {do ix^D = ixM^LLIM^D \}
          ! reset eps
          call eos_get_eps_one_grid(ps(igrid)%prim(ix^D, press_),ps(igrid)%prim(ix^D, rho_),ps(igrid)%prim(ix^D, eps_))
          ps(igrid)%prim(ix^D, eps_) = max( ps(igrid)%prim(ix^D, eps_), small_eps )

          ! fill atmosphere
          if (ps(igrid)%prim(ix^D,rho_) <= small_rho_thr ) then
             ps(igrid)%prim(ix^D,rho_)   = small_rho
             ps(igrid)%prim(ix^D,eps_)   = small_eps
             ps(igrid)%prim(ix^D,press_)   = small_press
             ps(igrid)%prim(ix^D,W_vel(1)) = 0.0d0
          end if
       {end do^D&\}
   
       ! update rest of the prim
       call phys_update_eos(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T,1:nprim))
    end do

    ! initialize metric
    if ( initialize_metric ) call cfc_metric_init()
  end subroutine ns_improve_initial_condition

  subroutine my_refine(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nprim), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    double precision             :: ratio(ixI^S)
    integer                      :: i_level, to_level
    double precision             :: phi_grav ! relativistic gravitational potential :=  1 - alp
    double precision             :: phi_grav_cut
    double precision, parameter  :: phi_grav_max = 0.2d0

    refine = 0
    coarsen = 0

    phi_grav = minval( w(ixO^S, alp_) )
    phi_grav = 1.0d0 - phi_grav
 
    to_level = refine_max_level
    phi_grav_cut = phi_grav_max 
    do i_level = refine_max_level-1, 1, -1
       phi_grav_cut = phi_grav_cut / 2.0d0
       if ( phi_grav < phi_grav_cut ) then
          to_level = i_level
       end if
    end do

    if ( level > to_level ) then
       refine = -1
       coarsen = 1
    else if ( level < to_level ) then
       refine = 1
       coarsen = -1
    end if

  end subroutine my_refine

  subroutine printlog()
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    logical              :: fileopen
    integer              :: i, iw, level
    double precision, allocatable     :: send_buffer(:), recive_buffer(:)
    double precision     :: volume_coverage(refine_max_level)
    integer              :: nx^D, nc, ncells, dit
    double precision     :: dtTimeLast, now, cellupdatesPerSecond
    double precision     :: activeBlocksPerCore, wctPerCodeTime, timeToFinish
    character(len=40)    :: fmt_string
    character(len=80)    :: filename
    character(len=2048)  :: line
    logical, save        :: opened  = .false.
    integer              :: amode, istatus(MPI_STATUS_SIZE)
    integer, parameter   :: my_unit = 20

    integer              :: igrid, iigrid
    double precision     :: dV(ixG^T)
    double precision     :: lfac(ixG^T)

    double precision     :: rho_max, rho_max_local
    double precision     :: alp_min, alp_min_local

    integer, parameter   :: total_volume = 1
    integer, parameter   :: total_mass = 2
    integer, parameter   :: n_var = 2

    ! find max rho
    rho_max_local = 0.0d0
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       rho_max_local = max(rho_max_local, maxval(ps(igrid)%prim(ixM^T,rho_)) )
    end do
    call MPI_ALLREDUCE(rho_max_local, rho_max, 1, mpi_double_precision, &
          MPI_MAX, icomm, ierrmpi)

    ! find min alpha
    alp_min_local = 1.0d0
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       alp_min_local = min(alp_min_local, minval(ps(igrid)%prim(ixM^T,alp_)) )
    end do
    call MPI_ALLREDUCE(alp_min_local, alp_min, 1, mpi_double_precision, &
          MPI_MIN, icomm, ierrmpi)

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! make sure all conserved variables are updated, maybe no needed
       call phys_to_conserved(ixG^LL,ixM^LL,ps(igrid)%cons(ixG^T,1:ncons),ps(igrid)%prim(ixG^T,1:nprim),ps(igrid)%x(ixG^T,1:ndim))

       ! calculate total volume
       dV(ixM^T) = ps(igrid)%prim(ixM^T,psi_)**6*ps(igrid)%dvolume(ixM^T)
       send_buffer(total_volume) = send_buffer(total_volume) + sum( dV(ixM^T) )

       ! calculate total rest mass
       send_buffer(total_mass) = send_buffer(total_mass) + &
          sum( ps(igrid)%cons(ixM^T,D_) * ps(igrid)%dvolume(ixM^T) )

    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    ! NOTE: 4 pi factor is not included
    recive_buffer = recive_buffer * 4.0d0 * dpi

    if (mype == 0) then

       ! On first entry, open the file and generate the header
       if (.not. opened) then

          filename = trim(base_filename) // ".log"

          ! Delete the log when not doing a restart run
          if (restart_from_file == undefined) then
             open(unit=my_unit,file=trim(filename),status='replace')
             close(my_unit, status='delete')
          end if

          amode    = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
          amode    = ior(amode,MPI_MODE_APPEND)

          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, &
               MPI_INFO_NULL, log_fh, ierrmpi)

          opened = .true.

          ! Start of file headern
          line = "it global_time dt rho_max alp_min M_rest"
          i = len_trim(line) + 2

          ! Only write header if not restarting
          if (restart_from_file == undefined) then
            call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'), &
                 len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
          end if
       end if

       ! Construct the line to be added to the log

       fmt_string = '(' // fmt_i // ',2' // fmt_r // ')'
       write(line, fmt_string) it, global_time, dt
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', 1, fmt_r // ')'
       write(line(i:), fmt_string) rho_max
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', 1, fmt_r // ')'
       write(line(i:), fmt_string) alp_min
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', n_var, fmt_r // ')'
       write(line(i:), fmt_string) recive_buffer(2:n_var)
       i = len_trim(line) + 2

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

    deallocate(send_buffer)
    deallocate(recive_buffer)


  end subroutine printlog

end module mod_usr
