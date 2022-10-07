module mod_usr
  use mod_physics
  use mod_grhd
  use mod_eos_hybrid
  use mod_eos_polytrope
  use mod_map_profile
  use mod_multigrid_coupling
  use mod_cfc
  
  implicit none

  logical, parameter  ::  initialize_metric = .True.
  logical, parameter  ::  evolve_metric = .True.
  !logical, parameter  ::  initialize_metric = .False.
  !logical, parameter  ::  evolve_metric = .False.

  logical :: bounce = .False.

contains

  subroutine usr_init()
    !use mod_global_parameters
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos
    use mod_cfc

    character*128  ::  profile_path = './profile/s15s7b2.short'
    integer  ::  Nr   ! number of grid of the profile
    double precision :: rho_index
    

    usr_init_one_grid => ccsn_init_one_grid
    usr_improve_initial_condition => ccsn_improve_initial_condition
    usr_print_log => printlog
    usr_refine_grid => my_refine
    usr_process_adv_global => my_poststep_analysis

    call set_coordinate_system("spherical")
    call grhd_activate()
    !call eos_hybrid_activate()
    call eos_polytrope_activate
    ! use atmosphere
    call eos_atmo_activate()

    !call initialize_gmunu()

    call mod_read_profile(profile_path)

    ! setup atmosphere density from the profile
    call eos_initialize_atmo(minval(prof%rho))
    !write(*,*) small_rho_thr, small_rho/rho_gf

    ! to use metric solver, we need to activate multigrid
    use_multigrid = .true.
    if (evolve_metric) then
       ! to use metric solver, we need to activate multigrid
       call cfc_solver_activate()
    end if

  end subroutine usr_init

  ! Initialize one grid
  subroutine ccsn_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    use mod_global_parameters
    use mod_eos

    integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nprim)

    integer  ::  ix1

    w = 0.0d0

    ! initally flat spacetime
    w(ix^S, alp_) = 1.0d0
    w(ix^S, psi_) = 1.0d0

    {do ix^D = ixM^LLIM^D \}
       call mod_map(w(ix1,rho_),x(ix1,r_),prof%rho,prof%radius,prof%Nr)
       call mod_map(w(ix1,W_vel(1)),x(ix1,r_),prof%vel,prof%radius,prof%Nr)
    {end do^D&\}


    {do ix^D = ixM^LLIM^D \}
       ! get_eps
       call eos_get_eps_one_grid(w(ix1, press_),w(ix1, rho_),w(ix1, eps_))
       call eos_get_pressure_one_grid(w(ix1, press_),w(ix1, rho_),w(ix1, eps_))
       call eos_get_cs2_one_grid(w(ix1, cs2_),w(ix1, rho_),w(ix1, eps_))
    {end do^D&\}

    ! fill atmosphere
    ! here we assume that the smallest rho in the profile is the atmosphere
    where (w(ix^S,rho_) <= small_rho_thr )
       w(ix^S,rho_)   = small_rho
       w(ix^S,W_vel(1)) = 0.0d0
       w(ix^S,press_)     = small_press
       w(ix^S,eps_)     = small_eps
    end where

    !perturbation
    !w(ix^S, veloc(1)) = 1.0d-8
    !where (x(ix^S,rho_) <= 1.0d1 )
       !w(ix^S,veloc(1)) = -5.0d-5 * dsin(dpi * x(ix^S,1)/1.0d1)
    !end where

  end subroutine ccsn_init_one_grid

  !> before the main loop, we improve the initial data here
  subroutine ccsn_improve_initial_condition()
    use mod_global_parameters
    use mod_physics
    use mod_eos
    use mod_cfc

    integer              :: ix^D
    integer              :: igrid, iigrid

    ! deallocate profile varibles to save memories
    deallocate(prof%radius)
    deallocate(prof%mass)
    deallocate(prof%rho)
    deallocate(prof%press)
    deallocate(prof%temp)
    deallocate(prof%eps)
    deallocate(prof%vel)
    deallocate(prof%ye)
    deallocate(prof%omega)

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
             !ps(igrid)%prim(ix^D,W_vel(2)) = 0.0d0
             !ps(igrid)%prim(ix^D,W_vel(3)) = 0.0d0
          end if
       {end do^D&\}
   
       ! update rest of the prim
       call phys_update_eos(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T,1:nprim))
    end do

    ! initialize metric
    !if (evolve_metric) then
       if ( initialize_metric ) call cfc_metric_init()
    !end if
  end subroutine ccsn_improve_initial_condition

  subroutine my_poststep_analysis(iit, qt)
    use mod_global_parameters
    use mod_eos
    use mod_eos_hybrid, only: rho_nuc
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt

    integer              :: igrid, iigrid
    double precision     :: rho_max, rho_max_local

    if (bounce) return

    rho_max_local = 0.0d0
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       rho_max_local = max(rho_max_local, maxval(ps(igrid)%prim(ixM^T,rho_)) )
    end do

    call MPI_ALLREDUCE(rho_max_local, rho_max, 1, mpi_double_precision, &
          MPI_MAX, icomm, ierrmpi)

    if (rho_max >= rho_nuc * rho_gf ) then
       bounce = .True.
       if (mype == 0) then
          write(*,'(a,i7,a,es12.4)') ' core bounce at it=',it,' global_time=',global_time
          !write(*,'(a,i7,a,i7,a,es12.4)') ' save a snapshot No.',&
          !     snapshotnext,' at it=',it,' global_time=',global_time

          open(unit=666,file='t_bounce.dat',status="unknown",&
               form='formatted',position="append")
          write(666,*) 'it=',it,' global_time=',global_time
          close(666)

          open(unit=666,file='savenow',status="unknown",&
               form='formatted',position="append")
          write(666,*) ''
          close(666)
       end if
    end if

  end subroutine my_poststep_analysis

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
    !   refine = -1
    !   coarsen = 1
    else if ( level < to_level ) then
    !   refine = 1
    !   coarsen = -1
    end if

    ratio(ixO^S) = dx(1, level) / x(ixO^S,1) 

    if ( any( ratio(ixO^S) > 1.0d-2 ) ) then
       ! the resolution is too low
       refine = 1
       coarsen = -1
    else if ( all( ratio(ixO^S) < 1.0d-2 ) ) then
       ! the resolution is too high
       refine = -1
       coarsen = 1
    end if
  end subroutine my_refine

  subroutine printlog()
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    logical              :: fileopen
    integer              :: i, iw, level
    double precision     :: wmean(0:nprim), wbuffer(0:nprim), total_volume
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
    double precision     :: tmp(1:2), mass_buffer, mass
    double precision     :: lfac(ixG^T)
    double precision     :: rho_max, rho_max_local
    double precision     :: alp_min, alp_min_local

    mass_buffer = 0.0d0

    ! calculate total rest mass
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       call phys_get_lfac2(ps(igrid), ixG^LL, ixM^LL, lfac(ixG^T))
       lfac(ixM^T) = dsqrt( lfac(ixM^T) )
       mass_buffer = mass_buffer + &
          sum(lfac(ixM^T) * ps(igrid)%prim(ixM^T,psi_)**6*ps(igrid)%prim(ixM^T,rho_)*ps(igrid)%dvolume(ixM^T))
    end do

    ! gather all the results!
    call MPI_ALLREDUCE(mass_buffer, mass, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

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
          line = "it global_time dt"
          i = len_trim(line) + 2
          write(line(i:),"(a,a)") trim("Mass"), " "
          write(line(i:),"(a,a)") trim("rho_max"), " "
          write(line(i:),"(a,a)") trim("alp_min"), " "

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
       write(line(i:), fmt_string) mass
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', 1, fmt_r // ')'
       write(line(i:), fmt_string) rho_max
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', 1, fmt_r // ')'
       write(line(i:), fmt_string) alp_min
       i = len_trim(line) + 2

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

  end subroutine printlog

end module mod_usr
