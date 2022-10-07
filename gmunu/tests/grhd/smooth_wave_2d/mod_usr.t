module mod_usr
  use mod_physics
  use mod_grhd
  use mod_eos_idealgas

  implicit none
  private

  double precision, parameter :: A = 2.0d-1
  double precision, parameter :: v0 = 2.0d-1
  double precision, parameter :: theta = dpi / 6.0d0

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid
    usr_print_log        => printlog

    call set_coordinate_system("Cartesian_2D")

    call grhd_activate()
    call eos_idealgas_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm_init_one_grid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics
    use mod_grhd
    use mod_eos

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim) ! primitive variables

    integer  ::  ix^D, idir, max_loc(1:ndir)

    {^IFONED call mpistop("This is a multi-D HD problem") }

    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0

    w(ixO^S,rho_) = 1.0d0 + A * dsin( 2.0d0 * dpi &
               * ( x(ixO^S,1) * dcos(theta) + x(ixO^S,2) * dsin(theta) ) )
    w(ixO^S,press_) = 1.0d0
    ! get eps
    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

    ! get Wv
    w(ixO^S,W_vel(1))   = v0 / dsqrt(1.0d0 - v0**2)

    ! update rest of the prim
    call phys_update_eos(ixI^L, ixO^L, w)

  end subroutine rm_init_one_grid

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
    double precision     :: rho_diff(ixG^T)
    double precision     :: rho_exact(ixG^T)

    integer, parameter   :: total_rho_exact = 1
    integer, parameter   :: L1_norm = 2
    integer, parameter   :: n_var = 2

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       rho_exact(ixM^T) = 1.0d0 + A * dsin( 2.0d0 * dpi &
               * ( (ps(igrid)%x(ixM^T,1) * dcos(theta) + ps(igrid)%x(ixM^T,2) * dsin(theta)) &
                   - ( v0 * dcos(theta) ) * global_time ) )
       rho_diff(ixM^T) = dabs( ps(igrid)%prim(ixM^T,rho_) - rho_exact(ixM^T) )

       send_buffer(total_rho_exact) = send_buffer(total_rho_exact) + sum(dabs(rho_exact(ixM^T)))
       send_buffer(L1_norm) = send_buffer(L1_norm) + sum(rho_diff(ixM^T))

    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    recive_buffer(L1_norm) = recive_buffer(L1_norm) / recive_buffer(total_rho_exact)

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
          line = "it global_time dt L1_norm"

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
