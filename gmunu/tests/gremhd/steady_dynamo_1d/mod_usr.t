module mod_usr
  use mod_physics
  use mod_gremhd
  use mod_eos_idealgas
  private

  public :: usr_init

  double precision, parameter :: resistivity = 0.1d0
  double precision, parameter :: xi_para = 0.5d0
  double precision, parameter :: wavenumber_k = 5.0d0

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid
    usr_get_resistivity  => my_get_resistivity
    usr_get_dynamo_coeff => my_get_dynamo_coeff
    usr_print_log => printlog

    call set_coordinate_system("Cartesian_1.75D")

    call gremhd_activate()
    call eos_idealgas_activate()
    ! setup atmosphere density from the profile
    !call eos_initialize_atmo(1.0d0)

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm_init_one_grid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics
    use mod_gremhd
    use mod_eos

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim) ! primitive variables

    integer  ::  ix^D

    w = 0.0d0

    ! no gauge effects
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(:)) = 0.0d0

    w(ixO^S,rho_)       = smalldouble
    w(ixO^S,press_)     = smalldouble
    w(ixO^S,W_vel(:))   = 0.0d0 

    w(ixO^S,Bvec(1))   = 0.0d0
    w(ixO^S,Bvec(2))   = 0.1d0 * dsin( wavenumber_k * x(ixO^S,1) )
    w(ixO^S,Bvec(3))   =-0.1d0 * dcos( wavenumber_k * x(ixO^S,1) )


    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       !call eos_get_pressure_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_cs2_one_grid(w(ix^D, cs2_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

  end subroutine rm_init_one_grid

  ! get resistivity
  subroutine my_get_resistivity(ixI^L,ixO^L,cons,x,eta)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(in) :: cons(ixI^S,1:ncons) 
    double precision, intent(out):: eta(ixI^S) 
    eta(ixO^S) = resistivity
  end subroutine my_get_resistivity

  ! get dynamo term
  subroutine my_get_dynamo_coeff(ixI^L,ixO^L,cons,x,xi)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(in) :: cons(ixI^S,1:ncons) 
    double precision, intent(out):: xi(ixI^S) 
    xi(ixO^S) = xi_para
  end subroutine my_get_dynamo_coeff

  subroutine printlog()
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_gremhd_phys_parameters
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

    integer              :: igrid, iigrid, idir
    double precision     :: lfac(ixG^T)

    double precision     :: Bvector(ixG^T,1:ndir)
    double precision     :: divb(ixG^T)

    integer, parameter   :: maxBx = 1
    integer, parameter   :: maxBy = 2
    integer, parameter   :: maxBz = 3
    integer, parameter   :: n_var = 3

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       send_buffer(maxBx) = max( send_buffer(maxBx), maxval(ps(igrid)%prim(ixM^T, Bvec(1))) )
       send_buffer(maxBy) = max( send_buffer(maxBy), maxval(ps(igrid)%prim(ixM^T, Bvec(2))) )
       send_buffer(maxBz) = max( send_buffer(maxBz), maxval(ps(igrid)%prim(ixM^T, Bvec(3))) )
    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_MAX, icomm, ierrmpi)

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
          line = "it global_time dt maxBx maxBy maxBz"

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
       write(line(i:), fmt_string) recive_buffer(1:n_var)
       i = len_trim(line) + 2

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

    deallocate(send_buffer)
    deallocate(recive_buffer)

  end subroutine printlog

end module mod_usr
