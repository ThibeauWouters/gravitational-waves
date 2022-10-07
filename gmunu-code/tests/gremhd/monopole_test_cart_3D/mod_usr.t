module mod_usr
  use mod_physics
  use mod_gremhd
  use mod_eos_idealgas

  implicit none

  ! aux output variables
  integer  :: divB_ = -1

  double precision, parameter            :: rg = 0.2d0

  double precision, parameter            :: resistivity = 1.0d-16
   
contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid
    usr_improve_initial_condition => rm_improve_initial_condition
    usr_process_adv_grid => cal_div_B
    usr_get_resistivity => my_get_resistivity
    usr_print_log => printlog

    call set_coordinate_system("Cartesian")
    call gremhd_activate()
    call eos_idealgas_activate()

    divB_ = var_set_auxvar('divB')

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
    double precision                :: v(ixI^S,1:ndir)
    double precision                :: r(ixI^S)

    integer  ::  ix^D, idir, max_loc(1:ndir)
    integer  ::  ixC^L

    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0
    w(ixO^S,beta(3)) = 0.0d0

    w(ixO^S,rho_)       = 1.0d0
    w(ixO^S,eps_)       = 1.0d-1
    w(ixO^S,press_)     = 2.0d0 / 3.0d0

    w(ixO^S,Bphi_)      = 0.0d0
    w(ixO^S,Bvec(1))    = 0.0d0
    w(ixO^S,Bvec(2))    = 0.0d0
    w(ixO^S,Bvec(3))    = 0.0d0

    r(ixO^S) = dsqrt( x(ixO^S,1)**2 + x(ixO^S,2)**2  + x(ixO^S,3)**2 )

    ! Gaussian monopoles
    where ( r(ixO^S) < rg )
       w(ixO^S,Bvec(1))    = dexp( - r(ixO^S)**2 / rg**2 ) - dexp(-1.0d0)
    end where

    w(ixO^S, divB_) = 0.0d0

    {do ix^D = ixO^LIM^D \}
       call eos_get_pressure_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_cs2_one_grid(w(ix^D, cs2_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
  end subroutine rm_init_one_grid
  
  !> before the main loop, we improve the initial data here
  subroutine rm_improve_initial_condition()
    use mod_global_parameters
    use mod_gremhd_phys_parameters
    integer              :: igrid, iigrid
    double precision                :: divb(ixG^T)
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       call gremhd_get_divb(ixG^LL, ixM^LL, ps(igrid), divb(ixG^T))
       ps(igrid)%prim(ixM^T, divB_) = ( divb(ixM^T) )
    end do
  end subroutine rm_improve_initial_condition

  !> calculate div B after the advance 
  subroutine cal_div_B(igrid,level,ixI^L,ixO^L,qt,s)
    use mod_global_parameters
    use mod_gremhd_phys_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s
    double precision                :: divb(ixI^S)
    call gremhd_get_divb(ixI^L, ixO^L, s, divb(ixI^S))
    s%prim(ixO^S, divB_) = ( divb(ixO^S) )
  end subroutine cal_div_B

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

    integer, parameter   :: total_mass = 1
    integer, parameter   :: total_divBdV = 2
    integer, parameter   :: total_divB2dV = 3
    integer, parameter   :: total_volume = 4
    integer, parameter   :: n_var = 4

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! calculate total rest mass
       call phys_get_lfac2(ixG^LL, ixM^LL, ps(igrid), lfac(ixG^T))
       lfac(ixM^T) = dsqrt( lfac(ixM^T) )
       send_buffer(total_mass) = send_buffer(total_mass) + &
          sum(lfac(ixM^T) * ps(igrid)%prim(ixM^T,psi_)**6*ps(igrid)%prim(ixM^T,rho_)*ps(igrid)%dvolume(ixM^T))

       ! calculate total volume
       send_buffer(total_volume) = send_buffer(total_volume) + sum(ps(igrid)%dvolume(ixM^T))

       ! calculate total divB dV
       send_buffer(total_divBdV) = send_buffer(total_divBdV) + &
          sum(dabs(ps(igrid)%prim(ixM^T, divB_)) * ps(igrid)%dvolume(ixM^T))

       ! calculate total divB^2 dV
       send_buffer(total_divB2dV) = send_buffer(total_divB2dV) + &
          sum(ps(igrid)%prim(ixM^T, divB_)**2 * ps(igrid)%dvolume(ixM^T))
    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    ! Volume average of divB
    recive_buffer(total_divBdV) = (recive_buffer(total_divBdV) / recive_buffer(total_volume))
    ! L2 norm of divB
    recive_buffer(total_divB2dV) = dsqrt( recive_buffer(total_divB2dV) / recive_buffer(total_volume) )

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
          line = "it global_time dt M0 avg_divB L2_divB"

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

       write(fmt_string, '(a,i0,a)') '(', n_var-1, fmt_r // ')'
       write(line(i:), fmt_string) recive_buffer(1:n_var-1)
       i = len_trim(line) + 2

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

    deallocate(send_buffer)
    deallocate(recive_buffer)

  end subroutine printlog

  ! get resistivity
  subroutine my_get_resistivity(ixI^L,ixO^L,cons,x,eta)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(in) :: cons(ixI^S,1:ncons) 
    double precision, intent(out):: eta(ixI^S) 
    eta(ixO^S) = resistivity
  end subroutine my_get_resistivity


end module mod_usr
