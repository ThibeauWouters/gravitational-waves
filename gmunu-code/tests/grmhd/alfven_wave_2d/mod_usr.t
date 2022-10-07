module mod_usr
  use mod_physics
  use mod_grmhd
  use mod_eos_idealgas

  implicit none
  private

  ! aux output variables
  integer, protected  :: divB_

  logical, parameter                     ::  initialize_mag_field = .false.

  double precision, parameter :: rho_0 = 1.0d0
  double precision, parameter :: P_0 = 5.0d-1
  double precision, parameter :: eta = 1.0d0
  double precision, parameter :: B_0 = 1.0d0

  double precision, parameter :: theta = 0.0d0 * dpi / 4.0d0  ! note: this is not working yet

  double precision            :: k ! wave vector 2pi / L_x
  double precision            :: va2_p

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid
    usr_improve_initial_condition => rm_improve_initial_condition
    usr_process_adv_grid => cal_div_B
    usr_print_log        => printlog

    call set_coordinate_system("Cartesian_2.5D")

    call grmhd_activate()
    call eos_idealgas_activate()
    divB_ = var_set_auxvar('divB')

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm_init_one_grid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics
    use mod_grmhd
    use mod_eos

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim) ! primitive variables

    double precision :: rho_h(ixI^S) ! rho(1+eps) + P
    double precision :: va2(ixI^S) ! Alfven speed
    double precision :: tmp(ixI^S)

    integer  ::  ix^D, idir

    {^IFONED call mpistop("This is a multi-D MHD problem") }

    k = 2.0d0 * dpi / dabs( xprobmax1 - xprobmin1 )
    !k = 1.0d0
    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0
    w(ixO^S,beta(3)) = 0.0d0

    w(ixO^S,rho_) = rho_0
    w(ixO^S,press_) = P_0
    ! get eps
    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

    ! get Alfven speed
    call eos_get_eps_one_grid(P_0,rho_0,va2_p) ! eps
    va2_p = rho_0 * ( 1.0d0 + va2_p ) + P_0   ! rho_h
    va2_p = 2.0d0 * B_0**2 / ( va2_p + B_0**2 * (1.0d0 + eta**2) ) ! tmp := 2*B_0 / (rho_h + B_0^2(1+eta^2))
    va2_p = va2_p / ( 1.0d0 + dsqrt( 1.0d0 - eta**2 * va2_p**2 ) ) 

    select case (iprob)

    case (0)
       ! small amplitude cp alfven wave test
       w(ixO^S,Bvec(1))   = B_0 * eta * dcos(theta) &
                          - B_0 * eta * dcos( k * ( x(ixO^S,1) * dcos(theta) - x(ixO^S,2) * dsin(theta) ) ) * dsin(theta)
       w(ixO^S,Bvec(2))   = B_0 * eta * dsin(theta) &
                          + B_0 * eta * dcos( k * ( x(ixO^S,1) * dcos(theta) - x(ixO^S,2) * dsin(theta) ) ) * dcos(theta)
       w(ixO^S,Bvec(3))   = B_0 * eta * dsin( k * ( x(ixO^S,1) * dcos(theta) - x(ixO^S,2) * dsin(theta) ) )  
       if (stagger_grid) stop

    case default
       error stop "Unknown iprob"
    end select

    ! get W
    tmp(ixO^S) = 1.0d0 / dsqrt(1.0d0 - eta**2 * va2_p)
    ! get Wv
    w(ixO^S,W_vel(1))   = tmp(ixO^S) * dsqrt(va2_p) * eta * dcos( k * ( x(ixO^S,1) * dcos(theta) - x(ixO^S,2) * dsin(theta) ) ) * dsin(theta) 
    w(ixO^S,W_vel(2))   = - tmp(ixO^S) * dsqrt(va2_p) * eta * dcos( k * ( x(ixO^S,1) * dcos(theta) - x(ixO^S,2) * dsin(theta) ) ) * dcos(theta) 
    w(ixO^S,W_vel(3))   = - tmp(ixO^S) * dsqrt(va2_p) * eta * dsin( k * ( x(ixO^S,1) * dcos(theta) - x(ixO^S,2) * dsin(theta) ) )

    w(ixO^S, divB_) = 0.0d0

    ! update rest of the prim
    call phys_update_eos(ixI^L, ixO^L, w)

  end subroutine rm_init_one_grid

  !> before the main loop, we improve the initial data here
  subroutine rm_improve_initial_condition()
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_physics
    use mod_eos

    integer              :: ix^D
    integer              :: igrid, iigrid
    double precision                :: v(ixG^T,1:ndir)
    double precision                :: lfac(ixG^T)
    double precision                :: r(ixG^T)

    ! initial cleaning is needed
    if ( initialize_mag_field ) then
       call phys_initial_clean_divb(it,global_time)
    end if
  end subroutine rm_improve_initial_condition

  !> calculate div B after the advance 
  subroutine cal_div_B(igrid,level,ixI^L,ixO^L,qt,s)
    use mod_global_parameters
    use mod_grmhd_phys_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s
    double precision                :: divb(ixI^S)
    call grmhd_get_divb(ixI^L, ixO^L, s, divb(ixI^S))
    s%prim(ixO^S, divB_) = dabs( divb(ixO^S) )
  end subroutine cal_div_B

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
    double precision     :: vz_diff(ixG^T)
    double precision     :: vz_exact(ixG^T)
    double precision     :: lfac(ixG^T)

    integer, parameter   :: total_vz_exact = 1
    integer, parameter   :: L1_norm = 2
    integer, parameter   :: n_var = 2

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       vz_exact(ixM^T) = - dsqrt(va2_p) * eta * dsin( k * ps(igrid)%x(ixM^T,1) )

       call phys_get_lfac2(ixG^LL, ixM^LL, ps(igrid), lfac(ixG^T))
       lfac(ixM^T) = dsqrt( lfac(ixM^T) )
       vz_diff(ixM^T) = dabs( ps(igrid)%prim(ixM^T, W_vel(3)) / lfac(ixM^T) - vz_exact(ixM^T) )

       send_buffer(total_vz_exact) = send_buffer(total_vz_exact) + sum(dabs(vz_exact(ixM^T)))
       send_buffer(L1_norm) = send_buffer(L1_norm) + sum(vz_diff(ixM^T))

    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    recive_buffer(L1_norm) = recive_buffer(L1_norm) / recive_buffer(total_vz_exact)

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
