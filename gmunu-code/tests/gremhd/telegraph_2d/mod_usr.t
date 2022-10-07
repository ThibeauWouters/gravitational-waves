module mod_usr
  use mod_physics
  use mod_gremhd
  use mod_eos_idealgas

  implicit none
  private

  ! aux output variables
  integer, protected  :: divB_ = -1
  integer, protected  :: divE_ = -1

  logical, parameter              :: initialize_mag_field = .true.

  double precision, parameter     :: sigma = 1.0d0
  double precision, parameter     :: B1 = 1.0d0

  double precision                :: mu
  double precision                :: k(1:3) ! wavevector

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid
    usr_improve_initial_condition => rm_improve_initial_condition
    usr_get_resistivity => my_get_resistivity
    usr_process_adv_grid => cal_div_B
    usr_print_log => printlog

    call set_coordinate_system("Cartesian_2.5D")

    call gremhd_activate()
    call eos_idealgas_activate()

    divB_ = var_set_auxvar('divB')
    divE_ = var_set_auxvar('divE')

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
    double precision                :: k_dot_x(ixI^S)
    double precision                :: E(ixI^S,1:ndir), B(ixI^S,1:ndir)
    double precision                :: kx, k2
    double precision                :: Lx, Ly
    double precision                :: alpha, sin_alpha, cos_alpha, tan_alpha

    integer  ::  ix^D, idir, max_loc(1:ndir)

    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0
    w(ixO^S,beta(3)) = 0.0d0

    Lx = dabs(xprobmax1 - xprobmin1)
    Ly = dabs(xprobmax2 - xprobmin2)

    tan_alpha = Lx / Ly
    sin_alpha = Lx / dsqrt( Ly**2 + Lx**2 )
    cos_alpha = Ly / dsqrt( Ly**2 + Lx**2 )

    kx = 2.0d0 * dpi / Lx
    k(1) = 1.0d0
    k(2) = tan_alpha
    k(3) = 0.0d0
    k(:) = kx * k(:)
    k2 = dot_product(k,k) 
    mu = dsqrt( k2 - 0.25 * sigma**2 )
    write(*,*) "T = 2pi / mu = ", 2.0d0 * dpi / mu

    k_dot_x(ixO^S) = 0.0d0
    do idir = 1, ndim
       k_dot_x(ixO^S) = k_dot_x(ixO^S) + k(idir) * x(ixO^S, idir)
    end do

    B = 0.0d0
    E = 0.0d0

    B(ixO^S, 3) = B1 * dexp( - 0.5d0 * sigma * global_time ) * dcos(k_dot_x(ixO^S))
    E(ixO^S, 2) = B1 * dexp( - 0.5d0 * sigma * global_time ) * &
                   ( mu / dsqrt(k2) * dcos(k_dot_x(ixO^S)) + 0.5d0 * sigma / dsqrt(k2) * dsin(k_dot_x(ixO^S)) )

    w(ixO^S, Evec(1)) = E(ixO^S,1) * cos_alpha - E(ixO^S,2) * sin_alpha
    w(ixO^S, Evec(2)) = E(ixO^S,2) * cos_alpha + E(ixO^S,1) * sin_alpha
    w(ixO^S, Evec(3)) = E(ixO^S,3)
    w(ixO^S, Bvec(1)) = B(ixO^S,1) * cos_alpha - B(ixO^S,2) * sin_alpha
    w(ixO^S, Bvec(2)) = B(ixO^S,2) * cos_alpha + B(ixO^S,1) * sin_alpha
    w(ixO^S, Bvec(3)) = B(ixO^S,3)

    w(ixO^S,rho_) = 1.0d12

    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_pressure_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_cs2_one_grid(w(ix^D, cs2_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}

  end subroutine rm_init_one_grid
  
  !> before the main loop, we improve the initial data here
  subroutine rm_improve_initial_condition()
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_physics
    use mod_gremhd_phys_parameters
    use mod_eos

    integer              :: ix^D
    integer              :: igrid, iigrid
    double precision                :: divb(ixG^T), dive(ixG^T)

    ! initial cleaning is needed
    if ( initialize_mag_field ) then
       call phys_initial_clean_divb(it,global_time)
    end if
    ! update divB
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       call gremhd_get_div(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T,Bvec(1:ndir)), divb(ixG^T), divb_4thorder )
       ps(igrid)%prim(ixM^T, divB_) = ( divb(ixM^T) )
       call gremhd_get_div(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T,Evec(1:ndir)), dive(ixG^T), dive_4thorder )
       ps(igrid)%prim(ixM^T, divE_) = ( dive(ixM^T) )
    end do
  end subroutine rm_improve_initial_condition

  !> calculate div B after the advance 
  subroutine cal_div_B(igrid,level,ixI^L,ixO^L,qt,s)
    use mod_global_parameters
    use mod_gremhd_phys_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s
    double precision                :: divb(ixI^S), dive(ixI^S)
    call gremhd_get_div(ixI^L, ixO^L, s%prim(ixI^S,Bvec(1:ndir)), divb(ixI^S), divb_4thorder )
    s%prim(ixO^S, divB_) = ( divb(ixO^S) )
    call gremhd_get_div(ixI^L, ixO^L, s%prim(ixI^S,Evec(1:ndir)), dive(ixI^S), dive_4thorder )
    s%prim(ixO^S, divE_) = ( dive(ixO^S) )
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

    double precision     :: Bvector(ixG^T,1:ndir)
    double precision     :: divb(ixG^T)
    double precision     :: k_dot_x(ixG^T)
    double precision     :: Bz_diff(ixG^T)
    double precision     :: Bz_exact(ixG^T)

    integer, parameter   :: total_volume = 1
    integer, parameter   :: L1_norm_Bz = 2
    integer, parameter   :: total_divBdV = 3
    integer, parameter   :: total_divB2dV = 4
    integer, parameter   :: n_var = 4

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       k_dot_x(ixM^T) = 0.0d0
       do idir = 1, ndim
          k_dot_x(ixM^T) = k_dot_x(ixM^T) + k(idir) * ps(igrid)%x(ixM^T, idir)
       end do
       Bz_exact(ixM^T) = B1 * dexp( - 0.5d0 * sigma * global_time ) * dcos(k_dot_x(ixM^T) - mu * global_time)
       Bz_diff(ixM^T) = dabs( ps(igrid)%prim(ixM^T, Bvec(3)) - Bz_exact(ixM^T) )

       send_buffer(L1_norm_Bz) = send_buffer(L1_norm_Bz) &
            + sum(Bz_diff(ixM^T) * ps(igrid)%dvolume(ixM^T))

       ! calculate total volume
       send_buffer(total_volume) = send_buffer(total_volume) + sum(ps(igrid)%dvolume(ixM^T))

       ! get divB
       do idir = 1, ndir
          Bvector(ixG^T,idir) = ps(igrid)%prim(ixG^T,psi_)**6 * ps(igrid)%prim(ixG^T, Bvec(idir))
       end do
       call gremhd_get_div(ixG^LL, ixM^LL, Bvector(ixG^T,1:ndir), divb(ixG^T), divb_4thorder )

       ! calculate total divB dV
       send_buffer(total_divBdV) = send_buffer(total_divBdV) + &
          sum(dabs(divB(ixM^T)) * ps(igrid)%dvolume(ixM^T))

       ! calculate total divB^2 dV
       send_buffer(total_divB2dV) = send_buffer(total_divB2dV) + &
          sum(divB(ixM^T)**2 * ps(igrid)%dvolume(ixM^T))
    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    ! Volume average of divB
    recive_buffer(total_divBdV) = (recive_buffer(total_divBdV) / recive_buffer(total_volume))
    ! L2 norm of divB
    recive_buffer(total_divB2dV) = dsqrt( recive_buffer(total_divB2dV) / recive_buffer(total_volume) )

    recive_buffer(L1_norm_Bz) = recive_buffer(L1_norm_Bz) / recive_buffer(total_volume)

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
          line = "it global_time dt L1_norm_Bz avg_divB L2_divB"

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
       write(line(i:), fmt_string) recive_buffer(2:n_var)
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
    eta(ixO^S) = 1.0d0 / sigma
  end subroutine my_get_resistivity

end module mod_usr
