module mod_usr
  use mod_physics
  use mod_grhd
  use mod_eos_idealgas
  use mod_XNS
  use mod_multigrid_coupling
  use mod_cfc

  implicit none
   
  logical, parameter  ::  initialize_metric = .True.
  logical, parameter  ::  evolve_metric = .True.
  !logical, parameter  ::  initialize_metric = .False.
  !logical, parameter  ::  evolve_metric = .False.

  logical, parameter  ::  read_id_after_refine = .True.

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos
    use mod_cfc

    !character*128  ::  profile_path = './BU8_640x64'
    character*128  ::  profile_path = './BU0_320x32'
    integer  ::  Nth, Nr   ! number of grid of the profile


    usr_init_one_grid => ns_init_one_grid
    usr_improve_initial_condition => ns_improve_initial_condition
    usr_refine_grid => my_refine
    usr_print_log => printlog

    call set_coordinate_system("spherical")
    call grhd_activate()
    ! use atmosphere
    call eos_atmo_activate()
    call eos_idealgas_activate()

    call initialize_gmunu()

    call mod_XNS_read_profile(profile_path)

    if ( prof%Nth == 1) then 
       stop "! the profile is 1D."
    else
       if ( ndir < 2 ) stop "! Dimension = 1 but the profile is 2D."
       !if ( prof%Nth .ne. domain_nx2) stop "! Nth of the profile must be the same as domain_nx2!"
    end if

    ! setup atmosphere density from the profile
    call eos_initialize_atmo(maxval(prof%rho))

    if (evolve_metric) then
       ! to use metric solver, we need to activate multigrid
       call cfc_solver_activate()
    end if

  end subroutine usr_init

  ! Initialize one grid
  subroutine ns_init_one_grid(ixG^L,ix^L,w,x)
    use mod_global_parameters
    use mod_eos
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nprim)

    double precision                :: v(ixG^S,1:ndir)
    double precision                :: lfac(ixG^S)
    double precision                :: gamma(ixG^S,1:3,1:3)
    integer  ::  ix^D, idir

    {^IFONED call mpistop("This is a multi-D HD problem") }

    w = 0.0d0
    v = 0.0d0

    {do ix^D = ix^LIM^D \} 
       call mod_XNS_map_2D(w(ix^D,rho_),x(ix^D,r_), x(ix^D,theta_),prof%rho(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,press_),x(ix^D,r_), x(ix^D,theta_),prof%press(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,alp_),x(ix^D,r_), x(ix^D,theta_),prof%alp(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,psi_),x(ix^D,r_), x(ix^D,theta_),prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,beta(3)),x(ix^D,r_), x(ix^D,theta_),prof%beta3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(v(ix^D,3),x(ix^D,r_), x(ix^D,theta_),prof%v3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
    {end do^D&\} 

    ! fill atmosphere
    where (w(ix^S,rho_) <= small_rho_thr )
       w(ix^S,rho_)   = small_rho
       w(ix^S,eps_)   = small_eps
       w(ix^S,press_)   = small_press
       v(ix^S,3) = 0.0d0
    end where

    ! get the metric
    call get_gamma_ij_hat(x(ixG^S, 1:ndim), ixG^L, ix^L, gamma(ixG^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ix^S,idir,idir) = gamma(ix^S,idir,idir) * w(ix^S, psi_)**4 
    end do

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
    w(ix^S,W_vel(2))   = v(ix^S, 2) * lfac(ix^S)
    w(ix^S,W_vel(3))   = v(ix^S, 3) * lfac(ix^S)

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

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       if (read_id_after_refine) then
          call usr_init_one_grid(ixG^LL, ixM^LL, ps(igrid)%prim, ps(igrid)%x)
       end if

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
             ps(igrid)%prim(ix^D,W_vel(2)) = 0.0d0
             ps(igrid)%prim(ix^D,W_vel(3)) = 0.0d0
          end if
       {end do^D&\}
   
       ! update rest of the prim
       call phys_update_eos(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T,1:nprim))
    end do

    ! initialize metric
    if (evolve_metric) then
       if ( initialize_metric ) call cfc_metric_init()
    end if

    ! deallocate profile varibles to save memories
    call mod_XNS_deallocate_var()

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

    if ( minval(x(ixO^S,1)) <= 5.0d-1 ) then
       ratio(ixO^S) = dx(2, level) * x(ixO^S,1) / dx(1, level)
       ratio(ixO^S) = ( ratio(ixO^S) - 1.0d0 )
       if ( any( ratio(ixO^S) <= 0.0d0 ) ) then
          refine = -1
          coarsen = 1
       else
          !refine = 1
          !coarsen = -1
       end if
    end if
  end subroutine my_refine

  subroutine printlog()
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_grhd_phys_parameters

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
    double precision     :: lfac(ixG^T), dV(ixG^T), h(ixG^T)
    double precision     :: tmp_1(ixG^T), tmp_2(ixG^T), tmp_3(ixG^T)
    double precision     :: y22_real(ixG^T), y22_img(ixG^T)
    double precision     :: e1_real(ixG^T), e1_img(ixG^T)
    double precision     :: e2_real(ixG^T), e2_img(ixG^T)
    double precision     :: e3_real(ixG^T), e3_img(ixG^T)
    double precision     :: e4_real(ixG^T), e4_img(ixG^T)

    double precision     :: rho_max, rho_max_local
    double precision     :: alp_min, alp_min_local

    integer              :: total_volume = -1
    integer              :: total_mass = -1
    integer              :: J_rot = -1
    integer              :: T_rot = -1
    integer              :: total_I11dot = -1
    integer              :: total_I12dot = -1
    integer              :: total_I13dot = -1
    integer              :: total_I22dot = -1
    integer              :: total_I23dot = -1
    integer              :: total_I33dot = -1
    integer              :: c1_real = -1
    integer              :: c1_img = -1
    integer              :: c2_real = -1
    integer              :: c2_img = -1
    integer              :: c3_real = -1
    integer              :: c3_img = -1
    integer              :: c4_real = -1
    integer              :: c4_img = -1
    integer              :: rho22_real = -1
    integer              :: rho22_img = -1
    integer              :: n_var = -1

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

    ! initialize variables
    n_var = 0
    total_volume = add_var()
    total_mass = add_var()
    c1_real = add_var()
    c1_img = add_var()
    c2_real = add_var()
    c2_img = add_var()
    c3_real = add_var()
    c3_img = add_var()
    c4_real = add_var()
    c4_img = add_var()
    rho22_real = add_var()
    rho22_img = add_var()
    J_rot = add_var()
    T_rot = add_var()
    total_I11dot = add_var()
    total_I12dot = add_var()
    total_I13dot = add_var()
    total_I22dot = add_var()
    total_I23dot = add_var()
    total_I33dot = add_var()

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! calculate total volume
       dV(ixM^T) = ps(igrid)%prim(ixM^T,psi_)**6*ps(igrid)%dvolume(ixM^T)
       send_buffer(total_volume) = send_buffer(total_volume) + sum( dV(ixM^T) )

       ! calculate total rest mass
       call phys_get_lfac2(ixG^LL, ixM^LL, ps(igrid), lfac(ixG^T))
       lfac(ixM^T) = dsqrt( lfac(ixM^T) )
       send_buffer(total_mass) = send_buffer(total_mass) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) )

       ! calculate real part of total rho in m=1 azimuthal decomposition 
       e1_real(ixM^T) = dcos(1.0d0*ps(igrid)%x(ixM^T,3))
       send_buffer(c1_real) = send_buffer(c1_real) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e1_real(ixM^T) )

       ! calculate img part of total rho in m=1 azimuthal decomposition 
       e1_img(ixM^T) = dsin(1.0d0*ps(igrid)%x(ixM^T,3))
       send_buffer(c1_img) = send_buffer(c1_img) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e1_img(ixM^T) )


       ! calculate real part of total rho in m=2 azimuthal decomposition 
       e2_real(ixM^T) = dcos(2.0d0*ps(igrid)%x(ixM^T,3))
       send_buffer(c2_real) = send_buffer(c2_real) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e2_real(ixM^T) )

       ! calculate img part of total rho in m=2 azimuthal decomposition 
       e2_img(ixM^T) = dsin(2.0d0*ps(igrid)%x(ixM^T,3))
       send_buffer(c2_img) = send_buffer(c2_img) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e2_img(ixM^T) )


       ! calculate real part of total rho in m=3 azimuthal decomposition 
       e3_real(ixM^T) = dcos(3.0d0*ps(igrid)%x(ixM^T,3))
       send_buffer(c3_real) = send_buffer(c3_real) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e3_real(ixM^T) )

       ! calculate img part of total rho in m=3 azimuthal decomposition 
       e3_img(ixM^T) = dsin(3.0d0*ps(igrid)%x(ixM^T,3))
       send_buffer(c3_img) = send_buffer(c3_img) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e3_img(ixM^T) )


       ! calculate real part of total rho in m=4 azimuthal decomposition 
       e4_real(ixM^T) = dcos(4.0d0*ps(igrid)%x(ixM^T,3))
       send_buffer(c4_real) = send_buffer(c4_real) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e4_real(ixM^T) )

       ! calculate img part of total rho in m=4 azimuthal decomposition 
       e4_img(ixM^T) = dsin(4.0d0*ps(igrid)%x(ixM^T,3))
       send_buffer(c4_img) = send_buffer(c4_img) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e4_img(ixM^T) )
                                      
    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

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
          line = trim(line) // " c1_real c1_img c2_real c2_img"
          line = trim(line) // " c3_real c3_img c4_real c4_img"
          line = trim(line) // " rho22_real rho22_img J_rot T_rot"
          line = trim(line) // " I11_dot I12_dot I13_dot I22_dot I23_dot I33_dot"

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

       write(fmt_string, '(a,i0,a)') '(', n_var-1, fmt_r // ')'
       write(line(i:), fmt_string) recive_buffer(2:n_var)
       i = len_trim(line) + 2

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

    deallocate(send_buffer)
    deallocate(recive_buffer)
    
    contains

       function add_var() result (i_var)
          integer :: i_var
          n_var = n_var + 1
          i_var = n_var
       end function add_var

  end subroutine printlog

end module mod_usr
