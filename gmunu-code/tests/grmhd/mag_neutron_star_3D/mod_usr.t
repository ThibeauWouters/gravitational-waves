module mod_usr
  use mod_physics
  use mod_grmhd
  use mod_eos_polytrope
  use mod_XNS
  use mod_multigrid_coupling
  use mod_cfc

  implicit none

  ! aux output variables
  integer  :: divB_ = -1
   
  logical, parameter  ::  initialize_mag_field = .False.
  logical, parameter  ::  initialize_metric = .False.
  logical, parameter  ::  evolve_metric = .False.

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos
    use mod_cfc

    character*128  ::  profile_path = '/shared/s2/chi-kit.cheong/Mag_NS_projects/initial_profile/T1K1'
    integer  ::  Nth, Nr   ! number of grid of the profile


    usr_init_one_grid => ns_init_one_grid
    usr_before_main_loop => ns_before_main_loop
    usr_process_adv_grid => cal_div_B
    usr_refine_grid => my_refine
    usr_print_log => printlog

    divB_ = var_set_auxvar('divB')

    call set_coordinate_system("spherical")
    call grmhd_activate()
    ! use atmosphere
    call eos_atmo_activate()
    call eos_polytrope_activate()

    call initialize_gmunu()

      call mod_XNS_read_profile(profile_path)

      if ( prof%Nth == 1) then 
         stop "! the profile is 1D."
      else
         if ( ndir < 2 ) stop "! Dimension = 1 but the profile is 2D."
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
    integer                         :: ix^D, idir

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
       call mod_XNS_map_2D(w(ix^D,Bvec(3)),x(ix^D,r_), x(ix^D,theta_),prof%b3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
    {end do^D&\} 

    ! fill atmosphere
    where (w(ix^S,rho_) <= small_rho_thr )
       !w(ix^S,rho_)   = small_rho
       !w(ix^S,eps_)   = small_eps
       !w(ix^S,press_)   = small_press
       !v(ix^S,1) = 0.0d0
       !v(ix^S,2) = 0.0d0
       !v(ix^S,3) = 0.0d0
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

    w(ix^S, divB_) = smalldouble

    ! reset eps
    {do ix^D = ix^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
    ! update rest of the prim
    call phys_update_eos(ixG^L, ix^L, w)
  end subroutine ns_init_one_grid
  
  !> deallocate profile varibles to save memories
  subroutine ns_before_main_loop()
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_physics
    use mod_eos
    use mod_cfc
    use mod_grmhd_phys_divb_mg

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
             ps(igrid)%prim(ix^D,W_vel(2)) = 0.0d0
             ps(igrid)%prim(ix^D,W_vel(3)) = 0.0d0
          end if
       {end do^D&\}
   
       ! update rest of the prim
       call phys_update_eos(ixG^LL, ixM^LL, ps(igrid)%prim(ixG^T,1:nprim))
    end do

    if (evolve_metric) then
       if ( initialize_metric ) call cfc_metric_init()
    end if
    if ( use_multigrid .and. initialize_mag_field ) then
       call grmhd_clean_divb_multigrid(clean_prim=.True.)
    end if
  end subroutine ns_before_main_loop

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

    ! core handling
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

  !> calculate div B after the advance 
  subroutine cal_div_B(igrid,level,ixI^L,ixO^L,qt,prim,x)
    use mod_global_parameters
    use mod_grmhd_phys_parameters

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: prim(ixI^S,1:nprim)

    double precision                :: b_vector(ixI^S,1:ndir)
    double precision                :: divb(ixI^S), Bvec2(ixI^S)
    integer                         :: idir

    do idir = 1, ndir
       b_vector(ixI^S,idir) = prim(ixI^S, psi_)**6 * prim(ixI^S, Bvec(idir))
    end do

    call grmhd_get_divb(ixI^L, ixO^L, b_vector(ixI^S,1:ndir), divb(ixI^S), divb_4thorder )
    prim(ixO^S, divB_) = dabs( divb(ixO^S)  )
  end subroutine cal_div_B

  subroutine printlog()
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_grmhd_phys_parameters
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
       call phys_get_lfac2(ps(igrid), ixG^LL, ixM^LL, lfac(ixG^T))
       lfac(ixM^T) = dsqrt( lfac(ixM^T) )
       send_buffer(total_mass) = send_buffer(total_mass) + &
          sum(lfac(ixM^T) * ps(igrid)%prim(ixM^T,psi_)**6*ps(igrid)%prim(ixM^T,rho_)*ps(igrid)%dvolume(ixM^T))

       ! calculate total volume
       send_buffer(total_volume) = send_buffer(total_volume) + sum(ps(igrid)%dvolume(ixM^T))

       ! get divB
       do idir = 1, ndir
          Bvector(ixG^T,idir) = ps(igrid)%prim(ixG^T,psi_)**6 * ps(igrid)%prim(ixG^T, Bvec(idir))
       end do
       call grmhd_get_divb(ixG^LL, ixM^LL, Bvector(ixG^T,1:ndir), divb(ixG^T), divb_4thorder )

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

end module mod_usr
