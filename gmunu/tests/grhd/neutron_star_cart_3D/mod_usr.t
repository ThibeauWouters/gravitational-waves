module mod_usr
  use mod_physics
  use mod_grhd
  use mod_eos_idealgas
  use mod_XNS
  use mod_multigrid_coupling
  use mod_cfc

  implicit none
   
  logical, parameter  ::  initialize_metric = .True.
  !logical, parameter  ::  initialize_metric = .False.

  ! if enforcing z_symm during initialisation
  logical, parameter  ::  z_symm = .True.
  logical, parameter  ::  read_id_after_refine = .True.
  logical, parameter  ::  perturbation = .True.

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos
    use mod_cfc

    character*128  ::  profile_path = 'BU0_320x32'

    usr_init_one_grid => ns_init_one_grid
    usr_improve_initial_condition => ns_improve_initial_condition
    usr_refine_grid => my_refine
    usr_print_log => printlog

    call set_coordinate_system("Cartesian")
    call grhd_activate()
    ! use atmosphere
    call eos_atmo_activate()
    call eos_idealgas_activate()

    call initialize_gmunu()

    call mod_XNS_read_profile(profile_path)
  
    ! setup atmosphere density from the profile
    call eos_initialize_atmo(maxval(prof%rho))

    ! to use metric solver, we need to activate multigrid
    call cfc_solver_activate()

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
    double precision                :: vphi(ixG^S), betaphi(ixG^S)
    double precision                :: lfac(ixG^S)
    double precision                :: gamma(ixG^S,1:3,1:3)
    integer                         :: ix^D, idir
    double precision                :: rad, theta, rad_xy
    double precision                :: sin_phi, cos_phi
    double precision                :: sin_theta, cos_theta

    {^IFONED call mpistop("This is a multi-D HD problem") }

    w = 0.0d0
    v = 0.0d0

    {do ix^D = ix^LIM^D \} 
       rad = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
       theta = acos( x(ix^D, 3) / rad )
       rad_xy = dsqrt( sum(x(ix^D,1:2)**2) ) 
       sin_phi = x(ix^D, 2) / rad_xy
       cos_phi = x(ix^D, 1) / rad_xy
       sin_theta = rad_xy / rad 

       if ( z_symm .and. theta > 0.5d0 * dpi ) theta = dpi - theta

       call mod_XNS_map_2D(w(ix^D,rho_),rad, theta,prof%rho(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,press_),rad, theta,prof%press(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,alp_),rad, theta,prof%alp(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,psi_),rad, theta,prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(betaphi(ix^D),rad, theta,prof%beta3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(vphi(ix^D),rad, theta,prof%v3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)

       vphi(ix^D) = vphi(ix^D) * rad * sin_theta
       betaphi(ix^D) = betaphi(ix^D) * rad * sin_theta

       v(ix^D, 1) = - sin_phi * vphi(ix^D)
       v(ix^D, 2) =   cos_phi * vphi(ix^D)
       w(ix^D, beta(1)) = - sin_phi * betaphi(ix^D)
       w(ix^D, beta(2)) =   cos_phi * betaphi(ix^D)
    {end do^D&\} 

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
    if ( initialize_metric ) call cfc_metric_init()

    ! perturbation is added once everything is initialised
    if (perturbation) then
       do iigrid = 1, igridstail
          igrid = igrids(iigrid)
          call init_perturbation(ixG^LL, ixM^LL, ps(igrid)%prim, ps(igrid)%x)
       end do
    end if

    ! deallocate profile varibles to save memories
    call mod_XNS_deallocate_var()

  end subroutine ns_improve_initial_condition

  ! Initialize one grid
  subroutine init_perturbation(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_eos
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim)

    double precision                :: v(ixI^S,1:ndir)
    double precision                :: vphi(ixI^S), betaphi(ixI^S)
    double precision                :: lfac(ixI^S)
    double precision                :: gamma(ixI^S,1:3,1:3)
    integer                         :: ix^D, idir
    double precision                :: rad, theta, phi, rad_xy, theta_tmp
    double precision                :: sin_phi, cos_phi
    double precision                :: sin_theta, cos_theta
    double precision                :: pre_factor
    double precision                :: r_star, delta(0:2) !for perturbation
    double precision                :: dvtheta(ixI^S)

    delta = 0.0d0

    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_)**4 
    end do

    ! add perturbation on velocities

    ! first, we need to factor out W for current Wv to get velocities
    lfac(ixO^S) = 1.0d0
    do idir = 1, ndir
       lfac(ixO^S) = lfac(ixO^S) + gamma(ixO^S,idir,idir) * w(ixO^S, W_vel(idir))**2
    end do
    lfac(ixO^S) = dsqrt( 1.0d0 / lfac(ixO^S) )

    ! set veloc -> W_vel
    do idir = 1, ndir
       v(ixO^S, idir) = w(ixO^S,W_vel(idir)) / lfac(ixO^S)
    end do

    delta(0) = 0.0d-3
    delta(1) = 1.0d-2
    delta(2) = 0.0d-3

    {do ix^D = ixO^LIM^D \}           
       if ( w(ix^D, rho_) <= small_rho_thr ) cycle
       rad = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
       theta = dacos( x(ix^D, 3) / rad )
       rad_xy = dsqrt( sum(x(ix^D,1:2)**2) )  
       sin_phi = x(ix^D, 2) / rad_xy
       cos_phi = x(ix^D, 1) / rad_xy
       sin_theta = rad_xy / rad 
       cos_theta = x(ix^D, 3) / rad 

       theta_tmp = theta
       if ( z_symm .and. theta_tmp > 0.5d0 * dpi ) theta_tmp = dpi - theta_tmp
       call mod_XNS_map_1D(r_star,theta_tmp,prof%r_star(:),prof%theta,prof%Nth) ! firstly, get r_star(theta)

       if ( rad < r_star ) then
          pre_factor = dsin( dpi * rad / r_star )
       else
          pre_factor = 0.0d0
       end if

       ! for l=2 mode, here we perturbe v^theta               
       dvtheta(ix^D) = delta(1) / (rad*w(ix^D, psi_)**2)**2 * pre_factor &
                       * sin_theta * cos_theta

       dvtheta(ix^D) = dvtheta(ix^D) * rad 
                   
       v(ix^D,1) = v(ix^D,1) + dvtheta(ix^D) * cos_theta * cos_phi 
       v(ix^D,2) = v(ix^D,2) + dvtheta(ix^D) * cos_theta * sin_phi          
       v(ix^D,3) = v(ix^D,3) - dvtheta(ix^D) * sin_theta

    {end do^D&\} 

    ! get W
    lfac(ixO^S) = 1.0d0
    ! calculate v^2 first
    do idir = 1, ndir
       lfac(ixO^S) = lfac(ixO^S) - gamma(ixO^S,idir,idir) * v(ixO^S, idir)**2
    end do
    ! Calculate the Lorentz factor from velocity
    lfac(ixO^S) = dsqrt( 1.0d0 / lfac(ixO^S) )

    ! set veloc -> W_vel
    w(ixO^S,W_vel(1))   = v(ixO^S, 1) * lfac(ixO^S)
    w(ixO^S,W_vel(2))   = v(ixO^S, 2) * lfac(ixO^S)
    w(ixO^S,W_vel(3))   = v(ixO^S, 3) * lfac(ixO^S)

    ! the following part is needed only when the perturbation 
    ! is applied on rho, press or eps,
    ! reset eps
    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
    ! update rest of the prim
    call phys_update_eos(ixI^L, ixO^L, w)

  end subroutine init_perturbation

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
    double precision     :: rsq(ixG^T)
    double precision     :: y22_real(ixG^T), y22_img(ixG^T)
    double precision     :: phi(ixG^T)
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

       ! calculate phi 
       phi(ixM^T) = datan2(ps(igrid)%x(ixM^T,2) , ps(igrid)%x(ixM^T,1))

       ! calculate real part of total rho in m=1 azimuthal decomposition 
       e1_real(ixM^T) = dcos(1.0d0*phi(ixM^T))
       send_buffer(c1_real) = send_buffer(c1_real) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e1_real(ixM^T) )

       ! calculate img part of total rho in m=1 azimuthal decomposition 
       e1_img(ixM^T) = dsin(1.0d0*phi(ixM^T))
       send_buffer(c1_img) = send_buffer(c1_img) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e1_img(ixM^T) )


       ! calculate real part of total rho in m=2 azimuthal decomposition 
       e2_real(ixM^T) = dcos(2.0d0*phi(ixM^T))
       send_buffer(c2_real) = send_buffer(c2_real) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e2_real(ixM^T) )

       ! calculate img part of total rho in m=2 azimuthal decomposition 
       e2_img(ixM^T) = dsin(2.0d0*phi(ixM^T))
       send_buffer(c2_img) = send_buffer(c2_img) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e2_img(ixM^T) )


       ! calculate real part of total rho in m=3 azimuthal decomposition 
       e3_real(ixM^T) = dcos(3.0d0*phi(ixM^T))
       send_buffer(c3_real) = send_buffer(c3_real) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e3_real(ixM^T) )

       ! calculate img part of total rho in m=3 azimuthal decomposition 
       e3_img(ixM^T) = dsin(3.0d0*phi(ixM^T))
       send_buffer(c3_img) = send_buffer(c3_img) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e3_img(ixM^T) )


       ! calculate real part of total rho in m=4 azimuthal decomposition 
       e4_real(ixM^T) = dcos(4.0d0*phi(ixM^T))
       send_buffer(c4_real) = send_buffer(c4_real) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e4_real(ixM^T) )

       ! calculate img part of total rho in m=4 azimuthal decomposition 
       e4_img(ixM^T) = dsin(4.0d0*phi(ixM^T))
       send_buffer(c4_img) = send_buffer(c4_img) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * &
          e4_img(ixM^T) )
                                      
       ! calculate real part of total rho in Y_22 projection
       rsq(ixM^T) = (ps(igrid)%x(ixM^T,1))**2 + (ps(igrid)%x(ixM^T,2))**2 + (ps(igrid)%x(ixM^T,3))**2
       y22_real(ixM^T) = ((ps(igrid)%x(ixM^T,1))**2 - (ps(igrid)%x(ixM^T,2))**2)/rsq(ixM^T)
       send_buffer(rho22_real) = send_buffer(rho22_real) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * y22_real(ixM^T) )

       ! calculate imaginary part of total rho in Y_22 projection
       y22_img(ixM^T) = -2*ps(igrid)%x(ixM^T,1) * ps(igrid)%x(ixM^T,2)/rsq(ixM^T)
       send_buffer(rho22_img) = send_buffer(rho22_img) + &
          sum( lfac(ixM^T) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) * y22_img(ixM^T) )


       ! calculate angular momentum 
       h(ixM^T) = 1.0d0 + ps(igrid)%prim(ixM^T,eps_) + ps(igrid)%prim(ixM^T,press_) / ps(igrid)%prim(ixM^T,rho_)
       tmp_1(ixM^T) = - ps(igrid)%cons(ixM^T,D_) * h(ixM^T) * ps(igrid)%x(ixM^T, 2) * &
              ( ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,W_vel(1)) - ps(igrid)%prim(ixM^T,beta(1))*lfac(ixM^T) )
       tmp_1(ixM^T) = tmp_1(ixM^T) + ps(igrid)%cons(ixM^T,D_) * h(ixM^T) * ps(igrid)%x(ixM^T, 1) * &
              ( ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,W_vel(2)) - ps(igrid)%prim(ixM^T,beta(2))*lfac(ixM^T) )
       tmp_1(ixM^T) = ps(igrid)%prim(ixM^T,psi_)**4 * tmp_1(ixM^T)
       send_buffer(J_rot) = send_buffer(J_rot) + &
          sum( tmp_1(ixM^T) * ps(igrid)%dvolume(ixM^T) )

       ! calculate rotational energy
       tmp_2(ixM^T) = ps(igrid)%x(ixM^T, 1) * ( ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,W_vel(2))/lfac(ixM^T) - ps(igrid)%prim(ixM^T,beta(2)) )
       tmp_2(ixM^T) = tmp_2(ixM^T) - &
        ps(igrid)%x(ixM^T, 2) * ( ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,W_vel(1))/lfac(ixM^T) - ps(igrid)%prim(ixM^T,beta(1)) )
       tmp_2(ixM^T) = tmp_2(ixM^T) / (ps(igrid)%x(ixM^T, 1)**2 + ps(igrid)%x(ixM^T, 2)**2) ! here tmp_2 is used to store Omega 
       send_buffer(T_rot) = send_buffer(T_rot) + &
           0.5d0 * sum( tmp_1(ixM^T) * tmp_2(ixM^T) * ps(igrid)%dvolume(ixM^T))

       ! calculate I11_dot
       send_buffer(total_I11dot) = send_buffer(total_I11dot) + &
          sum(  ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T,W_vel(1))*ps(igrid)%x(ixM^T,1) + ps(igrid)%prim(ixM^T,W_vel(1))*ps(igrid)%x(ixM^T,1) ) )

       ! calculate I12_dot
       send_buffer(total_I12dot) = send_buffer(total_I12dot) + &
          sum(  ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T,W_vel(1))*ps(igrid)%x(ixM^T,2) + ps(igrid)%prim(ixM^T,W_vel(2))*ps(igrid)%x(ixM^T,1) ) )

       ! calculate I13_dot
       send_buffer(total_I13dot) = send_buffer(total_I13dot) + &
          sum(  ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T,W_vel(1))*ps(igrid)%x(ixM^T,3) + ps(igrid)%prim(ixM^T,W_vel(3))*ps(igrid)%x(ixM^T,1) ) )

       ! calculate I22_dot
       send_buffer(total_I22dot) = send_buffer(total_I22dot) + &
          sum(  ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T,W_vel(2))*ps(igrid)%x(ixM^T,2) + ps(igrid)%prim(ixM^T,W_vel(2))*ps(igrid)%x(ixM^T,2) ) )

       ! calculate I23_dot
       send_buffer(total_I23dot) = send_buffer(total_I23dot) + &
          sum(  ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T,W_vel(2))*ps(igrid)%x(ixM^T,3) + ps(igrid)%prim(ixM^T,W_vel(3))*ps(igrid)%x(ixM^T,2) ) )

       ! calculate I33_dot
       send_buffer(total_I33dot) = send_buffer(total_I33dot) + &
          sum(  ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T,W_vel(3))*ps(igrid)%x(ixM^T,3) + ps(igrid)%prim(ixM^T,W_vel(3))*ps(igrid)%x(ixM^T,3) ) )
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
