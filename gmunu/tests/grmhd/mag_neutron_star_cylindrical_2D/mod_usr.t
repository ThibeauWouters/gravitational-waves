module mod_usr
  use mod_physics
  use mod_grmhd
  use mod_eos_idealgas
  use mod_XNS
  use mod_multigrid_coupling
  use mod_cfc

  implicit none

  ! aux output variables
  integer, protected  :: divB_ = -1
   
  logical, parameter  ::  init_b_from_vector_pot = .True.
  logical, parameter  ::  initialize_mag_field = .True.
  logical, parameter  ::  initialize_metric = .True.

  !logical, parameter  ::  init_b_from_vector_pot = .False.
  !logical, parameter  ::  initialize_mag_field = .False.
  !logical, parameter  ::  initialize_metric = .False.

  ! if enforcing z_symm during initialisation
  logical, parameter  ::  z_symm = .True.
  logical, parameter  ::  read_id_after_refine = .True.

  ! if perturb
  !logical, parameter  ::  perturbation = .True.
  logical, parameter  ::  perturbation = .False.

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos
    use mod_cfc

    character*128  ::  profile_path = '/users/chi-kit.cheong/simulations'

    profile_path = trim(profile_path) // '/initial_profile/BU0'

    usr_init_one_grid => ns_init_one_grid
    usr_improve_initial_condition => ns_improve_initial_condition
    if (init_b_from_vector_pot) usr_init_vector_potential => init_vec_pot_usr
    usr_process_adv_grid => update_aux_vars
    usr_refine_grid => my_refine
    usr_print_log => printlog

    call set_coordinate_system("cylindrical_2.5D")
    call grmhd_activate()
    ! use atmosphere
    call eos_atmo_activate()
    call eos_idealgas_activate()

    divB_ = var_set_auxvar('divB')

    call initialize_gmunu()

    call mod_XNS_read_profile(profile_path)
    if ( prof%Nth == 1) then 
       stop "! the profile is 1D."
    else
       if ( ndir < 2 ) stop "! Dimension = 1 but the profile is 2D."
    end if

    ! setup atmosphere density from the profile
    call eos_initialize_atmo(maxval(prof%rho))

    ! to use metric solver, we need to activate multigrid
    call cfc_solver_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine ns_init_one_grid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_eos
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim)

    double precision                :: v(ixI^S,1:ndir)
    double precision                :: vphi(ixI^S), betaphi(ixI^S)
    double precision                :: Br(ixI^S), Btheta(ixI^S), Bphi(ixI^S)
    double precision                :: lfac(ixI^S), gamma(ixI^S,1:3,1:3)
    double precision                :: rad, theta
    double precision                :: sin_phi, cos_phi
    double precision                :: sin_theta, cos_theta
    double precision                :: psi_tmp
    double precision                :: xi(ixI^S,1:ndim)
    integer                         :: ix^D, idir, ixC^L

    {^IFONED call mpistop("This is a multi-D HD problem") }

    w = 0.0d0
    v = 0.0d0

    {do ix^D = ixO^LIM^D \} 
       rad = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
       theta = dacos( x(ix^D, 2) / rad )
       if ( z_symm .and. theta > 0.5d0 * dpi ) theta = dpi - theta
       call mod_XNS_map_2D(w(ix^D,rho_),rad, theta,prof%rho(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,press_),rad, theta,prof%press(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,alp_),rad, theta,prof%alp(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,psi_),rad, theta,prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,beta(3)),rad, theta,prof%beta3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(v(ix^D, 3),rad, theta,prof%v3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D, Bvec(3)),rad, theta,prof%B3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
    {end do^D&\} 

    ! initialise B1 and B2 fields
    if (.not. init_b_from_vector_pot) then
       {do ix^D = ixO^LIM^D \} 
          rad = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
          theta = dacos( x(ix^D, 2) / rad )
          cos_theta = x(ix^D, 2) / rad 
          sin_theta = x(ix^D, 1) / rad 

          ! note: no z symm for pol B field
          call mod_XNS_map_2D(Br(ix^D),rad, theta,prof%B1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
          call mod_XNS_map_2D(Btheta(ix^D),rad, theta,prof%B2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
   
          ! rewrite them in orthnormal form
          !Br(ix^D) = Br(ix^D)
          Btheta(ix^D) = Btheta(ix^D) * rad 
   
          w(ix^D, Bvec(1)) = Br(ix^D) * sin_theta + Btheta(ix^D) * cos_theta
          w(ix^D, Bvec(2)) = Br(ix^D) * cos_theta - Btheta(ix^D) * sin_theta

       {end do^D&\} 
       ! if stagger grid CT is used
       if (stagger_grid) then
          block%conss = 0.0d0
          do idir=1,ndim
            ixCmin^D=ixImin^D;
            ixCmax^D=ixImax^D-kr(idir,^D);
            xi = x
            xi(ixI^S,idir)=xi(ixI^S,idir)+0.5d0*block%dx(ixI^S,idir)
            ! cell-face B_idir
            {do ix^D = ixC^LIM^D \} 
              rad = dsqrt( xi(ix^D,1)**2 + xi(ix^D,2)**2 )
              theta = dacos( xi(ix^D, 2) / rad )
              cos_theta = x(ix^D, 2) / rad 
              sin_theta = x(ix^D, 1) / rad 
              call mod_XNS_map_2D(psi_tmp, rad, theta, prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
              call mod_XNS_map_2D(Br(ix^D),rad, theta, prof%B1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
              call mod_XNS_map_2D(Btheta(ix^D),rad, theta, prof%B2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
              Btheta(ix^D) = Btheta(ix^D) * rad
              select case (idir)
              case (1)
                 block%conss(ix^D, idir) = Br(ix^D) * sin_theta + Btheta(ix^D) * cos_theta
              case (2)
                 block%conss(ix^D, idir) = Br(ix^D) * cos_theta - Btheta(ix^D) * sin_theta
              end select
              ! conformally rescale
              block%conss(ix^D, idir) = psi_tmp**6 * block%conss(ix^D, idir) 
            {end do^D&\} 
          end do
          ! although B field at cell centre are filled, we still update it based on B field at cell face
          call phys_face_to_center(ixO^L, block) 
       end if

    else
       call grmhd_b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block)
    end if

    ! fill atmosphere
    where (w(ixO^S,rho_) <= small_rho_thr )
       w(ixO^S,rho_)   = small_rho
       w(ixO^S,eps_)   = small_eps
       w(ixO^S,press_)   = small_press
       v(ixO^S,1) = 0.0d0
       v(ixO^S,2) = 0.0d0
       v(ixO^S,3) = 0.0d0
    end where

    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_)**4 
    end do

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

    ! reset eps
    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
    ! update rest of the prim
    call phys_update_eos(ixI^L, ixO^L, w)
  end subroutine ns_init_one_grid
  
  !> deallocate profile varibles to save memories
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
          ! Geometry subroutines expect this to be set
          block => ps(igrid)
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
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

    if (initialize_metric ) call cfc_metric_init()
    if ( initialize_mag_field ) &
           call phys_initial_clean_divb(it,global_time)

    ! update aux variables
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       call update_aux_vars(igrid,ps(igrid)%level,ixG^LL,ixM^LL,global_time, ps(igrid)) 
    end do

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
    double precision                :: rad, theta, theta_tmp
    double precision                :: sin_theta, cos_theta
    double precision                :: pre_factor
    double precision                :: r_star, delta(0:2) !for perturbation
    double precision                :: dvr(ixI^S), dvtheta(ixI^S)

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

    delta(0) = 1.0d-2
    delta(1) = 0.0d-2
    delta(2) = 0.0d-2

    {do ix^D = ixO^LIM^D \}           
       if ( w(ix^D, rho_) <= small_rho_thr ) cycle
       rad   = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
       theta = dacos( x(ix^D, 2) / rad )
       cos_theta = x(ix^D, 2) / rad 
       sin_theta = x(ix^D, 1) / rad 

       theta_tmp = theta
       if ( z_symm .and. theta_tmp > 0.5d0 * dpi ) theta_tmp = dpi - theta_tmp
       call mod_XNS_map_1D(r_star,theta_tmp,prof%r_star(:),prof%theta,prof%Nth) ! firstly, get r_star(theta)

       if ( rad < r_star ) then
          pre_factor = dsin( dpi * rad / r_star )
       else
          pre_factor = 0.0d0
       end if

       ! for l=0 mode, here we perturbe v^r         
       dvr(ix^D) = delta(0) / w(ix^D, psi_)**4 * pre_factor

       ! for l=2 mode, here we perturbe v^theta       
       dvtheta(ix^D) = delta(1) / (rad*w(ix^D, psi_)**2)**2 * pre_factor &
                       * sin_theta * cos_theta
   
       ! for l=4 mode, here we perturbe v^x(ix^D,theta_)
       dvtheta(ix^D) = dvtheta(ix^D) &
                     + delta(2) / (rad*w(ix^D, psi_)**2)**2 * pre_factor &
                       * sin_theta * cos_theta &
                       * (3.0d0 - 7.0d0 * cos_theta**2 )

       dvtheta(ix^D) = dvtheta(ix^D) * rad 

       v(ix^D,1) = v(ix^D,1) + dvr(ix^D) * sin_theta + dvtheta(ix^D) * cos_theta
       v(ix^D,2) = v(ix^D,2) + dvr(ix^D) * cos_theta - dvtheta(ix^D) * sin_theta

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

  subroutine init_vec_pot_usr(ixI^L, ixO^L, x, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixO^L,idir
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)
    double precision                   :: rad, theta
    integer                            :: ix^D
    A = 0.d0
    if (idir/=3) return ! only Aphi exist
    {do ix^D = ixO^LIM^D \} 
       if ( x(ix^D,1) < tiny(0.0d0) ) cycle
       rad = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
       theta = dacos( x(ix^D, 2) / rad )
       call mod_XNS_map_2D(A(ix^D),rad, theta,prof%A3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
    {end do^D&\} 
  end subroutine init_vec_pot_usr

  subroutine update_aux_vars(igrid,level,ixI^L,ixO^L,qt,s)
    use mod_global_parameters
    use mod_grmhd_phys_parameters
    use mod_eos
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s
    double precision                :: divb(ixI^S)
    integer                         :: ix^D
    !> calculate div B after the advance 
    call grmhd_get_divb(ixI^L, ixO^L, s, divb(ixI^S))
    s%prim(ixO^S, divB_) = dabs( divb(ixO^S) )
    ! fixme: do we really need this?
    !> when polytrope eos is used, P is function of rho only, without eps.
    !> we need to reset eps based on updated P and rho every steps
    if (eos_type == polytrope) then
       {do ix^D = ixI^LIM^D \}
          ! reset eps
          call eos_get_eps_one_grid(ps(igrid)%prim(ix^D, press_),ps(igrid)%prim(ix^D, rho_),ps(igrid)%prim(ix^D, eps_))
          ps(igrid)%prim(ix^D, eps_) = max( ps(igrid)%prim(ix^D, eps_), small_eps )
       {end do^D&\}
       ! update rest of the prim
       call phys_update_eos(ixI^L, ixO^L, ps(igrid)%prim(ixI^S,1:nprim))
    end if
  end subroutine update_aux_vars

  subroutine my_refine(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nprim), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    double precision             :: rho_max
    double precision, parameter  :: rho_cut = 1.0d-5
    double precision             :: ratio(ixI^S)
    integer                      :: i_level, to_level
    double precision             :: phi_grav ! relativistic gravitational potential :=  1 - alp
    double precision             :: phi_grav_cut
    double precision, parameter  :: phi_grav_max = 0.2d0

    refine = 0
    coarsen = 0

    ! check if the density is almost/already
    ! above rho_cut
    rho_max = maxval( w(ixO^S, rho_) )
    if (rho_max >= rho_cut) then
       ! always at the finest level
       refine = 1
       coarsen = -1
       return
    end if

    ! check if this is the most outer block
    if ( (maxval(x(ixO^S,1))+dx(1,level)) >= xprobmax1 .or. &
         (maxval(x(ixO^S,2))+dx(2,level)) >= xprobmax2 .or. &
         ((minval(x(ixO^S,2))-dx(2,level)) <= xprobmin2 .and. xprobmin2 < 0.0d0 ) &
         ) then
       ! always at the lowest level
       refine = -1
       coarsen = 1
       return
    end if

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

    ! if the current level is higher than the maximum allowed level, 
    ! we strictly limit the level.
    ! However, we do not restrict if the level is too low.
    ! The goal here is just put the upper bound of the level.

    ! if refine_criterion == 0, purely depents on user define level,
    ! then do the following
    else if ( refine_criterion==0 .and. level < to_level ) then
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
    double precision     :: divB(ixG^T)

    double precision     :: rho_max, rho_max_local
    double precision     :: alp_min, alp_min_local

    integer, parameter   :: total_volume = 1
    integer, parameter   :: total_mass = 2
    integer, parameter   :: J_rot = 3
    integer, parameter   :: T_rot = 4
    integer, parameter   :: total_EB_tor = 5
    integer, parameter   :: total_EB_pol = 6
    integer, parameter   :: total_divBdV = 7
    integer, parameter   :: total_I11dot = 8
    integer, parameter   :: total_I12dot = 9
    integer, parameter   :: total_I22dot = 10
    integer, parameter   :: n_var = 10

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

       ! calculate angular momentum 
       send_buffer(J_rot) = send_buffer(J_rot) + &
          sum( ps(igrid)%cons(ixM^T,mom(3)) * ps(igrid)%dvolume(ixM^T) )

       ! calculate kinetic energy (related to v phi only!) 
       call phys_get_lfac2(ixG^LL, ixM^LL, ps(igrid), lfac(ixG^T))
       lfac(ixM^T) = dsqrt( lfac(ixM^T) )
       send_buffer(T_rot) = send_buffer(T_rot) + &
          0.5d0 * sum( ps(igrid)%prim(ixM^T,W_vel(3))/lfac(ixM^T) * ps(igrid)%cons(ixM^T,mom(3)) * ps(igrid)%dvolume(ixM^T) )

       ! get divB
       call grmhd_get_divb(ixG^LL, ixM^LL, ps(igrid), divb(ixG^T))

       ! calculate total divB dV
       send_buffer(total_divBdV) = send_buffer(total_divBdV) + &
          sum(dabs(divB(ixM^T)) * dV(ixM^T))

       ! calculate total EB_tor
       send_buffer(total_EB_tor) = send_buffer(total_EB_tor) + &
          0.5d0 * sum( ps(igrid)%x(ixM^T,1)**2 * ps(igrid)%prim(ixM^T,psi_)**4 &
                * ps(igrid)%prim(ixM^T, Bvec(3))**2 * dV(ixM^T) )

       ! calculate total EB_pol
       send_buffer(total_EB_pol) = send_buffer(total_EB_pol) + &
          0.5d0 * sum( ps(igrid)%prim(ixM^T,psi_)**4 &
                * (ps(igrid)%prim(ixM^T, Bvec(1))**2 + ps(igrid)%prim(ixM^T, Bvec(2))**2)* dV(ixM^T) )

       ! calculate I11_dot
       send_buffer(total_I11dot) = send_buffer(total_I11dot) + &
          sum(  ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T,W_vel(1))*ps(igrid)%x(ixM^T,1) + ps(igrid)%prim(ixM^T,W_vel(1))*ps(igrid)%x(ixM^T,1) ) )

       ! calculate I12_dot
       send_buffer(total_I12dot) = send_buffer(total_I12dot) + &
          sum(  ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T,W_vel(1))*ps(igrid)%x(ixM^T,2) + ps(igrid)%prim(ixM^T,W_vel(2))*ps(igrid)%x(ixM^T,1) ) )

       ! calculate I22_dot
       send_buffer(total_I22dot) = send_buffer(total_I22dot) + &
          sum(  ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,rho_) * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T,W_vel(2))*ps(igrid)%x(ixM^T,2) + ps(igrid)%prim(ixM^T,W_vel(2))*ps(igrid)%x(ixM^T,2) ) )

    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    ! NOTE: 2pi factor is not included
    !recive_buffer = recive_buffer * 2.0d0 * dpi

    ! Volume average of divB
    recive_buffer(total_divBdV) = (recive_buffer(total_divBdV) / recive_buffer(total_volume))

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
          line = trim(line) //  " J_rot T_rot EB_tor EB_pol" 
          line = trim(line) // " avg_divB"
          line = trim(line) // " I11_dot I12_dot I22_dot"

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
