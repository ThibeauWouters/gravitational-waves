module mod_usr
  use mod_physics
  use mod_grmhd
  use mod_eos_idealgas
  use mod_XNS
  use mod_multigrid_coupling
  use mod_cfc

  implicit none
  private

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
  logical, parameter  ::  exclude_atmo = .true.  ! exclude atmosphere when post process

  public              :: usr_init

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos
    use mod_cfc

    character*128  ::  profile_path = '/users/chi-kit.cheong/simulations'
    !character*128  ::  profile_path = '/scratch/s1/TjonnieLi/patrick/simulations'
    profile_path = trim(profile_path) // '/initial_profile/BU0'

    usr_init_one_grid => ns_init_one_grid
    usr_improve_initial_condition => ns_improve_initial_condition
    if (init_b_from_vector_pot) usr_init_vector_potential => init_vec_pot_usr
    usr_process_adv_grid => update_aux_vars
    usr_refine_grid => my_refine
    usr_print_log => printlog

    call set_coordinate_system("Cartesian")
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

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim)

    double precision                :: v(ixI^S,1:ndir)
    double precision                :: vphi(ixI^S), betaphi(ixI^S)
    double precision                :: Br(ixI^S), Btheta(ixI^S), Bphi(ixI^S)
    double precision                :: lfac(ixI^S), gamma(ixI^S,1:3,1:3)
    double precision                :: xi(ixI^S,1:^ND)
    integer                         :: ix^D, idir, ixC^L
    double precision                :: rad, rad_xy, theta
    double precision                :: sin_theta, cos_theta
    double precision                :: sin_phi, cos_phi
    double precision                :: psi_tmp

    {^IFONED call mpistop("This is a multi-D HD problem") }

    w = 0.0d0
    v = 0.0d0

    {do ix^D = ixO^LIM^D \} 
       rad = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
       theta = dacos( x(ix^D, 3) / rad )
       if ( z_symm .and. theta > 0.5d0 * dpi ) theta = dpi - theta
       call mod_XNS_map_2D(w(ix^D,rho_),rad, theta,prof%rho(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,press_),rad, theta,prof%press(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,alp_),rad, theta,prof%alp(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,psi_),rad, theta,prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(betaphi(ix^D),rad, theta,prof%beta3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(vphi(ix^D),rad, theta,prof%v3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)

       v(ix^D, 1) = - x(ix^D, 2) * vphi(ix^D)
       v(ix^D, 2) =   x(ix^D, 1) * vphi(ix^D)
       w(ix^D, beta(1)) = - x(ix^D, 2) * betaphi(ix^D)
       w(ix^D, beta(2)) =   x(ix^D, 1) * betaphi(ix^D)
    {end do^D&\} 

    ! initialise B fields
    if (.not. init_b_from_vector_pot) then
       {do ix^D = ixO^LIM^D \} 
          rad = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
          rad_xy = dsqrt( sum(x(ix^D,1:2)**2) ) 
          theta = dacos( x(ix^D, 3) / rad )
          call mod_XNS_map_2D(Br(ix^D),rad, theta,prof%B1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
          call mod_XNS_map_2D(Btheta(ix^D),rad, theta,prof%B2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
          call mod_XNS_map_2D(Bphi(ix^D),rad, theta,prof%B3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
   
          Br(ix^D) = Br(ix^D)
          Btheta(ix^D) = Btheta(ix^D) * rad
          Bphi(ix^D) = Bphi(ix^D) * rad_xy
   
          sin_theta = rad_xy / rad
          cos_theta = x(ix^D, 3) / rad
          sin_phi = x(ix^D, 2) / rad_xy
          cos_phi = x(ix^D, 1) / rad_xy
   
          w(ix^D, Bvec(1)) = Br(ix^D) * sin_theta * cos_phi + Btheta(ix^D) * cos_theta * cos_phi - Bphi(ix^D) * sin_phi
          w(ix^D, Bvec(2)) = Br(ix^D) * sin_theta * sin_phi + Btheta(ix^D) * cos_theta * sin_phi + Bphi(ix^D) * cos_phi
          w(ix^D, Bvec(3)) = Br(ix^D) * cos_theta           - Btheta(ix^D) * sin_theta 
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
               rad = dsqrt( sum(xi(ix^D,1:ndim)**2) ) 
               rad_xy = dsqrt( sum(xi(ix^D,1:2)**2) ) 
               theta = dacos( xi(ix^D, 3) / rad )
               call mod_XNS_map_2D(psi_tmp,rad,theta,prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
               call mod_XNS_map_2D(Br(ix^D),rad, theta,prof%B1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
               call mod_XNS_map_2D(Btheta(ix^D),rad, theta,prof%B2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
               call mod_XNS_map_2D(Bphi(ix^D),rad, theta,prof%B3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
   
               Br(ix^D) = Br(ix^D)
               Btheta(ix^D) = Btheta(ix^D) * rad
               Bphi(ix^D) = Bphi(ix^D) * rad_xy

               sin_theta = rad_xy / rad
               cos_theta = xi(ix^D, 3) / rad
               sin_phi = xi(ix^D, 2) / rad_xy
               cos_phi = xi(ix^D, 1) / rad_xy
   
               select case (idir)
               case (1)
                  block%conss(ix^D,idir) = Br(ix^D) * sin_theta * cos_phi + Btheta(ix^D) * cos_theta * cos_phi - Bphi(ix^D) * sin_phi
               case (2)
                  block%conss(ix^D,idir) = Br(ix^D) * sin_theta * sin_phi + Btheta(ix^D) * cos_theta * sin_phi + Bphi(ix^D) * cos_phi
               case (3)
                  block%conss(ix^D,idir) = Br(ix^D) * cos_theta           - Btheta(ix^D) * sin_theta 
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

    w(ixO^S, divB_) = 0.d0

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

    ! update aux and conserved variables
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       call phys_to_conserved(ixG^LL,ixM^LL,ps(igrid)%cons(ixG^T,1:ncons),ps(igrid)%prim(ixG^T,1:nprim),ps(igrid)%x(ixG^T,1:ndim))
       call update_aux_vars(igrid,ps(igrid)%level,ixG^LL,ixM^LL,global_time, ps(igrid)) 
    end do

    ! deallocate profile varibles to save memories
    call mod_XNS_deallocate_var()

  end subroutine ns_improve_initial_condition

  subroutine init_vec_pot_usr(ixI^L, ixO^L, x, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixO^L,idir
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)
    double precision                   :: rad, rad_xy_sqr, theta
    double precision                   :: p_thr, press
    integer                            :: ix^D
    double precision, parameter        :: Ab = 2.823d0 !4.0d-2 
    A = 0.d0
    if (idir==3) return ! only A_phi exist
    ! press_thr
    p_thr = 4.0d-2 * maxval(prof%press)
    {do ix^D = ixO^LIM^D \} 
       rad = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
       !rad_xy_sqr = sum(x(ix^D,1:2)**2)
       theta = dacos( x(ix^D, 3) / rad )
       call mod_XNS_map_2D(press,rad,theta,prof%press(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       A(ix^D) = Ab * max( (press - p_thr), 0.0d0 ) ! Aphi, without 1/rad_xy_sqr factor
       select case ( idir )
       case (1)
          A(ix^D) = - A(ix^D) * x(ix^D, 2)
       case (2)
          A(ix^D) =   A(ix^D) * x(ix^D, 1)
       end select
    {end do^D&\} 
  end subroutine init_vec_pot_usr

  ! vector potential from XNS
  subroutine init_vec_pot_xns(ixI^L, ixO^L, x, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixO^L,idir
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)
    double precision                   :: rad, rad_xy_sqr, theta
    integer                            :: ix^D
    A = 0.d0
    if (idir==3) return ! only A_phi exist
    {do ix^D = ixO^LIM^D \} 
       rad_xy_sqr = sum(x(ix^D,1:2)**2)
       if ( rad_xy_sqr < tiny(0.0d0) ) cycle
       rad = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
       theta = dacos( x(ix^D, 3) / rad )
       call mod_XNS_map_2D(A(ix^D),rad,theta,prof%A3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       select case ( idir )
       case (1)
          A(ix^D) = - A(ix^D) * x(ix^D, 2) / rad_xy_sqr
       case (2)
          A(ix^D) =   A(ix^D) * x(ix^D, 1) / rad_xy_sqr
       end select
    {end do^D&\} 
  end subroutine init_vec_pot_xns

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
         (minval(x(ixO^S,1))-dx(1,level)) <= xprobmin1 .or. &
         (maxval(x(ixO^S,2))+dx(2,level)) >= xprobmax2 .or. &
         (minval(x(ixO^S,2))-dx(2,level)) <= xprobmin2 .or. &
         (maxval(x(ixO^S,3))+dx(3,level)) >= xprobmax3 .or. &
         ((minval(x(ixO^S,3))-dx(3,level)) <= xprobmin3 .and. abs(xprobmin3) <= smalldouble) &
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

  subroutine update_aux_vars(igrid,level,ixI^L,ixO^L,qt,s)
    use mod_global_parameters
    use mod_physics
    use mod_eos
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s
    double precision                :: divb(ixI^S)
    integer                         :: ix^D
    !> calculate div B after the advance 
    call phys_get_divb(ixI^L, ixO^L, s, divb(ixI^S))
    s%prim(ixO^S, divB_) = dabs( divb(ixO^S) )

    !> when polytrope eos is used, P is function of rho only, without eps.
    !> we need to reset eps based on updated P and rho every steps
    if (eos_type == polytrope) then
       {do ix^D = ixI^LIM^D \}
          ! reset eps
          call eos_get_eps_one_grid(ps(igrid)%prim(ix^D, press_),ps(igrid)%prim(ix^D, rho_),ps(igrid)%prim(ix^D, eps_))
          ps(igrid)%prim(ix^D, eps_) = max( ps(igrid)%prim(ix^D, eps_), small_eps )
       {end do^D&\}
   
       ! update rest of the prim
       !call phys_update_eos(ixI^L, ixO^L, ps(igrid)%prim(ixI^S,1:nprim))
    end if
  end subroutine update_aux_vars

  subroutine printlog()
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_eos
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

    integer              :: igrid, iigrid, idir, jdir, kdir
    double precision     :: h(ixG^T)
    double precision     :: lfac(ixG^T), dV(ixG^T), dV_hat(ixG^T)
    double precision     :: phi(ixG^T)
    double precision     :: tmp_1(ixG^T), tmp_2(ixG^T)

    double precision     :: rho_max, rho_max_local
    double precision     :: alp_min, alp_min_local
    double precision     :: B_max, B_max_local

    integer              :: total_volume = -1
    integer              :: rest_mass = -1
    integer              :: proper_mass = -1
    integer              :: J_rot = -1
    integer              :: T_rot = -1
    integer              :: E_kin = -1
    integer              :: E_tot = -1
    integer              :: total_divBdV = -1
    integer              :: total_EB_tor = -1
    integer              :: total_EB = -1
    integer              :: total_EEM = -1
    integer              :: total_B = -1
    integer              :: Iij_dot(1:6) ! I11, I12, I13, I22, I23, I33
    integer              :: cm_rho(1:4,1:2) 
    integer              :: cm_Bz(1:4,1:2) 
    integer              :: n_var = -1

    ! initialise array
    Iij_dot = -1
    cm_rho = -1

    ! initialise local variables
    B_max_local = 0.0d0
    rho_max_local = 0.0d0
    alp_min_local = 1.0d0

    ! initialize variables, output will follow this order
    n_var = 0
    total_volume = add_var()
    rest_mass    = add_var()
    proper_mass  = add_var()
    total_divBdV = add_var()
    total_B      = add_var()
    total_EB_tor = add_var()
    total_EB     = add_var()
    total_EEM    = add_var()
    J_rot        = add_var()
    T_rot        = add_var()
    E_kin        = add_var()
    E_tot        = add_var()
    do idir = 1, 4; do jdir = 1, 2
       cm_rho(idir,jdir) = add_var()
    end do; end do
    do idir = 1, 4; do jdir = 1, 2
       cm_Bz(idir,jdir) = add_var()
    end do; end do

    ! Iij_dot
    do idir = 1, 6
       Iij_dot(idir) = add_var()
    end do

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       ! find max rho
       rho_max_local = max(rho_max_local, maxval(ps(igrid)%prim(ixM^T,rho_)) )
       ! find min alpha
       alp_min_local = min(alp_min_local, minval(ps(igrid)%prim(ixM^T,alp_)) )

       ! make sure all conserved variables are updated, maybe no needed
       !call phys_to_conserved(ixG^LL,ixM^LL,ps(igrid)%cons(ixG^T,1:ncons),ps(igrid)%prim(ixG^T,1:nprim),ps(igrid)%x(ixG^T,1:ndim))

       call phys_get_lfac2(ixG^LL, ixM^LL, ps(igrid), lfac(ixG^T))
       lfac(ixM^T) = dsqrt( lfac(ixM^T) )

       ! calculate angle phi for coordinate transformation
       phi(ixM^T) = datan2( ps(igrid)%x(ixM^T, 2) , ps(igrid)%x(ixM^T, 1) )

       dV_hat(ixM^T) = ps(igrid)%dvolume(ixM^T)
       ! exclude atmosphere part
       if (exclude_atmo) then
          where ( ps(igrid)%prim(ixM^T, rho_) <= small_rho_thr ) 
             dV_hat(ixM^T) = 0.0d0
          end where
       end if

       ! calculate cell volume
       dV(ixM^T) = ps(igrid)%prim(ixM^T,psi_)**6 * dV_hat(ixM^T)
       ! calculate total volume
       send_buffer(total_volume) = send_buffer(total_volume) + sum( dV(ixM^T) )

       ! local rest mass
       tmp_2(ixM^T) = ps(igrid)%cons(ixM^T,D_) * dV_hat(ixM^T) 
       ! calculate total rest mass
       send_buffer(rest_mass) = send_buffer(rest_mass) + sum( tmp_2(ixM^T) )
       ! calculate total proper mass
       send_buffer(proper_mass) = send_buffer(proper_mass) + sum( tmp_2(ixM^T) * (1.0d0 + ps(igrid)%prim(ixM^T, eps_)) )

       ! calculate m=idir azimuthal decomposition 
       do idir = 1, 4; do jdir = 1, 2
          if (jdir == 1) then ! real part
             tmp_1(ixM^T) = dcos( dble(idir) * phi(ixM^T) )
          else ! img part
             tmp_1(ixM^T) = dsin( dble(idir) * phi(ixM^T) )
          end if
          send_buffer(cm_rho(idir,jdir)) = send_buffer(cm_rho(idir,jdir)) &
              + sum( tmp_2(ixM^T) * tmp_1(ixM^T) )
          send_buffer(cm_Bz(idir,jdir)) = send_buffer(cm_Bz(idir,jdir)) &
              + sum( ps(igrid)%prim(ixM^T, Bvec(3)) * tmp_1(ixM^T) * dV(ixM^T) )
       end do; end do

       ! calculate Iij_dot
       kdir = 1
       do idir = 1, 3
          do jdir = idir, 3
             send_buffer(Iij_dot(kdir)) = send_buffer(Iij_dot(kdir)) + &
                sum(  tmp_2(ixM^T) / lfac(ixM^T) &
                      * (  ps(igrid)%prim(ixM^T,W_vel(idir))*ps(igrid)%x(ixM^T,jdir) &
                         + ps(igrid)%prim(ixM^T,W_vel(jdir))*ps(igrid)%x(ixM^T,idir) ) )
             kdir = kdir + 1
          end do
       end do

       ! calculate angular momentum 
       h(ixM^T) = 1.0d0 + ps(igrid)%prim(ixM^T,eps_) + ps(igrid)%prim(ixM^T,press_) / ps(igrid)%prim(ixM^T,rho_)
       tmp_1(ixM^T) = - ps(igrid)%cons(ixM^T,D_) * h(ixM^T) * ps(igrid)%x(ixM^T, 2) * &
              ( ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,W_vel(1)) - ps(igrid)%prim(ixM^T,beta(1))*lfac(ixM^T) )
       tmp_1(ixM^T) = tmp_1(ixM^T) + ps(igrid)%cons(ixM^T,D_) * h(ixM^T) * ps(igrid)%x(ixM^T, 1) * &
              ( ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,W_vel(2)) - ps(igrid)%prim(ixM^T,beta(2))*lfac(ixM^T) )
       tmp_1(ixM^T) = ps(igrid)%prim(ixM^T,psi_)**4 * tmp_1(ixM^T)
       send_buffer(J_rot) = send_buffer(J_rot) + sum( tmp_1(ixM^T) * dV_hat(ixM^T) )

       ! calculate rotational energy
       tmp_2(ixM^T) = ps(igrid)%x(ixM^T, 1) * ( ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,W_vel(2))/lfac(ixM^T) - ps(igrid)%prim(ixM^T,beta(2)) )
       tmp_2(ixM^T) = tmp_2(ixM^T) - &
           ps(igrid)%x(ixM^T, 2) * ( ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,W_vel(1))/lfac(ixM^T) - ps(igrid)%prim(ixM^T,beta(1)) )
       tmp_2(ixM^T) = tmp_2(ixM^T) / (ps(igrid)%x(ixM^T, 1)**2 + ps(igrid)%x(ixM^T, 2)**2) ! here tmp_2 is used to store Omega 
       send_buffer(T_rot) = send_buffer(T_rot) + &
           0.5d0 * sum( tmp_1(ixM^T) * tmp_2(ixM^T) * dV_hat(ixM^T))

       ! calculate kin energy
       ! u_i v^i
       tmp_1(ixM^T) = ps(igrid)%prim(ixM^T,psi_)**4 / lfac(ixM^T) &
                * ( ps(igrid)%prim(ixM^T, W_vel(1))**2 &
                  + ps(igrid)%prim(ixM^T, W_vel(2))**2 &
                  + ps(igrid)%prim(ixM^T, W_vel(3))**2 )
       tmp_1(ixM^T) = ps(igrid)%cons(ixM^T,D_) * h(ixM^T) * tmp_1(ixM^T)
       send_buffer(E_kin) = send_buffer(E_kin) + 0.5d0 * sum( tmp_1(ixM^T) * dV_hat(ixM^T))

       tmp_1(ixM^T) = ps(igrid)%cons(ixM^T,tau_) + ps(igrid)%cons(ixM^T,D_)
       send_buffer(E_tot) = send_buffer(E_tot) + sum( tmp_1(ixM^T) * dV_hat(ixM^T))

       ! calculate total divB dV
       send_buffer(total_divBdV) = send_buffer(total_divBdV) + &
          sum(dabs(ps(igrid)%prim(ixM^T, divB_)) * dV(ixM^T))

       ! calculate local B2, store in tmp1
       tmp_1(ixM^T) = ps(igrid)%prim(ixM^T,psi_)**4 &
                * ( ps(igrid)%prim(ixM^T, Bvec(1))**2 &
                  + ps(igrid)%prim(ixM^T, Bvec(2))**2 &
                  + ps(igrid)%prim(ixM^T, Bvec(3))**2 )

       ! calculate total magnetic energy
       send_buffer(total_EB) = send_buffer(total_EB) &
          + 0.5d0 * sum( ps(igrid)%prim(ixM^T,alp_) * tmp_1(ixM^T) * dV(ixM^T) )

       ! calculate EEM
       send_buffer(total_EEM) = send_buffer(total_EEM) &
          + 0.5d0 * sum( ps(igrid)%prim(ixM^T,alp_) * tmp_1(ixM^T)*(2.0d0 - 1.0d0/lfac(ixM^T)**2) * dV(ixM^T) )
       ! B dot v
       tmp_2(ixM^T) = ps(igrid)%prim(ixM^T,psi_)**4 / lfac(ixM^T) &
                * ( ps(igrid)%prim(ixM^T, Bvec(1))*ps(igrid)%prim(ixM^T, W_vel(1)) &
                  + ps(igrid)%prim(ixM^T, Bvec(2))*ps(igrid)%prim(ixM^T, W_vel(2)) &
                  + ps(igrid)%prim(ixM^T, Bvec(3))*ps(igrid)%prim(ixM^T, W_vel(3)) )
       send_buffer(total_EEM) = send_buffer(total_EEM) &
          - 0.5d0 * sum( ps(igrid)%prim(ixM^T,alp_) * tmp_2(ixM^T)**2 * dV(ixM^T) )

       ! calculate local |B|, store in tmp1
       tmp_1(ixM^T) = dsqrt(tmp_1(ixM^T))

       ! find max |B|
       B_max_local = max(B_max_local, maxval(tmp_1(ixM^T)) )

       ! calculate total |B|
       send_buffer(total_B) = send_buffer(total_B) + sum( tmp_1(ixM^T) * dV(ixM^T) )

       ! calculate total EB_tor, store in tmp_1
       tmp_1(ixM^T) = ( - ps(igrid)%prim(ixM^T, Bvec(1)) * dsin(phi(ixM^T)) &
                        + ps(igrid)%prim(ixM^T, Bvec(2)) * dcos(phi(ixM^T)) )
       send_buffer(total_EB_tor) = send_buffer(total_EB_tor) + &
          0.5d0 * sum( ps(igrid)%prim(ixM^T,psi_)**4 * dV(ixM^T) * tmp_1(ixM^T)**2 )

    end do

    ! gather all the results!
    call MPI_ALLREDUCE(rho_max_local, rho_max, 1, mpi_double_precision, &
          MPI_MAX, icomm, ierrmpi)
    call MPI_ALLREDUCE(alp_min_local, alp_min, 1, mpi_double_precision, &
          MPI_MIN, icomm, ierrmpi)
    call MPI_ALLREDUCE(B_max_local, B_max, 1, mpi_double_precision, &
          MPI_MAX, icomm, ierrmpi)
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    ! NOTE: 2 factor is not included
    !recive_buffer = recive_buffer * 2.0d0

    ! Volume average of divB
    recive_buffer(total_divBdV) = (recive_buffer(total_divBdV) / recive_buffer(total_volume))
    ! Volume average of |B|
    recive_buffer(total_B) = (recive_buffer(total_B) / recive_buffer(total_volume))

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
          line = "it global_time dt rho_max alp_min B_max" 
          line = trim(line) // " Vol_tot M_rest M_proper"
          line = trim(line) // " divB_avg B_avg EB_tor EB_tot E_EM"
          line = trim(line) // " J_rot T_rot E_kin E_tot"
          line = trim(line) // " c1_real c1_img c2_real c2_img"
          line = trim(line) // " c3_real c3_img c4_real c4_img"
          line = trim(line) // " c1_Bz_real c1_Bz_img c2_Bz_real c2_Bz_img"
          line = trim(line) // " c3_Bz_real c3_Bz_img c4_Bz_real c4_Bz_img"
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

       write(fmt_string, '(a,i0,a)') '(', 1, fmt_r // ')'
       write(line(i:), fmt_string) B_max
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', n_var, fmt_r // ')'
       write(line(i:), fmt_string) recive_buffer(1:n_var)
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
