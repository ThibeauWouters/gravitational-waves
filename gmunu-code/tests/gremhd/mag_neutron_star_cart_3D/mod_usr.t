module mod_usr
  use mod_physics
  use mod_gremhd
  use mod_eos_idealgas
  use mod_XNS
  use mod_multigrid_coupling
  use mod_cfc

  implicit none

  ! aux output variables
  integer, protected  :: divB_ = -1
  integer, protected  :: divE_ = -1
   
  logical, parameter  ::  init_b_from_vector_pot = .True.
  logical, parameter  ::  initialize_mag_field = .True.
  logical, parameter  ::  initialize_metric = .True.
  logical, parameter  ::  evolve_metric = .True.

  !logical, parameter  ::  init_b_from_vector_pot = .False.
  !logical, parameter  ::  initialize_mag_field = .False.
  !logical, parameter  ::  initialize_metric = .False.
  !logical, parameter  ::  evolve_metric = .False.

  ! conducivity in the center of the star
  double precision, parameter  ::  sigma_0 = 1.0d6

  logical, parameter  ::  read_id_after_refine = .True.

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
    integer  ::  Nth, Nr   ! number of grid of the profile

    profile_path = trim(profile_path) // '/initial_profile/mag_AU0_1.40_p_ref_large_domain'
    !profile_path = trim(profile_path) // '/initial_profile/mag_ns_pol_weak_large_domain'

    usr_init_one_grid => ns_init_one_grid
    usr_improve_initial_condition => ns_improve_initial_condition
    if (init_b_from_vector_pot) usr_init_vector_potential => init_vec_pot_usr
    usr_get_resistivity => my_get_resistivity
    usr_process_adv_grid => update_aux_vars
    usr_refine_grid => my_refine
    usr_print_log => printlog

    call set_coordinate_system("Cartesian")
    call gremhd_activate()
    ! use atmosphere
    call eos_atmo_activate()
    call eos_idealgas_activate()

    divB_ = var_set_auxvar('divB')
    divE_ = var_set_auxvar('divE')

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
    double precision                :: Er(ixI^S), Etheta(ixI^S), Ephi(ixI^S)
    double precision                :: lfac(ixI^S), gamma(ixI^S,1:3,1:3)
    double precision                :: xi(ixI^S,1:^ND)
    integer                         :: ix^D, idir, ixC^L
    double precision                :: rad, rad_xy, theta
    double precision                :: sin_phi, cos_phi
    double precision                :: sin_theta, cos_theta
    double precision                :: psi_tmp

    {^IFONED call mpistop("This is a multi-D HD problem") }

    w = 0.0d0
    v = 0.0d0

    {do ix^D = ixO^LIM^D \} 
       rad_xy = dsqrt( sum(x(ix^D,1:2)**2) ) 
       rad = dsqrt( sum(x(ix^D,1:ndim)**2) ) 
       theta = acos( x(ix^D, 3) / rad )
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

       call mod_XNS_map_2D(Er(ix^D),rad, theta,prof%E_1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(Etheta(ix^D),rad, theta,prof%E_2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(Ephi(ix^D),rad, theta,prof%E_3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)

       sin_theta = rad_xy / rad
       cos_theta = x(ix^D, 3) / rad
       sin_phi = x(ix^D, 2) / rad_xy
       cos_phi = x(ix^D, 1) / rad_xy
   
       ! fixme: E is E_i
       !Er(ix^D) = Er(ix^D)
       Etheta(ix^D) = Etheta(ix^D) * rad
       Ephi(ix^D) = Ephi(ix^D) * rad_xy

       w(ix^D, Evec(1)) = Er(ix^D) * sin_theta * cos_phi + Etheta(ix^D) * cos_theta * cos_phi - sin_phi * Ephi(ix^D)
       w(ix^D, Evec(2)) = Er(ix^D) * sin_theta * sin_phi + Etheta(ix^D) * cos_theta * sin_phi + cos_phi * Ephi(ix^D)
       w(ix^D, Evec(3)) = Er(ix^D) * cos_theta - Etheta(ix^D) * sin_theta
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
   
          !Br(ix^D) = Br(ix^D)
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
   
               !Br(ix^D) = Br(ix^D)
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
       call gremhd_b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block)
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
  
  !> before the main loop, we improve the initial data here
  subroutine ns_improve_initial_condition()
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_physics
    use mod_eos
    use mod_cfc
    use mod_gremhd_phys_divb_mg

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

    if (evolve_metric .and. initialize_metric ) &
           call cfc_metric_init()
    if ( initialize_mag_field ) &
           call phys_initial_clean_divb(it,global_time)

    ! update aux variables
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       block => ps(igrid) 
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
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
  end subroutine init_vec_pot_usr

  ! get resistivity
  subroutine my_get_resistivity(ixI^L,ixO^L,cons,x,eta)
    use mod_global_parameters
    use mod_eos
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(in) :: cons(ixI^S,1:ncons) 
    double precision, intent(out):: eta(ixI^S) 

    double precision             :: tmp
    double precision             :: sigma(ixI^S)
    integer                      :: ix^D

    integer, parameter           :: conductivity_formular = 2
    double precision, parameter  :: D_tol = 1.0d-2
    double precision, parameter  :: D_rel_over_D_atmo = 1.0d2

    select case (conductivity_formular)
    case (1)
    ! method 1: directly assign values
       eta(ixO^S) = 1.0d0 / sigma_0 
    case (2)
    ! method 2: PRD.88.044020
       sigma(ixI^S) = smalldouble !epsilon(0.0d0)
       sigma(ixO^S) = 1.0d0 - small_D / cons(ixO^S, D_) 
       where ( sigma(ixO^S) <= smalldouble )
          sigma(ixO^S) = smalldouble
       else where
          sigma(ixO^S) = sigma_0 * sigma(ixO^S)**2
       end where
       eta(ixI^S) = 1.0d0 / sigma(ixI^S) 
    case (3)
    ! method 3: PRD.92.084064
       sigma(ixI^S) = epsilon(0.0d0)
       sigma(ixO^S) = dexp( 2.0d0 * D_tol * ( cons(ixO^S, D_) / small_D - D_rel_over_D_atmo ) )
       sigma(ixO^S) = 1.0d0 - 2.0d0 / ( 1.0d0 + sigma(ixO^S) )
       !sigma(ixO^S) = sigma_0 * max( sigma(ixO^S), 0.0d0 )
       sigma(ixO^S) = sigma_0 * sigma(ixO^S)
       where ( sigma(ixI^S) <= 0.0d0 )
          eta(ixI^S) = huge(sigma_0) ! largest value in double
       else where ( sigma(ixI^S) > sigma_0 )
          eta(ixI^S) = 1.0d0 / sigma_0
       else where
          eta(ixI^S) = 1.0d0 / sigma(ixI^S) 
       end where
    end select
  end subroutine my_get_resistivity

  subroutine update_aux_vars(igrid,level,ixI^L,ixO^L,qt,s)
    use mod_global_parameters
    use mod_gremhd_phys_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s
    double precision                :: divb(ixI^S)
    !> calculate div B after the advance 
    call gremhd_get_divb(ixI^L, ixO^L, s, divb(ixI^S))
    s%prim(ixO^S, divB_) = dabs( divb(ixO^S) )
  end subroutine update_aux_vars

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
    double precision     :: h(ixG^T)
    double precision     :: lfac(ixG^T),dV(ixG^T)
    double precision     :: phi(ixG^T)
    double precision     :: tmp_1(ixG^T), tmp_2(ixG^T), tmp_3(ixG^T)
    double precision     :: rsq(ixG^T)
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
    integer              :: total_divBdV = -1
    integer              :: total_EB_tor = -1
    integer              :: total_EB = -1
    integer              :: total_EE = -1
    integer              :: total_I11dot = -1
    integer              :: total_I12dot = -1
    integer              :: total_I13dot = -1
    integer              :: total_I22dot = -1
    integer              :: total_I23dot = -1
    integer              :: total_I33dot = -1
    integer              :: total_B = -1
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
    total_divBdV = add_var()
    total_B = add_var()
    total_EB_tor = add_var()
    total_EB = add_var()
    total_EE = add_var()
    J_rot = add_var()
    T_rot = add_var()
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

       ! make sure all conserved variables are updated, maybe no needed
       call phys_to_conserved(ixG^LL,ixM^LL,ps(igrid)%cons(ixG^T,1:ncons),ps(igrid)%prim(ixG^T,1:nprim),ps(igrid)%x(ixG^T,1:ndim))

       call phys_get_lfac2(ixG^LL, ixM^LL, ps(igrid), lfac(ixG^T))
       lfac(ixM^T) = dsqrt( lfac(ixM^T) )
       ! calculate angle phi for coordinate transformation
       phi(ixM^T) = datan2( ps(igrid)%x(ixM^T, 2) , ps(igrid)%x(ixM^T, 1) )

       ! calculate cell volume
       dV(ixM^T) = ps(igrid)%prim(ixM^T,psi_)**6*ps(igrid)%dvolume(ixM^T)
       ! calculate total volume
       send_buffer(total_volume) = send_buffer(total_volume) + sum( dV(ixM^T) )

       ! calculate total rest mass
       send_buffer(total_mass) = send_buffer(total_mass) + &
          sum( ps(igrid)%cons(ixM^T,D_) * ps(igrid)%dvolume(ixM^T) )

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

       tmp_2(ixM^T) = ps(igrid)%x(ixM^T, 1) * ( ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,W_vel(2))/lfac(ixM^T) - ps(igrid)%prim(ixM^T,beta(2)) )
       tmp_2(ixM^T) = tmp_2(ixM^T) - &
        ps(igrid)%x(ixM^T, 2) * ( ps(igrid)%prim(ixM^T,alp_) * ps(igrid)%prim(ixM^T,W_vel(1))/lfac(ixM^T) - ps(igrid)%prim(ixM^T,beta(1)) )
       ! calculate rotational energy
       tmp_2(ixM^T) = tmp_2(ixM^T) / (ps(igrid)%x(ixM^T, 1)**2 + ps(igrid)%x(ixM^T, 2)**2) ! here tmp_2 is used to store Omega 
       send_buffer(T_rot) = send_buffer(T_rot) + &
           0.5d0 * sum( tmp_1(ixM^T) * tmp_2(ixM^T) * ps(igrid)%dvolume(ixM^T))


       ! calculate total divB dV
       send_buffer(total_divBdV) = send_buffer(total_divBdV) + &
          sum(dabs(ps(igrid)%prim(ixM^T, divB_)) * dV(ixM^T))

       ! calculate total |B|
       send_buffer(total_B) = send_buffer(total_B) + &
          sum( ps(igrid)%prim(ixM^T,psi_)**2 * dV(ixM^T) &
                * dsqrt( ps(igrid)%prim(ixM^T, Bvec(1))**2 &
                  + ps(igrid)%prim(ixM^T, Bvec(2))**2 &
                  + ps(igrid)%prim(ixM^T, Bvec(3))**2 )&
                    )

       ! calculate total EB_tor
       tmp_1(ixM^T) = ( - ps(igrid)%prim(ixM^T, Bvec(1)) * dsin(phi(ixM^T)) &
                        + ps(igrid)%prim(ixM^T, Bvec(2)) * dcos(phi(ixM^T)) )
       send_buffer(total_EB_tor) = send_buffer(total_EB_tor) + &
          0.5d0 * sum( ps(igrid)%prim(ixM^T,psi_)**4 * dV(ixM^T) * tmp_1(ixM^T)**2 )

       ! calculate total EB
       send_buffer(total_EB) = send_buffer(total_EB) + &
          0.5d0 * sum( ps(igrid)%prim(ixM^T,psi_)**4 * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T, Bvec(1))**2 &
                  + ps(igrid)%prim(ixM^T, Bvec(2))**2 &
                  + ps(igrid)%prim(ixM^T, Bvec(3))**2 )&
                    )

       ! calculate total EE
       send_buffer(total_EE) = send_buffer(total_EE) + &
          0.5d0 * sum( ps(igrid)%prim(ixM^T,psi_)**4 * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T, Evec(1))**2 &
                  + ps(igrid)%prim(ixM^T, Evec(2))**2 &
                  + ps(igrid)%prim(ixM^T, Evec(3))**2 )&
                    )

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
          line = "it global_time dt rho_max alp_min Vol_tot M_rest"
          line = trim(line) // " divB_avg B_avg EB_tor EB_tot EE_tot"
          line = trim(line) // " J_rot T_rot"
          line = trim(line) // " c1_real c1_img c2_real c2_img"
          line = trim(line) // " c3_real c3_img c4_real c4_img"
          line = trim(line) // " rho22_real rho22_img"
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
