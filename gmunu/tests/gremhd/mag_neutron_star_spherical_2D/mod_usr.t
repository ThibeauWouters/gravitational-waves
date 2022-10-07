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
   
  !logical, parameter  ::  initialize_mag_field = .False.
  !logical, parameter  ::  initialize_metric = .False.
  !logical, parameter  ::  evolve_metric = .False.
  logical, parameter  ::  initialize_mag_field = .True.
  logical, parameter  ::  initialize_metric = .True.
  logical, parameter  ::  evolve_metric = .True.

  ! conducivity in the center of the star
  double precision, parameter  ::  sigma_0 = 1.0d6

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos
    use mod_cfc

    !character*128  ::  profile_path = '/users/chi-kit.cheong/simulations/initial_profile/T1K1'
    character*128  ::  profile_path = '/users/chi-kit.cheong/simulations/initial_profile/mag_ns_pol_weak'
    integer  ::  Nth, Nr   ! number of grid of the profile


    usr_init_one_grid => ns_init_one_grid
    usr_improve_initial_condition => ns_improve_initial_condition
    usr_get_resistivity => my_get_resistivity
    usr_process_adv_grid => update_aux_vars
    usr_refine_grid => my_refine
    usr_print_log => printlog

    call set_coordinate_system("spherical_2.5D")
    call gremhd_activate()
    ! use atmosphere
    call eos_atmo_activate()
    call eos_idealgas_activate()

    divB_ = var_set_auxvar('divB')
    divE_ = var_set_auxvar('divE')

    call initialize_gmunu()

    call mod_XNS_read_profile(profile_path)

    if ( prof%Nth == 1) then 
       call mpistop("! the profile is 1D.")
    else
       if ( ndir < 2 ) &
          call mpistop("! Dimension = 1 but the profile is 2D.")
    end if

    ! setup atmosphere density from the profile
    call eos_initialize_atmo(maxval(prof%rho))

    if (initialize_metric .or. evolve_metric) then
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
    double precision                :: E_i(ixI^S,1:ndir)
    double precision                :: lfac(ixI^S), gamma(ixI^S,1:3,1:3)
    double precision                :: psi_tmp
    double precision                :: xi(ixI^S,1:ndim)
    integer                         :: ix^D, idir, ixC^L
    double precision                :: r_star, delta ! for perturbation

    {^IFONED call mpistop("This is a multi-D HD problem") }

    w = 0.0d0
    v = 0.0d0
    w(ixO^S, divB_) = smalldouble

    {do ix^D = ixO^LIM^D \} 
       call mod_XNS_map_2D(w(ix^D,rho_),x(ix^D,r_), x(ix^D,theta_),prof%rho(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,press_),x(ix^D,r_), x(ix^D,theta_),prof%press(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,alp_),x(ix^D,r_), x(ix^D,theta_),prof%alp(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,psi_),x(ix^D,r_), x(ix^D,theta_),prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,beta(3)),x(ix^D,r_), x(ix^D,theta_),prof%beta3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(v(ix^D,3),x(ix^D,r_), x(ix^D,theta_),prof%v3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       
       call mod_XNS_map_2D(w(ix^D,Bvec(3)),x(ix^D,r_), x(ix^D,theta_),prof%B3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,Bvec(2)),x(ix^D,r_), x(ix^D,theta_),prof%B2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(w(ix^D,Bvec(1)),x(ix^D,r_), x(ix^D,theta_),prof%B1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)

       call mod_XNS_map_2D(E_i(ix^D,3),x(ix^D,r_), x(ix^D,theta_),prof%E_3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(E_i(ix^D,2),x(ix^D,r_), x(ix^D,theta_),prof%E_2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
       call mod_XNS_map_2D(E_i(ix^D,1),x(ix^D,r_), x(ix^D,theta_),prof%E_1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
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
           call mod_XNS_map_2D( psi_tmp ,xi(ix^D,r_), xi(ix^D,theta_),prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
           select case (idir)
           case (1)
              call mod_XNS_map_2D( block%conss(ix^D, idir) ,xi(ix^D,r_), xi(ix^D,theta_),prof%B1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
           case (2)
              call mod_XNS_map_2D( block%conss(ix^D, idir) ,xi(ix^D,r_), xi(ix^D,theta_),prof%B2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
           end select
           ! conformally rescale
           block%conss(ix^D, idir) = psi_tmp**6 * block%conss(ix^D, idir) 
         {end do^D&\} 
       end do

       ! although B field at cell centre are filled, we still update it based on B field at cell face
       call phys_face_to_center(ixO^L, block) 
    end if
 
    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_)**4 
    end do

    do idir = 1, ndir
       w(ixO^S,Evec(idir)) = E_i(ixO^S, idir) / gamma(ixO^S,idir,idir)
    end do

    ! add perturbation on velocities
    delta = 0.0d-3
    {do ix^D = ixO^LIM^D \}
          call mod_XNS_map_1D(r_star,x(ix^D,theta_),prof%r_star(:),prof%theta,prof%Nth) ! firstly, get r_star(theta)
          ! here we perturbe v^r as an example
          v(ix^D,1) = delta / gamma(ix^D, 1, 1) * dsin( dpi * x(ix^D,r_) / r_star )
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

    ! reset eps
    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {end do^D&\}
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
  end subroutine ns_improve_initial_condition

  ! get resistivity
  subroutine my_get_resistivity(ixI^L,ixO^L,cons,x,eta)
    use mod_global_parameters
    use mod_eos
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(in) :: cons(ixI^S,1:ncons) 
    double precision, intent(out):: eta(ixI^S) 

    double precision             :: tmp
    integer                      :: ix^D

    {do ix^D = ixM^LLIM^D \}
       tmp = 1.0d0 - small_rho / max( cons(ix^D,D_), small_D )
       tmp = min( max( smalldouble, tmp ), 1.0d0 )
       tmp = sigma_0 * tmp**2
       eta(ix^D) = min( 1.0d0 / tmp, sigma_0 )
       !eta(ix^D) = 1.0d0 / tmp
    {end do^D&\}
  end subroutine my_get_resistivity

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

    if ( minval(x(ixO^S,1)) <= 1.0d0 ) then
       ratio(ixO^S) = dx(2, level) * x(ixO^S,1) / dx(1, level)
       ratio(ixO^S) = ( ratio(ixO^S) - 1.0d0 )
       if ( any( ratio(ixO^S) <= smalldouble ) ) then
          refine = -1
          coarsen = 1
       else
          !refine = 1
          !coarsen = -1
       end if
    end if
  end subroutine my_refine

  subroutine update_aux_vars(igrid,level,ixI^L,ixO^L,qt,s)
    use mod_global_parameters
    use mod_gremhd_phys_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s
    double precision                :: divb(ixI^S), dive(ixI^S)
    !> calculate div B after the advance 
    call gremhd_get_divb(ixI^L, ixO^L, s, divb(ixI^S))
    s%prim(ixO^S, divB_) = dabs( divb(ixO^S) )
    ! also the div E
    call gremhd_get_div(ixI^L, ixO^L, s%prim(ixI^S,Evec(1:ndir)), dive(ixI^S), dive_4thorder )
    s%prim(ixO^S, divE_) = ( dive(ixO^S) )
  end subroutine update_aux_vars

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
    double precision     :: dV(ixG^T)
    double precision     :: lfac(ixG^T)
    double precision     :: h(ixG^T)

    double precision     :: rho_max, rho_max_local
    double precision     :: alp_min, alp_min_local

    integer              :: total_volume = -1
    integer              :: total_mass = -1
    integer              :: J_rot = -1
    integer              :: T_rot = -1
    integer              :: total_B = -1
    integer              :: total_EB_tor = -1
    integer              :: total_EB_pol = -1
    integer              :: total_EE = -1
    integer              :: total_divBdV = -1
    integer              :: total_divEdV = -1
    integer              :: total_I11dot = -1
    integer              :: total_I12dot = -1
    integer              :: total_I22dot = -1
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
    total_divEdV = add_var()
    total_B = add_var()
    total_EB_tor = add_var()
    total_EB_pol = add_var()
    total_EE = add_var()
    J_rot = add_var()
    T_rot = add_var()
    total_I11dot = add_var()
    total_I12dot = add_var()
    total_I22dot = add_var()

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

       ! calculate total volume
       dV(ixM^T) = ps(igrid)%prim(ixM^T,psi_)**6*ps(igrid)%dvolume(ixM^T)
       send_buffer(total_volume) = send_buffer(total_volume) + sum( dV(ixM^T) )

       ! calculate total rest mass
       send_buffer(total_mass) = send_buffer(total_mass) + &
          sum( ps(igrid)%cons(ixM^T,D_) * ps(igrid)%dvolume(ixM^T) )

       ! fixme: may not correct
       ! calculate angular momentum 
       send_buffer(J_rot) = send_buffer(J_rot) + &
          sum( ps(igrid)%cons(ixM^T,mom(3)) * ps(igrid)%dvolume(ixM^T) )

       ! calculate kinetic energy (related to v phi only!) 
       send_buffer(T_rot) = send_buffer(T_rot) + &
          0.5d0 * sum( ps(igrid)%prim(ixM^T,W_vel(3))/lfac(ixM^T) * ps(igrid)%cons(ixM^T,mom(3)) * ps(igrid)%dvolume(ixM^T) )

       ! calculate total |B|
       send_buffer(total_B) = send_buffer(total_B) + &
          sum( ps(igrid)%prim(ixM^T,psi_)**2 * dV(ixM^T) &
                * dsqrt( ps(igrid)%prim(ixM^T, Bvec(1))**2 &
                  + ps(igrid)%x(ixM^T,1)**2 * ps(igrid)%prim(ixM^T, Bvec(2))**2 &
                  + ps(igrid)%x(ixM^T,1)**2 * dsin(ps(igrid)%x(ixM^T,2))**2 * ps(igrid)%prim(ixM^T, Bvec(3))**2 )&
                    )

       ! calculate total divB dV
       send_buffer(total_divBdV) = send_buffer(total_divBdV) + &
          sum(dabs(ps(igrid)%prim(ixM^T, divB_)) * dV(ixM^T))

       ! calculate total divE dV
       send_buffer(total_divEdV) = send_buffer(total_divEdV) + &
          sum(dabs(ps(igrid)%prim(ixM^T, divE_)) * dV(ixM^T))

       ! calculate total EB_tor
       send_buffer(total_EB_tor) = send_buffer(total_EB_tor) + &
          0.5d0 * sum( ps(igrid)%x(ixM^T,1)**2 * dsin(ps(igrid)%x(ixM^T,2))**2 * ps(igrid)%prim(ixM^T,psi_)**4 &
                * ps(igrid)%prim(ixM^T, Bvec(3))**2 * dV(ixM^T) )

       ! calculate total EB_pol
       send_buffer(total_EB_pol) = send_buffer(total_EB_pol) + &
          0.5d0 * sum( ps(igrid)%prim(ixM^T,psi_)**4 &
                * (ps(igrid)%prim(ixM^T, Bvec(1))**2 &
                  + ps(igrid)%x(ixM^T,1)**2 * ps(igrid)%prim(ixM^T, Bvec(2))**2) * dV(ixM^T) )

       ! calculate total EE
       send_buffer(total_EE) = send_buffer(total_EE) + &
          0.5d0 * sum( ps(igrid)%prim(ixM^T,psi_)**4 * dV(ixM^T) &
                * ( ps(igrid)%prim(ixM^T, Evec(1))**2 &
                  + ps(igrid)%x(ixM^T,1)**2 * ps(igrid)%prim(ixM^T, Evec(2))**2 &
                  + ps(igrid)%x(ixM^T,1)**2 * dsin(ps(igrid)%x(ixM^T,2))**2 * ps(igrid)%prim(ixM^T, Evec(3))**2 )&
                    )

       ! fixme: may not correct
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

    ! Volume average of divB
    recive_buffer(total_divBdV) = (recive_buffer(total_divBdV) / recive_buffer(total_volume))
    ! Volume average of divE
    recive_buffer(total_divEdV) = (recive_buffer(total_divEdV) / recive_buffer(total_volume))
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
          line = trim(line) // " divB_avg divE_avg B_avg"
          line = trim(line) //  " EB_tor EB_pol EE" 
          line = trim(line) //  " J_rot T_rot" 
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
