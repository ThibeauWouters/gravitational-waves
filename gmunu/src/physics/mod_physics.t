!> This module defines the procedures of a physics module. It contains function
!> pointers for the various supported routines. An actual physics module has to
!> set these pointers to its implementation of these routines.
module mod_physics
  use mod_global_parameters, only: name_len, max_nvar
  use mod_physics_hllc
  use mod_physics_roe
  use mod_physics_ppm

  implicit none
  public

  !> String describing the physics type of the simulation
  character(len=name_len) :: physics_type = ""

  !> To use wider stencils in flux calculations. A value of 1 will extend it by
  !> one cell in both directions, in any dimension
  integer :: phys_wider_stencil = 0

  !> Whether the physics routines require diagonal ghost cells, for example for
  !> computing a curl.
  logical :: phys_req_diagonal = .true.

  !> Array per direction per variable, which can be used to specify that certain
  !> fluxes have to be treated differently
  integer, allocatable :: flux_type(:, :)

  !-------------------------------------------------------------------!
  !          Primitive variables (-1 if not present)
  !-------------------------------------------------------------------!

  ! --------- Hydro variables --------- !
  !> Index of the density (in the prim array)
  integer :: rho_ = -1
  !> Indices of the W * velocities
  integer, allocatable :: W_vel(:)
  !> Indices of the velocities
  integer, allocatable :: veloc(:)
  !> Index of the pressure
  integer :: press_ = -1
  !> Index of the specific internal energy (Note that in GRHD depend directly on eps)
  integer :: eps_ = -1
  !> Index of the sound speed square
  integer :: cs2_ = -1

  ! --------- MHD variables --------- !
  !> Indices of the E field
  integer, allocatable :: Evec(:)
  !> Indices of the B field
  integer, allocatable :: Bvec(:)
  !> Index of the scalar field when using GLM
  integer :: Bphi_ = -1
  integer :: Ephi_ = -1

  ! --------- Radiation Hydro variables --------- !
  !> Indices of the Eddington-factor
  integer :: nu_zeta = -1
  !> Indices of the neutrino energy density
  integer :: nu_E = -1
  !> Indices of the neutrino momentum density
  integer, allocatable :: nu_F(:)
  !> Indices of the neutrino momentum density over energy density
  integer, allocatable :: nu_F_over_E(:)

  ! --------- Metric variables --------- !
  !> Index of the lapse function
  integer :: alp_ = -1
  !> Index of the conformal factor
  integer :: psi_ = -1
  !> Indices of the shift vector 
  integer, allocatable :: beta(:)
  !> Indices of the vector potential
  integer, allocatable :: vecX(:)

  !> Indices of the debug
  integer, allocatable :: debug_vars(:)

  !-------------------------------------------------------------------!
  !          Conserved variables (-1 if not present)
  !-------------------------------------------------------------------!
  !> Indicate where is the Density
  integer :: D_ = -1
  !> Indices of the momentum density
  integer, allocatable :: mom(:)
  !> Indicate where is the Energy
  integer :: tau_ = -1
  !> Indices of the Econs field
  integer, allocatable :: Econs(:)
  !> Indices of the Bcons field
  integer, allocatable :: Bcons(:)
  !> Index of the scalar field when using GLM
  integer :: Bphi_cons_ = -1
  integer :: Ephi_cons_ = -1
  !> Indices of the conserved neutrino energy density
  integer :: nu_Econs = -1
  !> Indices of the conserved neutrino momentum density
  integer, allocatable :: nu_Fcons(:)
  !-------------------------------------------------------------------!

  !> omitted
  integer, parameter   :: flux_nul = -1
  !> Indicates a normal flux
  integer, parameter   :: flux_default        = 0
  !> Indicates the flux should be treated with tvdlf
  integer, parameter   :: flux_tvdlf          = 1
  !> Indicates dissipation should be omitted
  integer, parameter   :: flux_no_dissipation = 2

  !> Type for special methods defined per variable
  type iw_methods
    integer :: test
    !> If this is set, use the routine as a capacity function when adding fluxes
    procedure(sub_get_var), pointer, nopass :: inv_capacity => null()
  end type iw_methods

  !> Special methods defined per variable
  type(iw_methods) :: phys_iw_methods(max_nvar)

  ! hydrodynamics
  procedure(sub_check_params), pointer       :: phys_check_params           => null()
  procedure(sub_update_eos), pointer         :: phys_update_eos             => null()
  procedure(sub_convert), pointer            :: phys_to_conserved           => null()
  procedure(sub_convert), pointer            :: phys_to_primitive           => null()
  procedure(sub_modify_wLR), pointer         :: phys_modify_wLR             => null()
  procedure(sub_modify_flux), pointer        :: phys_modify_flux            => null()
  procedure(sub_get_var_from_state), pointer :: phys_get_lfac2              => null()
  procedure(sub_get_var_from_state), pointer :: phys_get_tilde_S            => null()
  procedure(sub_get_var), pointer            :: phys_get_csound2            => null()
  procedure(sub_get_cmax), pointer           :: phys_get_cmax               => null()
  procedure(sub_get_cbounds), pointer        :: phys_get_cbounds            => null()
  procedure(sub_get_cbounds_and_ct_velocity), pointer    :: phys_get_cbounds_and_ct_velocity => null()
  procedure(sub_get_flux), pointer           :: phys_get_flux               => null()
  procedure(sub_get_dt), pointer             :: phys_get_dt                 => null()
  procedure(sub_add_source_geom), pointer    :: phys_add_source_geom        => null()
  procedure(sub_add_source_geom_fv), pointer :: phys_add_source_geom_fv     => null()
  procedure(sub_add_source), pointer         :: phys_add_source             => null()
  procedure(sub_global_source), pointer      :: phys_global_source_after    => null()
  procedure(sub_special_advance), pointer    :: phys_special_advance        => null()
  procedure(sub_check_prim), pointer         :: phys_check_prim             => null()
  procedure(sub_boundary_adjust), pointer    :: phys_boundary_adjust        => null()
  procedure(sub_write_info), pointer         :: phys_write_info             => null()
  procedure(sub_small_values), pointer       :: phys_handle_small_values    => null()
  procedure(sub_get_var_from_state), pointer :: phys_get_divb               => null()
  procedure(sub_get_var), pointer            :: phys_get_dive               => null()
  procedure(sub_update_faces), pointer       :: phys_update_faces           => null()
  procedure(sub_face_to_center), pointer     :: phys_face_to_center         => null()
  procedure(sub_implicit_update), pointer    :: phys_implicit_update        => null()
  procedure(sub_evaluate_implicit),pointer   :: phys_evaluate_implicit      => null()

  procedure(process_global), pointer         :: phys_clean_divb             => null()
  procedure(process_global), pointer         :: phys_initial_clean_divb     => null()

  ! Called at the beginning of every time step (after determining dt)
  procedure(process_grid), pointer    :: phys_process_grid     => null()
  procedure(process_global), pointer  :: phys_process_global   => null()
  ! Called every time step just after advance (with w^(n+1), it^n, global_time^n)
  procedure(process_grid), pointer   :: phys_process_adv_grid   => null()
  procedure(process_global), pointer :: phys_process_adv_global => null()

  abstract interface

     subroutine sub_check_params
     end subroutine sub_check_params

     subroutine sub_boundary_adjust(igrid,psb)
       use mod_global_parameters
       integer, intent(in) :: igrid
       type(state), target :: psb(max_blocks)
     end subroutine sub_boundary_adjust

     subroutine sub_convert(ixI^L, ixO^L, cons, prim, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout) :: prim(ixI^S, nprim)
       double precision, intent(inout) :: cons(ixI^S, ncons)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
     end subroutine sub_convert

     subroutine sub_update_eos(ixI^L, ixO^L, prim)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout) :: prim(ixI^S, nprim)
     end subroutine sub_update_eos

     subroutine sub_modify_wLR(ixI^L, ixO^L, qt, consL, consR, primL, primR, s, idir)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idir
       double precision, intent(in)    :: qt
       double precision, intent(inout) :: consL(ixI^S,1:ncons), consR(ixI^S,1:ncons)
       double precision, intent(inout) :: primL(ixI^S,1:nprim), primR(ixI^S,1:nprim)
       type(state)                     :: s
     end subroutine sub_modify_wLR

     subroutine sub_modify_flux(ixI^L, ixO^L, consL, consR, primL, primR, xi, s, idir, f)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idir
       double precision, intent(inout) :: consL(ixI^S,1:ncons), consR(ixI^S,1:ncons)
       double precision, intent(inout) :: primL(ixI^S,1:nprim), primR(ixI^S,1:nprim)
       double precision, intent(in)    :: xi(ixI^S, 1:^ND)
       type(state)                     :: s
       double precision, intent(inout) :: f(ixI^S, ncons)
     end subroutine sub_modify_flux

     subroutine sub_get_cmax(prim, x, ixI^L, ixO^L, idim, cmax)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: prim(ixI^S, nprim), x(ixI^S, 1:^ND)
       double precision, intent(inout) :: cmax(ixI^S)
     end subroutine sub_get_cmax

     subroutine sub_get_cbounds(consL, consR, primL, primR, x, ixI^L, ixO^L, idim, cmax, cmin)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: consL(ixI^S, ncons), consR(ixI^S, ncons)
       double precision, intent(in)    :: primL(ixI^S, nprim), primR(ixI^S, nprim)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
       double precision, intent(inout) :: cmax(ixI^S)
       double precision, intent(inout), optional :: cmin(ixI^S)
     end subroutine sub_get_cbounds

     subroutine sub_get_cbounds_and_ct_velocity(consL, consR, primL, primR, x, ixI^L, ixO^L, idim, vcts, cmax, cmin)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: consL(ixI^S, ncons), consR(ixI^S, ncons)
       double precision, intent(in)    :: primL(ixI^S, nprim), primR(ixI^S, nprim)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
       type(ct_velocity), intent(inout):: vcts
       double precision, intent(inout) :: cmax(ixI^S)
       double precision, intent(inout), optional :: cmin(ixI^S)
     end subroutine sub_get_cbounds_and_ct_velocity

     subroutine sub_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: wC(ixI^S, 1:ncons)
       double precision, intent(in)    :: w(ixI^S, 1:nprim)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
       double precision, intent(out)   :: f(ixI^S, ncons)
     end subroutine sub_get_flux

     subroutine sub_add_source_geom(qdt, ixI^L, ixO^L, consCT, primCT, cons, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: qdt, x(ixI^S, 1:^ND)
       double precision, intent(in)    :: consCT(ixI^S, 1:ncons), primCT(ixI^S, 1:nprim)
       double precision, intent(inout) :: cons(ixI^S, 1:ncons)
     end subroutine sub_add_source_geom

     subroutine sub_add_source_geom_fv(qdt, ixI^L, ixO^L, consL, consR, primL, primR, consCT, primCT, cons, x, idims)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idims
       double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
       double precision, intent(in)    :: consL(ixI^S, 1:ncons), primL(ixI^S, 1:nprim)
       double precision, intent(in)    :: consR(ixI^S, 1:ncons), primR(ixI^S, 1:nprim)
       double precision, intent(in)    :: consCT(ixI^S, 1:ncons), primCT(ixI^S, 1:nprim)
       double precision, intent(inout) :: cons(ixI^S, 1:ncons)
     end subroutine sub_add_source_geom_fv

     subroutine sub_add_source(qdt, ixI^L, ixO^L, primCT, cons, x, &
          qsourcesplit, active)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: qdt
       double precision, intent(in)    :: primCT(ixI^S, 1:nprim), x(ixI^S, 1:ndim)
       double precision, intent(inout) :: cons(ixI^S, 1:ncons)
       logical, intent(in)             :: qsourcesplit
       logical, intent(inout)          :: active !< Needs to be set to true when active
     end subroutine sub_add_source

     !> Add global source terms on complete domain (potentially implicit)
     subroutine sub_global_source(qdt, qt, active)
       double precision, intent(in) :: qdt    !< Current time step
       double precision, intent(in) :: qt     !< Current time
       logical, intent(inout)       :: active !< Output if the source is active
     end subroutine sub_global_source

     !> Add special advance in each advect step
     subroutine sub_special_advance(qt, psa)
       use mod_global_parameters
       double precision, intent(in) :: qt     !< Current time
       type(state), target :: psa(max_blocks) !< Compute based on this state
     end subroutine sub_special_advance

     subroutine sub_get_dt(prim, ixI^L, ixO^L, dtnew, dx^D, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
       double precision, intent(in)    :: prim(ixI^S, 1:nprim)
       double precision, intent(inout) :: dtnew
     end subroutine sub_get_dt

     subroutine sub_check_prim(ixI^L, ixO^L, prim, prim_flag)
       use mod_global_parameters
       integer, intent(in)          :: ixI^L, ixO^L
       double precision, intent(in) :: prim(ixI^S,nprim)
       integer, intent(inout)       :: prim_flag(ixG^T)
     end subroutine sub_check_prim

     subroutine sub_write_info(file_handle)
       integer, intent(in) :: file_handle
     end subroutine sub_write_info

     subroutine sub_small_values(prim, x, ixI^L, ixO^L, update_eos, subname)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L,ixO^L
       double precision, intent(inout) :: prim(ixI^S,1:nprim)
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       character(len=*), intent(in)    :: subname
       logical, intent(in)             :: update_eos
     end subroutine sub_small_values

     subroutine sub_get_var_from_state(ixI^L, ixO^L, ps_in, out)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L, ixO^L
       type(state), intent(in)       :: ps_in
       double precision, intent(out) :: out(ixI^S)
     end subroutine sub_get_var_from_state

     subroutine sub_get_var(ixI^L, ixO^L, prim, out)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L, ixO^L
       double precision, intent(in)  :: prim(ixI^S, nprim)
       double precision, intent(out) :: out(ixI^S)
     end subroutine sub_get_var

     subroutine sub_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s,vcts)
       use mod_global_parameters
       integer, intent(in)                :: ixI^L, ixO^L
       double precision, intent(in)       :: qt,qdt
       ! cell-center primitive variables
       double precision, intent(in)       :: wprim(ixI^S,1:nprim)
       ! velocity structure
       type(state)                        :: sCT, s
       type(ct_velocity)                  :: vcts
       double precision, intent(in)       :: fC(ixI^S,1:ncons,1:ndim)
       double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)
     end subroutine sub_update_faces

     subroutine sub_face_to_center(ixO^L,s)
       use mod_global_parameters
       integer, intent(in)                :: ixO^L
       type(state)                        :: s
     end subroutine sub_face_to_center

     !> for processing after the advance (PIC-MHD, e.g.)
     subroutine process_grid(igrid,level,ixI^L,ixO^L,qt,w,x)
       use mod_global_parameters
       integer, intent(in)             :: igrid,level,ixI^L,ixO^L
       double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,1:nprim)
     end subroutine process_grid
 
     !> for processing after the advance (PIC-MHD, e.g.)
     subroutine process_global(iit,qt)
       use mod_global_parameters
       integer, intent(in)          :: iit
       double precision, intent(in) :: qt
     end subroutine process_global

     subroutine sub_evaluate_implicit(qtC,psa)
       use mod_global_parameters
       type(state), target :: psa(max_blocks)   
       double precision, intent(in) :: qtC      
     end subroutine sub_evaluate_implicit

     subroutine sub_implicit_update(dtfactor,qdt,qtC,psa,psb)
       use mod_global_parameters
       type(state), target :: psa(max_blocks)   
       type(state), target :: psb(max_blocks)   
       double precision, intent(in) :: qdt
       double precision, intent(in) :: qtC      
       double precision, intent(in) :: dtfactor 
     end subroutine sub_implicit_update

  end interface

contains

  subroutine phys_check()
    use mod_global_parameters, only: nprim, ndir

    use mod_physics_hllc, only: phys_hllc_check
    use mod_physics_roe, only: phys_roe_check
    use mod_physics_ppm, only: phys_ppm_check

    if (physics_type == "") call mpistop("Error: no physics module loaded")

    call phys_hllc_check()
    call phys_roe_check()
    call phys_ppm_check()

    ! Checks whether the required physics methods have been defined
    if (.not. associated(phys_check_params)) &
         phys_check_params => dummy_check_params

    if (.not. associated(phys_update_eos)) then
         phys_update_eos => dummy_update_eos 
         write(*,*) "WARNING: phys_update_eos is not loaded"
    end if

    if (.not. associated(phys_to_conserved)) &
         call mpistop("Error: phys_to_conserved not defined")

    if (.not. associated(phys_to_primitive)) &
         call mpistop("Error: phys_to_primitive not defined")

    if (.not. associated(phys_modify_wLR)) &
         phys_modify_wLR => dummy_modify_wLR

    if (.not. associated(phys_modify_flux)) &
         phys_modify_flux => dummy_modify_flux

    if (.not. associated(phys_get_cmax)) &
         call mpistop("Error: no phys_get_cmax not defined")

    if (.not. associated(phys_get_cbounds)) &
         call mpistop("Error: no phys_get_cbounds not defined")

    if (.not. associated(phys_get_flux)) &
         call mpistop("Error: no phys_get_flux not defined")

    if (.not. associated(phys_get_dt)) &
         call mpistop("Error: no phys_get_dt not defined")

    if (.not. associated(phys_add_source_geom)) &
         phys_add_source_geom => dummy_add_source_geom

    if (.not. associated(phys_add_source_geom_fv)) &
         phys_add_source_geom_fv => dummy_add_source_geom_fv

    if (.not. associated(phys_add_source)) &
         phys_add_source => dummy_add_source

    if (.not. associated(phys_check_prim)) &
         phys_check_prim => dummy_check_prim

    if (.not. associated(phys_boundary_adjust)) &
         phys_boundary_adjust => dummy_boundary_adjust

    if (.not. associated(phys_write_info)) &
         phys_write_info => dummy_write_info

    if (.not. associated(phys_handle_small_values)) &
         phys_handle_small_values => dummy_small_values

    if (.not. associated(phys_update_faces)) &
         phys_update_faces => dummy_update_faces

    if (.not. associated(phys_face_to_center)) &
         phys_face_to_center => dummy_face_to_center

  end subroutine phys_check

  subroutine dummy_init_params
  end subroutine dummy_init_params

  subroutine dummy_check_params
  end subroutine dummy_check_params

  subroutine dummy_convert(ixI^L, ixO^L, cons, prim, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: prim(ixI^S, nprim)
    double precision, intent(inout) :: cons(ixI^S, ncons)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
  end subroutine dummy_convert

  subroutine dummy_update_eos(ixI^L, ixO^L, prim)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: prim(ixI^S,1:nprim)
  end subroutine dummy_update_eos

  subroutine dummy_modify_wLR(ixI^L, ixO^L, qt, consL, consR, primL, primR, s, idir)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: consL(ixI^S,1:ncons), consR(ixI^S,1:ncons)
    double precision, intent(inout) :: primL(ixI^S,1:nprim), primR(ixI^S,1:nprim)
    type(state)                     :: s
  end subroutine dummy_modify_wLR

  subroutine dummy_modify_flux(ixI^L, ixO^L, consL, consR, primL, primR, xi, s, idir, f)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(inout) :: consL(ixI^S,1:ncons), consR(ixI^S,1:ncons)
    double precision, intent(inout) :: primL(ixI^S,1:nprim), primR(ixI^S,1:nprim)
    double precision, intent(in)    :: xi(ixI^S, 1:^ND)
    type(state)                     :: s
    double precision, intent(inout) :: f(ixI^S, ncons)
  end subroutine dummy_modify_flux

  subroutine dummy_add_source_geom(qdt, ixI^L, ixO^L, consCT, primCT, cons, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: consCT(ixI^S, 1:ncons), primCT(ixI^S, 1:nprim)
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
  end subroutine dummy_add_source_geom

  subroutine dummy_add_source_geom_fv(qdt, ixI^L, ixO^L, consL, consR, primL, primR, consCT, primCT, cons, x, idims)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idims
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(in)    :: consL(ixI^S, 1:ncons), primL(ixI^S, 1:nprim)
    double precision, intent(in)    :: consR(ixI^S, 1:ncons), primR(ixI^S, 1:nprim)
    double precision, intent(in)    :: consCT(ixI^S, 1:ncons), primCT(ixI^S, 1:nprim)
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
  end subroutine dummy_add_source_geom_fv

  subroutine dummy_add_source(qdt, ixI^L, ixO^L, primCT, cons, x, &
       qsourcesplit, active)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: primCT(ixI^S, 1:nprim), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active
    ! Don't have to set active, since it starts as .false.
  end subroutine dummy_add_source

  subroutine dummy_check_prim(ixI^L, ixO^L, prim, w_flag)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: prim(ixI^S,nprim)
    integer, intent(inout)       :: w_flag(ixG^T)

    w_flag(ixO^S) = 0             ! All okay
  end subroutine dummy_check_prim

  subroutine dummy_boundary_adjust(igrid,psb)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state), target :: psb(max_blocks)
  end subroutine dummy_boundary_adjust

  subroutine dummy_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh !< File handle
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    ! Number of physics parameters
    integer, parameter                  :: n_par = 0

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)
  end subroutine dummy_write_info

  subroutine dummy_small_values(prim, x, ixI^L, ixO^L, update_eos, subname)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: prim(ixI^S,1:nprim)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    logical, intent(in)             :: update_eos
  end subroutine dummy_small_values

  subroutine dummy_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixI^S,1:nprim)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixI^S,1:ncons,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)
  end subroutine dummy_update_faces

  subroutine dummy_face_to_center(ixO^L,s)
    use mod_global_parameters
    integer, intent(in)                :: ixO^L
    type(state)                        :: s
  end subroutine dummy_face_to_center

  subroutine dummy_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)   
    double precision, intent(in) :: qtC      
    integer :: iigrid, igrid

    ! Just copy in nul state
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       psa(igrid)%prim = 0.0d0*psa(igrid)%prim
       if(stagger_grid) psa(igrid)%conss = 0.0d0*psa(igrid)%conss
    end do
    !$OMP END PARALLEL DO

  end subroutine dummy_evaluate_implicit

  subroutine dummy_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)   
    type(state), target :: psb(max_blocks)   
    double precision, intent(in) :: qdt      
    double precision, intent(in) :: qtC      
    double precision, intent(in) :: dtfactor 
    integer :: iigrid, igrid

    ! Just copy in psb state when using the scheme without implicit part
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       psa(igrid)%cons = psb(igrid)%cons
       if(stagger_grid) psa(igrid)%conss = psb(igrid)%conss
    end do
    !$OMP END PARALLEL DO

  end subroutine dummy_implicit_update


end module mod_physics
