!> General relativistic radiation hydrodynamics physics module (only in CFC now)
module mod_grrhd_phys
  use mod_physics
  use mod_grrhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grrhd_phys_init

contains

  !> Read this module's parameters from a file
  subroutine grrhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grrhd_list/ evolve_hydro, &
                           use_GR, tolerance, iter_max, &
                           use_radiation

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grrhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine grrhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine grrhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = 1.0d0!grrhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine grrhd_write_info


  !> Initialize the module
  subroutine grrhd_phys_init()
    use mod_global_parameters
    use mod_grrhd_phys_convert
    use mod_grrhd_phys_flux
    use mod_grrhd_phys_add_source

    integer :: itr, idir

    call grrhd_read_params(par_files)

    physics_type = "grrhd"
  
    ! Determine primitive variables that needed to be reconstructed
    rho_ = var_set_primvar('rho', var_type=hydro_var)
    allocate(W_vel(ndir))
    do idir = 1, ndir
       W_vel(idir) = var_set_primvar('W_vel',idir, var_type=hydro_var)
    end do
    eps_ = var_set_primvar('eps', var_type=hydro_var,need_bc=.false.)
    if ( use_radiation ) then
       nu_E = var_set_primvar('nu_E', var_type=hydro_var,need_bc=.false.)
       allocate(nu_F_over_E(ndir))
       do idir = 1, ndir
          nu_F_over_E(idir) = var_set_primvar('nu_F_over_E',idir, var_type=hydro_var,need_bc=.false.)
       end do
    end if

    ! Determine primitive variables that needed NOT to be reconstructed
    press_ = var_set_primvar('press', var_type=hydro_var,need_bc=.false.,need_rec=.false.)
    cs2_ = var_set_primvar('cs2', var_type=hydro_var,need_bc=.false.,need_rec=.false.)
    if ( use_radiation ) then
       ! fixme: nu_zeta actually needed to be reconstructed, 
       ! but normally will not use the same scheme as in hydro.
       ! Most of the time PC is used(first order scheme), so I put it here
       nu_zeta = var_set_primvar('nu_zeta', var_type=hydro_var,need_bc=.false.)
    end if


    ! Set index of metric variables, it has to be reconstructed so use fluxvar
    alp_ = var_set_primvar("alp", var_type=metric_var, need_bc=.False.,need_rec=.false.)
    psi_ = var_set_primvar("psi", var_type=metric_var, need_bc=.False.,need_rec=.false.)
    allocate(beta(ndir))
    do idir=1,ndir
       beta(idir) = var_set_primvar("beta", idir, var_type=metric_var, need_bc=.False.,need_rec=.false.)
    enddo
    allocate(vecX(ndir))
    do idir=1,ndir
       vecX(idir) = var_set_primvar("X", idir, var_type=metric_var, need_bc=.False.,need_rec=.false.)
    enddo

    ! Determine flux (cons) variables
    D_ = var_set_consvar('D', var_type=hydro_var)
    allocate(mom(ndir))
    do idir = 1, ndir
       mom(idir) = var_set_consvar('mom',idir, var_type=hydro_var)
    end do
    tau_ = var_set_consvar('tau', var_type=hydro_var)
    if ( use_radiation ) then
       nu_Econs = var_set_consvar('nu_Econs', var_type=hydro_var)
       allocate(nu_Fcons(ndir))
       do idir = 1, ndir
          nu_Fcons(idir) = var_set_consvar('nu_Fcons',idir, var_type=hydro_var)
       end do
    end if

    phys_get_dt              => grrhd_get_dt
    phys_get_lfac2           => grrhd_get_lfac2
    phys_get_tilde_S         => grrhd_get_tilde_S
    phys_get_csound2         => grrhd_get_csound2
    phys_get_cmax            => grrhd_get_cmax
    phys_update_eos          => grrhd_update_eos
    phys_check_params        => grrhd_check_params
    phys_check_prim          => grrhd_check_prim
    phys_write_info          => grrhd_write_info
    phys_handle_small_values => grrhd_handle_small_values

    call grrhd_phys_convert_init()
    call grrhd_phys_flux_init()
    call grrhd_phys_add_source_init()
    if ( use_radiation ) then
       phys_modify_wLR            => grrhd_modify_wLR
    end if

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, ncons))
    else if (any(shape(flux_type) /= [ndir, ncons])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    flux_type = flux_default

    if ( .not.evolve_hydro ) then
       flux_type(:, D_) = flux_nul
       flux_type(:, tau_) = flux_nul
       flux_type(:, mom(1:ndir)) = flux_nul
    end if

    nvector      = 3 ! No. vector vars
    if ( use_radiation ) then
       nvector = nvector + 1
    end if
    allocate(iw_vector(nvector))
    iw_vector(1) = W_vel(1) - 1
    iw_vector(2) = beta(1) - 1
    iw_vector(3) = vecX(1) - 1
    if ( use_radiation ) then
       iw_vector(4) = nu_F_over_E(1) - 1
    end if

  end subroutine grrhd_phys_init

  subroutine grrhd_check_params
    use mod_global_parameters
    use mod_eos
    use mod_usr_methods

    ! check whether eos module has been loaded
    call eos_check
    if ( use_radiation ) then
       if (.not.use_imex_scheme ) then
          call mpistop("Error: currently only IMEX schemes are supported when use radiation.")
       end if

       if(.not. associated(usr_get_opacities)) then
          write(*,*) "Error: usr_get_opacities must be defined in GRRHD"
          !call mpistop("Error: usr_get_opacities must be defined in GRRHD")
       end if
    end if

  end subroutine grrhd_check_params

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine grrhd_update_eos(ixI^L, ixO^L, prim)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)                       :: ixI^L, ixO^L
    double precision, intent(inout)           :: prim(ixI^S, 1:nprim)
    integer                                   :: ix^D

    ! rho and eps are given, update the rest of the primitive variables
    {do ix^D = ixO^LIM^D \}
       call grrhd_update_eos_one_point(prim(ix^D,:))
    {enddo^D&\}
  end subroutine grrhd_update_eos

  !> Calculate cmax_idim within ixO^L
  subroutine grrhd_get_cmax(prim, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    integer                                   :: idir, ix^D
    double precision                          :: lambda(ixI^S,1:2)
    call grrhd_get_lambda(ixI^L, ixO^L, idim, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambda(ixI^S,1:2))
    cmax(ixO^S) = max(abs(lambda(ixO^S,1)), abs(lambda(ixO^S,2)))

  end subroutine grrhd_get_cmax

  subroutine grrhd_get_csound2(ixI^L,ixO^L,prim,csound2)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: prim(ixI^S,nprim)
    double precision, intent(out)   :: csound2(ixI^S)
    integer                         :: idir, ix^D
    ! Note: the input are the prim variables
    {do ix^D = ixO^LIM^D \}
       ! update cs2
       call eos_get_cs2_one_grid(csound2(ix^D),prim(ix^D, rho_),prim(ix^D, eps_))
       ! strictly require cs2 is physical
       csound2(ix^D) = max( min( 1.0d0, csound2(ix^D) ) , 0.0d0)
    {enddo^D&\}
  end subroutine grrhd_get_csound2

  subroutine grrhd_get_lfac2(ixI^L,ixO^L,ps_in,lfac2)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixI^L, ixO^L
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: lfac2(ixI^S)
    integer                         :: idir, ix^D
    double precision                :: gamma(ixI^S,1:3,1:3)

    associate(prim=>ps_in%prim,x=>ps_in%x)
    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:^ND), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
    end do
    
    lfac2(ixO^S) = 1.0d0 
    ! calculate W^2 = 1 + W^2 * v^2
    do idir = 1, ndir
       lfac2(ixO^S) = lfac2(ixO^S) + gamma(ixO^S,idir,idir) * prim(ixO^S, W_vel(idir))**2
    end do
    end associate
  end subroutine grrhd_get_lfac2

  ! get psi^6 (\gamma_{ij}S^{ij}), which is part of the source terms in cfc_alp
  subroutine grrhd_get_tilde_S(ixI^L,ixO^L,ps_in,tilde_S)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixI^L, ixO^L
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: tilde_S(ixI^S)

    double precision                :: h(ixI^S)   
    double precision                :: gamma(ixI^S,1:3,1:3)
    integer                         :: idir

    associate(prim=>ps_in%prim,x=>ps_in%x)
    call grrhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                h=h(ixI^S) )

    ! nota that S = S_i v^i + 3 p
    tilde_S(ixO^S) = 3.0d0 * prim(ixO^S, press_)
    do idir = 1, ndir
       tilde_S(ixO^S) = tilde_S(ixO^S) &
        + ( ( prim(ixO^S, rho_) * h(ixO^S) ) * prim(ixO^S, W_vel(idir)) ) & 
                 * gamma(ixO^S,idir,idir) * prim(ixO^S, W_vel(idir))

    end do
    tilde_S(ixO^S) = tilde_S(ixO^S) * prim(ixO^S, psi_)**6
    end associate
  end subroutine grrhd_get_tilde_S

  subroutine grrhd_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nprim)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

  end subroutine grrhd_get_dt

  !> Returns 0 in argument flag where values are ok
  subroutine grrhd_check_prim(ixI^L, ixO^L, prim, flag)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: prim(ixI^S, 1:nprim)
    integer, intent(inout)       :: flag(ixI^S)
    flag(ixO^S) = 0
    ! this subroutine is used only in MP5, so we care about rho and eps only.
    where(prim(ixO^S, rho_) < smalldouble ) 
       flag(ixO^S) = rho_
    else where(prim(ixO^S, eps_) < smalldouble ) 
       flag(ixO^S) = eps_
    end where
  end subroutine grrhd_check_prim

  subroutine grrhd_modify_wLR(ixI^L, ixO^L, qt, consL, consR, primL, primR, s, idir)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: consL(ixI^S,1:ncons), consR(ixI^S,1:ncons)
    double precision, intent(inout) :: primL(ixI^S,1:nprim), primR(ixI^S,1:nprim)
    type(state)                     :: s

    call grrhd_get_nu_zeta(ixI^L, ixO^L, primL, s%x)
    call grrhd_get_nu_zeta(ixI^L, ixO^L, primR, s%x)
  end subroutine grrhd_modify_wLR

end module mod_grrhd_phys
