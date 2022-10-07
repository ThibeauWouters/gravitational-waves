!> General relativistic hydrodynamics physics module (only in CFC now)
module mod_gremhd_phys
  use mod_physics
  use mod_gremhd_phys_parameters

  implicit none
  private


  ! Public methods
  public :: gremhd_phys_init
  public :: gremhd_get_div
  public :: gremhd_b_from_vector_potential
  public :: gremhd_get_divb

contains

  !> Read this module's parameters from a file
  subroutine gremhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n
    character(len=name_len)      :: type_divb_fix = "", type_dive_fix = ""

    namelist /gremhd_list/ evolve_hydro, evolve_EM, &
                           use_GR, tolerance, iter_max, lfac_max, &
                           rescale_f, tol_im, iter_max_im, &
                           elliptic_cleaning, &
                           type_dive_fix, dive_4thorder, &
                           type_divb_fix, divb_4thorder, &
                           divB_glm_kappa, divE_glm_kappa, &
                           divB_mg_n_cycle, divB_mg_redblack, divB_mg_tol, divB_mg_it_max, &
                           type_ct, face_to_center_order

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, gremhd_list, end=111)
111    close(unitpar)
    end do

    if (mype == 0) then
       if ( .not. evolve_hydro ) write(*,*) "Hydro is not evolving"
       if ( .not. evolve_EM ) write(*,*) "EM is not evolving"
    end if

    select case (type_divb_fix)
    case ('none')
       type_divb = none
    case ('multigrid')
       type_divb = multigrid
    case ('glm')
       type_divb = glm
    case ('ct')
       type_divb = ct
    case default
       if (mype == 0) write(*,*) "Warning: type_divb_fix is not specificied, now we assume type_divb = none!"
       type_divb = none
    end select

    select case (type_dive_fix)
    case ('none')
       type_dive = none
    case ('glm')
       type_dive = glm
    case default
       if (mype == 0) write(*,*) "Warning: type_dive_fix is not specificied, now we assume type_dive = none!"
       type_dive = none
    end select

  end subroutine gremhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine gremhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = 1.0d0!gremhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine gremhd_write_info


  !> Initialize the module
  subroutine gremhd_phys_init()
    use mod_global_parameters
    use mod_cfc, only: reconstruct_cfc
    use mod_gremhd_phys_convert
    use mod_gremhd_phys_flux
    use mod_gremhd_phys_add_source
    use mod_gremhd_phys_implicit_update
    use mod_gremhd_phys_divb_mg
    use mod_gremhd_phys_divb_ct

    integer :: itr, idir

    call gremhd_read_params(par_files)

    physics_type = "gremhd"
  
    ! Determine primitive variables 
    rho_ = var_set_primvar('rho', var_type=hydro_var)
    allocate(W_vel(ndir))
    do idir = 1, ndir
       W_vel(idir) = var_set_primvar('W_vel',idir, var_type=hydro_var)
    end do
    allocate(Evec(ndir))
    do idir = 1, ndir
       Evec(idir) = var_set_primvar('Evec',idir, var_type=hydro_var)
    end do

    allocate(Bvec(ndir))
    ! in the normal cases, B needs to be reconstructed
    ! B is stored in stagger grid, no needed to be reconstructed
    do idir = 1, ndir
       Bvec(idir) = var_set_primvar('Bvec',idir, var_type=hydro_var,need_rec=.not.stagger_grid)
    end do

    eps_ = var_set_primvar('eps', var_type=hydro_var,need_bc=.false.)
    if (type_divb == glm) then
       ! introduce a scalar field Phi for divergence cleaning
       Bphi_ = var_set_primvar('Bphi', var_type=hydro_var,need_bc=.false.)
    end if
    if (type_dive == glm) then
       ! introduce a scalar field Phi for divergence cleaning
       Ephi_ = var_set_primvar('Ephi', var_type=hydro_var,need_bc=.false.)
    end if

    press_ = var_set_primvar('press', var_type=hydro_var,need_bc=.false.,need_rec=.false.)
    cs2_ = var_set_primvar('cs2', var_type=hydro_var,need_bc=.false.,need_rec=.false.)

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

    ! Determine hydro flux variables
    D_ = var_set_consvar('D', var_type=hydro_var)
    allocate(mom(ndir))
    do idir = 1, ndir
       mom(idir) = var_set_consvar('mom',idir, var_type=hydro_var)
    end do
    tau_ = var_set_consvar('tau', var_type=hydro_var)
    allocate(Econs(ndir))
    do idir = 1, ndir
       Econs(idir) = var_set_consvar('Econs',idir, var_type=hydro_var)
    end do
    allocate(Bcons(ndir))
    do idir = 1, ndir
       Bcons(idir) = var_set_consvar('Bcons',idir, var_type=hydro_var)
    end do
    if (type_divb == glm) then
       ! introduce a scalar field Phi for divergence cleaning
       Bphi_cons_ = var_set_consvar('Bphi_cons', var_type=hydro_var)
    end if
    if (type_dive == glm) then
       ! introduce a scalar field Phi for divergence cleaning
       Ephi_cons_ = var_set_consvar('Ephi_cons', var_type=hydro_var)
    end if

    phys_get_dt              => gremhd_get_dt
    phys_get_lfac2           => gremhd_get_lfac2
    phys_get_tilde_S         => gremhd_get_tilde_S
    phys_get_csound2         => gremhd_get_csound2
    phys_get_cmax            => gremhd_get_cmax
    phys_get_divb            => gremhd_get_divb
    phys_update_eos          => gremhd_update_eos
    phys_check_params        => gremhd_check_params
    phys_check_prim          => gremhd_check_prim
    phys_write_info          => gremhd_write_info
    phys_handle_small_values => gremhd_handle_small_values

    call gremhd_phys_convert_init()
    call gremhd_phys_flux_init()
    call gremhd_phys_add_source_init()
    ! fixme: if not evolve_EM, divB might have some problem
    if ( evolve_EM ) call gremhd_phys_implicit_update_init()

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, ncons))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, ncons])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .False.
    select case (type_divb)
    case (glm)
       ! nothing
    case (multigrid)
       elliptic_cleaning = .True.
    case (ct)
       call gremhd_phys_divb_ct_init()
       ! diagonal ghost cell is needed when use CT
       phys_req_diagonal = .true.
       ! flux of B field is handled with CT, no need to evolve with other hydro
       flux_type(:,Bcons(1:ndim)) = flux_nul
       ! if use ct, the tol should able to go below double precision
       divB_mg_tol = min(1.0d-16, divB_mg_tol)
    case default
       ! nothing
    end select

    if ( .not.evolve_hydro ) then
       flux_type(:, D_) = flux_nul
       flux_type(:, tau_) = flux_nul
       flux_type(:, mom(1:ndir)) = flux_nul
    end if
    if ( .not.evolve_EM ) then
       flux_type(:, Econs(1:ndir)) = flux_nul
       flux_type(:, Bcons(1:ndir)) = flux_nul
    end if

    if (elliptic_cleaning) then
       call gremhd_phys_divb_mg_init()
       ! setup div cleaning here. fixme: maybe need to move it our later
       phys_clean_divb => gremhd_clean_divb_multigrid
    end if

    ! type of inital clean divb, currently we only have mg
    phys_initial_clean_divb => gremhd_initial_clean_divb_multigrid

    nvector      = 5 ! No. vector vars
    allocate(iw_vector(1:nvector))
    iw_vector(1) = W_vel(1) - 1
    iw_vector(2) = beta(1) - 1 
    iw_vector(3) = vecX(1) - 1
    iw_vector(4) = Evec(1) - 1
    iw_vector(5) = Bvec(1) - 1

  end subroutine gremhd_phys_init

  subroutine gremhd_check_params
    use mod_global_parameters
    use mod_eos
    use mod_usr_methods

    ! check whether eos module has been loaded
    call eos_check

    if (evolve_EM) then
       if (.not.use_imex_scheme ) then
          call mpistop("Error: currently only IMEX schemes are supported in gremhd.")
       end if
       
       if(.not. associated(usr_get_resistivity)) then
          call mpistop("Error: usr_get_resistivity must be defined in GRRMHD")
       end if
    end if
 
    !fixme: add more check on ct

  end subroutine gremhd_check_params

  !> Update the eos (prim) variables, p, temp, entropy, cs2
  subroutine gremhd_update_eos(ixI^L, ixO^L, prim)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L
    double precision, intent(inout)           :: prim(ixI^S, 1:nprim)
    integer                                   :: ix^D
    ! rho and eps are given, update the rest of the primitive variables
    {do ix^D = ixO^LIM^D \}
       call gremhd_update_eos_one_point(prim(ix^D, 1:nprim))
    {enddo^D&\}
  end subroutine gremhd_update_eos

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine gremhd_get_cmax(prim, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: prim(ixI^S, nprim), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    integer                                   :: idir, ix^D
    double precision                          :: lambda(ixI^S,1:2)
    call gremhd_get_lambda(ixI^L, ixO^L, idim, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambda(ixI^S,1:2))
    cmax(ixO^S) = max(abs(lambda(ixO^S,1)), abs(lambda(ixO^S,2)))
  end subroutine gremhd_get_cmax

  subroutine gremhd_get_csound2(ixI^L,ixO^L,prim,csound2)
    use mod_global_parameters
    use mod_eos
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: prim(ixI^S,1:nprim)
    double precision, intent(out)   :: csound2(ixI^S)
    integer                         :: idir, ix^D
    ! Note: the input are the prim variables
    {do ix^D = ixO^LIM^D \}
       ! update cs2
       call eos_get_cs2_one_grid(csound2(ix^D),prim(ix^D, rho_),prim(ix^D, eps_))
       ! strictly require cs2 is physical
       if ( (csound2(ix^D) > 1.0d0).or.(csound2(ix^D) < 0.0d0) ) then
          csound2(ix^D) = 0.0d0
       end if
    {enddo^D&\}
  end subroutine gremhd_get_csound2

  ! get psi^6 (\gamma_{ij}S^{ij}), which is part of the source terms in cfc_alp
  subroutine gremhd_get_tilde_S(ixI^L,ixO^L,ps_in,tilde_S)
    use mod_global_parameters
    use mod_geometry, only: get_gamma_ij_hat

    integer, intent(in)             :: ixI^L, ixO^L
    type(state), intent(in)         :: ps_in
    double precision, intent(out)   :: tilde_S(ixI^S)

    double precision                :: E2(ixI^S), B2(ixI^S)
    double precision                :: v2(ixI^S)
    double precision                :: h(ixI^S)    ! enthalpy
    double precision                :: lfac2(ixI^S) ! Lorentz factor
    double precision                :: gamma(ixI^S,1:3,1:3)
    integer                         :: idir

    associate(prim=>ps_in%prim,x=>ps_in%x)
    call gremhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), lfac2=lfac2(ixI^S), &
                E2=E2(ixI^S), B2=B2(ixI^S), &
                h=h(ixI^S) )

    
    tilde_S(ixO^S) = prim(ixO^S, rho_) * h(ixO^S) * ( lfac2(ixO^S) - 1.0d0 ) & ! W2v2 + 1 = lfac2
                  + 3.0d0 * prim(ixO^S, press_) &
                  + 0.5d0 * ( E2(ixO^S) + B2(ixO^S) )
    tilde_S(ixO^S) = tilde_S(ixO^S) * prim(ixO^S, psi_)**6
    end associate
  end subroutine gremhd_get_tilde_S

  subroutine gremhd_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nprim)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

  end subroutine gremhd_get_dt

  !> Returns 0 in argument flag where values are ok
  subroutine gremhd_check_prim(ixI^L, ixO^L, prim, flag)
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
  end subroutine gremhd_check_prim

  !> calculate magnetic field from vector potential
  subroutine gremhd_b_from_vector_potential(ixIs^L, ixI^L, ixO^L, s)
    use mod_global_parameters
    use mod_usr_methods, only: usr_init_vector_potential
    use mod_geometry

    integer, intent(in)                :: ixIs^L, ixI^L, ixO^L
    type(state), intent(inout)         :: s

    integer                            :: ixC^L, hxC^L, hxO^L, idim1, idim2, idir, jdir

    double precision, allocatable      :: A(:^D&,:), dA(:^D&,:,:)
    double precision, allocatable      :: sqrt_gamma_hat(:^D&)
    double precision, allocatable      :: xC(:^D&,:),xCC(:^D&,:)
    double precision, allocatable      :: circ(:^D&,:)

    associate(prim=>s%prim, cons=>s%cons, x=>s%x, conss=>s%conss)

    if (stagger_grid) then
       allocate(A(ixIs^S,1:3))
       allocate(xC(ixIs^S,1:ndim), xCC(ixIs^S,1:ndim))
       allocate(circ(ixIs^S,1:ndim))
       ! extend one layer of cell center locations in xCC
       xCC=0.d0
       xCC(ixI^S,1:ndim)=x(ixI^S,1:ndim)
       {
       xCC(ixIsmin^D^D%ixI^S,1:ndim) = x(ixImin^D^D%ixI^S,1:ndim)
       xCC(ixIsmin^D^D%ixIs^S,^D) = x({ixImin^DD,},^D) - block%dx({ixImin^DD,},^D)
       \}
       {^IFTHREED
       xCC(ixImin1:ixImax1,ixIsmin2,ixIsmin3,1)=x(ixImin1:ixImax1,ixImin2,ixImin3,1)
       xCC(ixIsmin1,ixImin2:ixImax2,ixIsmin3,2)=x(ixImin1,ixImin2:ixImax2,ixImin3,2)
       xCC(ixIsmin1,ixIsmin2,ixImin3:ixImax3,3)=x(ixImin1,ixImin2,ixImin3:ixImax3,3)
       }
   
       A=0.0d0
       do idir=7-2*ndim,3
         ixCmax^D=ixOmax^D;
         ixCmin^D=ixOmin^D-1+kr(idir,^D);
         do idim1=1,ndim
           ! Get edge coordinates
           if (idim1/=idir) then
             xC(ixC^S,idim1)=xCC(ixC^S,idim1) + 0.5d0 * s%dx(ixC^S,idim1)
           else
             xC(ixC^S,idim1)=xCC(ixC^S,idim1)
           end if
         end do
         ! Initialise vector potential at the edge
         call usr_init_vector_potential(ixIs^L, ixC^L, xC, A(ixIs^S,idir), idir)
         A(ixC^S,idir) = A(ixC^S,idir) * s%dsC(ixC^S,idir)
       end do
   
       ! Set NaN to zero (can happen e.g. on axis):
       where( isnan(A(ixG^T,1:ndir)) )
          A(ixG^T,1:ndir)=0.0d0
       end where
   
       ! Take the curl of the vector potential 
       circ=0.0d0
       ! Calculate circulation on each face
       do idim1=1,ndim ! Coordinate perpendicular to face 
         ixCmax^D=ixOmax^D;
         ixCmin^D=ixOmin^D-kr(idim1,^D);
         do idim2=1,ndim
           do idir=7-2*ndim,3 ! Direction of line integral
             if(lvc(idim1,idim2,idir)==0) cycle
             ! Assemble indices
             hxC^L=ixC^L-kr(idim2,^D);
             ! Add line integrals in direction idir
             circ(ixC^S,idim1) = circ(ixC^S,idim1) &
                              + lvc(idim1,idim2,idir) * &
                               ( A(ixC^S,idir) - A(hxC^S,idir) )
           end do
         end do
         ! Divide by the area of the face to get B
         where ( s%surfaceC(ixC^S,idim1) < tiny(0.0d0) )
           circ(ixC^S,idim1) = 0.0d0
         elsewhere
           circ(ixC^S,idim1) = circ(ixC^S,idim1) / s%surfaceC(ixC^S,idim1)
         end where
         ! store in Bconss
         conss(ixC^S,idim1) = circ(ixC^S,idim1)
       end do

       ! update cell center magnetic field from Bconss
       call phys_face_to_center(ixO^L, s)

       deallocate(A)
       deallocate(xC, xCC)
       deallocate(circ)
    else
    ! if not stagger grid

       allocate(A(ixI^S,1:3), dA(ixI^S,1:3,1:3))
       allocate(sqrt_gamma_hat(ixI^S))

       ! Initialise vector potential at the edge
       A=0.0d0
       do idir=1, ndir
         call usr_init_vector_potential(ixI^L, ixI^L, x, A(ixI^S, idir), idir)
       end do

       dA = 0.0d0
       do idir=1, ndim
          do jdir=1, ndir
             call partial_d( A(ixI^S, jdir) ,ixI^L, ixO^L, idir, dA(ixI^S,jdir,idir) )
          end do
       end do

       call get_sqrt_gamma_hat(x, ixI^L, ixO^L, sqrt_gamma_hat(ixI^S))
       do idir=1, ndim
          idim1=mod(idir,ndir)+1   ! next direction
          idim2=mod(idir+1,ndir)+1 ! next next direction
          cons(ixO^S, Bcons(idir)) = ( dA(ixO^S,idim2,idim1) -  dA(ixO^S,idim1,idim2) ) / sqrt_gamma_hat(ixO^S)
          prim(ixO^S, Bvec(idir)) = cons(ixO^S, Bcons(idir)) / prim(ixO^S, psi_)**6
       end do

       deallocate(A, dA)
       deallocate(sqrt_gamma_hat)

    end if

    end associate

  end subroutine gremhd_b_from_vector_potential

end module mod_gremhd_phys
