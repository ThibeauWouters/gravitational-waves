module mod_grmhd_phys_divb_mg
  use mod_physics
  use mod_grmhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grmhd_phys_divb_mg_init
  public :: grmhd_clean_divb_multigrid
  public :: grmhd_initial_clean_divb_multigrid

contains

  !> Initialize the module
  subroutine grmhd_phys_divb_mg_init()
    use mod_global_parameters
    use_multigrid = .true.
  end subroutine grmhd_phys_divb_mg_init

  subroutine grmhd_initial_clean_divb_multigrid(iit,qt)
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_ghostcells_update, only: getbc
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer                      :: iigrid, igrid, idir
    double precision             :: divB_mg_tol_tmp
    if (.not. use_multigrid ) then
       call mg_setup_multigrid()
    end if
    divB_mg_tol_tmp = divB_mg_tol
    divB_mg_tol = -1.0d0 ! initially, we want divB as small as possible
    if (mype==0) then
      write(*, '(A,ES9.2,A)') ' Start cleaning B field'
    end if

    ! update Bcon based on the current psi and Bvec
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! up to ndim is enough
       do idir = 1, ndim
          ps(igrid)%cons(ixM^T,Bcons(idir)) = &
             ps(igrid)%prim(ixM^T,Bvec(idir)) * ps(igrid)%prim(ixM^T, psi_)**6 
       end do
    end do
    !$OMP END PARALLEL DO

    ! make sure psi are updated in the ghost cells
    call getbc(global_time,0.0d0,ps,1,nprim,phys_req_diagonal)

    call grmhd_clean_divb_mg()
    if (mype==0) then
      write(*, '(A,ES9.2,A)') ' Cleaned B field'
    end if
    divB_mg_tol = divB_mg_tol_tmp
  end subroutine grmhd_initial_clean_divb_multigrid

  subroutine grmhd_clean_divb_multigrid(iit,qt)
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    call grmhd_clean_divb_mg()
  end subroutine grmhd_clean_divb_multigrid

  ! elliptic cleaning on Bcons
  subroutine grmhd_clean_divb_mg()
    use mod_ghostcells_update, only: getbc
    use mod_forest
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_geometry

    integer                      :: iigrid, igrid, id
    integer                      :: n, nc, lvl, ix^L, idir, ixC^L
    type(tree_node), pointer     :: pnode
    double precision             :: tmp(ixG^T), grad(ixG^T, ndir)
    double precision             :: res = huge(1.0d0)
    double precision             :: res_old = huge(1.0d0)
    double precision, parameter  :: residual_reduction = 1.d-10
    double precision             :: max_divb

    ! update Bvec based on the current psi and Bcons
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! up to ndim is enough
       do idir = 1, ndim
          ps(igrid)%prim(ixM^T,Bvec(idir)) = &
             ps(igrid)%cons(ixM^T,Bcons(idir)) / ps(igrid)%prim(ixM^T, psi_)**6 
       end do
    end do
    !$OMP END PARALLEL DO
    call getbc(global_time,0.0d0,ps,Bvec(1),Bvec(ndim),phys_req_diagonal)

    mg%operator_type = mg_laplacian
    if (divB_mg_redblack) then
       mg%smoother_type = mg_smoother_gsrb
    else
       mg%smoother_type = mg_smoother_gs
    end if
    mg%n_cycle_up   = divB_mg_n_cycle(1)
    mg%n_cycle_down = divB_mg_n_cycle(2)
    call mg_set_methods(mg)

    ! Set boundary conditions
    do n = 1, 2*ndim
       idir = (n+1)/2
       select case (typeboundary(Bvec(idir), n))
       case ('symm')
          ! d/dx B = 0, take phi = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('asymm')
          ! B = 0, so grad(phi) = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('cont','noinflow')
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('special')
          ! Assume Dirichlet boundary conditions, derivative zero
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('periodic')
          ! Nothing to do here
       case default
          !print *, "divb_multigrid warning: unknown b.c.: ", &
          !     trim(typeboundary(Bvec(idir), n))
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       end select
    end do

    if (divB_4thorder) then
       ix^L=ixM^LL^LADD2;
    else
       ix^L=ixM^LL^LADD1;
    end if

    max_divb = 0.0d0

    ! Store divergence of B as right-hand side
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       call grmhd_get_divb(ixG^LL, ixM^LL, ps(igrid), tmp(ixG^T))

       if ( coordinate /= cartesian ) then
          ! rescale the rhs: rhs = rhs * r**2
          tmp(ixM^T) = tmp(ixM^T) * ps(igrid)%x(ixM^T, r_)**2
          {^NOONED
          if ( coordinate == spherical ) then
             ! rescale the rhs: rhs_final = rhs * r**2 sin_theta**2
             tmp(ixM^T) = tmp(ixM^T) * dsin(ps(igrid)%x(ixM^T, 2))**2 
          end if
          }
       end if

       ! divB as the RHS
       mg%boxes(id)%cc({1:nc}, mg_irhs) = tmp(ixM^T)

       !if (stagger_grid) &
       !   max_divb = max(max_divb, maxval(abs(tmp(ixM^T))))
       !write(*,*) max_divb
    end do

    ! fixme: maybe we could do laplician first and see if res < tol
    ! Solve laplacian(phi) = divB
    do n = 1, divB_mg_it_max
       res_old = res
       call mg_fas_fmg(mg, n>1, max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       !if (mype == 0) write(*,*) 'divB cleaning :', n, res
       if (res < divB_mg_tol) then
          exit
       end if
    end do

    if (mype == 0 .and. n > divB_mg_it_max .and. divB_mg_tol > 0.0d0) then
       print *, "divb_multigrid warning: not fully converged"
       print *, "current amplitude of divb: ", res
       print *, "multigrid smallest grid: ", mg%domain_size_lvl(:, mg%lowest_lvl)
       print *, "note: smallest grid ideally has <= 8 cells"
       print *, "multigrid dx/dy/dz ratio: ", mg%dr(:, 1)/mg%dr(1, 1)
       print *, "note: dx/dy/dz should be similar"
    end if

    ix^L=ixM^LL^LADD1;

    ! Correct the magnetic field
    do iigrid = 1, igridstail
       igrid = igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       ! Compute the gradient of phi, we need ghost zones as well
       tmp(ix^S) = mg%boxes(id)%cc({:,}, mg_iphi)

       if(stagger_grid) then
         do idir =1, ndim
           ixCmin^D=ixMlo^D-kr(idir,^D);
           ixCmax^D=ixMhi^D;
           call gradientx(tmp,ps(igrid)%x,ixG^LL,ixC^L,idir,grad(ixG^T,idir),.false.)
           ! Apply the correction B* = B - gradient(phi)
           ps(igrid)%conss(ixC^S,idir)=ps(igrid)%conss(ixC^S,idir)-grad(ixC^S,idir)
         end do
         ! change cell-center magnetic field Bcons
         call phys_face_to_center(ixM^LL,ps(igrid))
       else
         do idir = 1, ndim
            call gradient(tmp(ixG^T),ixG^LL,ixM^LL,idir,grad(ixG^T, idir))
            ! Apply the correction B* = B - gradient(phi)
            ps(igrid)%cons(ixM^T, Bcons(idir)) = ps(igrid)%cons(ixM^T, Bcons(idir)) &
                                - grad(ixM^T, idir)
            ps(igrid)%prim(ixM^T, Bvec(idir)) = ps(igrid)%prim(ixM^T, Bvec(idir)) &
                                - grad(ixM^T, idir) / ps(igrid)%prim(ixM^T, psi_)**6 
         end do
       end if
    end do
  end subroutine grmhd_clean_divb_mg

end module mod_grmhd_phys_divb_mg
