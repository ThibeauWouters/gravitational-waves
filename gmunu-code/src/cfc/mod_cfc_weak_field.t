module mod_cfc_weak_field
  use mod_cfc_parameters
  implicit none
  private

  ! public methods
  public :: cfc_solve_newtonian_phi

  contains

  subroutine cfc_phi_solver_init()
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_forest
    use mod_physics
    use mod_geometry

    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    double precision               :: src1(ixG^T)
    integer                        :: ix^L

    ix^L=ixM^LL^LADD1;
    mg%operator_type = mg_cfc_psi
    call mg_set_methods(mg)

    select case (coordinate)
    case (cartesian)
       do nb=1, 2*ndim
          select case (typeboundary(psi_, nb)) 
          case ('cont','noinflow')
             ! outer boundary, use Robin B. C.
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_robin
             mg%bc(nb, mg_iphi)%bc_value = 1.0d0
          case ('symm')
             ! inner boundary, use neumann
             mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(nb, mg_iphi)%bc_value = 0.0d0
          case default
             call mpistop("Error: The boundary conditions of psi/alp are not set correctly.")
          end select
       end do
    case (cylindrical)
       ! Neumann boundary condition for r=0
       mg%bc(1, mg_iphi)%bc_type = mg_bc_neumann
       mg%bc(1, mg_iphi)%bc_value = 0.0d0
       ! Robin boundary condition for r=rmax
       mg%bc(2, mg_iphi)%bc_type = mg_bc_robin
       mg%bc(2, mg_iphi)%bc_value = 1.0d0
       {^NOONED
       ! B.C. for z=zmin
       if( abs(xprobmin2) < smalldouble ) then
          mg%bc(3, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(3, mg_iphi)%bc_value = 0.0d0
       else
          mg%bc(3, mg_iphi)%bc_type = mg_bc_robin
          mg%bc(3, mg_iphi)%bc_value = 1.0d0
       end if
       ! B.C. for z = zmax
       mg%bc(4, mg_iphi)%bc_type = mg_bc_robin
       mg%bc(4, mg_iphi)%bc_value = 1.0d0
       }
       {^IFTHREED
       ! B.C. for phi
       if ( poleB(1,1) ) then
          ! pi-periodic conditions is applied at r=0
          ! nothing to do here
       else
          ! Neumann boundary condition for phi = phi_min and phi_max
          mg%bc(5:6, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(5:6, mg_iphi)%bc_value = 0.0d0
       end if
       }
    case (spherical)
       ! Neumann boundary condition for r=0
       mg%bc(1, mg_iphi)%bc_type = mg_bc_neumann
       mg%bc(1, mg_iphi)%bc_value = 0.0d0
       ! Robin boundary condition for r=rmax
       mg%bc(2, mg_iphi)%bc_type = mg_bc_robin
       mg%bc(2, mg_iphi)%bc_value = 1.0d0
       {^NOONED
       ! Neumann boundary condition for every boundary
       do nb=3, mg_num_neighbors
          mg%bc(nb, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(nb, mg_iphi)%bc_value = 0.0d0
       end do
       }
       {^IFTHREED
       ! B.C. for phi
       if ( poleB(1,2) .or. poleB(2,2) ) then
          ! pi-periodic conditions is applied at northpole or southpole
          ! nothing to do here
       else
          ! Neumann boundary condition for phi = phi_min and phi_max
          mg%bc(5:6, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(5:6, mg_iphi)%bc_value = 0.0d0
       end if
       }
    end select

    ! copy the data into MG solver, and set up the source terms
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! layer of ghost cells on grid leaves are not included here

       ! source term 
       mg%boxes(id)%cc({1:nc}, mg_irhs) = 4.0d0 * dpi * ps(igrid)%prim(ixM^T, rho_) 
       ! fixme: work out a more general case
       if ( coordinate /= cartesian ) then
          ! rescale the rhs: rhs = rhs * r**2
          mg%boxes(id)%cc({1:nc}, mg_irhs) = mg%boxes(id)%cc({1:nc}, mg_irhs) &
               * ps(igrid)%x(ixM^T, 1)**2 
          {^NOONED
          stop "error in cfc weak field: not ready yet"
          }
       end if
    end do
  end subroutine cfc_phi_solver_init

  subroutine cfc_solve_newtonian_phi()
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_forest
    use mod_physics

    integer                        :: id, iigrid, igrid
    double precision               :: res = huge(1.0d0)
    integer                        :: nc, lvl, mg_it, nb
    type(tree_node), pointer       :: pnode
    integer                        :: ix^L

    call cfc_phi_solver_init()

    do mg_it = 1, cfc_it_max
       call mg_fas_fmg(mg, (mg_it>1), max_res=res)
       !call mg_fas_vcycle(mg, max_res=res)
       !write(*,*) 'solving psi:', mg_it, res
       if (res <= cfc_tol(2)) exit
       if (mod(mg_it,cfc_print)==0) then
          if (mype==0) write(*,*) 'solving newtonian phi:', it, mg_it, res
       end if
    end do
   
    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge newtonian phi."
          write(*,*) "! it=", it,"N = ", mg_it, " Res = ", res
          write(*,*) "! Maybe the cfc_tol is too small or cfc_it_max is too small"
       end if
    end if

    ! copy the data from MG solver
    ix^L=ixM^LL^LADD1;
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       ps(igrid)%prim(ix^S, psi_) = 1.0d0 - 0.5d0 * mg%boxes(id)%cc({:,}, mg_iphi)
       ps(igrid)%prim(ix^S, alp_) = 1.0d0 + mg%boxes(id)%cc({:,}, mg_iphi)
       ps(igrid)%prim(ix^S, vecX(1:ndir)) = 0.0d0
       ps(igrid)%prim(ix^S, beta(1:ndir)) = 0.0d0
    end do

  end subroutine cfc_solve_newtonian_phi

end module mod_cfc_weak_field
