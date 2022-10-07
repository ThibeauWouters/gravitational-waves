module mod_cfc_alp
  use mod_cfc_parameters
  implicit none
  private

  ! public methods
  public :: cfc_solve_alp

  contains

  subroutine cfc_solve_alp(A2)
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_forest
    use mod_physics
    use mod_eos, only: small_rho_thr

    double precision, intent(in)   :: A2(ixG^T,1:igridstail)
    integer                        :: id, iigrid, igrid
    double precision               :: res = huge(1.0d0)
    integer                        :: nc, lvl, mg_it
    type(tree_node), pointer       :: pnode
    double precision               :: tilde_S(ixG^T)
    double precision               :: src1(ixG^T)
    integer                        :: idir
    integer                        :: ix^L

    ix^L=ixM^LL^LADD1;

    mg%operator_type = mg_cfc_alp
    call mg_set_methods(mg)

    ! Same boundary conditions as psi, so we will not do it again here

    ! copy the data into MG solver
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! initial guess of (alp*psi-1) 
       mg%boxes(id)%cc({:,}, mg_iphi) = (ps(igrid)%prim(ix^S, alp_) * ps(igrid)%prim(ix^S, psi_) - 1.0d0)

       ! source term 1
       call phys_get_tilde_S(ixG^LL, ixM^LL, ps(igrid), tilde_S(ixG^T))
       src1(ixM^T) = - 2.0d0 * dpi / ps(igrid)%prim(ixM^T, psi_)**2 &
           * ( ps(igrid)%cons(ixM^T, tau_) + ps(igrid)%cons(ixM^T, D_) &
               + 2.0d0 * tilde_S(ixM^T) ) &
           - 7.0d0/8.0d0 / ps(igrid)%prim(ixM^T, psi_)**8 &
           * ( A2(ixM^T, iigrid) )

       mg%boxes(id)%cc({1:nc}, mg_itmp1) = src1(ixM^T)

       ! source term 2 (save psi), note that extra layer is needed for non-spherical cases
       mg%boxes(id)%cc({:,}, mg_itmp2) = ps(igrid)%prim(ix^S, psi_)
       ! No RHS 
       mg%boxes(id)%cc({1:nc}, mg_irhs) = 0.0d0
    end do

    do mg_it = 1, cfc_it_max
       call mg_fas_fmg(mg, (cfc_initial_guess .or. mg_it > 1), max_res=res)
       !write(*,*) 'solving alp:', mg_it, res
       if (res <= cfc_tol(1)) exit
       if ( cfc_test_alp .or. (mod(mg_it,cfc_print)==0) ) then
          if (mype==0) write(*,*) 'solving alp:', it, mg_it, res
       end if
    end do
   
    if (mg_it >= cfc_it_max) then
       if (mype==0) then
          write(*,*) "! Warning: Fail to converge alp."
          write(*,*) "! it=", it,"N = ", mg_it, " Res = ", res
          write(*,*) "! Maybe the cfc_tol is too small or cfc_it_max is too small"
       end if
    end if

    ! copy the data from MG solver
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)
       ! Include one layer of ghost cells on grid leaves
       ! the extra layer is needed, as we need to workout d/dr (alp/psi**6) !
       ps(igrid)%prim(ix^S, alp_) = ( mg%boxes(id)%cc({:,}, mg_iphi) + 1.0d0 ) &
                                  / ps(igrid)%prim(ix^S, psi_)
    end do

  end subroutine cfc_solve_alp

end module mod_cfc_alp
