module mod_gremhd_phys_implicit_update
  use mod_physics
  use mod_gremhd_phys_parameters

  implicit none
  private
  
  logical :: debug_flag = .False.

  ! Public methods
  public :: gremhd_phys_implicit_update_init

contains

  !> Initialize the module
  subroutine gremhd_phys_implicit_update_init()
    phys_implicit_update     => gremhd_implicit_update
    phys_evaluate_implicit   => gremhd_evaluate_implicit
  end subroutine gremhd_phys_implicit_update_init

  !> Evaluate Implicit part in place, i.e. psa==>F_im(psa)
  subroutine gremhd_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)   !< Compute implicit part from this state and update it
    double precision, intent(in) :: qtC      !< Both states psa and psb at this time level

    integer :: iigrid, igrid

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! primitive variables W_vel are needed in stiff source term
       call phys_to_primitive(ixG^LL,ixM^LL,psa(igrid)%cons(ixG^T,1:ncons),&
                        psa(igrid)%prim(ixG^T,1:nprim),psa(igrid)%x(ixG^T,1:ndim))
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call gremhd_get_stiff_source(qtC,ixG^LL,ixM^LL, &
        psa(igrid)%cons(ixG^T,1:ncons),psa(igrid)%prim(ixG^T,1:nprim),psa(igrid)%x(ixG^T,1:ndim))
    end do
    !$OMP END PARALLEL DO

  end subroutine gremhd_evaluate_implicit

  !> Implicit solve of psa=psb+dtfactor*dt*F_im(psa)
  subroutine gremhd_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)   !< Compute implicit part from this state and update it
    type(state), target :: psb(max_blocks)   !< Will be unchanged, as on entry
    double precision, intent(in) :: qdt      !< overall time step dt
    double precision, intent(in) :: qtC      !< Both states psa and psb at this time level
    double precision, intent(in) :: dtfactor !< Advance psa=psb+dtfactor*qdt*F_im(psa)

    integer :: iigrid, igrid

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       psa(igrid)%cons = psb(igrid)%cons
       if (stagger_grid) psa(igrid)%conss = psb(igrid)%conss
       call gremhd_implicit_update_grid(dtfactor,qdt,qtC,ixG^LL,ixM^LL, &
        psa(igrid)%cons(ixG^T,1:ncons),psa(igrid)%prim(ixG^T,1:nprim),psa(igrid)%x(ixG^T,1:ndim))
    end do
    !$OMP END PARALLEL DO
  end subroutine gremhd_implicit_update

  subroutine gremhd_get_stiff_source(qtC,ixI^L, ixO^L, cons, prim, x)
    use mod_global_parameters
    use mod_usr_methods
    double precision, intent(in)    :: qtC      !< Both states psa and psb at this time level
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
    double precision, intent(in)    :: prim(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision                :: lfac(ixI^S) ! Lorentz factor
    double precision                :: psi6(ixI^S) ! conformal factor **6
    double precision                :: eta(ixI^S) ! resistivity
    double precision                :: xi(ixI^S)  ! proper dynamo action term
    double precision                :: qE_dot_v(ixI^S), qB_dot_v(ixI^S) 
    double precision                :: gamma(ixI^S,1:3,1:3)
    double precision                :: sqrt_gamma(ixI^S)
    double precision                :: qJi(ixI^S,1:3) ! conformally transformed Current, but only the stiff part
    integer                         :: idir, jdir, kdir

    call usr_get_resistivity(ixI^L, ixO^L, cons(ixI^S, 1:ncons), x(ixI^S, 1:ndim), eta(ixI^S))
    call gremhd_get_intermediate_variables(ixI^L, ixO^L, prim(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), sqrt_gamma=sqrt_gamma(ixI^S), &
                lfac=lfac(ixI^S) )

    qE_dot_v(ixO^S) = 0.0d0
    do idir = 1, ndir
       qE_dot_v(ixO^S) = qE_dot_v(ixO^S) &
          + cons(ixO^S, Econs(idir)) * prim(ixO^S, W_vel(idir)) * gamma(ixO^S, idir, idir)
    end do
    qE_dot_v(ixO^S) = qE_dot_v(ixO^S) / lfac(ixO^S)
    qB_dot_v(ixO^S) = 0.0d0
    do idir = 1, ndir
       qB_dot_v(ixO^S) = qB_dot_v(ixO^S) &
          + cons(ixO^S, Bcons(idir)) * prim(ixO^S, W_vel(idir)) * gamma(ixO^S, idir, idir)
    end do
    qB_dot_v(ixO^S) = qB_dot_v(ixO^S) / lfac(ixO^S)

    ! dynamo part
    if (associated(usr_get_dynamo_coeff)) then
       call usr_get_dynamo_coeff(ixI^L, ixO^L, cons(ixI^S, 1:ncons), x(ixI^S, 1:ndim), xi(ixI^S))
       do idir = 1, ndir
          qJi(ixO^S, idir) = lfac(ixO^S) * cons(ixO^S, Bcons(idir)) - qB_dot_v(ixO^S) * prim(ixO^S, W_vel(idir))
          do jdir = 1, ndir
             do kdir = 1, ndir
                qJi(ixO^S, idir) = qJi(ixO^S, idir) &
                  - lvc(idir, jdir, kdir) / gamma(ixO^S, idir, idir) * sqrt_gamma(ixO^S) &
                    * prim(ixO^S, W_vel(jdir)) * cons(ixO^S, Econs(kdir))
             end do
          end do
       end do
       do idir = 1, ndir
          qJi(ixO^S, idir) = -xi(ixO^S) * qJi(ixO^S, idir)
       end do
    else
       qJi(ixO^S, :) = 0.0d0
    end if

    do idir = 1, ndir
       qJi(ixO^S, idir) = qJi(ixO^S, idir) &
             + lfac(ixO^S) * cons(ixO^S, Econs(idir)) - qE_dot_v(ixO^S) * prim(ixO^S, W_vel(idir))
       do jdir = 1, ndir
          do kdir = 1, ndir
             qJi(ixO^S, idir) = qJi(ixO^S, idir) &
               + lvc(idir, jdir, kdir) / gamma(ixO^S, idir, idir) * sqrt_gamma(ixO^S) &
                 * prim(ixO^S, W_vel(jdir)) * cons(ixO^S, Bcons(kdir))
          end do
       end do
    end do

    ! finally, we store the source terms directly in the cons
    cons(ixO^S, 1:ncons) = 0.0d0
    do idir = 1, ndir
       cons(ixO^S, Econs(idir)) = - prim(ixO^S, alp_) / eta(ixO^S) * qJi(ixO^S, idir)
    end do
  end subroutine gremhd_get_stiff_source

  subroutine gremhd_implicit_update_grid(dtfactor,qdt,qtC,ixI^L, ixO^L, cons, prim, x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    use mod_eos
    use mod_eos_idealgas, only: eos_gamma
    use mod_rootfinding
    double precision, intent(in)    :: qdt      !< overall time step dt
    double precision, intent(in)    :: qtC      !< Both states psa and psb at this time level
    double precision, intent(in)    :: dtfactor !< Advance psa=psb+dtfactor*qdt*F_im(psa)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
    double precision, intent(inout) :: prim(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision                :: one_over_eos_gamma_1 = 0.25d0  !< (eos_gamma-1)/eos_gamma  
    double precision                :: sigma_L = 0.0d0  !< this is defined as alp * dt
    double precision                :: eta_tilde = 0.0d0  !< this is defined as eta / ( alp * dt )

    integer                         :: idir, jdir, kdir
    double precision                :: psi6(ixI^S)
    integer                         :: ix^D 
    double precision                :: gamma(ixI^S,1:3,1:3) ! metric gamma_ij
    double precision                :: eta(ixI^S) ! resistivity
    double precision                :: xi(ixI^S)  ! proper dynamo action term
    double precision                :: sqrt_gamma(ixI^S)

    ! variables needed for the root founding
    integer                         :: error_code
    double precision                :: cons_tmp(ixI^S, 1:ncons)
    double precision                :: Ei(1:ndir)
    double precision                :: ui_hat(1:ndir)
    double precision                :: Jac(1:ndir, 1:ndir)
    double precision                :: lfac
   
    one_over_eos_gamma_1 = (eos_gamma - 1.0d0) / eos_gamma
    !one_over_eos_gamma_1 = (atmo_gamma - 1.0d0) / atmo_gamma

    call usr_get_resistivity(ixI^L, ixO^L, cons(ixI^S, 1:ncons), x(ixI^S, 1:ndim), eta(ixI^S))
    if(associated(usr_get_dynamo_coeff)) then
       call usr_get_dynamo_coeff(ixI^L, ixO^L, cons(ixI^S, 1:ncons), x(ixI^S, 1:ndim), xi(ixI^S))
    else 
       ! no dynamo effect
       xi(ixO^S) = 0.0d0
    end if

    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * prim(ixO^S, psi_)**4 
    end do
    psi6(ixO^S) = prim(ixO^S, psi_)**6 
    call get_sqrt_gamma_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, sqrt_gamma(ixI^S))
    sqrt_gamma(ixO^S) = psi6(ixO^S) * sqrt_gamma(ixO^S)
 
    ! conformal transformation back to normal conserved variables
    cons_tmp(ixO^S, D_) = cons(ixO^S, D_) / psi6(ixO^S)
    cons_tmp(ixO^S, tau_) = cons(ixO^S, tau_) / psi6(ixO^S)
    do idir = 1, ndir
       cons_tmp(ixO^S, mom(idir)) = cons(ixO^S, mom(idir)) / psi6(ixO^S)
       cons_tmp(ixO^S, Econs(idir)) = cons(ixO^S, Econs(idir)) / psi6(ixO^S)
       cons_tmp(ixO^S, Bcons(idir)) = cons(ixO^S, Bcons(idir)) / psi6(ixO^S)
    end do

    {do ix^D = ixO^LIM^D \}
      ! get sigma_L
      sigma_L = dtfactor * qdt * prim(ix^D, alp_) ! sigma_L = alpha * delta_t
      eta_tilde = eta(ix^D) / sigma_L
      if ( (cons_tmp(ix^D, D_) <= small_D) .or. &
          (cons_tmp(ix^D, tau_) <= small_tau) .or. &
          (eta(ix^D) >= 1.0d15) ) then 
         ! note: EM part is contained in this tau, however this relation should be held as well.
         ! atmosphere handling 
         ! skip the primitive recovery
         ui_hat(1:ndir) = 0.0d0
      else
         ! initial guess of ui
         do idir = 1, ndir
            ui_hat(idir) = prim(ix^D, W_vel(idir))
         end do
         call rootfinding_global_multid_newton_raphson(ui_hat, tol_im, iter_max_im, &
            error_code, all_func_of_u, fmin_of_f)
         !call rootfinding_multid_newton_raphson(ui_hat, tol_im, iter_max_im, &
         !     error_code, all_func_of_u)
         ! check the solution
         select case (error_code)
         !case (0) ! root is found
         !   write(*,*) "z= ", z
         case (-1) ! nothing happened
            call mpistop("have you ever attemp to find the root in con2prim?")
         case (1) ! failed to find the root
            if (debug_flag) then
               write(*,*) 'location x =',x(ix^D, 1:ndim)
               write(*,*) 'rho =', prim(ix^D, rho_), small_rho
               write(*,*) 'D = ', cons(ix^D, D_), cons_tmp(ix^D, D_), small_rho_thr
               write(*,*) 'D_atmo/D = ', small_D/cons(ix^D, D_)
               write(*,*) 'Wv = ', prim(ix^D, W_vel(:))
               write(*,*) 'ui = ', ui_hat(1:ndir)
               !write(*,*) 'f(Wv) = ', func_of_u(prim(ix^D, W_vel(:)))
               !write(*,*) 'f = ', func_of_u(ui_hat)
               write(*,*) 'eta =',eta(ix^D), eta_tilde
               write(*,*) 'sigma =',1.0d0/eta(ix^D)
               write(*,*) 'cons = ', cons(ix^D, 1:ncons)
               call mpistop("Error: Failed to find the root in gremhd_implicit_update")
            end if
            ! treat it as atmosphere
            ui_hat = 0.0d0
         case (2) ! z is NaN
            call mpistop("Error: the root is NaN in gremhd_implicit_update")
         case (3) ! reached local minimum
            if (debug_flag) then
               write(*,*) 'location x =',x(ix^D, 1:ndim)
               write(*,*) 'rho =', prim(ix^D, rho_), small_rho
               write(*,*) 'D = ', cons(ix^D, D_), cons_tmp(ix^D, D_), small_rho_thr
               write(*,*) 'D_atmo/D = ', small_D/cons(ix^D, D_)
               write(*,*) 'Wv = ', prim(ix^D, W_vel(:))
               write(*,*) 'ui = ', ui_hat(1:ndir)
               !write(*,*) 'f(Wv) = ', func_of_u(prim(ix^D, W_vel(:)))
               !write(*,*) 'f = ', func_of_u(ui_hat)
               write(*,*) 'eta =',eta(ix^D), eta_tilde
               write(*,*) 'sigma =',1.0d0/eta(ix^D)
               write(*,*) 'cons = ', cons(ix^D, 1:ncons)
               call mpistop("local minimum in gremhd_implicit_update")
            end if
         case (4) ! singular jacobian
            call mpistop("Error: singular jacobian in gremhd_implicit_update")
         case (5) ! slope is larger or equals to zero
            call mpistop("Error: slope is larger or equals to zero in gremhd_implicit_update")
         end select
      end if

      ! update conserved E field from the roots
      call get_vars(ui_hat, W_out=lfac, E_out=Ei)
      ! check if the resulting rho is lower then atmo
      if ( cons_tmp(ix^D, D_) / lfac <= small_rho_thr ) then
         ui_hat = 0.0d0 ! reset ui
         ! recalculate E field
         call get_vars(ui_hat, E_out=Ei)
      end if

      do idir = 1, ndir
         cons_tmp(ix^D, Econs(idir)) = Ei(idir)
         ! save the solution into prim, as an initial guess for other con2prim if needed
         prim(ix^D, W_vel(idir)) = ui_hat(idir)
      end do
    {enddo^D&\}

    do idir = 1, ndir
       cons(ixO^S, Econs(idir)) = cons_tmp(ixO^S, Econs(idir)) * psi6(ixO^S) 
    end do

    contains

       !> func f(ui) for finding the root
       subroutine get_vars(ui_in, W_out, E_out, f_out, fjac_out)
          implicit none
          double precision, intent(in)    :: ui_in(1:ndir)
          double precision, intent(out), optional   :: W_out
          double precision, intent(out), optional   :: E_out(1:ndir)
          double precision, intent(out), optional   :: f_out(1:ndir)
          double precision, intent(out), optional   :: fjac_out(1:ndir, 1:ndir)

          double precision                :: ui(1:ndir)

          double precision                :: E_dot_u, B_dot_u
          double precision                :: coeff_A(0:5)
          double precision                :: coeff_A_dot(0:5)

          integer                         :: idir,jdir,kdir,adir,bdir
          double precision                :: dhDdu(1:ndir)
          double precision                :: E_field(1:ndir)
          double precision                :: S_i_prime(1:ndir)
          double precision                :: dEdu(1:ndir,1:ndir)
          double precision                :: fluid_energy
          double precision                :: W_hat, hD
          
          W_hat = 1.0d0
          do idir = 1, ndir
             W_hat = W_hat + ui_in(idir)**2 * gamma(ix^D, idir, idir)
          end do
          W_hat = dsqrt( W_hat )
          if ( present(W_out) ) W_out = W_hat
    
          ui(:) = ui_in(:)

          ! get E field from u
          coeff_A(0) = W_hat + eta_tilde + xi(ix^D)**2 * ( W_hat**2 - 1.0d0 ) / ( W_hat + eta_tilde ) 
          coeff_A(1) = ( 1.0d0 + xi(ix^D)**2 * W_hat / ( W_hat + eta_tilde ) ) * eta_tilde / ( 1.0d0 + W_hat * eta_tilde )
          coeff_A(2) = - xi(ix^D) * eta_tilde / ( W_hat + eta_tilde )
          coeff_A(3) = xi(ix^D) * ( 1.0d0 + eta_tilde * W_hat ) / ( W_hat + eta_tilde )
          coeff_A(4) = xi(ix^D) * ( 1.0d0 - eta_tilde**2 + xi(ix^D)**2 ) / ( W_hat + eta_tilde ) / ( 1.0d0 + W_hat * eta_tilde )
          coeff_A(5) = - 1.0d0 - ( xi(ix^D)**2 * W_hat ) / ( W_hat + eta_tilde )

          E_dot_u = 0.0d0
          B_dot_u = 0.0d0
          do idir = 1, ndir
             E_dot_u = E_dot_u + cons_tmp(ix^D, Econs(idir)) * ui(idir) * gamma(ix^D, idir, idir)
             B_dot_u = B_dot_u + cons_tmp(ix^D, Bcons(idir)) * ui(idir) * gamma(ix^D, idir, idir)
          end do
          
          do idir = 1, ndir
             E_field(idir) = eta_tilde * cons_tmp(ix^D, Econs(idir)) &
                        + coeff_A(1) * E_dot_u * ui(idir) &
                        + coeff_A(3) * cons_tmp(ix^D, Bcons(idir)) &
                        + coeff_A(4) * B_dot_u * ui(idir) 

             jdir=mod(idir,ndir)+1   ! 'Next' direction
             kdir=mod(idir+1,ndir)+1 ! 'Next Next' direction

             ! gamma^{ii} (\eta_{ijk} u^j B^k) terms
             E_field(idir) = E_field(idir) &
                     + sqrt_gamma(ix^D) / gamma(ix^D, idir, idir) * ( &
                       ui(jdir) * ( coeff_A(2) * cons_tmp(ix^D, Econs(kdir)) + coeff_A(5) * cons_tmp(ix^D, Bcons(kdir)) ) & 
                       - ui(kdir) * ( coeff_A(2) * cons_tmp(ix^D, Econs(jdir)) + coeff_A(5) * cons_tmp(ix^D, Bcons(jdir)) ) & 
                       )
!             do jdir = 1, ndir
!                do kdir = 1, ndir
!                   E_field(idir) = E_field(idir) &
!                     + lvc(idir,jdir,kdir) / sqrt_gamma(ix^D) * ui(jdir) * gamma(ix^D, jdir, jdir) * gamma(ix^D, kdir, kdir) &
!                     * ( coeff_A(2) * cons_tmp(ix^D, Econs(kdir)) + coeff_A(5) * cons_tmp(ix^D, Bcons(kdir)) ) 
!                end do
!             end do
          end do
          E_field = E_field / coeff_A(0)
          if (present(E_out)) E_out = E_field

          if (present(f_out).or.present(fjac_out)) then
             fluid_energy = 0.0d0
             ! get EM energy first
             do idir = 1, ndir
                fluid_energy = fluid_energy + gamma(ix^D, idir, idir) * ( E_field(idir)**2 + cons_tmp(ix^D, Bcons(idir))**2 )
             end do
             ! fluid_energy = total energy (tau+D) - EM energy
             fluid_energy = cons_tmp(ix^D, tau_) - 0.5d0 * fluid_energy + cons_tmp(ix^D, D_)
   
             ! modified enthalpy
             hD = ( W_hat * fluid_energy - one_over_eos_gamma_1 * cons_tmp(ix^D, D_) ) &
                  / ( W_hat**2 - one_over_eos_gamma_1 ) 

             ! S_i_prime = S_i - EM part
             do idir = 1, ndir
                S_i_prime(idir) = cons_tmp(ix^D, mom(idir))
                jdir=mod(idir,ndir)+1   ! 'Next' direction
                kdir=mod(idir+1,ndir)+1 ! 'Next Next' direction
                S_i_prime(idir) = S_i_prime(idir) - sqrt_gamma(ix^D) * &
                      ( E_field(jdir) * cons_tmp(ix^D, Bcons(kdir)) & 
                       -E_field(kdir) * cons_tmp(ix^D, Bcons(jdir)) )
!                do jdir = 1, ndir
!                   do kdir = 1, ndir
!                      S_i_prime(idir) = S_i_prime(idir) & 
!                        - sqrt_gamma(ix^D) * lvc(idir,jdir,kdir) * E_field(jdir) * cons_tmp(ix^D, Bcons(kdir))
!                   end do
!                end do
             end do
   
             if (present(f_out)) then
                if (rescale_f) then
                   do idir = 1, ndir
                      f_out(idir) = gamma(ix^D, idir, idir) * ui(idir) - S_i_prime(idir) / hD
                   end do
                else
                   do idir = 1, ndir
                      f_out(idir) = gamma(ix^D, idir, idir) * ui(idir) * hD - S_i_prime(idir)
                   end do
                end if
             end if
             
             if (present(fjac_out)) then 
                ! get dEdu
                coeff_A_dot(0) = 1.0d0 &
                    + xi(ix^D)**2 * ( 1.0d0 + W_hat**2 + 2.0d0 * eta_tilde * W_hat ) &
                      / ( W_hat + eta_tilde )**2
                coeff_A_dot(1) = - eta_tilde**2 / ( 1.0d0 + W_hat * eta_tilde )**2 &
                                * ( 1.0d0 + xi(ix^D)**2 * ( W_hat**2 - 1.0d0 ) / ( W_hat + eta_tilde )**2 )
                coeff_A_dot(2) = xi(ix^D) * eta_tilde / ( W_hat + eta_tilde )**2
                coeff_A_dot(3) = xi(ix^D) * ( eta_tilde**2 - 1.0d0 ) / ( W_hat + eta_tilde )**2
                coeff_A_dot(4) = - xi(ix^D) * ( 1.0d0 - eta_tilde**2 + xi(ix^D)**2 ) &
                                * ( 1.0d0 + eta_tilde**2 + 2.0d0 * eta_tilde * W_hat) &
                               / ( W_hat + eta_tilde )**2 / ( 1.0d0 + W_hat * eta_tilde )**2 
                coeff_A_dot(5) = - xi(ix^D)**2 * eta_tilde / ( W_hat + eta_tilde )**2
   
                dEdu(:,:) = 0.0d0
                do adir = 1, ndir
                   do bdir = 1, ndir
      
                      dEdu(adir, bdir) = ( - coeff_A_dot(0) * E_field(adir) + coeff_A_dot(3) * cons_tmp(ix^D, Bcons(adir)) ) * ui(bdir) * gamma(ix^D, bdir, bdir) &
                        + W_hat * ui(adir) * ( coeff_A(1) * cons_tmp(ix^D, Econs(bdir)) + coeff_A(4) * cons_tmp(ix^D, Bcons(bdir)) ) * gamma(ix^D, bdir, bdir) &
                        + E_dot_u * ( coeff_A(1) * W_hat * kr(adir,bdir) + coeff_A_dot(1) * ui(adir) * ui(bdir) * gamma(ix^D, bdir, bdir) ) &
                        + B_dot_u * ( coeff_A(4) * W_hat * kr(adir,bdir) + coeff_A_dot(4) * ui(adir) * ui(bdir) * gamma(ix^D, bdir, bdir) )

                      jdir=mod(adir,ndir)+1   ! 'Next' direction
                      kdir=mod(adir+1,ndir)+1 ! 'Next Next' direction
                      dEdu(adir, bdir) = dEdu(adir, bdir) &
                         + sqrt_gamma(ix^D) / gamma(ix^D, adir, adir)  &
                          * ( ( cons_tmp(ix^D, Bcons(kdir)) * ( coeff_A(5) * W_hat * kr(jdir, bdir) &
                              + coeff_A_dot(5) * ui(jdir) * ui(bdir) * gamma(ix^D, bdir, bdir)) &
                             + cons_tmp(ix^D, Econs(kdir)) * ( coeff_A(2) * W_hat * kr(jdir, bdir) &
                              + coeff_A_dot(2) * ui(jdir) * ui(bdir) * gamma(ix^D, bdir, bdir)) ) &
                            - ( cons_tmp(ix^D, Bcons(jdir)) * ( coeff_A(5) * W_hat * kr(kdir, bdir) &
                              + coeff_A_dot(5) * ui(kdir) * ui(bdir) * gamma(ix^D, bdir, bdir)) &
                             + cons_tmp(ix^D, Econs(jdir)) * ( coeff_A(2) * W_hat * kr(kdir, bdir) &
                              + coeff_A_dot(2) * ui(kdir) * ui(bdir) * gamma(ix^D, bdir, bdir)) ) )

!                      do jdir = 1, ndir
!                         do kdir = 1, ndir
!                            dEdu(adir, bdir) = dEdu(adir, bdir) &
!                                   + lvc(adir,jdir,kdir) / sqrt_gamma(ix^D) * gamma(ix^D, kdir, kdir)  &
!                                   * ( cons_tmp(ix^D, Bcons(kdir)) * ( coeff_A(5) * W_hat * gamma(ix^D, jdir, bdir) &
!                                      + coeff_A_dot(5) * ui(jdir) * ui(bdir) * gamma(ix^D, jdir, jdir) * gamma(ix^D, bdir, bdir)) &
!                                     + cons_tmp(ix^D, Econs(kdir)) * ( coeff_A(2) * W_hat * gamma(ix^D, jdir, bdir) &
!                                      + coeff_A_dot(2) * ui(jdir) * ui(bdir) * gamma(ix^D, jdir, jdir) * gamma(ix^D, bdir, bdir)) )
!                         end do
!                      end do
      
                   end do
                end do
                dEdu(:,:) = dEdu(:,:) / coeff_A(0) / W_hat
       
                do jdir = 1, ndir
                   dhDdu(jdir) = ( fluid_energy / W_hat - 2.0d0 * hD ) * ui(jdir) * gamma(ix^D, jdir, jdir)
                   do idir = 1, ndir
                      dhDdu(jdir) = dhDdu(jdir) &
                            - W_hat * E_field(idir) * gamma(ix^D, idir, idir) * dEdu(idir, jdir)
                   end do
                end do
                dhDdu = dhDdu / ( W_hat**2 - one_over_eos_gamma_1 )
   
                ! finally, get Jac
                if (rescale_f) then
                   do idir = 1, ndir
                      do jdir = 1, ndir
      
                         fjac_out(idir, jdir) = gamma(ix^D, idir, jdir) &
                                   + S_i_prime(idir) / hD**2 * dhDdu(jdir)
   
                         adir=mod(idir,ndir)+1   ! 'Next' direction
                         bdir=mod(idir+1,ndir)+1 ! 'Next Next' direction
                         fjac_out(idir, jdir) = fjac_out(idir, jdir) &
                               + sqrt_gamma(ix^D) / hD * &
                                 ( dEdu(adir, jdir) * cons_tmp(ix^D, Bcons(bdir)) &
                                 - dEdu(bdir, jdir) * cons_tmp(ix^D, Bcons(adir)) )
!                         do adir = 1, ndir
!                            do bdir = 1, ndir
!                               fjac_out(idir, jdir) = fjac_out(idir, jdir) &
!                                   + sqrt_gamma(ix^D) / hD * lvc(idir,adir,bdir) * dEdu(adir, jdir) * cons_tmp(ix^D, Bcons(bdir))
!                            end do
!                         end do
   
                      end do
                   end do
                else
                   do idir = 1, ndir
                      do jdir = 1, ndir
      
                         fjac_out(idir, jdir) = gamma(ix^D, idir, jdir) * hD &
                                   + gamma(ix^D, idir, idir) * ui(idir) * dhDdu(jdir)

                         adir=mod(idir,ndir)+1   ! 'Next' direction
                         bdir=mod(idir+1,ndir)+1 ! 'Next Next' direction
                         fjac_out(idir, jdir) = fjac_out(idir, jdir) &
                               + sqrt_gamma(ix^D) * &
                                 ( dEdu(adir, jdir) * cons_tmp(ix^D, Bcons(bdir)) &
                                 - dEdu(bdir, jdir) * cons_tmp(ix^D, Bcons(adir)) )
   
!                         do adir = 1, ndir
!                            do bdir = 1, ndir
!                               fjac_out(idir, jdir) = fjac_out(idir, jdir) &
!                                   + sqrt_gamma(ix^D) * lvc(idir,adir,bdir) * dEdu(adir, jdir) * cons_tmp(ix^D, Bcons(bdir))
!                            end do
!                         end do
   
                      end do
                   end do
                end if


             end if
          end if
       end subroutine get_vars

       !> master function f(ui) and its Jacobian for finding the root
       subroutine all_func_of_u(ui,f_of_u,Jac)
          implicit none
          double precision, dimension(:), intent(in)     :: ui
          double precision, dimension(size(ui))          :: f_of_u
          double precision, dimension(size(ui), size(ui)):: Jac
          call get_vars(ui, f_out=f_of_u, fjac_out=Jac)
       end subroutine all_func_of_u

       function fmin_of_f(ui) result (f)
          implicit none
          double precision, dimension(:), intent(in)    :: ui
          double precision                              :: f
          double precision, dimension(1:size(ui))       :: fvec
          integer                                       :: idir
          call get_vars(ui, f_out=fvec)
          f = 0.0d0
          do idir=1, ndir
             f = f + fvec(idir)**2 / gamma(ix^D, idir, idir)
          end do
          !f = dot_product(fvec(:),fvec(:))
          f = 0.5d0 * f 
       end function fmin_of_f
  end subroutine gremhd_implicit_update_grid

end module mod_gremhd_phys_implicit_update
