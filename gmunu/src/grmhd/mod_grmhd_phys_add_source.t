module mod_grmhd_phys_add_source
  use mod_physics
  use mod_grmhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: grmhd_phys_add_source_init

contains

  !> Initialize the module
  subroutine grmhd_phys_add_source_init()
    use mod_global_parameters
    use mod_geometry

    integer :: itr, idir

    if (coordinate /= cartesian) then
       ! we need geom source terms
       phys_add_source_geom     => grmhd_add_source_geom
    end if
    phys_add_source          => grmhd_add_source
  end subroutine grmhd_phys_add_source_init

  subroutine get_g_up_munu(g, alp, beta, gamma, ixI^L, ixO^L)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: alp(ixI^S), beta(ixI^S, 1:3)
    double precision, intent(in)    :: gamma(ixI^S, 1:3, 1:3)
    double precision, intent(out)   :: g(ixI^S, 0:3, 0:3)
    integer                         :: mu, nu
    g(ixO^S, 0:3,0:3) = 0.0d0
    g(ixO^S, 0,0) = - 1.0d0 / alp(ixO^S)**2
    do mu = 1, 3
       g(ixO^S, 0, mu) = beta(ixO^S, mu) / alp(ixO^S)**2
       g(ixO^S, mu, 0) = g(ixO^S, 0, mu)  
    end do
    ! this is actually \gamma^{ij}
    do mu = 1, 3
       g(ixO^S, mu, mu) = 1.0d0 / gamma(ixO^S, mu, mu)
    end do
    ! g^{ij} = \gamma^{ij} - \beta^i \beta^j / \alpha^2
    do mu = 1, 3
       do nu = 1, 3
          g(ixO^S, mu, nu) = g(ixO^S, mu, nu) - &
                             beta(ixO^S, mu) * beta(ixO^S, nu) / alp(ixO^S)**2
       end do
    end do
  end subroutine get_g_up_munu

  !> Add gravitational source terms to w
  subroutine grmhd_add_source(qdt, ixI^L, ixO^L, primCT, cons, x, qsourcesplit, active)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(in)    :: primCT(ixI^S, 1:nprim)
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active !< Needs to be set to true when active

    ! MHD variables needed when using GLM
    double precision                :: grad_Bphi(ixI^S,1:ndir)

    {^NOONED
    double precision                :: cot_theta(ixI^S)
    }
    double precision                :: alppsi6(ixI^S) 
    double precision                :: add_source(ixI^S,1:ncons)
    double precision                :: htot(ixI^S) ! enthalpy
    double precision                :: Ptot(ixI^S) ! enthalpy
    double precision                :: lfac(ixI^S) ! Lorentz factor
    ! metric variables
    double precision                :: alp(ixI^S)
    double precision                :: psi(ixI^S), psi4(ixI^S), psi6(ixI^S)
    double precision                :: betai(ixI^S,1:3)
    ! dervitaves of the metric variables
    double precision                :: dalp(ixI^S,1:3), dpsi(ixI^S,1:3)
    double precision                :: dbeta(ixI^S,1:3,1:3)

    ! 3-metric gamma_ij
    double precision                :: gamma(ixI^S,1:3,1:3)
    ! extrinsic curvature K_ij
    double precision                :: K_ij(ixI^S,1:3,1:3)
    ! 4-metric g^{\mu\nu}
    double precision                :: g(ixI^S,0:3,0:3)
    ! dervitaves of the 3-metric gamma_{jk},i
    double precision                :: Dgamma(ixI^S,1:3,1:3,1:3)
    ! Christoffel symbols of the reference 3-metric gamma_hat
    double precision                :: christoffel(ixI^S,1:3,1:3,1:3)
    ! covariant dervitaves of beta D_i beta^k
    double precision                :: D_beta(ixI^S,1:3,1:3)
    ! energy-momentum tensor T^{\mu\nu}
    double precision                :: T(ixI^S,0:3,0:3)
    ! 4-velocity u^{\mu}
    double precision                :: u(ixI^S,0:3)
    double precision                :: bmu(ixI^S,0:3)
    integer                         :: iw,idir,jdir,kdir

    !-----------------------------------------------------------------------
    ! Part 0: For the split source part
    !-----------------------------------------------------------------------
    ! NO split source terms at the moment
    if (qsourcesplit) return
    
    add_source = 0.0d0

    ! fixme: maybe use pointers?
    ! initialize the metric
    psi(ixI^S) = primCT(ixI^S,psi_)
    psi4(ixI^S) = psi(ixI^S)**4
    psi6(ixI^S) = psi(ixI^S)**6 !fixme: really need this?
    alp(ixI^S) = primCT(ixI^S,alp_)
    betai(ixI^S,1:3) = 0.0d0
    do idir = 1, ndir
       betai(ixI^S,idir) = primCT(ixI^S, beta(idir))
    end do
    ! alp * psi**6
    alppsi6(ixI^S) = alp(ixI^S) * psi(ixI^S)**6
    ! volume averaged Christoffel symbols
    Christoffel(ixI^S,1:3,1:3,1:3) = block%christoffel(ixI^S,1:3,1:3,1:3)

    call grmhd_get_intermediate_variables(ixI^L, ixO^L, primCT(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                bmu=bmu(ixI^S,0:3), &
                lfac=lfac(ixI^S), &
                htot=htot(ixI^S), &
                Ptot=Ptot(ixI^S) )
    
    !-----------------------------------------------------------------------
    ! Part 1: normal source terms 
    !-----------------------------------------------------------------------
   
    if ( evolve_EM .and. type_divb == divb_glm ) then
       add_source(ixO^S, Bphi_cons_) = - alppsi6(ixO^S) * primCT(ixO^S, Bphi_) * divB_glm_kappa
       do idir = 1, ndim
          call gradient( primCT(ixI^S, Bphi_),ixI^L,ixO^L,idir,grad_Bphi(ixI^S,idir) )
          add_source(ixO^S, Bcons(idir)) = - alppsi6(ixO^S) * grad_Bphi(ixO^S, idir) / gamma(ixO^S, idir, idir)
       end do

       ! add the source terms to Bphi and Bcons
       cons(ixO^S, Bphi_cons_) = cons(ixO^S, Bphi_cons_) + qdt * add_source(ixO^S, Bphi_cons_)
       do idir = 1, ndim
          cons(ixO^S, Bcons(idir)) =  cons(ixO^S, Bcons(idir)) + qdt * add_source(ixO^S, Bcons(idir))
       end do
    end if
    
    !-----------------------------------------------------------------------
    ! Part 2: gravitational source terms, only needed when use_GR = .true.
    !-----------------------------------------------------------------------
    if ( .not. use_GR ) return

    ! calculate derivitives of the metric variables
    dalp(ixI^S,1:3) = 0.0d0
    dbeta(ixI^S,1:3,1:3) = 0.0d0
    dpsi(ixI^S,1:3) = 0.0d0
    do idir = 1, ndim
       call partial_d( alp(ixI^S) ,ixI^L,ixO^L,idir,dalp(ixI^S,idir) )
       call partial_d( psi(ixI^S) ,ixI^L,ixO^L,idir,dpsi(ixI^S,idir) )
       dpsi(ixO^S,idir) = dpsi(ixO^S,idir) / psi(ixO^S)
       do jdir = 1, ndir
          call partial_d( betai(ixI^S,jdir) ,ixI^L,ixO^L,idir,dbeta(ixI^S,jdir,idir) )
       end do
    end do

    ! covariant derivative of beta: partial_i beta^k + Gamma^k_{ij} beta^j
    D_beta(ixO^S, 1:3, 1:3) = dbeta(ixO^S, 1:3, 1:3)
    if ( coordinate /= cartesian ) then
       do idir = 1, 3
          do kdir = 1, 3 
             do jdir = 1, 3 
                D_beta(ixO^S, kdir, idir) = D_beta(ixO^S, kdir, idir) &
                                          + Christoffel(ixO^S, kdir,idir,jdir) * betai(ixO^S,jdir)
             end do
          end do
       end do
    end if

    ! partial dervitaves of metric D_i gamma_jk
    Dgamma(ixI^S,1:3,1:3,1:3) = 0.0d0
    do kdir = 1, ndim
       do idir = 1, 3
          Dgamma(ixO^S,idir,idir,kdir) = 4.0d0 * gamma(ixO^S,idir,idir) * dpsi(ixO^S,kdir)
       end do
    end do

    ! Note that g(mu,nu) here is g^{\mu\nu}
    call get_g_up_munu(g(ixI^S,0:3,0:3), alp(ixI^S), betai(ixI^S,1:3), gamma(ixI^S,1:3,1:3), ixI^L, ixO^L)

    ! Calculate the 4-velocity
    u(ixO^S,0:3) = 0.0d0
    u(ixO^S,0) = lfac(ixO^S) / alp(ixO^S)
    do idir = 1, ndir
       u(ixO^S, idir) = primCT(ixO^S, W_vel(idir)) - lfac(ixO^S) * betai(ixO^S,idir) / alp(ixO^S)
    end do

    ! energy-momentum tensor T^{\mu\nu}
    do idir = 0,3
       do jdir = 0,3
          T(ixO^S,idir,jdir) = primCT(ixO^S, rho_) * htot(ixO^S) * u(ixO^S, idir) * u(ixO^S, jdir) &
                               + Ptot(ixO^S) * g(ixO^S,idir,jdir) &
                               - bmu(ixO^S, idir) * bmu(ixO^S, jdir)
       end do
    end do

    ! Now calculate the source terms for HD (mom, tau) in the compact form
    do idir = 1, ndir
       do jdir = 1, 3 
          do kdir = 1, 3 
             add_source(ixO^S, mom(idir)) = add_source(ixO^S, mom(idir)) &
                                       + 0.5d0 * Dgamma(ixO^S, jdir,kdir,idir) * ( &
                                        T(ixO^S,0,0)*betai(ixO^S,kdir)*betai(ixO^S,jdir) &
                                         + 2.0d0*T(ixO^S,0,jdir)*betai(ixO^S,kdir) + T(ixO^S, jdir,kdir) )
          end do
       end do
    end do

    do idir = 1, ndir
       add_source(ixO^S, mom(idir)) = add_source(ixO^S, mom(idir)) &
                                 - T(ixO^S,0,0)*alp(ixO^S)*dalp(ixO^S, idir) 
       ! T^{0j} * gamma_jk * D_i beta^k
       do kdir = 1, 3 
          do jdir = 1, 3 
             add_source(ixO^S, mom(idir)) = add_source(ixO^S, mom(idir)) &
                                    + T(ixO^S,0,jdir) * gamma(ixO^S,jdir,kdir) * D_beta(ixO^S,kdir,idir)
          end do
       end do
    end do

    ! To compute the tau source term ( add_source(:,:,:,tau_) ), we need to work out the extrinsic curvature K_{ij}
    ! fixme: make get_K_ij as a stand alone subroutine
    K_ij(ixO^S,1:3,1:3) = 0.0d0

    select case (coordinate)
    case (cartesian)
       K_ij(ixO^S,1,1) = gamma(ixO^S,1,1)/(3.0d0*alp(ixO^S)) &
                    * ( 2.0d0*dbeta(ixO^S,1,1) & 
                      {^NOONED - dbeta(ixO^S,2,2) } &
                      {^IFTHREED - dbeta(ixO^S,3,3) } )
   
       K_ij(ixO^S,2,2) = gamma(ixO^S,2,2)/(3.0d0*alp(ixO^S)) &
                    * ( - dbeta(ixO^S,1,1) &
                      {^NOONED + 2.0d0*dbeta(ixO^S,2,2) } &
                      {^IFTHREED - dbeta(ixO^S,3,3) } )
   
       K_ij(ixO^S,3,3) = gamma(ixO^S,3,3)/(3.0d0*alp(ixO^S)) &
                    * ( - dbeta(ixO^S,1,1) &
                      {^NOONED - dbeta(ixO^S,2,2) } &
                      {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
    case (cylindrical)
       K_ij(ixO^S,1,1) = gamma(ixO^S,1,1)/(3.0d0*alp(ixO^S)) &
                    * ( 2.0d0*dbeta(ixO^S,1,1) - betai(ixO^S,1)/x(ixO^S,r_) & 
                      {^NOONED - dbeta(ixO^S,2,2) } &
                      {^IFTHREED - dbeta(ixO^S,3,3) } )
   
       K_ij(ixO^S,2,2) = gamma(ixO^S,2,2)/(3.0d0*alp(ixO^S)) &
                    * ( - dbeta(ixO^S,1,1) - betai(ixO^S,1)/x(ixO^S,r_) &
                      {^NOONED + 2.0d0 * dbeta(ixO^S,2,2) } &
                      {^IFTHREED - dbeta(ixO^S,3,3) } )
   
       K_ij(ixO^S,3,3) = gamma(ixO^S,3,3)/(3.0d0*alp(ixO^S)) &
                    * ( - dbeta(ixO^S,1,1) + 2.0d0 * betai(ixO^S,1)/x(ixO^S,r_) &
                      {^NOONED - dbeta(ixO^S,2,2) } &
                      {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
    case (spherical)
       {^NOONED
       cot_theta(ixO^S) = dcos(x(ixO^S,theta_))/dsin(x(ixO^S,theta_))
       }

       K_ij(ixO^S,1,1) = gamma(ixO^S,1,1)/(3.0d0*alp(ixO^S)) &
                    * ( 2.0d0*dbeta(ixO^S,1,1) - 2.0d0*betai(ixO^S,1)/x(ixO^S,r_) & 
                      {^NOONED - dbeta(ixO^S,2,2) - betai(ixO^S,2)*cot_theta(ixO^S) } &
                      {^IFTHREED - dbeta(ixO^S,3,3) } )
   
       K_ij(ixO^S,2,2) = gamma(ixO^S,2,2)/(3.0d0*alp(ixO^S)) &
                    * ( -dbeta(ixO^S,1,1) + betai(ixO^S,1)/x(ixO^S,r_) &
                      {^NOONED + 2.0d0*dbeta(ixO^S,2,2) - betai(ixO^S,2)*cot_theta(ixO^S) } &
                      {^IFTHREED - dbeta(ixO^S,3,3) } )
   
       K_ij(ixO^S,3,3) = gamma(ixO^S,3,3)/(3.0d0*alp(ixO^S)) &
                    * ( -dbeta(ixO^S,1,1) + betai(ixO^S,1)/x(ixO^S,r_) &
                      {^NOONED - dbeta(ixO^S,2,2) + 2.0d0*betai(ixO^S,2)*cot_theta(ixO^S) } &
                      {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
    end select
    {^NOONED
    K_ij(ixO^S,1,2) = ( dbeta(ixO^S,1,2)*gamma(ixO^S,1,1) + dbeta(ixO^S,2,1)*gamma(ixO^S,2,2) ) / (2.0d0*alp(ixO^S))
    K_ij(ixO^S,1,3) = ( dbeta(ixO^S,3,1)*gamma(ixO^S,3,3) + dbeta(ixO^S,1,3)*gamma(ixO^S,1,1) ) / (2.0d0*alp(ixO^S))
    K_ij(ixO^S,2,3) = ( dbeta(ixO^S,3,2)*gamma(ixO^S,3,3) + dbeta(ixO^S,2,3)*gamma(ixO^S,2,2) ) / (2.0d0*alp(ixO^S))
    ! K_ji=K_ij
    do idir=1,2
       do jdir=idir+1,3
          K_ij(ixO^S,jdir,idir) = K_ij(ixO^S,idir,jdir)
       end do
    end do
    }

    ! Now work out the source term of tau
    do idir=1,3
       add_source(ixO^S,tau_) = add_source(ixO^S,tau_) &
                  + T(ixO^S,0,0) * ( -betai(ixO^S,idir)*dalp(ixO^S,idir) ) &
                  + T(ixO^S,0,idir) * ( -dalp(ixO^S,idir) ) 
    end do

    do idir=1,3
       do jdir=1,3
          add_source(ixO^S,tau_) = add_source(ixO^S,tau_) &
                    + T(ixO^S,0,idir)*( 2.0d0*betai(ixO^S,jdir)*K_ij(ixO^S,idir,jdir) ) &
                    + T(ixO^S,0,0) * ( K_ij(ixO^S,idir,jdir)*betai(ixO^S,idir)*betai(ixO^S,jdir) ) &
                    + T(ixO^S,idir,jdir) * K_ij(ixO^S,idir,jdir)
       end do
    end do

    if ( evolve_hydro ) then
       cons(ixO^S, tau_) = cons(ixO^S, tau_) + qdt * alppsi6(ixO^S) * add_source(ixO^S, tau_)
       do idir =1, ndir
         cons(ixO^S, mom(idir)) = cons(ixO^S, mom(idir)) + qdt * alppsi6(ixO^S) * add_source(ixO^S,mom(idir))
       end do
    end if ! end hydro source

    ! Now work out the source term of Bphi if glm is being used
    if ( evolve_EM .and. type_divb == divb_glm ) then
       add_source(ixO^S,Bphi_cons_) = 0.0d0 ! reset add source for Bphi
       do idir=1, ndir
          add_source(ixO^S, Bphi_cons_) = add_source(ixO^S, Bphi_cons_) &
                   + K_ij(ixO^S,idir,idir) / gamma(ixO^S,idir,idir) 
       end do
       add_source(ixO^S, Bphi_cons_) = alp(ixO^S) * primCT(ixO^S, Bphi_) * add_source(ixO^S,Bphi_cons_) 

       do idir=1, ndim ! be careful! 
          add_source(ixO^S,Bphi_cons_) = add_source(ixO^S,Bphi_cons_) &
                       + primCT(ixO^S, Bvec(idir)) * dalp(ixO^S, idir)
       end do

       ! add the source terms on cons
       cons(ixO^S, Bphi_cons_) = cons(ixO^S, Bphi_cons_) + qdt * psi6(ixO^S) * add_source(ixO^S, Bphi_cons_)
    end if

  end subroutine grmhd_add_source

  !> Add geometrical source terms to w
  subroutine grmhd_add_source_geom(qdt, ixI^L, ixO^L, consCT, primCT, cons, x)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(in)    :: consCT(ixI^S, 1:ncons), primCT(ixI^S, 1:nprim)
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)

    integer                         :: iw,idir,jdir,kdir, h1x^L{^NOONED, h2x^L}
    double precision                :: source(ixI^S)
    double precision                :: fluxCT(ixI^S,1:3,1:3)
    double precision                :: alppsi6(ixI^S)
    double precision                :: B_dot_v(ixI^S)
    double precision                :: v_hat(ixI^S,1:ndir)
    double precision                :: gamma(ixI^S,1:3,1:3) ! metric gamma_ij
    double precision                :: B_over_W(ixI^S, 1:ndir) 
    double precision                :: Ptot(ixI^S) ! total pressure
    double precision                :: lfac(ixI^S) ! Lorentz factor

    if ( .not. evolve_hydro ) return

    call grmhd_get_intermediate_variables(ixI^L, ixO^L, primCT(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                lfac=lfac(ixI^S), &
                v_hat=v_hat(ixI^S,1:ndir), &
                B_dot_v=B_dot_v(ixI^S), &
                Ptot=Ptot(ixI^S) )

    ! alp * psi**6
    alppsi6(ixO^S) = primCT(ixO^S,alp_) * primCT(ixO^S,psi_)**6

    ! B^i / W
    do idir = 1, ndir
       B_over_W(ixO^S, idir) = primCT(ixO^S, Bvec(idir)) / lfac(ixO^S)
    end do

    ! Only Momentum flux f^i_j is needed only
    fluxCT = 0.0d0
    do jdir = 1, ndir
       do idir = 1, ndir
         fluxCT(ixO^S, idir, jdir) = v_hat(ixO^S,idir) * consCT(ixO^S, mom(jdir)) &
                       - alppsi6(ixO^S) * ( B_over_W(ixO^S, jdir) + B_dot_v(ixO^S) * primCT(ixO^S, W_vel(jdir)) ) &
                       * gamma(ixO^S, jdir, jdir) * B_over_W(ixO^S, idir)
       end do
    end do
    ! and the pressure terms
    do jdir = 1, 3
       fluxCT(ixO^S, jdir,jdir) = fluxCT(ixO^S, jdir, jdir) + alppsi6(ixO^S) * Ptot(ixO^S)
    end do

    select case (coordinate)
    case (cylindrical)
       ! geo_source[r] = Gamma^phi_{r phi}*f^phi_phi
       source(ixO^S) = block%christoffel(ixO^S,3,1,3) * fluxCT(ixO^S,3,3)
       cons(ixO^S, mom(r_)) = cons(ixO^S, mom(r_)) + qdt * source(ixO^S)

    case (spherical)
       ! geo_source[r] = Gamma^2_{12}*f^2_2 + Gamma^3_{13}*f^3_3
       source(ixO^S) = block%christoffel(ixO^S,2,1,2) * fluxCT(ixO^S,2,2) & 
                      +block%christoffel(ixO^S,3,1,3) * fluxCT(ixO^S,3,3) 

       cons(ixO^S, mom(r_)) = cons(ixO^S, mom(r_)) + qdt * source(ixO^S)

       {^NOONED
       ! geo_source[theta] = Gamma^1_{22}*f^2_1 + Gamma^2_{21}*f^1_2 + Gamma^3_{23}*f^3_3
       source(ixO^S) = block%christoffel(ixO^S,1,2,2) * fluxCT(ixO^S,2,1) & 
                        +block%christoffel(ixO^S,2,2,1) * fluxCT(ixO^S,1,2) &
                        +block%christoffel(ixO^S,3,2,3) * fluxCT(ixO^S,3,3) 

       cons(ixO^S, mom(2)) = cons(ixO^S, mom(2)) + qdt * source(ixO^S)
       ! geo_source[phi] = 0
       }
    end select
  end subroutine grmhd_add_source_geom

end module mod_grmhd_phys_add_source
