module mod_gremhd_phys_add_source
  use mod_physics
  use mod_gremhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: gremhd_phys_add_source_init

contains

  !> Initialize the module
  subroutine gremhd_phys_add_source_init()
    use mod_global_parameters
    use mod_geometry

    integer :: itr, idir

    if (coordinate /= cartesian) then
       ! we need geom source terms
       phys_add_source_geom     => gremhd_add_source_geom
    end if
    phys_add_source          => gremhd_add_source
  end subroutine gremhd_phys_add_source_init

  subroutine get_g_munu(g_munu, alp, beta, gamma, ixI^L, ixO^L)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: alp(ixI^S), beta(ixI^S, 1:3)
    double precision, intent(in)    :: gamma(ixI^S, 1:3, 1:3)
    double precision, intent(out)   :: g_munu(ixI^S, 0:3, 0:3)
    integer                         :: mu, nu
    g_munu(ixO^S, 0:3,0:3) = 0.0d0
    g_munu(ixO^S, 0,0) = - alp(ixO^S)**2
    do mu = 1, 3
       g_munu(ixO^S, 0, 0) = g_munu(ixO^S, 0, 0) + beta(ixO^S, mu)**2 * gamma(ixO^S, mu, mu)
    end do
    do mu = 1, 3
       g_munu(ixO^S, 0, mu) = beta(ixO^S, mu) * gamma(ixO^S, mu, mu)
       g_munu(ixO^S, mu, 0) = g_munu(ixO^S, 0, mu)  
    end do
    ! this is actually \gamma_{ij}
    do mu = 1, 3
       g_munu(ixO^S, mu, mu) = gamma(ixO^S, mu, mu)
    end do
  end subroutine get_g_munu

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
  subroutine gremhd_add_source(qdt, ixI^L, ixO^L, primCT, cons, x, qsourcesplit, active)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(in)    :: primCT(ixI^S, 1:nprim)
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active !< Needs to be set to true when active

    ! MHD variables needed when using GLM
    double precision                :: grad_phi(ixI^S,1:ndir)

    {^NOONED
    double precision                :: cot_theta(ixI^S)
    }
    double precision                :: alppsi6(ixI^S) 
    double precision                :: add_source(ixI^S,1:ncons)
    double precision                :: h(ixI^S) ! enthalpy
    double precision                :: lfac(ixI^S) ! Lorentz factor
    ! metric variables
    double precision                :: alp(ixI^S)
    double precision                :: psi(ixI^S), psi4(ixI^S), psi6(ixI^S)
    double precision                :: betai(ixI^S,1:3)
    ! dervitaves of the metric variables
    double precision                :: dalp(ixI^S,1:3), dpsi(ixI^S,1:3)
    double precision                :: dbeta(ixI^S,1:3,1:3)

    ! sqrt(gamma_ij)
    double precision                :: sqrt_gamma(ixI^S)
    ! 3-metric gamma_ij
    double precision                :: gamma(ixI^S,1:3,1:3)
    ! extrinsic curvature K_ij
    double precision                :: K_ij(ixI^S,1:3,1:3)
    ! 4-metric g_{\mu\nu}
    double precision                :: g_munu(ixI^S,0:3,0:3)
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
    ! Maxwell tensor F^{\mu\nu}
    double precision                :: F(ixI^S,0:3,0:3)
    ! square of Maxwell tensor F^{\mu\nu}F_{\mu\nu}
    double precision                :: F2(ixI^S)
    ! 4-velocity n^{\mu} of the Eulerian observer
    double precision                :: nmu(ixI^S,0:ndir), n_mu(ixI^S,0:ndir)
    ! 4-velocity u^{\mu}
    double precision                :: u(ixI^S,0:3)
    ! 4-EM fields E_mu and B_mu
    double precision                :: Bmu(ixI^S,0:3), Emu(ixI^S,0:3)
    double precision                :: B_mu(ixI^S,0:3), E_mu(ixI^S,0:3)
    double precision                :: B2(ixI^S), E2(ixI^S)

    double precision                :: v_hat(ixI^S,1:3)
    double precision                :: qrho_e(ixI^S) ! conformally rescaled charge density
    integer                         :: idir,jdir,kdir,ldir

    !-----------------------------------------------------------------------
    ! Part 0: For the split source part
    !-----------------------------------------------------------------------
    ! NO split source terms at the moment
    if (qsourcesplit) return
    
    add_source(ixO^S, 1:ncons) = 0.0d0

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
    
    call gremhd_get_intermediate_variables(ixI^L, ixO^L, primCT(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                sqrt_gamma=sqrt_gamma(ixI^S), gamma=gamma(ixI^S,1:3,1:3), &
                lfac=lfac(ixI^S), v_hat=v_hat(ixI^S,1:ndir), &
                E2=E2(ixI^S), B2=B2(ixI^S), qrho_e=qrho_e(ixI^S), &
                h=h(ixI^S))

    !-----------------------------------------------------------------------
    ! Part 1: normal source terms 
    !-----------------------------------------------------------------------

    if ( evolve_EM ) then
   
       if ( type_divb == glm ) then
          ! add kappa source terms to Bphi
          add_source(ixO^S, Bphi_cons_) = - alppsi6(ixO^S) * primCT(ixO^S, Bphi_) * divB_glm_kappa
          cons(ixO^S, Bphi_cons_) = cons(ixO^S, Bphi_cons_) + qdt * add_source(ixO^S, Bphi_cons_)

          ! add the source terms to Bcons
          do idir = 1, ndim
             call gradient( primCT(ixI^S, Bphi_),ixI^L,ixO^L,idir,grad_phi(ixI^S,idir) )
             add_source(ixO^S, Bcons(idir)) = - alppsi6(ixO^S) * grad_phi(ixO^S, idir) / gamma(ixO^S, idir, idir)
             cons(ixO^S, Bcons(idir)) = cons(ixO^S, Bcons(idir)) + qdt * add_source(ixO^S, Bcons(idir))
          end do
       end if
   
       if ( type_dive == glm ) then
          ! add kappa source terms to Ephi
          add_source(ixO^S, Ephi_cons_) = primCT(ixO^S, alp_) * qrho_e(ixO^S) & 
                   - alppsi6(ixO^S) * primCT(ixO^S, Ephi_) * divE_glm_kappa
          cons(ixO^S, Ephi_cons_) = cons(ixO^S, Ephi_cons_) + qdt * add_source(ixO^S, Ephi_cons_)

          ! add the source terms to Econs
          do idir = 1, ndim
             call gradient( primCT(ixI^S, Ephi_),ixI^L,ixO^L,idir,grad_phi(ixI^S,idir) )
             add_source(ixO^S, Econs(idir)) = - alppsi6(ixO^S) * grad_phi(ixO^S, idir) / gamma(ixO^S, idir, idir)
             cons(ixO^S, Econs(idir)) =  cons(ixO^S, Econs(idir)) + qdt * add_source(ixO^S, Econs(idir))
          end do
       end if

       ! this is the explict part of the current
       do idir = 1, ndir
          cons(ixO^S, Econs(idir)) = cons(ixO^S, Econs(idir)) - qdt * qrho_e(ixO^S) * v_hat(ixO^S, idir)
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
    call get_g_munu(g_munu(ixI^S,0:3,0:3), alp(ixI^S), betai(ixI^S,1:3), gamma(ixI^S,1:3,1:3), ixI^L, ixO^L)

    ! Calculate the 4-velocity of the Eulerian observer
    nmu(ixO^S,0) = 1.0d0 / alp(ixO^S)
    n_mu(ixO^S,0) = - alp(ixO^S)
    do idir = 1, ndir
       nmu(ixO^S, idir) = - betai(ixO^S, idir) / alp(ixO^S)
       n_mu(ixO^S, idir) = 0.0d0
    end do

    ! Calculate the 4-EM fields
    Emu(ixO^S, 0) = 0.0d0
    do idir = 1, ndir
       Emu(ixO^S, idir) = primCT(ixO^S, Evec(idir))
    end do
    B_mu(ixO^S, 0:ndir) = 0.0d0
    do idir = 0, ndir
       ! jdir is not start from 0 since B^0 = 0
       do jdir = 1, ndir
          B_mu(ixO^S, idir) = B_mu(ixO^S, idir) + primCT(ixO^S, Bvec(jdir)) * g_munu(ixO^S, idir, jdir)
       end do
    end do

    ! Maxwell tensor F^{\mu\nu}
    do idir = 0,3
       do jdir = 0,3

          F(ixO^S, idir, jdir) = nmu(ixO^S, idir) * Emu(ixO^S, jdir) - nmu(ixO^S, jdir) * Emu(ixO^S, idir)
          do kdir = 0,3
             do ldir = 0,3
                F(ixO^S, idir, jdir) = F(ixO^S, idir, jdir) &
                    - four_lvc(idir, jdir, kdir, ldir) * B_mu(ixO^S, ldir) * n_mu(ixO^S, kdir) &
                       / ( alp(ixO^S) * sqrt_gamma(ixO^S) ) 
             end do
          end do

       end do
    end do

    ! F^{\mu\nu} * F_{\mu\nu}
    F2(ixO^S) = 2.0d0 * ( B2(ixO^S) - E2(ixO^S) )

    ! Calculate the 4-velocity
    u(ixO^S,0:3) = 0.0d0
    u(ixO^S,0) = lfac(ixO^S) / alp(ixO^S)
    do idir = 1, ndir
       u(ixO^S, idir) = primCT(ixO^S, W_vel(idir)) - lfac(ixO^S) * betai(ixO^S, idir) / alp(ixO^S) 
    end do

    ! total energy-momentum tensorT^{\mu\nu}
    do idir = 0,3
       do jdir = 0,3

          T(ixO^S,idir,jdir) = primCT(ixO^S, rho_) * h(ixO^S) * u(ixO^S, idir) * u(ixO^S, jdir) &
                               + ( primCT(ixO^S, press_) - 0.25d0 * F2(ixO^S) ) * g(ixO^S, idir, jdir)
          do kdir = 0,3
             do ldir = 0,3
                T(ixO^S, idir, jdir) = T(ixO^S, idir, jdir) &
                      + g_munu(ixO^S, kdir, ldir) * F(ixO^S, idir, kdir) * F(ixO^S, jdir, ldir)
             end do
          end do

       enddo
    enddo

    ! Now calculate the source terms for HD (mom, tau) in the compact form
    do idir = 1, ndir
       do jdir = 1, 3 
          do kdir = 1, 3 
             add_source(ixO^S, mom(idir)) = add_source(ixO^S, mom(idir)) &
                                       + 0.5d0 * Dgamma(ixO^S, jdir, kdir, idir) * ( &
                                        T(ixO^S, 0, 0) * betai(ixO^S, kdir) * betai(ixO^S, jdir) &
                                         + 2.0d0 * T(ixO^S, 0, jdir) * betai(ixO^S, kdir) + T(ixO^S, jdir, kdir) )
          end do
       end do
    end do

    do idir = 1, ndir
       add_source(ixO^S, mom(idir)) = add_source(ixO^S, mom(idir)) &
                                 - T(ixO^S, 0, 0) * alp(ixO^S) * dalp(ixO^S, idir) 
       ! T^{0j} * gamma_jk * D_i beta^k
       do kdir = 1, 3 
          do jdir = 1, 3 
             add_source(ixO^S, mom(idir)) = add_source(ixO^S, mom(idir)) &
                                    + T(ixO^S, 0, jdir) * gamma(ixO^S, jdir, kdir) * D_beta(ixO^S, kdir, idir)
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
                  + T(ixO^S,0,0)*( -betai(ixO^S,idir)*dalp(ixO^S,idir) ) &
                  + T(ixO^S,0,idir)*( -dalp(ixO^S,idir) )

       do jdir=1,3
          add_source(ixO^S,tau_) = add_source(ixO^S,tau_) &
                     + T(ixO^S,0,idir)*( 2.0d0*betai(ixO^S,jdir)*K_ij(ixO^S,idir,jdir) )
       enddo

    enddo

    do idir=1,3
       do jdir=1,3
          add_source(ixO^S,tau_) = add_source(ixO^S,tau_) &
                    + T(ixO^S,0,0) * ( K_ij(ixO^S,idir,jdir)*betai(ixO^S,idir)*betai(ixO^S,jdir) ) &
                    + T(ixO^S,idir,jdir) * K_ij(ixO^S,idir,jdir)
       enddo
    enddo

    if ( evolve_hydro ) then
       cons(ixO^S, tau_) = cons(ixO^S, tau_) + qdt * alppsi6(ixO^S) * add_source(ixO^S, tau_)
       do idir =1, ndir
         cons(ixO^S, mom(idir)) = cons(ixO^S, mom(idir)) + qdt * alppsi6(ixO^S) * add_source(ixO^S,mom(idir))
       end do
    end if

    ! Now work out the source term of Bphi if glm is being used
    if ( evolve_EM ) then
       if ( type_divb == glm ) then
          add_source(ixO^S,Bphi_cons_) = 0.0d0 ! reset add source for Bphi
          do idir=1, ndir
             add_source(ixO^S, Bphi_cons_) = add_source(ixO^S, Bphi_cons_) &
                      + K_ij(ixO^S,idir,idir) / gamma(ixO^S,idir,idir) 
          enddo
          add_source(ixO^S, Bphi_cons_) = alp(ixO^S) * primCT(ixO^S, Bphi_) * add_source(ixO^S,Bphi_cons_) 
   
          do idir=1, ndim ! be careful! 
             add_source(ixO^S,Bphi_cons_) = add_source(ixO^S,Bphi_cons_) &
                          + primCT(ixO^S, Bvec(idir)) * dalp(ixO^S, idir)
          enddo
   
          ! add the source terms on cons
          cons(ixO^S, Bphi_cons_) = cons(ixO^S, Bphi_cons_) +  qdt * psi6(ixO^S) * add_source(ixO^S, Bphi_cons_)
       end if

       if ( type_dive == glm ) then
          add_source(ixO^S,Ephi_cons_) = 0.0d0 ! reset add source for Ephi
          do idir=1, ndir
             add_source(ixO^S, Ephi_cons_) = add_source(ixO^S, Ephi_cons_) &
                      + K_ij(ixO^S,idir,idir) / gamma(ixO^S,idir,idir) 
          enddo
          add_source(ixO^S, Ephi_cons_) = alp(ixO^S) * primCT(ixO^S, Ephi_) * add_source(ixO^S,Ephi_cons_) 
   
          do idir=1, ndim ! be careful! 
             add_source(ixO^S,Ephi_cons_) = add_source(ixO^S,Ephi_cons_) &
                          + primCT(ixO^S, Evec(idir)) * dalp(ixO^S, idir)
          enddo
   
          ! add the source terms on cons
          cons(ixO^S, Ephi_cons_) = cons(ixO^S, Ephi_cons_) +  qdt * psi6(ixO^S) * add_source(ixO^S, Ephi_cons_)
       end if
    end if

  end subroutine gremhd_add_source

  !> Add geometrical source terms to w
  subroutine gremhd_add_source_geom(qdt, ixI^L, ixO^L, consCT, primCT, cons, x)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(in)    :: consCT(ixI^S, 1:ncons), primCT(ixI^S, 1:nprim)
    double precision, intent(inout) :: cons(ixI^S, 1:ncons)

    integer                         :: idir,jdir,kdir, h1x^L{^NOONED, h2x^L}
    double precision                :: source(ixI^S)
    double precision                :: fluxCT(ixI^S,1:3,1:3)
    double precision                :: alppsi6(ixI^S)
    double precision                :: h(ixI^S)
    double precision                :: E2(ixI^S), B2(ixI^S) 
    double precision                :: gamma(ixI^S,1:3,1:3) ! metric gamma_ij
    double precision                :: lfac(ixI^S) ! Lorentz factor

    if ( .not. evolve_hydro ) return

    call gremhd_get_intermediate_variables(ixI^L, ixO^L, primCT(ixI^S, 1:nprim), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                h=h(ixI^S), E2=E2(ixI^S), B2=B2(ixI^S) )

    ! alp * psi**6
    alppsi6(ixO^S) = primCT(ixO^S,alp_) * primCT(ixO^S,psi_)**6

    ! Only Momentum flux f^i_j is needed only
    fluxCT = 0.0d0
    do jdir = 1, ndir
       do idir = 1, ndir
         fluxCT(ixO^S, idir, jdir) = - primCT(ixO^S, beta(idir)) * consCT(ixO^S, mom(jdir)) &
                       + alppsi6(ixO^S) * gamma(ixO^S, jdir, jdir) * ( primCT(ixO^S, rho_) * h(ixO^S) * &
                         primCT(ixO^S, W_vel(idir)) * primCT(ixO^S, W_vel(jdir)) &
                       - primCT(ixO^S, Evec(idir)) * primCT(ixO^S, Evec(jdir)) &
                       - primCT(ixO^S, Bvec(idir)) * primCT(ixO^S, Bvec(jdir)) )
       end do
    end do
    ! and the pressure terms
    do jdir = 1, 3
       fluxCT(ixO^S, jdir, jdir) = fluxCT(ixO^S, jdir, jdir) &
               + alppsi6(ixO^S) * ( primCT(ixO^S, press_) + 0.5d0 * ( E2(ixO^S) + B2(ixO^S) )   )
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

       cons(ixO^S, mom(theta_)) = cons(ixO^S, mom(theta_)) + qdt * source(ixO^S)
       ! geo_source[phi] = 0
       }
    end select
  end subroutine gremhd_add_source_geom

end module mod_gremhd_phys_add_source
