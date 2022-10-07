module mod_gremhd_phys_divb_ct
  use mod_physics
  use mod_gremhd_phys_parameters

  implicit none
  private

  ! Public methods
  public :: gremhd_phys_divb_ct_init

contains

  !> Initialize the module
  subroutine gremhd_phys_divb_ct_init()
    use mod_global_parameters
    stagger_grid = .true.
    nconss=ndim
    phys_get_cbounds_and_ct_velocity => gremhd_get_cbounds_and_ct_velocity
    phys_modify_wLR                  => gremhd_modify_wLR
    phys_face_to_center              => gremhd_face_to_center
    phys_update_faces                => gremhd_update_faces
  end subroutine gremhd_phys_divb_ct_init

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  !> Also prepare velocities for ct methods
  subroutine gremhd_get_cbounds_and_ct_velocity(consL, consR, primL, primR, x, ixI^L, ixO^L, idim, vcts, cmax, cmin)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative left and right status
    double precision, intent(in)    :: consL(ixI^S, 1:ncons), consR(ixI^S, 1:ncons)
    ! primitive left and right status
    double precision, intent(in)    :: primL(ixI^S, 1:nprim), primR(ixI^S, 1:nprim)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    type(ct_velocity), intent(inout):: vcts
    double precision, intent(inout) :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision, dimension(ixI^S,1:2) :: lambdaL, lambdaR
    double precision, dimension(ixI^S,1:ndir) :: vhatL, vhatR
    double precision, dimension(ixI^S,1:2) :: tmp_c
    double precision                       :: gamma_hat(ixI^S,1:3,1:3)
    double precision, dimension(ixI^S)     :: sqrt_gamma_hat, inv_sqrt_gamma_hat
    integer                                :: idir,idimE,idimN

      call gremhd_get_lambda(ixI^L, ixO^L, idim, primL(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambdaL(ixI^S,1:2), v_hat=vhatL(ixI^S,1:ndir))
      call gremhd_get_lambda(ixI^L, ixO^L, idim, primR(ixI^S, 1:nprim), x(ixI^S, 1:ndim), lambdaR(ixI^S,1:2), v_hat=vhatR(ixI^S,1:ndir))
  
      tmp_c(ixO^S,1)=max(0.0d0, lambdaL(ixO^S,1), lambdaR(ixO^S,1) )
      tmp_c(ixO^S,2)=min(0.0d0, lambdaL(ixO^S,2), lambdaR(ixO^S,2) ) 
  
      if(present(cmin)) then
        cmax(ixO^S) = tmp_c(ixO^S,1)
        cmin(ixO^S) = tmp_c(ixO^S,2)
      else
        cmax(ixO^S) = max(abs(tmp_c(ixO^S,1)), abs(tmp_c(ixO^S,2)))
      end if

      ! calculate velocities related to different CT schemes
      select case(type_ct)
      case('average')
      case('uct_hll')
        if(.not.allocated(vcts%alpELC)) then
          allocate(vcts%cbarmin(ixI^S,1:ndim),vcts%cbarmax(ixI^S,1:ndim)) 
          allocate(vcts%alpELC(ixI^S,1:ndir),vcts%alpERC(ixI^S,1:ndir))
          allocate(vcts%betaC(ixI^S,1:ndir,2))
        end if
        ! Store magnitude of characteristics
        if(present(cmin)) then
          vcts%cbarmin(ixO^S,idim)=max(-cmin(ixO^S),zero)
          vcts%cbarmax(ixO^S,idim)=max( cmax(ixO^S),zero)
        else
          vcts%cbarmax(ixO^S,idim)=max( cmax(ixO^S),zero)
          vcts%cbarmin(ixO^S,idim)=vcts%cbarmax(ixO^S,idim)
        end if

        idimN=mod(idim,ndir)+1 ! 'Next' direction
        idimE=mod(idim+1,ndir)+1 ! Electric field direction

        ! Store betas
        ! note that in general betaL are the same as betaR
        vcts%betaC(ixO^S,idim,1) = 0.5d0 * (primL(ixO^S,beta(idimN)) + primR(ixO^S,beta(idimN)))
        vcts%betaC(ixO^S,idim,2) = 0.5d0 * (primL(ixO^S,beta(idimE)) + primR(ixO^S,beta(idimE)))

        call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_hat(ixI^S,1:3,1:3))
        call get_sqrt_gamma_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, sqrt_gamma_hat(ixI^S))
        where ( sqrt_gamma_hat(ixO^S) > smalldouble )
           inv_sqrt_gamma_hat(ixO^S) = 1.0d0 / sqrt_gamma_hat(ixO^S)
        else where
           inv_sqrt_gamma_hat(ixO^S) = 0.0d0
        end where

        ! Store alpE over sqrt_gamma_hat
        vcts%alpELC(ixO^S,idim)=primL(ixO^S, alp_) * primL(ixO^S,Evec(idimE)) &
                * primL(ixO^S, psi_)**4 * gamma_hat(ixO^S, idimE, idimE) * inv_sqrt_gamma_hat(ixO^S) 
        vcts%alpERC(ixO^S,idim)=primR(ixO^S, alp_) * primR(ixO^S,Evec(idimE)) &
                * primR(ixO^S, psi_)**4 * gamma_hat(ixO^S, idimE, idimE) * inv_sqrt_gamma_hat(ixO^S)

      case default
        call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
      end select

  end subroutine gremhd_get_cbounds_and_ct_velocity

  subroutine gremhd_modify_wLR(ixI^L,ixO^L,qt,consL,consR,primL,primR,s,idir)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: consL(ixI^S,1:ncons), consR(ixI^S,1:ncons)
    double precision, intent(inout) :: primL(ixI^S,1:nprim), primR(ixI^S,1:nprim)
    type(state)                     :: s
    if(stagger_grid) then
      consL(ixO^S,Bcons(idir)) = s%conss(ixO^S,idir)
      consR(ixO^S,Bcons(idir)) = s%conss(ixO^S,idir)
      primL(ixO^S,Bvec(idir))  = s%conss(ixO^S,idir) / primL(ixO^S, psi_)**6
      primR(ixO^S,Bvec(idir))  = s%conss(ixO^S,idir) / primR(ixO^S, psi_)**6
    end if
  end subroutine gremhd_modify_wLR

  !> calculate cell-center values from face-center values
  subroutine gremhd_face_to_center(ixO^L,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixO^L
    type(state)                        :: s

    integer                            :: fxO^L, gxO^L, hxO^L, jxO^L, kxO^L, idim

    associate(prim=>s%prim, cons=>s%cons, conss=>s%conss)

    select case( face_to_center_order )
    case (2)
       ! calculate cell-center values from face-center values in 2nd order
       do idim=1,ndim
         ! Displace index to the left
         ! Even if ixI^L is the full size of the prim arrays, this is ok
         ! because the staggered arrays have an additional place to the left.
         hxO^L=ixO^L-kr(idim,^D);
         ! Interpolate to cell barycentre using arithmetic average
         ! This might be done better later, to make the method less diffusive.
         s%cons(ixO^S,Bcons(idim))=0.5d0/s%surface(ixO^S,idim)*&
           (conss(ixO^S,idim)*s%surfaceC(ixO^S,idim)&
           +conss(hxO^S,idim)*s%surfaceC(hxO^S,idim))
       end do
   
    case (4)
       ! calculate cell-center values from face-center values in 4th order
       do idim=1,ndim
         gxO^L=ixO^L-2*kr(idim,^D);
         hxO^L=ixO^L-kr(idim,^D);
         jxO^L=ixO^L+kr(idim,^D);
   
         ! Interpolate to cell barycentre using fourth order central formula
         s%cons(ixO^S,Bcons(idim))=(0.0625d0/s%surface(ixO^S,idim))*&
                ( -conss(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
            +9.0d0*conss(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
            +9.0d0*conss(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
                  -conss(jxO^S,idim)*s%surfaceC(jxO^S,idim) )
       end do
   
    case (6)
       ! calculate cell-center values from face-center values in 6th order
       do idim=1,ndim
         fxO^L=ixO^L-3*kr(idim,^D);
         gxO^L=ixO^L-2*kr(idim,^D);
         hxO^L=ixO^L-kr(idim,^D);
         jxO^L=ixO^L+kr(idim,^D);
         kxO^L=ixO^L+2*kr(idim,^D);
   
         ! Interpolate to cell barycentre using sixth order central formula
         s%cons(ixO^S,Bcons(idim))=(0.00390625d0/s%surface(ixO^S,idim))* &
            (  +3.0d0*conss(fxO^S,idim)*s%surfaceC(fxO^S,idim) &
              -25.0d0*conss(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
             +150.0d0*conss(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
             +150.0d0*conss(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
              -25.0d0*conss(jxO^S,idim)*s%surfaceC(jxO^S,idim) &
               +3.0d0*conss(kxO^S,idim)*s%surfaceC(kxO^S,idim) )
       end do
    end select

    ! also update cell-center Bvec
    do idim = 1, ndim
       prim(ixO^S, Bvec(idim)) = cons(ixO^S, Bcons(idim)) / prim(ixO^S, psi_)**6 
    end do

    end associate

  end subroutine gremhd_face_to_center

  subroutine gremhd_update_faces(ixI^L,ixO^L,qt,qdt,prim,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: prim(ixI^S,1:nprim)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixI^S,1:ncons,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    integer                            :: hxC^L,ixC^L
    integer                            :: idim1,idim2,idir
    double precision                   :: circ(ixI^S,1:ndim)

    ! note that fC here is actually fdA
    select case(type_ct)
    case('average')
      call update_faces_average(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s)
    !case('uct_contact')
    !  call update_faces_contact(ixI^L,ixO^L,qt,qdt,prim,fC,fE,sCT,s,vcts)
    case('uct_hll')
      call update_faces_hll(ixI^L,ixO^L,qt,qdt,fE,sCT,s,vcts)
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

    ! Calculate circulation on each face: interal(fE dot dl)
    circ(ixI^S,1:ndim)=0.0d0
    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1) &
                  + lvc(idim1,idim2,idir) * ( fE(ixC^S,idir) - fE(hxC^S,idir) )
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      !where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
      where(s%surfaceC(ixC^S,idim1) > tiny(0.0d0))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=0.0d0
      end where
      ! Time update cell-face magnetic field component
      s%conss(ixC^S,idim1)=s%conss(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

  end subroutine gremhd_update_faces

  !> get electric field though averaging neighors to update faces in CT
  subroutine update_faces_average(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixI^S,1:ncons,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    integer                            :: hxC^L,ixC^L,jxC^L !,ixCm^L
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(x=>s%x)
    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D-1;

    fE=0.0d0

    do idim1=1,ndim 
      iwdim1 = Bcons(idim1)
      do idim2=1,ndim
        iwdim2 = Bcons(idim2)
        do idir=7-2*ndim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ! Assemble indices
            jxC^L=ixC^L+kr(idim1,^D);
            hxC^L=ixC^L+kr(idim2,^D);

            ! Interpolate to edges
            fE(ixC^S,idir) = 0.25d0 * ( &
                  (fC(ixC^S,iwdim1,idim2)+fC(jxC^S,iwdim1,idim2)) &
                 -(fC(ixC^S,iwdim2,idim1)+fC(hxC^S,iwdim2,idim1)) )
                 ! (fC(ixC^S,iwdim1,idim2)+fC(jxC^S,iwdim1,idim2))/dxlevel(idim1) &
                 !-(fC(ixC^S,iwdim2,idim1)+fC(hxC^S,iwdim2,idim1))/dxlevel(idim2)  )

            fE(ixC^S,idir)=qdt * s%dsC(ixC^S,idir) * fE(ixC^S,idir)
            !fE(ixC^S,idir)=qdt * fE(ixC^S,idir)
          end if

        end do
      end do
    end do
    end associate
  end subroutine update_faces_average

  !> update faces using UCT contact mode by Del Zanna, L., Zanotti, O., Bucciantini, N., & Londrillo, P. in Astronomy & Astrophysics , 473, 11 (2007)
  subroutine update_faces_hll(ixI^L,ixO^L,qt,qdt,fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_geometry

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts

    double precision                   :: betaL(ixI^S,2), betaR(ixI^S,2)
    double precision                   :: alpELL(ixI^S), alpELR(ixI^S)
    double precision                   :: alpERL(ixI^S), alpERR(ixI^S)
    double precision                   :: EhatLL(ixI^S), EhatLR(ixI^S)
    double precision                   :: EhatRL(ixI^S), EhatRR(ixI^S)
    double precision                   :: btilL(s%ixGs^S,ndim)
    double precision                   :: btilR(s%ixGs^S,ndim)
    double precision                   :: cp(ixI^S,2)
    double precision                   :: cm(ixI^S,2)
    double precision                   :: sqrt_gamma_hat(s%ixGs^S)
    double precision                   :: xC(s%ixGs^S,1:ndim),xCC(s%ixGs^S,1:ndim)
    integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,ixCm^L
    integer                            :: idim1,idim2,idir

    associate(bfacesCT=>sCT%conss,x=>s%x,&
              alpELC=>vcts%alpELC,alpERC=>vcts%alpERC,&
              betaC=>vcts%betaC,&
              cbarmin=>vcts%cbarmin,cbarmax=>vcts%cbarmax  )

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! Loop over components of electric field

    ! idir: electric field component we need to calculate
    ! idim1: directions in which we already performed the reconstruction
    ! idim2: directions in which we perform the reconstruction

    fE=zero

    ! extend one layer of cell center locations in xCC
    xCC=0.d0
    xCC(ixI^S,1:ndim)=x(ixI^S,1:ndim)
    {
    xCC(s%ixGsmin^D^D%ixI^S,1:ndim) = x(ixImin^D^D%ixI^S,1:ndim)
    xCC(s%ixGsmin^D^D%s%ixGs^S,^D) = x({ixImin^DD,},^D) - block%dx({ixImin^DD,},^D)
    \}
    {^IFTHREED
    xCC(ixImin1:ixImax1,s%ixGsmin2,s%ixGsmin3,1)=x(ixImin1:ixImax1,ixImin2,ixImin3,1)
    xCC(s%ixGsmin1,ixImin2:ixImax2,s%ixGsmin3,2)=x(ixImin1,ixImin2:ixImax2,ixImin3,2)
    xCC(s%ixGsmin1,s%ixGsmin2,ixImin3:ixImax3,3)=x(ixImin1,ixImin2,ixImin3:ixImax3,3)
    }

    do idir=7-2*ndim,3
      ! Indices
      ! idir: electric field component
      ! idim1: one surface
      ! idim2: the other surface
      ! cyclic permutation: idim1,idim2,idir=1,2,3
      ! Velocity components on the surface
      ! follow cyclic premutations:
      ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-1+kr(idir,^D);

      ! Get edge coordinates and sqrt_gamma_hat
      do idim1=1,ndim
        if (idim1/=idir) then
          xC(ixC^S,idim1)=xCC(ixC^S,idim1) + 0.5d0 * s%dx(ixC^S,idim1)
        else
          xC(ixC^S,idim1)=xCC(ixC^S,idim1)
        end if
      end do
      call get_sqrt_gamma_hat(xC(s%ixGs^S, 1:ndim), s%ixGs^L, ixC^L, sqrt_gamma_hat(s%ixGs^S))

      ! Set indices and directions
      idim1=mod(idir,ndir)+1
      idim2=mod(idir+1,ndir)+1

      jxC^L=ixC^L+kr(idim1,^D);
      ixCp^L=ixC^L+kr(idim2,^D);

      ! Reconstruct electric fields
      call reconstruct(ixI^L,ixC^L,idim2, alpELC(ixI^S,idim1),&
               alpELL(ixI^S),alpELR(ixI^S))
      call reconstruct(ixI^L,ixC^L,idim2, alpERC(ixI^S,idim1),&
               alpERL(ixI^S),alpERR(ixI^S))

      ! Reconstruct betas
      call reconstruct(ixI^L,ixC^L,idim2,betaC(ixI^S,idim1,1),&
               betaL(ixI^S,1),betaR(ixI^S,1)) 
      call reconstruct(ixI^L,ixC^L,idim2,betaC(ixI^S,idim1,2),&
               betaL(ixI^S,2),betaR(ixI^S,2)) 

      ! Reconstruct magnetic fields
      ! Eventhough the arrays are larger, reconstruct works with
      ! the limits ixG.
      call reconstruct(ixI^L,ixC^L,idim2,bfacesCT(ixI^S,idim1),&
               btilL(ixI^S,idim1),btilR(ixI^S,idim1))
      call reconstruct(ixI^L,ixC^L,idim1,bfacesCT(ixI^S,idim2),&
               btilL(ixI^S,idim2),btilR(ixI^S,idim2))

      ! compute electric field Ehat
      EhatLL(ixC^S) = alpELL(ixC^S) + sqrt_gamma_hat(ixC^S) * ( &
                       + betaL(ixC^S,1)*btilL(ixC^S,idim2) &
                       - betaL(ixC^S,2)*btilL(ixC^S,idim1) )
      EhatLR(ixC^S) = alpELR(ixC^S) + sqrt_gamma_hat(ixC^S) * ( &
                       + betaL(ixC^S,1)*btilL(ixC^S,idim2) &
                       - betaR(ixC^S,2)*btilR(ixC^S,idim1) )
      EhatRL(ixC^S) = alpERL(ixC^S) + sqrt_gamma_hat(ixC^S) * ( &
                       + betaR(ixC^S,1)*btilR(ixC^S,idim2) &
                       - betaL(ixC^S,2)*btilL(ixC^S,idim1) )
      EhatRR(ixC^S) = alpERR(ixC^S) + sqrt_gamma_hat(ixC^S) * ( &
                       + betaR(ixC^S,1)*btilR(ixC^S,idim2) &
                       - betaR(ixC^S,2)*btilR(ixC^S,idim1) )

      ! Take the maximum characteristic
      cm(ixC^S,1)=max(cbarmin(ixCp^S,idim1),cbarmin(ixC^S,idim1))
      cp(ixC^S,1)=max(cbarmax(ixCp^S,idim1),cbarmax(ixC^S,idim1))

      cm(ixC^S,2)=max(cbarmin(jxC^S,idim2),cbarmin(ixC^S,idim2))
      cp(ixC^S,2)=max(cbarmax(jxC^S,idim2),cbarmax(ixC^S,idim2))
     
      ! Calculate eletric field over sqrt_gamma_hat
      ! fixme: possible 1/0
      fE(ixC^S,idir)=( cp(ixC^S,1)*cp(ixC^S,2)*EhatLL(ixC^S) &
                     + cp(ixC^S,1)*cm(ixC^S,2)*EhatLR(ixC^S) &
                     + cm(ixC^S,1)*cp(ixC^S,2)*EhatRL(ixC^S) &
                     + cm(ixC^S,1)*cm(ixC^S,2)*EhatRR(ixC^S) ) &
                       / (cp(ixC^S,1)+cm(ixC^S,1)) / (cp(ixC^S,2)+cm(ixC^S,2)) &
                     + ( cp(ixC^S,1)*cm(ixC^S,1) * (btilR(ixC^S,idim2) - btilL(ixC^S,idim2)) ) &
                     / (cp(ixC^S,1)+cm(ixC^S,1)) &
                     - ( cp(ixC^S,2)*cm(ixC^S,2) * (btilR(ixC^S,idim1) - btilL(ixC^S,idim1)) ) &
                     / (cp(ixC^S,2)+cm(ixC^S,2)) 

      fE(ixC^S,idir) = qdt * s%dsC(ixC^S,idir) * fE(ixC^S,idir)

    end do
    end associate
  end subroutine update_faces_hll

end module mod_gremhd_phys_divb_ct
