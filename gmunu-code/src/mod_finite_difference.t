!> Module with finite difference methods for fluxes
module mod_finite_difference

  implicit none
  private

  public :: fd

contains

  subroutine fd(qdt,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,snew,fC,fE,dx^D,x)
    use mod_physics
    use mod_source, only: addsource2
    use mod_finite_volume, only: reconstruct_LR
    use mod_global_parameters
    use mod_usr_methods

    double precision, intent(in)                                     :: qdt, qtC, qt, dx^D
    integer, intent(in)                                              :: ixI^L, ixO^L, idims^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in)            :: x

    type(state)                                                      :: sCT, snew
    double precision, dimension(ixI^S,1:ncons,1:ndim), intent(out)  :: fC
    double precision, dimension(ixI^S,7-2*ndim:3)                    :: fE

    double precision, dimension(ixI^S,1:ncons)                      :: fCT
    double precision, dimension(ixI^S,1:nprim)                          :: fm, fp, fmR, fpL, wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:ncons) :: consL, consR
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nprim) :: primL, primR
    double precision, dimension(ixI^S)      :: cmaxC
    double precision, dimension(ixI^S)      :: cminC
    double precision, dimension(1:ndim)     :: dxinv, dxdim
    logical                                                          :: transport
    integer                                                          :: idims, iw, ixC^L, ix^L, hxO^L

    type(ct_velocity) :: vcts

    associate(primCT=>sCT%prim, &
              consCT=>sCT%cons, &
              cons_new=>snew%cons)

    fC=0.d0

    ^D&dxinv(^D)=-qdt/dx^D;
    ^D&dxdim(^D)=dx^D;
    do idims= idims^LIM

       ! Get fluxes for the whole grid (mesh+nghostcells)
       {^D& ixmin^D = ixOmin^D - nghostcells * kr(idims,^D)\}
       {^D& ixmax^D = ixOmax^D + nghostcells * kr(idims,^D)\}

       hxO^L=ixO^L-kr(idims,^D);

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D); ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
         ixmax^D=ixmax^D+nghostcells-nghostcells*kr(idims,^D); ixmin^D=ixmin^D-nghostcells+nghostcells*kr(idims,^D);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
       end if

       call phys_get_flux(consCT,primCT,x,ixG^LL,ix^L,idims,fCT)

       do iw=1,ncons
          ! Lax-Friedrich splitting:
          fp(ix^S,iw) = half * (fCT(ix^S,iw) + tvdlfeps * cmax_global * consCT(ix^S,iw))
          fm(ix^S,iw) = half * (fCT(ix^S,iw) - tvdlfeps * cmax_global * consCT(ix^S,iw))
       end do ! iw loop

       ! now do the reconstruction of fp and fm:
       call reconstructL(ixI^L,ixC^L,idims,fp,fpL)
       call reconstructR(ixI^L,ixC^L,idims,fm,fmR)

       fC(ixC^S,1:ncons,idims) = fpL(ixC^S,1:ncons) + fmR(ixC^S,1:ncons)
       !if(associated(usr_set_flux)) call usr_set_flux(ixI^L,ixC^L,qt,wLC,wRC,wLp,wRp,sCT,idims,fC)

       if(stagger_grid) then
         ! apply limited reconstruction for left and right status at cell interfaces
         !call reconstruct_LR(typelimiter, nhydro_lo, nhydro_lo+nreconstruct-1, &
         !      ixI^L,ixC^L,ixC^L,idims,primCT,consL,consR,primL,primR,x,dxdim(idims))
         call phys_get_cbounds(consL,consR,primL,primR,x,ixI^L,ixC^L,idims,cmaxC,cminC)
       end if

    end do !idims loop

    if(stagger_grid) call phys_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,snew,vcts)

    do idims= idims^LIM
       hxO^L=ixO^L-kr(idims,^D);
       do iw = nc_hydro_lo,nc_hydro_hi
          if (slab_uniform) then
             fC(ixI^S,iw,idims) = dxinv(idims) * fC(ixI^S,iw,idims)
             cons_new(ixO^S,iw)=cons_new(ixO^S,iw)+(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
          else
             fC(ixI^S,iw,idims)=-qdt*fC(ixI^S,iw,idims)*block%surfaceC(ixI^S,idims)
             cons_new(ixO^S,iw)=cons_new(ixO^S,iw)+ &
                  (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))/block%dvolume(ixO^S)
          end if
       end do ! iw loop
    end do ! Next idims

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,consCT,primCT,cons_new,x)

    if(stagger_grid) call phys_face_to_center(ixO^L,snew)

    ! check and optionally correct unphysical values
    !if(fix_small_values) then
    !   call phys_handle_small_values(.false.,cons_new,x,ixI^L,ixO^L,'fd')
    !endif

    !call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
    !     ixI^L,ixO^L,1,nprim,qtC,consCT,qt,cons_new,x,.false.)
    end associate

  end subroutine fd

  subroutine reconstructL(ixI^L,iL^L,idims,w,wLC)
    use mod_global_parameters
    !use mod_mp5
    use mod_limiter

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nprim)

    double precision, intent(out)   :: wLC(ixI^S,1:nprim) 

    double precision                :: ldw(ixI^S), dwC(ixI^S)
    integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, iw
    double precision                :: a2max

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterL(1,ncons,ixI^L,iL^L,idims,w,wLC)
!    case (limiter_weno5)
!       call WENO5limiterL(ixI^L,iL^L,idims,w,wLC,1)
!    case (limiter_weno5nm)
!       call WENO5NMlimiterL(ixI^L,iL^L,idims,w,wLC,1)
!    case (limiter_wenoz5)
!       call WENO5limiterL(ixI^L,iL^L,idims,w,wLC,2)
!    case (limiter_wenoz5nm)
!       call WENO5NMlimiterL(ixI^L,iL^L,idims,w,wLC,2)
!    case (limiter_wenozp5)
!       call WENO5limiterL(ixI^L,iL^L,idims,w,wLC,3)
!    case (limiter_wenozp5nm)
!       call WENO5NMlimiterL(ixI^L,iL^L,idims,w,wLC,3)
    case default 

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);

       wLC(kxC^S,1:ncons) = w(kxC^S,1:ncons)

       jxR^L=iL^L+kr(idims,^D);

       ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=1,ncons
          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

          call dwlimiter2(dwC,ixI^L,ixC^L,idims,type_limiter(block%level),ldw)

          wLC(iL^S,iw)=wLC(iL^S,iw)+half*ldw(iL^S)
       end do
    end select

  end subroutine reconstructL

  subroutine reconstructR(ixI^L,iL^L,idims,w,wRC)
    use mod_global_parameters
    !use mod_mp5
    use mod_limiter

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nprim)

    double precision, intent(out)   :: wRC(ixI^S,1:nprim) 

    double precision                :: rdw(ixI^S), dwC(ixI^S)
    integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, kxR^L, iw
    double precision                :: a2max

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterR(1,ncons,ixI^L,iL^L,idims,w,wRC)
!    case (limiter_weno5)
!       call WENO5limiterR(ixI^L,iL^L,idims,w,wRC,1)
!    case (limiter_weno5nm)
!       call WENO5NMlimiterR(ixI^L,iL^L,idims,w,wRC,1)
!    case (limiter_wenoz5)
!       call WENO5limiterR(ixI^L,iL^L,idims,w,wRC,2)
!    case (limiter_wenoz5nm)
!       call WENO5NMlimiterR(ixI^L,iL^L,idims,w,wRC,2)
!    case (limiter_wenozp5)
!       call WENO5limiterR(ixI^L,iL^L,idims,w,wRC,3)
!    case (limiter_wenozp5nm)
!       call WENO5NMlimiterR(ixI^L,iL^L,idims,w,wRC,3)
    case default 

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);

       wRC(kxC^S,1:ncons)=w(kxR^S,1:ncons)

       jxR^L=iL^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=1,ncons
          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

          call dwlimiter2(dwC,ixI^L,ixC^L,idims,type_limiter(block%level),rdw)

          wRC(iL^S,iw)=wRC(iL^S,iw)-half*rdw(jxR^S)
       end do
    end select

  end subroutine reconstructR

end module mod_finite_difference
