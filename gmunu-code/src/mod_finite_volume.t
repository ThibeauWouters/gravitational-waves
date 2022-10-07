!> Module with finite volume methods for fluxes
module mod_finite_volume
  implicit none
  private

  public :: finite_volume
  public :: reconstruct_LR

contains

  !> finite volume method
  subroutine finite_volume(method,qdt,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,snew,sold,fC,fE,dx^D,xbar)
    use mod_physics
    use mod_global_parameters
    use mod_limiter
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
    use mod_usr_methods
    use mod_cfc

    integer, intent(in)                                   :: method
    double precision, intent(in)                          :: qdt, qtC, qt, dx^D
    integer, intent(in)                                   :: ixI^L, ixO^L, idims^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in) :: xbar
    type(state)                                           :: sCT, snew, sold
    double precision, dimension(ixI^S,1:ncons,1:ndim)     :: fC
    double precision, dimension(ixI^S,7-2*ndim:3)         :: fE

    ! cell-face location coordinates
    double precision, dimension(ixI^S,1:ndim)             :: xi
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:ncons)            :: consL, consR
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nprim)            :: primL, primR
    double precision, dimension(ixI^S,1:ncons)            :: fLC, fRC
    double precision, dimension(ixI^S)                    :: cmaxC, cminC
    double precision, dimension(ixI^S)                    :: inv_volume
    double precision, dimension(1:ndim)                   :: dxinv, dxdim
    integer, dimension(ixI^S)                             :: patchf
    integer                                               :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L

    ! when CT is being used
    type(ct_velocity)                                     :: vcts

    ! variables that related to positivity preserving limiter
    ! low order flux
    double precision, dimension(:^D&,:,:), allocatable    :: fC_low
    !double precision, parameter                           :: epsD = smalldouble, epstau = smalldouble
    double precision, parameter                           :: epsD = 0.0d0, epstau = 0.0d0 

    associate(x=>sCT%x, primCT=>sCT%prim, consCT=>sCT%cons, &
              cons_new=>snew%cons, cons_old=>sold%cons)

    if (positivity_preserving) then
       allocate(fC_low(ixI^S,1:ncons,1:ndim))
       fC_low=0.0d0
    end if

    fC=0.0d0

    ! The flux calculation contracts by one in the idims direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ix^L=ixO^L;
    do idims= idims^LIM
       ix^L=ix^L^LADD2*kr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in fv : Nonconforming input limits")

    ^D&dxinv(^D)=-qdt/dx^D;
    ^D&dxdim(^D)=dx^D;

    do idims= idims^LIM

       ! Assemble indices
       hxO^L=ixO^L-kr(idims,^D);
       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);

       if(stagger_grid) then
          ! ct needs all transverse cells
          ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D); ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
       end if

       ! Determine stencil size, indices to reconstruct to
       {ixCRmin^D = max(ixCmin^D - phys_wider_stencil,ixGlo^D)\}
       {ixCRmax^D = min(ixCmax^D + phys_wider_stencil,ixGhi^D)\}

       ! get cell-face coordinates from cell-center
       xi(ixI^S,1:ndim) = x(ixI^S,1:ndim)
       xi(ixI^S,idims)  = xi(ixI^S,idims)+0.5d0*sCT%dx(ixI^S,idims)

       ! primR and primL are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idims they are
       ! shifted half a cell towards the 'higher' direction.
       primR(kxC^S,1:nprim)=primCT(kxR^S,1:nprim)
       primL(kxC^S,1:nprim)=primCT(kxC^S,1:nprim)

       ! low order TVD-Lax-Friedrich flux first.
       if (positivity_preserving) then
         ! fixme: maybe this part can be faster
         !primR(kxC^S,1:nprim)=primCT(kxR^S,1:nprim)
         !primL(kxC^S,1:nprim)=primCT(kxC^S,1:nprim)
         consR(kxC^S,1:ncons)=consCT(kxR^S,1:ncons)
         consL(kxC^S,1:ncons)=consCT(kxC^S,1:ncons)
         ! special modification of left and right status before flux evaluation
         call phys_modify_wLR(ixI^L,ixCR^L,qt,consL,consR,primL,primR,sCT,idims)
         call phys_get_flux(consL,primL,xi,ixI^L,ixC^L,idims,fLC)
         call phys_get_flux(consR,primR,xi,ixI^L,ixC^L,idims,fRC)
         call phys_get_cbounds(consL,consR,primL,primR,xi,ixI^L,ixC^L,idims,cmaxC)
         call get_Riemann_flux_tvdlf()
         if ( slab_uniform ) then
           do iw=nc_hydro_lo,nc_hydro_hi
             fC_low(ixC^S,iw,idims)=fC(ixC^S,iw,idims)
           end do
         else 
           do iw=nc_hydro_lo,nc_hydro_hi
             fC_low(ixC^S,iw,idims)=fC(ixC^S,iw,idims)*block%surfaceC(ixC^S,idims)
           end do
         end if
       end if ! end pp limiter

       !primR(kxC^S,nmetric_lo:nmetric_hi)=primCT(kxR^S,nmetric_lo:nmetric_hi)
       !primL(kxC^S,nmetric_lo:nmetric_hi)=primCT(kxC^S,nmetric_lo:nmetric_hi)
       ! get cell-face metric
       ! fixme: maybe this part can be faster
       if ( .not. reconstruct_cfc ) then
          call cfc_metric_interpolation(ixI^L,ixCR^L,idims,primCT,xbar,primL,xi)
          call cfc_metric_interpolation(ixI^L,ixCR^L,idims,primCT,xbar,primR,xi)
          !primR(kxC^S,nmetric_lo:nmetric_hi)=primL(kxR^S,nmetric_lo:nmetric_hi)
       end if

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(type_limiter(block%level), &
         ixI^L,ixCR^L,ixCR^L,idims,primCT,consL,consR,primL,primR,xbar,xi,dxdim(idims))

       ! special modification of left and right status before flux evaluation
       call phys_modify_wLR(ixI^L,ixCR^L,qt,consL,consR,primL,primR,sCT,idims)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(consL,primL,xi,ixI^L,ixC^L,idims,fLC)
       call phys_get_flux(consR,primR,xi,ixI^L,ixC^L,idims,fRC)

       ! estimating bounds for the minimum and maximum signal velocities
       if(method==fs_tvdlf.or.method==fs_tvdmu) then
         if (.not.stagger_grid) then
           call phys_get_cbounds(consL,consR,primL,primR,xi,ixI^L,ixC^L,idims,cmaxC)
         else
           call phys_get_cbounds_and_ct_velocity(consL,consR,primL,primR,xi,ixI^L,ixC^L,idims,vcts,cmaxC)
         end if
       else
         if (.not.stagger_grid) then
           call phys_get_cbounds(consL,consR,primL,primR,xi,ixI^L,ixC^L,idims,cmaxC,cminC)
         else
           call phys_get_cbounds_and_ct_velocity(consL,consR,primL,primR,xi,ixI^L,ixC^L,idims,vcts,cmaxC,cminC)
         end if
       end if

       ! use approximate Riemann solver to get flux at interfaces
       select case(method)
       case(fs_hll)
         call get_Riemann_flux_hll()
       case(fs_tvdlf)
         call get_Riemann_flux_tvdlf()
       case(fs_tvdmu)
         call get_Riemann_flux_tvdmu()
       case default
         call mpistop('unkown Riemann flux in finite volume')
       end select

       ! special modification of flux
       ! fixme: for M1, dont use it at the moment
       !call phys_modify_flux(ixI^L,ixCR^L,consL,consR,primL,primR,xi,sCT,idims,fC)

       if (.not.slab_uniform) then
          do iw=nc_hydro_lo,nc_hydro_hi
            fC(ixC^S,iw,idims)=fC(ixC^S,iw,idims)*block%surfaceC(ixC^S,idims)
          end do
       end if

       if (.not.slab.and.idimsmin==1) &
            call phys_add_source_geom_fv(qdt,ixI^L,ixO^L,consL,consR,primL,primR,consCT,primCT,cons_new,xbar,idims)

       ! If use positivity preserving limiter, with fC and fC_low, work out the modify the flux
       if (positivity_preserving) call positivity_preserving_limiter()

    end do ! Next idims

    if(associated(usr_set_flux)) then
       do idims= idims^LIM
          call usr_set_flux(ixI^L,ixC^L,idims,fC,xi)
       end do ! Next idims
    end if

    if (stagger_grid) call phys_update_faces(ixI^L,ixO^L,qt,qdt,primCT,fC,fE,sCT,snew,vcts)

    do idims= idims^LIM
       hxO^L=ixO^L-kr(idims,^D);

       do iw=nc_hydro_lo,nc_hydro_hi
         if ( flux_type(idims, iw) == flux_nul ) then
            fC(ixI^S,iw,idims)=0.0d0
         end if
       end do

       ! Multiply the fluxes by -dt/dx since Flux fixing expects this
       if (slab_uniform) then
          do iw = nc_hydro_lo,nc_hydro_hi
            fC(ixI^S,iw,idims)=dxinv(idims)*fC(ixI^S,iw,idims)

            if (associated(phys_iw_methods(iw)%inv_capacity)) then
              call phys_iw_methods(iw)%inv_capacity(ixI^L, ixO^L, cons_new, inv_volume)
              cons_new(ixO^S,iw)=cons_new(ixO^S,iw) + inv_volume * &
                   (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
            else
              cons_new(ixO^S,iw)=cons_new(ixO^S,iw) &
                   + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
            end if

          end do
        else
          if (associated(phys_iw_methods(iw)%inv_capacity)) then
            call phys_iw_methods(iw)%inv_capacity(ixI^L, ixO^L, cons_new, inv_volume)
          else
            inv_volume(ixO^S) = 1.0d0
          end if
          inv_volume(ixO^S) = inv_volume(ixO^S)/block%dvolume(ixO^S)

          do iw=nc_hydro_lo,nc_hydro_hi
            fC(ixI^S,iw,idims)=-qdt*fC(ixI^S,iw,idims)
            cons_new(ixO^S,iw)=cons_new(ixO^S,iw) &
                  + inv_volume(ixO^S) * (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims)) 
          enddo
       end if

       ! For the MUSCL scheme apply the characteristic based limiter
       if (method==fs_tvdmu) &
            call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,consL,consR,cons_new,xbar,fC,dx^D)

    end do ! Next idims

    if (.not.slab.and.idimsmin==1) &
         call phys_add_source_geom(qdt,ixI^L,ixO^L,consCT,primCT,cons_new,x)

    if (stagger_grid) call phys_face_to_center(ixO^L,snew)

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nprim,qtC,primCT,qt,cons_new,xbar,.false.)

    if (positivity_preserving) then
       deallocate(fC_low)
    end if

  end associate

  contains

    subroutine get_Riemann_flux_tvdmu()
      do iw=nc_hydro_lo,nc_hydro_hi
         ! To save memory we use fLC to store (F_L+F_R)/2=0.5d0*(fLC+fRC)
         fLC(ixC^S, iw)=0.5d0*(fLC(ixC^S, iw)+fRC(ixC^S, iw))
         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf()
      double precision :: fac(ixC^S)

      fac = -0.5d0*tvdlfeps*cmaxC(ixC^S)

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=nc_hydro_lo,nc_hydro_hi

         ! To save memory we use fLC to store (F_L+F_R)/2=0.5d0*(fLC+fRC)
         fLC(ixC^S, iw)=0.5d0*(fLC(ixC^S, iw)+fRC(ixC^S, iw))

         ! Add TVDLF dissipation to the flux
         if (flux_type(idims, iw) /= flux_no_dissipation) then
            fLC(ixC^S, iw)=fLC(ixC^S, iw) + fac*(consR(ixC^S,iw)-consL(ixC^S,iw))
         end if

         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
      end do ! Next iw
    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll()

      double precision :: fac(ixC^S), div(ixC^S)

      where(cminC(ixC^S) >= zero)
        patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S) <= zero)
        patchf(ixC^S) =  2
      elsewhere
        patchf(ixC^S) =  1
        fac(ixC^S) = tvdlfeps*cminC(ixC^S)*cmaxC(ixC^S)
        div(ixC^S) = 1/(cmaxC(ixC^S)-cminC(ixC^S))
      endwhere

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=nc_hydro_lo,nc_hydro_hi
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixC^S, iw) = 0.5d0*(fLC(ixC^S, iw) + fRC(ixC^S, iw) &
                 -tvdlfeps*max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                 (consR(ixC^S,iw)-consL(ixC^S,iw)))
         else
            where(patchf(ixC^S)==1)
               ! Add hll dissipation to the flux
               fLC(ixC^S, iw) = (cmaxC(ixC^S)*fLC(ixC^S, iw)-cminC(ixC^S) * fRC(ixC^S, iw) &
                    +fac(ixC^S)*(consR(ixC^S,iw)-consL(ixC^S,iw))) * div(ixC^S)
            elsewhere(patchf(ixC^S)== 2)
               fLC(ixC^S, iw)=fRC(ixC^S, iw)
            elsewhere(patchf(ixC^S)==-2)
               fLC(ixC^S, iw)=fLC(ixC^S, iw)
            endwhere
         endif

         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
      end do ! Next iw
    end subroutine get_Riemann_flux_hll

    subroutine positivity_preserving_limiter()
      implicit none
      double precision, dimension(ixI^S)   :: inv_lambda
      double precision, dimension(ixI^S)   :: theta, theta_tau

      if (slab_uniform) then
         inv_lambda(ixI^S) = dxdim(idims)/(qdt)
      else
         inv_lambda(ixI^S) = block%dvolume(ixI^S)/(qdt)
      end if

      call get_theta(ixI^L,ixC^L,idims,epsD,inv_lambda(ixI^S),sCT%cons(ixI^S,D_),fC_low(ixI^S,D_,idims),fC(ixI^S,D_,idims),theta)
      call get_theta(ixI^L,ixC^L,idims,epstau,inv_lambda(ixI^S),sCT%cons(ixI^S,tau_),fC_low(ixI^S,tau_,idims),fC(ixI^S,tau_,idims),theta_tau)

      theta(ixC^S) = min(theta(ixC^S),theta_tau(ixC^S))

      do iw = nc_hydro_lo,nc_hydro_hi
         ! note that pp limiter cannot act on stagger grid
         if ( stagger_grid ) then
            if ( (iw>=Bcons(1)) .and. (iw<=Bcons(ndim)) ) cycle
         end if
         fC(ixC^S,iw,idims) = theta(ixC^S)*fC(ixC^S,iw,idims) &
              + (1.0d0-theta(ixC^S))*fC_low(ixC^S,iw,idims)
      end do ! Next iw
    end subroutine positivity_preserving_limiter

  end subroutine finite_volume

  !> Determine the upwinded consL(ixL) and consR(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(typelimiter_in,ixI^L,ixL^L,ixR^L,idims, &
         prim,consL,consR,primL,primR,xbar,xi,dxdim)
    use mod_physics
    use mod_eos
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: typelimiter_in
    integer, intent(in) :: ixI^L, ixL^L, ixR^L, idims
    double precision, intent(in) :: dxdim
    double precision, dimension(ixI^S,1:nprim) :: prim
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:ncons) :: consL, consR
    ! left and right constructed status in primitive form
    double precision, dimension(ixI^S,1:nprim) :: primL, primR
    double precision, dimension(ixI^S,1:ndim) :: xbar,xi

    integer            :: jxR^L, ixC^L, jxC^L, iw
    double precision   :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)
    logical            :: fattening = .False.

    do iw=1, nprim
       if (.not. var_reconstruct(iw) ) cycle

       select case (typelimiter_in)
       case (limiter_venk)
          call venklimiter(ixI^L,ixL^L,idims,dxdim,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw)) 
       case (limiter_mp5)
          call MP5limiter(ixI^L,ixL^L,idims,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw))
       case (limiter_weno3)
          call WENO3limiter(ixI^L,ixL^L,idims,dxdim,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),1)
       case (limiter_wenoyc3)
          call WENO3limiter(ixI^L,ixL^L,idims,dxdim,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),2)
       case (limiter_weno5)
          call WENO5limiter(ixI^L,ixL^L,idims,dxdim,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),1)
       case (limiter_weno5nm)
          call WENO5NMlimiter(ixI^L,ixL^L,idims,dxdim,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),1)
       case (limiter_wenoz5)
          call WENO5limiter(ixI^L,ixL^L,idims,dxdim,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),2)
       case (limiter_wenoz5nm)
          call WENO5NMlimiter(ixI^L,ixL^L,idims,dxdim,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),2)
       case (limiter_wenozp5)
          call WENO5limiter(ixI^L,ixL^L,idims,dxdim,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),3)
       case (limiter_wenozp5nm)
          call WENO5NMlimiter(ixI^L,ixL^L,idims,dxdim,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),3)
       case (limiter_weno5cu6)
          call WENO5CU6limiter(ixI^L,ixL^L,idims,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw))
       case (limiter_teno5ad)
          call TENO5ADlimiter(ixI^L,ixL^L,idims,dxdim,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw))
       case (limiter_weno7)
          call WENO7limiter(ixI^L,ixL^L,idims,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),1)
       case (limiter_mpweno7)
          call WENO7limiter(ixI^L,ixL^L,idims,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),2)
       case (limiter_exeno7)
          call exENO7limiter(ixI^L,ixL^L,idims,prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw))
       case (limiter_ppm)
          fattening = .True. ! fixme: remove this
          ! Since PPM uses the ordinary grid-index:
          ixCmin^D=ixLmin^D+kr(^D,idims);
          ixCmax^D=ixLmax^D;
          ! our fattening is only available for hydro
          call PPMlimiter(ixI^L,ixC^L,idims,prim(ixI^S,iw),prim(ixI^S,iw),primL(ixI^S,iw),primR(ixI^S,iw),fattening)
       case default
          jxR^L=ixR^L+kr(idims,^D);
          ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idims,^D);
          jxC^L=ixC^L+kr(idims,^D);
   
          if (loglimit(iw)) then
             prim(ixCmin^D:jxCmax^D,iw)=dlog10(prim(ixCmin^D:jxCmax^D,iw))
             primL(ixL^S,iw)=dlog10(primL(ixL^S,iw))
             primR(ixR^S,iw)=dlog10(primR(ixR^S,iw))
          end if

          ! original version
          !dwC(ixC^S)=( prim(jxC^S,iw)-prim(ixC^S,iw) )
          ! limit flux from left and/or right
          !call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter_in,ldw,rdw)
          !primL(ixL^S,iw)=primL(ixL^S,iw)+0.5d0*ldw(ixL^S)
          !primR(ixR^S,iw)=primR(ixR^S,iw)-0.5d0*rdw(jxR^S)

          ! new version
          !dxi(ixC^S) = dxdim!xi(jxC^S,idims) - xi(ixC^S,idims)
          dwC(ixC^S)=( prim(jxC^S,iw)-prim(ixC^S,iw) ) &
                      / (xbar(jxC^S,idims) - xbar(ixC^S,idims)) &
                       * dxdim

          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter_in,ldw,rdw)
          primL(ixL^S,iw)=prim(ixL^S,iw)+ldw(ixL^S)*(xi(ixL^S,idims)-xbar(ixL^S,idims))/dxdim
          primR(ixR^S,iw)=prim(jxR^S,iw)+rdw(jxR^S)*(xi(ixR^S,idims)-xbar(jxR^S,idims))/dxdim

          if (loglimit(iw)) then
             prim(ixCmin^D:jxCmax^D,iw)=10.0d0**prim(ixCmin^D:jxCmax^D,iw)
             primL(ixL^S,iw)=10.0d0**primL(ixL^S,iw)
             primR(ixR^S,iw)=10.0d0**primR(ixR^S,iw)
          end if
       end select

    end do

    ! reconstructing hydro variables, 
    if ( fix_small_values ) then
       call phys_handle_small_values(primL,xi,ixI^L,ixL^L, .False., 'reconstruction L')
       call phys_handle_small_values(primR,xi,ixI^L,ixR^L, .False., 'reconstruction R')
    end if

    ! we need to update rest of the prim variables
    call phys_update_eos(ixI^L,ixL^L,primL)
    call phys_update_eos(ixI^L,ixR^L,primR)
   
    call phys_to_conserved(ixI^L,ixL^L,consL,primL,xi)
    call phys_to_conserved(ixI^L,ixR^L,consR,primR,xi)
  end subroutine reconstruct_LR

  subroutine get_theta(ixI^L,ixC^L,idims,eps,inv_lambda,u,flow,fhigh,theta)
    use mod_global_parameters, only: kr, smalldouble
    integer, intent(in)                             :: ixI^L, ixC^L, idims
    double precision, intent(in)                    :: eps
    double precision, dimension(ixI^S), intent(in)  :: inv_lambda, u
    double precision, dimension(ixI^S), intent(in)  :: flow, fhigh
    double precision, dimension(ixI^S), intent(out) :: theta

    integer                                         :: ixCp^L, ixOp^L
    double precision, dimension(ixI^S)              :: tmp, thp, thm
    double precision, dimension(ixI^S)              :: diff_fdA

    ! Note: here we assume that u( i=0 ) is given
    ixCp^L=ixC^L+kr(idims,^D);
    ixOpmin^D=ixCmin^D; ixOpmax^D=ixCpmax^D;
    
    thm(ixC^S) = 1.0d0
    thp(ixC^S) = 1.0d0
    
    tmp(ixOp^S) = 0.5d0*inv_lambda(ixOp^S)*(u(ixOp^S)-eps)

    diff_fdA(ixC^S) = -flow(ixC^S) + fhigh(ixC^S)
    where (diff_fdA(ixC^S) == 0.0d0)
       diff_fdA(ixC^S) = smalldouble ! avoid flow = fhight case
    end where
    
    where (tmp(ixC^S) < fhigh(ixC^S))
       thm(ixC^S) = tmp(ixC^S) - flow(ixC^S)
       thm(ixC^S) = thm(ixC^S) / (diff_fdA(ixC^S))
    end where
    
    where (tmp(ixCp^S) < -fhigh(ixC^S))
       thp(ixC^S) = - tmp(ixCp^S) - flow(ixC^S)
       thp(ixC^S) = thp(ixC^S) / (diff_fdA(ixC^S))
    end where

    theta(ixC^S) = min(thm(ixC^S),thp(ixC^S))
    theta(ixC^S) = min(max(theta(ixC^S),0.0d0),1.0d0)
  end subroutine get_theta

end module mod_finite_volume
