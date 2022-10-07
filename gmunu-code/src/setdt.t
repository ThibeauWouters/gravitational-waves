!>setdt  - set dt for all levels between levmin and levmax. 
!>         dtpar>0  --> use fixed dtpar for all level
!>         dtpar<=0 --> determine CFL limited timestep 
subroutine setdt()
  use mod_global_parameters
  use mod_physics
  use mod_usr_methods, only: usr_get_dt

  integer :: iigrid, igrid, ncycle, ncycle2, ifile, idim
  double precision :: dtnew, qdtnew, dtmin_mype, factor, dx^D, dxmin^D

  double precision :: dtmax, dxmin, cmax_mype, v(ixG^T)

  if (dtpar<=zero) then
     dtmin_mype=bigdouble
     cmax_mype = zero
     !$OMP PARALLEL DO PRIVATE(igrid,qdtnew,dtnew,dx^D) REDUCTION(min:dtmin_mype) REDUCTION(max:cmax_mype)
     do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
        dtnew=bigdouble
        dx^D=rnode(rpdx^D_,igrid);
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        block=>ps(igrid)

        call getdt_courant(ps(igrid)%prim,ixG^LL,ixM^LL,qdtnew,dx^D,ps(igrid)%xbar)
        dtnew=min(dtnew,qdtnew)

        call phys_get_dt(ps(igrid)%prim,ixG^LL,ixM^LL,qdtnew,dx^D,ps(igrid)%xbar)
        dtnew=min(dtnew,qdtnew)

        if (associated(usr_get_dt)) then
           call usr_get_dt(ps(igrid)%prim,ixG^LL,ixM^LL,qdtnew,dx^D,ps(igrid)%xbar)
        end if

        dtnew          = min(dtnew,qdtnew)
        dtmin_mype     = min(dtmin_mype,dtnew)
     end do
     !$OMP END PARALLEL DO
  else
     dtmin_mype=dtpar
  end if

  if (dtmin_mype<dtmin) then
     write(unitterm,*)"Error: Time step too small!", dtmin_mype
     write(unitterm,*)"on processor:", mype, "at time:", global_time," step:", it
     write(unitterm,*)"Lower limit of time step:", dtmin
     crash=.true.
  end if

  if (slowsteps>it-it_init+1) then
     factor=one-(one-dble(it-it_init+1)/dble(slowsteps))**2
     dtmin_mype=dtmin_mype*factor
  end if

  if(final_dt_reduction)then
     !if (dtmin_mype>time_max-global_time) then
     !   write(unitterm,*)"WARNING final timestep artificially reduced!"
     !   write(unitterm,*)"on processor:", mype, "at time:", global_time," step:", it
     !endif
     if(time_max-global_time<=dtmin) then
        !write(unitterm,*)'Forcing to leave timeloop as time is reached!'
        final_dt_exit=.true.
     endif
     dtmin_mype=min(dtmin_mype,time_max-global_time)
  end if

  if (dtpar<=zero) then
     call MPI_ALLREDUCE(dtmin_mype,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                        icomm,ierrmpi)
  else
     dt=dtmin_mype
  end if

  if(any(dtsave(1:nfile)<bigdouble).or.any(tsave(isavet(1:nfile),1:nfile)<bigdouble))then
     dtmax = minval(ceiling(global_time/dtsave(1:nfile))*dtsave(1:nfile))-global_time
     do ifile=1,nfile
        dtmax = min(tsave(isavet(ifile),ifile)-global_time,dtmax)
     end do
     if(dtmax > smalldouble)then 
       dt=min(dt,dtmax)
     else
       ! dtmax=0 means dtsave is divisible by global_time
       dt=min(dt,minval(dtsave(1:nfile)))
     end if      
  end if

  if(mype==0) then
    if(any(dtsave(1:nfile)<dt)) then
      write(unitterm,*) 'Warning: timesteps: ',dt,' exceeding output intervals ', dtsave(1:nfile)
    endif
  endif   

  ! global Lax-Friedrich finite difference flux splitting needs fastest wave-speed
  ! so does GLM: 
  if(need_global_cmax) call MPI_ALLREDUCE(cmax_mype,cmax_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)

  contains

    !> compute CFL limited dt (for variable time stepping)
    subroutine getdt_courant(w,ixI^L,ixO^L,dtnew,dx^D,x)
      use mod_global_parameters
      use mod_physics, only: phys_get_cmax
      
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in)    :: dx^D
      double precision, intent(inout) :: w(ixI^S,1:nprim), dtnew
      
      integer :: idims
      double precision :: courantmax, dxinv(1:ndim), courantmaxtot, courantmaxtots
      double precision :: cmax(ixI^S), cmaxtot(ixI^S), tmp(ixI^S)

      dtnew=bigdouble
      
      courantmax=zero
      courantmaxtot=zero
      courantmaxtots=zero
      
      ^D&dxinv(^D)=one/dx^D;
      
      cmaxtot(ixO^S)=zero
      
      do idims=1,ndim
        call phys_get_cmax(w,x,ixI^L,ixO^L,idims,cmax)
        if(need_global_cmax) cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
        if(slab_uniform) then
          cmaxtot(ixO^S)=cmaxtot(ixO^S)+cmax(ixO^S)*dxinv(idims)
          courantmax=max(courantmax,maxval(cmax(ixO^S)*dxinv(idims)))
        else
          tmp(ixO^S)=cmax(ixO^S)/block%ds(ixO^S,idims)
          cmaxtot(ixO^S)=cmaxtot(ixO^S)+tmp(ixO^S)
          courantmax=max(courantmax,maxval(tmp(ixO^S)))
        end if
        courantmaxtot=courantmaxtot+courantmax
      end do
      
      select case (typecourant)
      case ('minimum')
         ! courantmax='max(c/dx)'
         if (courantmax>smalldouble)     dtnew=min(dtnew,courantpar/courantmax)
      case ('summax')
         ! courantmaxtot='summed max(c/dx)'
         if (courantmaxtot>smalldouble)  dtnew=min(dtnew,courantpar/courantmaxtot)
      case ('maxsum')
         ! courantmaxtots='max(summed c/dx)'
         courantmaxtots=max(courantmaxtots,maxval(cmaxtot(ixO^S)))
         if (courantmaxtots>smalldouble) dtnew=min(dtnew,courantpar/courantmaxtots)
      case default
         write(unitterm,*)'Unknown typecourant=',typecourant
         call mpistop("Error from getdt_courant: no such typecourant!")
      end select
      
    end subroutine getdt_courant

end subroutine setdt
