!> fill ghost cells at a physical boundary
subroutine bc_phys(iside,idims,time,qdt,s,ixG^L,ixB^L)
  use mod_usr_methods, only: usr_special_bc
  use mod_bc_data, only: bc_data_set
  use mod_global_parameters
  use mod_physics

  integer, intent(in) :: iside, idims, ixG^L,ixB^L
  double precision, intent(in) :: time,qdt
  type(state), intent(inout) :: s
  double precision :: wtmp(ixG^S,1:nprim)

  integer :: idir, is
  integer :: ixOs^L,hxO^L,jxO^L
  double precision :: Q(ixG^S),Qp(ixG^S) 
  integer :: iw, iB, ix^D, ixO^L, ixM^L, nghostcellsi,iib^D
  logical  :: isphysbound
  integer :: vel_index

  associate(x=>s%x,w=>s%prim,ws=>s%conss)
  select case (idims)
  {case (^D)
     if (iside==2) then
        ! maximal boundary
        iB=2*^D
        ixOmin^DD=ixBmax^D+1-nghostcells^D%ixOmin^DD=ixBmin^DD;
        ixOmax^DD=ixBmax^DD;
        ! cont/symm/asymm types
        do iw=1,nprim
           select case (typeboundary(iw,iB))
           case ("symm")
              w(ixO^S,iw) = w(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S,iw)
           case ("asymm")
              w(ixO^S,iw) =-w(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S,iw)
           case ("cont")
              do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
              end do
           case("noinflow")
              ! noinflow applies only on velocities, for other vectors, use cont
              ! Note that although ^D = ndim is not ndir, in case of ndir > ndim, 
              ! the lastest velocity has no corresponding boundary
              ! so apply cont should be enough.
              if ( allocated(veloc) ) then
                 vel_index = veloc(^D)
              else if ( allocated(W_vel) ) then
                 vel_index = W_vel(^D)
              else
                 call mpistop("something got wrong, no veloc is allocated")
              end if

              if ( iw==vel_index )then
                do ix^D=ixOmin^D,ixOmax^D
                    w(ix^D^D%ixO^S,iw) = max(w(ixOmin^D-1^D%ixO^S,iw),zero)
                end do
              else
                do ix^D=ixOmin^D,ixOmax^D
                    w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
                end do
              end if
           case ("special", "bc_data")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("character")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("aperiodic")
              !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixO^S,iw) = - w(ixO^S,iw)
           case ("periodic")
  !            call mpistop("periodic bc info should come from neighbors")
           case default
              write (unitterm,*) "Undefined boundarytype ",typeboundary(iw,iB), &
                 "for variable iw=",iw," and side iB=",iB
           end select
        end do
        if(stagger_grid) then
          do idir=1,nconss
          ! At this stage, extrapolation is applied only to the tangential components
            if(idir==^D) cycle 
            ixOsmax^DD=ixOmax^DD;
            ixOsmin^DD=ixOmin^DD-kr(^DD,idir);
            select case(typeboundary(Bvec(idir),iB))
            case ("symm")
              ws(ixOs^S,idir) = ws(ixOsmin^D-1:ixOsmin^D-nghostcells:-1^D%ixOs^S,idir)
            case ("asymm")
              ws(ixOs^S,idir) =-ws(ixOsmin^D-1:ixOsmin^D-nghostcells:-1^D%ixOs^S,idir)
            case ("cont","noinflow")
              do ix^D=ixOsmin^D,ixOsmax^D
                 ws(ix^D^D%ixOs^S,idir) = ws(ixOsmin^D-1^D%ixOs^S,idir)
              end do
            case ("periodic")
            case ("special")
               ! skip it here, do AFTER all normal type boundaries are set
            case default
              write (unitterm,*) "Undefined boundarytype ",typeboundary(iw,iB), &
                 "for variable iw=",iw," and side iB=",iB
            end select
          end do
          ! Now that the tangential components are set,
          ! fill the normal components using a prescription for the divergence.
          ! This prescription is given by the typeB for the normal component.
          do idir=1,nconss
            ! Consider only normal direction
            if (idir/=^D) cycle
            ixOs^L=ixO^L;
            hxO^L=ixO^L-nghostcells*kr(^DD,^D);
            ! Calculate divergence and partial divergence
            call phys_get_divb(ixG^L,hxO^L,s,Q)
            select case(typeboundary(Bvec(idir),iB))
            case("symm")
              ws(ixOs^S,idir)=zero
              do ix^D=0,nghostcells-1
                call phys_get_divb(ixG^L,ixO^L,s,Qp)
                ws(ixOsmin^D+ix^D^D%ixOs^S,idir)=&
                  (Q(hxOmax^D-ix^D^D%hxO^S)*s%dvolume(hxOmax^D-ix^D^D%hxO^S)&
                 -Qp(ixOmin^D+ix^D^D%ixO^S)*s%dvolume(ixOmin^D+ix^D^D%ixO^S))&
                  /s%surfaceC(ixOsmin^D+ix^D^D%ixOs^S,^D)
              end do
            case("asymm")
              ws(ixOs^S,idir)=zero
              do ix^D=0,nghostcells-1
                call phys_get_divb(ixG^L,ixO^L,s,Qp)
                ws(ixOsmin^D+ix^D^D%ixOs^S,idir)=&
                 (-Q(hxOmax^D-ix^D^D%hxO^S)*s%dvolume(hxOmax^D-ix^D^D%hxO^S)&
                 -Qp(ixOmin^D+ix^D^D%ixO^S)*s%dvolume(ixOmin^D+ix^D^D%ixO^S))&
                  /s%surfaceC(ixOsmin^D+ix^D^D%ixOs^S,^D)
              end do
            case("cont")
              ws(ixOs^S,idir)=zero
              do ix^D=0,nghostcells-1
                call phys_get_divb(ixG^L,ixO^L,s,Qp)
                ws(ixOsmin^D+ix^D^D%ixOs^S,idir)=&
                  (Q(hxOmax^D^D%hxO^S)*s%dvolume(hxOmax^D^D%hxO^S)&
                 -Qp(ixOmin^D+ix^D^D%ixO^S)*s%dvolume(ixOmin^D+ix^D^D%ixO^S))&
                  /s%surfaceC(ixOsmin^D+ix^D^D%ixOs^S,^D)
              end do
            case("periodic")
            end select
          end do
          ! Fill cell averages
          !call phys_face_to_center(ixO^L,s)
        end if
     else
        ! minimal boundary
        iB=2*^D-1
        ixOmin^DD=ixBmin^DD;
        ixOmax^DD=ixBmin^D-1+nghostcells^D%ixOmax^DD=ixBmax^DD;
        ! cont/symm/asymm types
        do iw=1,nprim
           select case (typeboundary(iw,iB))
           case ("symm")
              w(ixO^S,iw) = w(ixOmax^D+nghostcells:ixOmax^D+1:-1^D%ixO^S,iw)
           case ("asymm")
              w(ixO^S,iw) =-w(ixOmax^D+nghostcells:ixOmax^D+1:-1^D%ixO^S,iw)
           case ("cont")
              do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
              end do
           case("noinflow")
              if ( allocated(veloc) ) then
                 vel_index = veloc(^D)
              else if ( allocated(W_vel) ) then
                 vel_index = W_vel(^D)
              else
                 call mpistop("something got wrong, no veloc is allocated")
              end if

              if ( iw==vel_index )then
                 do ix^D=ixOmin^D,ixOmax^D
                   w(ix^D^D%ixO^S,iw) = min(w(ixOmax^D+1^D%ixO^S,iw),zero)
                 end do
              else
                 do ix^D=ixOmin^D,ixOmax^D
                   w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
                 end do
              end if
           case ("special", "bc_data")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("character")
              ! skip it here, do AFTER all normal type boundaries are set
           case ("aperiodic")
              !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixO^S,iw) = - w(ixO^S,iw)
           case ("periodic")
  !            call mpistop("periodic bc info should come from neighbors")
           case default
              write (unitterm,*) "Undefined boundarytype ",typeboundary(iw,iB), &
                 "for variable iw=",iw," and side iB=",iB
           end select
        end do
        if(stagger_grid) then
          do idir=1,nconss
          ! At this stage, extrapolation is applied only to the tangential components
            if(idir==^D) cycle 
            ixOsmax^DD=ixOmax^DD;
            ixOsmin^DD=ixOmin^DD-kr(^DD,idir);
            select case(typeboundary(Bvec(idir),iB))
            case ("symm")
              ws(ixOs^S,idir) = ws(ixOsmax^D+nghostcells:ixOsmax^D+1:-1^D%ixOs^S,idir)
            case ("asymm")
              ws(ixOs^S,idir) =-ws(ixOsmax^D+nghostcells:ixOsmax^D+1:-1^D%ixOs^S,idir)
            case ("cont","noinflow")
              do ix^D=ixOsmin^D,ixOsmax^D
                 ws(ix^D^D%ixOs^S,idir) = ws(ixOsmax^D+1^D%ixOs^S,idir)
              end do
            case ("periodic")
            case ("special")
               ! skip it here, do AFTER all normal type boundaries are set
            case default
              write (unitterm,*) "Undefined boundarytype ",typeboundary(iw,iB), &
                 "for variable iw=",iw," and side iB=",iB
            end select
          end do
          ! Now that the tangential components are set,
          ! fill the normal components using a prescription for the divergence.
          ! This prescription is given by the typeB for the normal component.
          do idir=1,nconss
            ! Consider only normal direction
            if (idir/=^D) cycle
            ixOs^L=ixO^L-kr(^DD,^D);
            jxO^L=ixO^L+nghostcells*kr(^DD,^D);
            ! Calculate divergence and partial divergence
            call phys_get_divb(ixG^L,jxO^L,s,Q)
            select case(typeboundary(Bvec(idir),iB))
            case("symm")
              ws(ixOs^S,idir)=zero
              do ix^D=0,nghostcells-1
                call phys_get_divb(ixG^L,ixO^L,s,Qp)
                ws(ixOsmax^D-ix^D^D%ixOs^S,idir)=&
                 -(Q(jxOmin^D+ix^D^D%jxO^S)*s%dvolume(jxOmin^D+ix^D^D%jxO^S)&
                 -Qp(ixOmax^D-ix^D^D%ixO^S)*s%dvolume(ixOmax^D-ix^D^D%ixO^S))&
                  /s%surfaceC(ixOsmax^D-ix^D^D%ixOs^S,^D)
              end do
            case("asymm")
              ws(ixOs^S,idir)=zero
              do ix^D=0,nghostcells-1
                call phys_get_divb(ixG^L,ixO^L,s,Qp)
                ws(ixOsmax^D-ix^D^D%ixOs^S,idir)=&
                 -(-Q(jxOmin^D+ix^D^D%jxO^S)*s%dvolume(jxOmin^D+ix^D^D%jxO^S)&
                 -Qp(ixOmax^D-ix^D^D%ixO^S)*s%dvolume(ixOmax^D-ix^D^D%ixO^S))&
                  /s%surfaceC(ixOsmax^D-ix^D^D%ixOs^S,^D)
              end do
            case("cont")
              ws(ixOs^S,idir)=zero
              do ix^D=0,nghostcells-1
                call phys_get_divb(ixG^L,ixO^L,s,Qp)
                ws(ixOsmax^D-ix^D^D%ixOs^S,idir)=&
                 -(Q(jxOmin^D^D%jxO^S)*s%dvolume(jxOmin^D^D%jxO^S)&
                 -Qp(ixOmax^D-ix^D^D%ixO^S)*s%dvolume(ixOmax^D-ix^D^D%ixO^S))&
                  /s%surfaceC(ixOsmax^D-ix^D^D%ixOs^S,^D)
              end do
            case("periodic")
            end select
          end do
          ! Fill cell averages
          !call phys_face_to_center(ixO^L,s)
        end if
     end if \}
  end select

  ! do user defined special boundary conditions
  if (any(typeboundary(1:nprim,iB)=="special")) then
     if (.not. associated(usr_special_bc)) &
          call mpistop("usr_special_bc not defined")
     call usr_special_bc(time,ixG^L,ixO^L,iB,w,x)
  end if

  ! fill boundary conditions from external data vtk files
  if (any(typeboundary(1:nprim,iB)=="bc_data")) then
     call bc_data_set(time,ixG^L,ixO^L,iB,w,x)
  end if

  {#IFDEF EVOLVINGBOUNDARY
  if (any(typeboundary(1:nprim,iB)=="character")) then
    ixM^L=ixM^LL;
    if(ixGmax1==ixGhi1) then
      nghostcellsi=nghostcells
    else
      nghostcellsi=ceiling(nghostcells*0.5d0)
    end if
    select case (idims)
    {case (^D)
       if (iside==2) then
          ! maximal boundary
          ixOmin^DD=ixGmax^D+1-nghostcellsi^D%ixOmin^DD=ixBmin^DD;
          ixOmax^DD=ixBmax^DD;
          if(all(w(ixO^S,1:nprim)==0.d0)) then
            do ix^D=ixOmin^D,ixOmax^D
               w(ix^D^D%ixO^S,1:nprim) = w(ixOmin^D-1^D%ixO^S,1:nprim)
            end do
          end if
          if(qdt>0.d0.and.ixGmax^D==ixGhi^D) then
            ixOmin^DD=ixOmin^D^D%ixOmin^DD=ixMmin^DD;
            ixOmax^DD=ixOmax^D^D%ixOmax^DD=ixMmax^DD;
            wtmp(ixG^S,1:nprim)=pso(block%igrid)%prim(ixG^S,1:nprim)
            call characteristic_project(idims,iside,ixG^L,ixO^L,wtmp,x,dxlevel,qdt)
            w(ixO^S,1:nprim)=wtmp(ixO^S,1:nprim)
          end if
       else
          ! minimal boundary
          ixOmin^DD=ixBmin^DD;
          ixOmax^DD=ixGmin^D-1+nghostcellsi^D%ixOmax^DD=ixBmax^DD;
          if(all(w(ixO^S,1:nprim)==0.d0)) then
            do ix^D=ixOmin^D,ixOmax^D
               w(ix^D^D%ixO^S,1:nprim) = w(ixOmax^D+1^D%ixO^S,1:nprim)
            end do
          end if
          if(qdt>0.d0.and.ixGmax^D==ixGhi^D) then
            ixOmin^DD=ixOmin^D^D%ixOmin^DD=ixMmin^DD;
            ixOmax^DD=ixOmax^D^D%ixOmax^DD=ixMmax^DD;
            wtmp(ixG^S,1:nprim)=pso(block%igrid)%prim(ixG^S,1:nprim)
            call characteristic_project(idims,iside,ixG^L,ixO^L,wtmp,x,dxlevel,qdt)
            w(ixO^S,1:nprim)=wtmp(ixO^S,1:nprim)
          end if
       end if \}
    end select
    if(ixGmax1==ixGhi1) then
      call identifyphysbound(block%igrid,isphysbound,iib^D)   
      if(iib1==-1.and.iib2==-1) then
        do ix2=nghostcells,1,-1 
          do ix1=nghostcells,1,-1 
            w(ix^D,1:nprim)=(w(ix1+1,ix2+1,1:nprim)+w(ix1+1,ix2,1:nprim)+w(ix1,ix2+1,1:nprim))/3.d0
          end do
        end do
      end if
      if(iib1== 1.and.iib2==-1) then
        do ix2=nghostcells,1,-1 
          do ix1=ixMmax1+1,ixGmax1
            w(ix^D,1:nprim)=(w(ix1-1,ix2+1,1:nprim)+w(ix1-1,ix2,1:nprim)+w(ix1,ix2+1,1:nprim))/3.d0
          end do
        end do
      end if
      if(iib1==-1.and.iib2== 1) then
        do ix2=ixMmax2+1,ixGmax2
          do ix1=nghostcells,1,-1 
            w(ix^D,1:nprim)=(w(ix1+1,ix2-1,1:nprim)+w(ix1+1,ix2,1:nprim)+w(ix1,ix2-1,1:nprim))/3.d0
          end do
        end do
      end if
      if(iib1== 1.and.iib2== 1) then
        do ix2=ixMmax2+1,ixGmax2
          do ix1=ixMmax1+1,ixGmax1
            w(ix^D,1:nprim)=(w(ix1-1,ix2-1,1:nprim)+w(ix1-1,ix2,1:nprim)+w(ix1,ix2-1,1:nprim))/3.d0
          end do
        end do
      end if
    end if
  end if
  }
  !end do

  ! Fill cell averages
  if (stagger_grid) call phys_face_to_center(ixO^L,s)

  end associate
end subroutine bc_phys

!> fill inner boundary values
subroutine getintbc(time,ixG^L)
  use mod_usr_methods, only: usr_internal_bc
  use mod_global_parameters

  double precision, intent(in)   :: time
  integer, intent(in)            :: ixG^L

  integer :: iigrid, igrid, ixO^L

  ixO^L=ixG^L^LSUBnghostcells;

  !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
     block=>ps(igrid)
     ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

     if (associated(usr_internal_bc)) then
        call usr_internal_bc(node(plevel_,igrid),time,ixG^L,ixO^L,ps(igrid)%prim,ps(igrid)%x)
     end if
  end do
  !$OMP END PARALLEL DO

end subroutine getintbc
