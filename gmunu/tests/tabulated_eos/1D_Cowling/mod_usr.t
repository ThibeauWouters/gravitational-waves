module mod_usr
  use mod_physics
  use mod_grhd
  use mod_eos_tabulated
  use mod_XNS

  implicit none

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    character*128  ::  profile_path = './BU0_1D_640'
    integer  ::  Nth, Nr   ! number of grid of the profile
    double precision :: xrho, xtemp, xye, xenr


    usr_init_one_grid => ns_init_one_grid
    usr_before_main_loop => ns_before_main_loop
    usr_print_log => printlog

    call set_coordinate_system("spherical")
    call grhd_activate()
    call eos_tabulated_activate()
    ! use atmosphere
    call eos_atmo_activate()

    call initialize_gmunu()

    call mod_XNS_read_profile(profile_path)
    if ( prof%Nth == 1) then  
    ! This is a 1D profile.
       if ( ndir .ne. 1 ) stop "! Dimension = 2 but the profile is 1D."
    else
    ! This is a 2D profile.
       if ( ndir < 2 ) stop "! Dimension = 1 but the profile is 2D."
    endif

    ! setup atmosphere density from the profile
    call eos_initialize_atmo(maxval(prof%rho))

  xrho = 10.0d0**1.474994d1 * 1.61930347d-18
  xtemp = 63.0d0
  xye = 0.2660725d0
  xenr = 0.0d0

  call eos_get_eps_one_grid(0.0d0, xrho, xenr, xtemp, xye)

  write(6,*) "######################################"
  write(6,"(1P10E15.6)") xrho,xtemp,xye
  write(6,"(1P10E15.6)") xenr
  write(6,*) "######################################"

  call eos_get_temp_one_grid(xrho, xenr, xtemp, xye)

  write(6,*) "######################################"
  write(6,"(1P10E15.6)") xrho,xtemp,xye
  write(6,"(1P10E15.6)") xenr
  write(6,*) "######################################"

  write(6,*) "######################################"
  write(6,*) "######################################"
  write(6,*) "######################################"
  xenr = xenr + 100.0d0 * 1.11265006d-21
  xtemp = 120.0d0
  call eos_get_temp_one_grid(xrho, xenr, xtemp, xye)

  write(6,*) "######################################"
  write(6,"(1P10E15.6)") xrho,xtemp,xye
  write(6,"(1P10E15.6)") xenr
  write(6,*) "######################################"

    stop

  end subroutine usr_init

  ! Initialize one grid
  subroutine ns_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    use mod_global_parameters
    use mod_grhd
    use mod_eos

    integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nprim)

    integer  ::  ix1

    w = 0.0d0

    do ix1 =ixmin1, ixmax1 
       call mod_XNS_map(w(ix1,rho_),x(ix1,r_),prof%rho(:,1),prof%radius,prof%Nr)
       call mod_XNS_map(w(ix1,press_),x(ix1,r_),prof%press(:,1),prof%radius,prof%Nr)
       call mod_XNS_map(w(ix1,alp_),x(ix1,r_),prof%alp(:,1),prof%radius,prof%Nr)
       call mod_XNS_map(w(ix1,psi_),x(ix1,r_),prof%psi(:,1),prof%radius,prof%Nr)
    end do
    w(ixmin1:ixmax1,beta(1))   = 0.0d0


    ! fill atmosphere
    ! here we assume that the smallest rho in the profile is the atmosphere
    !where (w(ixmin1:ixmax1,rho_) <= minval(prof%rho) )
    !   w(ixmin1:ixmax1,rho_)   = 0.0d0
    !end where
    where (w(ixmin1:ixmax1,rho_) <= small_rho_thr )
       w(ixmin1:ixmax1,rho_)   = small_rho
       w(ixmin1:ixmax1,veloc(1)) = 0.0d0
       w(ixmin1:ixmax1,press_)     = small_press
       w(ixmin1:ixmax1,eps_)     = small_eps
    end where

    do ix1 =ixmin1, ixmax1 
       ! get_eps
       call eos_get_eps_one_grid(w(ix1, press_),w(ix1, rho_),w(ix1, eps_))
       call eos_get_pressure_one_grid(w(ix1, press_),w(ix1, rho_),w(ix1, eps_))
    end do

  end subroutine ns_init_one_grid
  
  !> deallocate profile varibles to save memories
  subroutine ns_before_main_loop()
    deallocate(prof%radius)
    deallocate(prof%theta)
    deallocate(prof%rho,prof%press,prof%v3)
    deallocate(prof%psi,prof%alp,prof%beta3)
  end subroutine ns_before_main_loop

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)
    use mod_eos, only: small_rho, small_press,small_eps
    ! special boundary types, user defined
    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nprim)
    integer :: ixI^L

    ixImin^DD=ixBmin^DD;
    ixImax^DD=ixBmin^D-1+nghostcells^D%ixImax^DD=ixBmax^DD;
    ! fix rho = small_rho at the outter boundary:
    w(ixImin1:ixImax1,rho_) = small_rho
    w(ixImin1:ixImax1,press_)   = small_press
    w(ixImin1:ixImax1,eps_)   = small_eps
    !do ix2=ixImin2,ixImax2
    !   w(ixImin1:ixImax1,ix2,rho_) = w(ixImin1:ixImax1,ixImax2+1,rho_) 
    !   w(ixImin1:ixImax1,ix2,press_)   = w(ixImin1:ixImax1,ixImax2+1,press_) 
    !end do

    !where(dabs(x(ixI^S,1))<1.0d0)
    !   w(ixI^S,rho_)   = 1.0d0
    !   w(ixI^S,veloc(1))=zero
    !   w(ixI^S,veloc(2)) = 0.995d0
    !   w(ixI^S,press_)     = 5.0d-3
    !else where
       ! Reflective:
       !   w(ixI^S,rho_) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,rho_) 
       !   w(ixI^S,e_) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,e_) 
       !   w(ixI^S,veloc(1)) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,veloc(1))
       !   w(ixI^S,veloc(2)) =-w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,veloc(2))
    !end where

  end subroutine specialbound_usr

  subroutine printlog()

    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters

    logical              :: fileopen
    integer              :: i, iw, level
    double precision     :: wmean(1:nprim), total_volume
    double precision     :: volume_coverage(refine_max_level)
    integer              :: nx^D, nc, ncells, dit
    double precision     :: dtTimeLast, now, cellupdatesPerSecond
    double precision     :: activeBlocksPerCore, wctPerCodeTime, timeToFinish
    character(len=40)    :: fmt_string
    character(len=80)    :: filename
    character(len=2048)  :: line
    logical, save        :: opened  = .false.
    integer              :: amode, istatus(MPI_STATUS_SIZE)
    integer, parameter   :: my_unit = 20

    ! Compute the volume-average of w**1 = w
    call get_volume_average(1, wmean, total_volume)

    ! Compute the volume coverage
    call get_volume_coverage(volume_coverage)
    wmean(1:nprim) = ps(1)%prim(ixMlo^D,1:nprim)

    if (mype == 0) then

       ! On first entry, open the file and generate the header
       if (.not. opened) then

          filename = trim(base_filename) // ".log"

          ! Delete the log when not doing a restart run
          if (restart_from_file == undefined) then
             open(unit=my_unit,file=trim(filename),status='replace')
             close(my_unit, status='delete')
          end if

          amode    = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
          amode    = ior(amode,MPI_MODE_APPEND)

          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, &
               MPI_INFO_NULL, log_fh, ierrmpi)

          opened = .true.

          ! Start of file headern
          line = "it global_time dt"
          do level=1,nprim
             i = len_trim(line) + 2
             write(line(i:),"(a,a)") trim(prim_names(level)), " "
          end do

          ! Only write header if not restarting
          if (restart_from_file == undefined) then
            call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'), &
                 len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
          end if
       end if

       ! Construct the line to be added to the log

       fmt_string = '(' // fmt_i // ',2' // fmt_r // ')'
       write(line, fmt_string) it, global_time, dt
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', nprim, fmt_r // ')'
       write(line(i:), fmt_string) wmean(1:nprim)
       i = len_trim(line) + 2

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

  end subroutine printlog

end module mod_usr
