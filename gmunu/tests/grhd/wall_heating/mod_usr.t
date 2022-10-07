module mod_usr
  use mod_physics
  use mod_grhd
  use mod_eos_idealgas

  implicit none
  private
  ! alpha = 0 for Cartesian
  ! alpha = 1 for cylindrical
  ! alpha = 2 for spherical
  integer, parameter :: alpha = 1
 
  double precision, parameter :: W_in = 1.0d1
  double precision, parameter :: rho_0 = 1.0d0
  double precision            ::  v_in, v_sh, sigma

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid
    usr_special_bc       => specialbound_usr
    usr_print_log        => printlog
    
    select case (alpha)
    case (0)
       call set_coordinate_system("Cartesian_1D")
    case (1)
       call set_coordinate_system("cylindrical")
    case (2)
       call set_coordinate_system("spherical")
    case default
       call mpistop("alpha can only be 0, 1 or 2.")
    end select

    call grhd_activate()
    call eos_idealgas_activate()

    v_in = dsqrt(1.0d0 - 1.0d0 / W_in**2)
    ! shock velocity
    v_sh = ( eos_gamma - 1.0d0 ) * W_in * dabs(v_in) / ( W_in + 1.0d0 )
    ! compression ratio
    sigma = (eos_gamma+1.0d0)/(eos_gamma-1.0d0) &
            +(eos_gamma*(W_in-1.0d0))/(eos_gamma-1.0d0)


    ! setup atmosphere density from the profile
    !call eos_initialize_atmo(1.0d0)
  end subroutine usr_init

  ! Initialize one grid
  subroutine rm_init_one_grid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics
    use mod_grhd
    use mod_eos

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim) ! primitive variables

    integer                         ::  ix^D, idir

    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0

    where (dabs(x(ixO^S,1)) >= v_sh * global_time)
       w(ixO^S,rho_)       = rho_0 * ( 1.0d0 + dabs(v_in / x(ixO^S,1)) * global_time )**alpha 
       w(ixO^S,W_vel(1))   = - sign(1.0d0,x(ixO^S,1)) * dsqrt( W_in**2 - 1.0d0)
    else where
       w(ixO^S,rho_)       = sigma * rho_0 * ( 1.0d0 + dabs(v_in) / v_sh )**alpha
       w(ixO^S,W_vel(1))   = 0.0d0
    end where
    w(ixO^S,press_)     = 1.0d-7 * W_in * (eos_gamma-1.0d0) * w(ixO^S, rho_) 

    ! reset eps
    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
    ! update rest of the prim
    call phys_update_eos(ixI^L, ixO^L, w)

  end subroutine rm_init_one_grid

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)
    use mod_eos

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nprim)
    integer :: ixO^L, ix^D

    if (iB == ismax^D) then
       ixOmin^DD=ixBmax^D+1-nghostcells^D%ixOmin^DD=ixBmin^DD;
       ixOmax^DD=ixBmax^DD;
    else
       ixOmin^DD=ixBmin^DD;
       ixOmax^DD=ixBmin^D-1+nghostcells^D%ixOmax^DD=ixBmax^DD;
    end if
   
    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0

    where (dabs(x(ixO^S,1)) >= v_sh * global_time)
       w(ixO^S,rho_)       = rho_0 * ( 1.0d0 + dabs(v_in / x(ixO^S,1) ) * global_time )**alpha 
       w(ixO^S,W_vel(1))   = - sign(1.0d0, x(ixO^S,1)) * dsqrt( W_in**2 - 1.0d0)
    else where
       w(ixO^S,rho_)       = sigma * rho_0 * ( 1.0d0 + dabs(v_in) / v_sh )**alpha
       w(ixO^S,W_vel(1))   = 0.0d0
    end where

    where (x(ixG^S,1) < v_sh * global_time)
       !w(ixG^S,W_vel(1))   = 0.0d0
       !w(ixG^S,rho_)       = sigma * rho_0 * ( 1.0d0 + dabs(v_in) / v_sh )**alpha
    end where

    w(ixO^S,press_)     = 1.0d-7 * W_in * (eos_gamma-1.0d0) * w(ixO^S, rho_) 
    ! dont forget to update these primitive variables
    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
    ! update rest of the prim
    call phys_update_eos(ixG^L, ixO^L, w)

  end subroutine specialbound_usr

  subroutine printlog()
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    logical              :: fileopen
    integer              :: i, iw, level
    double precision, allocatable     :: send_buffer(:), recive_buffer(:)
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

    integer, parameter   :: total_volume = 1
    integer, parameter   :: total_L1dV = 2
    integer, parameter   :: total_L2dV = 3
    integer, parameter   :: total_soldV = 4
    integer, parameter   :: total_sol2dV = 5
    integer, parameter   :: n_var = 5

    double precision     :: rho_sol(ixG^T), rho_diff(ixG^T)
    integer              :: igrid, iigrid
    integer              :: ix^D

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! calculate total volume
       send_buffer(total_volume) = send_buffer(total_volume) + sum(ps(igrid)%dvolume(ixM^T))

       ! calculate the difference from analytic solution of rho
       rho_diff = 0.0d0
       rho_sol = 0.0d0
       {do ix^D = ixM^LL \}
          rho_sol(ix^D) = analytic_sol_of_rho(ps(igrid)%x(ix^D,1),global_time)
       {enddo^D&\}
       ! calculate total soldV
       send_buffer(total_soldV) = send_buffer(total_soldV) + &
          sum( rho_sol(ixM^T) * ps(igrid)%dvolume(ixM^T))
       ! calculate total sol2dV
       send_buffer(total_sol2dV) = send_buffer(total_sol2dV) + &
          sum( rho_sol(ixM^T)**2 * ps(igrid)%dvolume(ixM^T))

       rho_diff(ixM^T) = dabs( ps(igrid)%prim(ixM^T,rho_) - rho_sol(ixM^T) )
       ! calculate total L1dV
       send_buffer(total_L1dV) = send_buffer(total_L1dV) + &
          sum( rho_diff(ixM^T) * ps(igrid)%dvolume(ixM^T))
       ! calculate total L2dV
       send_buffer(total_L2dV) = send_buffer(total_L2dV) + &
          sum( rho_diff(ixM^T)**2 * ps(igrid)%dvolume(ixM^T))
    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    ! Volume average of L1 and L2 norm
    recive_buffer(total_L1dV) = (recive_buffer(total_L1dV) / recive_buffer(total_soldV))
    recive_buffer(total_L2dV) = dsqrt(recive_buffer(total_L2dV) / recive_buffer(total_sol2dV))

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
          line = "it global_time dt total_V L1 L2"

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

       write(fmt_string, '(a,i0,a)') '(', n_var, fmt_r // ')'
       write(line(i:), fmt_string) recive_buffer(1:n_var)
       i = len_trim(line) + 2

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

    deallocate(send_buffer)
    deallocate(recive_buffer)

  end subroutine printlog
 
  function analytic_sol_of_rho(x,t)
     implicit none
     double precision :: analytic_sol_of_rho
     double precision, intent(in) :: x,t
     if ( x >= v_sh * t) then
        analytic_sol_of_rho = rho_0 * ( 1.0d0 + dabs(v_in) / x * t )**alpha 
     else
        analytic_sol_of_rho = rho_0 * ( 1.0d0 + dabs(v_in) / v_sh )**alpha 
     end if
  end function analytic_sol_of_rho

end module mod_usr
