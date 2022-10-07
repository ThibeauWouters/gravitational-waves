module mod_usr
  use mod_physics
  use mod_gremhd
  use mod_eos_idealgas

  implicit none
  private

  ! aux output variables
  integer, protected  :: divB_ = -1
  integer, protected  :: divE_ = -1

  logical, parameter                     :: initialize_mag_field = .true.
  double precision, parameter            :: resistivity = 1.0d-3

  double precision, parameter            :: q0 = 7.0d-1
  double precision, parameter            :: p0 = 1.0d-1

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid
    usr_special_bc       => specialbound_usr
    usr_improve_initial_condition => rm_improve_initial_condition
    usr_get_resistivity => my_get_resistivity
    usr_process_adv_grid => update_aux_vars
    usr_print_log => printlog

    call set_coordinate_system("Cartesian_2.5D")
    call gremhd_activate()
    call eos_idealgas_activate()

    divB_ = var_set_auxvar('divB')
    divE_ = var_set_auxvar('divE')

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm_init_one_grid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics
    use mod_gremhd
    use mod_eos

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim) ! primitive variables
    double precision                :: v(ixI^S,1:ndir)
    double precision                :: lfac(ixI^S)
    double precision                :: r2(ixI^S)!, cos_theta(ixI^S), sin_theta(ixI^S)
    double precision                :: Er(ixI^S), vphi(ixI^S)

    integer  ::  ix^D, idir, max_loc(1:ndir)

    w = 0.0d0
    v = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0
    w(ixO^S,beta(3)) = 0.0d0

    r2(ixO^S) = ( x(ixO^S,1)**2 + x(ixO^S,2)**2 )
    !cos_theta(ixO^S) = x(ixO^S,1) / r(ixO^S)
    !sin_theta(ixO^S) = x(ixO^S,2) / r(ixO^S)
    !Er(ixO^S) = 0.5d0 * q0 * r(ixO^S) / (r(ixO^S)**2 + 1.0d0)
    !vphi(ixO^S) = - 0.5d0 * q0 * r(ixO^S) / dsqrt( (r(ixO^S)**2 + 1.0d0)**2 - 0.25d0 * q0**2 )
    !w(ixO^S, Evec(1)) = Er(ixO^S) * cos_theta(ixO^S)
    !w(ixO^S, Evec(2)) = Er(ixO^S) * sin_theta(ixO^S)
    !v(ixO^S, 1) = - vphi(ixO^S) * sin_theta(ixO^S)
    !v(ixO^S, 2) =   vphi(ixO^S) * cos_theta(ixO^S)

    w(ixO^S, Evec(1)) = 0.5d0 * q0 * x(ixO^S,1) / (r2(ixO^S) + 1.0d0)
    w(ixO^S, Evec(2)) = 0.5d0 * q0 * x(ixO^S,2) / (r2(ixO^S) + 1.0d0)
    v(ixO^S, 1) =   0.5d0 * q0 * x(ixO^S,2) / dsqrt( (r2(ixO^S) + 1.0d0)**2 - 0.25d0 * q0**2 )
    v(ixO^S, 2) = - 0.5d0 * q0 * x(ixO^S,1) / dsqrt( (r2(ixO^S) + 1.0d0)**2 - 0.25d0 * q0**2 )

    w(ixO^S,rho_) = 1.0d0
    w(ixO^S,press_) = - w(ixO^S,rho_) * ( eos_gamma - 1.0d0 ) / eos_gamma &
                      + ( p0 + w(ixO^S,rho_) * ( eos_gamma - 1.0d0 ) / eos_gamma ) * &
                        ( (4.0d0*r2(ixO^S) + 4.0d0 - q0**2) / (r2(ixO^S) + 1.0d0) / (4.0d0-q0**2) )**(0.5d0*eos_gamma/(eos_gamma-1.0d0))

    w(ixO^S, Bvec(3))    = dsqrt( (r2(ixO^S) + 1.0d0)**2 - 0.25d0 * q0**2 ) / (r2(ixO^S) + 1.0d0)

    ! get W
    lfac(ixO^S) = 1.0d0
    ! calculate 1 - v^2 first
    do idir = 1, ndir
       lfac(ixO^S) = lfac(ixO^S) - v(ixO^S, idir)**2
    end do
    ! Calculate the Lorentz factor from velocity
    lfac(ixO^S) = dsqrt( 1.0d0 / lfac(ixO^S) )

    ! set veloc -> W_vel
    w(ixO^S,W_vel(1))   = v(ixO^S, 1) * lfac(ixO^S)
    w(ixO^S,W_vel(2))   = v(ixO^S, 2) * lfac(ixO^S)
    w(ixO^S,W_vel(3))   = v(ixO^S, 3) * lfac(ixO^S)

    w(ixO^S, divB_) = 0.0d0
    w(ixO^S, divE_) = q0 / (r2(ixO^S) + 1.0d0)**2

    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
    call phys_update_eos(ixI^L, ixO^L, w)

  end subroutine rm_init_one_grid
  
  !> before the main loop, we improve the initial data here
  subroutine rm_improve_initial_condition()
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_physics
    use mod_eos

    integer              :: ix^D
    integer              :: igrid, iigrid

    ! initial cleaning is needed
    if ( initialize_mag_field ) then
       call phys_initial_clean_divb(it,global_time)
    end if

    ! update aux variables
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       block => ps(igrid) 
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       call update_aux_vars(igrid,ps(igrid)%level,ixG^LL,ixM^LL,global_time, ps(igrid)) 
    end do
  end subroutine rm_improve_initial_condition

  subroutine update_aux_vars(igrid,level,ixI^L,ixO^L,qt,s)
    use mod_global_parameters
    use mod_gremhd_phys_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s
    double precision                :: divb(ixI^S), dive(ixI^S)
    !> calculate div B after the advance 
    call gremhd_get_divb(ixI^L, ixO^L, s, divb(ixI^S))
    s%prim(ixO^S, divB_) = dabs( divb(ixO^S) )
    ! also the div E
    call gremhd_get_div(ixI^L, ixO^L, s%prim(ixI^S,Evec(1:ndir)), dive(ixI^S), dive_4thorder )
    s%prim(ixO^S, divE_) = ( dive(ixO^S) )
  end subroutine update_aux_vars

  subroutine printlog()
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_gremhd_phys_parameters
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

    integer              :: igrid, iigrid, idir
    double precision     :: lfac(ixG^T)

    double precision     :: Bvector(ixG^T,1:ndir)
    double precision     :: Einit(ixG^T,1:ndir)
    double precision     :: divb(ixG^T)
    double precision     :: r2(ixG^T)
    double precision     :: q_diff(ixG^T), q_exact(ixG^T), q_init(ixG^T)
    double precision     :: Bz_diff(ixG^T)
    double precision     :: Bz_exact(ixG^T)

    integer              :: total_volume = -1
    integer              :: total_q = -1
    integer              :: L1_norm_q = -1
    integer              :: L1_norm_Bz = -1
    integer              :: total_divBdV = -1
    integer              :: total_divB2dV = -1
    integer              :: n_var = -1

    ! initialize variables
    n_var = 0
    total_volume = add_var()
    total_q = add_var()
    L1_norm_q = add_var()
    L1_norm_Bz = add_var()
    total_divBdV = add_var()
    total_divB2dV = add_var()

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! calculate L1_q
       r2(ixG^T) = ( ps(igrid)%x(ixG^T,1)**2 + ps(igrid)%x(ixG^T,2)**2 )
       Einit(ixG^T, 1) = 0.5d0 * q0 * ps(igrid)%x(ixG^T,1) / (r2(ixG^T) + 1.0d0)
       Einit(ixG^T, 2) = 0.5d0 * q0 * ps(igrid)%x(ixG^T,2) / (r2(ixG^T) + 1.0d0)
       Einit(ixG^T, 3) = 0.0d0
       call gremhd_get_div(ixG^LL, ixM^LL, Einit(ixG^T,1:ndir), q_init(ixG^T), dive_4thorder )
       !q_exact(ixM^T) = q0 / ( r2(ixM^T) + 1.0d0 )**2
       !q_diff(ixM^T) = dabs( 1.0d0 - ps(igrid)%prim(ixM^T, divE_) / q_exact(ixM^T)  )
       !q_diff(ixM^T) = dabs( 1.0d0 - ps(igrid)%prim(ixM^T, divE_) / q_init(ixM^T)  )
       q_diff(ixM^T) = dabs( ps(igrid)%prim(ixM^T, divE_) - q_init(ixM^T)  )
       send_buffer(L1_norm_q) = send_buffer(L1_norm_q) + sum(q_diff(ixM^T) * ps(igrid)%dvolume(ixM^T) )

       ! calculate L1_Bz
       Bz_exact(ixM^T) = dsqrt( (r2(ixM^T) + 1.0d0)**2 - 0.25d0 * q0**2 ) / (r2(ixM^T) + 1.0d0)
       Bz_diff(ixM^T) = dabs( ps(igrid)%prim(ixM^T, Bvec(3)) - Bz_exact(ixM^T) )
       send_buffer(L1_norm_Bz) = send_buffer(L1_norm_Bz) + sum(Bz_diff(ixM^T) * ps(igrid)%dvolume(ixM^T) )

       ! calculate total volume
       send_buffer(total_volume) = send_buffer(total_volume) + sum(ps(igrid)%dvolume(ixM^T))
       ! calculate total charge
       send_buffer(total_q) = send_buffer(total_q) &
                  + sum( ps(igrid)%prim(ixM^T, divE_) * ps(igrid)%dvolume(ixM^T) )

       ! calculate total divB dV
       send_buffer(total_divBdV) = send_buffer(total_divBdV) + &
          sum( dabs(ps(igrid)%prim(ixM^T, divB_)) * ps(igrid)%dvolume(ixM^T) )

       ! calculate total divB^2 dV
       send_buffer(total_divB2dV) = send_buffer(total_divB2dV) + &
          sum( dabs(ps(igrid)%prim(ixM^T, divB_))**2 * ps(igrid)%dvolume(ixM^T) )
    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    ! Volume average of divB
    recive_buffer(total_divBdV) = (recive_buffer(total_divBdV) / recive_buffer(total_volume))
    ! L2 norm of divB
    recive_buffer(total_divB2dV) = dsqrt( recive_buffer(total_divB2dV) / recive_buffer(total_volume) )

    !recive_buffer(total_q) = recive_buffer(total_q) ! / recive_buffer(total_volume)

    recive_buffer(L1_norm_q) = recive_buffer(L1_norm_q) &
           / recive_buffer(total_volume) !/ recive_buffer(total_q)
    recive_buffer(L1_norm_Bz) = recive_buffer(L1_norm_Bz) / recive_buffer(total_volume)

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
          line = "it global_time dt q_tot L1_norm_q L1_norm_Bz avg_divB L2_divB"

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

       write(fmt_string, '(a,i0,a)') '(', n_var-1, fmt_r // ')'
       write(line(i:), fmt_string) recive_buffer(2:n_var)
       i = len_trim(line) + 2

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

    deallocate(send_buffer)
    deallocate(recive_buffer)
    
    contains

       function add_var() result (i_var)
          integer :: i_var
          n_var = n_var + 1
          i_var = n_var
       end function add_var

  end subroutine printlog

  ! get resistivity
  subroutine my_get_resistivity(ixI^L,ixO^L,cons,x,eta)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(in) :: cons(ixI^S,1:ncons) 
    double precision, intent(out):: eta(ixI^S) 
    eta(ixO^S) = resistivity
  end subroutine my_get_resistivity

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)
    use mod_eos
    ! special boundary types, user defined
    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nprim)

    double precision                :: v(ixG^S,1:ndir)
    double precision                :: lfac(ixG^S)
    double precision                :: r2(ixG^S) !, cos_theta(ixG^S), sin_theta(ixG^S)
    double precision                :: Er(ixG^S), vphi(ixG^S)
    integer :: ixO^L, ix^D, idir
   
    ! no gauge effect
    w(ixB^S,psi_) = 1.0d0
    w(ixB^S,alp_) = 1.0d0
    w(ixB^S,beta(1)) = 0.0d0
    w(ixB^S,beta(2)) = 0.0d0
    w(ixB^S,beta(3)) = 0.0d0

    r2(ixB^S) = ( x(ixB^S,1)**2 + x(ixB^S,2)**2 )

    !Er(ixB^S) = 0.5d0 * q0 * r(ixB^S) / (r(ixB^S)**2 + 1.0d0)
    !vphi(ixB^S) = - 0.5d0 * q0 * r(ixB^S) / dsqrt( (r(ixB^S)**2 + 1.0d0)**2 - 0.25d0 * q0**2 )

    w(ixB^S,rho_) = 1.0d0
    w(ixB^S,press_) = - w(ixB^S,rho_) * ( eos_gamma - 1.0d0 ) / eos_gamma &
                      + ( p0 + w(ixB^S,rho_) * ( eos_gamma - 1.0d0 ) / eos_gamma ) * &
                        ( (4.0d0*r2(ixB^S) + 4.0d0 - q0**2) / (r2(ixB^S) + 1.0d0) / (4.0d0-q0**2) )**(0.5d0*eos_gamma/(eos_gamma-1.0d0))

    w(ixB^S, Bvec(3))    = dsqrt( (r2(ixB^S) + 1.0d0)**2 - 0.25d0 * q0**2 ) / (r2(ixB^S) + 1.0d0)

    w(ixB^S, Evec(1)) = 0.5d0 * q0 * x(ixB^S,1) / (r2(ixB^S) + 1.0d0)
    w(ixB^S, Evec(2)) = 0.5d0 * q0 * x(ixB^S,2) / (r2(ixB^S) + 1.0d0)
    w(ixB^S, Evec(3)) = 0.0d0

    v(ixB^S, 1) =   0.5d0 * q0 * x(ixB^S,2) / dsqrt( (r2(ixB^S) + 1.0d0)**2 - 0.25d0 * q0**2 )
    v(ixB^S, 2) = - 0.5d0 * q0 * x(ixB^S,1) / dsqrt( (r2(ixB^S) + 1.0d0)**2 - 0.25d0 * q0**2 )
    v(ixB^S, 3) = 0.0d0

    ! get W
    lfac(ixB^S) = 1.0d0
    ! calculate v^2 first
    do idir = 1, ndir
       lfac(ixB^S) = lfac(ixB^S) - v(ixB^S, idir)**2
    end do
    ! Calculate the Lorentz factor from velocity
    lfac(ixB^S) = dsqrt( 1.0d0 / lfac(ixB^S) )

    ! set veloc -> W_vel
    w(ixB^S,W_vel(1))   = v(ixB^S, 1) * lfac(ixB^S)
    w(ixB^S,W_vel(2))   = v(ixB^S, 2) * lfac(ixB^S)
    w(ixB^S,W_vel(3))   = v(ixB^S, 3) * lfac(ixB^S)

    w(ixB^S, divB_) = 0.0d0
    w(ixB^S, divE_) = q0 / (r2(ixB^S) + 1.0d0)**2

    ! dont forget to update these primitive variables
    {do ix^D = ixB^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
    ! update rest of the prim
    call phys_update_eos(ixG^L, ixB^L, w)

  end subroutine specialbound_usr

end module mod_usr
