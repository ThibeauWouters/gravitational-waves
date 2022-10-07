module mod_usr
  use mod_physics
  use mod_grmhd
  use mod_eos_idealgas

  implicit none
  private

  ! aux output variables
  integer, protected  :: divB_ = -1
  integer, protected  :: B2_ = -1
  integer, protected  :: mag_press_ = -1
  integer, protected  :: gamma_ = -1

  logical, parameter                     :: initialize_mag_field = .false.
  double precision, parameter            :: r_loop = 0.3d0
  double precision, parameter            :: A_loop = 1.0d-3
  double precision, parameter            :: X_0 = -1.0d0
  double precision, parameter            :: Y_0 = -0.5d0
  

  public :: usr_init

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables
    use mod_eos

    usr_init_one_grid    => rm_init_one_grid
    usr_improve_initial_condition => rm_improve_initial_condition
    usr_refine_grid => my_refine
    usr_process_adv_grid => update_aux_vars
    usr_print_log => printlog

    call set_coordinate_system("polar_2.5D")

    call grmhd_activate()
    call eos_idealgas_activate()

    divB_ = var_set_auxvar('divB')
    B2_ = var_set_auxvar('B2')
    mag_press_ = var_set_auxvar('mag_press')
    gamma_ = var_set_auxvar('gamma')

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm_init_one_grid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics
    use mod_grmhd
    use mod_eos

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nprim) ! primitive variables
    double precision                :: v(ixI^S,1:ndir)
    double precision                :: lfac(ixI^S)
    double precision                :: equ
    double precision                :: Bx, By
    double precision                :: xi(ixI^S,1:ndim)

    integer  ::  ix^D, idir, ixC^L

    w = 0.0d0

    ! no gauge effect
    w(ixO^S,psi_) = 1.0d0
    w(ixO^S,alp_) = 1.0d0
    w(ixO^S,beta(1)) = 0.0d0
    w(ixO^S,beta(2)) = 0.0d0
    w(ixO^S,beta(3)) = 0.0d0

    w(ixO^S,rho_)       = 1.0d0
    w(ixO^S,press_)     = 3.0d0

    ! in cart coord.
    v(ixO^S, 1)   =  0.2d0
    v(ixO^S, 2)   =  1.0d-1 ! / 2.4d1
    lfac(ixO^S) = v(ixO^S,1)**2 + v(ixO^S,2)**2
    lfac(ixO^S) = 1.0d0 / dsqrt( 1.0d0 - lfac(ixO^S) )

    ! transform in cylindrical
    w(ixO^S,W_vel(1))   =   v(ixO^S, 1) * dcos( x(ixO^S,2) ) + v(ixO^S, 2) * dsin( x(ixO^S,2) )
    w(ixO^S,W_vel(2))   = - v(ixO^S, 1) * dsin( x(ixO^S,2) ) + v(ixO^S, 2) * dcos( x(ixO^S,2) )

    w(ixO^S,W_vel(1))   = w(ixO^S,W_vel(1)) * lfac(ixO^S)
    w(ixO^S,W_vel(2))   = w(ixO^S,W_vel(2)) * lfac(ixO^S) / x(ixO^S,1)
    w(ixO^S,W_vel(3))   = 0.0d0

    {do ix^D = ixO^LIM^D \}
       equ = x(ix^D,1)**2 - 2 * x(ix^D,1) * ( X_0 * dcos(x(ix^D,2)) + Y_0 * dsin(x(ix^D,2)) )&
             - r_loop**2 + X_0**2 + Y_0**2
       if ( equ <= 0.0d0 ) then
          equ = x(ix^D,1)**2 - 2 * x(ix^D,1) * ( X_0 * dcos(x(ix^D,2)) + Y_0 * dsin(x(ix^D,2)) ) + X_0**2 + Y_0**2
          equ = dsqrt( dabs(equ) )
          Bx = - A_loop / equ * ( x(ix^D,1) * dsin(x(ix^D,2)) - Y_0 )
          By =   A_loop / equ * ( x(ix^D,1) * dcos(x(ix^D,2)) - X_0 )
          w(ix^D,Bvec(1))    =  Bx * dcos(x(ix^D,2)) + By * dsin(x(ix^D,2))
          w(ix^D,Bvec(2))    =  By * dcos(x(ix^D,2)) - Bx * dsin(x(ix^D,2))
          w(ix^D,Bvec(2))    =  w(ix^D,Bvec(2)) / x(ix^D,1)
          !w(ix^D,Bvec(1))    =   A_loop / equ * ( - X_0 * dsin(x(ix^D,2)) + Y_0 * dcos(x(ix^D,2)) )
          !w(ix^D,Bvec(2))    =   A_loop / equ * ( x(ix^D,1) - X_0 * dcos(x(ix^D,2)) - Y_0 * dsin(x(ix^D,2)) ) / x(ix^D,1)
       end if
    {enddo^D&\}

    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_pressure_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_cs2_one_grid(w(ix^D, cs2_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
  end subroutine rm_init_one_grid
  
  !> before the main loop, we improve the initial data here
  subroutine rm_improve_initial_condition()
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_physics
    use mod_eos

    integer              :: ix^D
    integer              :: igrid, iigrid
    double precision                :: v(ixG^T,1:ndir)
    double precision                :: lfac(ixG^T)
    double precision                :: divb(ixG^T)

    ! initial cleaning is needed
    if ( initialize_mag_field ) then
       call phys_initial_clean_divb(it,global_time)
    end if
    ! update divB
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       call update_aux_vars(igrid,1,ixG^LL,ixM^LL,global_time, ps(igrid)) ! level is not used here, set as 1
    end do
  end subroutine rm_improve_initial_condition

  !> calculate div B after the advance 
  subroutine update_aux_vars(igrid,level,ixI^L,ixO^L,qt,s)
    use mod_global_parameters
    use mod_grmhd_phys_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s
    double precision                :: lfac2(ixI^S)
    double precision                :: tmp(ixI^S)

    call phys_get_lfac2(ixI^L, ixO^L, s, lfac2(ixI^S))
    s%prim(ixO^S, gamma_) = dsqrt( lfac2(ixO^S) )

    call grmhd_get_divb(ixI^L, ixO^L, s, tmp(ixI^S))
    s%prim(ixO^S, divB_) = dabs( tmp(ixO^S) )

    tmp(ixO^S) =    s%prim(ixO^S, Bvec(1))**2 &
                  + s%x(ixO^S, 1)**2 * s%prim(ixO^S, Bvec(2))**2 &
                  + s%prim(ixO^S, Bvec(3))**2 
    s%prim(ixO^S, B2_) = dabs( tmp(ixO^S) ) ! B2

    tmp(ixO^S) =    s%prim(ixO^S, Bvec(1))*s%prim(ixO^S, W_vel(1)) &
                  + s%x(ixO^S, 1)**2 * s%prim(ixO^S, Bvec(2))*s%prim(ixO^S, W_vel(2)) &
                  + s%prim(ixO^S, Bvec(3))*s%prim(ixO^S, W_vel(3))
    s%prim(ixO^S, mag_press_) = 0.5d0 * ( s%prim(ixO^S, B2_) + tmp(ixO^S)**2 ) / lfac2(ixO^S)
 
  end subroutine update_aux_vars

  subroutine my_refine(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nprim), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    integer                      :: i_level, to_level
    double precision             :: ratio(ixI^S)
    double precision             :: B2(ixI^S)
    double precision, parameter  :: B2_min = 1.0d-10

    refine = 0
    coarsen = 0

    B2(ixO^S) = w(ixO^S, Bvec(1))**2 + w(ixO^S, Bvec(1))**2

    if ( any( B2(ixO^S) > B2_min ) ) then
       refine = 1
       coarsen = -1
    else
       refine = -1
       coarsen = 1
    end if

    if ( minval(x(ixO^S,1)) <= 1.0d-1 ) then
       ratio(ixO^S) = dx(2, level) * x(ixO^S,1) / dx(1, level)
       ratio(ixO^S) = ( ratio(ixO^S) - 1.0d0 )
       if ( any( ratio(ixO^S) <= 0.0d0 ) ) then
          refine = -1
          coarsen = 1
       else
          !refine = 1
          !coarsen = -1
       end if
    end if
  end subroutine my_refine

  subroutine printlog()
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_grmhd_phys_parameters
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
    double precision     :: divb(ixG^T)

    integer, parameter   :: total_mass = 1
    integer, parameter   :: total_divBdV = 2
    integer, parameter   :: total_divB2dV = 3
    integer, parameter   :: total_volume = 4
    integer, parameter   :: n_var = 4

    allocate(send_buffer(1:n_var))
    allocate(recive_buffer(1:n_var))
    send_buffer = 0.0d0
    recive_buffer = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! calculate total rest mass
       call phys_get_lfac2(ixG^LL, ixM^LL, ps(igrid), lfac(ixG^T))
       lfac(ixM^T) = dsqrt( lfac(ixM^T) )
       send_buffer(total_mass) = send_buffer(total_mass) + &
          sum(lfac(ixM^T) * ps(igrid)%prim(ixM^T,psi_)**6*ps(igrid)%prim(ixM^T,rho_)*ps(igrid)%dvolume(ixM^T))

       ! calculate total volume
       send_buffer(total_volume) = send_buffer(total_volume) + sum(ps(igrid)%dvolume(ixM^T))

       ! get divB
       call grmhd_get_divb(ixG^LL, ixM^LL, ps(igrid), divb(ixG^T))

       ! calculate total divB dV
       send_buffer(total_divBdV) = send_buffer(total_divBdV) + &
          sum(dabs(divB(ixM^T)) * ps(igrid)%dvolume(ixM^T))

       ! calculate total divB^2 dV
       send_buffer(total_divB2dV) = send_buffer(total_divB2dV) + &
          sum(divB(ixM^T)**2 * ps(igrid)%dvolume(ixM^T))
    end do

    ! gather all the results!
    call MPI_ALLREDUCE(send_buffer, recive_buffer, n_var, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    ! Volume average of divB
    recive_buffer(total_divBdV) = (recive_buffer(total_divBdV) / recive_buffer(total_volume))
    ! L2 norm of divB
    recive_buffer(total_divB2dV) = dsqrt( recive_buffer(total_divB2dV) / recive_buffer(total_volume) )

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
          line = "it global_time dt M0 avg_divB L2_divB"

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
       write(line(i:), fmt_string) recive_buffer(1:n_var-1)
       i = len_trim(line) + 2

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

    deallocate(send_buffer)
    deallocate(recive_buffer)

  end subroutine printlog

end module mod_usr

!@ARTICLE{2011ApJS..193....6B,
!       author = {{Beckwith}, Kris and {Stone}, James M.},
!        title = "{A Second-order Godunov Method for Multi-dimensional
!Relativistic Magnetohydrodynamics}",
!      journal = {\apjs},
!     keywords = {magnetohydrodynamics: MHD, methods: numerical, relativistic
!processes, Astrophysics - High Energy Astrophysical Phenomena, General
!Relativity and Quantum Cosmology, Physics - Computational Physics},
!         year = 2011,
!        month = mar,
!       volume = {193},
!       number = {1},
!          eid = {6},
!        pages = {6},
!          doi = {10.1088/0067-0049/193/1/6},
!archivePrefix = {arXiv},
!       eprint = {1101.3573},
! primaryClass = {astro-ph.HE},
!       adsurl = {https://ui.adsabs.harvard.edu/abs/2011ApJS..193....6B},
!      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
!}


