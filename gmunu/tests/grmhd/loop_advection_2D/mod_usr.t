module mod_usr
  use mod_physics
  use mod_grmhd
  use mod_eos_idealgas

  implicit none
  private

  ! aux output variables
  integer, protected  :: divB_ = -1

  logical, parameter                     :: init_b_from_vector_pot = .true.
  logical, parameter                     :: initialize_mag_field = .true.
  double precision, parameter            :: r_loop = 0.3d0
  double precision, parameter            :: A_loop = 1.0d-3
  

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
    usr_process_adv_grid => cal_div_B
    usr_print_log => printlog
    if (init_b_from_vector_pot) usr_init_vector_potential => init_vec_pot_usr

    call set_coordinate_system("Cartesian_2.5D")

    call grmhd_activate()
    call eos_idealgas_activate()

    divB_ = var_set_auxvar('divB')

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
    double precision                :: r(ixI^S)
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

    v(ixO^S, 1)   =  0.2d0
    v(ixO^S, 2)   =  1.0d-1 ! / 2.4d1

    lfac(ixO^S) = v(ixO^S,1)**2 + v(ixO^S,2)**2
    lfac(ixO^S) = 1.0d0 / dsqrt( 1.0d0 - lfac(ixO^S) )

    w(ixO^S,W_vel(1))   = v(ixO^S, 1) * lfac(ixO^S)
    w(ixO^S,W_vel(2))   = v(ixO^S, 2) * lfac(ixO^S)
    w(ixO^S,W_vel(3))   = 0.0d0

    if (.not. init_b_from_vector_pot) then
       r(ixO^S) = dsqrt( x(ixO^S,1)**2 + x(ixO^S,2)**2 )
       where ( r(ixO^S) <= r_loop )
          w(ixO^S,Bvec(1))    = -A_loop * x(ixO^S, 2) / r(ixO^S)
          w(ixO^S,Bvec(2))    =  A_loop * x(ixO^S, 1) / r(ixO^S)
          w(ixO^S,Bvec(3))    = 0.0d0
       else where
          w(ixO^S,Bvec(1))    = 0.0d0
          w(ixO^S,Bvec(2))    = 0.0d0
          w(ixO^S,Bvec(3))    = 0.0d0
       end where
       if (stagger_grid) then
          block%conss = 0.0d0
          do idir=1,ndim
            ixCmin^D=ixImin^D;
            ixCmax^D=ixImax^D-kr(idir,^D);
            xi = x
            xi(ixI^S,idir)=xi(ixI^S,idir)+0.5d0*block%dx(ixI^S,idir)
            ! cell-face B_idir
            r(ixC^S) = dsqrt( xi(ixC^S,1)**2 + xi(ixC^S,2)**2 )
            where ( r(ixC^S) <= r_loop )
               block%conss(ixC^S,idir) = A_loop * xi(ixC^S, ndim-idir+1) / r(ixC^S)
            else where
               block%conss(ixC^S,idir) = 0.0d0
            end where
          end do
          block%conss(ixI^S,1) = - block%conss(ixI^S,1)
          call phys_face_to_center(ixO^L, block)
       end if
    else
       call grmhd_b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block)
    end if

    !if (stagger_grid) &
    !   call phys_face_to_center(ixO^L, block)

    w(ixO^S, divB_) = 0.0d0

    {do ix^D = ixO^LIM^D \}
       call eos_get_eps_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_pressure_one_grid(w(ix^D, press_),w(ix^D, rho_),w(ix^D, eps_))
       call eos_get_cs2_one_grid(w(ix^D, cs2_),w(ix^D, rho_),w(ix^D, eps_))
    {enddo^D&\}
  end subroutine rm_init_one_grid

  subroutine init_vec_pot_usr(ixI^L, ixO^L, x, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixO^L,idir
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)
    double precision                   :: r(ixI^S)
    A = 0.d0
    if (idir==3) then
      r(ixO^S) = dsqrt( x(ixO^S,1)**2 + x(ixO^S,2)**2 )
      where( r(ixO^S) <= r_loop )
        A(ixO^S) = A_loop * ( r_loop - r(ixO^S) )
      end where
    end if
  end subroutine init_vec_pot_usr
  
  !> before the main loop, we improve the initial data here
  subroutine rm_improve_initial_condition()
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters
    use mod_physics
    use mod_eos

    integer              :: ix^D
    integer              :: igrid, iigrid
    double precision                :: lfac(ixG^T)
    double precision                :: divb(ixG^T)

    ! initial cleaning is needed
    if ( initialize_mag_field ) then
       call phys_initial_clean_divb(it,global_time)
    end if
    ! update divB
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       call grmhd_get_divb(ixG^LL, ixM^LL, ps(igrid), divb(ixG^T))
       ps(igrid)%prim(ixM^T, divB_) = ( divb(ixM^T) )
    end do
  end subroutine rm_improve_initial_condition

  !> calculate div B after the advance 
  subroutine cal_div_B(igrid,level,ixI^L,ixO^L,qt,s)
    use mod_global_parameters
    use mod_grmhd_phys_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s
    double precision                :: divb(ixI^S)
    call grmhd_get_divb(ixI^L, ixO^L, s, divb(ixI^S))
    s%prim(ixO^S, divB_) = dabs( divb(ixO^S) )
  end subroutine cal_div_B

  subroutine my_refine(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nprim), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    integer                      :: i_level, to_level
    double precision                :: B2(ixI^S)
    double precision, parameter     :: B2_min = 1.0d-10

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


