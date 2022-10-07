module mod_variables
  use mod_basic_types

  implicit none
  public

  integer, parameter:: hydro_var = 1
  integer, parameter:: metric_var = 3

  !> Total number of flux (conservative) variables
  integer           :: ncons = 0

  !> Number of conserved hydro variables
  integer           :: nc_hydro = 0
  !> index of min conserved hydro variables
  integer           :: nc_hydro_lo = 0
  !> index of max conserved hydro variables
  integer           :: nc_hydro_hi = 0

  !> Number of primitive variables which need user to specify boundary type
  integer           :: nprimbc = 0

  !> Number of primitive variables
  integer           :: nprim = 0

  !> Number of primitive hydro variables
  integer           :: nhydro = 0
  !> index of min primitive hydro variables
  integer           :: nhydro_lo = 0
  !> index of max primitive hydro variables
  integer           :: nhydro_hi = 0

  !> Number of metric variables
  integer           :: nmetric = 0
  !> index of min metric variables
  integer           :: nmetric_lo = 0
  !> index of max metric variables
  integer           :: nmetric_hi = 0

  !> Number of stagger conservative variables
  integer           :: nconss = 0

  !> Number of vector variables (used for writing output)
  integer           :: nvector = 0

  !> Indices of vector variables
  integer, dimension(:), allocatable :: iw_vector


  !> Number of auxiliary variables
  integer :: nwaux

  !> Number of auxiliary variables that are only included in output
  integer :: nwauxio

  !> Number of extra variables
  integer :: nwextra

  !> Maximum number of variables
  integer, parameter :: max_nvar = 50

  logical            :: var_reconstruct(max_nvar)

  !> Primitive variable names
  character(len=name_len) :: prim_names(max_nvar)

  !> Conservative variable names
  character(len=name_len) :: cons_names(max_nvar)

contains

  !> Set primitive variable
  function var_set_primvar(name_prim, ix, var_type, need_bc, need_rec) result(iw)
    character(len=*), intent(in)  :: name_prim !< Primitive name
    integer, intent(in), optional :: ix        !< Optional index (to make var1, var2, ...)
    integer, intent(in), optional :: var_type  !< variable type
    logical, intent(in), optional :: need_bc   !< Require boundary condition (default: true)
    logical, intent(in), optional :: need_rec  !< Require reconstruction for hydro only (default: true)

    integer                       :: iw, var_t
    logical                       :: flag

    ! total number of primitive variables
    nprim  = nprim + 1
    iw     = nprim

    flag = .true.
    if (present(need_bc)) flag = need_bc
    if (flag) nprimbc = nprimbc + 1

    if (.not. present(ix)) then
      prim_names(nprim) = name_prim
    else
      write(prim_names(nprim),"(A,I0)") name_prim, ix
    end if

    var_t = -1
    if (present(var_type)) var_t = var_type

    select case (var_t)
    case (hydro_var)
       if ( nhydro_lo == 0 ) nhydro_lo = iw
       nhydro = nhydro + 1
       nhydro_hi = iw

       flag = .true.
       if (present(need_rec)) flag = need_rec
       var_reconstruct(iw) = flag
    case (metric_var)
       if ( nmetric_lo == 0 ) nmetric_lo = iw
       nmetric = nmetric + 1
       nmetric_hi = iw
    case default
       ! nothing to do here, this is an aux variable
    end select
  end function var_set_primvar

  !> Set conservative variable
  function var_set_consvar(name_cons, ix, var_type) result(iw)
    character(len=*), intent(in)  :: name_cons !< Conservative name
    integer, intent(in), optional :: ix        !< Optional index (to make var1, var2, ...)
    integer, intent(in), optional :: var_type  !< variable type
    integer                       :: iw, var_t

    ncons  = ncons + 1
    iw     = ncons

    if (.not. present(ix)) then
      cons_names(ncons) = name_cons
    else
      write(cons_names(ncons),"(A,I0)") name_cons, ix
    end if

    var_t = -1
    if (present(var_type)) var_t = var_type

    select case (var_t)
    case (hydro_var)
       if ( nc_hydro_lo == 0 ) nc_hydro_lo = iw
       nc_hydro = nc_hydro + 1
       nc_hydro_hi = iw
    case default
       ! nothing to do here, this is an aux variable
    end select
  end function var_set_consvar

  !> add aux variable
  function var_set_auxvar(name_prim, ix) result(iw)
    character(len=*), intent(in)  :: name_prim !< Primitive name
    integer, intent(in), optional :: ix        !< Optional index (to make var1, var2, ...)
    integer                       :: iw
    character(len=2048)           :: error_message
    
    if ( nprim <= 0 ) then
       error_message = "cannot add aux var " // trim(name_prim)
       error_message = trim(error_message) // "add aux var after physics initialization!"
       call mpistop(error_message)
    end if

    ! total number of primitive variables
    nprim  = nprim + 1
    iw     = nprim
    if (.not. present(ix)) then
      prim_names(nprim) = name_prim
    else
      write(prim_names(nprim),"(A,I0)") name_prim, ix
    end if
  end function var_set_auxvar


end module mod_variables
