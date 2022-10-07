module mod_map_profile
  implicit none
  private

  double precision, parameter :: ggrav = 6.673d-8
  double precision, parameter :: clite = 2.99792458d10
  double precision, parameter :: clite_g = 2.99792458d10

  double precision, parameter :: rho_gf = 1.61930347d-18
  double precision, parameter :: press_gf = 1.80171810d-39
  double precision, parameter :: eps_gf = 1.11265006d-21
  double precision, parameter :: time_gf = 2.03001708d+05
  double precision, parameter :: mass_gf = 5.02765209d-34
  double precision, parameter :: length_gf = 6.77140812d-06
  double precision, parameter :: energy_gf = 5.59424238d-55
  double precision, parameter :: lum_gf = 2.7556091d-60

  double precision, parameter :: mev_to_erg = 1.60217733d-6
  double precision, parameter :: erg_to_mev = 6.24150636d5
  double precision, parameter :: amu_cgs = 1.66053873d-24
  double precision, parameter :: massn_cgs = 1.674927211d-24
  double precision, parameter :: amu_mev = 931.49432d0
  double precision, parameter :: kb_erg = 1.380658d-16
  double precision, parameter :: kb_mev = 8.61738568d-11
  double precision, parameter :: temp_mev_to_kelvin = 1.1604447522806d10
  double precision, parameter :: planck = 6.626176d-27
  double precision, parameter :: avo = 6.0221367d23
  double precision, parameter :: hbarc_mevcm = 1.97326966d-11
  double precision, parameter  :: kboltz_cgs = 1.380662d-16

  type profile_t
     integer                       :: Nr
     double precision, allocatable :: radius(:)
     double precision, allocatable :: mass(:)
     double precision, allocatable :: rho(:)
     double precision, allocatable :: press(:)
     double precision, allocatable :: temp(:)
     double precision, allocatable :: eps(:)
     double precision, allocatable :: vel(:)
     double precision, allocatable :: ye(:)
     double precision, allocatable :: omega(:)
  end type profile_t

  ! save the profle
  type(profile_t) :: prof

  public :: profile_t
  public :: prof
  public :: mod_read_profile
  public :: mod_map

contains

subroutine mod_read_profile(profile_path)
  character*(*) :: profile_path
  character*128 :: filename
  double precision :: minrho, buffer

  integer :: Nr
  integer :: i, ibuffer

  open(666,file=trim(profile_path),status='unknown', & 
       form='formatted',action='read')
  read(666,*) Nr

  prof%Nr = Nr

  allocate(prof%radius(Nr))
  allocate(prof%mass(Nr))
  allocate(prof%rho(Nr))
  allocate(prof%press(Nr))
  allocate(prof%temp(Nr))
  allocate(prof%eps(Nr))
  allocate(prof%vel(Nr))
  allocate(prof%ye(Nr))
  allocate(prof%omega(Nr))
  
  ! read profile
  do i = 1, Nr
     read(666,*) ibuffer, &
            prof%mass(i), &
            prof%radius(i), &
            prof%temp(i), &
            prof%rho(i), &
            prof%vel(i), &
            prof%ye(i), &
            buffer
            !prof%press(i), &
            !prof%eps(i), &
            !prof%omega(i), &
  end do

  close(666)

  ! go to c=G=Msun=1
  prof%radius = prof%radius * length_gf
  prof%rho = prof%rho * rho_gf
  prof%vel = prof%vel / clite
  prof%omega = prof%omega / time_gf
  prof%temp = prof%temp * kboltz_cgs / mev_to_erg

  !minrho = minval(prof%rho)
  !do i=1,Nr
     !if(prof%rho(i,j)<=minrho) prof%rho(i,j) = 0.0d0
  !enddo
end subroutine mod_read_profile

! **************************************************************
subroutine mod_linterp(x1,x2,y1,y2,x,y)

! perform linear interpolation      
  implicit none

  real*8 slope,x1,x2,y1,y2,x,y

  if (x2.lt.x1) then
     stop "Error in linterp!"
  endif

  slope = (y2 - y1) / (x2 - x1)

  y = slope*(x-x1) + y1
 
end subroutine  mod_linterp

! ***************************************************************

subroutine mod_find_index(zones,array,goal,upper_index,lower_index)
  
! bisection search
  implicit none
  
  integer zones,i
  real*8 array(zones)
  real*8 goal
  integer middle_index,upper_index,lower_index

  lower_index = 1
  upper_index = zones
  
  do while ( (upper_index - lower_index) .gt. 1 )
     middle_index = (lower_index + upper_index) * 0.5
     if ( (goal .ge. array(lower_index)) &
          .and. (goal .le. array(middle_index)) ) then
        upper_index = middle_index
     else
        if ( (goal .ge. array(middle_index)) &
             .and. (goal .le. array(upper_index)) ) then
           lower_index = middle_index
        endif
     endif
  enddo
      
end subroutine mod_find_index

! ******************************************************************

subroutine mod_map(point_value,point_radius0,parray,pradius,zones)

  implicit none
  
  real*8 point_value, point_radius, point_radius0
  real*8 pradius(*), parray(*)
  integer zones
  integer upper_index, lower_index

  point_radius = abs(point_radius0)
  
  if (point_radius .ge. pradius(1) .and. & 
       point_radius .lt. pradius(zones) )  then
     
     call mod_find_index(zones,pradius,point_radius, &
          upper_index,lower_index)
     
     call mod_linterp( pradius(lower_index),pradius(upper_index), &
          parray(lower_index), parray(upper_index),  & 
          point_radius, point_value )

  else if (point_radius .lt. pradius(1)) then
     ! linear extrapolation
     call mod_linterp(pradius(1),pradius(2), & 
          parray(1),parray(2),point_radius,point_value)

  else if (point_radius .gt. pradius(zones)) then
     ! linear extrapolation
     call mod_linterp(pradius(zones-1),pradius(zones), & 
          parray(zones-1),parray(zones),point_radius,point_value)
  endif
  
  
end subroutine mod_map

end module mod_map_profile

