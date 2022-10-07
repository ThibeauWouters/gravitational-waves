module mod_XNS
  implicit none
  private

  type XNS_profile
     integer                       :: Nr, Nth
     double precision, allocatable :: radius(:), theta(:)
     double precision, allocatable :: rho(:,:), press(:,:), v3(:,:)
     double precision, allocatable :: psi(:,:), alp(:,:), beta3(:,:)
  end type XNS_profile

  ! save the profle
  type(XNS_profile) :: prof

  public :: XNS_profile
  public :: prof
  public :: mod_XNS_read_profile
  !public :: mod_XNS_get_dim
  !public :: mod_XNS_get_grid
  !public :: mod_XNS_get_profile
  public :: mod_XNS_map

contains

subroutine mod_XNS_read_profile(profile_path)

  character*(*) profile_path
  character*128 :: filename
  double precision :: minrho

  integer :: Nth, Nr
  integer :: i,j

  call mod_XNS_get_dim(profile_path, Nr, Nth)

  prof%Nr = Nr
  prof%Nth = Nth
  allocate(prof%radius(Nr))
  allocate(prof%theta(Nth))
  allocate(prof%rho(Nr,Nth),prof%press(Nr,Nth),prof%v3(Nr,Nth))
  allocate(prof%psi(Nr,Nth),prof%alp(Nr,Nth),prof%beta3(Nr,Nth))
  
  ! read profile
  call mod_XNS_get_grid(profile_path, Nr, Nth, prof%radius, prof%theta)
  call mod_XNS_get_profile(profile_path, Nr, Nth, &
                             prof%rho,prof%press,prof%v3,prof%psi,prof%alp,prof%beta3)


  minrho = minval(prof%rho)

  do j=1,Nth
     do i=1,Nr
        if(prof%rho(i,j)<=minrho) prof%rho(i,j) = 0.0d0
     enddo
  enddo
end subroutine mod_XNS_read_profile

subroutine mod_XNS_get_dim(lprofile_name, Nr, Nth)

  character*(*) lprofile_name
  character*128 :: filename

  integer, intent(out) :: Nth, Nr
  
  !--------------------------------------------------------!       
  ! read Grid.dat profile      
  !--------------------------------------------------------!       
  filename = trim(adjustl(lprofile_name))//"/Grid.dat"
  open(666,file=trim(filename),status='unknown', & 
       form='formatted',action='read')
  read(666,*) Nth, Nr

  close(666)

end subroutine mod_XNS_get_dim


subroutine mod_XNS_get_grid(lprofile_name, Nr, Nth, radius, theta)

  character*(*) lprofile_name
  character*128 :: filename

  integer, intent(inout) :: Nth, Nr
  double precision, intent(out) :: radius(Nr), theta(Nth)
  
  integer i,j
  
  !--------------------------------------------------------!       
  ! read Grid.dat profile      
  !--------------------------------------------------------!       
  filename = trim(adjustl(lprofile_name))//"/Grid.dat"
  open(666,file=trim(filename),status='unknown', & 
       form='formatted',action='read')
  read(666,*) Nth, Nr

  do i=1, Nth
     read(666,*) theta(i)
  enddo

  do i=1, Nr
     read(666,*) radius(i)
  enddo

  close(666)

end subroutine mod_XNS_get_grid

subroutine mod_XNS_get_profile(lprofile_name, Nr, Nth, rho, press, v3, psi, alp, beta3)

  character*(*) lprofile_name
  character*128 :: filename

  integer, intent(inout) :: Nth, Nr
  double precision, intent(out) :: rho(Nr,Nth),press(Nr,Nth),v3(Nr,Nth)
  double precision, intent(out) :: psi(Nr,Nth),alp(Nr,Nth),beta3(Nr,Nth)

  integer i,j
  
  double precision :: vtot2, discrim, rho_min
  double precision :: dx,dtheta,dphi, phi_max

  !write(*,*) "! Now read the Hydroeq.dat profile."
  !--------------------------------------------------------!       
  ! read Hydroeq.dat profile      
  !--------------------------------------------------------!       
  filename = trim(adjustl(lprofile_name))//"/Hydroeq.dat"
  open(666,file=trim(filename),status='unknown', & 
       form='formatted',action='read')
  read(666,*) Nth, Nr

  do j=1,Nth
     do i=1,Nr
        read(666,*) rho(i,j),press(i,j),psi(i,j), &
                    v3(i,j),alp(i,j),beta3(i,j)
     enddo
  enddo

  close(666)

end subroutine mod_XNS_get_profile

! **************************************************************
subroutine mod_XNS_linterp(x1,x2,y1,y2,x,y)

! perform linear interpolation      
  implicit none

  real*8 slope,x1,x2,y1,y2,x,y

  if (x2.lt.x1) then
     stop "Error in linterp!"
  endif

  slope = (y2 - y1) / (x2 - x1)

  y = slope*(x-x1) + y1
 
end subroutine  mod_XNS_linterp

! ***************************************************************

subroutine mod_XNS_find_index(zones,array,goal,upper_index,lower_index)
  
! bisection search
  implicit none
  
  integer zones,i
  real*8 array(*)
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
      
end subroutine mod_XNS_find_index

! ******************************************************************

subroutine mod_XNS_map(point_value,point_radius0,parray,pradius,zones)

  implicit none
  
  real*8 point_value, point_radius, point_radius0
  real*8 pradius(*), parray(*)
  integer zones
  integer upper_index, lower_index

  point_radius = abs(point_radius0)
  
  if (point_radius .ge. pradius(1) .and. & 
       point_radius .lt. pradius(zones) )  then
     
     call mod_XNS_find_index(zones,pradius,point_radius, &
          upper_index,lower_index)
     
     call mod_XNS_linterp( pradius(lower_index),pradius(upper_index), &
          parray(lower_index), parray(upper_index),  & 
          point_radius, point_value )

  else if (point_radius .lt. pradius(1)) then
     ! linear extrapolation
     call mod_XNS_linterp(pradius(1),pradius(2), & 
          parray(1),parray(2),point_radius,point_value)

  else if (point_radius .gt. pradius(zones)) then
     ! linear extrapolation
     call mod_XNS_linterp(pradius(zones-1),pradius(zones), & 
          parray(zones-1),parray(zones),point_radius,point_value)
  endif
  
  
end subroutine mod_XNS_map

end module mod_XNS

