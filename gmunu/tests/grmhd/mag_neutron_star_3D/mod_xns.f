module mod_XNS
  implicit none
  private

  type XNS_profile
     integer                       :: Nr, Nth
     double precision, allocatable :: radius(:), theta(:)
     double precision, allocatable :: rho(:,:), press(:,:), v3(:,:)
     double precision, allocatable :: psi(:,:), alp(:,:), beta3(:,:)
     double precision, allocatable :: B3(:,:)
  end type XNS_profile

  ! save the profle
  type(XNS_profile) :: prof

  public :: XNS_profile
  public :: prof
  public :: mod_XNS_read_profile
  public :: mod_XNS_map_1D
  public :: mod_XNS_map_2D

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
  allocate(prof%B3(Nr,Nth))
  
  ! read profile
  call mod_XNS_get_grid(profile_path, Nr, Nth, prof%radius, prof%theta)
  call mod_XNS_get_profile(profile_path, Nr, Nth, &
                             prof%rho,prof%press,prof%v3,prof%psi,prof%alp,prof%beta3,prof%b3)


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

subroutine mod_XNS_get_profile(lprofile_name, Nr, Nth, rho, press, v3, psi, alp, beta3, b3)

  character*(*) lprofile_name
  character*128 :: filename

  integer, intent(inout) :: Nth, Nr
  double precision, intent(out) :: rho(Nr,Nth),press(Nr,Nth),v3(Nr,Nth),b3(Nr,Nth)
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
                    v3(i,j),alp(i,j),beta3(i,j),b3(i,j)
     enddo
  enddo

  close(666)

end subroutine mod_XNS_get_profile

subroutine mod_XNS_map_1D(point_value, point_radius, parray, pradius ,nr)

  implicit none
  
  double precision, intent(out)  :: point_value
  double precision, intent(in)   :: point_radius
  double precision, dimension(1:nr), intent(in) :: parray
  double precision, dimension(1:nr), intent(in) :: pradius
  integer, intent(in)   :: nr 

  integer :: ir
  integer :: pm1r
  double precision   :: dr
  double precision   :: del_r

  ! find the closest points
  ir = minloc( dabs(pradius(2:nr-1)-point_radius), dim = 1 )
  del_r = point_radius - pradius(ir)
  if (del_r > 0.0d0) then
     pm1r = 1
  else
     pm1r = -1
  end if
 
  dr = ( pradius(ir + pm1r) - pradius(ir) )
  point_value = parray(ir) &
                + del_r/dr * (parray(ir + pm1r) - parray(ir))

end subroutine mod_XNS_map_1D

subroutine mod_XNS_map_2D(point_value, point_radius, point_theta,&
                             parray, &
                             pradius ,ptheta ,nr, ntheta)

  implicit none
  
  double precision, intent(out)  :: point_value
  double precision, intent(in)   :: point_radius, point_theta
  double precision, dimension(1:nr,1:ntheta), intent(in) :: parray
  double precision, dimension(1:nr), intent(in) :: pradius
  double precision, dimension(1:ntheta), intent(in) :: ptheta
  integer, intent(in)   :: nr, ntheta

  integer :: ir, itheta
  integer :: pm1r, pm1theta
  double precision   :: dr, dtheta
  double precision   :: del_r, del_theta

  double precision   :: fh(4), a(4)

  ! find the closest points
  ir = minloc( dabs(pradius(2:nr-1)-point_radius), dim = 1 )
  itheta = minloc( dabs(ptheta(2:ntheta-1)-point_theta), dim = 1 )

  del_r = point_radius - pradius(ir)
  del_theta = point_theta - ptheta(itheta)

  if (del_r >= 0.0d0) then
     pm1r = 1
  else
     pm1r = -1
  end if
  if (del_theta >= 0.0d0) then
     pm1theta = 1
  else
     pm1theta = -1
  end if
 
  dr = ( pradius(ir + pm1r) - pradius(ir) )
  dtheta = ( ptheta(itheta + pm1theta) - ptheta(itheta) )

  fh(1) = parray(ir, itheta)
  fh(2) = parray(ir + pm1r, itheta)
  fh(3) = parray(ir, itheta + pm1theta)
  fh(4) = parray(ir + pm1r, itheta + pm1theta)

  a(1) = fh(1)
  a(2) = (fh(2)-fh(1))/dr
  a(3) = (fh(3)-fh(1))/dtheta
  a(4) = (fh(4)-fh(2)-fh(3)+fh(1))/dtheta/dr

  point_value  = a(1) &
                 + a(2) * del_r &
                 + a(3) * del_theta &
                 + a(4) * del_theta * del_r
end subroutine mod_XNS_map_2D

end module mod_XNS

