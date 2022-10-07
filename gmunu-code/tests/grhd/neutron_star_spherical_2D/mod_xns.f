module mod_XNS
  implicit none
  private

  type XNS_profile
     integer                       :: Nr, Nth
     double precision, allocatable :: radius(:), theta(:)
     double precision, allocatable :: rho(:,:), press(:,:), v3(:,:)
     double precision, allocatable :: psi(:,:), alp(:,:), beta3(:,:)
     double precision, allocatable :: B3(:,:), B1(:,:), B2(:,:)
     double precision, allocatable :: J3(:,:), J1(:,:), J2(:,:)
     double precision, allocatable :: A3(:,:), At(:,:)
     double precision, allocatable :: E_3(:,:), E_1(:,:), E_2(:,:)
     double precision, allocatable :: r_star(:)
  end type XNS_profile

  ! save the profle
  type(XNS_profile) :: prof

  public :: XNS_profile
  public :: prof
  public :: mod_XNS_allocate_var
  public :: mod_XNS_deallocate_var
  public :: mod_XNS_read_profile
  public :: mod_XNS_map_1D
  public :: mod_XNS_map_2D

contains

subroutine mod_XNS_allocate_var(Nr, Nth)
  integer, intent(in) :: Nth, Nr
  allocate(prof%radius(Nr))
  allocate(prof%theta(Nth))
  allocate(prof%rho(Nr,Nth),prof%press(Nr,Nth),prof%v3(Nr,Nth))
  allocate(prof%psi(Nr,Nth),prof%alp(Nr,Nth),prof%beta3(Nr,Nth))
  allocate(prof%B3(Nr,Nth),prof%B1(Nr,Nth),prof%B2(Nr,Nth))
  allocate(prof%A3(Nr,Nth),prof%At(Nr,Nth))
  allocate(prof%E_3(Nr,Nth),prof%E_1(Nr,Nth),prof%E_2(Nr,Nth))
  allocate(prof%J3(Nr,Nth),prof%J1(Nr,Nth),prof%J2(Nr,Nth))
  allocate(prof%r_star(Nth))
end subroutine mod_XNS_allocate_var

subroutine mod_XNS_deallocate_var
  deallocate(prof%radius)
  deallocate(prof%theta)
  deallocate(prof%rho,prof%press,prof%v3)
  deallocate(prof%psi,prof%alp,prof%beta3)
  deallocate(prof%B3,prof%B1,prof%B2)
  deallocate(prof%A3,prof%At)
  deallocate(prof%E_3,prof%E_1,prof%E_2)
  deallocate(prof%J3,prof%J1,prof%J2)
  deallocate(prof%r_star)
end subroutine mod_XNS_deallocate_var

subroutine mod_XNS_read_profile(profile_path)

  character*(*) profile_path
  character*128 :: filename
  double precision :: minrho

  integer :: Nth, Nr
  integer :: i,j

  call mod_XNS_get_dim(profile_path, Nr, Nth)

  prof%Nr = Nr
  prof%Nth = Nth
  call mod_XNS_allocate_var(Nr, Nth)

  ! read profile
  call mod_XNS_get_grid(profile_path, Nr, Nth, prof%radius, prof%theta)
  call mod_XNS_get_profile(profile_path, Nr, Nth, &
                             prof%rho,prof%press,prof%v3,prof%psi,prof%alp,prof%beta3, &
                             prof%B3, prof%B2, prof%B1, &
                             prof%A3, &
                             prof%E_3, prof%E_2, prof%E_1, &
                             prof%At, &
                             prof%J3, prof%J2, prof%J1 &
                             )


  minrho = minval(prof%rho)

  prof%r_star(:) = HUGE(0.0d0) ! initiality set r_star at infinity

  do j=1,Nth
     do i=1,Nr
        ! set atmosphere
        if(prof%rho(i,j)<=minrho) then
           prof%rho(i,j) = TINY(0.0d0)
           prof%r_star(j) = min( prof%r_star(j), prof%radius(i) )
        end if
        ! A_tilde = A / ( r sin(theta) ), so transform it back
        prof%A3(i,j) = prof%A3(i,j) * prof%radius(i) * dsin(prof%theta(j))
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

subroutine mod_XNS_get_profile(lprofile_name, Nr, Nth, rho, press, v3, &
               psi, alp, beta3, &
               B3, B2, B1, A3, E_3, E_2, E_1, At, J3, J2, J1)

  character*(*) lprofile_name
  character*128 :: filename

  integer, intent(inout) :: Nth, Nr
  double precision, intent(out) :: rho(Nr,Nth),press(Nr,Nth),v3(Nr,Nth)
  double precision, intent(out) :: psi(Nr,Nth),alp(Nr,Nth),beta3(Nr,Nth)
  double precision, intent(out) :: B3(Nr,Nth), B1(Nr,Nth), B2(Nr,Nth)
  double precision, intent(out) :: A3(Nr,Nth), At(Nr,Nth)
  double precision, intent(out) :: J3(Nr,Nth), J1(Nr,Nth), J2(Nr,Nth)
  double precision, intent(out) :: E_3(Nr,Nth), E_1(Nr,Nth), E_2(Nr,Nth)

  integer i,j
  
  double precision :: vtot2, discrim, rho_min
  double precision :: dx,dtheta,dphi, phi_max
  
  logical :: mag_exist = .False.

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

  ! Here, we need to check if B field exists
  filename = trim(adjustl(lprofile_name))//"/Hydroeq_mag.dat"
  inquire(file=filename, exist=mag_exist)
  if (mag_exist) then
     !--------------------------------------------------------!       
     ! read Hydroeq_mag.dat profile      
     !--------------------------------------------------------!       
     !filename = trim(adjustl(lprofile_name))//"/Hydroeq_mag.dat"
     open(666,file=trim(filename),status='unknown', & 
          form='formatted',action='read')
     read(666,*) Nth, Nr
   
     do j=1,Nth
        do i=1,Nr
           read(666,*) B3(i,j), B1(i,j), B2(i,j), &
                       A3(i,j), &
                       E_3(i,j), E_1(i,j), E_2(i,j), &
                       At(i,j), &
                       J3(i,j), J1(i,j), J2(i,j)
        enddo
     enddo
   
     close(666)
  end if

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
  ir = min( max(ir, 2), nr-1 )
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
  ! make sure the index falls into correct range
  ir = min( max(ir, 2), nr-1 )
  itheta = min( max(itheta, 2), ntheta-1 )

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

