!> Module for eos
module mod_eos_tabulated
  use mod_eos
  use mod_eos_tabulated_parameters

  implicit none
  public

contains

  !> Read this module's parameters from a file
  subroutine eos_tabulated_read_params(files)
    use mod_eos
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /eos_tabulated_list/ eos_path, eos_precision

    do n = 1, size(files)
       ! Try to read in the namelists. They can be absent or in a different
       ! order, since we rewind before each read.
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, eos_tabulated_list, end=111)
111    close(unitpar)
    end do

  end subroutine eos_tabulated_read_params

  subroutine eos_tabulated_activate()
    use mod_global_parameters
    use mod_eos
    use mod_eos_readtable

    call eos_tabulated_read_params(par_files)
    eos_type = tabulated

    call readtable(eos_path)

    eos_get_pressure_one_grid         => tabulated_get_pressure_one_grid
    eos_get_eps_one_grid              => tabulated_get_eps_one_grid
    eos_get_cs2_one_grid              => tabulated_get_cs2_one_grid
    eos_get_temp_one_grid             => tabulated_get_temp_one_grid
    eos_get_beta_eqm_ye_one_grid      => tabulated_get_beta_eqm_ye_one_grid

    eos_eps_get_all_one_grid          => tabulated_eps_get_all_one_grid
    eos_temp_get_all_one_grid         => tabulated_temp_get_all_one_grid

    eos_get_eps_range                 => tabulated_get_eps_range

  end subroutine eos_tabulated_activate

  !> Get Pressure based on rho, eps, ye 
  subroutine tabulated_get_pressure_one_grid(prs,rho,eps,temp,ye)
    !use mod_eos
    use mod_interpolation
    implicit none
    
    double precision, intent(inout) :: prs
    double precision, intent(in)    :: rho
    double precision, intent(in)    :: eps
    double precision, intent(inout), optional :: temp
    double precision, intent(in), optional    :: ye

    double precision             :: log_rho, log_temp, log_eps, temp_in

    if (rho<small_rho_thr) then
       call atmo_get_pressure_one_grid(prs,rho,eps)
       return
    end if

    if (.not.present(ye))   call mpistop("nuc_eos get_press: need input ye!")

    if (rho>eos_rhomax)     call mpistop("nuc_eos get_press: rho > rhomax")
    if (rho<eos_rhomin)     call mpistop("nuc_eos get_press: rho < rhomin")
    if (ye>eos_yemax)       call mpistop("nuc_eos get_press: ye > yemax")
    if (ye<eos_yemin)       call mpistop("nuc_eos get_press: ye < yemin")
    if (eps>eos_epsmax)     call mpistop("nuc_eos get_press: eps > epsmax")
    if (eps<eos_epsmin)     call mpistop("nuc_eos get_press: eps < epsmin")

    ! if temp present, use it as an initial guess 
    if (present(temp)) temp_in = temp

    ! limit the initial guess of temperature
    temp_in = max(eos_tempmin, min(eos_tempmax, temp_in))

    log_rho  = dlog10(rho)
    log_eps  = dlog10(eps + energy_shift)
    log_temp = dlog10(temp_in)

    call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
    if (present(temp)) temp = 10.0d0**log_temp

    call intep3d(log_rho, log_temp, ye, &
           prs, eos_tables(:,:,:,i_logpress), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

    prs = 10.0d0**prs

  end subroutine tabulated_get_pressure_one_grid

  !> Get specific internal energy from (rho, temp, ye)
  subroutine tabulated_get_eps_one_grid(prs,rho,eps,temp,ye)
    !use mod_eos
    use mod_interpolation

    implicit none
    
    double precision, intent(in) :: prs
    double precision, intent(in) :: rho
    double precision, intent(in), optional :: temp, ye
    double precision, intent(inout) :: eps

    double precision             :: log_rho, log_temp

    if (.not.present(temp)) call mpistop("nuc_eos get_eps: need input temp!")
    if (.not.present(ye))   call mpistop("nuc_eos get_eps: need input ye!")

    if (rho<small_rho_thr) then
       call atmo_get_eps_one_grid(prs,rho,eps)
       return
    end if

    if (rho>eos_rhomax)     call mpistop("nuc_eos get_eps: rho > rhomax")
    if (rho<eos_rhomin)     call mpistop("nuc_eos get_eps: rho < rhomin")
    if (temp>eos_tempmax)   call mpistop("nuc_eos get_eps: temp > tempmax")
    if (temp<eos_tempmin)   call mpistop("nuc_eos get_eps: temp < tempmin")
    if (ye>eos_yemax)       call mpistop("nuc_eos get_eps: ye > yemax")
    if (ye<eos_yemin)       call mpistop("nuc_eos get_eps: ye < yemin")

    log_rho  = dlog10(rho)
    log_temp = dlog10(temp)

    call intep3d(log_rho, log_temp, ye, &
           eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

    eps = 10.0d0**eps - energy_shift

  end subroutine tabulated_get_eps_one_grid

  !> Get cs2 from (rho, eps, ye)
  subroutine tabulated_get_cs2_one_grid(cs2,rho,eps,temp,ye)
    !use mod_eos
    use mod_interpolation

    implicit none
    
    double precision, intent(inout) :: cs2
    double precision, intent(in) :: rho
    double precision, intent(in) :: eps
    double precision, intent(in), optional :: temp, ye

    double precision             :: log_rho, log_eps, log_temp

    if (.not.present(temp)) call mpistop("nuc_eos get_cs2: need input temp!")
    if (.not.present(ye))   call mpistop("nuc_eos get_cs2: need input ye!")

    if (rho<small_rho_thr) then
       call atmo_get_cs2_one_grid(cs2,rho,eps)
       return
    end if

    if (rho>eos_rhomax)     call mpistop("nuc_eos get_cs2: rho > rhomax")
    if (rho<eos_rhomin)     call mpistop("nuc_eos get_cs2: rho < rhomin")
    if (temp>eos_tempmax)   call mpistop("nuc_eos get_cs2: temp > tempmax")
    if (temp<eos_tempmin)   call mpistop("nuc_eos get_cs2: temp < tempmin")
    if (ye>eos_yemax)       call mpistop("nuc_eos get_cs2: ye > yemax")
    if (ye<eos_yemin)       call mpistop("nuc_eos get_cs2: ye < yemin")
    if (eps>eos_epsmax)     call mpistop("nuc_eos get_cs2: eps > epsmax")
    if (eps<eos_epsmin)     call mpistop("nuc_eos get_cs2: eps < epsmin")

    log_rho  = dlog10(rho)
    log_temp = dlog10(temp)
    log_eps  = dlog10(eps + energy_shift)

    call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)

    call intep3d(log_rho, log_temp, ye, &
           cs2, eos_tables(:,:,:,i_cs2), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)

  end subroutine tabulated_get_cs2_one_grid

  !> Get Temperature from (rho, eps, ye)
  subroutine tabulated_get_temp_one_grid(rho,eps,temp,ye)
    !use mod_eos
    use mod_interpolation

    implicit none
    
    double precision, intent(in)    :: rho, eps, ye
    double precision, intent(inout) :: temp

    double precision                :: log_rho, log_temp, log_eps

    if (rho<small_rho_thr) then
       temp = atmo_temp
       return
    end if

    if (rho>eos_rhomax)   call mpistop("nuc_eos get_temp: rho > rhomax"  ) 
    if (rho<eos_rhomin)   call mpistop("nuc_eos get_temp: rho < rhomin"  ) 
    if (ye>eos_yemax)     call mpistop("nuc_eos get_temp: ye > yemax"    ) 
    if (ye<eos_yemin)     call mpistop("nuc_eos get_temp: ye < yemin"    ) 
    if (eps>eos_epsmax)   call mpistop("nuc_eos get_temp: eps > epsmax"  ) 
    if (eps<eos_epsmin)   call mpistop("nuc_eos get_temp: eps < epsmin"  ) 

    ! limit the initial guess of temperature
    temp = max(eos_tempmin, min(eos_tempmax, temp))

    log_rho = dlog10(rho)
    log_eps = dlog10(eps + energy_shift)
    log_temp = dlog10(temp)

    call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)

    temp = 10.0d0**log_temp

    ! check the results
    if (temp>eos_tempmax) call mpistop("nuc_eos get_temp: temp > tempmax") 
    if (temp<eos_tempmin) call mpistop("nuc_eos get_temp: temp < tempmin") 

  end subroutine tabulated_get_temp_one_grid

  !> Get all tab eos var from (rho, eps, ye)
  subroutine tabulated_eps_get_all_one_grid(rho,eps,ye,temp,prs,ent,cs2,dedt,&
                 dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu)
    use mod_interpolation
    implicit none

    double precision, intent(inout) :: rho, eps, ye
    double precision, intent(inout), optional :: ent, prs, temp, cs2, dedt
    double precision, intent(inout), optional :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar
    double precision, intent(inout), optional :: mu_e,mu_n,mu_p,muhat,munu

    double precision             :: log_rho, log_temp, log_eps, temp_in

    double precision             :: ffx(nvars)

    if (rho<small_rho_thr) then
       rho = small_rho
       eps = small_eps
       ye = atmo_ye
       ! fixme: should include temp and rest of the variables as well
       if (present(prs))    call atmo_get_pressure_one_grid(prs,rho,eps)
       if (present(cs2))    call atmo_get_cs2_one_grid(cs2,rho,eps)
       if (present(ent))    ent  = 4.0d0
       if (present(temp))   temp = atmo_temp
       return
    end if

    if (rho>eos_rhomax) call mpistop("nuc_eos eps_get_all: rho > rhomax")
    if (ye>eos_yemax)   call mpistop("nuc_eos eps_get_all: ye  > yemax ")
    if (ye<eos_yemin)   call mpistop("nuc_eos eps_get_all: ye  < yemin ")
    if (eps>eos_epsmax) call mpistop("nuc_eos eps_get_all: eps > epsmax")
    if (eps<eos_epsmin) call mpistop("nuc_eos eps_get_all: eps < epsmin")
    
    temp_in = eos_tempmin
    if (present(temp)) temp_in = temp

    log_temp = dlog10(temp_in) 
    log_rho  = dlog10(rho)
    log_eps  = dlog10(eps + energy_shift)

    call tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)

    call intep3d_many(log_rho, log_temp, ye, &
                      ffx, eos_tables,  &
                      nrho, ntemp, nye, nvars,     &
                      logrho_table, logtemp_table, ye_table)

    if (present(temp))     temp = 10.d0**log_temp
    if (present(prs))      prs  = 10.d0**ffx(i_logpress)
    if (present(ent))      ent  = ffx(i_entropy)
    if (present(munu))     munu = ffx(i_munu)
    !if (ffx(i_cs2)> 1.0d0) write(*,*) "nuc_eos:cs2>1.0d0" !stop "nuc_eos: cs2 > 1.0d0"
    !if (ffx(i_cs2)< 0.0d0) write(*,*) "nuc_eos:cs2<1.0d0" !stop "nuc_eos: cs2 < 0.0d0"
    if (present(cs2))      cs2  = ffx(i_cs2)
    !  derivatives
    if (present(dedt))     dedt    = ffx(i_dedT)
    if (present(dpdrhoe))  dpdrhoe = ffx(i_dpdrhoe)
    if (present(dpderho))  dpderho = ffx(i_dpderho)
    !  chemical potentials
    if (present(muhat))    muhat   = ffx(i_muhat)
    if (present(mu_e))     mu_e    = ffx(i_mu_e)
    if (present(mu_p))     mu_p    = ffx(i_mu_p)
    if (present(mu_n))     mu_n    = ffx(i_mu_n)
    !  compositions
    if (present(xa))       xa      = ffx(i_xa)
    if (present(xh))       xh      = ffx(i_xh)
    if (present(xn))       xn      = ffx(i_xn)
    if (present(xp))       xp      = ffx(i_xp)
    if (present(abar))     abar    = ffx(i_abar)
    if (present(zbar))     zbar    = ffx(i_zbar)

  end subroutine tabulated_eps_get_all_one_grid

  !> Get all eos tab var from (rho, temp, ye)
  subroutine tabulated_temp_get_all_one_grid(rho,temp,ye,eps,prs,ent,cs2,dedt,&
                       dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,mu_e,mu_n,mu_p,muhat,munu)
    use mod_interpolation
    implicit none

    double precision, intent(in)    :: ye
    double precision, intent(inout) :: eps, rho, temp
    double precision, intent(inout), optional :: ent, prs, cs2, dedt
    double precision, intent(inout), optional :: dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar
    double precision, intent(inout), optional :: mu_e,mu_n,mu_p,muhat,munu

    double precision             :: log_rho, log_temp, log_eps
    double precision             :: ffx(nvars)

    if (rho<small_rho_thr) then
       ! fixme: should include rest of the variables as well
       if (present(prs))    call atmo_get_pressure_one_grid(prs,rho,eps)
       if (present(cs2))    call atmo_get_cs2_one_grid(cs2,rho,eps)
       if (present(ent)) ent = 4.0d0
       return
    end if

    if (rho>eos_rhomax) call mpistop("nuc_eos temp_get_all: rho > rhomax")
    if (ye>eos_yemax)   call mpistop("nuc_eos temp_get_all: ye  > yemax ")
    if (ye<eos_yemin)   call mpistop("nuc_eos temp_get_all: ye  < yemin ")
    if (eps>eos_epsmax) call mpistop("nuc_eos temp_get_all: eps > epsmax")
    if (eps<eos_epsmin) call mpistop("nuc_eos temp_get_all: eps < epsmin")

    log_temp = dlog10(temp)   
    log_rho  = dlog10(rho)

    call intep3d_many(log_rho, log_temp, ye, &
                      ffx, eos_tables,  &
                      nrho, ntemp, nye, nvars,     &
                      logrho_table, logtemp_table, ye_table)

    eps  = 10.0d0**ffx(i_logenergy) - energy_shift
    if (present(prs))      prs  = 10.d0**ffx(i_logpress)
    if (present(ent))      ent  = ffx(i_entropy)
    if (present(munu))     munu = ffx(i_munu)
    !if (ffx(i_cs2)> 1.0d0) write(*,*) "nuc_eos:cs2>1.0d0" !stop "nuc_eos: cs2 > 1.0d0"
    !if (ffx(i_cs2)< 0.0d0) write(*,*) "nuc_eos:cs2<1.0d0" !stop "nuc_eos: cs2 < 0.0d0"
    if (present(cs2))      cs2  = ffx(i_cs2)
    !  derivatives
    if (present(dedt))     dedt    = ffx(i_dedT)
    if (present(dpdrhoe))  dpdrhoe = ffx(i_dpdrhoe)
    if (present(dpderho))  dpderho = ffx(i_dpderho)
    !  chemical potentials
    if (present(muhat))    muhat   = ffx(i_muhat)
    if (present(mu_e))     mu_e    = ffx(i_mu_e)
    if (present(mu_p))     mu_p    = ffx(i_mu_p)
    if (present(mu_n))     mu_n    = ffx(i_mu_n)
    !  compositions
    if (present(xa))       xa      = ffx(i_xa)
    if (present(xh))       xh      = ffx(i_xh)
    if (present(xn))       xn      = ffx(i_xn)
    if (present(xp))       xp      = ffx(i_xp)
    if (present(abar))     abar    = ffx(i_abar)
    if (present(zbar))     zbar    = ffx(i_zbar)

  end subroutine tabulated_temp_get_all_one_grid

  subroutine tabulated_get_eps_range(rho, eps_min, eps_max, ye)
    use mod_interpolation
    implicit none

    double precision, intent(in) :: rho
    double precision, intent(in), optional :: ye
    double precision, intent(out) :: eps_max, eps_min
    double precision              :: log_rho, log_temp, eps

    if (.not.present(ye))   call mpistop("get_eps_range: Tabulate EOS need ye!")

    !find eps_min    
    log_rho  = dlog10(rho)
    log_temp = dlog10(eos_tempmin)
    
    call intep3d(log_rho, log_temp, ye, &
                eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
    eps_min = 10.d0**eps - energy_shift

    !find eps_max
    log_rho  = dlog10(rho)
    log_temp = dlog10(eos_tempmax)
    
    call intep3d(log_rho, log_temp, ye, &
                eps, eos_tables(:,:,:,i_logenergy), &
           nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
    eps_max = 10.d0**eps - energy_shift

    ! Debug
    !eps_min = eos_epsmin
    !eps_max = eos_epsmax

  end subroutine tabulated_get_eps_range  

  subroutine tabulated_logeps_to_logtemp(log_rho,log_eps,log_temp,ye)
    use mod_interpolation
    use mod_rootfinding

    implicit none
    
    double precision, intent(in)    :: log_rho, log_eps, ye
    double precision, intent(inout) :: log_temp

    double precision                :: logtemp_min, logtemp_max
    double precision                :: log_eps_local, log_eps_min, log_eps_max
    ! for root finding
    integer                         :: error_code

    error_code = -1
    log_eps_local = log_eps

    ! range of the root
    logtemp_min = logtemp_table(1)
    logtemp_max = logtemp_table(ntemp)
    !log_temp = max(logtemp_min, min(logtemp_max, log_temp))

    ! check if the initial guess is close enough
    !if (dabs(func_eps_of_temp(log_temp)) < eos_precision * dabs(log_eps_local) ) return

    ! get logtemp from logeps
    !call rootfinding_illinois(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
    call rootfinding_brent(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
    select case (error_code)
    !case (0) ! root is found
    !   write(*,*) "z= ", z
    case (-1) ! nothing happened
       call mpistop("have you ever attemp to find the root in tabulated eos?")
    case (2) ! z is NaN
       call mpistop("NaN")
    case (1,3) ! if log_temp cannot be found or is not bracketed
       ! adjust the range of log_eps, find the root again
       error_code = -1
       call eos_get_eps_range(10.0d0**log_rho, log_eps_min, log_eps_max, ye)
       log_eps_max = dlog10(log_eps_max+energy_shift)
       log_eps_min = dlog10(log_eps_min+energy_shift)
       log_eps_local = max(log_eps_min, min(log_eps_max, log_eps_local))
       ! using the bound values of logeps to find logtemp again
       !call rootfinding_illinois(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
       call rootfinding_brent(log_temp, logtemp_min, logtemp_max, eos_precision, eos_iter_max, error_code, func_eps_of_temp)
       select case (error_code)
       !case (0) ! root is found
       !   write(*,*) "z= ", z
       case (-1) ! nothing happened
          call mpistop("forgot to find the root in tabulated eos after adjusting eps range?")
       case (1) ! failed to find the root
          call mpistop("Fail to find the root in tabulated eos after adjusting eps range")
       case (2) ! z is NaN
          call mpistop("NaN after adjusting eps range")
       case (3) ! z is not bracketed
          !write(*,*)  "log_eps, logeps_min, logeps_max, rho, ye "
          !write(*,*)  log_eps_local, log_eps_min, log_eps_max, 10**log_rho, ye
          !write(*,*) log_temp, logtemp_min, logtemp_max, func_eps_of_temp(logtemp_min), func_eps_of_temp(logtemp_max)
          !write(*,*)  "---------------------------"
          ! If still not bracketed, probably because they are out of the range
          ! fixme: maybe there are better ways to handle this case
          if (log_eps_local <= log_eps_min) then
             log_temp = logtemp_min
          else if (log_eps_local >= log_eps_max) then
             log_temp = logtemp_max
          else
             call mpistop("the root is not bracketed in tabulated eos after adjusting eps range")
          end if
          !do error_code = 1,ntemp
          !   write(*,*) logtemp_table(error_code), func_eps_of_temp(logtemp_table(error_code))
          !end do
          !call mpistop("the root is not bracketed in tabulated eos after adjusting eps range")
       end select
    end select

    contains
       double precision function func_eps_of_temp(log_temp_in)
         implicit none
         double precision, intent(in) :: log_temp_in
         call intep3d(log_rho, log_temp_in, ye, & 
                  func_eps_of_temp, eos_tables(:,:,:,i_logenergy), &
                   nrho, ntemp, nye, logrho_table, logtemp_table, ye_table)
         func_eps_of_temp = 1.0d0 - func_eps_of_temp / log_eps_local 
       end function
  end subroutine tabulated_logeps_to_logtemp

  subroutine tabulated_get_beta_eqm_ye_one_grid(rho,temp,ye)
    use mod_interpolation
    use mod_rootfinding
    implicit none

    double precision, intent(in)              :: rho, temp
    double precision, intent(inout)           :: ye

    double precision                          :: log_rho, log_temp, log_eps, ye_min, ye_max
    integer                                   :: error_code = -1

    ! given rho and temp, get Ye which satisfies beta equilibrium
    if (rho>eos_rhomax)     call mpistop("nuc_eos beta_eqm: rho > rhomax")
    if (rho<eos_rhomin)     call mpistop("nuc_eos beta_eqm: rho < rhomin")
    if (temp>eos_tempmax)   call mpistop("nuc_eos beta_eqm: temp > tempmax")
    if (temp<eos_tempmin)   call mpistop("nuc_eos beta_eqm: temp < tempmin")

    ! Beta equilibrium requires that \mu_n =  \mu_p +\mu_e -\mu_nu
    ! Since we assume that \mu_nu should be negligible we demand
    ! \mu_n-\mu_p-\mu_e = \mu_hat =0

    log_rho  = dlog10(rho)
    log_temp = dlog10(temp)

    ye_min = ye_table(1)
    ye_max = ye_table(nye)

    ! get eps from temperature
    call rootfinding_brent(ye, ye_min, ye_max, eos_precision, eos_iter_max, error_code, func_munu_of_ye)
    select case (error_code)
    !case (0) ! root is found
    !   write(*,*) "z= ", z
    case (-1) ! nothing happened
       call mpistop("have you ever attemp to find the root in tabulated eos?")
    case (1) ! failed to find the root
       call mpistop("Fail to find the root in beta_eq")
    case (2) ! z is NaN
       call mpistop("NaN in beta_eq")
    case (3) ! z is not bracketed
       call mpistop("the root is not bracketed in beta_eq")
    end select

    contains
       double precision function func_munu_of_ye(ye_in)
         implicit none
         double precision, intent(in) :: ye_in
         double precision             :: ffx(nvars)
         call intep3d_many(log_rho, log_temp, ye_in, &
                      ffx, eos_tables,  &
                      nrho, ntemp, nye, nvars,     &
                      logrho_table, logtemp_table, ye_table)
         !f(Ye) =  mu_e(Ye) + mu_p(Ye) - mu_n(Ye)  !ignored the restmass difference
         func_munu_of_ye = ffx(i_mu_e) + ffx(i_mu_p) - ffx(i_mu_n)
       end function
  end subroutine tabulated_get_beta_eqm_ye_one_grid
 
end module mod_eos_tabulated
