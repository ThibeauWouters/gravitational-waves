module mod_eos_readtable
  use mod_eos
  use mod_eos_tabulated_parameters

  implicit none
  public

contains

   subroutine readtable(eos_filename)
   ! This routine reads the table and initializes
   ! all variables in the module. 
     use mod_global_parameters
     use mod_eos, only: &
         eos_rhomin,  eos_rhomax, &
         eos_yemin,   eos_yemax,  &
         eos_tempmin, eos_tempmax, &
         eos_epsmin,  eos_epsmax

     use hdf5 
   
     implicit none
   
     character(*), intent(in)      :: eos_filename
   
     character(len=100) :: message
   
   ! HDF5 vars
     integer(HID_T) :: file_id,dset_id,dspace_id
     integer(HSIZE_T) :: dims1(1), dims3(3)
     logical :: dataset_there
     integer :: error,rank,accerr
     integer :: i,j,k
   
     double precision :: amu_cgs_andi
     double precision :: buffer1,buffer2,buffer3,buffer4
     double precision, allocatable      :: enthalpy(:,:,:)
     accerr=0
   
     if (mype==0) write(*,*) "Reading stellarcollapse.org EOS Table"
   
     call h5open_f(error)
   
     call h5fopen_f (trim(adjustl(eos_filename)), H5F_ACC_RDONLY_F, file_id, error)
   
     if (mype==0) write(6,*) trim(adjustl(eos_filename))
   
   ! read scalars
     dims1(1)=1
     call h5dopen_f(file_id, "pointsrho", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nrho, dims1, error)
     call h5dclose_f(dset_id,error)
   
     if(error.ne.0) then
        stop "Could not read EOS table file"
     endif
   
     dims1(1)=1
     call h5dopen_f(file_id, "pointstemp", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntemp, dims1, error)
     call h5dclose_f(dset_id,error)
   
     if(error.ne.0) then
        stop "Could not read EOS table file"
     endif
   
     dims1(1)=1
     call h5dopen_f(file_id, "pointsye", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nye, dims1, error)
     call h5dclose_f(dset_id,error)
   
     if(error.ne.0) then
        stop "Could not read EOS table file"
     endif
   
     ! Check if we have one of the new SRO tables or if this is one of
     ! the old O'Connor & Ott tables on stellarcollapse.org
     ! The latter tables have the speed of sound not divided by the
     ! specific enthalpy
   
     ! first check if we have the relevant dataset
     call h5lexists_f(file_id,"have_rel_cs2",dataset_there,error)
     if(dataset_there) then
        dims1(1)=1
        call h5dopen_f(file_id, "have_rel_cs2", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, have_rel_cs2, dims1, error)
        call h5dclose_f(dset_id,error)
        if(error.ne.0) then
           stop "Could not read EOS table file"
        endif
     else
        have_rel_cs2 = 0
     endif
     
     if (mype==0) write(message,"(a25,i5,i5,i5)") "We have nrho ntemp nye: ", nrho,ntemp,nye
     if (mype==0) write(*,*) message
   
     allocate(eos_tables(nrho,ntemp,nye,nvars))
   
     ! index variable mapping:
   
     dims3(1)=nrho
     dims3(2)=ntemp
     dims3(3)=nye
     call h5dopen_f(file_id, "logpress", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_logpress), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
     call h5dopen_f(file_id, "logenergy", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_logenergy), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
     call h5dopen_f(file_id, "entropy", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_entropy), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
     call h5dopen_f(file_id, "munu", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_munu), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
     call h5dopen_f(file_id, "cs2", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_cs2), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
     call h5dopen_f(file_id, "dedt", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_dedt), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
     call h5dopen_f(file_id, "dpdrhoe", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_dpdrhoe), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
     call h5dopen_f(file_id, "dpderho", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_dpderho), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
   ! chemical potentials
     call h5dopen_f(file_id, "muhat", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_muhat), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     call h5dopen_f(file_id, "mu_e", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_mu_e), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     call h5dopen_f(file_id, "mu_p", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_mu_p), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     call h5dopen_f(file_id, "mu_n", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_mu_n), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
   ! compositions
     call h5dopen_f(file_id, "Xa", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_xa), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     call h5dopen_f(file_id, "Xh", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_xh), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     call h5dopen_f(file_id, "Xn", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_xn), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     call h5dopen_f(file_id, "Xp", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_xp), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
   
   ! average nucleus
     call h5dopen_f(file_id, "Abar", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_abar), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     call h5dopen_f(file_id, "Zbar", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_zbar), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
   ! Gamma
     call h5dopen_f(file_id, "gamma", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eos_tables(:,:,:,i_gamma), dims3, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     allocate(logrho_table(nrho))
     dims1(1)=nrho
     call h5dopen_f(file_id, "logrho", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logrho_table, dims1, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     allocate(logtemp_table(ntemp))
     dims1(1)=ntemp
     call h5dopen_f(file_id, "logtemp", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logtemp_table, dims1, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     allocate(ye_table(nye))
     dims1(1)=nye
     call h5dopen_f(file_id, "ye", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ye_table, dims1, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
   
     call h5dopen_f(file_id, "energy_shift", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, energy_shift, dims1, error)
     call h5dclose_f(dset_id,error)
     energy_shift = energy_shift * eps_gf ! to code unit
     accerr=accerr+error
   
     if(accerr.ne.0) then
       stop "Problem reading EOS table file"
     endif
   
   
     call h5fclose_f (file_id,error)
   
     call h5close_f (error)

     ! transform the table into code unit
     logrho_table(1:nrho) = logrho_table(1:nrho) + log10(rho_gf)
     eos_tables(1:nrho,1:ntemp,1:nye,i_logpress)  = eos_tables(1:nrho,1:ntemp,1:nye,i_logpress) + log10(press_gf)
     eos_tables(1:nrho,1:ntemp,1:nye,i_logenergy) = eos_tables(1:nrho,1:ntemp,1:nye,i_logenergy) + log10(eps_gf)
     eos_tables(1:nrho,1:ntemp,1:nye,i_cs2)       = eos_tables(1:nrho,1:ntemp,1:nye,i_cs2)/(clight)**2
     eos_tables(1:nrho,1:ntemp,1:nye,i_dpderho)   = eos_tables(1:nrho,1:ntemp,1:nye,i_dpderho) + log10(press_gf) - log10(eps_gf)
     eos_tables(1:nrho,1:ntemp,1:nye,i_dpdrhoe)   = eos_tables(1:nrho,1:ntemp,1:nye,i_dpdrhoe) + log10(press_gf) - log10(rho_gf)

     ! calculate enthalpy
     allocate(enthalpy(nrho,ntemp,nye))
     do k=1,nye
        do j=1,ntemp
           do i=1,nrho
              enthalpy(i,j,k) = (1.0d0 + 10.0d0**eos_tables(i,j,k,i_logenergy) - energy_shift + &
                   10.0d0**eos_tables(i,j,k,i_logpress)/10.0d0**logrho_table(i))
           enddo
        enddo
     enddo
   
     ! set min-max values:
     eos_rhomin = 10.0d0**logrho_table(1)
     eos_rhomax = 10.0d0**logrho_table(nrho)
   
     eos_yemin = ye_table(1)
     eos_yemax = ye_table(nye)
   
     eos_tempmin = 10.0d0**logtemp_table(1)
     eos_tempmax = 10.0d0**logtemp_table(ntemp)
   
     eos_epsmin = 10.0d0**minval(eos_tables(:,:,:,i_logenergy)) - energy_shift
     eos_epsmax = 10.0d0**maxval(eos_tables(:,:,:,i_logenergy)) - energy_shift

     eos_hmin = minval(enthalpy(:,:,:))
 
     ! correct legacy table speed of sound
     ! (divide by relativistic enthalpy, fix units)
     if(have_rel_cs2 == 0) &
        eos_tables(1:nrho,1:ntemp,1:nye,i_cs2) = eos_tables(1:nrho,1:ntemp,1:nye,i_cs2) / enthalpy(1:nrho,1:ntemp,1:nye)
   
     deallocate(enthalpy)
   
     if (mype==0) write(6,*) "Done reading eos tables"
   
   end subroutine readtable
   
end module 
