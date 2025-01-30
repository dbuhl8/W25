! This tool reads a netCDF simulation data field and returns the uw and mean u

! The program needs 1+ arguments:
! The name(s) of the simdat file

program netcdfextract
  use netcdf
  implicit none

  character(len=*),Parameter   :: errstring = 'Usage: netcdf2raw infile'
  character(len=100) :: netcdf_file
  character(len=100) :: arg
  character(len=34)  :: fname
  integer :: argnum,arglength,ierr
  real    :: t,dt
  integer :: ncid
  integer :: nfields
  logical :: notemp = .false.,nochem = .false.
  logical :: noux = .false.,nouy = .false.,nouz = .false.
  integer :: Temp_varid,Chem_varid
  integer :: ux_varid,uy_varid,uz_varid
  integer :: t_varid,dt_varid,timestep_varid
  integer :: x_dimid,y_dimid,z_dimid,time_dimid
  integer :: B_therm_varid,B_comp_varid,D_visc_varid,D_therm_varid,D_comp_varid
  integer :: S_therm_varid,S_comp_varid
  integer :: Gammax_varid,Gammay_varid,Gammaz_varid
  real    :: B_therm,B_comp,D_visc,D_therm,D_comp,S_therm,S_comp
  real    :: Gammax,Gammay,Gammaz
  real    :: dx,dy,dz
  integer :: Nx,Ny,Nz,nsteps,nstep
  integer :: starts(4),counts(4)
  integer :: n,k
  integer :: timestep,firststep,laststep
  integer :: counter ! Db
  character(len=12) :: cdummy

  real, dimension(:,:,:), allocatable :: ux,uy,uz,wz,temp,tdisp,tot_energy,mdisp
  real, dimension(:,:), allocatable :: ux_bar,uy_bar,uz_bar,uz_rms,mdisp_bar
  real, dimension(:,:), allocatable :: growth_rate
  real, dimension(:,:), allocatable :: tdisp_bar, barotropic_ER,lz_bar,vortz_bar
  real :: invRo=2 ! this needs to be changed for each simulation it is ran on
  complex :: c1 = (1.0, 0.0), ci = (0.0, 1.0)
  integer :: i1,j1,k1,ncount, num_files

  open(55,file='vort_plot.dat')

  argnum = COMMAND_ARGUMENT_COUNT() !Fortran 2003, but most compilers know 
  
  ! Dante - modified the command argument logic and stucture so that this may
  ! work on more than 1 simdat file

  if (argnum.ne.1) then 
    do num_files = 1, argnum
      CALL GET_COMMAND_ARGUMENT(num_files,netcdf_file)
      write(*,'(2a)') "Processing ",trim(netcdf_file)
      write(*,'(a,3I5)') "Looking for timestep",nstep

      call check( nf90_open(netcdf_file, nf90_nowrite, ncid) )
      call check( nf90_inq_varid(ncid,'t',t_varid) ) 
      call check( nf90_inq_varid(ncid,'timestep',timestep_varid) ) 
      ! Looks to see if the field exists, modifies flag if not.  
      call check_for_field(ux_varid,'ux',noux)
      call check_for_field(uy_varid,'uy',nouy)
      call check_for_field(uz_varid,'uz',nouz)
      call check_for_field(temp_varid,'Temp',notemp)
      ! Reads in X, Y, and TIME, as well as their dimensions.  
      call check(nf90_inq_dimid(ncid,"X",x_dimid) )
      call check(nf90_inquire_dimension(ncid,x_dimid,len=Nx) )
      call check(nf90_inq_dimid(ncid,"Y",y_dimid) )
      call check(nf90_inquire_dimension(ncid,y_dimid,len=Ny) )
      call check(nf90_inq_dimid(ncid,"Z",z_dimid) )
      call check(nf90_inquire_dimension(ncid,z_dimid,len=Nz) )
      call check(nf90_inq_dimid(ncid,"TIME",time_dimid) )
      call check(nf90_inquire_dimension(ncid,time_dimid,len=nsteps) )

      Write(*,'(a,i4,a,i4,a,i4)') "Nx =",Nx," Ny =",Ny," Nz =",Nz
      Write(*,'(a,i7)') "Number of saved timesteps =",nsteps
      ! Reads in the input parameters 
      call check(nf90_inq_varid(ncid,"B_therm", B_therm_varid) )
      call check(nf90_get_var(ncid,B_therm_varid,B_therm) )
      call check(nf90_inq_varid(ncid,"B_comp", B_comp_varid) )
      call check(nf90_get_var(ncid,B_comp_varid,B_comp) )
      call check(nf90_inq_varid(ncid,"D_visc", D_visc_varid) )
      call check(nf90_get_var(ncid,D_visc_varid,D_visc) )
      call check(nf90_inq_varid(ncid,"D_therm", D_therm_varid) )
      call check(nf90_get_var(ncid,D_therm_varid,D_therm) )
      call check(nf90_inq_varid(ncid,"D_comp", D_comp_varid) )
      call check(nf90_get_var(ncid,D_comp_varid,D_comp) )
      call check(nf90_inq_varid(ncid,"S_therm", S_therm_varid) )
      call check(nf90_get_var(ncid,S_therm_varid,S_therm) )
      call check(nf90_inq_varid(ncid,"S_comp", S_comp_varid) )
      call check(nf90_get_var(ncid,S_comp_varid,S_comp) )
      call check(nf90_inq_varid(ncid,"Gammax", Gammax_varid) )
      call check(nf90_get_var(ncid,Gammax_varid,Gammax) )
      call check(nf90_inq_varid(ncid,"Gammay", Gammay_varid) )
      call check(nf90_get_var(ncid,Gammay_varid,Gammay) )
      call check(nf90_inq_varid(ncid,"Gammaz", Gammaz_varid) )
      call check(nf90_get_var(ncid,Gammaz_varid,Gammaz) )
      print*,"halli"

      dx = Gammax/Nx
      dy = Gammay/Ny
      dz = Gammaz/Nz
      ! Allocate arrays for ux, uy, uz
      allocate(ux(0:Nx-1,0:Ny-1,0:Nz-1))
      allocate(uy(0:Nx-1,0:Ny-1,0:Nz-1))
      allocate(uz(0:Nx-1,0:Ny-1,0:Nz-1))
      allocate(wz(0:Nx-1,0:Ny-1,0:Nz-1))
      allocate(temp(0:Nx-1,0:Ny-1,0:Nz-1))
      allocate(tdisp(0:Nx-1,0:Ny-1,0:Nz-1))
      allocate(mdisp(0:Nx-1,0:Ny-1,0:Nz-1))
      allocate(tot_energy(0:Nx-1,0:Ny-1,0:Nz-1))

      allocate(vortz_bar(0:Nx-1,0:Ny-1))
      allocate(uz_bar(0:Nx-1,0:Ny-1))
      allocate(ux_bar(0:Nx-1,0:Ny-1))
      allocate(uy_bar(0:Nx-1,0:Ny-1))
      allocate(uz_rms(0:Nx-1,0:Ny-1))
      allocate(tdisp_bar(0:Nx-1,0:Ny-1))
      allocate(mdisp_bar(0:Nx-1,0:Ny-1))
      allocate(lz_bar(0:Nx-1,0:Ny-1))
      allocate(barotropic_ER(0:Nx-1,0:Ny-1))
      allocate(growth_rate(0:Nx-1,0:Ny-1))

      starts = (/ 1, 1, 1, 1 /)
      counts = (/ Nx, Ny, Nz, 1 /)
      ! Read the netcdf file, find the correct timestep.
      n = 0
      timestep = 0 
      write(55, "(3I16, 2F16.7)") Nx, Ny, nsteps, dx, dy
      do while (n.le.nsteps-1)
         n=n+1
         call check(nf90_get_var(ncid,timestep_varid,timestep,(/ n /)))
         print*,timestep,n
         !if(n.eq.1) firststep = timestep
        !end do
        !if(timestep.lt.nstep) then
        !laststep = timestep
              !print*,'Target step is not in this file range (',firststep,',',laststep,')'
        !endif
       
        ! Once found, read the time, and the ux, uy, uz vectors. 
        call check(nf90_get_var(ncid,t_varid,t,(/ n /) ) )   
        starts(4) = n
        if (.not. noux) then
           call check(nf90_get_var(ncid,ux_varid,ux,starts,counts) )
        endif
        if (.not. nouy) then
           call check(nf90_get_var(ncid,uy_varid,uy,starts,counts) )
        endif
        if (.not. nouz) then
           call check(nf90_get_var(ncid,uz_varid,uz,starts,counts) )
        endif
        if (.not. notemp) then
           call check(nf90_get_var(ncid,temp_varid,temp,starts,counts) )
        endif

        vortz_bar = 0.
        wz = 0.
        tdisp = 0.
        mdisp = 0.
        tot_energy = 0.
        ux_bar = 0.
        uy_bar = 0.
        uz_bar = 0.
        uz_rms = 0.
        tdisp_bar = 0.
        mdisp_bar = 0.
        lz_bar = 0.
        barotropic_ER = 0.
        growth_rate = 0.
      
        do k1=1,Nz-2
           do j1=1,Ny-2
              do i1=1,Nx-2
                  wz(i1,j1,k1) = (uy(i1+1,j1,k1)-uy(i1-1,j1,k1))/(2*dx) -(ux(i1,j1+1,k1)-ux(i1,j1-1,k1))/(2*dy)
                  tdisp(i1,j1,k1) = ((temp(i1+1,j1,k1)-temp(i1-1,j1,k1))/(2*dx))**2 +&
                                    ((temp(i1,j1+1,k1)-temp(i1,j1-1,k1))/(2*dy))**2 +&
                                    ((temp(i1,j1,k1+1)-temp(i1,j1,k1-1))/(2*dz))**2
                  tot_energy(i1,j1,k1) = ux(i1,j1,k1)**2 +&
                                         uy(i1,j1,k1)**2 +&
                                         uz(i1,j1,k1)**2
                  ! this is quite long.... sorry
                  ! |grad u|^2
                  mdisp(i1,j1,k1) = ((ux(i1+1,j1,k1)-ux(i1-1,j1,k1))/(2*dx))**2+&
                                    ((ux(i1,j1+1,k1)-ux(i1,j1-1,k1))/(2*dy))**2+&
                                    ((ux(i1,j1,k1+1)-ux(i1,j1,k1-1))/(2*dz))**2+&
                                    ((uy(i1+1,j1,k1)-uy(i1-1,j1,k1))/(2*dx))**2+&
                                    ((uy(i1,j1+1,k1)-uy(i1,j1-1,k1))/(2*dy))**2+&
                                    ((uy(i1,j1,k1+1)-uy(i1,j1,k1-1))/(2*dz))**2+&
                                    ((uz(i1+1,j1,k1)-uz(i1-1,j1,k1))/(2*dx))**2+&
                                    ((uz(i1,j1+1,k1)-uz(i1,j1-1,k1))/(2*dy))**2+&
                                    ((uz(i1,j1,k1+1)-uz(i1,j1,k1-1))/(2*dz))**2
              enddo
           enddo
        enddo

        !vertical averaging
        do k1 = 1, Nz-2
          vortz_bar = vortz_bar + wz(:,:,k1)
          uz_bar = uz_bar + uz(:,:,k1)
          uy_bar = uy_bar + uy(:,:,k1)
          ux_bar = ux_bar + ux(:,:,k1)
          uz_rms = uz_rms + uz(:,:,k1)**2
          tdisp_bar = tdisp_bar + tdisp(:,:,k1)
          mdisp_bar = mdisp_bar + mdisp(:,:,k1)
          barotropic_ER = barotropic_ER + tot_energy(:,:,k1)
        end do

        vortz_bar = vortz_bar/(Nz-2)
        uz_bar = uz_bar /(Nz-2)
        uy_bar = uy_bar /(Nz-2)
        ux_bar = ux_bar /(Nz-2)
        uz_rms = sqrt(uz_rms /(Nz-2))
        tdisp_bar = tdisp_bar/(Nz-2)
        mdisp_bar = mdisp_bar/(Nz-2)
        barotropic_ER = barotropic_ER/(Nz-2)
        barotropic_ER = (ux_bar**2 + uy_bar**2 + uz_bar**2)/barotropic_ER

        lz_bar = acf3D(uz, Nz, dz, Nx, Ny)
        
        !compute real part of the largest eigenvalue (Hattori & Hirota 2023)
        do j1 = 1, Ny-2
          do i1 = 1, Nx-2
            growth_rate(i1,j1) = -(ux_bar(i1+1,j1)-ux_bar(i1-1,j1))/(4*dx) - &
              (uy_bar(i1,j1+1)-uy_bar(i1,j1-1))/(4*dy) + real(sqrt(c1*(&
              ((ux_bar(i1+1,j1)-ux_bar(i1-1,j1))/(2*dx)&
              +(uy_bar(i1,j1+1)-uy_bar(i1,j1-1))/(2*dy))**2 - 4*(&
              -(invRo**2) + invRo*((ux_bar(i1,j1+1)-ux_bar(i1,j1-1))/(2*dy)&
                +(uy_bar(i1+1,j1)-uy_bar(i1-1,j1))/(2*dx)) + &
              (ux_bar(i1+1,j1)-ux_bar(i1-1,j1))/(2*dx)&
              *(uy_bar(i1,j1+1)-uy_bar(i1,j1-1))/(2*dy) + &
              -(uy_bar(i1+1,j1)-uy_bar(i1-1,j1))/(2*dx)*&
              (ux_bar(i1,j1+1)-ux_bar(i1,j1-1))/(2*dy))&
              )))/2.
          end do 
        end do 

        do j1 = 1, Ny-2
          do i1 = 1, Nx-2
            write(55,"(2I16, 11E16.7)") i1, j1, t, ux_bar(i1,j1),&
              uy_bar(i1,j1), uz_bar(i1,j1), vortz_bar(i1,j1), tdisp_bar(i1,j1),&
              uz_rms(i1,j1), lz_bar(i1,j1), barotropic_ER(i1,j1),&
              mdisp_bar(i1,j1), growth_rate(i1,j1)
          end do 
        end do 
        write(55,*)
        write(55,*)
      enddo
      deallocate(ux,uy,uz)
      deallocate(wz,tdisp,temp,tot_energy,mdisp)
      deallocate(ux_bar,uy_bar,uz_bar,tdisp_bar,uz_rms,mdisp_bar,growth_rate)
      deallocate(vortz_bar,lz_bar,barotropic_ER)
    end do
  else ! if only one simdat file is given
    CALL GET_COMMAND_ARGUMENT(1,netcdf_file)
    write(*,'(2a)') "Processing ",trim(netcdf_file)
    write(*,'(a,3I5)') "Looking for timestep",nstep

    call check( nf90_open(netcdf_file, nf90_nowrite, ncid) )
    call check( nf90_inq_varid(ncid,'t',t_varid) ) 
    call check( nf90_inq_varid(ncid,'timestep',timestep_varid) ) 
    ! Looks to see if the field exists, modifies flag if not.  
    call check_for_field(ux_varid,'ux',noux)
    call check_for_field(uy_varid,'uy',nouy)
    call check_for_field(uz_varid,'uz',nouz)
    call check_for_field(temp_varid,'Temp',notemp)
    ! Reads in X, Y, and TIME, as well as their dimensions.  
    call check(nf90_inq_dimid(ncid,"X",x_dimid) )
    call check(nf90_inquire_dimension(ncid,x_dimid,len=Nx) )
    call check(nf90_inq_dimid(ncid,"Y",y_dimid) )
    call check(nf90_inquire_dimension(ncid,y_dimid,len=Ny) )
    call check(nf90_inq_dimid(ncid,"Z",z_dimid) )
    call check(nf90_inquire_dimension(ncid,z_dimid,len=Nz) )
    call check(nf90_inq_dimid(ncid,"TIME",time_dimid) )
    call check(nf90_inquire_dimension(ncid,time_dimid,len=nsteps) )

    Write(*,'(a,i4,a,i4,a,i4)') "Nx =",Nx," Ny =",Ny," Nz =",Nz
    Write(*,'(a,i7)') "Number of saved timesteps =",nsteps
    ! Reads in the input parameters 
    call check(nf90_inq_varid(ncid,"B_therm", B_therm_varid) )
    call check(nf90_get_var(ncid,B_therm_varid,B_therm) )
    call check(nf90_inq_varid(ncid,"B_comp", B_comp_varid) )
    call check(nf90_get_var(ncid,B_comp_varid,B_comp) )
    call check(nf90_inq_varid(ncid,"D_visc", D_visc_varid) )
    call check(nf90_get_var(ncid,D_visc_varid,D_visc) )
    call check(nf90_inq_varid(ncid,"D_therm", D_therm_varid) )
    call check(nf90_get_var(ncid,D_therm_varid,D_therm) )
    call check(nf90_inq_varid(ncid,"D_comp", D_comp_varid) )
    call check(nf90_get_var(ncid,D_comp_varid,D_comp) )
    call check(nf90_inq_varid(ncid,"S_therm", S_therm_varid) )
    call check(nf90_get_var(ncid,S_therm_varid,S_therm) )
    call check(nf90_inq_varid(ncid,"S_comp", S_comp_varid) )
    call check(nf90_get_var(ncid,S_comp_varid,S_comp) )
    call check(nf90_inq_varid(ncid,"Gammax", Gammax_varid) )
    call check(nf90_get_var(ncid,Gammax_varid,Gammax) )
    call check(nf90_inq_varid(ncid,"Gammay", Gammay_varid) )
    call check(nf90_get_var(ncid,Gammay_varid,Gammay) )
    call check(nf90_inq_varid(ncid,"Gammaz", Gammaz_varid) )
    call check(nf90_get_var(ncid,Gammaz_varid,Gammaz) )
    print*,"halli"


    dx = Gammax/Nx
    dy = Gammay/Ny
    dz = Gammaz/Nz
    ! Allocate arrays for ux, uy, uz
    allocate(ux(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(uy(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(uz(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(wz(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(temp(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(tdisp(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(mdisp(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(tot_energy(0:Nx-1,0:Ny-1,0:Nz-1))

    allocate(vortz_bar(0:Nx-1,0:Ny-1))
    allocate(uz_bar(0:Nx-1,0:Ny-1))
    allocate(ux_bar(0:Nx-1,0:Ny-1))
    allocate(uy_bar(0:Nx-1,0:Ny-1))
    allocate(uz_rms(0:Nx-1,0:Ny-1))
    allocate(tdisp_bar(0:Nx-1,0:Ny-1))
    allocate(mdisp_bar(0:Nx-1,0:Ny-1))
    allocate(lz_bar(0:Nx-1,0:Ny-1))
    allocate(barotropic_ER(0:Nx-1,0:Ny-1))
    allocate(growth_rate(0:Nx-1,0:Ny-1))

    starts = (/ 1, 1, 1, 1 /)
    counts = (/ Nx, Ny, Nz, 1 /)
    ! Read the netcdf file, find the correct timestep.
    n = 0
    timestep = 0 
    write(55, "(3I16, 2F16.7)") Nx, Ny, nsteps, dx, dy
    do while (n.le.nsteps-1)
       n=n+1
       call check(nf90_get_var(ncid,timestep_varid,timestep,(/ n /)))
       print*,timestep,n
       !if(n.eq.1) firststep = timestep
      !end do
      !if(timestep.lt.nstep) then
      !laststep = timestep
            !print*,'Target step is not in this file range (',firststep,',',laststep,')'
      !endif
     
      ! Once found, read the time, and the ux, uy, uz vectors. 
      call check(nf90_get_var(ncid,t_varid,t,(/ n /) ) )   
      starts(4) = n
      if (.not. noux) then
         call check(nf90_get_var(ncid,ux_varid,ux,starts,counts) )
      endif
      if (.not. nouy) then
         call check(nf90_get_var(ncid,uy_varid,uy,starts,counts) )
      endif
      if (.not. nouz) then
         call check(nf90_get_var(ncid,uz_varid,uz,starts,counts) )
      endif
      if (.not. notemp) then
         call check(nf90_get_var(ncid,temp_varid,temp,starts,counts) )
      endif

      vortz_bar = 0.
      wz = 0.
      tdisp = 0.
      mdisp = 0.
      tot_energy = 0.
      ux_bar = 0.
      uy_bar = 0.
      uz_bar = 0.
      uz_rms = 0.
      tdisp_bar = 0.
      mdisp_bar = 0.
      lz_bar = 0.
      barotropic_ER = 0.
    
      !Calculat/ the horizontal average of uw and of u for each z, wrte it to file
      do k1=1,Nz-2
         do j1=1,Ny-2
            do i1=1,Nx-2
                wz(i1,j1,k1) = (uy(i1+1,j1,k1)-uy(i1-1,j1,k1))/(2*dx) -(ux(i1,j1+1,k1)-ux(i1,j1-1,k1))/(2*dy)
                tdisp(i1,j1,k1) = ((temp(i1+1,j1,k1)-temp(i1-1,j1,k1))/(2*dx))**2 +&
                                  ((temp(i1,j1+1,k1)-temp(i1,j1-1,k1))/(2*dy))**2 +&
                                  ((temp(i1,j1,k1+1)-temp(i1,j1,k1-1))/(2*dz))**2
                tot_energy(i1,j1,k1) = ux(i1,j1,k1)**2 +&
                                       uy(i1,j1,k1)**2 +&
                                       uz(i1,j1,k1)**2
                mdisp(i1,j1,k1) = ((ux(i1+1,j1,k1)-ux(i1-1,j1,k1))/(2*dx))**2+&
                                  ((ux(i1,j1+1,k1)-ux(i1,j1-1,k1))/(2*dy))**2+&
                                  ((ux(i1,j1,k1+1)-ux(i1,j1,k1-1))/(2*dz))**2+&
                                  ((uy(i1+1,j1,k1)-uy(i1-1,j1,k1))/(2*dx))**2+&
                                  ((uy(i1,j1+1,k1)-uy(i1,j1-1,k1))/(2*dy))**2+&
                                  ((uy(i1,j1,k1+1)-uy(i1,j1,k1-1))/(2*dz))**2+&
                                  ((uz(i1+1,j1,k1)-uz(i1-1,j1,k1))/(2*dx))**2+&
                                  ((uz(i1,j1+1,k1)-uz(i1,j1-1,k1))/(2*dy))**2+&
                                  ((uz(i1,j1,k1+1)-uz(i1,j1,k1-1))/(2*dz))**2
            enddo
         enddo
      enddo

      !vertical averaging
      do k1 = 1, Nz-2
        vortz_bar = vortz_bar + wz(:,:,k1)
        uz_bar = uz_bar + uz(:,:,k1)
        uy_bar = uy_bar + uy(:,:,k1)
        ux_bar = ux_bar + ux(:,:,k1)
        uz_rms = uz_rms + uz(:,:,k1)**2
        tdisp_bar = tdisp_bar + tdisp(:,:,k1)
        mdisp_bar = mdisp_bar + mdisp(:,:,k1)
        barotropic_ER = barotropic_ER + tot_energy(:,:,k1)
      end do

      vortz_bar = vortz_bar/(Nz-2)
      uz_bar = uz_bar /(Nz-2)
      uy_bar = uy_bar /(Nz-2)
      ux_bar = ux_bar /(Nz-2)
      uz_rms = sqrt(uz_rms /(Nz-2))
      tdisp_bar = tdisp_bar/(Nz-2)
      mdisp_bar = mdisp_bar/(Nz-2)
      barotropic_ER = barotropic_ER/(Nz-2)
      barotropic_ER = (ux_bar**2 + uy_bar**2 + uz_bar**2)/barotropic_ER

      lz_bar = acf3D(uz, Nz, dz, Nx, Ny)

      do j1 = 1, Ny-2
        do i1 = 1, Nx-2
          growth_rate(i1,j1) = -(ux_bar(i1+1,j1)-ux_bar(i1-1,j1))/(4*dx) - &
              (uy_bar(i1,j1+1)-uy_bar(i1,j1-1))/(4*dy) + real(sqrt(c1*(&
              ((ux_bar(i1+1,j1)-ux_bar(i1-1,j1))/(2*dx)&
              +(uy_bar(i1,j1+1)-uy_bar(i1,j1-1))/(2*dy))**2 - 4*(&
              -(invRo**2) + invRo*((ux_bar(i1,j1+1)-ux_bar(i1,j1-1))/(2*dy)&
                +(uy_bar(i1+1,j1)-uy_bar(i1-1,j1))/(2*dx)) + &
              (ux_bar(i1+1,j1)-ux_bar(i1-1,j1))/(2*dx)&
              *(uy_bar(i1,j1+1)-uy_bar(i1,j1-1))/(2*dy) + &
              -(uy_bar(i1+1,j1)-uy_bar(i1-1,j1))/(2*dx)*&
              (ux_bar(i1,j1+1)-ux_bar(i1,j1-1))/(2*dy)))))/2.
        end do 
      end do 

      do j1 = 1, Ny-2
        do i1 = 1, Nx-2
          write(55,"(2I16, 12E16.7)") i1, j1, t, ux_bar(i1,j1),&
            uy_bar(i1,j1), uz_bar(i1,j1), vortz_bar(i1,j1), tdisp_bar(i1,j1),&
            uz_rms(i1,j1), lz_bar(i1,j1), barotropic_ER(i1,j1),&
            mdisp_bar(i1,j1), growth_rate(i1,j1)
        end do 
      end do 
      write(55,*)
      write(55,*)
    enddo
    deallocate(ux,uy,uz)
    deallocate(wz,tdisp,temp,tot_energy,mdisp)
    deallocate(ux_bar,uy_bar,uz_bar,tdisp_bar,uz_rms,mdisp_bar)
    deallocate(vortz_bar,lz_bar,barotropic_ER,growth_rate)
  end if
  
  close(55)
  close(20)
  contains
    subroutine check_for_field(varid,name,nofield)
      implicit none
      integer :: varid
      character (len=*) :: name
      logical :: nofield
      ierr=nf90_inq_varid(ncid,name,varid) 
      if (ierr == nf90_noerr) then 
         nfields = nfields+1
         nofield = .false.
      else
         if (ierr.eq.nf90_enotvar) then 
            nofield = .true.
            varid = 0 
            write(*,'(3a)') 'No ',name,' found!'
         else
            call check(ierr)
         end if
      endif
    end subroutine check_for_field

    subroutine check(status)
      integer, intent ( in) :: status
      
      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if
    end subroutine check
    
    function acf3D(field, Np, dn, N1, N2) result(val)
      implicit none
      real :: field(:,:,:)
      real :: dn
      integer :: Np, N1, N2, sk1, si1, sj1, l
      real, dimension(N1, N2) :: num, dem, val
      real, dimension(int(Np/2),N1,N2) :: acf
  
      do l = 1, int(Np/2)
        num = 0.
        dem = 0.
        do sk1 = 0, Np-l-1 
           num = num + field(:,:,sk1)*field(:,:,sk1+l)
           dem = dem + field(:,:,sk1)**2
        end do 
        acf(l, :,:) = num/dem
      end do 
      do sj1 = 1, N2
        do si1 = 1, N1
          do sk1 = 1, int(Np/2)
            if(acf(sk1,si1,sj1) .lt. 0.5) then
              val(si1,sj1) = sk1*dn +&
                dn*(0.5-acf(sk1,si1,sj1))/(acf(sk1,si1,sj1)-acf(sk1-1,si1,sj1))
              exit
            end if
          end do 
        end do 
      end do 
    end function acf3D
end program netcdfextract
