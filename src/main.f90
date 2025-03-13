program TECO
    use datatypes
    use driver
    use io_mod
    use mcmc
    use mcmc_mod
    use mod_spinup
    use sensitivity
    
    implicit none
    integer :: count_mode 
    integer num_args, ierr
    type(site_data_type) :: st
    ! logical :: do_sensitive

    ! do_sensitive = .False.
    do_sa = .False.
    print *, ""
    write(*,*) "# -----------------------------------------"
    ! get the count of command
    num_args = COMMAND_ARGUMENT_COUNT()
    ! check if have the file
    if (num_args /= 1) then
        write(*,*) "Usage: ./my_program <conf_nml_file>"
        stop
    end if
    ! read the command
    call GET_COMMAND_ARGUMENT(1, conf_nml_file, ierr)
    conf_nml_file = adjustl(trim("configs/"))//adjustl(trim(conf_nml_file))
    write(*,*) "Reading the file of ", conf_nml_file

    print *, ""
    write(*,*) "# -----------------------------------------"
    call read_configs_nml(conf_nml_file)
    ! check the three mode: do_simu; do_mcmc; do_spinup
    count_mode = 0
    if (do_simu)   count_mode = count_mode + 1
    if (do_mcmc)   count_mode = count_mode + 1
    if (do_spinup) count_mode = count_mode + 1

    if (count_mode .eq. 0) then
        print *, "# You must choose a run mode."
        print *, "# Please make sure one of the three modes is True in file of TECO_model_configs.nml"
        print *, "#    *do_simu; *do_mcmc; *do_spinup"
        print *, ""
        if (do_sa) then
            print *, "# Running the sensitive analysis ... "
            continue
        else
            stop
        endif
    elseif (count_mode .gt. 1) then
        print *, "# You can only select one mode out of the three, please check your file."
        print *, "# Please check the file of TECO_model_configs.nml"
        print *, "#    *do_simu; *do_mcmc; *do_spinup"
        print *, ""
        stop
    else
        continue
    endif

    write(*,*) "# Start to run the case of """, adjustl(trim(case_name)), """"
    ! update the in-out path and create the revelent ouput paths
    call createNewCase() 

    call get_forcingdata()                      ! read forcing data
    ! nHours  = nforcing                          
    ! nDays   = int(nHours/24.)
    ! nYears  = int(nforcing/(365*24))
    ! nMonths = nYears*12

    ! print *, allocated(outVars_y%gpp), nYears
    if (.not. do_snow) call get_snowdepth()
    if (do_restart)then
        ! call read_restart(restartfile)     ! this module in "writeOutput2nc.f90"
        ! call initialize_with_restart()
    endif

    ! read the initilized values and parameters
    params_nml_file = adjustl(trim("configs/"))//adjustl(trim(params_nml_file))
    call read_nml_params_initValues(params_nml_file)
    allocate(st%sp(npft))
    ! initialize the output
    allocate(outVars_h%sp(npft))
    allocate(outVars_d%sp(npft))
    allocate(outVars_m%sp(npft))
    allocate(outVars_y%sp(npft))

    if(do_simu)then
        print *, "# Start to run simulation mode."
        ! assign the output variables
        ! if(do_out_hr)  call assign_outVars(tot_outVars_h, nHours,  npft)
        ! if(do_out_day) call assign_outVars(tot_outVars_d, nDays,   npft)
        ! if(do_out_mon) call assign_outVars(tot_outVars_m, nMonths, npft)
        ! if(do_out_yr)  call assign_outVars(tot_outVars_y, nYears,  npft)
        ! initilize the model 
        ! call initilize(file_site_params, files_pft_params, st)
        call initialize_teco(st)!, .False.)

        ! Start to run TECO model
        call teco_simu(st, .True.)            ! run simulation
        ! write the output data
! #ifdef USE_NETCDF
!         if(do_out_hr)  call write_outputs_nc(outDir_h, outVars_h, nHours, "hourly") 
!         if(do_out_day) call write_outputs_nc(outDir_d, outVars_d, nDays,  "daily") 
!         if(do_out_mon) call write_outputs_nc(outDir_m, outVars_m, nMonths,"monthly") 
! #endif
        ! print*, outvars_d%sp(1)%gpp
        ! if(do_out_hr)  call write_outputs_csv(outDir_csv, outVars_h, nHours, "hourly") 
        ! if(do_out_day) call write_outputs_csv(outDir_csv, outVars_d, nDays,  "daily") 
        ! if(do_out_mon) call write_outputs_csv(outDir_csv, outVars_m, nMonths,"monthly") 

        if(do_out_hr)  call deallocate_results(outVars_h, npft)
        if(do_out_day) call deallocate_results(outVars_d, npft)
        if(do_out_mon) call deallocate_results(outVars_m, npft)
        if(do_out_yr)  call deallocate_results(outVars_y, npft)
        
    elseif(do_spinup)then
        print*, "# Start to run MCMC mode."
        call initialize_teco(st)!, .False.)
        call init_spinup_variables()    ! initilize the spin-up variables
        call run_spinup(st)               ! run spin-up loops
        call write_spinup_res()         ! write the results of SPIN-UP
        call write_restart()            ! write the result file
        call deallo_spinup_variables()  ! deallocate the variables of SPIN-UP
    elseif(do_mcmc) then
        print*, "# Start to run MCMC mode."
        call initialize_teco(st)!, .False.)
        call init_mcmc(mcmc_configfile, st)  ! initilize the MCMC 
        call run_mcmc(st)                 ! run MCMC
        call deallocate_mcmc()          ! deallocate the MCMC variables 
    elseif(do_sa) then
        print*, "# Start to run sensitivity analysis ..."
        call initialize_teco(st)!, .False.)
        call init_sa(mcmc_configfile, st)  ! initilize the MCMC used to do sensitivity analysis
        call run_sa(st)
        ! call deallocate_sa()
        call deallocate_mcmc()
    endif
    
    if(allocated(outVars_h%sp)) deallocate(outVars_h%sp)
    if(allocated(outVars_d%sp)) deallocate(outVars_d%sp)
    if(allocated(outVars_m%sp)) deallocate(outVars_m%sp)
    if(allocated(outVars_y%sp)) deallocate(outVars_y%sp)
    if(allocated(sp_names))     deallocate(sp_names)
    if(allocated(st%sp))        deallocate(st%sp)

    call deallocate_date_type()
end program TECO


subroutine createNewCase()
    use datatypes
    ! use mcmc_functions
    ! use mod_mcmc
    implicit none
    ! create a new case to run the TECO model
    !   * create the output path

    print *, "# Update and create the output dirs"

    ! update the full path of input file
    climfile        = adjustl(trim(inDir))//"/"//adjustl(trim(climfile))       ! climate file name
    snowdepthfile   = adjustl(trim(inDir))//"/"//adjustl(trim(snowdepthfile))  ! snow depthfile
    restartfile     = adjustl(trim(inDir))//"/"//adjustl(trim(restartfile))    ! restartfile
    watertablefile  = adjustl(trim(inDir))//"/"//adjustl(trim(watertablefile)) ! Jian: maybe used when not run soil_physical
    ! check the inputfile
    call check_inputfile(climfile, "climate file")
    if (.not. do_snow)    call check_inputfile(snowdepthfile,  "The file of snowdepth")
    if (do_restart)       call check_inputfile(restartfile,    "The file of restart")
    if (.not. do_soilphy) call check_inputfile(watertablefile, "The file of water table")
    ! create the outdir
    call CreateFolder(adjustl(trim(outdir)))

    ! update and create the output dir of case
    outdir_case = adjustl(trim(outdir))//"/"//adjustl(trim(case_name))
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
    outdir_case = adjustl(trim(outdir))//"\"//adjustl(trim(case_name))
#endif
    call CreateFolder(adjustl(trim(outdir_case)))

    ! if (index(achar(1), 'W') == 1) then
    !     print *, "This is Windows platform"
    ! else
    !     print *, "This is not Windows platform"
    ! end if

    ! Check if the platform is Windows
! #if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
!   WRITE(*, *) 'Running on Windows.'
! #else
!   WRITE(*, *) 'Running on a non-Windows platform.'
! #endif

    ! update and create the output dir of each format outputs
#ifdef USE_NETCDF
    outDir_nc  = adjustl(trim(outdir_case))//"/"//adjustl(trim(outDir_nc))
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
    outDir_nc  = adjustl(trim(outdir_case))//"\"//adjustl(trim(outDir_nc))
#endif
    call CreateFolder(adjustl(trim(outDir_nc)))
#endif
    outDir_csv = adjustl(trim(outdir_case))//"/"//adjustl(trim(outDir_csv))
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
    outDir_csv = adjustl(trim(outdir_case))//"\"//adjustl(trim(outDir_csv))
#endif
    call CreateFolder(adjustl(trim(outDir_csv)))
    ! print*, outDir_csv

    ! update and create the output for each time frequency of nc-format outputs
    if (do_out_hr) then
#ifdef USE_NETCDF
        outDir_h = adjustl(trim(outdir_nc))//"/Hourly"
        call CreateFolder(adjustl(trim(outDir_h)))
#endif
    endif
    if (do_out_day) then
        outDir_d = adjustl(trim(outdir_nc))//"/Daily"
        call CreateFolder(adjustl(trim(outDir_d)))
    endif
    if (do_out_mon) then
        outDir_m = adjustl(trim(outdir_nc))//"/Monthly"
        call CreateFolder(adjustl(trim(outDir_m)))
    endif
    if (do_out_yr) then
        outDir_y = adjustl(trim(outdir_nc))//"/Yearly"
        call CreateFolder(adjustl(trim(outDir_y)))
    endif

    if (do_spinup)then
        outDir_spinup = adjustl(trim(outdir_case))//"/"//adjustl(trim(outDir_spinup))
        call CreateFolder(adjustl(trim(outDir_spinup)))
        outfile_restart = adjustl(trim(outDir_spinup))//"/restart.nc"
    endif

    if (do_mcmc)then
        outDir_mcmc = adjustl(trim(outdir_case))//"/"//adjustl(trim(outDir_mcmc))
        call CreateFolder(adjustl(trim(outDir_mcmc)))
        ! if (do_mc_out_hr)then
            outDir_mcmc_h = adjustl(trim(outDir_mcmc))//"/results_mcmc_Hourly"
            call CreateFolder(adjustl(trim(outDir_mcmc_h)))
        ! endif
        ! if (do_mc_out_day) then
            outDir_mcmc_d = adjustl(trim(outDir_mcmc))//"/results_mcmc_Daily"
            call CreateFolder(adjustl(trim(outDir_mcmc_d)))
        ! endif
        ! if (do_mc_out_mon) then
            outDir_mcmc_m = adjustl(trim(outDir_mcmc))//"/results_mcmc_Monthly"
            call CreateFolder(adjustl(trim(outDir_mcmc_m)))
        ! endif
    endif

    if(do_restart)then
        restartfile = adjustl(trim(inDir))//adjustl(trim(restartfile))
    endif
end subroutine createNewCase

subroutine CreateFolder(path_new)
    implicit none
    character(len=*), INTENT(in) :: path_new
    character (len=:), allocatable :: cmdChar
    logical :: dirExists
    ! ----------------------------------------------------
    allocate(character(len=6+len(path_new)) :: cmdChar)
    cmdChar = "mkdir "//path_new
    inquire( file=trim(path_new)//'\.', exist=dirExists )  ! Works with gfortran, but not ifort
    ! inquire( directory=newDirPath, exist=dirExists )         ! Works with ifort, but not gfortran
    print *, cmdChar
    if (.not. dirExists) call system(cmdChar)
    deallocate(cmdChar)
end subroutine CreateFolder

subroutine check_inputfile(filepath, whatfile)
    implicit none
    character(*), intent(in) :: filepath, whatfile
    logical :: file_exists

    inquire(file=filepath, exist=file_exists)
    
    if (file_exists) then
        print *, "# ", whatfile," exists: "
        print *, "#     ", filepath
    else
        print *, "# ", whatfile," does not exist: "
        print *, "#     ", filepath
        stop
    end if
end subroutine check_inputfile

function replace_slash(s) result(r)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: r
    integer :: i
    do i = 1, len(s)
      if (s(i:i) == "/") then
        r(i:i) = "\"
      else
        r(i:i) = s(i:i)
      end if
    end do
  end function replace_slash
