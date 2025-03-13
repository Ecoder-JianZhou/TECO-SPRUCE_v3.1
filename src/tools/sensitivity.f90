module sensitivity
    use datatypes
    use driver
    use mcmc_mod
    implicit none
    integer :: npar4DA
    
contains

    subroutine init_sa(in_mcmc_configfile, st)
        implicit none
        ! *** in_mcmc_configfile : array of files for vegetation parameters
        ! ***              st    : return the st data
        ! *** put values to the variables for MCMC
        character(*), intent(in) :: in_mcmc_configfile
        type(site_data_type), intent(inout) :: st
        real(8), allocatable :: temp_parmin(:), temp_parmax(:), temp_parval(:)
        integer, allocatable :: temp_idx_st(:), temp_idx_sp(:,:)
        integer :: ipft, ipar_st, ipar_sp
        integer, allocatable :: mark_npar(:)

        ! read the nml file of MCMC configs (eg. TECO_MCMC_configs.nml)
        call readMCMC_configs_NML(in_mcmc_configfile)    ! update mc_params: st parameters and sp parameters
        ! handle the parameters for MCMC
        allocate(temp_parmin(npar), temp_parmax(npar), temp_parval(npar))  ! allocate the temporary parmin value
        ! allocate(temp_paridx(npar))  ! mark the index of parameters for MCMC

        allocate(temp_idx_st(nsite_params))
        allocate(temp_idx_sp(npft, nspec_params))
        npar4DA = 0        ! How many parameters need to be optimized
        npar4st = 0
        allocate(mark_npar(1+npft))
        do ipar_st = 1, nsite_params
            if(mc_params%st(ipar_st)%parmax .ne.  mc_params%st(ipar_st)%parmin)then
                npar4DA = npar4DA + 1
                npar4st = npar4st + 1
                temp_parmin(npar4DA) = mc_params%st(ipar_st)%parmin
                temp_parmax(npar4DA) = mc_params%st(ipar_st)%parmax
                temp_parval(npar4DA) = mc_params%st(ipar_st)%parval
                temp_idx_st(npar4st) = ipar_st
                print*,npar4DA,mc_params%st(ipar_st)%parname, mc_params%st(ipar_st)%parval, &
                    mc_params%st(ipar_st)%parmin, mc_params%st(ipar_st)%parmax
            endif
            ! print*,ipar_st, mc_params%st(ipar_st)%parname, mc_params%st(ipar_st)%parval, &
            ! mc_params%st(ipar_st)%parmin, mc_params%st(ipar_st)%parmax
        enddo
        mark_npar(1) = npar4st
        do ipft = 1, npft
            npar4sp = 0
            do ipar_sp = 1, nspec_params
                if(mc_params%sp(ipft)%sp_params(ipar_sp)%parmax .gt.  mc_params%sp(ipft)%sp_params(ipar_sp)%parmin)then
                    npar4DA = npar4DA + 1
                    npar4sp = npar4sp + 1
                    temp_parmin(npar4DA) = mc_params%sp(ipft)%sp_params(ipar_sp)%parmin
                    temp_parmax(npar4DA) = mc_params%sp(ipft)%sp_params(ipar_sp)%parmax
                    temp_parval(npar4DA) = mc_params%sp(ipft)%sp_params(ipar_sp)%parval
                    temp_idx_sp(ipft, npar4sp) = ipar_sp
                    print*,npar4DA,mc_params%sp(ipft)%sp_params(ipar_sp)%parname, mc_params%sp(ipft)%sp_params(ipar_sp)%parval, &
                        mc_params%sp(ipft)%sp_params(ipar_sp)%parmin, mc_params%sp(ipft)%sp_params(ipar_sp)%parmax
                endif
                ! print*,mc_params%sp(ipft)%sp_params(ipar_sp)%parname, mc_params%sp(ipft)%sp_params(ipar_sp)%parval, &
                !         mc_params%sp(ipft)%sp_params(ipar_sp)%parmin, mc_params%sp(ipft)%sp_params(ipar_sp)%parmax
            enddo
            mark_npar(1+ipft) = npar4sp
        enddo

        allocate(mc_DApar%DAparmin(npar4DA),    mc_DApar%DAparmax(npar4DA), &
                 mc_DApar%DApar(npar4DA),       mc_DApar%DApar_old(npar4DA), &
                 mc_DApar%DAparidx_st(npar4st), mc_DApar%DAparidx_sp(npft, npar4sp))

        mc_DApar%DAparmin    = temp_parmin(:npar4DA)
        mc_DApar%DAparmax    = temp_parmax(:npar4DA)
        mc_DApar%DApar       = temp_parval(:npar4DA)
        mc_DApar%DApar_old   = mc_DApar%DApar         ! mark as old parameters
        mc_DApar%DAparidx_st = temp_idx_st
        mc_DApar%DAparidx_sp = temp_idx_sp

        deallocate(temp_parmin, temp_parmax, temp_parval, temp_idx_st, temp_idx_sp)
    end subroutine init_sa

    subroutine run_sa(st)
        implicit none
        type(site_data_type), intent(inout) :: st
        integer :: ipar, ipar4sa
        character(10) :: str_i_simu, str_i_par

        do_out_hr   = .True.
        do_out_day  = .False.
        do_out_mon  = .True.
        do_out_yr   = .False.
        do ipar = 1, npar4DA
            ! select 100 parameters range
            write(str_i_par, "(I0.3)") ipar
            do ipar4sa = 1, 100
                print*, ipar, ipar4sa
                mc_DApar%DApar(ipar) = mc_DApar%DAparmin(ipar) + ipar4sa*(mc_DApar%DAparmax(ipar) - mc_DApar%DAparmin(ipar))/100
                write(str_i_simu, "(I0.3)") ipar4sa
                mc_str_n = adjustl(trim(str_i_par))//"-"//adjustl(trim(str_i_simu))
                call mc_update_mc_params()
                call mc_update_params4simu()
                call initialize_teco(st)!, .True.)
                call teco_simu(st, .True.)            ! run the model
            enddo
        enddo
    end subroutine run_sa
    
end module sensitivity