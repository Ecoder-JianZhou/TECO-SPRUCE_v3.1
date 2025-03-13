module mcmc
    use driver
    use datatypes
    use mcmc_mod
    use MCMC_outputs

    implicit none

    integer ipar, covexist, npar4DA

    real(8) :: fact_rejet
    real(8) J_last(20), J_new(20), accept_rate, J_show_old, J_show_new, delta_scale, delta_scale_min, delta_scale_max
    integer new, reject
    logical do_cov2createNewPars, do_cov
    integer, allocatable :: mark_npar(:)
    integer :: before_upgrade, after_upgrade

    contains
    
    subroutine init_mcmc(in_mcmc_configfile, st)
    ! *** in_mcmc_configfile : array of files for vegetation parameters
    ! ***              st    : return the st data
    ! *** put values to the variables for MCMC
        implicit none
        character(*), intent(in) :: in_mcmc_configfile
        type(site_data_type), intent(inout) :: st
        real(8), allocatable :: temp_parmin(:), temp_parmax(:), temp_parval(:)
        integer, allocatable :: temp_idx_st(:), temp_idx_sp(:,:)
        integer :: ipft, ipar_st, ipar_sp

        ! read the nml file of MCMC configs (eg. TECO_MCMC_configs.nml)
        call readMCMC_configs_NML(in_mcmc_configfile)    ! update mc_params: st parameters and sp parameters
        ! read the observational data
        call readObsData() ! return a type array of vars4MCMC

        ! initilize the parameters and initial values in TECO model
        ! allocate(npar4DA(npft))         ! How many parameters need to be optimized
        ! allocate(mc_parvals(npft))      ! set: parval, parmin, parmax
        ! allocate(mc_in_params(npft))    ! 
        ! allocate(mc_init_params(npft))

        ! if(allocated(vegn%allSp)) then
        !     do ipft = 1, npft
        !         allocate(mc_parvals(ipft)%parval(npar), mc_parvals(ipft)%parmin(npar), mc_parvals(ipft)%parmax(npar))
        !         call readParamNml(adjustl(trim("configs/"//adjustl(trim(files_vegn_params(ipft))))), &
        !             mc_in_params(ipft), mc_init_params(ipft), &
        !             mc_parvals(ipft)%parval, mc_parvals(ipft)%parmin, mc_parvals(ipft)%parmax)
        !         call initilize_site(mc_in_params(ipft), mc_init_params(ipft))  ! Jian: this version not separate the site parameters and pft parameters
        !         call initilize_spec(vegn%allSp(ipft), mc_in_params(ipft), mc_init_params(ipft))
        !     enddo
        ! endif ! finish the initilize parameters and the mc_parvals

        
        ! handle the parameters for MCMC
        ! allocate(mc_DApar(npft)) 
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

        ! give some values to the parameters for MCMC
        covexist      = 0
        fact_rejet    = 2.4/sqrt(real(npar4DA))

        ! record
        allocate(mc_DApar%coefhistory(ncov, npar4DA))
        ! create the coefnorm for generating the new parameters
        allocate(mc_DApar%coefnorm(npar4DA)) 
        allocate(mc_DApar%coefac(npar4DA))
        do ipar = 1, npar4DA
            mc_DApar%coefnorm(ipar) = 0.5
            mc_DApar%coefac(ipar)   = mc_DApar%coefnorm(ipar)
        enddo
        allocate(mc_DApar%gamnew(npar4DA, npar4DA))
        J_last = 9000000.0
        ! init the outputs
        call init_mcmc_outputs(npar4DA)
        do_cov = .True.
        do_cov2createNewPars = .False.
    end subroutine init_mcmc

    subroutine run_mcmc(st)
        implicit none
        type(site_data_type), intent(inout) :: st
        integer temp_upgraded, ipft, mark4scale, ishow, nonaccept
        real(8) rand, init_scale, rand_scale
        
        print *, "# Start to run mcmc ..."
        call generate_newPar()
        ! print *, "here 1..."
        ! call generate_rand_newpar()
        upgraded = 0
        new = 0
        mark4scale = 0
        init_scale = search_scale
        nonaccept  = 0
        delta_scale = 0.1
        delta_scale_min = 0.01
        delta_scale_max = 1
        ! print *, "here 2..."
        do iDAsimu = 1, nDAsimu
            ! write(*,*) iDAsimu, "/", nDAsimu, J_last, J_new, upgraded, accept_rate
            ! print *, "here 3..."
            write(*,*) iDAsimu, "/", nDAsimu,  J_show_old, J_show_new, upgraded, accept_rate
            ! do ishow = 1, 20
            !     write(*,*) iDAsimu, "/", nDAsimu, upgraded, ishow, J_last(ishow), J_new(ishow), J_last(ishow) - J_new(ishow)
            ! enddo
    !         write(*,*) iDAsimu, "/", nDAsimu, J_last(1),"/", J_new(1),";", J_last(2),"/", J_new(2),";",&
    ! & J_last(3),"/", J_new(3),";", J_last(4),"/", J_new(4),";", J_last(5),"/", J_new(5),";",upgraded, accept_rate
            call mcmc_functions_init()  ! initialize the mc_itime ... variables
                
            ! print*,"before:", mc_DApar%DApar
            ! generate parameters 
            call generate_newPar()
            ! print*,"after: ",  mc_DApar%DApar  
            ! update the parameters
            call mc_update_mc_params()
            call mc_update_params4simu()
            
            ! call renewMDpars(mc_parvals(ipft)%parval, mc_in_params(ipft))          ! call update parameters in TECO model
            
            ! call initialize()           ! initialize the TECO model 
            call initialize_teco(st)!, .True.)
            
            ! finish ! initialize the TECO model
            call teco_simu(st, .False.)            ! run the model
            if(ISNAN(st%Rsoilab3)) then
                cycle
            endif 
            ! print*, "test: ", st%sp(3)%GLmax, st%sp(3)%GRmax, st%sp(3)%Gsmax
            temp_upgraded = upgraded
            ! call costFuncObs()          ! calculate the cost between observations and simulations
            call costFuncObs_old()
            ! call costFuncObs_update()
            ! call costFuncObs_new()

            ! if upgraded is updated in costFuncObs
            if (upgraded .gt. temp_upgraded) then
                new =  new + 1  ! new is for what?
                if (covexist .eq. 1)then
                    mc_DApar%coefac = mc_DApar%coefnorm                           ! coefac is old parameter sets? coef is the new one; coefnorm 
                    mc_DApar%coefhistory(new, :) = mc_DApar%coefnorm              ! coefhistory used to create new coef matrix
                else
                    mc_DApar%DApar_old = mc_DApar%DApar
                    do ipar = 1, npar4DA
                        mc_DApar%coefnorm(ipar) = (mc_DApar%DApar(ipar)    - mc_DApar%DAparmin(ipar))&
                                                / (mc_DApar%DAparmax(ipar) - mc_DApar%DAparmin(ipar))
                    enddo
                endif
                ! do ipft = 1, npft
                    mc_DApar%coefhistory(new, :) = mc_DApar%coefnorm 
                ! enddo
                if(new .ge. ncov)new=0
                ! update the parameters sets
                ! do ipft = 1, npft
                ! print*, "size: ",  size(arr_params_set%tot_paramsets), nDAsimu,npar4DA
                    arr_params_set%tot_paramsets(upgraded,:) = mc_DApar%DApar
                ! enddo
                mark4scale = 0
                search_scale = init_scale
                nonaccept = 0
            else
                reject = reject + 1
                mark4scale = mark4scale + 1
                ! if(.not. do_cov2createNewPars)then
                !     if (mark4scale > 20) then 
                !         if (nonaccept < 20) then
                !             ! call random_number(rand_scale)
                !             ! search_scale = max(search_scale, min(rand_scale, 0.5))
                !             search_scale = AMIN1(search_scale*1.1, 0.5)
                !             print*, "increase search scale: ", search_scale
                !             nonaccept = nonaccept + 1
                !         ! else:
                !             mark4scale = 0
                !         endif
                !     endif
                ! ! nonaccept = nonaccept + 1
                ! endif
                ! if (accept_rate < 0.1) then 
                !     delta_scale = amax1(AMIN1(delta_scale*0.7, delta_scale_max), delta_scale_min) 
                !     if(do_cov2createNewPars) then
                !         fact_rejet = amin1(amax1(fact_rejet*(1-0.3), 1.), 2.4)
                !         print*, "changed fact_rejet < 0.1: ", fact_rejet
                !     endif
                !     print*, "chanaged delta_scale and fact_reject: ", delta_scale, fact_rejet
                ! else
                !     delta_scale = amax1(AMIN1(delta_scale*1.3, delta_scale_max), delta_scale_min)
                !     if(do_cov2createNewPars) then
                !         fact_rejet = amin1(amax1(fact_rejet*(1.3), 1.), 2.4)
                !         print*, "changed fact_rejet > 0.1: ", fact_rejet
                !     endif
                ! endif
            endif

            ! updates of the covariance matrix
            if(do_cov)then
                if (.not. do_cov2createNewPars .and. mod(upgraded, ncov).eq.0 .and. upgraded .ne. 0)then
                    do_cov2createNewPars = .True.
                    mc_DApar%coefac   = mc_DApar%coefnorm          ! coefnorm: normized values between min and max values
                    call varcov(mc_DApar%coefhistory, mc_DApar%gamnew, npar4DA, ncov) !
                    if(.not.(all(mc_DApar%gamnew==0.)))then
                        mc_DApar%gamma = mc_DApar%gamnew
                        call racine_mat(mc_DApar%gamma, mc_DApar%gamnew, npar4DA)
                        mc_DApar%gamma = mc_DApar%gamnew
                    endif
                endif
                    

                if(mod(upgraded, ncov).eq.0 .and. covexist.eq.1 .and. upgraded .ne. 0)then
                    call varcov(mc_DApar%coefhistory, mc_DApar%gamnew, npar4DA, ncov)
                    if(.not.(all(mc_DApar%gamnew==0.)))then
                        mc_DApar%gamma = mc_DApar%gamnew
                        call racine_mat(mc_DApar%gamma, mc_DApar%gamnew, npar4DA)
                        mc_DApar%gamma = mc_DApar%gamnew
                    endif
                endif
            endif
        enddo

        ! summary
        call mcmc_param_outputs(upgraded, npar4DA, st)!, mc_DApar)  ! update the variable of 
    end subroutine run_mcmc

    subroutine generate_newPar()
        ! This subroutine is used to generate the new parameters to run MCMC
        ! Based on the Shuang's code, it need to use the coef to generate the new parameters.
        implicit none
        ! real(8), intent(in) :: par_old(:), par_min(:), par_max(:)
        ! real(8), intent(inout) :: par_new(:) 
        integer igenPar, parflag, ipft
        real(8) rand_harvest, rand

        call random_seed()

        ! DApar_old = DApar                   ! mark as old parameters 
        if (do_cov2createNewPars) then
            parflag = 1                 ! mark
            do while(parflag .gt. 0)    ! create the new coefnorm
                ! create the new coefnorm based on old one of coefac
                call gengaussvect(fact_rejet*mc_DApar%gamma, mc_DApar%coefac, &
                        mc_DApar%coefnorm, npar4DA)          ! generate the new cov parameters
                parflag = 0
                do igenPar = 1, npar4DA                                                ! check the cov 
                    if(mc_DApar%coefnorm(igenPar).lt.0. .or. mc_DApar%coefnorm(igenPar).gt.1.)then
                        parflag=parflag+1
                        ! write(*,*)'out of range',igenPar
                    endif
                enddo
            enddo
            ! create the new parameters from 
            do ipar = 1, npar4DA
                mc_DApar%DApar(ipar) = mc_DApar%DAparmin(ipar) + &
                    mc_DApar%coefnorm(ipar) * (mc_DApar%DAparmax(ipar)-mc_DApar%DAparmin(ipar))
            enddo
        else ! do not run cov to create new parameters, just random selections
            do igenPar = 1, npar4DA     ! for each parameters
999             continue
                call random_number(rand_harvest)    
                rand = rand_harvest - 0.5           ! create a random number in [-0.5, 0.5]
                if((mc_DApar%DApar_old(igenPar) .gt. mc_DApar%DAparmax(igenPar)) &
                    &   .or. (mc_DApar%DApar_old(igenPar) .lt. mc_DApar%DAparmin(igenPar))) then 
                        mc_DApar%DApar_old(igenPar) = rand_harvest*(mc_DApar%DAparmax(igenPar) - mc_DApar%DAparmin(igenPar))
                endif
                if(mc_DApar%DAparmax(igenPar) .eq. mc_DApar%DAparmin(igenPar)) then
                    mc_DApar%DApar_old(igenPar) = mc_DApar%DAparmax(igenPar)
                endif
                mc_DApar%DApar(igenPar) = mc_DApar%DApar_old(igenPar) + &
                    rand*(mc_DApar%DAparmax(igenPar) - mc_DApar%DAparmin(igenPar)) * search_scale   ! create new parameter
                if((mc_DApar%DApar(igenPar) .gt. mc_DApar%DAparmax(igenPar)) &
                    &   .or. (mc_DApar%DApar(igenPar) .lt. mc_DApar%DAparmin(igenPar))) then 
                    ! print*, "here 1: ",  igenPar,mc_DApar%DApar(igenPar), mc_DApar%DAparmin(igenPar), mc_DApar%DAparmax(igenPar)
                        ! stop
                    goto 999                  ! judge the range of new parameter
                endif
            enddo
        endif
        return   ! mainly return the DApar, meanwhile update the coefnorm
    end subroutine generate_newPar

    subroutine costFuncObs_update()
        implicit none
        real(8) J_cost, delta_J(12), cs_rand
        integer :: ipft, iaccep(12), i, iupdata, nupdata
        
        J_new = 0
        nupdata = 12

        !! cLeaf + cStem:Tree
        if(vars4MCMC%cPlant_tree_y%existOrNot)then
            call CalculateCost(vars4MCMC%cPlant_tree_y%mdData(:,4), vars4MCMC%cPlant_tree_y%obsData(:,4),&
                 vars4MCMC%cPlant_tree_y%obsData(:,5), J_cost)
            J_new(1) = J_new(1) + J_cost
        endif

        ! ! ANPP_Tree_y
        if(vars4MCMC%ANPP_Tree_y%existOrNot)then
            call CalculateCost(vars4MCMC%ANPP_Tree_y%mdData(:,4), vars4MCMC%ANPP_Tree_y%obsData(:,4),&
                 vars4MCMC%ANPP_Tree_y%obsData(:,5), J_cost)
            J_new(2) = J_new(2) + J_cost
        endif

        ! leaf_mass_shrub_y
        if(vars4MCMC%leaf_mass_shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%leaf_mass_shrub_y%mdData(:,4), vars4MCMC%leaf_mass_shrub_y%obsData(:,4),&
                 vars4MCMC%leaf_mass_shrub_y%obsData(:,5), J_cost)
            ! print*, vars4MCMC%leaf_mass_shrub_y%mdData(:,4), vars4MCMC%leaf_mass_shrub_y%obsData(:,4)
            J_new(3) = J_new(3) + J_cost*10
        endif

        ! stem_mass_shrub_y
        if(vars4MCMC%stem_mass_shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%stem_mass_shrub_y%mdData(:,4), vars4MCMC%stem_mass_shrub_y%obsData(:,4),&
                 vars4MCMC%stem_mass_shrub_y%obsData(:,5), J_cost)
            J_new(4) = J_new(4) + J_cost
        endif

        ! ! ! ANPP_Shrub_y
        if(vars4MCMC%ANPP_Shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%ANPP_Shrub_y%mdData(:,4), vars4MCMC%ANPP_Shrub_y%obsData(:,4),&
                 vars4MCMC%ANPP_Shrub_y%obsData(:,5), J_cost)
            J_new(5) = J_new(5) + J_cost ! test *10 to better constrain
        endif

        ! BNPP_y  ! tree + shrub
        if(vars4MCMC%BNPP_y%existOrNot)then
            call CalculateCost(vars4MCMC%BNPP_y%mdData(:,4), vars4MCMC%BNPP_y%obsData(:,4),&
                 vars4MCMC%BNPP_y%obsData(:,5), J_cost)
            J_new(6) = J_new(6) + J_cost*10
        endif

        ! C plant sphagnum
        if(vars4MCMC%cPlant_sphag_y%existOrNot)then
            call CalculateCost(vars4MCMC%cPlant_sphag_y%mdData(:,4), vars4MCMC%cPlant_sphag_y%obsData(:,4),&
                 vars4MCMC%cPlant_sphag_y%obsData(:,5), J_cost)
            J_new(7) = J_new(7) + J_cost
        endif

        ! ! ! NPP_sphag_y
        if(vars4MCMC%NPP_sphag_y%existOrNot)then
            call CalculateCost(vars4MCMC%NPP_sphag_y%mdData(:,4), vars4MCMC%NPP_sphag_y%obsData(:,4),&
                 vars4MCMC%NPP_sphag_y%obsData(:,5), J_cost)
            J_new(8) = J_new(8) + J_cost*10
        endif  

        if(vars4MCMC%gpp_d%existOrNot)then
            call CalculateCost(vars4MCMC%gpp_d%mdData(:,4), vars4MCMC%gpp_d%obsData(:,4),&
                 vars4MCMC%gpp_d%obsData(:,5), J_cost)
            J_new(9) = J_new(9) + J_cost!*1000
        endif

        if(vars4MCMC%er_h%existOrNot)then
            call CalculateCost(vars4MCMC%er_h%mdData(:,4), vars4MCMC%er_h%obsData(:,4),&
                 vars4MCMC%er_h%obsData(:,5), J_cost)
            J_new(10) = J_new(10) + J_cost!*1000
        endif

        ! nee_h         ! shrub + sphag.
        ! if(vars4MCMC%nee_h%existOrNot)then
        !     call CalculateCost(vars4MCMC%nee_h%mdData(:,4), vars4MCMC%nee_h%obsData(:,4),&
        !          vars4MCMC%nee_h%obsData(:,5), J_cost)
        !     J_new(11) = J_new(11) + J_cost!*1000
        ! endif

        if(vars4MCMC%cSoil_y%existOrNot)then
            call CalculateCost(vars4MCMC%cSoil_y%mdData(:,4), vars4MCMC%cSoil_y%obsData(:,4),&
                 vars4MCMC%cSoil_y%obsData(:,5), J_cost)
            J_new(12) = J_new(12) + J_cost/10
        endif

        ! ch4_h
        ! if(vars4MCMC%ch4_h%existOrNot)then
        !     call CalculateCost(vars4MCMC%ch4_h%mdData(:,4), vars4MCMC%ch4_h%obsData(:,4),&
        !          vars4MCMC%ch4_h%obsData(:,5), J_cost)
        !     J_new(13) = J_new(13) + J_cost
        ! endif
        ! =====================================================================================

        

        ! ! ! ! ! LAI_d         ! tree  + Shrub
        ! if(vars4MCMC%LAI_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%LAI_d%mdData(:,4), vars4MCMC%LAI_d%obsData(:,4),&
        !          vars4MCMC%LAI_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif

        ! ! ! leaf_resp_shrub_d 
        ! ! if(vars4MCMC%leaf_resp_shrub_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%leaf_resp_shrub_d%mdData(:,4), vars4MCMC%leaf_resp_shrub_d%obsData(:,4),&
        ! !          vars4MCMC%leaf_resp_shrub_d%obsData(:,5), J_cost)
        ! !     J_new(2) = J_new(2) + J_cost
        ! ! endif
        ! ! ! print*, "J_new13: ", J_new

        ! ! ! leaf_resp_tree_d 
        ! ! if(vars4MCMC%leaf_resp_tree_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%leaf_resp_tree_d%mdData(:,4), vars4MCMC%leaf_resp_tree_d%obsData(:,4),&
        ! !          vars4MCMC%leaf_resp_tree_d%obsData(:,5), J_cost)
        ! !     J_new(1) = J_new(1) + J_cost
        ! ! endif

         ! ! ! photo_shrub_d 
        ! ! if(vars4MCMC%photo_shrub_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%photo_shrub_d%mdData(:,4), vars4MCMC%photo_shrub_d%obsData(:,4),&
        ! !          vars4MCMC%photo_shrub_d%obsData(:,5), J_cost)
        ! !     J_new(2) = J_new(2) + J_cost
        ! ! endif
        ! ! ! ! print*, "J_new18: ", J_new

        ! ! ! photo_tree_d 
        ! ! if(vars4MCMC%photo_tree_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%photo_tree_d%mdData(:,4), vars4MCMC%photo_tree_d%obsData(:,4),&
        ! !          vars4MCMC%photo_tree_d%obsData(:,5), J_cost)
        ! !     J_new(1) = J_new(1) + J_cost
        ! ! endif

        ! if(vars4MCMC%er_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%er_d%mdData(:,4), vars4MCMC%er_d%obsData(:,4),&
        !          vars4MCMC%er_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif

        ! nee_d         ! Shrub + sphag.
        ! if(vars4MCMC%nee_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%nee_d%mdData(:,4), vars4MCMC%nee_d%obsData(:,4),&
        !          vars4MCMC%nee_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif


        ! ! ! ch4_d 
        ! ! if(vars4MCMC%ch4_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%ch4_d%mdData(:,4), vars4MCMC%ch4_d%obsData(:,4),&
        ! !          vars4MCMC%ch4_d%obsData(:,5), J_cost)
        ! !     J_new(5) = J_new(5) + J_cost
        ! ! endif
        ! ! ! print*, "J_new15: ", J_new
        
        ! ! ! print*, "J_new16: ", J_new
        
        ! ! ! CN_shag_d 
        ! ! if(vars4MCMC%CN_shag_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%CN_shag_d%mdData(:,4), vars4MCMC%CN_shag_d%obsData(:,4),&
        ! !          vars4MCMC%CN_shag_d%obsData(:,5), J_cost)
        ! !     J_new(3) = J_new(3) + J_cost
        ! ! endif
        ! ! ! print*, "J_new17: ", J_new

       
        ! ! print*, "J_new19: ", J_new
        ! ! ------------------------------------------------------------------------------------
        ! ! write(*,*) "here2",J_new
        iaccep = 0
        nupdata = 12 
        do i = 1,12
            ! if (iDAsimu .eq. i*nDAsimu/nupdata+1) then
            !    if(i<12) J_last(i+1) = J_new(i+1)
            ! endif
            if(J_new(i) .eq. 0) then ! no data is available
                delta_J(i) = -0.1
            else
                delta_J(i) = J_new(i) - J_last(i)
            endif

            delta_J(i) = delta_J(i)

            call random_number(cs_rand)
            if(AMIN1(1.0, exp(-delta_J(i))) .gt. cs_rand)then
                iaccep(i) = iaccep(i) + 1
                if (iDAsimu < 0.8*nDAsimu)then
                    call update_params(i)
                endif
            endif
        enddo
        ! if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
        !     ! if (AMIN1(1.0, exp(-delta_J(5))) .gt. cs_rand) then
        !         ! if (sum(iaccep(1:5)) .gt. 4)then
        !             upgraded = upgraded + 1
        !             J_last = J_new
        !         ! endif
        !     ! endif
        ! endif
        ! iupdata = 12
        
        ! print*, iDAsimu, nDAsimu, iupdata
        ! if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
        !     upgraded = upgraded + 1
        !     J_last = J_new
        ! endif


        if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
        ! if(AMIN1(1.0, exp(-sum(delta_J(1:2)))) .gt. cs_rand)then
            ! if (AMIN1(1.0, exp(-sum(delta_J(3:5)))) .gt. cs_rand)then
            !     if(AMIN1(1.0, exp(-delta_J(6))) .gt. cs_rand)then
            !         if(AMIN1(1.0, exp(-sum(delta_J(7:8)))) .gt. cs_rand)then
                        if(AMIN1(1.0, exp(-sum(delta_J(9:11)))) .gt. cs_rand)then
            !                 if(AMIN1(1.0, exp(-delta_J(12))) .gt. cs_rand)then
            upgraded = upgraded + 1
            J_last = J_new
            !                 endif
                        endif
            !         endif
            !     endif
            ! endif
        endif
            
        J_show_new = sum(J_new)
        J_show_old = sum(J_last) 
        accept_rate = real(upgraded)/real(iDAsimu)
    end subroutine costFuncObs_update

    subroutine update_params(i)
        implicit none
        integer, intent(in) :: i
        if ((i == 1) .or. (i == 2)) call generate_newPar_ipft(26, 49)
        if ((i == 3) .or. (i == 4) .or. (i==5)) then
            call generate_newPar_ipft(50, 73)
        endif
        if ((i == 7) .or. (i==8)) call generate_newPar_ipft(74, 97)
    end subroutine update_params

    subroutine generate_newPar_ipft(istart, iend)
        ! This subroutine is used to generate the new parameters to run MCMC
        ! Based on the Shuang's code, it need to use the coef to generate the new parameters.
        implicit none
        ! real(8), intent(in) :: par_old(:), par_min(:), par_max(:)
        ! real(8), intent(inout) :: par_new(:) 
        integer igenPar, parflag, ipft, istart, iend
        real(8) rand_harvest, rand

        call random_seed()

        ! DApar_old = DApar                   ! mark as old parameters 
        if (do_cov2createNewPars) then
            parflag = 1                 ! mark
            do while(parflag .gt. 0)    ! create the new coefnorm
                ! create the new coefnorm based on old one of coefac
                call gengaussvect(fact_rejet*mc_DApar%gamma, mc_DApar%coefac, &
                        mc_DApar%coefnorm, npar4DA)          ! generate the new cov parameters
                parflag = 0
                do igenPar = istart, iend                                                ! check the cov 
                    if(mc_DApar%coefnorm(igenPar).lt.0. .or. mc_DApar%coefnorm(igenPar).gt.1.)then
                        parflag=parflag+1
                        ! write(*,*)'out of range',igenPar
                    endif
                enddo
            enddo
            ! create the new parameters from 
            do ipar = istart, iend
                mc_DApar%DApar(ipar) = mc_DApar%DAparmin(ipar) + &
                    mc_DApar%coefnorm(ipar) * (mc_DApar%DAparmax(ipar)-mc_DApar%DAparmin(ipar))
            enddo
        else ! do not run cov to create new parameters, just random selections
            do igenPar = istart, iend     ! for each parameters
999             continue
                call random_number(rand_harvest)    
                rand = rand_harvest - 0.5           ! create a random number in [-0.5, 0.5]
                if((mc_DApar%DApar_old(igenPar) .gt. mc_DApar%DAparmax(igenPar)) &
                    &   .or. (mc_DApar%DApar_old(igenPar) .lt. mc_DApar%DAparmin(igenPar))) then 
                        mc_DApar%DApar_old(igenPar) = rand_harvest*(mc_DApar%DAparmax(igenPar) - mc_DApar%DAparmin(igenPar))
                endif
                if(mc_DApar%DAparmax(igenPar) .eq. mc_DApar%DAparmin(igenPar)) then
                    mc_DApar%DApar_old(igenPar) = mc_DApar%DAparmax(igenPar)
                endif
                mc_DApar%DApar(igenPar) = mc_DApar%DApar_old(igenPar) + &
                    rand*(mc_DApar%DAparmax(igenPar) - mc_DApar%DAparmin(igenPar)) * search_scale   ! create new parameter
                if((mc_DApar%DApar(igenPar) .gt. mc_DApar%DAparmax(igenPar)) &
                    &   .or. (mc_DApar%DApar(igenPar) .lt. mc_DApar%DAparmin(igenPar))) then 
                    ! print*, "here 1: ",  igenPar,mc_DApar%DApar(igenPar), mc_DApar%DAparmin(igenPar), mc_DApar%DAparmax(igenPar)
                        ! stop
                    goto 999                  ! judge the range of new parameter
                endif
            enddo
        endif
        return   ! mainly return the DApar, meanwhile update the coefnorm
    end subroutine generate_newPar_ipft


    subroutine costFuncObs_old()
        implicit none
        real(8) J_cost, delta_J(20), cs_rand, delta_J_new
        integer :: ipft, iaccep(20), i, iupdata, nupdata, iiii
        
        J_new = 0
        nupdata = 12

        ! cLeaf + cStem:Tree
        if(vars4MCMC%cPlant_tree_y%existOrNot)then
            call CalculateCost(vars4MCMC%cPlant_tree_y%mdData(:,4), vars4MCMC%cPlant_tree_y%obsData(:,4),&
                 vars4MCMC%cPlant_tree_y%obsData(:,5), J_cost)
            J_new(1) = J_new(1) + J_cost*10000
        endif

        ! ! ANPP_Tree_y
        if(vars4MCMC%ANPP_Tree_y%existOrNot)then
            call CalculateCost(vars4MCMC%ANPP_Tree_y%mdData(:,4), vars4MCMC%ANPP_Tree_y%obsData(:,4),&
                 vars4MCMC%ANPP_Tree_y%obsData(:,5), J_cost)
            J_new(2) = J_new(2) + J_cost*1000
        endif

        ! leaf_mass_shrub_y
        if(vars4MCMC%leaf_mass_shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%leaf_mass_shrub_y%mdData(:,4), vars4MCMC%leaf_mass_shrub_y%obsData(:,4),&
                 vars4MCMC%leaf_mass_shrub_y%obsData(:,5), J_cost)
            ! print*, vars4MCMC%leaf_mass_shrub_y%mdData(:,4), vars4MCMC%leaf_mass_shrub_y%obsData(:,4)
            J_new(3) = J_new(3) + J_cost*10000
        endif

        ! stem_mass_shrub_y
        if(vars4MCMC%stem_mass_shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%stem_mass_shrub_y%mdData(:,4), vars4MCMC%stem_mass_shrub_y%obsData(:,4),&
                 vars4MCMC%stem_mass_shrub_y%obsData(:,5), J_cost)
            J_new(4) = J_new(4) + J_cost*10000
        endif

        ! ! ! ANPP_Shrub_y
        if(vars4MCMC%ANPP_Shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%ANPP_Shrub_y%mdData(:,4), vars4MCMC%ANPP_Shrub_y%obsData(:,4),&
                 vars4MCMC%ANPP_Shrub_y%obsData(:,5), J_cost)
            J_new(5) = J_new(5) + J_cost*10000 ! test *10 to better constrain
        endif

        ! BNPP_y  ! tree + shrub
        if(vars4MCMC%BNPP_y%existOrNot)then
            call CalculateCost(vars4MCMC%BNPP_y%mdData(:,4), vars4MCMC%BNPP_y%obsData(:,4)*1.3,&
                 vars4MCMC%BNPP_y%obsData(:,5), J_cost)
            J_new(6) = J_new(6) + J_cost*20000
        endif

        ! C plant sphagnum
        if(vars4MCMC%cPlant_sphag_y%existOrNot)then
            call CalculateCost(vars4MCMC%cPlant_sphag_y%mdData(:,4), vars4MCMC%cPlant_sphag_y%obsData(:,4),&
                 vars4MCMC%cPlant_sphag_y%obsData(:,5), J_cost)
            J_new(7) = J_new(7) + J_cost*10000
        endif

        ! ! ! NPP_sphag_y
        if(vars4MCMC%NPP_sphag_y%existOrNot)then
            call CalculateCost(vars4MCMC%NPP_sphag_y%mdData(:,4), vars4MCMC%NPP_sphag_y%obsData(:,4),&
                 vars4MCMC%NPP_sphag_y%obsData(:,5), J_cost)
            J_new(8) = J_new(8) + J_cost*20000
        endif  
        ! print*, vars4MCMC%NPP_sphag_y%mdData(:,4), vars4MCMC%NPP_sphag_y%obsData(:,4)
        ! print*, J_new(8), J_last(8)

        ! print*,"before:", J_new
        if(vars4MCMC%gpp_d%existOrNot)then
            call CalculateCost(vars4MCMC%gpp_d%mdData(:,4), vars4MCMC%gpp_d%obsData(:,4),&
                 vars4MCMC%gpp_d%obsData(:,5), J_cost)
            J_new(9) = J_new(9) + J_cost*10000
        endif
        ! print*, vars4MCMC%gpp_d%mdData(:,4), vars4MCMC%gpp_d%obsData(:,4)
        ! print*,"after:", J_new

        if(vars4MCMC%er_h%existOrNot)then
            call CalculateCost(vars4MCMC%er_h%mdData(:,4), vars4MCMC%er_h%obsData(:,4),&
                 vars4MCMC%er_h%obsData(:,5), J_cost)
            J_new(10) = J_new(10) + J_cost*50000
        endif

        ! nee_h         ! shrub + sphag.
        if(vars4MCMC%nee_h%existOrNot)then
            call CalculateCost(vars4MCMC%nee_h%mdData(:,4), vars4MCMC%nee_h%obsData(:,4),&
                 vars4MCMC%nee_h%obsData(:,5), J_cost)
            J_new(11) = J_new(11) + J_cost*10000
        endif

        if(vars4MCMC%cSoil_y%existOrNot)then
            call CalculateCost(vars4MCMC%cSoil_y%mdData(:,4), vars4MCMC%cSoil_y%obsData(:,4),&
                 vars4MCMC%cSoil_y%obsData(:,5), J_cost)
            J_new(12) = J_new(12) + J_cost*1000
        endif

        ! ch4_h
        if(vars4MCMC%ch4_h%existOrNot)then
            call CalculateCost(vars4MCMC%ch4_h%mdData(:,4)*1000000, vars4MCMC%ch4_h%obsData(:,4)*1000000,&
                 vars4MCMC%ch4_h%obsData(:,5), J_cost)
            J_new(13) = J_new(13) + J_cost*5
        endif

        if(vars4MCMC%rh_y%existOrNot)then
            call CalculateCost(vars4MCMC%rh_y%mdData(:,4), vars4MCMC%rh_y%obsData(:,4),&
                 vars4MCMC%rh_y%obsData(:,5), J_cost)
            J_new(14) = J_new(14) + J_cost*30000
        endif

        if(vars4MCMC%ch4_y%existOrNot)then
            call CalculateCost(vars4MCMC%ch4_y%mdData(:,4), vars4MCMC%ch4_y%obsData(:,4),&
                 vars4MCMC%ch4_y%obsData(:,5), J_cost)
            J_new(15) = J_new(15) + J_cost*100
        endif

        ! water table
        if(vars4MCMC%zwt_h%existOrNot)then
            call CalculateCost(vars4MCMC%zwt_h%mdData(:,4)*150, vars4MCMC%zwt_h%obsData(:,4)*150,&
                 vars4MCMC%zwt_h%obsData(:,5), J_cost)
            ! print*, size(vars4MCMC%zwt_h%mdData(:,4)), size(vars4MCMC%zwt_h%obsData(:,4))
!                  do iiii = 1, size(vars4MCMC%zwt_h%mdData(:,4))
! print*, vars4MCMC%zwt_h%mdData(iiii,1), vars4MCMC%zwt_h%mdData(iiii,2), &
! & vars4MCMC%zwt_h%mdData(iiii,3), vars4MCMC%zwt_h%mdData(iiii,4), &
! & vars4MCMC%zwt_h%obsData(iiii,1),vars4MCMC%zwt_h%obsData(iiii,2), &
! & vars4MCMC%zwt_h%obsData(iiii,3),vars4MCMC%zwt_h%obsData(iiii,4)
! enddo
            J_new(16) = J_new(16) + J_cost
        endif



        ! LAI_shrub_d        
        if(vars4MCMC%lai_shrub_d%existOrNot)then
            call CalculateCost(vars4MCMC%lai_shrub_d%mdData(:,4), vars4MCMC%lai_shrub_d%obsData(:,4),&
                 vars4MCMC%lai_shrub_d%obsData(:,5), J_cost)
            J_new(17) = J_new(17) + J_cost*500
        endif

        ! LAI_tree_d        
        if(vars4MCMC%lai_tree_d%existOrNot)then
            call CalculateCost(vars4MCMC%lai_tree_d%mdData(:,4), vars4MCMC%lai_tree_d%obsData(:,4),&
                 vars4MCMC%lai_tree_d%obsData(:,5), J_cost)
            J_new(18) = J_new(18) + J_cost*500
        endif

        ! photo_tree_h       
        if(vars4MCMC%photo_tree_h%existOrNot)then
            call CalculateCost(vars4MCMC%photo_tree_h%mdData(:,4), vars4MCMC%photo_tree_h%obsData(:,4),&
                 vars4MCMC%photo_tree_h%obsData(:,5), J_cost)
            J_new(19) = J_new(19) + J_cost*500
        endif

        ! photo_shrub_d        
        if(vars4MCMC%photo_shrub_d%existOrNot)then
            call CalculateCost(vars4MCMC%photo_shrub_d%mdData(:,4), vars4MCMC%photo_shrub_d%obsData(:,4),&
                 vars4MCMC%photo_shrub_d%obsData(:,5), J_cost)
            J_new(20) = J_new(20) + 0!J_cost
        endif
        ! =====================================================================================

        

        ! ! ! ! ! LAI_d         ! tree  + Shrub
        ! if(vars4MCMC%LAI_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%LAI_d%mdData(:,4), vars4MCMC%LAI_d%obsData(:,4),&
        !          vars4MCMC%LAI_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif

        ! ! ! leaf_resp_shrub_d 
        ! ! if(vars4MCMC%leaf_resp_shrub_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%leaf_resp_shrub_d%mdData(:,4), vars4MCMC%leaf_resp_shrub_d%obsData(:,4),&
        ! !          vars4MCMC%leaf_resp_shrub_d%obsData(:,5), J_cost)
        ! !     J_new(2) = J_new(2) + J_cost
        ! ! endif
        ! ! ! print*, "J_new13: ", J_new

        ! ! ! leaf_resp_tree_d 
        ! ! if(vars4MCMC%leaf_resp_tree_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%leaf_resp_tree_d%mdData(:,4), vars4MCMC%leaf_resp_tree_d%obsData(:,4),&
        ! !          vars4MCMC%leaf_resp_tree_d%obsData(:,5), J_cost)
        ! !     J_new(1) = J_new(1) + J_cost
        ! ! endif

         ! ! ! photo_shrub_d 
        ! ! if(vars4MCMC%photo_shrub_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%photo_shrub_d%mdData(:,4), vars4MCMC%photo_shrub_d%obsData(:,4),&
        ! !          vars4MCMC%photo_shrub_d%obsData(:,5), J_cost)
        ! !     J_new(2) = J_new(2) + J_cost
        ! ! endif
        ! ! ! ! print*, "J_new18: ", J_new

        ! ! ! photo_tree_d 
        ! ! if(vars4MCMC%photo_tree_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%photo_tree_d%mdData(:,4), vars4MCMC%photo_tree_d%obsData(:,4),&
        ! !          vars4MCMC%photo_tree_d%obsData(:,5), J_cost)
        ! !     J_new(1) = J_new(1) + J_cost
        ! ! endif

        ! if(vars4MCMC%er_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%er_d%mdData(:,4), vars4MCMC%er_d%obsData(:,4),&
        !          vars4MCMC%er_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif

        ! nee_d         ! Shrub + sphag.
        ! if(vars4MCMC%nee_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%nee_d%mdData(:,4), vars4MCMC%nee_d%obsData(:,4),&
        !          vars4MCMC%nee_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif


        ! ! ! ch4_d 
        ! ! if(vars4MCMC%ch4_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%ch4_d%mdData(:,4), vars4MCMC%ch4_d%obsData(:,4),&
        ! !          vars4MCMC%ch4_d%obsData(:,5), J_cost)
        ! !     J_new(5) = J_new(5) + J_cost
        ! ! endif
        ! ! ! print*, "J_new15: ", J_new
        
        ! ! ! print*, "J_new16: ", J_new
        
        ! ! ! CN_shag_d 
        ! ! if(vars4MCMC%CN_shag_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%CN_shag_d%mdData(:,4), vars4MCMC%CN_shag_d%obsData(:,4),&
        ! !          vars4MCMC%CN_shag_d%obsData(:,5), J_cost)
        ! !     J_new(3) = J_new(3) + J_cost
        ! ! endif
        ! ! ! print*, "J_new17: ", J_new

       
        ! ! print*, "J_new19: ", J_new
        ! ! ------------------------------------------------------------------------------------
        ! ! write(*,*) "here2",J_new
        iaccep = 0
        nupdata = 12 
        do i = 1,20
            ! if (iDAsimu .eq. i*nDAsimu/nupdata+1) then
            !    if(i<12) J_last(i+1) = J_new(i+1)
            ! endif
            if(J_new(i) .eq. 0) then ! no data is available
                delta_J(i) = -0.1
            else
                delta_J(i) = J_new(i) - J_last(i)
            endif

            delta_J(i) = delta_J(i)

            call random_number(cs_rand)
            if(AMIN1(1.0, exp(-delta_J(i))) .gt. cs_rand)then
                iaccep(i) = iaccep(i) + 1
            endif
        enddo
        ! if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
        !     ! if (AMIN1(1.0, exp(-delta_J(5))) .gt. cs_rand) then
        !         ! if (sum(iaccep(1:5)) .gt. 4)then
        !             upgraded = upgraded + 1
        !             J_last = J_new
        !         ! endif
        !     ! endif
        ! endif
        ! iupdata = 12
        
        ! print*, iDAsimu, nDAsimu, iupdata
        ! if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
        !     upgraded = upgraded + 1
        !     J_last = J_new
        ! endif
        delta_J_new = (sum(J_new(1:15)) - sum(J_last(1:15)))!/15 * delta_scale!0.05
        ! if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
        if(AMIN1(1.0, exp(-delta_J_new)) .gt. cs_rand)then
        ! if(AMIN1(1.0, exp(-sum(delta_J(1:2)))) .gt. cs_rand)then
            ! if (AMIN1(1.0, exp(-sum(delta_J(3:5)))) .gt. cs_rand)then
            !     if(AMIN1(1.0, exp(-delta_J(6))) .gt. cs_rand)then
            !         if(AMIN1(1.0, exp(-sum(delta_J(7:8)))) .gt. cs_rand)then
                        ! if(AMIN1(1.0, exp(-delta_J(9))) .gt. cs_rand)then
                        !     if(AMIN1(1.0, exp(-delta_J(10))) .gt. cs_rand)then
            !                 if(AMIN1(1.0, exp(-delta_J(12))) .gt. cs_rand)then
            upgraded = upgraded + 1
            J_last = J_new
                            ! endif
                        ! endif
            !         endif
            !     endif
            ! endif
        endif
        ! print*, J_new, J_last
            
        J_show_new = sum(J_new)
        J_show_old = sum(J_last) 
        accept_rate = real(upgraded)/real(iDAsimu)

    end subroutine costFuncObs_old

    subroutine costFuncObs_new()
        implicit none
        real(8) J_cost, delta_J(12), cs_rand
        integer :: ipft, iaccep(12), i, iupdata, nupdata
        
        J_new = 0

        !! cLeaf + cStem
        if(vars4MCMC%cPlant_tree_y%existOrNot)then
            call CalculateCost(vars4MCMC%cPlant_tree_y%mdData(:,4), vars4MCMC%cPlant_tree_y%obsData(:,4),&
                 vars4MCMC%cPlant_tree_y%obsData(:,5), J_cost)
            J_new(1) = J_new(1) + J_cost*10
        endif

        ! ! ANPP_Tree_y
        if(vars4MCMC%ANPP_Tree_y%existOrNot)then
            call CalculateCost(vars4MCMC%ANPP_Tree_y%mdData(:,4), vars4MCMC%ANPP_Tree_y%obsData(:,4),&
                 vars4MCMC%ANPP_Tree_y%obsData(:,5), J_cost)
            J_new(2) = J_new(2) + J_cost*10
        endif

        ! leaf_mass_shrub_y
        if(vars4MCMC%leaf_mass_shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%leaf_mass_shrub_y%mdData(:,4), vars4MCMC%leaf_mass_shrub_y%obsData(:,4),&
                 vars4MCMC%leaf_mass_shrub_y%obsData(:,5), J_cost)
            J_new(3) = J_new(3) + J_cost*10
        endif

        ! stem_mass_shrub_y
        if(vars4MCMC%stem_mass_shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%stem_mass_shrub_y%mdData(:,4), vars4MCMC%stem_mass_shrub_y%obsData(:,4),&
                 vars4MCMC%stem_mass_shrub_y%obsData(:,5), J_cost)
            J_new(4) = J_new(4) + J_cost*10
        endif

        ! ! ! ANPP_Shrub_y
        if(vars4MCMC%ANPP_Shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%ANPP_Shrub_y%mdData(:,4), vars4MCMC%ANPP_Shrub_y%obsData(:,4),&
                 vars4MCMC%ANPP_Shrub_y%obsData(:,5), J_cost)
            J_new(5) = J_new(5) + J_cost*10 ! test *10 to better constrain
        endif

        ! BNPP_y  ! tree + shrub
        if(vars4MCMC%BNPP_y%existOrNot)then
            call CalculateCost(vars4MCMC%BNPP_y%mdData(:,4), vars4MCMC%BNPP_y%obsData(:,4),&
                 vars4MCMC%BNPP_y%obsData(:,5), J_cost)
            J_new(6) = J_new(6) + J_cost*100
        endif

        ! C plant sphagnum
        if(vars4MCMC%cPlant_sphag_y%existOrNot)then
            call CalculateCost(vars4MCMC%cPlant_sphag_y%mdData(:,4), vars4MCMC%cPlant_sphag_y%obsData(:,4),&
                 vars4MCMC%cPlant_sphag_y%obsData(:,5), J_cost)
            J_new(7) = J_new(7) + J_cost*10
        endif

        ! ! ! NPP_sphag_y
        if(vars4MCMC%NPP_sphag_y%existOrNot)then
            call CalculateCost(vars4MCMC%NPP_sphag_y%mdData(:,4), vars4MCMC%NPP_sphag_y%obsData(:,4),&
                 vars4MCMC%NPP_sphag_y%obsData(:,5), J_cost)
            J_new(8) = J_new(8) + J_cost*10
        endif  

        if(vars4MCMC%gpp_d%existOrNot)then
            call CalculateCost(vars4MCMC%gpp_d%mdData(:,4), vars4MCMC%gpp_d%obsData(:,4),&
                 vars4MCMC%gpp_d%obsData(:,5), J_cost)
            J_new(9) = J_new(9) + J_cost*1000
        endif

        if(vars4MCMC%er_h%existOrNot)then
            call CalculateCost(vars4MCMC%er_h%mdData(:,4), vars4MCMC%er_h%obsData(:,4),&
                 vars4MCMC%er_h%obsData(:,5), J_cost)
            J_new(10) = J_new(10) + J_cost*1000
        endif

        ! nee_h         ! shrub + sphag.
        if(vars4MCMC%nee_h%existOrNot)then
            call CalculateCost(vars4MCMC%nee_h%mdData(:,4), vars4MCMC%nee_h%obsData(:,4),&
                 vars4MCMC%nee_h%obsData(:,5), J_cost)
            J_new(11) = J_new(11) + J_cost*1000
        endif

        if(vars4MCMC%cSoil_y%existOrNot)then
            call CalculateCost(vars4MCMC%cSoil_y%mdData(:,4), vars4MCMC%cSoil_y%obsData(:,4),&
                 vars4MCMC%cSoil_y%obsData(:,5), J_cost)
            J_new(12) = J_new(12) + J_cost
        endif

        ! ch4_h
        ! if(vars4MCMC%ch4_h%existOrNot)then
        !     call CalculateCost(vars4MCMC%ch4_h%mdData(:,4), vars4MCMC%ch4_h%obsData(:,4),&
        !          vars4MCMC%ch4_h%obsData(:,5), J_cost)
        !     J_new(13) = J_new(13) + J_cost
        ! endif
        ! =====================================================================================

        

        ! ! ! ! ! LAI_d         ! tree  + Shrub
        ! if(vars4MCMC%LAI_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%LAI_d%mdData(:,4), vars4MCMC%LAI_d%obsData(:,4),&
        !          vars4MCMC%LAI_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif

        ! ! ! leaf_resp_shrub_d 
        ! ! if(vars4MCMC%leaf_resp_shrub_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%leaf_resp_shrub_d%mdData(:,4), vars4MCMC%leaf_resp_shrub_d%obsData(:,4),&
        ! !          vars4MCMC%leaf_resp_shrub_d%obsData(:,5), J_cost)
        ! !     J_new(2) = J_new(2) + J_cost
        ! ! endif
        ! ! ! print*, "J_new13: ", J_new

        ! ! ! leaf_resp_tree_d 
        ! ! if(vars4MCMC%leaf_resp_tree_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%leaf_resp_tree_d%mdData(:,4), vars4MCMC%leaf_resp_tree_d%obsData(:,4),&
        ! !          vars4MCMC%leaf_resp_tree_d%obsData(:,5), J_cost)
        ! !     J_new(1) = J_new(1) + J_cost
        ! ! endif

         ! ! ! photo_shrub_d 
        ! ! if(vars4MCMC%photo_shrub_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%photo_shrub_d%mdData(:,4), vars4MCMC%photo_shrub_d%obsData(:,4),&
        ! !          vars4MCMC%photo_shrub_d%obsData(:,5), J_cost)
        ! !     J_new(2) = J_new(2) + J_cost
        ! ! endif
        ! ! ! ! print*, "J_new18: ", J_new

        ! ! ! photo_tree_d 
        ! ! if(vars4MCMC%photo_tree_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%photo_tree_d%mdData(:,4), vars4MCMC%photo_tree_d%obsData(:,4),&
        ! !          vars4MCMC%photo_tree_d%obsData(:,5), J_cost)
        ! !     J_new(1) = J_new(1) + J_cost
        ! ! endif

        ! if(vars4MCMC%er_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%er_d%mdData(:,4), vars4MCMC%er_d%obsData(:,4),&
        !          vars4MCMC%er_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif

        ! nee_d         ! Shrub + sphag.
        ! if(vars4MCMC%nee_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%nee_d%mdData(:,4), vars4MCMC%nee_d%obsData(:,4),&
        !          vars4MCMC%nee_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif


        ! ! ! ch4_d 
        ! ! if(vars4MCMC%ch4_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%ch4_d%mdData(:,4), vars4MCMC%ch4_d%obsData(:,4),&
        ! !          vars4MCMC%ch4_d%obsData(:,5), J_cost)
        ! !     J_new(5) = J_new(5) + J_cost
        ! ! endif
        ! ! ! print*, "J_new15: ", J_new
        
        ! ! ! print*, "J_new16: ", J_new
        
        ! ! ! CN_shag_d 
        ! ! if(vars4MCMC%CN_shag_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%CN_shag_d%mdData(:,4), vars4MCMC%CN_shag_d%obsData(:,4),&
        ! !          vars4MCMC%CN_shag_d%obsData(:,5), J_cost)
        ! !     J_new(3) = J_new(3) + J_cost
        ! ! endif
        ! ! ! print*, "J_new17: ", J_new

       
        ! ! print*, "J_new19: ", J_new
        ! ! ------------------------------------------------------------------------------------
        ! ! write(*,*) "here2",J_new
        iaccep = 0
        nupdata = 13
        if (iDAsimu .eq. 1*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 2*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 3*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 4*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 5*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 6*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 7*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 8*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 9*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 10*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 11*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. (nupdata-1)*nDAsimu/nupdata+1) J_last = J_new
        
        do i = 1,12
            ! if (iDAsimu .eq. i*nDAsimu/nupdata+1) then
            !    if(i<12) J_last(i+1) = J_new(i+1)
            ! endif
            if(J_new(i) .eq. 0) then ! no data is available
                delta_J(i) = -0.1
            else
                delta_J(i) = J_new(i) - J_last(i)
            endif

            delta_J(i) = delta_J(i)

            call random_number(cs_rand)
            if(AMIN1(1.0, exp(-delta_J(i))) .gt. cs_rand)then
                iaccep(i) = iaccep(i) + 1
            endif
        enddo
        ! if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
        !     ! if (AMIN1(1.0, exp(-delta_J(5))) .gt. cs_rand) then
        !         ! if (sum(iaccep(1:5)) .gt. 4)then
        !             upgraded = upgraded + 1
        !             J_last = J_new
        !         ! endif
        !     ! endif
        ! endif
        ! iupdata = 12
        
        if(iDAsimu < 1*nDAsimu/nupdata) then 
            ! change to: cplant_sphagnum J(7)
            iupdata = 1
            if(AMIN1(1.0, exp(-delta_J(7))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(7)
            J_show_old = J_last(7)  
        else if(iDAsimu < 2*nDAsimu/nupdata) then
            ! change to: NPP_sphagnum J(8)
            iupdata = 2
            if(AMIN1(1.0, exp(-sum(delta_J(7:8)))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = sum(J_new(7:8))
            J_show_old = sum(J_last(7:8))
        else if(iDAsimu < 3*nDAsimu/nupdata) then
            ! change to: cLeaf Shrub J(3)
            iupdata = 3
            if(AMIN1(1.0, exp(-delta_J(3))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(3)
            J_show_old = J_last(3)
        else if(iDAsimu < 4*nDAsimu/nupdata) then
            iupdata = 4
            if(AMIN1(1.0, exp(-delta_J(4))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(4)
            J_show_old = J_last(4)
        else if(iDAsimu < 5*nDAsimu/nupdata) then
            iupdata = 5
            if(AMIN1(1.0, exp(-sum(delta_J(3:5)))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = sum(J_new(3:5))
            J_show_old = sum(J_last(3:5))
        else if(iDAsimu < 6*nDAsimu/nupdata) then
            ! BNPP
            iupdata = 6
            if(AMIN1(1.0, exp(-delta_J(6))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(6)
            J_show_old = J_last(6)
        else if(iDAsimu < 7*nDAsimu/nupdata) then
            ! GPP
            iupdata = 7
            if(AMIN1(1.0, exp(-delta_J(9))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(9)
            J_show_old = J_last(9)
        else if(iDAsimu < 8*nDAsimu/nupdata) then
            ! ER
            iupdata = 8
            if(AMIN1(1.0, exp(-delta_J(10))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(10)
            J_show_old = J_last(10)
        else if(iDAsimu < 9*nDAsimu/nupdata) then
            ! NEE 
            iupdata = 9
            if(AMIN1(1.0, exp(-sum(delta_J(9:11)))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = sum(J_new(9:11))
            J_show_old = sum(J_last(9:11))
        else if(iDAsimu < 10*nDAsimu/nupdata) then
            ! cplant Tree
            iupdata = 10
            if(AMIN1(1.0, exp(-delta_J(1))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(1)
            J_show_old = J_last(1)
        else if(iDAsimu < 11*nDAsimu/nupdata) then
            ! anpp tree
            iupdata = 11
            if(AMIN1(1.0, exp(-sum(delta_J(1:2)))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = sum(J_new(1:2))
            J_show_old = sum(J_last(1:2))
        else if (iDAsimu < 12*nDAsimu/nupdata) then
            ! csoil
            iupdata = 12
            if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(12)
            J_show_old = J_last(12)
        else
            if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = sum(J_new)
            J_show_old = sum(J_last)
        endif
            
        accept_rate = real(upgraded)/real(iDAsimu)

        ! do i = 1, nupdata
        !     if (iDAsimu .eq. i*nDAsimu/nupdata) then
        !         if (i==1) then 
        !             before_upgrade = 1
        !         else
        !             before_upgrade = after_upgrade
        !         endif
        !         after_upgrade = upgraded
        !         call  update_min_max_val_new_new(i)      
        !     endif
        ! enddo
        if(iDAsimu .eq. 1) then 
            before_upgrade = 1
        else if(iDAsimu .eq. 2*nDAsimu/nupdata) then
            after_upgrade = upgraded
            call update_min_max_val_new_new(2)
        ! else if(iDAsimu .eq. 3*nDAsimu/nupdata)then
        !     before_upgrade = after_upgrade + 1
        !     after_upgrade  = upgraded
        !     call update_min_max_val_new_new(3)
        ! else if(iDAsimu .eq. 4*nDAsimu/nupdata)then
        !     before_upgrade = after_upgrade + 1
        !     after_upgrade  = upgraded
        !     call update_min_max_val_new_new(4)
        else if(iDAsimu .eq. 5*nDAsimu/nupdata)then
            before_upgrade = after_upgrade + 1
            after_upgrade  = upgraded
            call update_min_max_val_new_new(5)
        else if(iDAsimu .eq. 6*nDAsimu/nupdata)then
            before_upgrade = after_upgrade + 1
            after_upgrade  = upgraded
            call update_min_max_val_new_new(6)
        ! else if(iDAsimu .eq. 7*nDAsimu/nupdata)then
        !     before_upgrade = after_upgrade + 1
        !     after_upgrade  = upgraded
        !     call update_min_max_val_new_new(7)
        ! else if(iDAsimu .eq. 8*nDAsimu/nupdata)then
        !     before_upgrade = after_upgrade + 1
        !     after_upgrade  = upgraded
        !     call update_min_max_val_new_new(8)
        else if(iDAsimu .eq. 9*nDAsimu/nupdata)then
            before_upgrade = after_upgrade + 1
            after_upgrade  = upgraded
            call update_min_max_val_new_new(9)
        ! else if(iDAsimu .eq. 10*nDAsimu/nupdata)then
        !     before_upgrade = after_upgrade + 1
        !     after_upgrade  = upgraded
        !     call update_min_max_val_new_new(10)
        else if(iDAsimu .eq. 11*nDAsimu/nupdata)then
            before_upgrade = after_upgrade + 1
            after_upgrade  = upgraded
            call update_min_max_val_new_new(11)
        else if(iDAsimu .eq. 12*nDAsimu/nupdata)then
            before_upgrade = after_upgrade + 1
            after_upgrade  = upgraded
            call update_min_max_val_new_new(12)
        endif
        
    end subroutine costFuncObs_new

    subroutine update_min_max_val_new_new(i)
        implicit none
        integer :: i, igenPar
        integer :: n_chg_pars, nBuilt_in
        real(8) :: new_max, new_min
        n_chg_pars = after_upgrade - before_upgrade
        nBuilt_in  = int(0.1*n_chg_pars)

        if (i == 2) then ! cPlant_sphagnum + NPP sphagnum
            do igenPar = 80, 97
                call update_mc(nBuilt_in, igenPar)
            enddo
        else if(i == 5) then ! cleaf + cstem + anpp shrub
            call update_mc(nBuilt_in, 52)
            ! 56-57, 59-61, 64-70
            call update_mc(nBuilt_in, 56)
            call update_mc(nBuilt_in, 57)
            call update_mc(nBuilt_in, 59)
            call update_mc(nBuilt_in, 60)
            call update_mc(nBuilt_in, 61)
            do igenPar = 64, 70
                call update_mc(nBuilt_in, igenPar)
            enddo
        else if(i == 9) then ! gpp +er+nee
            ! 1-3, 5-14, 23-25, 28,32-48, 56-57, 64-73, 80-97 
            do igenPar = 1,3
                call update_mc(nBuilt_in, igenPar)
            enddo
            do igenPar = 5,14
                call update_mc(nBuilt_in, igenPar)
            enddo
            do igenPar = 23,25
                call update_mc(nBuilt_in, igenPar)
            enddo
            call update_mc(nBuilt_in, 28)
            do igenPar = 32,48
                call update_mc(nBuilt_in, igenPar)
            enddo
            do igenPar = 56,57
                call update_mc(nBuilt_in, igenPar)
            enddo
            do igenPar = 64,73
                call update_mc(nBuilt_in, igenPar)
            enddo
            do igenPar = 80,97
                call update_mc(nBuilt_in, igenPar)
            enddo
        else if(i == 11) then
            do igenPar = 26, 49
                call update_mc(nBuilt_in, igenPar)
            enddo
        endif

    end subroutine update_min_max_val_new_new

    subroutine costFuncObs()
        implicit none
        real(8) J_cost, delta_J(12), cs_rand
        integer :: ipft, iaccep(12), i, iupdata, nupdata
        
        J_new = 0

        !! cLeaf + cStem
        if(vars4MCMC%cPlant_tree_y%existOrNot)then
            call CalculateCost(vars4MCMC%cPlant_tree_y%mdData(:,4), vars4MCMC%cPlant_tree_y%obsData(:,4),&
                 vars4MCMC%cPlant_tree_y%obsData(:,5), J_cost)
            J_new(1) = J_new(1) + J_cost*10
        endif

        ! ! ANPP_Tree_y
        if(vars4MCMC%ANPP_Tree_y%existOrNot)then
            call CalculateCost(vars4MCMC%ANPP_Tree_y%mdData(:,4), vars4MCMC%ANPP_Tree_y%obsData(:,4),&
                 vars4MCMC%ANPP_Tree_y%obsData(:,5), J_cost)
            J_new(2) = J_new(2) + J_cost*10
        endif

        ! leaf_mass_shrub_y
        if(vars4MCMC%leaf_mass_shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%leaf_mass_shrub_y%mdData(:,4), vars4MCMC%leaf_mass_shrub_y%obsData(:,4),&
                 vars4MCMC%leaf_mass_shrub_y%obsData(:,5), J_cost)
            J_new(3) = J_new(3) + J_cost*10
        endif

        ! stem_mass_shrub_y
        if(vars4MCMC%stem_mass_shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%stem_mass_shrub_y%mdData(:,4), vars4MCMC%stem_mass_shrub_y%obsData(:,4),&
                 vars4MCMC%stem_mass_shrub_y%obsData(:,5), J_cost)
            J_new(4) = J_new(4) + J_cost*10
        endif

        ! ! ! ANPP_Shrub_y
        if(vars4MCMC%ANPP_Shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%ANPP_Shrub_y%mdData(:,4), vars4MCMC%ANPP_Shrub_y%obsData(:,4),&
                 vars4MCMC%ANPP_Shrub_y%obsData(:,5), J_cost)
            J_new(5) = J_new(5) + J_cost*10 ! test *10 to better constrain
        endif

        ! BNPP_y  ! tree + shrub
        if(vars4MCMC%BNPP_y%existOrNot)then
            call CalculateCost(vars4MCMC%BNPP_y%mdData(:,4), vars4MCMC%BNPP_y%obsData(:,4),&
                 vars4MCMC%BNPP_y%obsData(:,5), J_cost)
            J_new(6) = J_new(6) + J_cost*100
        endif

        ! C plant sphagnum
        if(vars4MCMC%cPlant_sphag_y%existOrNot)then
            call CalculateCost(vars4MCMC%cPlant_sphag_y%mdData(:,4), vars4MCMC%cPlant_sphag_y%obsData(:,4),&
                 vars4MCMC%cPlant_sphag_y%obsData(:,5), J_cost)
            J_new(7) = J_new(7) + J_cost*10
        endif

        ! ! ! NPP_sphag_y
        if(vars4MCMC%NPP_sphag_y%existOrNot)then
            call CalculateCost(vars4MCMC%NPP_sphag_y%mdData(:,4), vars4MCMC%NPP_sphag_y%obsData(:,4),&
                 vars4MCMC%NPP_sphag_y%obsData(:,5), J_cost)
            J_new(8) = J_new(8) + J_cost*10
        endif  

        if(vars4MCMC%gpp_d%existOrNot)then
            call CalculateCost(vars4MCMC%gpp_d%mdData(:,4), vars4MCMC%gpp_d%obsData(:,4),&
                 vars4MCMC%gpp_d%obsData(:,5), J_cost)
            J_new(9) = J_new(9) + J_cost*1000
        endif

        if(vars4MCMC%er_h%existOrNot)then
            call CalculateCost(vars4MCMC%er_h%mdData(:,4), vars4MCMC%er_h%obsData(:,4),&
                 vars4MCMC%er_h%obsData(:,5), J_cost)
            J_new(10) = J_new(10) + J_cost*1000
        endif

        ! nee_h         ! shrub + sphag.
        if(vars4MCMC%nee_h%existOrNot)then
            call CalculateCost(vars4MCMC%nee_h%mdData(:,4), vars4MCMC%nee_h%obsData(:,4),&
                 vars4MCMC%nee_h%obsData(:,5), J_cost)
            J_new(11) = J_new(11) + J_cost*1000
        endif

        if(vars4MCMC%cSoil_y%existOrNot)then
            call CalculateCost(vars4MCMC%cSoil_y%mdData(:,4), vars4MCMC%cSoil_y%obsData(:,4),&
                 vars4MCMC%cSoil_y%obsData(:,5), J_cost)
            J_new(12) = J_new(12) + J_cost
        endif

        ! ch4_h
        ! if(vars4MCMC%ch4_h%existOrNot)then
        !     call CalculateCost(vars4MCMC%ch4_h%mdData(:,4), vars4MCMC%ch4_h%obsData(:,4),&
        !          vars4MCMC%ch4_h%obsData(:,5), J_cost)
        !     J_new(13) = J_new(13) + J_cost
        ! endif
        ! =====================================================================================

        

        ! ! ! ! ! LAI_d         ! tree  + Shrub
        ! if(vars4MCMC%LAI_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%LAI_d%mdData(:,4), vars4MCMC%LAI_d%obsData(:,4),&
        !          vars4MCMC%LAI_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif

        ! ! ! leaf_resp_shrub_d 
        ! ! if(vars4MCMC%leaf_resp_shrub_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%leaf_resp_shrub_d%mdData(:,4), vars4MCMC%leaf_resp_shrub_d%obsData(:,4),&
        ! !          vars4MCMC%leaf_resp_shrub_d%obsData(:,5), J_cost)
        ! !     J_new(2) = J_new(2) + J_cost
        ! ! endif
        ! ! ! print*, "J_new13: ", J_new

        ! ! ! leaf_resp_tree_d 
        ! ! if(vars4MCMC%leaf_resp_tree_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%leaf_resp_tree_d%mdData(:,4), vars4MCMC%leaf_resp_tree_d%obsData(:,4),&
        ! !          vars4MCMC%leaf_resp_tree_d%obsData(:,5), J_cost)
        ! !     J_new(1) = J_new(1) + J_cost
        ! ! endif

         ! ! ! photo_shrub_d 
        ! ! if(vars4MCMC%photo_shrub_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%photo_shrub_d%mdData(:,4), vars4MCMC%photo_shrub_d%obsData(:,4),&
        ! !          vars4MCMC%photo_shrub_d%obsData(:,5), J_cost)
        ! !     J_new(2) = J_new(2) + J_cost
        ! ! endif
        ! ! ! ! print*, "J_new18: ", J_new

        ! ! ! photo_tree_d 
        ! ! if(vars4MCMC%photo_tree_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%photo_tree_d%mdData(:,4), vars4MCMC%photo_tree_d%obsData(:,4),&
        ! !          vars4MCMC%photo_tree_d%obsData(:,5), J_cost)
        ! !     J_new(1) = J_new(1) + J_cost
        ! ! endif

        ! if(vars4MCMC%er_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%er_d%mdData(:,4), vars4MCMC%er_d%obsData(:,4),&
        !          vars4MCMC%er_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif

        ! nee_d         ! Shrub + sphag.
        ! if(vars4MCMC%nee_d%existOrNot)then
        !     call CalculateCost(vars4MCMC%nee_d%mdData(:,4), vars4MCMC%nee_d%obsData(:,4),&
        !          vars4MCMC%nee_d%obsData(:,5), J_cost)
        !     J_new(4) = J_new(4) + J_cost
        ! endif


        ! ! ! ch4_d 
        ! ! if(vars4MCMC%ch4_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%ch4_d%mdData(:,4), vars4MCMC%ch4_d%obsData(:,4),&
        ! !          vars4MCMC%ch4_d%obsData(:,5), J_cost)
        ! !     J_new(5) = J_new(5) + J_cost
        ! ! endif
        ! ! ! print*, "J_new15: ", J_new
        
        ! ! ! print*, "J_new16: ", J_new
        
        ! ! ! CN_shag_d 
        ! ! if(vars4MCMC%CN_shag_d%existOrNot)then
        ! !     call CalculateCost(vars4MCMC%CN_shag_d%mdData(:,4), vars4MCMC%CN_shag_d%obsData(:,4),&
        ! !          vars4MCMC%CN_shag_d%obsData(:,5), J_cost)
        ! !     J_new(3) = J_new(3) + J_cost
        ! ! endif
        ! ! ! print*, "J_new17: ", J_new

       
        ! ! print*, "J_new19: ", J_new
        ! ! ------------------------------------------------------------------------------------
        ! ! write(*,*) "here2",J_new
        iaccep = 0
        nupdata = 13
        if (iDAsimu .eq. 1*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 2*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 3*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 4*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 5*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 6*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 7*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 8*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 9*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 10*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. 11*nDAsimu/nupdata+1) J_last = J_new
        if (iDAsimu .eq. (nupdata-1)*nDAsimu/nupdata+1) J_last = J_new
        
        do i = 1,12
            ! if (iDAsimu .eq. i*nDAsimu/nupdata+1) then
            !    if(i<12) J_last(i+1) = J_new(i+1)
            ! endif
            if(J_new(i) .eq. 0) then ! no data is available
                delta_J(i) = -0.1
            else
                delta_J(i) = J_new(i) - J_last(i)
            endif

            delta_J(i) = delta_J(i)

            call random_number(cs_rand)
            if(AMIN1(1.0, exp(-delta_J(i))) .gt. cs_rand)then
                iaccep(i) = iaccep(i) + 1
            endif
        enddo
        ! if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
        !     ! if (AMIN1(1.0, exp(-delta_J(5))) .gt. cs_rand) then
        !         ! if (sum(iaccep(1:5)) .gt. 4)then
        !             upgraded = upgraded + 1
        !             J_last = J_new
        !         ! endif
        !     ! endif
        ! endif
        ! iupdata = 12
        
        if(iDAsimu < 1*nDAsimu/nupdata) then 
            iupdata = 1
            if(AMIN1(1.0, exp(-sum(delta_J(1:2)))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = sum(J_new(1:2))
            J_show_old = sum(J_last(1:2))
        else if(iDAsimu < 2*nDAsimu/nupdata) then
            iupdata = 2
            if(AMIN1(1.0, exp(-delta_J(2))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(2)
            J_show_old = J_last(2)
        else if(iDAsimu < 3*nDAsimu/nupdata) then
            iupdata = 3
            if(AMIN1(1.0, exp(-delta_J(3))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(3)
            J_show_old = J_last(3)
        else if(iDAsimu < 4*nDAsimu/nupdata) then
            iupdata = 4
            if(AMIN1(1.0, exp(-delta_J(4))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(4)
            J_show_old = J_last(4)
        else if(iDAsimu < 5*nDAsimu/nupdata) then
            iupdata = 5
            if(AMIN1(1.0, exp(-delta_J(5))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(5)
            J_show_old = J_last(5)
        else if(iDAsimu < 6*nDAsimu/nupdata) then
            iupdata = 6
            if(AMIN1(1.0, exp(-delta_J(6))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(6)
            J_show_old = J_last(6)
        else if(iDAsimu < 7*nDAsimu/nupdata) then
            iupdata = 7
            if(AMIN1(1.0, exp(-delta_J(7))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(7)
            J_show_old = J_last(7)
        else if(iDAsimu < 8*nDAsimu/nupdata) then
            iupdata = 8
            if(AMIN1(1.0, exp(-sum(delta_J(7:8)))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = sum(J_new(7:8))
            J_show_old = sum(J_last(7:8))
        else if(iDAsimu < 9*nDAsimu/nupdata) then
            iupdata = 9
            if(AMIN1(1.0, exp(-delta_J(9))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(9)
            J_show_old = J_last(9)
        else if(iDAsimu < 10*nDAsimu/nupdata) then
            iupdata = 10
            if(AMIN1(1.0, exp(-delta_J(10))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(10)
            J_show_old = J_last(10)
        else if(iDAsimu < 11*nDAsimu/nupdata) then
            iupdata = 11
            if(AMIN1(1.0, exp(-sum(delta_J(9:11)))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = sum(J_new(9:11))
            J_show_old = sum(J_last(9:11))
        else if (iDAsimu < 12*nDAsimu/nupdata) then
            iupdata = 12
            if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = J_new(12)
            J_show_old = J_last(12)
        else
            if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
                upgraded = upgraded + 1
                J_last = J_new
            endif
            J_show_new = sum(J_new)
            J_show_old = sum(J_last)
        endif
        ! print*, iDAsimu, nDAsimu, iupdata
        ! if(AMIN1(1.0, exp(-delta_J(iupdata))) .gt. cs_rand)then
        !     upgraded = upgraded + 1
        !     J_last = J_new
        ! endif
            
        
        accept_rate = real(upgraded)/real(iDAsimu)

        do i = 1, nupdata
            if (iDAsimu .eq. i*nDAsimu/nupdata) then
                if (i==1) then 
                    before_upgrade = 1
                else
                    before_upgrade = after_upgrade
                endif
                after_upgrade = upgraded
                call  update_min_max_val_new(i)      
            endif
        enddo
    end subroutine costFuncObs

    subroutine update_min_max_val_new(i)
        implicit none
        integer :: i, igenPar
        integer :: n_chg_pars, nBuilt_in
        real(8) :: new_max, new_min
        n_chg_pars = after_upgrade - before_upgrade
        nBuilt_in  = int(0.7*n_chg_pars)
        ! print*, before_upgrade, after_upgrade, nBuilt_in
        ! print*, shape(arr_params_set%tot_paramsets)
        if (i == 1) then ! update tree
            call update_mc(nBuilt_in, 33)
            call update_mc(nBuilt_in, 34)
            call update_mc(nBuilt_in, 35)

            call update_mc(nBuilt_in, 45)
            call update_mc(nBuilt_in, 46)
            call update_mc(nBuilt_in, 47)
            ! call update_mc(nBuilt_in, 48)
            ! call update_mc(nBuilt_in, 49)
        else if(i == 2) then
            do igenPar = 26, 49
                call update_mc(nBuilt_in, igenPar)
            enddo
        else if (i == 3) then ! cleaf shrub
            call update_mc(nBuilt_in, 57)
            call update_mc(nBuilt_in, 69)
            ! do igenPar = 69,73
            !     call update_mc(nBuilt_in, igenPar)
            ! enddo
        else if (i == 4) then ! cStem shrub
            call update_mc(nBuilt_in, 58)
            call update_mc(nBuilt_in, 70)
            ! do igenPar = 69,73
            !     call update_mc(nBuilt_in, igenPar)
            ! enddo
        else if (i==5) then ! anpp shrub
            do igenPar = 50,73
                call update_mc(nBuilt_in, igenPar)
            enddo
        else if (i==6) then ! BNPP
            do igenPar = 26,73
                call update_mc(nBuilt_in, igenPar)
            enddo
        else if (i==7) then ! C plant sphagnum
            call update_mc(nBuilt_in, 93)
            call update_mc(nBuilt_in, 94)
            call update_mc(nBuilt_in, 95)
            call update_mc(nBuilt_in, 81)
            call update_mc(nBuilt_in, 82)
            call update_mc(nBuilt_in, 83)
            ! do igenPar = 93,97
            !     call update_mc(nBuilt_in, igenPar)
            ! enddo
        else if (i==8) then  ! NPP shagnum
            do igenPar = 74,97
                call update_mc(nBuilt_in, igenPar)
            enddo
        else if (i==9) then  ! gpp
            do igenPar = 50,56
                call update_mc(nBuilt_in, igenPar)
            enddo
            do igenPar = 64,66
                call update_mc(nBuilt_in, igenPar)
            enddo
            ! ----------------
            do igenPar = 74,80
                call update_mc(nBuilt_in, igenPar)
            enddo
            do igenPar = 88,90
                call update_mc(nBuilt_in, igenPar)
            enddo
        else if (i==10) then ! er
            do igenPar = 60,66
                call update_mc(nBuilt_in, igenPar)
            enddo
            do igenPar = 84,90
                call update_mc(nBuilt_in, igenPar)
            enddo
        else if (i==11) then ! nee
            do igenPar = 50,97
                call update_mc(nBuilt_in, igenPar)
            enddo
        else if (i==12) then ! soil
            do igenPar = 1,97
                call update_mc(nBuilt_in, igenPar)
            enddo
        endif

    end subroutine update_min_max_val_new

    subroutine update_mc(nBuilt_in, igenPar)
        implicit none
        integer :: nBuilt_in, igenPar
        real(8) :: new_max, new_min
        new_max = maxval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, igenPar))
        new_min = minval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, igenPar))
        mc_DApar%DAparmax(igenPar) = new_max                  
        mc_DApar%DAparmin(igenPar) = new_min
        if ((mc_DApar%DAparmax(igenPar) - mc_DApar%DAparmin(igenPar))<0.01)then
            mc_DApar%DAparmax(igenPar) = mc_DApar%DAparmin(igenPar)
        endif
    end subroutine update_mc

    subroutine update_min_max_val(i)
        implicit none
        integer :: i, igenPar
        integer :: n_chg_pars, nBuilt_in
        real(8) :: new_max, new_min
        n_chg_pars = after_upgrade - before_upgrade
        nBuilt_in  = int(0.5*n_chg_pars)
        ! print*, before_upgrade, after_upgrade, nBuilt_in
        ! print*, shape(arr_params_set%tot_paramsets)
        if (i < 3) then ! update tree
            do igenPar = 1, mark_npar(2)     ! for each parameters
                new_max = maxval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, mark_npar(1) + igenPar))
                new_min = minval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, mark_npar(1) + igenPar))
                mc_DApar%DAparmax(mark_npar(1) + igenPar) = new_max                  
                mc_DApar%DAparmin(mark_npar(1) + igenPar) = new_min                    
            enddo
            if ((mc_DApar%DAparmax(mark_npar(1) + igenPar) - mc_DApar%DAparmin(mark_npar(1) + igenPar))<0.01)then
                mc_DApar%DAparmax(mark_npar(1) + igenPar) = mc_DApar%DAparmin(mark_npar(1) + igenPar)
            endif
            ! call generate_rand_newpar()
        else if (i < 6) then  ! update shrub
            do igenPar = 1, mark_npar(3)     ! for each parameters
                mc_DApar%DAparmax(sum(mark_npar(1:2)) + igenPar) = &
                    maxval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, sum(mark_npar(1:2)) + igenPar))
                mc_DApar%DAparmin(sum(mark_npar(1:2)) + igenPar) = &
                    minval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, sum(mark_npar(1:2)) + igenPar))
            enddo
            ! call generate_rand_newpar()
        else if (i == 6) then
            do igenPar = 1, mark_npar(2) + mark_npar(3)     ! for each parameters
                mc_DApar%DAparmax(mark_npar(1) + igenPar) = &
                    maxval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, mark_npar(1) + igenPar))
                mc_DApar%DAparmin(mark_npar(1) + igenPar) = &
                    minval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, mark_npar(1) + igenPar))
            enddo
            ! call generate_rand_newpar()
        else if (i < 9) then
            do igenPar = 1, mark_npar(4)     ! for each parameters
                mc_DApar%DAparmax(sum(mark_npar(1:3)) + igenPar) = &
                    maxval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, sum(mark_npar(1:3)) + igenPar))
                mc_DApar%DAparmin(sum(mark_npar(1:3)) + igenPar) = &
                    minval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, sum(mark_npar(1:3)) + igenPar))
            enddo
            ! call generate_rand_newpar()
        else if (i<12) then
            do igenPar = 1, sum(mark_npar(3:4))     ! for each parameters
                mc_DApar%DAparmax(sum(mark_npar(1:2)) + igenPar) = &
                    maxval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, sum(mark_npar(1:2)) + igenPar))
                mc_DApar%DAparmin(sum(mark_npar(1:2)) + igenPar) = &
                    minval(arr_params_set%tot_paramsets(before_upgrade + nBuilt_in:after_upgrade, sum(mark_npar(1:2)) + igenPar))
            enddo
            ! call generate_rand_newpar()
        endif

    end subroutine update_min_max_val

    subroutine generate_rand_newpar()
        implicit none
        ! real(8), intent(in) :: par_old(:), par_min(:), par_max(:)
        ! real(8), intent(inout) :: par_new(:) 
        integer igenPar, parflag, ipft
        real(8) rand_harvest, rand

        call random_seed()
        ! print*, "test heare"
        do igenPar = 1, npar4DA     ! for each parameters
585         continue
            call random_number(rand_harvest)    
            rand = rand_harvest           ! create a random number in [-0.5, 0.5]
            mc_DApar%DApar_old(igenPar) = mc_DApar%DAparmin(igenPar) + &
                rand*(mc_DApar%DAparmax(igenPar) - mc_DApar%DAparmin(igenPar)) * search_scale   ! create new parameter
            print*, "test2: ", igenPar, mc_DApar%DApar(igenPar), mc_DApar%DAparmin(igenPar),mc_DApar%DAparmax(igenPar)
            if((mc_DApar%DApar_old(igenPar) .gt. mc_DApar%DAparmax(igenPar)) &
                &   .or. (mc_DApar%DApar_old(igenPar) .lt. mc_DApar%DAparmin(igenPar))) then 
            ! print*, "here 2: ",  igenPar,mc_DApar%DApar_old(igenPar), mc_DApar%DAparmin(igenPar), mc_DApar%DAparmax(igenPar)
                    ! stop
                goto 585                  ! judge the range of new parameter
            endif
        enddo
        return   ! mainly return the DApar, meanwhile update the coefnorm
    end subroutine generate_rand_newpar

    subroutine CalculateCost(datMod4MCMC, datObs4MCMC, stdObs4MCMC, JCost)
        ! calculate the cost of the observation and simulation, and update the number of updated
        implicit none
        real(8), intent(in) :: datMod4MCMC(:), datObs4MCMC(:), stdObs4MCMC(:)
        integer nLine, iLine, nCost
        real(8) JCost, dObsSimu, std4cal

        nLine = size(datObs4MCMC)
        nCost = 0
        JCost = 0.

        do iLine = 1, nLine
            if(datObs4MCMC(iLine) .gt. -999 .and. datMod4MCMC(iLine) .gt. -999)then
                nCost    = nCost + 1   
                dObsSimu = datMod4MCMC(iLine) - datObs4MCMC(iLine) 
                if (stdObs4MCMC(iLine) < 0.51) then 
                    std4cal = 0.5
                else
                    std4cal = stdObs4MCMC(iLine)
                    std4cal  = abs(0.2*datObs4MCMC(iLine))
                endif
                if (datObs4MCMC(iLine) < 1) then
                    std4cal = 0.5
                else
                    std4cal  = abs(0.2*datObs4MCMC(iLine))
                endif
                ! std4cal  = abs(0.2*datObs4MCMC(iLine))
                if (std4cal <=0) std4cal = 1
                JCost    = JCost + (dObsSimu*dObsSimu)/(2*std4cal*std4cal)
            endif
        enddo
        if(nCost .gt. 0) JCost=JCost/real(nCost)
        return ! JCost
    end subroutine CalculateCost

    subroutine racine_mat(M, Mrac,npara)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Square root of a matrix							  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer npara,i, nrot
        real(8) M(npara,npara),Mrac(npara,npara)
        real(8) valpr(npara),vectpr(npara,npara)
        Mrac=0.
        call jacobi(M,npara,npara,valpr,vectpr,nrot)
        do i=1,npara
            if(valpr(i).ge.0.) then
                Mrac(i,i)=sqrt(valpr(i))
            else
                print*, 'WARNING!!! Square root of the matrix is undefined.'
                print*, ' A negative eigenvalue has been set to zero - results may be wrong'
                Mrac=M
                return
            endif
        enddo
        Mrac=matmul(matmul(vectpr, Mrac),transpose(vectpr))

    end subroutine racine_mat

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Extraction of the eigenvalues and the eigenvectors !!
    !! of a matrix (Numerical Recipes)					  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE jacobi(a,n,np,d,v,nrot)
        INTEGER :: n,np,nrot
        real(8) :: a(np,np),d(np),v(np,np)
        INTEGER, PARAMETER :: NMAX=500
        INTEGER :: i,ip,iq,j
        real(8) :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
        
        do ip=1,n
            do iq=1,n
                v(ip,iq)=0.
            end do
            v(ip,ip)=1.
        end do
        
        do ip=1,n
            b(ip)=a(ip,ip)
            d(ip)=b(ip)
            z(ip)=0.
        end do
        
        nrot=0
        do i=1,50
            sm=0.
            do ip=1,n-1
                do iq=ip+1,n
                    sm=sm+abs(a(ip,iq))
                end do
            end do
            if(sm.eq.0.)return
            if(i.lt.4)then
                tresh=0.2*sm/n**2
            else
                tresh=0.
            endif
            do ip=1,n-1
                do iq=ip+1,n
                    g=100.*abs(a(ip,iq))
                    if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
                        a(ip,iq)=0.
                    else if(abs(a(ip,iq)).gt.tresh)then
                        h=d(iq)-d(ip)
                        if(abs(h)+g.eq.abs(h))then
                            t=a(ip,iq)/h
                        else
                            theta=0.5*h/a(ip,iq)
                            t=1./(abs(theta)+sqrt(1.+theta**2))
                            if(theta.lt.0.) then
                                t=-t
                            endif
                        endif
                        c=1./sqrt(1+t**2)
                        s=t*c
                        tau=s/(1.+c)
                        h=t*a(ip,iq)
                        z(ip)=z(ip)-h
                        z(iq)=z(iq)+h
                        d(ip)=d(ip)-h
                        d(iq)=d(iq)+h
                        a(ip,iq)=0.
                        do j=1,ip-1
                            g=a(j,ip)
                            h=a(j,iq)
                            a(j,ip)=g-s*(h+g*tau)
                            a(j,iq)=h+s*(g-h*tau)
                        end do
                        do j=ip+1,iq-1
                            g=a(ip,j)
                            h=a(j,iq)
                            a(ip,j)=g-s*(h+g*tau)
                            a(j,iq)=h+s*(g-h*tau)
                        end do
                        do j=iq+1,n
                            g=a(ip,j)
                            h=a(iq,j)
                            a(ip,j)=g-s*(h+g*tau)
                            a(iq,j)=h+s*(g-h*tau)
                        end do
                        do j=1,n
                            g=v(j,ip)
                            h=v(j,iq)
                            v(j,ip)=g-s*(h+g*tau)
                            v(j,iq)=h+s*(g-h*tau)
                        end do
                        nrot=nrot+1
                    endif
                end do
            end do
            do ip=1,n
                b(ip)=b(ip)+z(ip)
                d(ip)=b(ip)
                z(ip)=0.
            end do
        end do
        print*, 'too many iterations in jacobi'
        return
    END subroutine jacobi

    subroutine gengaussvect(gamma_racine,xold,xnew,npara)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Generation of a random vector from a multivariate  !!
    !! normal distribution with mean zero and covariance  !!
    !! matrix gamma.									  !!
    !! Beware!!! In order to improve the speed of the	  !!
    !! algorithms, the subroutine use the Square root	  !!
    !! matrix of gamma									  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer npara, i
        real(8) gamma_racine(npara,npara)
        real(8) x(npara),xold(npara),xnew(npara)
        
        do i=1,npara
            x(i)=rangauss(25)
        enddo
        
        x = matmul(gamma_racine, x)
        xnew = xold + x
    end subroutine gengaussvect

   real(8) function rangauss(idum)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Generation of a random number from a standard	  !!
    !! normal distribution. (Numerical Recipes)           !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer idum
        real(8) v1, v2, r, fac, gset
        real(8) r_num
        integer :: iset
        
        ! data iset/0/
        iset = 0
        if(iset==0) then
1           CALL random_number(r_num)
            v1=2.*r_num-1
            CALL random_number(r_num)
            v2=2.*r_num-1
            r=(v1)**2+(v2)**2
            if(r>=1) go to 1
            fac=sqrt(-2.*log(r)/r)
            gset=v1*fac
            rangauss=v2*fac
            iset=1
        else
            rangauss=gset
            iset=0
        end if
        return
    end function

    subroutine varcov(tab,varcovar,npara,ncov)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! variance matrix of a matrix of data				  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer npara,ncov
        real(8) tab(ncov,npara),tab2(ncov,npara)
        real(8) varcovar(npara,npara)
        
        call centre(tab,tab2,npara,ncov)
        
        varcovar = matmul(transpose(tab2), tab2)*(1./real(ncov))
        
    end subroutine varcov



    subroutine centre(mat,mat_out,npara,ncov)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Compute the centered matrix, ie. the matrix minus  !!
    !! the column means									  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer npara,i,ncov
        real(8) mat(ncov,npara),mat_out(ncov,npara)
        ! real(8) mean

        do i=1,npara
            mat_out(:,i) = mat(:,i) - mean(mat(:,i),ncov)
        enddo

    end subroutine centre

    real(8) function mean(tab,ncov)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! mean of a vector									  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer ncov, incov
        real(8) tab(ncov)
        real(8) mean_tt
        mean_tt=0.
        do incov=1,ncov
        mean_tt=mean_tt+tab(incov)/real(ncov)
        enddo
        mean=mean_tt
    End Function mean

    ! real(8) function mean(tab, ncov)
    !     integer ncov
    !     real(8) tab(ncov)
    ! end function mean
    
    subroutine deallocate_mcmc()
    ! deallocate some variables and summary the information of MCMC
        implicit none
        integer :: ipft
        
        

        ! if(allocated(mc_parvals)) then
        !     npft = size(mc_parvals)
        !     do ipft = 1, npft
        !         if(allocated(mc_parvals(ipft)%parval)) deallocate(mc_parvals(ipft)%parval)
        !         if(allocated(mc_parvals(ipft)%parmin)) deallocate(mc_parvals(ipft)%parmin)
        !         if(allocated(mc_parvals(ipft)%parmax)) deallocate(mc_parvals(ipft)%parmax)
        !     enddo
        !     deallocate(mc_parvals)
        ! endif

        ! if(allocated(MDparval)) deallocate(MDparval)
        ! if(allocated(mc_DApar)) then
            ! npft = size(mc_DApar)
            ! do ipft = 1, npft
                if(allocated(mc_DApar%DAparmin))  deallocate(mc_DApar%DAparmin)
                if(allocated(mc_DApar%DAparmax))  deallocate(mc_DApar%DAparmax)
                if(allocated(mc_DApar%DAparidx_st))  deallocate(mc_DApar%DAparidx_st)
                if(allocated(mc_DApar%DAparidx_sp)) deallocate(mc_DApar%DAparidx_sp)
                if(allocated(mc_DApar%DApar))     deallocate(mc_DApar%DApar)
                if(allocated(mc_DApar%DApar_old)) deallocate(mc_DApar%DApar_old)

                if(allocated(mc_DApar%coefhistory)) deallocate(mc_DApar%coefhistory)
                if(allocated(mc_DApar%coefnorm))    deallocate(mc_DApar%coefnorm)
                if(allocated(mc_DApar%coefac))      deallocate(mc_DApar%coefac)
        
                if(allocated(mc_DApar%gamnew))   deallocate(mc_DApar%gamnew)
            ! enddo
            ! deallocate(mc_DApar)
        ! endif
        

        if(allocated(vars4MCMC%ANPP_Shrub_y%obsData))  deallocate(vars4MCMC%ANPP_Shrub_y%obsData)
        if(allocated(vars4MCMC%ANPP_Tree_y%obsData))  deallocate(vars4MCMC%ANPP_Tree_y%obsData)
        if(allocated(vars4MCMC%NPP_sphag_y%obsData))  deallocate(vars4MCMC%NPP_sphag_y%obsData)
        if(allocated(vars4MCMC%BNPP_y%obsData))  deallocate(vars4MCMC%BNPP_y%obsData)        ! tree + shrub
        if(allocated(vars4MCMC%er_d%obsData))  deallocate(vars4MCMC%er_d%obsData)          ! shrub + sphag.
        if(allocated(vars4MCMC%er_h%obsData))  deallocate(vars4MCMC%er_h%obsData)          ! shrub + sphag.
        if(allocated(vars4MCMC%gpp_d%obsData))  deallocate(vars4MCMC%gpp_d%obsData)         ! Shrub + sphag.
        if(allocated(vars4MCMC%nee_d%obsData))  deallocate(vars4MCMC%nee_d%obsData)         ! Shrub + sphag.
        if(allocated(vars4MCMC%nee_h%obsData))  deallocate(vars4MCMC%nee_h%obsData)         ! shrub + sphag.
        if(allocated(vars4MCMC%LAI_d%obsData))  deallocate(vars4MCMC%LAI_d%obsData)         ! tree  + Shrub
        !
        if(allocated(vars4MCMC%leaf_mass_shrub_y%obsData))  deallocate(vars4MCMC%leaf_mass_shrub_y%obsData)
        if(allocated(vars4MCMC%stem_mass_shrub_y%obsData))  deallocate(vars4MCMC%stem_mass_shrub_y%obsData)
        if(allocated(vars4MCMC%leaf_resp_shrub_d%obsData))  deallocate(vars4MCMC%leaf_resp_shrub_d%obsData) 
        if(allocated(vars4MCMC%leaf_resp_tree_d%obsData))  deallocate(vars4MCMC%leaf_resp_tree_d%obsData) 
        ! methane
        if(allocated(vars4MCMC%ch4_d%obsData))  deallocate(vars4MCMC%ch4_d%obsData) 
        if(allocated(vars4MCMC%ch4_h%obsData))  deallocate(vars4MCMC%ch4_h%obsData) 
        ! 
        if(allocated(vars4MCMC%CN_shag_d%obsData))  deallocate(vars4MCMC%CN_shag_d%obsData) 
        if(allocated(vars4MCMC%photo_shrub_d%obsData))  deallocate(vars4MCMC%photo_shrub_d%obsData) 
        if(allocated(vars4MCMC%photo_tree_d%obsData))  deallocate(vars4MCMC%photo_tree_d%obsData) 
        ! ---------------------------------------------------------------------------------------------

        if(allocated(vars4MCMC%ANPP_Shrub_y%mdData))  deallocate(vars4MCMC%ANPP_Shrub_y%mdData)
        if(allocated(vars4MCMC%ANPP_Tree_y%mdData))  deallocate(vars4MCMC%ANPP_Tree_y%mdData)
        if(allocated(vars4MCMC%NPP_sphag_y%mdData))  deallocate(vars4MCMC%NPP_sphag_y%mdData)
        if(allocated(vars4MCMC%BNPP_y%mdData))  deallocate(vars4MCMC%BNPP_y%mdData)        ! tree + shrub
        if(allocated(vars4MCMC%er_d%mdData))  deallocate(vars4MCMC%er_d%mdData)          ! shrub + sphag.
        if(allocated(vars4MCMC%er_h%mdData))  deallocate(vars4MCMC%er_h%mdData)          ! shrub + sphag.
        if(allocated(vars4MCMC%gpp_d%mdData))  deallocate(vars4MCMC%gpp_d%mdData)         ! Shrub + sphag.
        if(allocated(vars4MCMC%nee_d%mdData))  deallocate(vars4MCMC%nee_d%mdData)         ! Shrub + sphag.
        if(allocated(vars4MCMC%nee_h%mdData))  deallocate(vars4MCMC%nee_h%mdData)         ! shrub + sphag.
        if(allocated(vars4MCMC%LAI_d%mdData))  deallocate(vars4MCMC%LAI_d%mdData)         ! tree  + Shrub
        !
        if(allocated(vars4MCMC%leaf_mass_shrub_y%mdData))  deallocate(vars4MCMC%leaf_mass_shrub_y%mdData)
        if(allocated(vars4MCMC%stem_mass_shrub_y%mdData))  deallocate(vars4MCMC%stem_mass_shrub_y%mdData)
        if(allocated(vars4MCMC%leaf_resp_shrub_d%mdData))  deallocate(vars4MCMC%leaf_resp_shrub_d%mdData) 
        if(allocated(vars4MCMC%leaf_resp_tree_d%mdData))  deallocate(vars4MCMC%leaf_resp_tree_d%mdData) 
        ! methane
        if(allocated(vars4MCMC%ch4_d%mdData))  deallocate(vars4MCMC%ch4_d%mdData) 
        if(allocated(vars4MCMC%ch4_h%mdData))  deallocate(vars4MCMC%ch4_h%mdData) 
        ! 
        if(allocated(vars4MCMC%CN_shag_d%mdData))  deallocate(vars4MCMC%CN_shag_d%mdData) 
        if(allocated(vars4MCMC%photo_shrub_d%mdData))  deallocate(vars4MCMC%photo_shrub_d%mdData) 
        if(allocated(vars4MCMC%photo_tree_d%mdData))  deallocate(vars4MCMC%photo_tree_d%mdData) 
        ! ---------------------------------------------------------------------------------------------

        if(allocated(parnames)) deallocate(parnames)
        ! if(allocated(npar4DA)) deallocate(npar4DA)
        ! if(allocated(fact_rejet)) deallocate(fact_rejet)

        ! in MCMC_outputs module
        ! do ipft = 1, npft
            if(allocated(arr_params_set%tot_paramsets)) deallocate(arr_params_set%tot_paramsets)
            if(allocated(arr_params_set%sel_paramsets)) deallocate(arr_params_set%sel_paramsets)
            if(allocated(arr_params_set%upg_paramsets)) deallocate(arr_params_set%upg_paramsets)
        ! enddo
        ! if(allocated(arr_params_set)) deallocate(arr_params_set)

        ! if (do_mc_out_hr)then
        !     call deallocate_mcmc_outs_type(sel_paramsets_outs_h)
        !     call deallocate_mcmc_outs_type(tot_paramsets_outs_h)
        ! endif
        if (do_mc_out_day)then
            call deallocate_mcmc_outs_type(sel_paramsets_outs_d)
            call deallocate_mcmc_outs_type(tot_paramsets_outs_d)
        endif
        ! if (do_mc_out_mon)then
        !     call deallocate_mcmc_outs_type(sel_paramsets_outs_m)
        !     call deallocate_mcmc_outs_type(tot_paramsets_outs_m)
        ! endif
        if(allocated(mark_npar)) deallocate(mark_npar)
    end subroutine deallocate_mcmc
end module mcmc
