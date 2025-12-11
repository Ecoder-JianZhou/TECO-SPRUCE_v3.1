module mcmc
    use datatypes
    use driver ! run teco simulation
    use mcmc_mod

    integer npar4DA ! number of parameters to be optimize
    integer upgraded, iDAsimu

    logical :: do_cov2create_new_params ! if do cov to create new parameters
    real(8), allocatable :: DAparVal(:)
    real(8), allocatable :: DAparMin(:)
    real(8), allocatable :: DAparMax(:)
    real(8), allocatable :: DAparOldVal(:)
    integer, allocatable :: idxStPar(:), idxFixStep(:)
    ! ------- divide different parameters for different observations ---------
    integer, allocatable :: idxSpecialPars_st(:), idxSpecialPars_sp(:,:)

    type(index_species_parameters_to_opt) :: idxSpPar(numPFT) ! type defined in MCMC_mod
    real(8), allocatable ::  search_scale(:) ! fixed step for parameters

    real(8) J_last(numObsFiles), J_new(numObsFiles)
    ! save parameters in MCMC processes
    real(8), allocatable :: totMcParamsets(:,:) ! upgraded, nParams
    ! mark paramters or observation after 100 cyles no paramters updated
    integer, allocatable :: mark_idx_obs(:) ! for each observation file, 1: updated, 0: not updated
    integer :: mark_simu, isave
    real(8) :: scale_sel_params
    real(8), allocatable :: mark_obs(:,:)

    integer, parameter :: istep_1 = 501,  istep_2 = 1001, istep_3 = 1501 ! all, tree, shrub
    integer, parameter :: istep_4 = 2001, istep_5 = 2501, istep_6 = 5501 ! sphag, soil, annual ch4+Rh
    integer, parameter :: istep_7 = 10501,istep_8 = 15501, istep_9 = 25501 ! annual Tree, Shrub, Sphag
    integer, parameter :: istep_10 = 50000 ! all annual parameters
    ! integer, parameter :: istep_1 = 11,  istep_2 = 21, istep_3 = 31
    ! integer, parameter :: istep_4 = 41, istep_5 = 51, istep_6 = 61

    ! integer, parameter :: istep_1 = 11,  istep_2 = 21, istep_3 = 31 ! all, tree, shrub
    ! integer, parameter :: istep_4 = 41, istep_5 = 51, istep_6 = 61 ! sphag, soil, annual ch4+Rh
    ! integer, parameter :: istep_7 = 71,istep_8 = 81, istep_9 = 91 ! annual Tree, Shrub, Sphag
    ! integer, parameter :: istep_10 = 100 ! all annual parameters

    contains

    subroutine run_mcmc(st)
        implicit none
        type(site_data_type), intent(inout) :: st
        integer :: temp_upgraded, i
        real(8) :: upgraded_rate, sum_j_new, sum_j_last, temp_new_parval(npar4DA)

        print*, "# Start to run MCMC ... "
        call random_seed()
        ! 1. generate new parameters
        ! call generate_new_parameters()
        allocate(mark_obs(100, numObsFiles))
        upgraded = 0
        mark_simu = 0
        scale_sel_params = 0.1
        ! 2. do cylces
        do iDAsimu = 1, 50000!mcset%nDAsimu
            upgraded_rate = (real(upgraded, 8) / real(iDAsimu, 8))!*100

            if(iDAsimu<istep_1)then
                sum_j_new = sum(J_new)
                sum_j_last = sum(J_last)
            elseif(iDAsimu < istep_2)then
                sum_j_new = sum(J_new(9:13))
                sum_j_last = sum(J_last(9:13))
            elseif(iDAsimu<istep_3)then
                sum_j_new = sum(J_new(2:6)) + J_new(20)
                sum_j_last = sum(J_last(2:6)) + J_last(20)
            elseif(iDAsimu<istep_4)then
                sum_j_new = sum(J_new(7:8))
                sum_j_last = sum(J_last(7:8))
            elseif(iDAsimu<istep_5)then
                sum_j_new = sum(J_new(14:15)) + J_new(19)
                sum_j_last = sum(J_last(14:15)) + J_last(19)
            elseif(iDAsimu < istep_6)then
                sum_j_new = J_new(15) + J_new(19)
                sum_j_last = J_last(15) + J_last(19)
            elseif(iDAsimu < istep_7)then
                sum_j_new = sum(J_new(9:13))
                sum_j_last = sum(J_last(9:13))
            elseif(iDAsimu < istep_8)then
                sum_j_new = sum(J_new(2:6)) + J_new(20)
                sum_j_last = sum(J_last(2:6)) + J_last(20)
            elseif(iDAsimu < istep_9)then
                sum_j_new = sum(J_new(7:8))
                sum_j_last = sum(J_last(7:8))
            elseif(iDAsimu < istep_10)then
                sum_j_new = sum(J_new)
                sum_j_last = sum(J_last)
            ! elseif(iDAsimu<istep_6)then
            !     sum_j_new = sum(J_new(2:4))+J_new(7)+sum(J_new(9:12))+sum(J_new(16:20)) 
            !     sum_j_last = sum(J_last(2:4))+J_last(7)+sum(J_last(9:12))+sum(J_last(16:20))
            else
                sum_j_new = sum(J_new)
                sum_j_last = sum(J_last)
            endif
            ! sum_j_new = sum(J_new(9:13))
            ! sum_j_last = sum(J_last(9:13))

            write(*, '(I6, A1, I6, 2X, F12.2, 2X, F12.2, 2X, F12.2, 2X, F12.2, 2X, I6, 2X, F6.2, 2X, F8.6, 2X, F8.6)') &
            iDAsimu, '/', mcset%nDAsimu, sum_j_new, sum_j_last,sum(J_new), sum(J_last), upgraded, upgraded_rate, search_scale(1), &
            search_scale(idxFixStep(1))
                ! initilize
            do i = 1, numObsFiles
                mcVarData(i)%mark_idx = 1
            end do

            ! update mcparams
            if(iDAsimu > 1) call generate_new_parameters()
            temp_new_parval = DAparVal
            DAparVal = DAparOldVal

            if(iDAsimu < istep_1)then ! all parameters without annual DA
                do_annual_mcmc = .False.
                DAparVal(1:25)    = temp_new_parval(1:25)
                DAparVal(42:66)   = temp_new_parval(42:66)
                DAparVal(107:131) = temp_new_parval(107:131)
                DAparVal(172:196) = temp_new_parval(172:196)
            elseif(iDAsimu < istep_2)then ! tree parameters without annual DA
                do_annual_mcmc = .False.
                ! DAparVal(1:25)    = temp_new_parval(1:25)
                DAparVal(42:66)   = temp_new_parval(42:66)
                ! DAparVal(107:131) = temp_new_parval(107:131)
                ! DAparVal(172:196) = temp_new_parval(172:196)
            elseif(iDAsimu < istep_3)then ! Shrub parameters without annual DA
                do_annual_mcmc = .False.
                ! DAparVal(1:25)    = temp_new_parval(1:25)
                ! DAparVal(42:66)   = temp_new_parval(42:66)
                DAparVal(107:131) = temp_new_parval(107:131)
                ! DAparVal(172:196) = temp_new_parval(172:196)
            elseif(iDAsimu < istep_4)then ! Sphag. parameters without annual DA
                do_annual_mcmc = .False.
                ! DAparVal(1:25)    = temp_new_parval(1:25)
                ! DAparVal(42:66)   = temp_new_parval(42:66)
                ! DAparVal(107:131) = temp_new_parval(107:131)
                DAparVal(172:196) = temp_new_parval(172:196)
            elseif(iDAsimu < istep_5)then ! soil parameters without annual DA
                do_annual_mcmc = .False.
                DAparVal(1:25)    = temp_new_parval(1:25)
                ! DAparVal(42:66)   = temp_new_parval(42:66)
                ! DAparVal(107:131) = temp_new_parval(107:131)
                ! DAparVal(172:196) = temp_new_parval(172:196)
            elseif(iDAsimu < istep_6)then ! start annual DA
                do_annual_mcmc = .True.
                DAparVal(26:41) = temp_new_parval(26:41)
            elseif(iDAsimu < istep_7)then
                do_annual_mcmc = .True.
                DAparVal(59:60)    = temp_new_parval(59:60)
                DAparVal(67:106)   = temp_new_parval(67:106)
            elseif(iDAsimu < istep_8)then
                do_annual_mcmc = .True.
                DAparVal(124:125) = temp_new_parval(124:125)
                DAparVal(132:171) = temp_new_parval(132:171)
            elseif(iDAsimu < istep_9)then
                do_annual_mcmc = .True.
                DAparVal(189:190) = temp_new_parval(189:190)
                DAparVal(197:236) = temp_new_parval(197:236)
            elseif(iDAsimu < istep_10)then ! all parameters with annual DA
                do_annual_mcmc = .True.
                DAparVal(26:41)   = temp_new_parval(26:41)
                
                DAparVal(59:60)   = temp_new_parval(59:60)
                DAparVal(67:106)  = temp_new_parval(67:106)

                DAparVal(124:125) = temp_new_parval(124:125)
                DAparVal(132:171) = temp_new_parval(132:171)

                DAparVal(189:190) = temp_new_parval(189:190)
                DAparVal(197:236) = temp_new_parval(197:236)
            else
                do_annual_mcmc = .True.
                DAparVal = temp_new_parval
            endif

            if(iDAsimu == istep_5-1) then   
                DAparOldVal(26:33)   = DAparOldVal(6)  ! Q10rh
                DAparOldVal(34:41)   = DAparOldVal(16) ! Q10pro             
                DAparOldVal(67:74)   = DAparOldVal(57) ! Tree_Entrpy
                DAparOldVal(75:90)   = 1.0             ! Tree s_Vea s_ved
                DAparOldVal(91:98)   = DAparOldVal(48) ! Tree Vcmax0
                DAparOldVal(99:106)  = DAparOldVal(52) ! Tree Q10
                DAparOldVal(132:139) = DAparOldVal(122) ! Shrub_Entrpy
                DAparOldVal(140:155) = 1.0             ! Shrub s_Vea s_ved
                DAparOldVal(156:163) = DAparOldVal(113) ! Shrub Vcmax0
                DAparOldVal(164:171) = DAparOldVal(117) ! Shrub Q10
                DAparOldVal(197:204) = DAparOldVal(187) ! Sphag_Entrpy
                DAparOldVal(205:220) = 1.0             ! Sphag s_Vea s_ved
                DAparOldVal(221:228) = DAparOldVal(178) ! Sphag Vcmax0
                DAparOldVal(229:236) = DAparOldVal(182) ! Sphag Q10
            endif


            call update_mcParams(DAparVal, idxStPar, idxSpPar)
            call update_simuParams() ! update parameters to simulate

            call initialize_teco(st)
            call teco_simu(st, .False.)

            temp_upgraded = upgraded
            call cal_cost_function() ! check if update upgraded
            ! if(iDAsimu < 101) then
            !     mark_obs(iDAsimu, :) = J_new
            ! else
            !     mark_obs(1:99, :) = mark_obs(2:100, :)
            !     mark_obs(100, :) = J_new
            ! endif

            ! mcset%search_scale = mcset%search_scale * (1.0 + 0.05*(upgraded_rate - 0.25))
            ! mcset%search_scale = max(0.0000001, min(0.5, mcset%search_scale))
            
            ! do i = 1, size(idxFixStep)
            !     search_scale(idxFixStep(i)) = scale_sel_params
            ! enddo
            if(iDAsimu == 1)       search_scale = mcset%search_scale
            if(iDAsimu == istep_1) search_scale = mcset%search_scale
            if(iDAsimu == istep_2) search_scale = mcset%search_scale
            if(iDAsimu == istep_3) search_scale = mcset%search_scale
            if(iDAsimu == istep_4) search_scale = mcset%search_scale
            if(iDAsimu == istep_5) search_scale = mcset%search_scale
            if(iDAsimu == istep_6) search_scale = mcset%search_scale
            if(iDAsimu == istep_7) search_scale = mcset%search_scale
            if(iDAsimu == istep_8) search_scale = mcset%search_scale
            if(iDAsimu == istep_9) search_scale = mcset%search_scale
            if(iDAsimu == istep_10) search_scale = mcset%search_scale

            search_scale = search_scale * (1.0 + 0.01*(upgraded_rate - 0.25))
            search_scale = max(0.0000001, min(0.5, search_scale))
            
            if(upgraded > temp_upgraded)then
                ! DAparOldVal = DAparVal
                ! totMcParamsets(upgraded,:) = DAparVal
                totMcParamsets(upgraded,:) = DAparOldVal
                ! mark_simu = 0
            endif

            
            if(mod(iDAsimu, 5000) == 0) then 
                isave = 0
                call mcmc_param_outputs(upgraded, st)
            endif

            if(upgraded > 5)then
                if(iDAsimu == istep_5-1) then
                    isave = 1
                    call mcmc_param_outputs(upgraded, st)
                elseif(iDAsimu == istep_6-1) then
                    isave = 2
                    call mcmc_param_outputs(upgraded, st)
                elseif(iDAsimu == istep_7-1) then
                    isave = 3
                    call mcmc_param_outputs(upgraded, st)
                elseif(iDAsimu == istep_8-1) then
                    isave = 4
                    call mcmc_param_outputs(upgraded, st)
                elseif(iDAsimu == istep_9-1) then
                    isave = 5
                    call mcmc_param_outputs(upgraded, st)
                endif
            endif


            ! mark_simu = mark_simu + 1
            
            ! if(iDAsimu > 5000) then
            !     if(mark_simu > 100) then
            !         scale_sel_params = scale_sel_params * 0.998
            !     else
            !         scale_sel_params = scale_sel_params * 1.002
            !     endif
            ! endif
            ! scale_sel_params = max(min(scale_sel_params, 0.2), 0.001)
            
        enddo
        call mcmc_param_outputs(upgraded, st)
        ! 3. outputs 
    end subroutine run_mcmc

    subroutine init_mcmc(in_mcConfNmlFile)
        implicit none
        character(*), intent(in) :: in_mcConfNmlFile
        integer :: i

        call read_mcmc_conf_nml(in_mcConfNmlFile)
        ! get paramters to optimize
        call filter_params_to_optimize()
        ! print*, "idxFixStep: ", idxFixStep
        allocate(search_scale(npar4DA))
        search_scale = mcset%search_scale
        do i = 1, size(idxFixStep)
            search_scale(idxFixStep(i)) = 0.02
        enddo

        J_last = 9000000.0
        ! initilize outputs
        allocate(totMcParamsets(mcset%nDAsimu, npar4DA))
    end subroutine init_mcmc

    subroutine filter_params_to_optimize()
        implicit none
        real(8), dimension(numAllParams) :: temp_parmin, temp_parmax, temp_parval
        integer :: temp_idx_st(numStParams)
        integer :: i, j, nSpPar4DA(numPFT), nStPar4DA, iFixStep
        type(index_species_parameters_to_opt) :: temp_idx_sp(numPFT)
        ! setting different parameters for different observations
        integer :: temp_sp_params(npft, numSpParams), idx_mark, idxPar_st

        character(50) :: arr_keyParams(15)
        allocate(idxFixStep(24))! 6 for site parameters, 9 * 2 for species parameters

        ! arr_keyParams=["Tau_F       ", "Tau_C       ", "Tau_Micro   ", "Tau_SlowSOM ", &
        !                "Tau_Passive ", "Q10rh       ", "Q10pro      ", "Vmaxfraction", &
        !                "SLAx        ", "GLmax       ", "GRmax       ", "Gsmax       ", &
        !                "Vcmax0      ", "Tau_Leaf    ", "Q10         ", "Rl0         ", &
        !                "JV          ", "fn2r        "]

        arr_keyParams=["Tau_F       ", "Tau_C       ", "Tau_Micro   ", &
                       "Q10rh       ", "Q10pro      ", "r_me        ",  &
                       "SLAx        ", "GLmax       ", "GRmax       ", "Gsmax       ", &
                       "Vcmax0      ", "Q10         ", "Rl0         ", &
                       "JV          ", "fn2r        "]!, "s_cLeaf     ", "s_cStem     "]

        npar4DA = 0
        idxPar_st  = 0
        iFixStep = 0
        print*, "# check all parameters to be optimized ..."
        do i = 1, numStParams
            if(mcParams%st(i)%parmax > mcParams%st(i)%parmin) then
                npar4DA = npar4DA + 1
                temp_parmin(npar4DA) = mcParams%st(i)%parmin
                temp_parmax(npar4DA) = mcParams%st(i)%parmax
                temp_parval(npar4DA) = mcParams%st(i)%parVal
                temp_idx_st(npar4DA) = i
                print*, npar4DA, trim(mcParams%st(i)%parName), mcParams%st(i)%parVal, &
                        mcParams%st(i)%parmin, mcParams%st(i)%parmax
                if(any(mcParams%st(i)%parName == arr_keyParams))then
                    iFixStep = iFixStep + 1
                    idxFixStep(iFixStep) = npar4DA
                    print*, iFixStep,"fixed step for ", trim(mcParams%st(i)%parName), " at ", npar4DA
                endif
            endif
        enddo

        nStPar4DA = npar4DA
        ! site-based parameters
        allocate(idxSpecialPars_st(nStPar4DA))
        idxPar_st = idxPar_st + 1
        idxSpecialPars_st(idxPar_st) = npar4DA
        do i = 1, nStPar4DA
            idxSpecialPars_st(i) = i
        enddo

        do i = 1, numPFT
            nSpPar4DA(i) = 0
            allocate(temp_idx_sp(i)%idx(numSpParams))
            idx_mark = 0
            do j = 1, numSpParams
                if(mcParams%sp(i)%var(j)%parmax > mcParams%sp(i)%var(j)%parmin)then
                    npar4DA = npar4DA + 1
                    nSpPar4DA(i) = nSpPar4DA(i) + 1
                    temp_parmin(npar4DA) = mcParams%sp(i)%var(j)%parmin
                    temp_parmax(npar4DA) = mcParams%sp(i)%var(j)%parmax
                    temp_parval(npar4DA) = mcParams%sp(i)%var(j)%parVal
                    temp_idx_sp(i)%idx(nSpPar4DA(i)) = j
                    print*, npar4DA, trim(mcParams%sp(i)%var(j)%parName), mcParams%sp(i)%var(j)%parVal, &
                        mcParams%sp(i)%var(j)%parmin, mcParams%sp(i)%var(j)%parmax
                    
                    idx_mark = idx_mark + 1
                    temp_sp_params(i, idx_mark) = npar4DA
                    if(i>1)then ! except tree
                        if(any(mcParams%sp(i)%var(j)%parName == arr_keyParams))then
                            iFixStep = iFixStep + 1
                            idxFixStep(iFixStep) = npar4DA
                            print*, iFixStep, "fixed step for ", trim(mcParams%sp(i)%var(j)%parName), " at ", npar4DA
                        endif
                    endif
                endif
            enddo
            allocate(idxSpPar(i)%idx(nSpPar4DA(i)))
            idxSpPar(i)%idx = temp_idx_sp(i)%idx(:nSpPar4DA(i))
            if(allocated(temp_idx_sp(i)%idx)) deallocate(temp_idx_sp(i)%idx)
        enddo

        allocate(idxSpecialPars_sp(3, nSpPar4DA(1)))
        idxSpecialPars_sp = temp_sp_params(:,:nSpPar4DA(1))

        allocate(DAparMax(npar4DA), DAparMin(npar4DA), DAparVal(npar4DA), DAparOldVal(npar4DA), idxStPar(nStPar4DA))

        DAparMin = temp_parmin(:npar4DA)
        DAparMax = temp_parmax(:npar4DA)
        DAparVal = temp_parval(:npar4DA)
        idxStPar = temp_idx_st(:nStPar4DA)
        DAparOldVal = DAparVal

    end subroutine filter_params_to_optimize

    subroutine generate_new_parameters()
        implicit NONE
        integer :: i
        real(8) :: rand_harvest, rand

        
        if (do_cov2create_new_params) then
            
        else
        ! do normal update parameters
            do i = 1, npar4DA
                do
                    call random_number(rand_harvest)
                    rand = rand_harvest - 0.5 
                    ! Ensure parameter is within bounds
                    if(DAparOldVal(i) > DAparMax(i) .or. DAparOldVal(i) < DAparMin(i))then
                        DAparOldVal(i) = DAparMin(i) + rand_harvest*(DAparMax(i)-DAparMin(i))
                    endif
                    if(DAparMin(i) .eq. DAparMax(i)) then
                        DAparOldVal(i) = DAparMax(i)
                    endif
                    ! generate new parameters
                    DAparVal(i) = DAparOldVal(i) + rand*(DAparMax(i)-DAparMin(i))*search_scale(i)
                    ! if out of bounds, retry with a new random number
                    if (DAparVal(i) >= DAparMin(i) .and. DAparVal(i) <= DAparMax(i)) then
                        exit
                    else
                        DAparVal(i) = DAparOldVal(i)
                    endif
                enddo
            end do
        endif

    end subroutine

    subroutine cal_cost_function()
        implicit none
        integer i, mark_upgrade
        real(8) :: JCost, delta_J(numObsFiles), cs_rand
        ! real(8) :: delta_J_norm(numObsFiles),mean_delta_J,std_delta_J
        real(8) :: mean_obs, std_obs

        J_new = 0.
        do i = 1, numObsFiles

            if(mcVarData(i)%existOrNot)then
                call calculate_cost(mcVarData(i)%obsData(:,4), mcVarData(i)%mdData(:,4),JCost)
                if(i == 2 .or. i == 5 .or. i==7)then
                    J_new(i) = J_new(i) + JCost*10000
                else
                    J_new(i) = J_new(i) + JCost*1000!*obsWt(i)
                endif
            end if
            if(mod(iDAsimu, 10000) == 0) then 
                print*, "test_", i,": ", J_new(i) - J_last(i), J_new(i)
                ! print*, "test_",i,": ",J_new(i)
            endif
        enddo

        ! print*, mcVarData(20)%mdData(2,:),mcVarData(20)%obsData(2,:)

        delta_J = J_new - J_last
        mark_upgrade = 0

        call random_number(cs_rand)

        if(iDAsimu < istep_1)then
            if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
                DAparOldVal(1:25)    = DAparVal(1:25)
                DAparOldVal(42:66)   = DAparVal(42:66)
                DAparOldVal(107:131) = DAparVal(107:131)
                DAparOldVal(172:196) = DAparVal(172:196) 
                mark_upgrade = 1
                J_last = J_new
            endif
        elseif(iDAsimu < istep_2)then
            if(AMIN1(1.0, exp(-sum(delta_J(9:13)))) .gt. cs_rand)then
                DAparOldVal(42:66)   = DAparVal(42:66)            
                mark_upgrade = 1
                J_last = J_new
            endif
        elseif(iDAsimu < istep_3)then
            if(AMIN1(1.0, exp(-(sum(delta_J(2:6))+delta_J(20)))) .gt. cs_rand)then
                DAparOldVal(107:131) = DAparVal(107:131)
                mark_upgrade = 1
                J_last = J_new
            endif
        elseif(iDAsimu < istep_4)then
            if(AMIN1(1.0, exp(-(sum(delta_J(7:8))))) .gt. cs_rand)then
                DAparOldVal(172:196) = DAparVal(172:196)
                mark_upgrade = 1
                J_last = J_new
            endif
        elseif(iDAsimu < istep_5)then
            if(AMIN1(1.0, exp(-(sum(delta_J(14:15))+delta_J(19)))) .gt. cs_rand)then
                DAparOldVal(1:25)   = DAparVal(1:25)
                ! DAparOldVal(1:41)   = DAparVal(1:41)
                mark_upgrade = 1
                J_last = J_new
            endif
        elseif(iDAsimu < istep_6)then
            if(AMIN1(1.0, exp(-(delta_J(15)+delta_J(19)))) .gt. cs_rand)then
                ! DAparOldVal(1:25)   = DAparVal(1:25)
                DAparOldVal(26:41)   = DAparVal(26:41)
                mark_upgrade = 1
                J_last = J_new
            endif
        elseif(iDAsimu < istep_7)then
            if(AMIN1(1.0, exp(-sum(delta_J(9:13)))) .gt. cs_rand)then
                DAparOldVal(59:60)   = DAparVal(59:60)
                DAparOldVal(67:106)  = DAparVal(67:106)
                mark_upgrade = 1
                J_last = J_new
            endif
        elseif(iDAsimu < istep_8)then
            ! if(AMIN1(1.0, exp(-sum(delta_J(2:6))+delta_J(20))) .gt. cs_rand)then
            if(AMIN1(1.0, exp(-(sum(delta_J(2:4))+delta_J(7)+sum(delta_J(9:12))+sum(delta_J(15:20))))) .gt. cs_rand)then
                DAparOldVal(124:125) = DAparVal(124:125)
                DAparOldVal(132:171) = DAparVal(132:171)
                mark_upgrade = 1
                J_last = J_new
            endif
        elseif(iDAsimu < istep_9)then
            if(AMIN1(1.0, exp(-(sum(delta_J(2:4))+delta_J(7)+sum(delta_J(9:12))+sum(delta_J(15:20))))) .gt. cs_rand)then
                DAparOldVal(189:190) = DAparVal(189:190)
                DAparOldVal(197:236) = DAparVal(197:236)
                mark_upgrade = 1
                J_last = J_new
            endif
        elseif(iDAsimu < istep_10)then
            ! if(AMIN1(1.0, exp(-(sum(delta_J(2:4))+delta_J(7)+sum(delta_J(9:12))+sum(delta_J(15:20))))) .gt. cs_rand)then
            if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
                DAparOldVal(26:41)   = DAparVal(26:41)
                
                DAparOldVal(59:60)   = DAparVal(59:60)
                DAparOldVal(67:106)  = DAparVal(67:106)

                DAparOldVal(124:125) = DAparVal(124:125)
                DAparOldVal(132:171) = DAparVal(132:171)

                DAparOldVal(189:190) = DAparVal(189:190)
                DAparOldVal(197:236) = DAparVal(197:236)
                mark_upgrade = 1
                J_last = J_new
            endif

        ! elseif(iDAsimu < istep_6)then
        !     if(AMIN1(1.0, exp(-(sum(delta_J(2:4))+delta_J(7)+sum(delta_J(9:12))+sum(delta_J(15:20))))) .gt. cs_rand)then
        !         ! DAparOldVal(26:41)   = DAparVal(26:41)
        !         ! DAparOldVal(67:106)  = DAparVal(67:106)
        !         DAparOldVal(132:171) = DAparVal(132:171)
        !         DAparOldVal(197:236) = DAparVal(197:236)
        !         mark_upgrade = 1
        !         J_last = J_new
        !     endif
        else
            if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
                DAparOldVal = DAparVal
                mark_upgrade = 1
                J_last = J_new
            endif
        endif

        ! if(iDAsimu < 5001)then
        !     if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
        !     ! if(AMIN1(1.0, exp(-delta_J(20))) .gt. cs_rand)then
        !         upgraded = upgraded + 1
        !         J_last   = J_new
        !     endif

        ! elseif(iDAsimu<15001)then
        !     if(AMIN1(1.0, exp(-sum(delta_J(9:13)))) .gt. cs_rand)then
        !         do i = 1,  size(idxSpecialPars_sp(1,:))
        !             DAparOldVal(idxSpecialPars_sp(1,i)) = DAparVal(idxSpecialPars_sp(1,i))
        !         enddo
        !         J_last = J_new
        !         mark_upgrade = 1
        !     endif
        ! elseif(iDAsimu<20001)then
        !     if(AMIN1(1.0, exp(-(sum(delta_J(2:6))+delta_J(20)))) .gt. cs_rand)then
        !         do i = 1,  size(idxSpecialPars_sp(2,:))
        !             DAparOldVal(idxSpecialPars_sp(2,i)) = DAparVal(idxSpecialPars_sp(2,i))
        !         enddo
        !         J_last = J_new
        !         mark_upgrade = 1
        !     endif
        ! elseif(iDAsimu<25001)then
        !     if(AMIN1(1.0, exp(-(sum(delta_J(7:8))))) .gt. cs_rand)then
        !         do i = 1,  size(idxSpecialPars_sp(3,:))
        !             DAparOldVal(idxSpecialPars_sp(3,i)) = DAparVal(idxSpecialPars_sp(3,i))
        !         enddo
        !         J_last = J_new
        !         mark_upgrade = 1
        !     endif
        ! elseif(iDAsimu<30001)then
        !     if(AMIN1(1.0, exp(-(sum(delta_J(14:15))+delta_J(19)))) .gt. cs_rand)then
        !         DAparOldVal(2:23) = DAparVal(2:23)
        !         J_last = J_new
        !         mark_upgrade = 1
        !     endif
        ! endif
        
        if(mark_upgrade == 1) upgraded = upgraded + 1

    end subroutine

    ! subroutine calculate_cost(datObs, datMd, stdObs, JCost)
    !     implicit none
    !     real(8), intent(in) :: datObs(:), datMd(:), stdObs(:)
    !     real(8) :: JCost, dObsSimu, std4cal
    !     integer :: nLine, iLine, nCost
    !     nLine = size(datObs)
    !     nCost = 0
    !     JCost = 0.
    !     do iLine = 1, nLine
    !         if(datObs(iLine)>-999. .and. datMd(iLine)>-999.)then
    !             nCost = nCost+1
    !             dObsSimu = datMd(iLine) - datObs(iLine)
    !             std4cal  = max(0.2 * abs(datObs(iLine)), 0.5)
    !             if (std4cal <= 0) std4cal = 1.0
    !             JCost = JCost + (dObsSimu*dObsSimu)/(2.*std4cal*std4cal)
    !         endif
    !     enddo
    !     if(nCost .gt. 0) JCost=JCost/real(nCost)
    ! end subroutine

    subroutine calculate_cost(datObs, datMd, JCost)
        implicit none
        real(8), intent(in) :: datObs(:), datMd(:)!, stdObs(:)
        real(8) :: JCost, obs_mean, numerator, denominator
        integer :: n, i
        n = size(datObs)
        obs_mean = sum(datObs)/n
        numerator = 0.0
        denominator = 0.0
        do i = 1, n
            if(datObs(i)>-999. .and. datMd(i)>-999.)then
                numerator = numerator + (datMd(i) - datObs(i))**2
                denominator = denominator + (datObs(i) - obs_mean)**2
            endif
        end do
        if (denominator > 0.0) then
            JCost = (numerator / denominator)
        elseif (size(datObs) .eq. 1) then
            Jcost = numerator/(datObs(1) - 83094)**2
        else
            JCost = 0
        endif
        ! if(nCost .gt. 0) JCost=JCost/real(nCost)
        return
    end subroutine calculate_cost


    subroutine mcmc_param_outputs(nUpgraded, st)
        implicit none
        integer, intent(in) :: nUpgraded!, npar4DA, idxStPar(:)
        ! type(index_species_parameters_to_opt), intent(in) :: idxSpPar(:)
        type(site_data_type), intent(inout) :: st

        integer :: npar4st, npar4sp
        integer :: nBuilt_in, ipar, inum, ipft, isimu!, iline
        ! character(250) :: outfile_mc_ParamSets
        character(4000) :: header_line
        integer, allocatable :: rand_number(:)
        real(8), allocatable :: nonBuiltParamsets(:,:), randSelParamsets(:,:)

        npar4st = size(idxStPar)
        ! keep unless 1 parameter set
        nBuilt_in = max(1, int(0.1 * nUpgraded))
        ! get csv header
        header_line = ""
        do ipar = 1, npar4st
            if(ipar .eq. 1) then 
                header_line=trim(mcParams%st(idxStPar(ipar))%parname)
            else
                header_line = trim(header_line) // "," // trim(mcParams%st(idxStPar(ipar))%parname)
            endif
        end do

        do ipft = 1, npft
            npar4sp = size(idxSpPar(ipft)%idx)
            do ipar = 1, npar4sp
                header_line = trim(header_line) // "," // &
                    trim(mcParams%sp(ipft)%var(idxSpPar(ipft)%idx(ipar))%parname) // "_" // sp_names(ipft)
            end do
        end do
        ! copy parameters set (drop built-in part)
        allocate(nonBuiltParamsets(nUpgraded - nBuilt_in, npar4DA))
        nonBuiltParamsets = totMcParamsets(nBuilt_in:nUpgraded, :)
        ! write parameter sets to csv
        call write_parameter_sets(trim(outDir_mcmc) // "/total_parameter_sets.csv", header_line, nonBuiltParamsets)
        ! random select 100 parameter set
        allocate(rand_number(mcset%nRand))
        call generate_random_numbers(1, nUpgraded - nBuilt_in, rand_number)
        allocate(randSelParamsets(mcset%nRand, npar4DA))
        do inum = 1, mcset%nRand
            randSelParamsets(inum, :) = nonBuiltParamsets(rand_number(inum), :)
        end do
        ! write selected parameter sets
        call write_parameter_sets(trim(outDir_mcmc) // "/sel_parameter_sets.csv", header_line, randSelParamsets)

        ! run teco model according to the select variables
        ! do_out_hr   = .True.
        ! do_out_day  = .True.
        ! do_out_mon  = .False.
        ! do_out_yr   = .False.
        do isimu = 1, 1!mcset%nRand
            write(mc_str_n, "(I0.3)") isave
            ! call update_mcParams(randSelParamsets(isimu, :), idxStPar, idxSpPar)
            call update_mcParams(totMcParamsets(nUpgraded, :), idxStPar, idxSpPar)
            call update_simuParams() ! update parameters to simulate

            call initialize_teco(st)
            call teco_simu(st, .True.)
        end do

        deallocate(rand_number)
        deallocate(nonBuiltParamsets)
        deallocate(randSelParamsets)
    end subroutine mcmc_param_outputs

    subroutine generate_random_numbers(min_value, max_value, res_rand)
        implicit none
        integer, dimension(:), intent(inout) :: res_rand
        integer, intent(in) :: min_value, max_value
        integer :: i, j, temp!, range_size, available_numbers
        integer, dimension(max_value - min_value + 1) :: all_numbers
        real(8) :: r

        ! initialize the random
        ! call random_seed()

        ! initilize all_numbers array
        do i = 1, size(all_numbers)
            all_numbers(i) = min_value - 1 + i
        end do

        ! using Fisher-Yates method
        do i = size(all_numbers), 2, -1
            call random_number(r)
            j = int(r * i) + 1
            temp = all_numbers(i)
            all_numbers(i) = all_numbers(j)
            all_numbers(j) = temp
        end do

        ! get the before N random number 
        res_rand = all_numbers(1:size(res_rand,1))
    end subroutine generate_random_numbers

    subroutine write_parameter_sets(filename, header, paramsets)
        implicit none
        character(*), intent(in) :: filename, header
        real(8), intent(in) :: paramsets(:, :)
        integer :: iline
    
        open(118, file=filename, status='replace')
        write(118, *) trim(header)
        do iline = 1, size(paramsets, 1)
            write(118, '(*(ES10.3,:,","))') paramsets(iline,:)
        end do
        close(118)
    end subroutine write_parameter_sets

    subroutine mcmc_vars_deallocate()
        implicit none
        integer i

        do i = 1, numPFT
            if (allocated(idxSpPar(i)%idx)) deallocate(idxSpPar(i)%idx)
        enddo
        if(allocated(totMcParamsets)) deallocate(totMcParamsets)
        if(allocated(DAparVal)) deallocate(DAparVal)
        if(allocated(DAparMin)) deallocate(DAparMin)
        if(allocated(DAparMax)) deallocate(DAparMax)
        if(allocated(DAparOldVal)) deallocate(DAparOldVal)
        if(allocated(idxStPar))    deallocate(idxStPar)
        if(allocated(idxFixStep))  deallocate(idxFixStep)
        if(allocated(search_scale)) deallocate(search_scale)
        if(allocated(mark_obs)) deallocate(mark_obs)
        if(allocated(idxSpecialPars_st)) deallocate(idxSpecialPars_st)
        if(allocated(idxSpecialPars_sp)) deallocate(idxSpecialPars_sp)
    end subroutine mcmc_vars_deallocate
    
end module mcmc