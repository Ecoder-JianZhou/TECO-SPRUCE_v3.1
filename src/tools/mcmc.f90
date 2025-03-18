module mcmc
    use datatypes
    use driver ! run teco simulation
    use mcmc_mod

    integer npar4DA ! number of parameters to be optimize
    integer upgraded

    logical :: do_cov2create_new_params ! if do cov to create new parameters
    real(8), allocatable :: DAparVal(:)
    real(8), allocatable :: DAparMin(:)
    real(8), allocatable :: DAparMax(:)
    real(8), allocatable :: DAparOldVal(:)
    integer, allocatable :: idxStPar(:)
    type(index_species_parameters_to_opt) :: idxSpPar(numPFT) ! type defined in MCMC_mod

    real(8) J_last(numObsFiles), J_new(numObsFiles)
    ! save parameters in MCMC processes
    real(8), allocatable :: totMcParamsets(:,:) ! upgraded, nParams

    contains

    subroutine run_mcmc(st)
        implicit none
        type(site_data_type), intent(inout) :: st
        integer :: temp_upgraded, iDAsimu, i
        real(8) :: upgraded_rate

        print*, "# Start to run MCMC ... "
        ! 1. generate new parameters
        call generate_new_parameters()
        upgraded = 0
        ! 2. do cylces
        do iDAsimu = 1, mcset%nDAsimu
            upgraded_rate = (real(upgraded, 8) / real(iDAsimu, 8))*100
            write(*, '(I6, A1, I6, 2X, F12.2, 2X, F12.2, 2X, I6, 2X, F6.2)') &
                iDAsimu, '/', mcset%nDAsimu, sum(J_new), sum(J_last), upgraded, upgraded_rate
                ! initilize
            do i = 1, numObsFiles
                mcVarData(i)%mark_idx = 1
            end do

            ! update mcparams
            call generate_new_parameters()

            call update_mcParams(DAparVal, idxStPar, idxSpPar)
            call update_simuParams() ! update parameters to simulate

            call initialize_teco(st)
            call teco_simu(st, .False.)

            temp_upgraded = upgraded
            call cal_cost_function() ! check if update upgraded
            if(upgraded > temp_upgraded)then
                DAparOldVal = DAparVal
                totMcParamsets(upgraded,:) = DAparVal
            endif
        enddo
        call mcmc_param_outputs(upgraded, st)
        ! 3. outputs 
    end subroutine run_mcmc

    subroutine init_mcmc(in_mcConfNmlFile)
        implicit none
        character(*), intent(in) :: in_mcConfNmlFile

        call read_mcmc_conf_nml(in_mcConfNmlFile)
        ! get paramters to optimize
        call filter_params_to_optimize()
        J_last = 9000000.0
        ! initilize outputs
        allocate(totMcParamsets(mcset%nDAsimu, npar4DA))
    end subroutine init_mcmc

    subroutine filter_params_to_optimize()
        implicit none
        real(8), dimension(numAllParams) :: temp_parmin, temp_parmax, temp_parval
        integer :: temp_idx_st(numStParams)
        integer :: i, j, nSpPar4DA(numPFT), nStPar4DA
        type(index_species_parameters_to_opt) :: temp_idx_sp(numPFT)

        npar4DA = 0
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
            endif
        enddo

        nStPar4DA = npar4DA

        do i = 1, numPFT
            nSpPar4DA(i) = 0
            allocate(temp_idx_sp(i)%idx(numSpParams))
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
                endif
            enddo
            allocate(idxSpPar(i)%idx(nSpPar4DA(i)))
            idxSpPar(i)%idx = temp_idx_sp(i)%idx(:nSpPar4DA(i))
            if(allocated(temp_idx_sp(i)%idx)) deallocate(temp_idx_sp(i)%idx)
        enddo

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

        call random_seed()
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
                    DAparVal(i) = DAparOldVal(i) + rand*(DAparMax(i)-DAparMin(i))*mcset%search_scale
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
        integer i
        real(8) :: JCost, delta_J(numObsFiles), cs_rand

        J_new = 0.
        do i = 1, numObsFiles
            if(mcVarData(i)%existOrNot)then
                call calculate_cost(mcVarData(i)%obsData(:,4), mcVarData(i)%mdData(:,4),JCost)
                J_new(i) = J_new(i) + JCost*100*obsWt(i)
            end if
            ! print*, "test_", i,": ", J_new(i) - J_last(i)
            ! print*, "test_",i,": ",J_new(i)
        enddo

        ! print*, mcVarData(7)%mdData(1,4),mcVarData(7)%obsData(1,4)

        delta_J = J_new - J_last
        call random_number(cs_rand)
        if(AMIN1(1.0, exp(-sum(delta_J))) .gt. cs_rand)then
            upgraded = upgraded + 1
            J_last   = J_new
        endif

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
            numerator = numerator + (datMd(i) - datObs(i))**2
            denominator = denominator + (datObs(i) - obs_mean)**2
        end do
        if (denominator > 0.0) then
            JCost = (numerator / denominator)
        else
            ! JCost = 0
            Jcost = numerator**0.5
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
        character(1200) :: header_line
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
        do_out_hr   = .True.
        do_out_day  = .True.
        do_out_mon  = .False.
        do_out_yr   = .False.
        do isimu = 1, 1!mcset%nRand
            write(mc_str_n, "(I0.3)") isimu
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
        call random_seed()

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
    end subroutine mcmc_vars_deallocate
    
end module mcmc