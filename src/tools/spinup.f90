module mod_spinup
    use datatypes
    use driver

    !-------------------------------------------------------------------------------------------------------
    integer iloop
    integer :: nloop = 5000

    type spec_spinup_vars
        real(8), allocatable :: sp_gpp_y(:)
        real(8), allocatable :: sp_npp_y(:)
        real(8), allocatable :: sp_ra_y(:)
        real(8), allocatable :: sp_cLeaf_y(:)
        real(8), allocatable :: sp_cStem_y(:)
        real(8), allocatable :: sp_cRoot_y(:)
        real(8), allocatable :: sp_nsc_y(:)
        ! real(8), allocatable :: sp_cOther_y(:)
        real(8), allocatable :: sp_nLeaf_y(:)
        real(8), allocatable :: sp_nStem_y(:)
        real(8), allocatable :: sp_nRoot_y(:)
        real(8), allocatable :: sp_nsn_y(:)
        ! real(8), allocatable :: sp_nOther_y(:)
        real(8), allocatable :: sp_tran_y(:)
        real(8), allocatable :: sp_lai_y(:)
    end type spec_spinup_vars

    type tot_spinup_vars
        type(spec_spinup_vars), allocatable :: sp(:)
        ! carbon fluxes (Kg C m-2 s-1)
        real(8), allocatable :: sp_gpp_y(:)
        real(8), allocatable :: sp_npp_y(:)
        real(8), allocatable :: sp_ra_y(:)
        real(8), allocatable :: sp_rh_y(:) 
        real(8), allocatable :: sp_wetlandCH4_y(:)
        real(8), allocatable :: sp_wetlandCH4prod_y(:)
        real(8), allocatable :: sp_wetlandCH4cons_y(:)
        ! Carbon Pools  (KgC m-2)
        real(8), allocatable :: sp_cLeaf_y(:)
        real(8), allocatable :: sp_cStem_y(:)
        real(8), allocatable :: sp_cRoot_y(:)
        real(8), allocatable :: sp_cOther_y(:)                             
        real(8), allocatable :: sp_cLitter_y(:)
        real(8), allocatable :: sp_cLitterCwd_y(:)                                         ! litter (excluding coarse Stemy debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse Stemy debris
        real(8), allocatable :: sp_cSoil_y(:)
        real(8), allocatable :: sp_cSoilFast_y(:)
        real(8), allocatable :: sp_cSoilSlow_y(:)
        real(8), allocatable :: sp_cSoilPassive_y(:)                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        real(8), allocatable :: sp_CH4_y(:, :)                                                          ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        real(8), allocatable :: sp_fBNF_y(:)
        real(8), allocatable :: sp_fN2O_y(:)
        real(8), allocatable :: sp_fNloss_y(:)
        real(8), allocatable :: sp_fNnetmin_y(:)
        real(8), allocatable :: sp_fNdep_y(:)                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        real(8), allocatable :: sp_nLeaf_y(:)
        real(8), allocatable :: sp_nStem_y(:)
        real(8), allocatable :: sp_nRoot_y(:)
        real(8), allocatable :: sp_nOther_y(:)
        real(8), allocatable :: sp_nLitter_y(:)
        real(8), allocatable :: sp_nLitterCwd_y(:)
        real(8), allocatable :: sp_nSoil_y(:)
        real(8), allocatable :: sp_nMineral_y(:)                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        real(8), allocatable :: sp_hfls_y(:)
        real(8), allocatable :: sp_hfss_y(:)                                                       ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        real(8), allocatable :: sp_ec_y(:)
        real(8), allocatable :: sp_tran_y(:)
        real(8), allocatable :: sp_es_y(:)                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        real(8), allocatable :: sp_hfsbl_y(:)                                                         ! Snow sublimation
        real(8), allocatable :: sp_mrro_y(:)
        real(8), allocatable :: sp_mrros_y(:)
        real(8), allocatable :: sp_mrrob_y(:)
        real(8), allocatable :: sp_lai_y(:)
        ! test
        ! real(8), allocatable :: sp_test_y(:, :)
    end type tot_spinup_vars
    
    type(tot_spinup_vars) :: sp_outvars

    contains
    subroutine run_spinup(st)
        implicit none
        type(site_data_type), intent(inout) :: st
        integer :: test_i
        write(*,*)"This is spinup: ", nloop
        do iloop = 1, nloop
            write(*,*) "iloop: ", iloop
            test_i = iloop
            ! if (iloop .eq. int(nloop/100))then
            ! call initialize_spinup(st)
            if(mod(iloop, 100) .eq. 0 ) then
                do_out_hr   = .False.
            !     do_out_day  = .True.
            !     do_out_mon  = .False.
            !     do_out_yr   = .False.
                write(mc_str_n, "(I0.3)") test_i
                call teco_simu(st, .True.)
            else
                call teco_simu(st, .False.)
            endif

            ! if (iloop > 145) then
            !     do_out_hr   = .False.
            !     do_out_day  = .True.
            !     do_out_mon  = .False.
            !     do_out_yr   = .False.
            !     write(mc_str_n, "(I0.3)") iloop
            !     call teco_simu(st, .True.)
            ! else
            !     call teco_simu(st, .False.)
            ! endif
            call update_spinup_values()
        enddo
        
    end subroutine run_spinup
    ! ---------------------------------------
    subroutine initialize_spinup(st)
        implicit none
        type(site_data_type), intent(inout) :: st
        ! st%sftmp = -0.
        ! st%Tsnow = -20.
        ! st%Twater = 0.
        ! st%Tice   = 0.
        ! print*, "st%sftmp: ", st%sftmp
        ! print*, "st%G: ", st%G
        ! print*, "st%Tsoill: ", st%Tsoill
        ! st%G      = 20.5

        ! st%Esoil  =0.5*st%G
        ! st%snow_dsim =0.575
        ! st%dcount=50.
        ! st%dcount_soil=50.
        ! st%ice_tw = 0.0
        ! st%Tsoill = (/ -0.09, 0.73, 1.3, 1.95, 2.3, 3., 4., 4.5, 5., 5.98/)
        ! ! st%ice    = ice=(/0.021, 0.0, 0., 0., 0.0, 0.0, 0.0, 0.0,    &
        ! ! &   0.0, 0.0/)
        ! st%Vp(1:3) = 0.
        ! st%Vp(4:6)=0.001
        ! st%Vp(7:10)=0.01
        ! st%bubble_methane_tot  = 0.
        ! st%CH4= (/0.000152,0.05,0.6,0.7,0.7,1.7,1.7,1.7,1.7,1.7/)
        ! st%Nbub = 100.
    end subroutine initialize_spinup
    !-----------------------------------------------------------------------
    subroutine update_spinup_values()
        implicit none 
        integer ipft
        do ipft = 1, npft
            sp_outvars%sp(ipft)%sp_gpp_y(iloop)    = outVars_y%sp(ipft)%gpp
            sp_outvars%sp(ipft)%sp_npp_y(iloop)    = outVars_y%sp(ipft)%npp
            sp_outvars%sp(ipft)%sp_ra_y(iloop)     = outVars_y%sp(ipft)%ra
            sp_outvars%sp(ipft)%sp_cLeaf_y(iloop)  = outVars_y%sp(ipft)%cLeaf
            sp_outvars%sp(ipft)%sp_cStem_y(iloop)  = outVars_y%sp(ipft)%cStem
            sp_outvars%sp(ipft)%sp_cRoot_y(iloop)  = outVars_y%sp(ipft)%cRoot
            sp_outvars%sp(ipft)%sp_nsc_y(iloop)    = outVars_y%sp(ipft)%nsc
            sp_outvars%sp(ipft)%sp_nLeaf_y(iloop)  = outVars_y%sp(ipft)%nLeaf
            sp_outvars%sp(ipft)%sp_nStem_y(iloop)  = outVars_y%sp(ipft)%nStem
            sp_outvars%sp(ipft)%sp_nRoot_y(iloop)  = outVars_y%sp(ipft)%nRoot
            sp_outvars%sp(ipft)%sp_nsn_y(iloop)    = outVars_y%sp(ipft)%nsn
            sp_outvars%sp(ipft)%sp_tran_y(iloop)   = outVars_y%sp(ipft)%tran
            sp_outvars%sp(ipft)%sp_lai_y(iloop)    = outVars_y%sp(ipft)%lai
        enddo
        sp_outvars%sp_gpp_y(iloop)            = outVars_y%gpp
        sp_outvars%sp_npp_y(iloop)            = outVars_y%npp
        sp_outvars%sp_ra_y(iloop)             = outVars_y%ra
        sp_outvars%sp_rh_y(iloop)             = outVars_y%rh
        sp_outvars%sp_wetlandCH4_y(iloop)     = outVars_y%wetlandCH4
        sp_outvars%sp_wetlandCH4prod_y(iloop) = outVars_y%wetlandCH4prod
        sp_outvars%sp_wetlandCH4cons_y(iloop) = outVars_y%wetlandCH4cons
        ! Carbon Pools  (KgC m-2)
        sp_outvars%sp_cLeaf_y(iloop)          = outVars_y%cLeaf
        sp_outvars%sp_cStem_y(iloop)          = outVars_y%cStem
        sp_outvars%sp_cRoot_y(iloop)          = outVars_y%cRoot
        sp_outvars%sp_cOther_y(iloop)         = outVars_y%cOther             
        sp_outvars%sp_cLitter_y(iloop)        = outVars_y%cLitter
        sp_outvars%sp_cLitterCwd_y(iloop)     = outVars_y%cLitterCwd                                         ! litter (excluding coarse Stemy debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse Stemy debris
        sp_outvars%sp_cSoil_y(iloop)          = outVars_y%cSoil
        sp_outvars%sp_cSoilFast_y(iloop)      = outVars_y%cSoilFast
        sp_outvars%sp_cSoilSlow_y(iloop)      = outVars_y%cSoilSlow
        sp_outvars%sp_cSoilPassive_y(iloop)   = outVars_y%cSoilPassive                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        sp_outvars%sp_CH4_y(iloop,:)         = outVars_y%CH4                                                          ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        sp_outvars%sp_fBNF_y(iloop)           = outVars_y%fBNF
        sp_outvars%sp_fN2O_y(iloop)           = outVars_y%fN2O
        sp_outvars%sp_fNloss_y(iloop)         = outVars_y%fNloss
        sp_outvars%sp_fNnetmin_y(iloop)       = outVars_y%fNnetmin
        sp_outvars%sp_fNdep_y(iloop)          = outVars_y%fNdep                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        sp_outvars%sp_nLeaf_y(iloop)          = outVars_y%nLeaf
        sp_outvars%sp_nStem_y(iloop)          = outVars_y%nStem
        sp_outvars%sp_nRoot_y(iloop)          = outVars_y%nRoot
        sp_outvars%sp_nOther_y(iloop)         = outVars_y%nOther
        sp_outvars%sp_nLitter_y(iloop)        = outVars_y%nLitter
        sp_outvars%sp_nLitterCwd_y(iloop)     = outVars_y%nLitterCwd
        sp_outvars%sp_nSoil_y(iloop)          = outVars_y%nSoil
        sp_outvars%sp_nMineral_y(iloop)       = outVars_y%nMineral                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        sp_outvars%sp_hfls_y(iloop)           = outVars_y%hfls
        sp_outvars%sp_hfss_y(iloop)           = outVars_y%hfss                                                       ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        sp_outvars%sp_ec_y(iloop)             = outVars_y%ec
        sp_outvars%sp_tran_y(iloop)           = outVars_y%tran
        sp_outvars%sp_es_y(iloop)             = outVars_y%es                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        sp_outvars%sp_hfsbl_y(iloop)          = outVars_y%hfsbl                                                         ! Snow sublimation
        sp_outvars%sp_mrro_y(iloop)           = outVars_y%mrro
        sp_outvars%sp_mrros_y(iloop)          = outVars_y%mrros
        sp_outvars%sp_mrrob_y(iloop)          = outVars_y%mrrob
        sp_outvars%sp_lai_y(iloop)            = outVars_y%lai
        
    end subroutine update_spinup_values

    !-----------------------------------------------------------------------
    subroutine init_spinup_variables()
        implicit none
        integer ipft
        allocate(sp_outvars%sp(npft))
        do ipft = 1, npft
            allocate(sp_outvars%sp(ipft)%sp_gpp_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_npp_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_ra_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_cLeaf_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_cStem_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_cRoot_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_nsc_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_nLeaf_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_nStem_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_nRoot_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_nsn_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_tran_y(nloop))
            allocate(sp_outvars%sp(ipft)%sp_lai_y(nloop))
        enddo
        allocate(sp_outvars%sp_gpp_y(nloop))
        allocate(sp_outvars%sp_npp_y(nloop))
        allocate(sp_outvars%sp_ra_y(nloop))
        allocate(sp_outvars%sp_rh_y(nloop)) 
        allocate(sp_outvars%sp_wetlandCH4_y(nloop))
        allocate(sp_outvars%sp_wetlandCH4prod_y(nloop))
        allocate(sp_outvars%sp_wetlandCH4cons_y(nloop))
        ! Carbon Pools  (KgC m-2)
        allocate(sp_outvars%sp_cLeaf_y(nloop))
        allocate(sp_outvars%sp_cStem_y(nloop))
        allocate(sp_outvars%sp_cRoot_y(nloop))
        allocate(sp_outvars%sp_cOther_y(nloop))                             
        allocate(sp_outvars%sp_cLitter_y(nloop))
        allocate(sp_outvars%sp_cLitterCwd_y(nloop))                                         ! litter (excluding coarse Stemy debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse Stemy debris
        allocate(sp_outvars%sp_cSoil_y(nloop))
        allocate(sp_outvars%sp_cSoilFast_y(nloop))
        allocate(sp_outvars%sp_cSoilSlow_y(nloop))
        allocate(sp_outvars%sp_cSoilPassive_y(nloop))                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        allocate(sp_outvars%sp_CH4_y(nloop,nlayers))                                                          ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        allocate(sp_outvars%sp_fBNF_y(nloop))
        allocate(sp_outvars%sp_fN2O_y(nloop))
        allocate(sp_outvars%sp_fNloss_y(nloop))
        allocate(sp_outvars%sp_fNnetmin_y(nloop))
        allocate(sp_outvars%sp_fNdep_y(nloop))                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        allocate(sp_outvars%sp_nLeaf_y(nloop))
        allocate(sp_outvars%sp_nStem_y(nloop))
        allocate(sp_outvars%sp_nRoot_y(nloop))
        allocate(sp_outvars%sp_nOther_y(nloop))
        allocate(sp_outvars%sp_nLitter_y(nloop))
        allocate(sp_outvars%sp_nLitterCwd_y(nloop))
        allocate(sp_outvars%sp_nSoil_y(nloop))
        allocate(sp_outvars%sp_nMineral_y(nloop))                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        allocate(sp_outvars%sp_hfls_y(nloop))
        allocate(sp_outvars%sp_hfss_y(nloop))                                                       ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        allocate(sp_outvars%sp_ec_y(nloop))
        allocate(sp_outvars%sp_tran_y(nloop))
        allocate(sp_outvars%sp_es_y(nloop))                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        allocate(sp_outvars%sp_hfsbl_y(nloop))                                                         ! Snow sublimation
        allocate(sp_outvars%sp_mrro_y(nloop))
        allocate(sp_outvars%sp_mrros_y(nloop))
        allocate(sp_outvars%sp_mrrob_y(nloop))
        allocate(sp_outvars%sp_lai_y(nloop))
        ! allocate(sp_test_y(nloop,9))
    end subroutine init_spinup_variables

    subroutine write_spinup_res()         ! write the results of SPIN-UP
        implicit none
        character(2000) :: header_csv
        character(300) :: csv_fileName
        integer :: unit, ipft,iostat,nformat,ilayer
        character(len=2500) :: format_string

        unit = 232
        csv_fileName = adjustl(trim(outDir_spinup))//"/TECO-SPRUCE_spinup.csv"
        open(newunit=unit, file=adjustl(trim(csv_fileName)), status='replace', action='write', iostat=iostat)
        ! Write header line
        header_csv = "isimu,"
        do ipft = 1, npft
            header_csv = adjustl(trim(header_csv))//"gpp_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"npp_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"ra_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"cLeaf_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"cStem_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"cRoot_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nsc_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nLeaf_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nStem_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nRoot_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nsn_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"tran_"//adjustl(trim(sp_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"lai_"//adjustl(trim(sp_names(ipft)))//"," 
        enddo
        header_csv = adjustl(trim(header_csv))//"gpp,npp,ra,rh,wetlandCH4,wetlandCH4prod,&
            wetlandCH4cons,cLeaf,cStem,cRoot,cOther,cLitter,cLitterCwd,cSoil,&
            cSoilFast,cSoilSlow,cSoilPassive,CH4_1,CH4_2,CH4_3,CH4_4,CH4_5,CH4_6,CH4_7,&
            CH4_8,CH4_9,CH4_10,fBNF,fN2O,fNloss,fNnetmin,fNdep,nLeaf,nStem,nRoot,nOther,&
            nLitter,nLitterCwd,nSoil,nMineral,hfls,hfss,ec,tran,es,hfsbl,&
            mrro,mrros,mrrob,lai"
        
        write(unit, *) adjustl(trim(header_csv))
        nformat       = 13*npft+50
        format_string = '((i4,",")' // repeat('(f15.4, ",")', nformat-1) // 'f15.4)'
        do iloop = 1, nloop
            write(unit, adjustl(trim(format_string)))iloop, &
                (sp_outvars%sp(ipft)%sp_gpp_y(iloop),     sp_outvars%sp(ipft)%sp_npp_y(iloop),     &    
                sp_outvars%sp(ipft)%sp_ra_y(iloop),       sp_outvars%sp(ipft)%sp_cLeaf_y(iloop),   & 
                sp_outvars%sp(ipft)%sp_cStem_y(iloop),    sp_outvars%sp(ipft)%sp_cRoot_y(iloop),   & 
                sp_outvars%sp(ipft)%sp_nsc_y(iloop),      sp_outvars%sp(ipft)%sp_nLeaf_y(iloop),   &  
                sp_outvars%sp(ipft)%sp_nStem_y(iloop),   sp_outvars%sp(ipft)%sp_nRoot_y(iloop),&
                sp_outvars%sp(ipft)%sp_nsn_y(iloop),   sp_outvars%sp(ipft)%sp_tran_y(iloop),& 
                sp_outvars%sp(ipft)%sp_lai_y(iloop),      ipft = 1, npft),& 
                sp_outvars%sp_gpp_y(iloop),        sp_outvars%sp_npp_y(iloop),            sp_outvars%sp_ra_y(iloop),      &
                sp_outvars%sp_rh_y(iloop),      sp_outvars%sp_wetlandCH4_y(iloop), sp_outvars%sp_wetlandCH4prod_y(iloop), &
                sp_outvars%sp_wetlandCH4cons_y(iloop),   sp_outvars%sp_cLeaf_y(iloop),      sp_outvars%sp_cStem_y(iloop),    &
                sp_outvars%sp_cRoot_y(iloop),          sp_outvars%sp_cOther_y(iloop),    sp_outvars%sp_cLitter_y(iloop),    &
                sp_outvars%sp_cLitterCwd_y(iloop),     sp_outvars%sp_cSoil_y(iloop),           &
                sp_outvars%sp_cSoilFast_y(iloop),  sp_outvars%sp_cSoilSlow_y(iloop),      sp_outvars%sp_cSoilPassive_y(iloop),  &
                (sp_outvars%sp_CH4_y(iloop,ilayer), ilayer = 1, nlayers),                   &
                sp_outvars%sp_fBNF_y(iloop),       sp_outvars%sp_fN2O_y(iloop),   sp_outvars%sp_fNloss_y(iloop),   &
                sp_outvars%sp_fNnetmin_y(iloop),   sp_outvars%sp_fNdep_y(iloop),   &
                sp_outvars%sp_nLeaf_y(iloop),      sp_outvars%sp_nStem_y(iloop),  sp_outvars%sp_nRoot_y(iloop),    &
                sp_outvars%sp_nOther_y(iloop),     sp_outvars%sp_nLitter_y(iloop), &
                sp_outvars%sp_nLitterCwd_y(iloop), sp_outvars%sp_nSoil_y(iloop),  sp_outvars%sp_nMineral_y(iloop), &
                sp_outvars%sp_hfls_y(iloop),       sp_outvars%sp_hfss_y(iloop),    &
                sp_outvars%sp_ec_y(iloop),         sp_outvars%sp_tran_y(iloop),   sp_outvars%sp_es_y(iloop),      &
                sp_outvars%sp_hfsbl_y(iloop),      sp_outvars%sp_mrro_y(iloop),   sp_outvars%sp_mrros_y(iloop),   &
                 sp_outvars%sp_mrrob_y(iloop),                       &     
                sp_outvars%sp_lai_y(iloop)
        enddo
        close(unit)
        return
    end subroutine write_spinup_res

    subroutine write_restart()            ! write the result file
        implicit none
        print*, "write_restart ..."
    end subroutine write_restart
    
    !-----------------------------------------------------------------------
    subroutine deallo_spinup_variables()
        implicit none
        integer ipft
        if(allocated(sp_outvars%sp))then
            do ipft = 1, npft
                if (allocated(sp_outvars%sp(ipft)%sp_gpp_y)) deallocate(sp_outvars%sp(ipft)%sp_gpp_y)
                if (allocated(sp_outvars%sp(ipft)%sp_npp_y))deallocate(sp_outvars%sp(ipft)%sp_npp_y)
                if (allocated(sp_outvars%sp(ipft)%sp_ra_y))deallocate(sp_outvars%sp(ipft)%sp_ra_y)
                if (allocated(sp_outvars%sp(ipft)%sp_cLeaf_y))deallocate(sp_outvars%sp(ipft)%sp_cLeaf_y)
                if (allocated(sp_outvars%sp(ipft)%sp_cStem_y))deallocate(sp_outvars%sp(ipft)%sp_cStem_y)
                if (allocated(sp_outvars%sp(ipft)%sp_cRoot_y))deallocate(sp_outvars%sp(ipft)%sp_cRoot_y)
                ! if (allocated(sp_outvars%sp(ipft)%sp_cOther_y))deallocate(sp_outvars%sp(ipft)%sp_cOther_y)
                if (allocated(sp_outvars%sp(ipft)%sp_nLeaf_y))deallocate(sp_outvars%sp(ipft)%sp_nLeaf_y)
                if (allocated(sp_outvars%sp(ipft)%sp_nStem_y))deallocate(sp_outvars%sp(ipft)%sp_nStem_y)
                if (allocated(sp_outvars%sp(ipft)%sp_nRoot_y))deallocate(sp_outvars%sp(ipft)%sp_nRoot_y)
                ! if (allocated(sp_outvars%sp(ipft)%sp_nOther_y))deallocate(sp_outvars%sp(ipft)%sp_nOther_y)
                if (allocated(sp_outvars%sp(ipft)%sp_lai_y))   deallocate(sp_outvars%sp(ipft)%sp_lai_y)
            enddo
            deallocate(sp_outvars%sp)
        endif
        if (allocated(sp_outvars%sp_gpp_y)) deallocate(sp_outvars%sp_gpp_y)
        if (allocated(sp_outvars%sp_npp_y)) deallocate(sp_outvars%sp_npp_y)
        if (allocated(sp_outvars%sp_ra_y)) deallocate(sp_outvars%sp_ra_y)
        if (allocated(sp_outvars%sp_rh_y)) deallocate(sp_outvars%sp_rh_y) 
        if (allocated(sp_outvars%sp_wetlandCH4_y)) deallocate(sp_outvars%sp_wetlandCH4_y)
        if (allocated(sp_outvars%sp_wetlandCH4prod_y)) deallocate(sp_outvars%sp_wetlandCH4prod_y)
        if (allocated(sp_outvars%sp_wetlandCH4cons_y)) deallocate(sp_outvars%sp_wetlandCH4cons_y)
        ! Carbon Pools  (KgC m-2)
        if (allocated(sp_outvars%sp_cLeaf_y)) deallocate(sp_outvars%sp_cLeaf_y)
        if (allocated(sp_outvars%sp_cStem_y)) deallocate(sp_outvars%sp_cStem_y)
        if (allocated(sp_outvars%sp_cRoot_y)) deallocate(sp_outvars%sp_cRoot_y)
        if (allocated(sp_outvars%sp_cOther_y)) deallocate(sp_outvars%sp_cOther_y)                             
        if (allocated(sp_outvars%sp_cLitter_y)) deallocate(sp_outvars%sp_cLitter_y)
        if (allocated(sp_outvars%sp_cLitterCwd_y)) deallocate(sp_outvars%sp_cLitterCwd_y)                                         ! litter (excluding coarse Stemy debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse Stemy debris
        if (allocated(sp_outvars%sp_cSoil_y)) deallocate(sp_outvars%sp_cSoil_y)
        if (allocated(sp_outvars%sp_cSoilFast_y)) deallocate(sp_outvars%sp_cSoilFast_y)
        if (allocated(sp_outvars%sp_cSoilSlow_y)) deallocate(sp_outvars%sp_cSoilSlow_y)
        if (allocated(sp_outvars%sp_cSoilPassive_y)) deallocate(sp_outvars%sp_cSoilPassive_y)                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        if (allocated(sp_outvars%sp_CH4_y)) deallocate(sp_outvars%sp_CH4_y)                                                          ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        if (allocated(sp_outvars%sp_fBNF_y)) deallocate(sp_outvars%sp_fBNF_y)
        if (allocated(sp_outvars%sp_fN2O_y)) deallocate(sp_outvars%sp_fN2O_y)
        if (allocated(sp_outvars%sp_fNloss_y)) deallocate(sp_outvars%sp_fNloss_y)
        if (allocated(sp_outvars%sp_fNnetmin_y)) deallocate(sp_outvars%sp_fNnetmin_y)
        if (allocated(sp_outvars%sp_fNdep_y)) deallocate(sp_outvars%sp_fNdep_y)                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        if (allocated(sp_outvars%sp_nLeaf_y)) deallocate(sp_outvars%sp_nLeaf_y)
        if (allocated(sp_outvars%sp_nStem_y)) deallocate(sp_outvars%sp_nStem_y)
        if (allocated(sp_outvars%sp_nRoot_y)) deallocate(sp_outvars%sp_nRoot_y)
        if (allocated(sp_outvars%sp_nOther_y)) deallocate(sp_outvars%sp_nOther_y)
        if (allocated(sp_outvars%sp_nLitter_y)) deallocate(sp_outvars%sp_nLitter_y)
        if (allocated(sp_outvars%sp_nLitterCwd_y)) deallocate(sp_outvars%sp_nLitterCwd_y)
        if (allocated(sp_outvars%sp_nSoil_y)) deallocate(sp_outvars%sp_nSoil_y)
        if (allocated(sp_outvars%sp_nMineral_y)) deallocate(sp_outvars%sp_nMineral_y)                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        if (allocated(sp_outvars%sp_hfls_y)) deallocate(sp_outvars%sp_hfls_y)
        if (allocated(sp_outvars%sp_hfss_y)) deallocate(sp_outvars%sp_hfss_y)                                                       ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        if (allocated(sp_outvars%sp_ec_y))    deallocate(sp_outvars%sp_ec_y)
        if (allocated(sp_outvars%sp_tran_y))  deallocate(sp_outvars%sp_tran_y)
        if (allocated(sp_outvars%sp_es_y))    deallocate(sp_outvars%sp_es_y)                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        if (allocated(sp_outvars%sp_hfsbl_y)) deallocate(sp_outvars%sp_hfsbl_y)                                                         ! Snow sublimation
        if (allocated(sp_outvars%sp_mrro_y))  deallocate(sp_outvars%sp_mrro_y)
        if (allocated(sp_outvars%sp_mrros_y)) deallocate(sp_outvars%sp_mrros_y)
        if (allocated(sp_outvars%sp_mrrob_y)) deallocate(sp_outvars%sp_mrrob_y)
        if (allocated(sp_outvars%sp_lai_y))   deallocate(sp_outvars%sp_lai_y)
        ! deallocate(sp_test_y)
    end subroutine deallo_spinup_variables

    !-----------------------------------------------------------------------
end module mod_spinup