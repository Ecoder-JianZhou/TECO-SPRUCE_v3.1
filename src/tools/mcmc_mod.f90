module mcmc_mod
    ! some functions that both driver module and MCMC module.
    use datatypes
    implicit none

    integer :: npar, nDAsimu, ncov, nRand,upgraded,iDAsimu
    real(8)    :: search_scale
    logical :: do_mc_out_hr, do_mc_out_day, do_mc_out_mon, do_mc_out_yr
    integer :: nsite_params, nspec_params, npar4st, npar4sp

    ! real(8) mc_lat, mc_Longitude, mc_wsmax, mc_wsmin
    ! real(8) mc_LAIMAX, mc_LAIMIN, mc_rdepth, mc_Rootmax, mc_Stemmax
    ! real(8) mc_SapR, mc_SapS, mc_SLA, mc_GLmax, mc_GRmax, mc_Gsmax, mc_stom_n
    ! real(8) mc_a1, mc_Ds0, mc_Vcmx0, mc_extkU, mc_xfang, mc_alpha
    ! real(8) mc_Tau_Leaf, mc_Tau_Stem, mc_Tau_Root, mc_Tau_F
    ! real(8) mc_Tau_C,  mc_Tau_Micro, mc_Tau_SlowSOM, mc_Tau_Passive
    ! real(8) mc_gddonset, mc_Q10, mc_Rl0, mc_Rs0, mc_Rr0
    ! real(8) mc_r_me, mc_Q10pro, mc_kCH4, mc_Omax, mc_CH4_thre
    ! real(8) mc_Tveg, mc_Tpro_me, mc_Toxi
    ! real(8) mc_f, mc_bubprob, mc_Vmaxfraction
    ! real(8) mc_Q10rh, mc_JV, mc_Entrpy
    ! real(8) mc_etaL, mc_etaW, mc_etaR
    ! real(8) mc_f_F2M, mc_f_C2M, mc_f_C2S, mc_f_M2S
    ! real(8) mc_f_M2P, mc_f_S2P, mc_f_S2M, mc_f_P2M

    type mc_tile_params
        character(50) :: parname
        real(8) :: parval
        real(8) :: parmin
        real(8) :: parmax
    end type mc_tile_params

    type mc_spec_params
        type(mc_tile_params), allocatable :: sp_params(:)
    end type mc_spec_params
    
    type mc_params_type
        type(mc_spec_params), allocatable :: sp(:)
        type(mc_tile_params), allocatable :: st(:)
    end type mc_params_type
    type(mc_params_type) :: mc_params

    ! type mc_in_params_type
    !     type(in_site_params) :: max_st_params
    !     type(in_site_params) :: min_st_params
    !     type(in_spec_params), allocatable :: max_sp_params(:)
    !     type(in_spec_params), allocatable :: min_sp_params(:)
    ! end type mc_in_params_type

    ! type(mc_in_params_type) :: mc_in_params

    ! type params_mcmc
    !     real(8) lat, Longitude, wsmax, wsmin
    !     real(8) LAIMAX, LAIMIN, rdepth, Rootmax, Stemmax
    !     real(8) SapR, SapS, SLA, GLmax, GRmax, Gsmax, stom_n
    !     real(8) a1, Ds0, Vcmx0, extkU, xfang, alpha
    !     real(8) Tau_Leaf, Tau_Stem, Tau_Root, Tau_F
    !     real(8) Tau_C,  Tau_Micro, Tau_SlowSOM, Tau_Passive
    !     real(8) gddonset, Q10, Rl0, Rs0, Rr0
    !     real(8) r_me, Q10pro, kCH4, Omax, CH4_thre
    !     real(8) Tveg, Tpro_me, Toxi
    !     real(8) f, bubprob, Vmaxfraction
    !     real(8) Q10rh, JV, Entrpy
    !     real(8) etaL, etaW, etaR
    !     real(8) f_F2M, f_C2M, f_C2S, f_M2S
    !     real(8) f_M2P, f_S2P, f_S2M, f_P2M
    ! end type params_mcmc

    type params_mcmc
        real(8), allocatable :: parval(:)
        real(8), allocatable :: parmin(:)
        real(8), allocatable :: parmax(:)
    end type params_mcmc
    type(params_mcmc), allocatable :: mc_parvals(:)

    type params_DApar
        real(8), allocatable :: DAparmin(:)
        real(8), allocatable :: DAparmax(:)
        real(8), allocatable :: DApar(:)
        real(8), allocatable :: DApar_old(:)
        ! integer, allocatable :: DAparidx(:)
        integer, allocatable :: DAparidx_st(:)
        integer, allocatable :: DAparidx_sp(:,:)
        real(8), allocatable :: gamma(:,:)
        real(8), allocatable :: gamnew(:,:)
        real(8), allocatable :: coefhistory(:,:)
        real(8), allocatable :: coefnorm(:)
        real(8), allocatable :: coefac(:)
    end type params_DApar
    type(params_DApar) :: mc_DApar          ! all variables for DATA ASSIMILATION

    ! type(nml_params_data_type),allocatable     :: mc_in_params(:)     !   in_params for MCMC
    ! type(nml_initValue_data_type), allocatable :: mc_init_params(:)   ! init_params for MCMC

    type params_sets
        real(8), allocatable :: tot_paramsets(:,:)
        real(8), allocatable :: upg_paramsets(:,:)
        real(8), allocatable :: sel_paramsets(:,:) 
    end type params_sets
    type(params_sets) :: arr_params_set

    ! type(nml_params_data_type) :: in_params, in_parval, in_parval_min, in_parval_max
    ! real(8), allocatable :: parval(:), parmin(:), parmax(:)
    character(20), allocatable :: parnames(:)

    ! observational file path
    character(500) :: obsfile_ANPP_Shrub_y
    character(500) :: obsfile_ANPP_Tree_y
    character(500) :: obsfile_NPP_sphag_y
    character(500) :: obsfile_BNPP_y        ! tree + shrub
    character(500) :: obsfile_er_d          ! shrub + sphag.
    character(500) :: obsfile_er_h          ! shrub + sphag.
    character(500) :: obsfile_gpp_d         ! Shrub + sphag.
    character(500) :: obsfile_nee_d         ! Shrub + sphag.
    character(500) :: obsfile_nee_h         ! shrub + sphag.
    character(500) :: obsfile_LAI_d         ! tree  + Shrub
    character(500) :: obsfile_rh_y
    !
    character(500) :: obsfile_leaf_mass_shrub_y
    character(500) :: obsfile_stem_mass_shrub_y
    character(500) :: obsfile_cPlant_tree_y
    character(500) :: obsfile_cPlant_sphag_y
    character(500) :: obsfile_cSoil_y
    character(500) :: obsfile_leaf_resp_shrub_d 
    character(500) :: obsfile_leaf_resp_tree_d 
    ! methane
    character(500) :: obsfile_ch4_d 
    character(500) :: obsfile_ch4_h 
    character(500) :: obsfile_ch4_y

    character(500) :: obsfile_watertable_h
    ! 
    character(500) :: obsfile_CN_shag_d 
    character(500) :: obsfile_photo_shrub_d 
    character(500) :: obsfile_photo_tree_d  
    !
    character(500) :: obsfile_LAI_tree_d 
    character(500) :: obsfile_LAI_shrub_d 
    character(500) :: obsfile_photo_tree_h 
    !
    character(500) :: obsfile_bnpp_tree_y 
    character(500) :: obsfile_bnpp_shrub_y 
    !
    character(500) :: obsfile_gpp_tree_y 

    ! variables for calculating the cost in MCMC processes
    type interCostVariable
        character(300) :: filepath
        logical :: existOrNot
        real(8), allocatable :: obsData(:,:)
        real(8), allocatable :: mdData(:,:)
        integer :: mc_itime
    end type interCostVariable

    type allCostVariables
    ! default variables, you can add the variable names here. (year, doy, hour, value, std.)
        ! carbon flux 
        type(interCostVariable) :: ANPP_Shrub_y
        type(interCostVariable) :: ANPP_Tree_y
        type(interCostVariable) :: NPP_sphag_y
        type(interCostVariable) :: BNPP_y        ! tree + shrub
        type(interCostVariable) :: er_d          ! shrub + sphag.
        type(interCostVariable) :: er_h          ! shrub + sphag.
        type(interCostVariable) :: gpp_d         ! Shrub + sphag.
        type(interCostVariable) :: nee_d         ! Shrub + sphag.
        type(interCostVariable) :: nee_h         ! shrub + sphag.
        type(interCostVariable) :: LAI_d         ! tree  + Shrub
        type(interCostVariable) :: rh_y
        !
        type(interCostVariable) :: leaf_mass_shrub_y
        type(interCostVariable) :: stem_mass_shrub_y
        type(interCostVariable) :: cPlant_tree_y
        type(interCostVariable) :: cPlant_sphag_y
        type(interCostVariable) :: cSoil_y
        type(interCostVariable) :: leaf_resp_shrub_d 
        type(interCostVariable) :: leaf_resp_tree_d 
        ! methane
        type(interCostVariable) :: ch4_d 
        type(interCostVariable) :: ch4_h 
        type(interCostVariable) :: ch4_y
        ! 
        type(interCostVariable) :: CN_shag_d 
        type(interCostVariable) :: photo_shrub_d 
        type(interCostVariable) :: photo_tree_d 
        type(interCostVariable) :: zwt_h
        !
        type(interCostVariable) :: lai_tree_d 
        type(interCostVariable) :: lai_shrub_d 
        type(interCostVariable) :: photo_tree_h
        !
        type(interCostVariable) :: bnpp_shrub_y 
        type(interCostVariable) :: bnpp_tree_y 
        !
        type(interCostVariable) :: gpp_tree_y ! just to keep tree gpp reasonable with CUE = 0.5
    end type allCostVariables

    type(allCostVariables) :: vars4MCMC      ! define a allCostVariables first

    ! variables for marking the cycle number
    integer mc_itime_gpp_d, mc_itime_nee_d, mc_itime_reco_d
    integer mc_itime_gpp_h, mc_itime_nee_h, mc_itime_reco_h
    integer mc_itime_ch4_h, mc_itime_cleaf, mc_itime_cStem
    integer mc_itime_anpp_y, mc_itime_bnpp_y, mc_itime_lai_h
    integer mc_itime_npp_y, mc_itime_reco_y
    integer mc_iyear,  mc_iday, mc_ihour

    type mcmc_spec_outvars_type
        character(50) :: spec_name
        ! carbon fluxes (Kg C m-2 s-1)
        real(8), allocatable :: gpp(:, :)
        real(8), allocatable :: nee(:, :)
        real(8), allocatable :: npp(:, :)
        real(8), allocatable :: nppLeaf(:, :)
        real(8), allocatable :: nppStem(:, :)
        real(8), allocatable :: nppRoot(:, :)
        real(8), allocatable :: nppOther(:, :)    ! According to SPRUCE-MIP, stem means above ground Stemy tissues which is different from Stem tissues.
        real(8), allocatable :: ra(:, :)
        real(8), allocatable :: raLeaf(:, :)
        real(8), allocatable :: raStem(:, :)
        real(8), allocatable :: raRoot(:, :)
        real(8), allocatable :: raOther(:, :)
        real(8), allocatable :: rMaint(:, :)
        real(8), allocatable :: rGrowth(:, :)
        real(8), allocatable :: nbp(:, :)
        ! Carbon Pools  (KgC m-2)
        real(8), allocatable :: cLeaf(:, :)
        real(8), allocatable :: cStem(:, :)
        real(8), allocatable :: cRoot(:, :)
        ! Nitrogen pools (kgN m-2)
        real(8), allocatable :: nLeaf(:, :)
        real(8), allocatable :: nStem(:, :)
        real(8), allocatable :: nRoot(:, :)
        ! real(8), allocatable :: nOther(:)
        ! water fluxes (kg m-2 s-1)
        real(8), allocatable :: tran(:, :)
        ! other
        real(8), allocatable :: lai(:, :)                     ! m2 m-2, Leaf area index
    end type mcmc_spec_outvars_type

    type mcmc_outVars_type
        type(mcmc_spec_outvars_type), allocatable :: sp(:)
        ! carbon fluxes (Kg C m-2 s-1)
        real(8), allocatable :: gpp(:, :)
        real(8), allocatable :: nee(:, :)
        real(8), allocatable :: npp(:, :)
        real(8), allocatable :: nppLeaf(:, :)
        real(8), allocatable :: nppStem(:, :)
        real(8), allocatable :: nppRoot(:, :)
        real(8), allocatable :: nppOther(:, :)           ! According to SPRUCE-MIP, stem means above ground Stemy tissues which is different from Stem tissues.
        real(8), allocatable :: ra(:, :)
        real(8), allocatable :: raLeaf(:, :)
        real(8), allocatable :: raStem(:, :)
        real(8), allocatable :: raRoot(:, :)
        real(8), allocatable :: raOther(:, :)
        real(8), allocatable :: rMaint(:, :)
        real(8), allocatable :: rGrowth(:, :)            ! maintenance respiration and growth respiration
        real(8), allocatable :: rh(:, :)
        real(8), allocatable :: nbp(:, :)                ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        real(8), allocatable :: wetlandCH4(:, :)
        real(8), allocatable :: wetlandCH4prod(:, :)
        real(8), allocatable :: wetlandCH4cons(:, :)     ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        real(8), allocatable :: cLeaf(:, :)
        real(8), allocatable :: cStem(:, :)
        real(8), allocatable :: cRoot(:, :)
        real(8), allocatable :: cOther(:, :)              ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        real(8), allocatable :: cLitter(:, :)
        real(8), allocatable :: cLitterCwd(:, :)          ! litter (excluding coarse Stemy debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse Stemy debris
        real(8), allocatable :: cSoil(:, :)
        real(8), allocatable :: cSoilLevels(:, :, :)
        real(8), allocatable :: cSoilFast(:, :)
        real(8), allocatable :: cSoilSlow(:, :)
        real(8), allocatable :: cSoilPassive(:, :)           ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        real(8), allocatable :: CH4(:, :, :)          ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        real(8), allocatable :: fBNF(:, :)
        real(8), allocatable :: fN2O(:, :)
        real(8), allocatable :: fNloss(:, :)
        real(8), allocatable :: fNnetmin(:, :)
        real(8), allocatable :: fNdep(:, :)                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        real(8), allocatable :: nLeaf(:, :)
        real(8), allocatable :: nStem(:, :)
        real(8), allocatable :: nRoot(:, :)
        real(8), allocatable :: nOther(:, :)
        real(8), allocatable :: nLitter(:, :)
        real(8), allocatable :: nLitterCwd(:, :)
        real(8), allocatable :: nSoil(:, :)
        real(8), allocatable :: nMineral(:, :)                ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        real(8), allocatable :: hfls(:, :)
        real(8), allocatable :: hfss(:, :)
        real(8), allocatable :: SWnet(:, :)
        real(8), allocatable :: LWnet(:, :)                   ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        real(8), allocatable :: ec(:, :)
        real(8), allocatable :: tran(:, :)
        real(8), allocatable :: es(:, :)                      ! Canopy evaporation; Canopy transpiration; Soil evaporation
        real(8), allocatable :: hfsbl(:, :)                   ! Snow sublimation
        real(8), allocatable :: mrro(:, :)
        real(8), allocatable :: mrros(:, :)
        real(8), allocatable :: mrrob(:, :)                   ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        real(8), allocatable :: mrso(:, :, :)           ! Kg m-2, soil moisture in each soil layer
        real(8), allocatable :: tsl(:, :, :)            ! K, soil temperature in each soil layer
        real(8), allocatable :: tsland(:, :)                  ! K, surface temperature
        real(8), allocatable :: wtd(:, :)                     ! m, Water table depth
        real(8), allocatable :: snd(:, :)                     ! m, Total snow depth
        real(8), allocatable :: lai(:, :)                     ! m2 m-2, Leaf area index            
    end type mcmc_outVars_type

    type(mcmc_outVars_type) sel_paramsets_outs_h
    type(mcmc_outVars_type) sel_paramsets_outs_d
    type(mcmc_outVars_type) sel_paramsets_outs_m
    ! total simulation outputs
    type(mcmc_outVars_type) tot_paramsets_outs_h
    type(mcmc_outVars_type) tot_paramsets_outs_d
    type(mcmc_outVars_type) tot_paramsets_outs_m

    contains

    subroutine mcmc_functions_init()
        implicit none
        vars4MCMC%ANPP_Shrub_y%mc_itime = 1
        vars4MCMC%ANPP_Tree_y%mc_itime  = 1
        vars4MCMC%NPP_sphag_y%mc_itime  = 1
        vars4MCMC%BNPP_y%mc_itime       = 1
        vars4MCMC%er_d%mc_itime         = 1
        vars4MCMC%er_h%mc_itime         = 1
        vars4MCMC%gpp_d%mc_itime        = 1
        vars4MCMC%nee_d%mc_itime        = 1
        vars4MCMC%nee_h%mc_itime        = 1
        vars4MCMC%LAI_d%mc_itime        = 1
        vars4MCMC%rh_y%mc_itime         = 1
        
        vars4MCMC%leaf_mass_shrub_y%mc_itime = 1
        vars4MCMC%stem_mass_shrub_y%mc_itime = 1
        vars4MCMC%cPlant_tree_y%mc_itime     = 1
        vars4MCMC%cPlant_sphag_y%mc_itime    = 1
        vars4MCMC%cSoil_y%mc_itime           = 1
        vars4MCMC%leaf_resp_shrub_d%mc_itime = 1
        vars4MCMC%leaf_resp_tree_d%mc_itime  = 1 

        vars4MCMC%ch4_d%mc_itime = 1
        vars4MCMC%ch4_h%mc_itime = 1
        vars4MCMC%ch4_y%mc_itime = 1

        vars4MCMC%CN_shag_d%mc_itime     = 1
        vars4MCMC%photo_shrub_d%mc_itime = 1
        vars4MCMC%photo_tree_d%mc_itime  = 1

        vars4MCMC%zwt_h%mc_itime  = 1

        vars4MCMC%lai_shrub_d%mc_itime  = 1
        vars4MCMC%lai_tree_d%mc_itime   = 1
        vars4MCMC%photo_tree_h%mc_itime = 1

        vars4MCMC%bnpp_tree_y%mc_itime  = 1
        vars4MCMC%bnpp_shrub_y%mc_itime   = 1

        vars4MCMC%gpp_tree_y%mc_itime  = 1

        mc_iyear = 1
        mc_iday  = 1
        mc_ihour = 1
    end subroutine mcmc_functions_init

    subroutine readMCMC_configs_NML(in_mcmc_configfile)
        implicit none
        character(*), intent(in) :: in_mcmc_configfile
        ! type(site_data_type), intent(inout) :: st
        integer io, ipft

        namelist /nml_mcmc_settings/ nDAsimu, search_scale, ncov, nRand, &
                do_mc_out_hr, do_mc_out_day, do_mc_out_mon, do_mc_out_yr
        
        namelist /nml_obsfiles/ obsfile_ANPP_Shrub_y, obsfile_ANPP_Tree_y, obsfile_NPP_sphag_y, &
            obsfile_BNPP_y, obsfile_er_d, obsfile_er_h, obsfile_gpp_d, obsfile_nee_d, &
            obsfile_nee_h, obsfile_LAI_d, obsfile_rh_y, obsfile_leaf_mass_shrub_y, obsfile_stem_mass_shrub_y, &
            obsfile_leaf_resp_shrub_d, obsfile_leaf_resp_tree_d, obsfile_ch4_d, obsfile_ch4_h, obsfile_ch4_y, & 
            obsfile_CN_shag_d, obsfile_photo_shrub_d, obsfile_photo_tree_d, &
            obsfile_cPlant_tree_y, obsfile_cSoil_y, obsfile_cPlant_sphag_y, obsfile_watertable_h, &
            obsfile_LAI_shrub_d, obsfile_LAI_tree_d, obsfile_photo_tree_h, &
            obsfile_bnpp_tree_y, obsfile_bnpp_shrub_y, obsfile_gpp_tree_y

        ! character(50) :: parname_1,  parname_2,  parname_3,  parname_4,  parname_5,  parname_6
        ! character(50) :: parname_7,  parname_8,  parname_9,  parname_10, parname_11, parname_12
        ! character(50) :: parname_13, parname_14, parname_15, parname_16, parname_17, parname_18
        ! character(50) :: parname_19, parname_20, parname_21, parname_22, parname_23, parname_24
        ! character(50) :: parname_25, parname_26, parname_27, parname_28, parname_29, parname_30
        ! character(50) :: parname_31

        real(8) :: wsmax(3), wsmin(3), extkU(3), Tau_F(3), Tau_C(3), Tau_Micro(3), Tau_SlowSOM(3), Tau_Passive(3) 
        real(8) :: Q10rh(3), etaW(3), f_F2M(3), f_C2M(3), f_C2S(3), f_M2S(3), f_M2P(3), f_S2P(3), f_S2M(3), f_P2M(3)
        real(8) :: r_me(3),  Q10pro(3), kCH4(3), Omax(3), CH4_thre(3), Tveg(3), Tpro_me(3), Toxi(3) 
        real(8) :: f(3),     bubprob(3), Vmaxfraction(3)
        real(8) :: f_fast(3), f_slow(3), s_soil(3) 
        real(8) :: par_shcap_snow(3), par_condu_snow(3), par_condu_b(3)
        real(8) :: par_fsub(3), par_rho_snow(3), par_decay_m(3)

        ! namelist /nml_site_params/ parname_1, wsmax, parname_2, wsmin, parname_3, extkU, &
        !     parname_4, Tau_F, parname_5, Tau_C, parname_6, Tau_Micro, parname_7, Tau_SlowSOM, &
        !     parname_8, Tau_Passive, parname_9, Q10rh, parname_10, f_F2M, parname_11, f_C2M, &
        !     parname_12, f_C2S, parname_13, f_M2S, parname_14, f_M2P, parname_15, f_S2P, &
        !     parname_16, f_S2M, parname_17, f_P2M, parname_18, r_me,  parname_19, Q10pro, & 
        !     parname_20, kCH4,  parname_21, Omax,  parname_22, CH4_thre, parname_23, Tveg, &
        !     parname_24, Tpro_me, parname_25, Toxi, parname_26, f, parname_27, bubprob, &
        !     parname_28, Vmaxfraction

        namelist /nml_site_params/ wsmax,  wsmin,  extkU, &
            Tau_F, Tau_C, Tau_Micro, Tau_SlowSOM, Tau_Passive, &
            Q10rh, etaW, f_F2M, f_C2M, f_C2S, f_M2S, f_M2P, f_S2P, f_S2M, f_P2M, &  
            r_me, Q10pro, kCH4, Omax, CH4_thre, &  
            Tveg, Tpro_me, Toxi, f, bubprob, Vmaxfraction, &
            f_fast, f_slow, s_soil, par_shcap_snow, par_condu_snow, par_condu_b, & !par_albedo_snow, &
            par_fsub, par_rho_snow, par_decay_m


        real(8) :: def_LAImax(max_npft),   min_LAImax(max_npft),   max_LAImax(max_npft)
        real(8) :: def_LAImin(max_npft),   min_LAImin(max_npft),   max_LAImin(max_npft)
        real(8) :: def_stom_n(max_npft),   min_stom_n(max_npft),   max_stom_n(max_npft)
        real(8) :: def_SapS(max_npft),     min_SapS(max_npft),     max_SapS(max_npft)
        real(8) :: def_SapR(max_npft),     min_SapR(max_npft),     max_SapR(max_npft)
        real(8) :: def_SLAx(max_npft),     min_SLAx(max_npft),     max_SLAx(max_npft)
        real(8) :: def_GLmax(max_npft),    min_GLmax(max_npft),    max_GLmax(max_npft)
        real(8) :: def_GRmax(max_npft),    min_GRmax(max_npft),    max_GRmax(max_npft)
        real(8) :: def_Gsmax(max_npft),    min_Gsmax(max_npft),    max_Gsmax(max_npft)
        real(8) :: def_alpha(max_npft),    min_alpha(max_npft),    max_alpha(max_npft) 
        real(8) :: def_Vcmax0(max_npft),   min_Vcmax0(max_npft),   max_Vcmax0(max_npft)
        real(8) :: def_Ds0(max_npft),      min_Ds0(max_npft),      max_Ds0(max_npft)
        real(8) :: def_xfang(max_npft),    min_xfang(max_npft),    max_xfang(max_npft)
        real(8) :: def_rdepth(max_npft),   min_rdepth(max_npft),   max_rdepth(max_npft)
        real(8) :: def_Rootmax(max_npft),  min_Rootmax(max_npft),  max_Rootmax(max_npft)
        real(8) :: def_Stemmax(max_npft),  min_Stemmax(max_npft),  max_Stemmax(max_npft)
        real(8) :: def_Tau_Leaf(max_npft), min_Tau_Leaf(max_npft), max_Tau_Leaf(max_npft)
        real(8) :: def_Tau_Stem(max_npft), min_Tau_Stem(max_npft), max_Tau_Stem(max_npft)
        real(8) :: def_Tau_Root(max_npft), min_Tau_Root(max_npft), max_Tau_Root(max_npft) 
        real(8) :: def_Q10(max_npft),      min_Q10(max_npft),      max_Q10(max_npft)
        real(8) :: def_Rl0(max_npft),      min_Rl0(max_npft),      max_Rl0(max_npft)
        real(8) :: def_Rs0(max_npft),      min_Rs0(max_npft),      max_Rs0(max_npft)
        real(8) :: def_Rr0(max_npft),      min_Rr0(max_npft),      max_Rr0(max_npft)
        real(8) :: def_JV(max_npft),       min_JV(max_npft),       max_JV(max_npft)
        real(8) :: def_Entrpy(max_npft),   min_Entrpy(max_npft),   max_Entrpy(max_npft)
        real(8) :: def_gddonset(max_npft), min_gddonset(max_npft), max_gddonset(max_npft)
        real(8) :: def_hmax(max_npft),     min_hmax(max_npft),     max_hmax(max_npft)
        real(8) :: def_hl0(max_npft),      min_hl0(max_npft),      max_hl0(max_npft)
        real(8) :: def_LAIMAX0(max_npft),  min_LAIMAX0(max_npft),  max_LAIMAX0(max_npft)
        real(8) :: def_la0(max_npft),      min_la0(max_npft),      max_la0(max_npft)
        real(8) :: def_fn2l(max_npft),     min_fn2l(max_npft),     max_fn2l(max_npft)
        real(8) :: def_fn2r(max_npft),     min_fn2r(max_npft),     max_fn2r(max_npft)
        ! s_cLeaf, s_cStem,   s_cRoot,  s_nsc
        real(8) :: def_s_cLeaf(max_npft), min_s_cLeaf(max_npft), max_s_cLeaf(max_npft)
        real(8) :: def_s_cStem(max_npft), min_s_cStem(max_npft), max_s_cStem(max_npft)
        real(8) :: def_s_cRoot(max_npft), min_s_cRoot(max_npft), max_s_cRoot(max_npft)
        real(8) :: def_s_nsc(max_npft),   min_s_nsc(max_npft),   max_s_nsc(max_npft)
        real(8) :: def_s_nsn(max_npft),   min_s_nsn(max_npft),   max_s_nsn(max_npft)

        ! namelist /nml_species_params/ parname_1, min_LAImax, max_LAImax, parname_2, min_LAImin, max_LAImin, &
        !     parname_3, min_stom_n, max_stom_n, parname_4, min_SapS, max_SapS, parname_5, min_SapR, max_SapR,&
        !     parname_6, min_SLAx, max_SLAx, parname_7, min_GLmax, max_GLmax, parname_8, min_GRmax, max_GRmax,&
        !     parname_9, min_Gsmax, max_Gsmax, parname_10, min_alpha, max_alpha, parname_11, min_Vcmax0, max_Vcmax0,&
        !     parname_12, min_Ds0, max_Ds0, parname_13, min_xfang, max_xfang, parname_14, min_rdepth, max_rdepth,&
        !     parname_15, min_Rootmax, max_Rootmax, parname_16, min_Stemmax, max_Stemmax, parname_17, min_Tau_Leaf,max_Tau_Leaf,&
        !     parname_18, min_Tau_Stem, max_Tau_Stem, parname_19, min_Tau_Root, max_Tau_Root, parname_20, min_Q10, max_Q10,&
        !     parname_21, min_Rl0, max_Rl0, parname_22, min_Rs0, max_Rs0, parname_23, min_Rr0, max_Rr0, &
        !     parname_24, min_JV, max_JV, parname_25, min_Entrpy, max_Entrpy, parname_26, min_gddonset, max_gddonset,&
        !     parname_27, min_hmax, max_hmax, parname_28, min_hl0, max_hl0, parname_29, min_LAIMAX0, max_LAIMAX0,&
        !     parname_30, min_la0, max_la0
        namelist /nml_species_params/ def_LAImax, min_LAImax, max_LAImax, def_LAImin, min_LAImin, max_LAImin, &
            def_stom_n, min_stom_n, max_stom_n, def_SapS, min_SapS, max_SapS, def_SapR, min_SapR, max_SapR,&
            def_SLAx, min_SLAx, max_SLAx, def_GLmax, min_GLmax, max_GLmax, def_GRmax, min_GRmax, max_GRmax,&
            def_Gsmax, min_Gsmax, max_Gsmax, def_alpha, min_alpha, max_alpha, def_Vcmax0, min_Vcmax0, max_Vcmax0,&
            def_Ds0, min_Ds0, max_Ds0, def_xfang, min_xfang, max_xfang, def_rdepth, min_rdepth, max_rdepth,&
            def_Rootmax, min_Rootmax, max_Rootmax, def_Stemmax, min_Stemmax, max_Stemmax, &
            def_Tau_Leaf, min_Tau_Leaf,max_Tau_Leaf, def_Tau_Stem, min_Tau_Stem, max_Tau_Stem, &
            def_Tau_Root, min_Tau_Root, max_Tau_Root, def_Q10, min_Q10, max_Q10,&
            def_Rl0, min_Rl0, max_Rl0, def_Rs0, min_Rs0, max_Rs0, def_Rr0, min_Rr0, max_Rr0, &
            def_JV, min_JV, max_JV, def_Entrpy, min_Entrpy, max_Entrpy, def_gddonset, min_gddonset, max_gddonset,&
            def_hmax, min_hmax, max_hmax, def_hl0, min_hl0, max_hl0, def_LAIMAX0, min_LAIMAX0, max_LAIMAX0,&
            def_la0, min_la0, max_la0, def_fn2l, min_fn2l, max_fn2l, def_fn2r, min_fn2r, max_fn2r,&
            def_s_cLeaf, min_s_cLeaf, max_s_cLeaf, def_s_cStem, min_s_cStem, max_s_cStem, &
            def_s_cRoot, min_s_cRoot, max_s_cRoot, def_s_nsc, min_s_nsc, max_s_nsc, def_s_nsn, min_s_nsn, max_s_nsn

        nsite_params = 29+3+6
        nspec_params = 32+5
        npar = nsite_params + npft * nspec_params

        ! allocate(mc_in_params%max_sp_params(npft))
        ! allocate(mc_in_params%min_sp_params(npft))

        allocate(mc_params%sp(npft))
        allocate(mc_params%st(nsite_params))
        
        do ipft = 1, npft
            if (allocated(mc_params%sp))then
                allocate(mc_params%sp(ipft)%sp_params(nspec_params))
            endif
        enddo

        print *, "mcmc_configfile: ", adjustl(trim("configs/"))//adjustl(trim(in_mcmc_configfile))
        open(145, file=adjustl(trim("configs/"))//adjustl(trim(in_mcmc_configfile)))
        read(145, nml=nml_mcmc_settings, iostat=io)
        read(145, nml=nml_obsfiles,      iostat=io)
        read(145, nml=nml_site_params,   iostat=io)
        ! update the values of site-based MCMC parameters
        ! call update_site_params(st, &
        !     &   st%lat, st%lon,  st%wsmax,  st%wsmin,  st%extkU, &
        !     &   st%Tau_F, st%Tau_C, st%Tau_Micro, st%Tau_SlowSOM, st%Tau_Passive, &
        !     &   st%Q10rh, st%etaW, st%f_F2M, st%f_C2M, st%f_C2S, st%f_M2S, st%f_M2P, &
        !     &   st%f_S2P, st%f_S2M, st%f_P2M, &  
        !     &   st%r_me, st%Q10pro, st%kCH4, st%Omax, st%CH4_thre, &  
        !     &   st%Tveg, st%Tpro_me, st%Toxi, st%f, st%bubprob, st%Vmaxfraction)

        call update_mc_tile_params(mc_params%st(1),  "wsmax",        wsmax(1),       wsmax(2),       wsmax(3))
        call update_mc_tile_params(mc_params%st(2),  "wsmin",        wsmin(1),       wsmin(2),       wsmin(3))
        call update_mc_tile_params(mc_params%st(3),  "extkU",        extkU(1),       extkU(2),       extkU(3))
        call update_mc_tile_params(mc_params%st(4),  "Tau_F",        Tau_F(1),       Tau_F(2),       Tau_F(3))
        call update_mc_tile_params(mc_params%st(5),  "Tau_C",        Tau_C(1),       Tau_C(2),       Tau_C(3))
        call update_mc_tile_params(mc_params%st(6),  "Tau_Micro",    Tau_Micro(1),   Tau_Micro(2),   Tau_Micro(3))
        call update_mc_tile_params(mc_params%st(7),  "Tau_SlowSOM",  Tau_SlowSOM(1), Tau_SlowSOM(2), Tau_SlowSOM(3))
        call update_mc_tile_params(mc_params%st(8),  "Tau_Passive",  Tau_Passive(1), Tau_Passive(2), Tau_Passive(3))
        call update_mc_tile_params(mc_params%st(9),  "Q10rh",        Q10rh(1),       Q10rh(2),       Q10rh(3))
        call update_mc_tile_params(mc_params%st(10),  "etaW",         etaW(1),        etaW(2),        etaW(3))
        call update_mc_tile_params(mc_params%st(11), "f_F2M",        f_F2M(1),       f_F2M(2),       f_F2M(3))
        call update_mc_tile_params(mc_params%st(12), "f_C2M",        f_C2M(1),       f_C2M(2),       f_C2M(3))
        call update_mc_tile_params(mc_params%st(13), "f_C2S",        f_C2S(1),       f_C2S(2),       f_C2S(3))
        call update_mc_tile_params(mc_params%st(14), "f_M2S",        f_M2S(1),       f_M2S(2),       f_M2S(3))
        call update_mc_tile_params(mc_params%st(15), "f_M2P",        f_M2P(1),       f_M2P(2),       f_M2P(3))
        call update_mc_tile_params(mc_params%st(16), "f_S2P",        f_S2P(1),       f_S2P(2),       f_S2P(3))
        call update_mc_tile_params(mc_params%st(17), "f_S2M",        f_S2M(1),       f_S2M(2),       f_S2M(3))
        call update_mc_tile_params(mc_params%st(18), "f_P2M",        f_P2M(1),       f_P2M(2),       f_P2M(3))
        call update_mc_tile_params(mc_params%st(19), "r_me",         r_me(1),        r_me(2),        r_me(3))
        call update_mc_tile_params(mc_params%st(20), "Q10pro",       Q10pro(1),      Q10pro(2),      Q10pro(3))
        call update_mc_tile_params(mc_params%st(21), "kCH4",         kCH4(1),        kCH4(2),        kCH4(3))
        call update_mc_tile_params(mc_params%st(22), "Omax",         Omax(1),        Omax(2),        Omax(3))
        call update_mc_tile_params(mc_params%st(23), "CH4_thre",     CH4_thre(1),    CH4_thre(2),    CH4_thre(3))
        call update_mc_tile_params(mc_params%st(24), "Tveg",         Tveg(1),        Tveg(2),        Tveg(3))
        call update_mc_tile_params(mc_params%st(25), "Tpro_me",      Tpro_me(1),     Tpro_me(2),     Tpro_me(3))
        call update_mc_tile_params(mc_params%st(26), "Toxi",         Toxi(1),        Toxi(2),        Toxi(3))
        call update_mc_tile_params(mc_params%st(27), "f",            f(1),           f(2),           f(3))
        call update_mc_tile_params(mc_params%st(28), "bubprob",      bubprob(1),     bubprob(2),     bubprob(3))
        call update_mc_tile_params(mc_params%st(29), "Vmaxfraction", Vmaxfraction(1),Vmaxfraction(2),Vmaxfraction(3))
        ! f_fast, f_slow, s_soil 
        call update_mc_tile_params(mc_params%st(30), "f_fast", f_fast(1),f_fast(2),f_fast(3))
        call update_mc_tile_params(mc_params%st(31), "f_slow", f_slow(1),f_slow(2),f_slow(3))
        call update_mc_tile_params(mc_params%st(32), "s_soil", s_soil(1),s_soil(2),s_soil(3))
        ! par_shcap_snow, par_condu_snow, par_condu_b, & !par_albedo_snow, &
            ! par_fsub, par_rho_snow, par_decay_m
        call update_mc_tile_params(mc_params%st(33), "par_shcap_snow", par_shcap_snow(1),par_shcap_snow(2),par_shcap_snow(3))
        call update_mc_tile_params(mc_params%st(34), "par_condu_snow", par_condu_snow(1),par_condu_snow(2),par_condu_snow(3))
        call update_mc_tile_params(mc_params%st(35), "par_condu_b",    par_condu_b(1),   par_condu_b(2),par_condu_b(3))
        call update_mc_tile_params(mc_params%st(36), "par_fsub",       par_fsub(1),par_fsub(2),par_fsub(3))
        call update_mc_tile_params(mc_params%st(37), "par_rho_snow",   par_rho_snow(1),par_rho_snow(2),par_rho_snow(3))
        call update_mc_tile_params(mc_params%st(38), "par_decay_m",    par_decay_m(1),par_decay_m(2),par_decay_m(3))
        ! update the values of spec-based MCMC parameters
        read(145, nml=nml_species_params,   iostat=io)
        do ipft = 1, npft
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(1),  "LAImax",   &
                &   def_LAImax(ipft),   min_LAImax(ipft),   max_LAImax(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(2),  "LAImin",   &
                &   def_LAImin(ipft),   min_LAImin(ipft),   max_LAImin(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(3),  "stom_n",   &
                &   def_stom_n(ipft),   min_stom_n(ipft),   max_stom_n(ipft)) 
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(4),  "SapS",     &
                &   def_SapS(ipft),     min_SapS(ipft),     max_SapS(ipft)) 
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(5),  "SapR",     &
                &   def_SapR(ipft),     min_SapR(ipft),     max_SapR(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(6),  "SLA",      &
                &   def_SLAx(ipft),     min_SLAx(ipft),     max_SLAx(ipft)) 
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(7),  "GLmax",    &
                &   def_GLmax(ipft),    min_GLmax(ipft),    max_GLmax(ipft)) 
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(8),  "GRmax",    &
                &   def_GRmax(ipft),    min_GRmax(ipft),    max_GRmax(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(9),  "Gsmax",    &
                &   def_Gsmax(ipft),    min_Gsmax(ipft),    max_Gsmax(ipft)) 
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(10), "alpha",    &
                &   def_alpha(ipft),    min_alpha(ipft),    max_alpha(ipft)) 
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(11), "Vcmax0",   &
                &   def_Vcmax0(ipft),   min_Vcmax0(ipft),   max_Vcmax0(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(12), "Ds0",      &
                &   def_Ds0(ipft),      min_Ds0(ipft),      max_Ds0(ipft)) 
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(13), "xfang",    &
                &   def_xfang(ipft),    min_xfang(ipft),    max_xfang(ipft)) 
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(14), "rdepth",   &
                &   def_rdepth(ipft),   min_rdepth(ipft),   max_rdepth(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(15), "Rootmax",  &
                &   def_Rootmax(ipft),  min_Rootmax(ipft),  max_Rootmax(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(16), "Stemmax",  &
                &   def_Stemmax(ipft),  min_Stemmax(ipft),  max_Stemmax(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(17), "Tau_Leaf", &
                &   def_Tau_Leaf(ipft), min_Tau_Leaf(ipft), max_Tau_Leaf(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(18), "Tau_Stem", &
                &   def_Tau_Stem(ipft), min_Tau_Stem(ipft), max_Tau_Stem(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(19), "Tau_Root", &
                &   def_Tau_Root(ipft), min_Tau_Root(ipft), max_Tau_Root(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(20), "Q10",      &
                &   def_Q10(ipft),      min_Q10(ipft),      max_Q10(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(21), "Rl0",      &
                &   def_Rl0(ipft),      min_Rl0(ipft),      max_Rl0(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(22), "Rs0",      &
                &   def_Rs0(ipft),      min_Rs0(ipft),      max_Rs0(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(23), "Rr0",      &
                &   def_Rr0(ipft),      min_Rr0(ipft),      max_Rr0(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(24), "JV",       &
                &   def_JV(ipft),       min_JV(ipft),       max_JV(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(25), "Entrpy",   &
                &   def_Entrpy(ipft),   min_Entrpy(ipft),   max_Entrpy(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(26), "gddonset", &
                &   def_gddonset(ipft), min_gddonset(ipft), max_gddonset(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(27), "hmax",     &
                &   def_hmax(ipft),     min_hmax(ipft),     max_hmax(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(28), "hl0",      &
                &   def_hl0(ipft),      min_hl0(ipft),      max_hl0(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(29), "LAIMAX0",  &
                &   def_LAIMAX0(ipft),  min_LAIMAX0(ipft),  max_LAIMAX0(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(30), "la0",      &
                &   def_la0(ipft),      min_la0(ipft),      max_la0(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(31), "fn2l",      &
                &   def_fn2l(ipft),      min_fn2l(ipft),      max_fn2l(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(32), "fn2r",      &
                &   def_fn2r(ipft),      min_fn2r(ipft),      max_fn2r(ipft))
            ! s_cLeaf, s_cStem,   s_cRoot,  s_nsc
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(33), "s_cLeaf",      &
                &   def_s_cLeaf(ipft),   min_s_cLeaf(ipft), max_s_cLeaf(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(34), "s_cStem",      &
                &   def_s_cStem(ipft),   min_s_cStem(ipft), max_s_cStem(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(35), "s_cRoot",      &
                &   def_s_cRoot(ipft),   min_s_cRoot(ipft), max_s_cRoot(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(36), "s_nsc",      &
                &   def_s_nsc(ipft),     min_s_nsc(ipft),   max_s_nsc(ipft))
            call update_mc_tile_params(mc_params%sp(ipft)%sp_params(37), "s_nsn",      &
                &   def_s_nsn(ipft),     min_s_nsn(ipft),   max_s_nsn(ipft))
        enddo
        close(145)
        ! give the filepath to each variable
        vars4MCMC%ANPP_Shrub_y%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_ANPP_Shrub_y))
        vars4MCMC%ANPP_Tree_y%filepath  = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_ANPP_Tree_y))
        vars4MCMC%NPP_sphag_y%filepath  = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_NPP_sphag_y))
        vars4MCMC%BNPP_y%filepath       = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_BNPP_y))        ! tree + shrub
        vars4MCMC%er_d%filepath         = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_er_d))          ! shrub + sphag.
        vars4MCMC%er_h%filepath         = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_er_h))          ! shrub + sphag.
        vars4MCMC%gpp_d%filepath        = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_gpp_d))         ! Shrub + sphag.
        vars4MCMC%nee_d%filepath        = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_nee_d))         ! Shrub + sphag.
        vars4MCMC%nee_h%filepath        = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_nee_h))         ! shrub + sphag.
        vars4MCMC%LAI_d%filepath        = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_LAI_d))         ! tree  + Shrub
        vars4MCMC%rh_y%filepath         = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_rh_y))
        !
        vars4MCMC%leaf_mass_shrub_y%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_leaf_mass_shrub_y))
        vars4MCMC%stem_mass_shrub_y%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_stem_mass_shrub_y))
        vars4MCMC%cPlant_tree_y%filepath     = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_cPlant_tree_y))
        vars4MCMC%cPlant_sphag_y%filepath     = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_cPlant_sphag_y))
        vars4MCMC%cSoil_y%filepath     = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_cSoil_y))
        vars4MCMC%leaf_resp_shrub_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_leaf_resp_shrub_d))
        vars4MCMC%leaf_resp_tree_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_leaf_resp_tree_d)) 
        ! methane
        vars4MCMC%ch4_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_ch4_d))
        vars4MCMC%ch4_h%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_ch4_h))
        vars4MCMC%ch4_y%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_ch4_y))
        ! 
        vars4MCMC%CN_shag_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_CN_shag_d))
        vars4MCMC%photo_shrub_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_photo_shrub_d)) 
        vars4MCMC%photo_tree_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_photo_tree_d))
        !
        vars4MCMC%zwt_h%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_watertable_h))
        !
        vars4MCMC%lai_shrub_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_LAI_shrub_d))
        vars4MCMC%lai_tree_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_LAI_tree_d))
        vars4MCMC%photo_tree_h%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_photo_tree_h))
        !
        vars4MCMC%bnpp_tree_y%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_bnpp_tree_y))
        vars4MCMC%bnpp_shrub_y%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_bnpp_shrub_y))
        !
        vars4MCMC%gpp_tree_y%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_gpp_tree_y))
    end subroutine readMCMC_configs_NML

    subroutine update_mc_tile_params(in_mc_par, parname, parval, parmin, parmax)
        implicit none
        type(mc_tile_params), intent(inout) :: in_mc_par
        character(*), intent(in) :: parname
        real(8), intent(in) :: parmin, parmax, parval
        in_mc_par%parname = parname
        in_mc_par%parval  = parval
        in_mc_par%parmin  = parmin
        in_mc_par%parmax  = parmax
        return
    end subroutine update_mc_tile_params

    subroutine readObsData()
        implicit none
        call readObsData_var(vars4MCMC%ANPP_Shrub_y)
        call readObsData_var(vars4MCMC%ANPP_Tree_y)
        call readObsData_var(vars4MCMC%NPP_sphag_y)
        call readObsData_var(vars4MCMC%BNPP_y)        ! tree + shrub
        call readObsData_var(vars4MCMC%er_d)          ! shrub + sphag.
        call readObsData_var(vars4MCMC%er_h)          ! shrub + sphag.
        call readObsData_var(vars4MCMC%gpp_d)         ! Shrub + sphag.
        call readObsData_var(vars4MCMC%nee_d)         ! Shrub + sphag.
        call readObsData_var(vars4MCMC%nee_h)         ! shrub + sphag.
        call readObsData_var(vars4MCMC%LAI_d)         ! tree  + Shrub
        call readObsData_var(vars4MCMC%rh_y)
        !
        call readObsData_var(vars4MCMC%leaf_mass_shrub_y)
        call readObsData_var(vars4MCMC%stem_mass_shrub_y)
        call readObsData_var(vars4MCMC%cPlant_tree_y)
        call readObsData_var(vars4MCMC%cPlant_sphag_y)
        call readObsData_var(vars4MCMC%cSoil_y)
        call readObsData_var(vars4MCMC%leaf_resp_shrub_d) 
        call readObsData_var(vars4MCMC%leaf_resp_tree_d) 
        ! methane
        call readObsData_var(vars4MCMC%ch4_d) 
        call readObsData_var(vars4MCMC%ch4_h) 
        call readObsData_var(vars4MCMC%ch4_y) 
        ! 
        call readObsData_var(vars4MCMC%CN_shag_d) 
        call readObsData_var(vars4MCMC%photo_shrub_d) 
        call readObsData_var(vars4MCMC%photo_tree_d) 
        !
        call readObsData_var(vars4MCMC%zwt_h)
        !
        call readObsData_var(vars4MCMC%lai_shrub_d)
        call readObsData_var(vars4MCMC%lai_tree_d)
        call readObsData_var(vars4MCMC%photo_tree_h)
        !
        call readObsData_var(vars4MCMC%bnpp_tree_y)
        call readObsData_var(vars4MCMC%bnpp_shrub_y)
        !
        call readObsData_var(vars4MCMC%gpp_tree_y)
    end subroutine readObsData

    subroutine readObsData_var(var_obsData)
        implicit none
        type(interCostVariable), intent(inout) :: var_obsData
        logical toExistOrNot
        integer toCountLines

        INQUIRE(FILE=var_obsData%filepath, EXIST=toExistOrNot)
        var_obsData%existOrNot = toExistOrNot
        if (var_obsData%existOrNot) then
            call ReadLineNumFromFile(var_obsData%filepath, toCountLines)
            allocate(var_obsData%obsData(toCountLines, 5))
            call ReadObsDataFromFile(var_obsData%filepath, toCountLines, var_obsData%obsData)
            allocate(var_obsData%mdData(toCountLines, 4))
        endif
        return
    end subroutine readObsData_var

    subroutine mc_update_mc_params()
        implicit none
        integer :: ipar, ipft, ipar4DA
        ipar4DA = 0
        do ipar = 1, npar4st
            ipar4DA = ipar4DA + 1
            call update_mc_tile_params(mc_params%st(mc_DApar%DAparidx_st(ipar)),&  
                &                      mc_params%st(mc_DApar%DAparidx_st(ipar))%parname, &
                &                      mc_DApar%DApar(ipar4DA), mc_DApar%DAparmin(ipar4DA), mc_DApar%DAparmax(ipar4DA))
        enddo
        do ipft = 1, npft
            do ipar = 1, npar4sp
                ipar4DA = ipar4DA + 1
                call update_mc_tile_params(mc_params%sp(ipft)%sp_params(mc_DApar%DAparidx_sp(ipft, ipar)),&  
                    &                      mc_params%sp(ipft)%sp_params(mc_DApar%DAparidx_sp(ipft, ipar))%parname, &
                    &                      mc_DApar%DApar(ipar4DA), mc_DApar%DAparmin(ipar4DA), mc_DApar%DAparmax(ipar4DA))
            enddo
        enddo
    end subroutine mc_update_mc_params

    subroutine mc_update_params4simu()
        implicit none
        integer :: imc_params, ipft
        ! mc_params -> in_params_vals
        do imc_params = 1, size(mc_params%st)
            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "wsmax") then
                ! print*, "check: ", mc_params%st(imc_params)%parname, adjustl(trim(mc_params%st(imc_params)%parname))
                in_params_vals%st_params%wsmax = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "wsmin") then    
                in_params_vals%st_params%wsmin = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "extkU") then      
                in_params_vals%st_params%extkU = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tau_F")  then     
                in_params_vals%st_params%Tau_F = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tau_C")  then     
                in_params_vals%st_params%Tau_C = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tau_Micro") then  
                in_params_vals%st_params%Tau_Micro = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tau_SlowSOM") then
                in_params_vals%st_params%Tau_SlowSOM = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tau_Passive") then
                in_params_vals%st_params%Tau_Passive = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Q10rh")  then     
                in_params_vals%st_params%Q10rh = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "etaW")  then      
                in_params_vals%st_params%etaW = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_F2M")  then     
                in_params_vals%st_params%f_F2M = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_C2M")  then     
                in_params_vals%st_params%f_C2M = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_C2S")  then      
                in_params_vals%st_params%f_C2S = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_M2S")  then       
                in_params_vals%st_params%f_M2S = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_M2P")  then       
                in_params_vals%st_params%f_M2P = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_S2P")  then       
                in_params_vals%st_params%f_S2P = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_S2M")  then       
                in_params_vals%st_params%f_S2M = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_P2M")  then       
                in_params_vals%st_params%f_P2M = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "r_me")  then         
                in_params_vals%st_params%r_me = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Q10pro")  then       
                in_params_vals%st_params%Q10pro = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "kCH4")  then         
                in_params_vals%st_params%kCH4 = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Omax")   then        
                in_params_vals%st_params%Omax = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "CH4_thre")  then     
                in_params_vals%st_params%CH4_thre = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tveg")  then         
                in_params_vals%st_params%Tveg = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tpro_me")  then      
                in_params_vals%st_params%Tpro_me = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Toxi")  then         
                in_params_vals%st_params%Toxi = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f")   then           
                in_params_vals%st_params%f = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "bubprob") then     
                in_params_vals%st_params%bubprob = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Vmaxfraction")  then 
                in_params_vals%st_params%Vmaxfraction = mc_params%st(imc_params)%parval
            endif

            ! f_fast, f_slow, s_soil 
            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_fast")  then 
                in_params_vals%st_params%f_fast = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_slow")  then 
                in_params_vals%st_params%f_slow = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "s_soil")  then 
                in_params_vals%st_params%s_soil = mc_params%st(imc_params)%parval
            endif

            ! par_shcap_snow, par_condu_snow, par_condu_b, & !par_albedo_snow, &
            ! par_fsub, par_rho_snow, par_decay_m
            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "par_shcap_snow")  then 
                in_params_vals%st_params%par_shcap_snow = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "par_condu_snow")  then 
                in_params_vals%st_params%par_condu_snow = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "par_condu_b")  then 
                in_params_vals%st_params%par_condu_b = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "par_fsub")  then 
                in_params_vals%st_params%par_fsub = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "par_rho_snow")  then 
                in_params_vals%st_params%par_rho_snow = mc_params%st(imc_params)%parval
            endif

            if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "par_decay_m")  then 
                in_params_vals%st_params%par_decay_m = mc_params%st(imc_params)%parval
            endif
        enddo

        do ipft = 1, npft
            do imc_params = 1, size(mc_params%sp(ipft)%sp_params)
                ! spec_name  = st%sp(ipft)%spec_name
                ! pft_weight = st%sp(ipft)%pft_weight
                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "stom_n") then
                    in_params_vals%sp_params(ipft)%stom_n = int(mc_params%sp(ipft)%sp_params(imc_params)%parval)
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)).eq. "LAImax") then
                    in_params_vals%sp_params(ipft)%LAImax = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "LAImin") then
                    in_params_vals%sp_params(ipft)%LAImin = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "SapS") then
                    in_params_vals%sp_params(ipft)%SapS   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "SapR") then 
                    in_params_vals%sp_params(ipft)%SapR   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "SLA") then   
                    in_params_vals%sp_params(ipft)%SLA   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "GLmax")  then 
                    in_params_vals%sp_params(ipft)%GLmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "GRmax") then  
                    in_params_vals%sp_params(ipft)%GRmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Gsmax") then   
                    in_params_vals%sp_params(ipft)%Gsmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "alpha")then
                    in_params_vals%sp_params(ipft)%alpha   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Vcmax0")then  
                    in_params_vals%sp_params(ipft)%Vcmax0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Ds0")then
                    in_params_vals%sp_params(ipft)%Ds0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "xfang")then    
                    in_params_vals%sp_params(ipft)%xfang   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "rdepth")then    
                    in_params_vals%sp_params(ipft)%rdepth   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Rootmax")then    
                    in_params_vals%sp_params(ipft)%Rootmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Stemmax")then    
                    in_params_vals%sp_params(ipft)%Stemmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Tau_Leaf") then   
                    in_params_vals%sp_params(ipft)%Tau_Leaf   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Tau_Stem")then    
                    in_params_vals%sp_params(ipft)%Tau_Stem   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif
                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Tau_Root") then   
                    in_params_vals%sp_params(ipft)%Tau_Root   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Q10") then    
                    in_params_vals%sp_params(ipft)%Q10   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Rl0") then   
                    in_params_vals%sp_params(ipft)%Rl0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Rs0") then   
                    in_params_vals%sp_params(ipft)%Rs0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Rr0") then    
                    in_params_vals%sp_params(ipft)%Rr0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif
                
                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "JV")  then  
                    in_params_vals%sp_params(ipft)%JV   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Entrpy") then   
                    in_params_vals%sp_params(ipft)%Entrpy   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "gddonset") then   
                    in_params_vals%sp_params(ipft)%gddonset   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "hmax")  then  
                    in_params_vals%sp_params(ipft)%hmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "hl0")  then  
                    in_params_vals%sp_params(ipft)%hl0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "LAIMAX0") then   
                    in_params_vals%sp_params(ipft)%LAIMAX0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "la0")  then  
                    in_params_vals%sp_params(ipft)%la0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "fn2l")  then  
                    in_params_vals%sp_params(ipft)%fn2l   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "fn2r")  then  
                    in_params_vals%sp_params(ipft)%fn2r   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                ! s_cLeaf, s_cStem,   s_cRoot,  s_nsc
                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "s_cLeaf")  then  
                    in_params_vals%sp_params(ipft)%s_cLeaf   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "s_cStem")  then  
                    in_params_vals%sp_params(ipft)%s_cStem   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "s_cRoot")  then  
                    in_params_vals%sp_params(ipft)%s_cRoot   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "s_nsc")  then  
                    in_params_vals%sp_params(ipft)%s_nsc   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif

                if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "s_nsn")  then  
                    in_params_vals%sp_params(ipft)%s_nsn   = mc_params%sp(ipft)%sp_params(imc_params)%parval
                endif
            enddo
        enddo
    end subroutine mc_update_params4simu

    ! subroutine mc_update_params4simu(st)
    !     implicit none
    !     type(site_data_type), intent(inout) :: st
    !     integer :: ipft, imc_params

    !     real(8) :: lat, lon,  wsmax,  wsmin,  extkU
    !     real(8) :: Tau_F, Tau_C, Tau_Micro, Tau_SlowSOM, Tau_Passive
    !     real(8) :: Q10rh, etaW, f_F2M, f_C2M, f_C2S, f_M2S, f_M2P, f_S2P, f_S2M, f_P2M  
    !     real(8) :: r_me, Q10pro, kCH4, Omax, CH4_thre  
    !     real(8) :: Tveg, Tpro_me, Toxi, f, bubprob, Vmaxfraction 

    !     integer :: stom_n
    !     character(50) :: spec_name
    !     real(8) :: pft_weight, LAImax, LAImin 
    !     real(8) :: SapS,     SapR,     SLA
    !     real(8) :: GLmax,    GRmax,    Gsmax
    !     real(8) :: alpha,    Vcmax0,   Ds0
    !     real(8) :: xfang,    rdepth,   Rootmax,  Stemmax
    !     real(8) :: Tau_Leaf, Tau_Stem, Tau_Root
    !     real(8) :: Q10,      Rl0,      Rs0,      Rr0
    !     real(8) :: JV,       Entrpy,   gddonset 
    !     real(8) ::  hmax,    hl0,      LAIMAX0,  la0

    !     lat = st%lat
    !     lon = st%lon
    !     do imc_params = 1, size(mc_params%st)
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "wsmax") then
    !             print*, "check: ", mc_params%st(imc_params)%parname, adjustl(trim(mc_params%st(imc_params)%parname))
    !             wsmax = mc_params%st(imc_params)%parval
    !         endif
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "wsmin")     wsmin = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "extkU")       extkU = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tau_F")       Tau_F = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tau_C")       Tau_C = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tau_Micro")   Tau_Micro = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tau_SlowSOM") Tau_SlowSOM = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tau_Passive") Tau_Passive = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Q10rh")       Q10rh = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "etaW")        etaW = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_F2M")       f_F2M = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_C2M")       f_C2M = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_C2S")       f_C2S = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_M2S")       f_M2S = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_M2P")       f_M2P = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_S2P")       f_S2P = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_S2M")       f_S2M = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f_P2M")       f_P2M = mc_params%st(imc_params)%parval

    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "r_me")         r_me = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Q10pro")       Q10pro = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "kCH4")         kCH4 = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Omax")         Omax = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "CH4_thre")     CH4_thre = mc_params%st(imc_params)%parval

    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tveg")         Tveg = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Tpro_me")      Tpro_me = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Toxi")         Toxi = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "f")            f = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "bubprob")      bubprob = mc_params%st(imc_params)%parval
    !         if(adjustl(trim(mc_params%st(imc_params)%parname)) .eq. "Vmaxfraction") Vmaxfraction = mc_params%st(imc_params)%parval
    !     enddo
  
    !     print *, lat, lon,  wsmax,  wsmin,  extkU, &
    !     & Tau_F, Tau_C, Tau_Micro, Tau_SlowSOM, Tau_Passive, &
    !     & Q10rh, etaW, f_F2M, f_C2M, f_C2S, f_M2S, f_M2P, f_S2P, f_S2M, f_P2M, &  
    !     & r_me, Q10pro, kCH4, Omax, CH4_thre, &  
    !     & Tveg, Tpro_me, Toxi, f, bubprob, Vmaxfraction
    !     print*, "update_Tauc: ", Tau_C
    !     call update_site_params(st, &
    !         & lat, lon,  wsmax,  wsmin,  extkU, &
    !         & Tau_F, Tau_C, Tau_Micro, Tau_SlowSOM, Tau_Passive, &
    !         & Q10rh, etaW, f_F2M, f_C2M, f_C2S, f_M2S, f_M2P, f_S2P, f_S2M, f_P2M, &  
    !         & r_me, Q10pro, kCH4, Omax, CH4_thre, &  
    !         & Tveg, Tpro_me, Toxi, f, bubprob, Vmaxfraction)

    !     ! update the species parameters


    !     do ipft = 1, npft
    !         do imc_params = 1, size(mc_params%sp(ipft)%sp_params)
    !             spec_name  = st%sp(ipft)%spec_name
    !             pft_weight = st%sp(ipft)%pft_weight

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "stom_n") then
    !                 stom_n = int(mc_params%sp(ipft)%sp_params(imc_params)%parval)
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)).eq. "LAImax") then
    !                 LAImax = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "LAImin") then
    !                 LAImin = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "SapS") then
    !                 SapS   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "SapR") then 
    !                 SapS   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "SLA") then   
    !                 SapS   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "GLmax")  then 
    !                 GLmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "GRmax") then  
    !                 GRmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Gsmax") then   
    !                 Gsmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "alpha")then
    !                 alpha   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Vcmax0")then  
    !                 Vcmax0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Ds0")then
    !                 Ds0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "xfang")then    
    !                 xfang   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "rdepth")then    
    !                 rdepth   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Rootmax")then    
    !                 Rootmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Stemmax")then    
    !                 Stemmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Tau_Leaf") then   
    !                 Tau_Leaf   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Tau_Stem")then    
    !                 Tau_Stem   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif
    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Tau_Root") then   
    !                 Tau_Root   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Q10") then    
    !                 Q10   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Rl0") then   
    !                 Rl0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Rs0") then   
    !                 Rs0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Rr0") then    
    !                 Rr0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif
                
    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "JV")  then  
    !                 JV   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "Entrpy") then   
    !                 Entrpy   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "gddonset") then   
    !                 gddonset   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "hmax")  then  
    !                 hmax   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif
    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "hl0")  then  
    !                 hl0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif
    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "LAIMAX0") then   
    !                 LAIMAX0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif

    !             if(adjustl(trim(mc_params%sp(ipft)%sp_params(imc_params)%parname)) .eq. "la0")  then  
    !                 la0   = mc_params%sp(ipft)%sp_params(imc_params)%parval
    !             endif
    !         enddo

    !         call update_spec_params(st%sp(ipft), &
    !             &   spec_name, stom_n,  pft_weight,  LAImax, LAImin, & 
    !             SapS,     SapR,      SLA, GLmax, GRmax,  Gsmax, &
    !             alpha,    Vcmax0,    Ds0,  xfang, rdepth, Rootmax,  Stemmax, &
    !             Tau_Leaf, Tau_Stem,  Tau_Root,    Q10,    Rl0,      Rs0,      Rr0, &
    !             JV,       Entrpy,    gddonset,    hmax,   hl0,      LAIMAX0,  la0)
    !     enddo
    ! end subroutine mc_update_params4simu

    ! subroutine renewMDpars(parval, re_in_params)
    !     implicit none
    !     real(8), intent(in) :: parval(:)
    !     type(nml_params_data_type), intent(inout) :: re_in_params

    !     re_in_params%lat         = parval(1)
    !     re_in_params%lon         = parval(2)
    !     re_in_params%wsmax       = parval(3)
    !     re_in_params%wsmin       = parval(4)                                            
    !     re_in_params%LAIMAX      = parval(5)
    !     re_in_params%LAIMIN      = parval(6)
    !     re_in_params%rdepth      = parval(7)
    !     re_in_params%Rootmax     = parval(8)
    !     re_in_params%Stemmax     = parval(9)                                    
    !     re_in_params%SapR        = parval(10)
    !     re_in_params%SapS        = parval(11)
    !     re_in_params%SLAx        = parval(12)
    !     re_in_params%GLmax       = parval(13)
    !     re_in_params%GRmax       = parval(14)
    !     re_in_params%Gsmax       = parval(15)
    !     re_in_params%stom_n      = parval(16)         
    !     re_in_params%a1          = parval(17)
    !     re_in_params%Ds0         = parval(18)
    !     re_in_params%Vcmax0      = parval(19)
    !     re_in_params%extkU       = parval(20)
    !     re_in_params%xfang       = parval(21)
    !     re_in_params%alpha       = parval(22)    
    !     re_in_params%Tau_Leaf    = parval(23)
    !     re_in_params%Tau_Stem    = parval(24)
    !     re_in_params%Tau_Root    = parval(25)
    !     re_in_params%Tau_F       = parval(26)
    !     re_in_params%Tau_C       = parval(27)
    !     re_in_params%Tau_Micro   = parval(28)
    !     re_in_params%Tau_SlowSOM = parval(29)
    !     re_in_params%Tau_Passive = parval(30)    
    !     re_in_params%gddonset    = parval(31)
    !     re_in_params%Q10         = parval(32)
    !     re_in_params%Rl0         = parval(33)     
    !     re_in_params%Rs0         = parval(34)    
    !     re_in_params%Rr0         = parval(35)                    
    !     re_in_params%r_me        = parval(36)
    !     re_in_params%Q10pro      = parval(37)
    !     re_in_params%kCH4        = parval(38)
    !     re_in_params%Omax         = parval(39)
    !     re_in_params%CH4_thre     = parval(40)
    !     re_in_params%Tveg         = parval(41)
    !     re_in_params%Tpro_me      = parval(42)
    !     re_in_params%Toxi         = parval(43)        
    !     re_in_params%f            = parval(44)
    !     re_in_params%bubprob      = parval(45)
    !     re_in_params%Vmaxfraction = parval(46)                                    
    !     re_in_params%Q10rh        = parval(47)
    !     re_in_params%JV           = parval(48)
    !     re_in_params%Entrpy       = parval(49)                    
    !     re_in_params%etaL         = parval(50)
    !     re_in_params%etaW         = parval(51)
    !     re_in_params%etaR         = parval(52)
    !     re_in_params%f_F2M        = parval(53)
    !     re_in_params%f_C2M        = parval(54)
    !     re_in_params%f_C2S        = parval(55)
    !     re_in_params%f_M2S        = parval(56)
    !     re_in_params%f_M2P        = parval(57)
    !     re_in_params%f_S2P        = parval(58)
    !     re_in_params%f_S2M        = parval(59)
    !     re_in_params%f_P2M        = parval(60)
    !     return
    ! end subroutine renewMDpars


    ! subroutine giveValues2var(filepath, existOrNot, data)
    !     implicit none
    !     character(500) filepath
    !     logical existOrNot
    !     real(8), allocatable :: data(:, :)
    !     integer count_lines

    !     INQUIRE(FILE=filepath, EXIST=existOrNot)
    !     if(existOrNot)then
    !         call ReadLineNumFromFile(filepath, count_lines)
    !         allocate(data(count_lines, 5))
    !         call ReadObsDataFromFile(filepath, count_lines, data)
    !     end if
    !     return
    ! end subroutine giveValues2var

    ! subroutine giveValues2par(arr_par)
    !     implicit none
    !     real(8), intent(inout) :: arr_par(:)

    !     arr_par(1)  = mc_lat
    !     arr_par(2)  = mc_Longitude 
    !     arr_par(3)  = mc_wsmax 
    !     arr_par(4)  = mc_wsmin                                                      
    !     arr_par(5)  = mc_LAIMAX
    !     arr_par(6)  = mc_LAIMIN    
    !     arr_par(7)  = mc_rdepth    
    !     arr_par(8)  = mc_Rootmax    
    !     arr_par(9)  = mc_Stemmax                                            
    !     arr_par(10) = mc_SapR    
    !     arr_par(11) = mc_SapS     
    !     arr_par(12) = mc_SLA        
    !     arr_par(13) = mc_GLmax    
    !     arr_par(14) = mc_GRmax    
    !     arr_par(15) = mc_Gsmax    
    !     arr_par(16) = mc_stom_n                                            
    !     arr_par(17) = mc_a1       
    !     arr_par(18) = mc_Ds0        
    !     arr_par(19) = mc_Vcmx0    
    !     arr_par(20) = mc_extkU    
    !     arr_par(21) = mc_xfang    
    !     arr_par(22) = mc_alpha                         
    !     arr_par(23) = mc_Tau_Leaf   
    !     arr_par(24) = mc_Tau_Stem   
    !     arr_par(25) = mc_Tau_Root   
    !     arr_par(26) = mc_Tau_F       
    !     arr_par(27) = mc_Tau_C       
    !     arr_par(28) = mc_Tau_Micro   
    !     arr_par(29) = mc_Tau_SlowSOM 
    !     arr_par(30) = mc_Tau_Passive                             
    !     arr_par(31) = mc_gddonset    
    !     arr_par(32) = mc_Q10         
    !     arr_par(33) = mc_Rl0        
    !     arr_par(34) = mc_Rs0        
    !     arr_par(35) = mc_Rr0                            
    !     arr_par(36) = mc_r_me   
    !     arr_par(37) = mc_Q10pro   
    !     arr_par(38) = mc_kCH4    
    !     arr_par(39) = mc_Omax   
    !     arr_par(40) = mc_CH4_thre 
    !     arr_par(41) = mc_Tveg  
    !     arr_par(42) = mc_Tpro_me 
    !     arr_par(43) = mc_Toxi               
    !     arr_par(44) = mc_f    
    !     arr_par(45) = mc_bubprob  
    !     arr_par(46) = mc_Vmaxfraction                                        
    !     arr_par(47) = mc_Q10rh  
    !     arr_par(48) = mc_JV   
    !     arr_par(49) = mc_Entrpy                                
    !     arr_par(50) = mc_etaL   
    !     arr_par(51) = mc_etaW  
    !     arr_par(52) = mc_etaR   
    !     arr_par(53) = mc_f_F2M   
    !     arr_par(54) = mc_f_C2M  
    !     arr_par(55) = mc_f_C2S 
    !     arr_par(56) = mc_f_M2S  
    !     arr_par(57) = mc_f_M2P 
    !     arr_par(58) = mc_f_S2P  
    !     arr_par(59) = mc_f_S2M  
    !     arr_par(60) = mc_f_P2M 

    ! end subroutine giveValues2par

    subroutine GetSimuData(get_iyear, get_iday, get_ihour, in_st, nHr, nDay, nMon, nYr)
        implicit none
        type(site_data_type), intent(in) :: in_st
        integer, intent(in) :: nHr, nDay, nMon, nYr

        integer get_iyear, get_iday, get_ihour
        integer i, ipft
        ! vars4MCMC%
        mc_iyear = get_iyear
        mc_iday  = get_iday
        mc_ihour = get_ihour + 1

        ! ------------------------------- with weight ---------------------------------------------------------------
        ! call GetSimuData_var(vars4MCMC%ANPP_Shrub_y,(outVars_y%sp(2)%nppLeaf + outVars_y%sp(2)%nppStem)*24*365)
        ! call GetSimuData_var(vars4MCMC%ANPP_Tree_y, (outVars_y%sp(1)%nppLeaf + outVars_y%sp(1)%nppStem)*24*365)
        ! call GetSimuData_var(vars4MCMC%NPP_sphag_y,  outVars_y%sp(3)%npp*24*365)
        ! call GetSimuData_var(vars4MCMC%BNPP_y,       &
        !     (in_st%sp(1)%pft_weight/(in_st%sp(1)%pft_weight + in_st%sp(2)%pft_weight)*outVars_y%sp(1)%nppRoot + &
        !      in_st%sp(2)%pft_weight/(in_st%sp(1)%pft_weight + in_st%sp(2)%pft_weight)*outVars_y%sp(2)%nppRoot)*24*365)        ! tree + shrub
        ! call GetSimuData_var(vars4MCMC%er_d,        &
        !     (in_st%sp(2)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_d%sp(2)%ra + &
        !      in_st%sp(3)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_d%sp(3)%ra + outVars_d%rh)*24)          ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%er_h,        &
        !     (in_st%sp(2)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_h%sp(2)%ra + &
        !      in_st%sp(3)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_h%sp(3)%ra + outVars_h%rh))          ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%gpp_d,       & ! change to hourly
        !     (in_st%sp(2)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_h%sp(2)%gpp + &
        !      in_st%sp(3)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_h%sp(3)%gpp))         ! Shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%nee_d,       &
        !     ((in_st%sp(2)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_d%sp(2)%ra  + &
        !       in_st%sp(3)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_d%sp(3)%ra  + outVars_d%rh)-&
        !      (in_st%sp(2)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_d%sp(2)%gpp + &
        !       in_st%sp(3)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_d%sp(3)%gpp))*24)         ! Shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%nee_h,       &
        !     ((in_st%sp(2)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_h%sp(2)%ra  + &
        !       in_st%sp(3)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_h%sp(3)%ra  + outVars_h%rh) -&
        !      (in_st%sp(2)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_h%sp(2)%gpp + &
        !       in_st%sp(3)%pft_weight/(in_st%sp(2)%pft_weight + in_st%sp(3)%pft_weight)*outVars_h%sp(3)%gpp)))         ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%LAI_d,       &
        !     (in_st%sp(1)%pft_weight/(in_st%sp(1)%pft_weight + in_st%sp(2)%pft_weight)*outVars_d%sp(1)%LAI + &
        !      in_st%sp(2)%pft_weight/(in_st%sp(1)%pft_weight + in_st%sp(2)%pft_weight)*outVars_d%sp(2)%LAI))         ! tree  + Shrub
        ! ! !
        ! call GetSimuData_var(vars4MCMC%leaf_mass_shrub_y, outVars_y%sp(2)%cleaf)
        ! call GetSimuData_var(vars4MCMC%stem_mass_shrub_y, outVars_y%sp(2)%cStem)

        ! call GetSimuData_var(vars4MCMC%cPlant_tree_y, outVars_d%sp(1)%cLeaf+outVars_d%sp(1)%cStem)
        ! call GetSimuData_var(vars4MCMC%cPlant_sphag_y, outVars_d%sp(3)%cLeaf+outVars_d%sp(3)%cStem + outVars_d%sp(3)%cRoot)
        ! call GetSimuData_var(vars4MCMC%cSoil_y, outVars_d%cSoilFast + outVars_d%cSoilSlow + outVars_d%cSoilPassive)

        ! call GetSimuData_var(vars4MCMC%leaf_resp_shrub_d, outVars_d%sp(2)%raLeaf*24) 
        ! call GetSimuData_var(vars4MCMC%leaf_resp_tree_d,  outVars_d%sp(1)%raLeaf*24) 
        ! ! methane
        ! call GetSimuData_var(vars4MCMC%ch4_d, outVars_d%wetlandCH4*24) 
        ! call GetSimuData_var(vars4MCMC%ch4_h, outVars_h%wetlandCH4) 
        ! ! 
        ! call GetSimuData_var(vars4MCMC%CN_shag_d,    ((outVars_d%sp(3)%cleaf + outVars_d%sp(3)%cStem + &
        !                                               outVars_d%sp(3)%cRoot)/(outVars_d%sp(3)%nLeaf + &
        !                                               outVars_d%sp(3)%nStem + outVars_d%sp(3)%nRoot))) 
        ! call GetSimuData_var(vars4MCMC%photo_shrub_d, outVars_d%sp(2)%gpp*24) 
        ! call GetSimuData_var(vars4MCMC%photo_tree_d,  outVars_d%sp(1)%gpp*24) 
        ! ----------------------------------------------------------------------------
        ! call GetSimuData_var(vars4MCMC%ANPP_Shrub_y,(outVars_y%sp(2)%nppLeaf + outVars_y%sp(2)%nppStem)*24*365)
        ! call GetSimuData_var(vars4MCMC%ANPP_Tree_y, (outVars_y%sp(1)%nppLeaf + outVars_y%sp(1)%nppStem)*24*365)
        ! call GetSimuData_var(vars4MCMC%NPP_sphag_y,  outVars_y%sp(3)%npp*24*365)
        ! call GetSimuData_var(vars4MCMC%BNPP_y,      (outVars_y%sp(1)%nppRoot + outVars_y%sp(2)%nppRoot)*24*365)        ! tree + shrub
        ! call GetSimuData_var(vars4MCMC%er_d,        (outVars_d%sp(2)%ra + outVars_d%sp(3)%ra + outVars_d%rh)*24)          ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%er_h,        (outVars_h%sp(2)%ra + outVars_h%sp(3)%ra + outVars_h%rh))          ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%gpp_d,       & ! change to hourly
        !     (outVars_h%sp(2)%gpp + outVars_h%sp(3)%gpp))         ! Shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%nee_d,       &
        !     ((outVars_d%sp(2)%ra  + outVars_d%sp(3)%ra  + outVars_d%rh)-&
        !      (outVars_d%sp(2)%gpp + outVars_d%sp(3)%gpp))*24)         ! Shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%nee_h,       &
        !     ((outVars_h%sp(2)%ra  + outVars_h%sp(3)%ra  + outVars_h%rh) -&
        !      (outVars_h%sp(2)%gpp + outVars_h%sp(3)%gpp)))         ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%LAI_d,       &
        !     (outVars_d%sp(1)%LAI + outVars_d%sp(2)%LAI))         ! tree  + Shrub
        ! ! !
        ! call GetSimuData_var(vars4MCMC%leaf_mass_shrub_y, outVars_y%sp(2)%cleaf)
        ! call GetSimuData_var(vars4MCMC%stem_mass_shrub_y, outVars_y%sp(2)%cStem)

        ! call GetSimuData_var(vars4MCMC%cPlant_tree_y, outVars_d%sp(1)%cLeaf+outVars_d%sp(1)%cStem)
        ! call GetSimuData_var(vars4MCMC%cPlant_sphag_y, outVars_d%sp(3)%cLeaf+outVars_d%sp(3)%cStem + outVars_d%sp(3)%cRoot)
        ! call GetSimuData_var(vars4MCMC%cSoil_y, outVars_d%cSoilFast + outVars_d%cSoilSlow + outVars_d%cSoilPassive)

        ! call GetSimuData_var(vars4MCMC%leaf_resp_shrub_d, outVars_d%sp(2)%raLeaf*24) 
        ! call GetSimuData_var(vars4MCMC%leaf_resp_tree_d,  outVars_d%sp(1)%raLeaf*24) 
        ! ! methane
        ! call GetSimuData_var(vars4MCMC%ch4_d, outVars_d%wetlandCH4*24) 
        ! call GetSimuData_var(vars4MCMC%ch4_h, outVars_h%wetlandCH4) 
        ! ! 
        ! call GetSimuData_var(vars4MCMC%CN_shag_d,    ((outVars_d%sp(3)%cleaf + outVars_d%sp(3)%cStem + &
        !                                               outVars_d%sp(3)%cRoot)/(outVars_d%sp(3)%nLeaf + &
        !                                               outVars_d%sp(3)%nStem + outVars_d%sp(3)%nRoot))) 
        ! call GetSimuData_var(vars4MCMC%photo_shrub_d, outVars_d%sp(2)%gpp*24) 
        ! call GetSimuData_var(vars4MCMC%photo_tree_d,  outVars_d%sp(1)%gpp*24) 
        ! ----------------------------------------------------------------------------
        ! call GetSimuData_var(vars4MCMC%ANPP_Shrub_y,(outVars_y%sp(2)%nppLeaf + outVars_y%sp(2)%nppStem)*24*365)
        ! call GetSimuData_var(vars4MCMC%ANPP_Tree_y, (outVars_y%sp(1)%nppLeaf + outVars_y%sp(1)%nppStem)*24*365)
        ! call GetSimuData_var(vars4MCMC%NPP_sphag_y,  outVars_y%sp(3)%npp*24*365)
        ! call GetSimuData_var(vars4MCMC%BNPP_y,      (outVars_y%sp(1)%nppRoot + outVars_y%sp(2)%nppRoot)*24*365)        ! tree + shrub
        ! call GetSimuData_var(vars4MCMC%er_d,        (outVars_d%sp(2)%ra + outVars_d%sp(3)%ra + outVars_d%rh)*24)          ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%er_h,        (outVars_h%sp(2)%ra + outVars_h%sp(3)%ra + outVars_h%rh))          ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%gpp_d,       & ! change to hourly
        !     (outVars_h%sp(2)%gpp + outVars_h%sp(3)%gpp))         ! Shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%nee_d,       &
        !     ((outVars_d%sp(2)%ra  + outVars_d%sp(3)%ra  + outVars_d%rh)-&
        !      (outVars_d%sp(2)%gpp + outVars_d%sp(3)%gpp))*24)         ! Shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%nee_h,       &
        !     ((outVars_h%sp(2)%ra  + outVars_h%sp(3)%ra  + outVars_h%rh) -&
        !      (outVars_h%sp(2)%gpp + outVars_h%sp(3)%gpp)))         ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%LAI_d,       &
        !     (outVars_d%sp(1)%LAI + outVars_d%sp(2)%LAI))         ! tree  + Shrub
        ! ! !
        ! call GetSimuData_var(vars4MCMC%leaf_mass_shrub_y, outVars_y%sp(2)%cleaf)
        ! call GetSimuData_var(vars4MCMC%stem_mass_shrub_y, outVars_y%sp(2)%cStem)

        ! call GetSimuData_var(vars4MCMC%cPlant_tree_y, outVars_d%sp(1)%cLeaf+outVars_d%sp(1)%cStem)
        ! call GetSimuData_var(vars4MCMC%cPlant_sphag_y, outVars_d%sp(3)%cLeaf+outVars_d%sp(3)%cStem + outVars_d%sp(3)%cRoot)
        ! call GetSimuData_var(vars4MCMC%cSoil_y, outVars_d%cSoilFast + outVars_d%cSoilSlow + outVars_d%cSoilPassive)

        ! call GetSimuData_var(vars4MCMC%leaf_resp_shrub_d, outVars_d%sp(2)%raLeaf*24) 
        ! call GetSimuData_var(vars4MCMC%leaf_resp_tree_d,  outVars_d%sp(1)%raLeaf*24) 
        ! ! methane
        ! call GetSimuData_var(vars4MCMC%ch4_d, outVars_d%wetlandCH4*24) 
        ! call GetSimuData_var(vars4MCMC%ch4_h, outVars_h%wetlandCH4) 
        ! ! 
        ! call GetSimuData_var(vars4MCMC%CN_shag_d,    ((outVars_d%sp(3)%cleaf + outVars_d%sp(3)%cStem + &
        !                                               outVars_d%sp(3)%cRoot)/(outVars_d%sp(3)%nLeaf + &
        !                                               outVars_d%sp(3)%nStem + outVars_d%sp(3)%nRoot))) 
        ! call GetSimuData_var(vars4MCMC%photo_shrub_d, outVars_d%sp(2)%gpp*24) 
        ! call GetSimuData_var(vars4MCMC%photo_tree_d,  outVars_d%sp(1)%gpp*24) 
        ! ---------------------------------weight the NPP ------------------------
        ! call GetSimuData_var(vars4MCMC%ANPP_Shrub_y,&
        !     in_st%sp(2)%pft_weight*(outVars_y%sp(2)%nppLeaf + outVars_y%sp(2)%nppStem)*24*365)
        ! call GetSimuData_var(vars4MCMC%ANPP_Tree_y, &
        !     in_st%sp(1)%pft_weight*(outVars_y%sp(1)%nppLeaf + outVars_y%sp(1)%nppStem)*24*365)
        ! call GetSimuData_var(vars4MCMC%NPP_sphag_y,  in_st%sp(3)%pft_weight*outVars_y%sp(3)%npp*24*365)
        ! call GetSimuData_var(vars4MCMC%BNPP_y,      &
        !     (in_st%sp(1)%pft_weight*outVars_y%sp(1)%nppRoot + in_st%sp(2)%pft_weight*outVars_y%sp(2)%nppRoot)*24*365)        ! tree + shrub
        ! call GetSimuData_var(vars4MCMC%er_d,        (outVars_d%sp(2)%ra + outVars_d%sp(3)%ra + outVars_d%rh)*24)          ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%er_h,        (outVars_h%sp(2)%ra + outVars_h%sp(3)%ra + outVars_h%rh))          ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%gpp_d,       & ! change to hourly
        !     (outVars_h%sp(2)%gpp + outVars_h%sp(3)%gpp))         ! Shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%nee_d,       &
        !     ((outVars_d%sp(2)%ra  + outVars_d%sp(3)%ra  + outVars_d%rh)-&
        !      (outVars_d%sp(2)%gpp + outVars_d%sp(3)%gpp))*24)         ! Shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%nee_h,       &
        !     ((outVars_h%sp(2)%ra  + outVars_h%sp(3)%ra  + outVars_h%rh) -&
        !      (outVars_h%sp(2)%gpp + outVars_h%sp(3)%gpp)))         ! shrub + sphag.
        ! call GetSimuData_var(vars4MCMC%LAI_d,       &
        !     (outVars_d%sp(1)%LAI + outVars_d%sp(2)%LAI))         ! tree  + Shrub
        ! ! !
        ! call GetSimuData_var(vars4MCMC%leaf_mass_shrub_y, in_st%sp(2)%pft_weight*outVars_y%sp(2)%cleaf)
        ! call GetSimuData_var(vars4MCMC%stem_mass_shrub_y, in_st%sp(2)%pft_weight*outVars_y%sp(2)%cStem)

        ! call GetSimuData_var(vars4MCMC%cPlant_tree_y,  in_st%sp(1)%pft_weight*(outVars_d%sp(1)%cLeaf+outVars_d%sp(1)%cStem))
        ! call GetSimuData_var(vars4MCMC%cPlant_sphag_y, &
        !     in_st%sp(3)%pft_weight*(outVars_d%sp(3)%cLeaf+outVars_d%sp(3)%cStem + outVars_d%sp(3)%cRoot))
        ! call GetSimuData_var(vars4MCMC%cSoil_y, outVars_d%cSoilFast + outVars_d%cSoilSlow + outVars_d%cSoilPassive)

        ! call GetSimuData_var(vars4MCMC%leaf_resp_shrub_d, outVars_d%sp(2)%raLeaf*24) 
        ! call GetSimuData_var(vars4MCMC%leaf_resp_tree_d,  outVars_d%sp(1)%raLeaf*24) 
        ! ! methane
        ! call GetSimuData_var(vars4MCMC%ch4_d, outVars_d%wetlandCH4*24) 
        ! call GetSimuData_var(vars4MCMC%ch4_h, outVars_h%wetlandCH4) 
        ! ! 
        ! call GetSimuData_var(vars4MCMC%CN_shag_d,    ((outVars_d%sp(3)%cleaf + outVars_d%sp(3)%cStem + &
        !                                               outVars_d%sp(3)%cRoot)/(outVars_d%sp(3)%nLeaf + &
        !                                               outVars_d%sp(3)%nStem + outVars_d%sp(3)%nRoot))) 
        ! call GetSimuData_var(vars4MCMC%photo_shrub_d, outVars_d%sp(2)%gpp*24) 
        ! call GetSimuData_var(vars4MCMC%photo_tree_d,  outVars_d%sp(1)%gpp*24) 
        ! ----------------------------------------------------------------------------
        ! call update_mcmc_tot_outputs(tot_paramsets_outs_d, outVars_d, nDay)

        ! add the weight to equal the observation data
        call GetSimuData_var(vars4MCMC%ANPP_Shrub_y,(outVars_y%sp(2)%nppLeaf + outVars_y%sp(2)%nppStem)*24*365*&
                                                     in_st%sp(2)%pft_weight)
        call GetSimuData_var(vars4MCMC%ANPP_Tree_y, (outVars_y%sp(1)%nppLeaf + outVars_y%sp(1)%nppStem)*24*365*&
                                                     in_st%sp(1)%pft_weight)
        call GetSimuData_var(vars4MCMC%NPP_sphag_y,  in_st%sp(3)%pft_weight*outVars_y%sp(3)%npp*24*365)
        call GetSimuData_var(vars4MCMC%BNPP_y,       &
            (in_st%sp(1)%pft_weight*outVars_y%sp(1)%nppRoot + in_st%sp(2)%pft_weight*outVars_y%sp(2)%nppRoot)*24*365)        ! tree + shrub
        call GetSimuData_var(vars4MCMC%er_d,        &
            (in_st%sp(2)%pft_weight*outVars_d%sp(2)%ra + in_st%sp(3)%pft_weight*outVars_d%sp(3)%ra + outVars_d%rh)*24)          ! shrub + sphag.
        call GetSimuData_var(vars4MCMC%er_h,        &
            (in_st%sp(2)%pft_weight*outVars_h%sp(2)%ra + in_st%sp(3)%pft_weight*outVars_h%sp(3)%ra + outVars_h%rh))          ! shrub + sphag.
        call GetSimuData_var(vars4MCMC%gpp_d,       & ! change to hourly
            (in_st%sp(2)%pft_weight*outVars_h%sp(2)%gpp + in_st%sp(3)%pft_weight*outVars_h%sp(3)%gpp))         ! Shrub + sphag.
        call GetSimuData_var(vars4MCMC%nee_d,       &
            ((in_st%sp(2)%pft_weight*outVars_d%sp(2)%ra  + in_st%sp(3)%pft_weight*outVars_d%sp(3)%ra  + outVars_d%rh)-&
             (in_st%sp(2)%pft_weight*outVars_d%sp(2)%gpp + in_st%sp(3)%pft_weight*outVars_d%sp(3)%gpp))*24)         ! Shrub + sphag.
        call GetSimuData_var(vars4MCMC%nee_h,       &
            ((in_st%sp(2)%pft_weight*outVars_h%sp(2)%ra  + in_st%sp(3)%pft_weight*outVars_h%sp(3)%ra  + outVars_h%rh) -&
             (in_st%sp(2)%pft_weight*outVars_h%sp(2)%gpp + in_st%sp(3)%pft_weight*outVars_h%sp(3)%gpp)))         ! shrub + sphag.
        call GetSimuData_var(vars4MCMC%LAI_d,       &
            (in_st%sp(1)%pft_weight*outVars_d%sp(1)%LAI  + in_st%sp(2)%pft_weight*outVars_d%sp(2)%LAI))         ! tree  + Shrub

        call GetSimuData_var(vars4MCMC%rh_y,       outVars_y%rh*24*365)         ! tree  + Shrub
        
        ! call GetSimuData_var(vars4MCMC%leaf_mass_shrub_y, in_st%sp(2)%pft_weight*outVars_y%sp(2)%cleaf)
        ! call GetSimuData_var(vars4MCMC%stem_mass_shrub_y, in_st%sp(2)%pft_weight*outVars_y%sp(2)%cStem)

        call GetSimuData_var(vars4MCMC%leaf_mass_shrub_y, in_st%sp(2)%pft_weight*outVars_d%sp(2)%cleaf)
        call GetSimuData_var(vars4MCMC%stem_mass_shrub_y, in_st%sp(2)%pft_weight*outVars_d%sp(2)%cStem)

        call GetSimuData_var(vars4MCMC%cPlant_tree_y, in_st%sp(1)%pft_weight*(outVars_d%sp(1)%cLeaf+outVars_d%sp(1)%cStem))
        call GetSimuData_var(vars4MCMC%cPlant_sphag_y, &
            in_st%sp(3)%pft_weight*(outVars_d%sp(3)%cLeaf+outVars_d%sp(3)%cStem + outVars_d%sp(3)%cRoot))
        call GetSimuData_var(vars4MCMC%cSoil_y, outVars_d%cSoilFast + outVars_d%cSoilSlow + outVars_d%cSoilPassive)

        call GetSimuData_var(vars4MCMC%leaf_resp_shrub_d, outVars_d%sp(2)%raLeaf*24) 
        call GetSimuData_var(vars4MCMC%leaf_resp_tree_d,  outVars_d%sp(1)%raLeaf*24) 
        ! methane
        call GetSimuData_var(vars4MCMC%ch4_d, outVars_d%wetlandCH4*24) 
        call GetSimuData_var(vars4MCMC%ch4_h, outVars_h%wetlandCH4) 
        call GetSimuData_var(vars4MCMC%ch4_y, outVars_y%wetlandCH4*24*365) 
        ! 
        call GetSimuData_var(vars4MCMC%CN_shag_d,    ((outVars_d%sp(3)%cleaf + outVars_d%sp(3)%cStem + &
                                                      outVars_d%sp(3)%cRoot)/(outVars_d%sp(3)%nLeaf + &
                                                      outVars_d%sp(3)%nStem + outVars_d%sp(3)%nRoot))) 
        call GetSimuData_var(vars4MCMC%photo_shrub_d, outVars_d%sp(2)%gpp*24) 
        call GetSimuData_var(vars4MCMC%photo_tree_d,  outVars_d%sp(1)%gpp*24) 

        call GetSimuData_var(vars4MCMC%zwt_h,  outVars_h%wtd) 

        call GetSimuData_var(vars4MCMC%lai_shrub_d,  outVars_d%sp(2)%lai) 
        call GetSimuData_var(vars4MCMC%lai_tree_d,  outVars_d%sp(1)%lai) 
        call GetSimuData_var(vars4MCMC%photo_tree_h, outVars_h%sp(1)%Aleaf(3)) 

        call GetSimuData_var(vars4MCMC%bnpp_tree_y,  in_st%sp(1)%pft_weight*outVars_y%sp(1)%nppRoot*24*365) 
        call GetSimuData_var(vars4MCMC%bnpp_shrub_y, in_st%sp(2)%pft_weight*outVars_y%sp(2)%nppRoot*24*365) 

        call GetSimuData_var(vars4MCMC%gpp_tree_y,  in_st%sp(1)%pft_weight*outVars_y%sp(1)%gpp*24*365) 

    end subroutine GetSimuData

    subroutine GetSimuData_var(var_obsData, var_mdData)
        implicit none
        type(interCostVariable), intent(inout) :: var_obsData
        real(8), intent(in) :: var_mdData

        if(var_obsData%existOrNot)then  ! if the observation file is existed
            if(var_obsData%mc_itime <= size(var_obsData%obsData, dim=1))then  ! still have observation data not being matched
                do while(var_obsData%obsData(var_obsData%mc_itime, 1) .lt. forcing(1)%year) ! some observation is beyond the range of simulation
                    var_obsData%mdData(var_obsData%mc_itime, 4) = -9999
                    var_obsData%mc_itime = var_obsData%mc_itime + 1
                enddo

                if(var_obsData%obsData(var_obsData%mc_itime, 1) .eq. mc_iyear .and. &
                   var_obsData%obsData(var_obsData%mc_itime, 2) .eq. mc_iday  .and. &
                   var_obsData%obsData(var_obsData%mc_itime, 3) .eq. mc_ihour) then
                        var_obsData%mdData(var_obsData%mc_itime, 1) = mc_iyear
                        var_obsData%mdData(var_obsData%mc_itime, 2) = mc_iday
                        var_obsData%mdData(var_obsData%mc_itime, 3) = mc_ihour
                        var_obsData%mdData(var_obsData%mc_itime, 4) = var_mdData
                        var_obsData%mc_itime = var_obsData%mc_itime + 1
                endif
            endif
        endif
    end subroutine GetSimuData_var

    subroutine update_mcmc_tot_outputs(tot_mcmc_outputs, simu_outputs, inTime)
        implicit none
        type(mcmc_outVars_type), intent(inout) :: tot_mcmc_outputs 
        type(outvars_data_type), intent(in)    :: simu_outputs
        integer, intent(in) :: inTime
        integer :: ipft

        do ipft = 1, npft
            tot_mcmc_outputs%sp(ipft)%gpp(iDAsimu, inTime)   = simu_outputs%sp(ipft)%gpp
            tot_mcmc_outputs%sp(ipft)%nee(iDAsimu, inTime)   = simu_outputs%sp(ipft)%nee
            tot_mcmc_outputs%sp(ipft)%npp(iDAsimu, inTime)   = simu_outputs%sp(ipft)%npp
            tot_mcmc_outputs%sp(ipft)%nppLeaf(iDAsimu, inTime)   = simu_outputs%sp(ipft)%nppLeaf
            tot_mcmc_outputs%sp(ipft)%nppStem(iDAsimu, inTime)   = simu_outputs%sp(ipft)%nppStem
            tot_mcmc_outputs%sp(ipft)%nppStem(iDAsimu, inTime)   = simu_outputs%sp(ipft)%nppStem
            tot_mcmc_outputs%sp(ipft)%nppRoot(iDAsimu, inTime)   = simu_outputs%sp(ipft)%nppRoot
            tot_mcmc_outputs%sp(ipft)%nppOther(iDAsimu, inTime)   = simu_outputs%sp(ipft)%nppOther    ! According to SPRUCE-MIP, stem means above ground Stemy tissues which is different from Stem tissues.
            tot_mcmc_outputs%sp(ipft)%ra(iDAsimu, inTime)   = simu_outputs%sp(ipft)%ra
            tot_mcmc_outputs%sp(ipft)%raLeaf(iDAsimu, inTime)   = simu_outputs%sp(ipft)%raLeaf
            tot_mcmc_outputs%sp(ipft)%raStem(iDAsimu, inTime)   = simu_outputs%sp(ipft)%raStem
            tot_mcmc_outputs%sp(ipft)%raRoot(iDAsimu, inTime)   = simu_outputs%sp(ipft)%raRoot
            tot_mcmc_outputs%sp(ipft)%raOther(iDAsimu, inTime)   = simu_outputs%sp(ipft)%raOther
            tot_mcmc_outputs%sp(ipft)%rMaint(iDAsimu, inTime)   = simu_outputs%sp(ipft)%rMaint
            tot_mcmc_outputs%sp(ipft)%rGrowth(iDAsimu, inTime)   = simu_outputs%sp(ipft)%rGrowth
            tot_mcmc_outputs%sp(ipft)%nbp(iDAsimu, inTime)   = simu_outputs%sp(ipft)%nbp
            ! Carbon Pools  (KgC m-2)
            tot_mcmc_outputs%sp(ipft)%cLeaf(iDAsimu, inTime)   = simu_outputs%sp(ipft)%cLeaf
            tot_mcmc_outputs%sp(ipft)%cStem(iDAsimu, inTime)   = simu_outputs%sp(ipft)%cStem
            tot_mcmc_outputs%sp(ipft)%cRoot(iDAsimu, inTime)   = simu_outputs%sp(ipft)%cRoot
            ! Nitrogen pools (kgN m-2)
            tot_mcmc_outputs%sp(ipft)%nLeaf(iDAsimu, inTime)   = simu_outputs%sp(ipft)%nLeaf
            tot_mcmc_outputs%sp(ipft)%nStem(iDAsimu, inTime)   = simu_outputs%sp(ipft)%nStem
            tot_mcmc_outputs%sp(ipft)%nRoot(iDAsimu, inTime)   = simu_outputs%sp(ipft)%nRoot
            ! tot_mcmc_outputs%sp(ipft)%nOther(:)
            ! water fluxes (kg m-2 s-1)
            tot_mcmc_outputs%sp(ipft)%tran(iDAsimu, inTime)   = simu_outputs%sp(ipft)%tran
            ! other
            tot_mcmc_outputs%sp(ipft)%lai(iDAsimu, inTime)   = simu_outputs%sp(ipft)%lai
        enddo

        tot_mcmc_outputs%gpp(iDAsimu, inTime)            = simu_outputs%gpp    
        tot_mcmc_outputs%nee(iDAsimu, inTime)            = simu_outputs%nee
        tot_mcmc_outputs%npp(iDAsimu, inTime)            = simu_outputs%npp
        tot_mcmc_outputs%nppLeaf(iDAsimu, inTime)        = simu_outputs%nppLeaf
        tot_mcmc_outputs%nppStem(iDAsimu, inTime)        = simu_outputs%nppStem
        tot_mcmc_outputs%nppStem(iDAsimu, inTime)        = simu_outputs%nppStem
        tot_mcmc_outputs%nppRoot(iDAsimu, inTime)        = simu_outputs%nppRoot
        tot_mcmc_outputs%nppOther(iDAsimu, inTime)       = simu_outputs%nppOther
        tot_mcmc_outputs%ra(iDAsimu, inTime)             = simu_outputs%ra
        tot_mcmc_outputs%raLeaf(iDAsimu, inTime)         = simu_outputs%raLeaf
        tot_mcmc_outputs%raStem(iDAsimu, inTime)         = simu_outputs%raStem
        tot_mcmc_outputs%raRoot(iDAsimu, inTime)         = simu_outputs%raRoot
        tot_mcmc_outputs%raOther(iDAsimu, inTime)        = simu_outputs%raOther
        tot_mcmc_outputs%rMaint(iDAsimu, inTime)         = simu_outputs%rMaint
        tot_mcmc_outputs%rGrowth(iDAsimu, inTime)        = simu_outputs%rGrowth
        tot_mcmc_outputs%rh(iDAsimu, inTime)             = simu_outputs%rh
        tot_mcmc_outputs%nbp(iDAsimu, inTime)            = simu_outputs%nbp
        tot_mcmc_outputs%wetlandCH4(iDAsimu, inTime)     = simu_outputs%wetlandCH4
        tot_mcmc_outputs%wetlandCH4prod(iDAsimu, inTime) = simu_outputs%wetlandCH4prod
        tot_mcmc_outputs%wetlandCH4cons(iDAsimu, inTime) = simu_outputs%wetlandCH4cons
        ! Carbon Pools  (KgC m-2)
        tot_mcmc_outputs%cLeaf(iDAsimu, inTime)          = simu_outputs%cLeaf
        tot_mcmc_outputs%cStem(iDAsimu, inTime)          = simu_outputs%cStem
        tot_mcmc_outputs%cRoot(iDAsimu, inTime)          = simu_outputs%cRoot
        tot_mcmc_outputs%cOther(iDAsimu, inTime)         = simu_outputs%cOther
        tot_mcmc_outputs%cLitter(iDAsimu, inTime)        = simu_outputs%cLitter
        tot_mcmc_outputs%cLitterCwd(iDAsimu, inTime)     = simu_outputs%cLitterCwd
        tot_mcmc_outputs%cSoil(iDAsimu, inTime)          = simu_outputs%cSoil
        tot_mcmc_outputs%cSoilLevels(iDAsimu, inTime, :) = simu_outputs%cSoilLevels
        tot_mcmc_outputs%cSoilFast(iDAsimu, inTime)      = simu_outputs%cSoilFast
        tot_mcmc_outputs%cSoilSlow(iDAsimu, inTime)      = simu_outputs%cSoilSlow
        tot_mcmc_outputs%cSoilPassive(iDAsimu, inTime)   = simu_outputs%cSoilPassive
        tot_mcmc_outputs%CH4(iDAsimu, inTime, :)         = simu_outputs%CH4
        ! Nitrogen fluxes (kgN m-2 s-1)
        tot_mcmc_outputs%fBNF(iDAsimu, inTime)           = simu_outputs%fBNF
        tot_mcmc_outputs%fN2O(iDAsimu, inTime)           = simu_outputs%fN2O
        tot_mcmc_outputs%fNloss(iDAsimu, inTime)         = simu_outputs%fNloss
        tot_mcmc_outputs%fNnetmin(iDAsimu, inTime)       = simu_outputs%fNnetmin
        tot_mcmc_outputs%fNdep(iDAsimu, inTime)          = simu_outputs%fNdep
        ! Nitrogen pools (kgN m-2)
        tot_mcmc_outputs%nLeaf(iDAsimu, inTime)          = simu_outputs%nLeaf
        tot_mcmc_outputs%nStem(iDAsimu, inTime)          = simu_outputs%nStem
        tot_mcmc_outputs%nRoot(iDAsimu, inTime)          = simu_outputs%nRoot
        tot_mcmc_outputs%nOther(iDAsimu, inTime)         = simu_outputs%nOther
        tot_mcmc_outputs%nLitter(iDAsimu, inTime)        = simu_outputs%nLitter
        tot_mcmc_outputs%nLitterCwd(iDAsimu, inTime)     = simu_outputs%nLitterCwd
        tot_mcmc_outputs%nSoil(iDAsimu, inTime)          = simu_outputs%nSoil
        tot_mcmc_outputs%nMineral(iDAsimu, inTime)       = simu_outputs%nMineral
        ! energy fluxes (W m-2)
        tot_mcmc_outputs%hfls(iDAsimu, inTime)           = simu_outputs%hfls
        tot_mcmc_outputs%hfss(iDAsimu, inTime)           = simu_outputs%hfss
        tot_mcmc_outputs%SWnet(iDAsimu, inTime)          = simu_outputs%SWnet
        tot_mcmc_outputs%LWnet(iDAsimu, inTime)          = simu_outputs%LWnet
        ! water fluxes (kg m-2 s-1)
        tot_mcmc_outputs%ec(iDAsimu, inTime)             = simu_outputs%ec
        tot_mcmc_outputs%tran(iDAsimu, inTime)           = simu_outputs%tran
        tot_mcmc_outputs%es(iDAsimu, inTime)             = simu_outputs%es
        tot_mcmc_outputs%hfsbl(iDAsimu, inTime)          = simu_outputs%hfsbl
        tot_mcmc_outputs%mrro(iDAsimu, inTime)           = simu_outputs%mrro
        tot_mcmc_outputs%mrros(iDAsimu, inTime)          = simu_outputs%mrros
        tot_mcmc_outputs%mrrob(iDAsimu, inTime)          = simu_outputs%mrrob
        ! other
        tot_mcmc_outputs%mrso(iDAsimu, inTime, :)        = simu_outputs%mrso
        tot_mcmc_outputs%tsl(iDAsimu, inTime, :)         = simu_outputs%tsl
        tot_mcmc_outputs%tsland(iDAsimu, inTime)         = simu_outputs%tsland
        tot_mcmc_outputs%wtd(iDAsimu, inTime)            = simu_outputs%wtd
        tot_mcmc_outputs%snd(iDAsimu, inTime)            = simu_outputs%snd
        tot_mcmc_outputs%lai(iDAsimu, inTime)            = simu_outputs%lai
        return
    end subroutine update_mcmc_tot_outputs


    ! subroutine ReadLineNumFromFile(filepath, count_lines)
    !     implicit none
    !     character(len=*), intent(in) :: filepath
    !     character(len=100) header, line
    !     integer STAT, count_lines

    !     open(38, file=trim(filepath), status="old", action="read", iostat=STAT) ! open file
    !     read(38, '(a100)') header           ! read the header of the file
    !     count_lines = 0                     ! initilize the count_lines
    !     do while(.TRUE.)
    !         read(38, *, iostat=STAT) line   ! read each line
    !         if(STAT .ne. 0) exit            ! until the end of the file
    !         count_lines = count_lines + 1   ! recording the count of the lines
    !     enddo
    !     return
    ! end subroutine ReadLineNumFromFile

    subroutine ReadObsDataFromFile(filepath, count_lines, resData)
        ! Jian: note that this subroutine is used to read the observational data. 
        ! The observational file must be .txt format, and with 5 columns: year, doy, hour, value, std.
        implicit none
        character(len=*), intent(in) :: filepath
        character(len=100) header
        integer STAT, count_lines, iline, n
        real(8) resData(count_lines, 5), readData(5) ! 5 colunms: year, doy, hour, value, std.

        OPEN(34, FILE=trim(filepath), status='old', ACTION='read', IOSTAT=STAT) ! open file
        read(34, '(a100)') header
        iline = 1
        do
            read(34,*,iostat=STAT, end=567) (readData(n), n = 1, 5)
            if(STAT .ne. 0) exit
            resData(iline, :) = readData
            iline = iline + 1
        end do
567     continue
        close(34)
        return
    end subroutine ReadObsDataFromFile

end module mcmc_mod