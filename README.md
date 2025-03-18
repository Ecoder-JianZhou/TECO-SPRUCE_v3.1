# TECO-SPRUCE_v3.1
develop TECO-SPRUCE_v3.1

# To do
1. modify the PAR in forcing data in SPRUCE site (W/m2)
2. Tsoil and zwt use the observation
3. CH4 adds the varying POX

# Done
1. copy all code and settings from TECO-SPRUCE_v3.0
    revise the structure of TECO-SPRUCE to teco and mcmc
        TECO core code provides two interfaces: 1. parameter update; 2. outputs with date
    
    remove the nc-format outputs


Tsoil and zwt in forcing data


# add new parameter to do data assimilation
    parameter.nml: pox = 0.5
    mcmc_conf.nml: 
        st(39)%parname = 'pox'
        st(39)%parval = 0.5
        st(39)%parmin = 0.3
        st(39)%parmax = 0.9
    datatype.f90:
        site_data_type: 
            real(8) :: pox
        in_site_params:
            real(8) :: pox
        call update_site_params( ...
            in_params_vals%st_params%pox)
        subroutine update_site_params(...pox)
            real(8), intent(in) :: pox

    io_mod.f90
        subroutine read_nml_params_initValues(in_paramsNmlFile)
            real(8) :: pox 
            namelist /nml_site_params/ pox

    mcmc_mod.f90
        case ("pox");         in_params_vals%st_params%pox         = mcParams%st(i)%parval

