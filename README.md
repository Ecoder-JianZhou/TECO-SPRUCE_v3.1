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