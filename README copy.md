# TECO-SPRUCE_v3.0
The code is TECO-SPRUCE_v3.0, which is developed based on TECO-SPRUCE_v2.0 (https://github.com/ecomary/teco_spruce_v2_0). Code is used in the paper of "The fate of peatland carbon interactively determined by elevated carbon dioxide and warming".

author: Jian Zhou 
e-mail: jz964@cornell.edu

Note: this model cannot run on CodeOcean due to the absence of netcdf.mod. Compiling it within the Docker environment is challenging because the output file must be directed to the results folder, but intermediate files are generated during the compilation process. To resolve this, you will need to download the code to your local computer and run it there.

You also can find the code in github: https://github.com/Ecoder-JianZhou/TECO-SPRUCE_v3.0

Files:
  - configs: the settings of each simulation.
  - inputs: forcing data and observation in each experimental plot.
  - outputs: the outputs of each simulation.
  - run_scripts: Compilation script and execution script.
  - src: source code files of TECO-SPRUCE_v3.0
  - README.txt
  - run_pretreat_mcmc.sh
  - run_pretreat_simu.sh
  - run_treat_mcmc.sh
  - run_treat_simu_noacc.sh
  - run_treat_simulation.sh

Environment:
  - Linux system 
  - gfortran
  - netcdf

Estimated time:
  - Running 50,000 iterations of a data assimilation process (MCMC) on a personal computer (i7-12700) requires approximately 20 hours.
  - A single simulation can be completed in one minute.

How to run:

1. Simulate the treatment experiments
    To run all experiments once:
        ./run_scripts/treat_simulation/make_desktop_all.sh
        ./run_treat_simulation.sh
    To run each experiment (P04 as an example):
        ./run_scripts/treat_simulation/make_desktop_P04.sh
        ./run_scripts/treat_simulation/run_desktop_P04.sh

2. Run the simulations for the treatment experiments without accounting for plant community changes. Use the optimized parameters from Plot 06 for all plots.
    To run all experiments once:
        ./run_scripts/treat_simu_noacc/make_desktop_all.sh
        ./run_treat_simu_noacc.sh
    To run each experiment (P04 as an example):
        ./run_scripts/treat_simu_noacc/make_desktop_P04.sh
        ./run_scripts/treat_simu_noacc/run_desktop_P04.sh

3. To obtain parameters optimized using individual plot observations.
    To run all experiments once:
        ./run_scripts/treat_mcmc/make_desktop_all.sh
        ./run_treat_mcmc.sh
    To run each experiment (P04 as an example):
        ./run_scripts/treat_mcmc/make_desktop_P04.sh
        ./run_scripts/treat_mcmc/run_desktop_P04.sh

4. To derive the initial states of each plot based on individual plot observations.
    To run all experiments once:
        ./run_scripts/pretreat_mcmc/make_desktop_all.sh
        ./run_pretreat_mcmc.sh
    To run each experiment (P04 as an example):
        ./run_scripts/pretreat_mcmc/make_desktop_P04.sh
        ./run_scripts/pretreat_mcmc/run_desktop_P04.sh

5. To simulate the pre-treatment experiments
    To run all experiments once:
        ./run_scripts/pretreat_simu/make_desktop_all.sh
        ./run_pretreat_simu.sh
    To run each experiment (P04 as an example):
        ./run_scripts/pretreat_simu/make_desktop_P04.sh
        ./run_scripts/pretreat_simu/run_desktop_P04.sh

The results will be found in "outputs"