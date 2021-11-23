# CLM5PDAF

This is a custom README for the CLM5PDAF coupling as presented in 

* Strebel, L., Bogena, H., Vereecken, H., and Hendricks Franssen, H.-J.: Coupling the Community Land Model version 5.0 to the parallel data assimilation framework PDAF: Description and applications, Geosci. Model Dev. Discuss. (preprint), [https://doi.org/10.5194/gmd-2021-38](https://doi.org/10.5194/gmd-2021-38), in review, 2021.


The code in this repository is a branch of the Terrestrial System Modeling Platform (TSMP or [TerrSysMP](https://www.terrsysmp.org) ) see README_TSMP.md for more details.

To run CLM5PDAF, a system where CLM5, PDAF, and TSMP can be compiled an run independently is necessary.

TSMP requires the following file structure:

- tsmp
    - bldsva
    - clm5_0
    - pdaf1_1

The main developements of this branch can be found in:

* tsmp/bldsva/intf_DA/pdaf1_1/model/enkf_clm_mod_5.F90
* tsmp/bldsva/intf_DA/pdaf1_1/model/enkf_clm_5.F90
* tsmp/bldsva/intf_DA/pdaf1_1/tsmp/clm5_0/cime_comp_mod.F90
* tsmp/bldsva/intf_DA/pdaf1_1/tsmp/clm5_0/seq_comm_mct.F90
* tsmp/bldsva/intf_DA/pdaf1_1/framework/pdaf_terrsysmp.F90

If CLM5 and PDAF are compilable and your system is configured for TSMP, then CLM5PDAF can be built with:

./build_tsmp.ksh -v 4.4.0MCTPDAF -c clm -m *YOURMACHINENAME* -O Intel

Once built the executable tsmp-pdaf can be found at tsmp/bin/*YOURMACHINENAME*_4.4.0MCTPDAF_clm .
To generate the namelists for CLM5 in ensemble mode a helper python script called "create_ensemble_namelists.py" can be found in this repository. Of course the paths to the input files for CLM5 have to be adjusted before it can be used on a different system.

Also an example configuration file "enkfpf.par" can be found in this repo.

Once the namelists are created and customized the CLM5PDAF can be executed like this:

./tsmp-pdaf -n_modeltasks 96 -screen 3 -filtertype 2 -subtype 1 -delt_obs 1 -rms_obs 0.02 -obs_filename observation_files/WU/2009/SM_CLM

The observation files used in this study can be found in the WU.tar archive file.

After a successfull execution the python helper script collect_all_ensemble_stats.py can be used to collect the history files of each ensemble member into an ensemble mean, min, max, and variance collection for each variable.

Such a collection can then for example be used to create the figures of the study. The python helper script "create_plots.py" can be customized to achieve this.
