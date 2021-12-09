# Simulation_Case_20170913
This set of test-particle simulation codes is used to reproduce the observations for pitch-angle distribution of protons fluxes from Van Allen Probe B during the magnetic dip event on 13 September 2017. The basic algorithm of this set of codes is to use the Runge-Kutta method ("ode45" function in MATLAB) to calculate the motion of the particle's guiding center (Northrop, 1963). 

Cite the code: [![DOI](https://zenodo.org/badge/435510084.svg)](https://zenodo.org/badge/latestdoi/435510084)

Requirements
------
This code should be compatible with Windows, Mac, and Linux operating systems, with MATLAB installed.

The demo results are generated by MATLAB R2020a, in Windows 10 system.

The central arithmetic is the "ode45" function to solve the differential equations, so the requirement for the MATLAB version is very low.

To obtain the background proton flux without injection, the AP8max model is used, which could be accessed from the standard IRBEM library (https://github.com/PRBEM/IRBEM.git).

Code Description
------
**Trace_magdip_proton_rbsp_case.m**: the main function to simulate the perpendicular-moving proton fluxes at virtual spacecraft location (which follows the equatorial projection of the Van Allen Probe B orbit). The electromagnetic fields include the background fields and magnetic dip model.

**Trace_magdip_proton_rbsp_case_smpa.m**: the main function to simulate the bouncing proton fluxes at virtual spacecraft location (which follows the equatorial projection of the Van Allen Probe B orbit). The electromagnetic fields only include the background fields.

**Plot_rbsp_flux_case.m**: function to plot the comparison between Van Allen Probe B observations and simulation results. 

**Trace_magdip_proton_rbsp_case_plot.m**: main function to trace two perpendicular-moving protons observed at the dip center backward in time to the injection region, and to plot their drift trajectories. 

**Traj_3DAGCT_forward_deepen_time.m**: main computational process. This file utilizes the "ode45" function in MATLAB to calculate the motion of the particle's guiding center.

**Calc_B_deepen.m & Calc_E_deepen.m**: calculate the magnetic field and electric field during the computational process.

**odetpbar.m & textprogressbar.m**: visualize the progress of the ode solver.

Data Preparation
-----------
To run this set of code, first one should put the ephemeris datasets and HOPE measurement files from Van Allen Probe B during this event into the "/data" folder. Specially, the ephemeris datasets are "rbspb_def_MagEphem_TS04D_20170912_v1.0.0.h5" and "rbspb_def_MagEphem_TS04D_20170913_v1.0.0.h5". The HOPE-Level 3 datasets are "rbspb_rel04_ect-hope-pa-l3_20170912_v7.4.0.cdf" and "rbspb_rel04_ect-hope-pa-l3_20170913_v7.2.0.cdf".

Run Instructions
-----
After data preparation, first run **Trace_magdip_proton_rbsp_case.m**, and this step will save the results in the "/simulation_result" folder, see Line 228 `save('simulation_result\sim_90_X5.5_Lm_full.mat');`. (Run time: several hours.)

Then run **Trace_magdip_proton_rbsp_case_smpa.m**, and this step will also save the results in the "/simulation_result" folder, see Line 229 `save('simulation_result\sim_27_X5.5_Lm_full.mat');`. (Run time: several hours.)

Then run **Plot_rbsp_flux_case.m**, you will get the comparison between the energy spectra from satellite observations and simulation results. (Run time: less than one minute.)

Run **Trace_magdip_proton_rbsp_case_plot.m**, the drift trajectories of injected perpendicular-moving protons will be displayed. (Run time: less than one minute.)

Demo Output
-----
In the "/simulation_result" folder, you will find **flux_energy_spectrum.png** for comparison between the satellite observations and simulation results, and **trajectory_1.png**, **trajectory_2.png**, **trajectory_3.png** for the drift trajectories of perpendicular-moving protons. You can find the corresponding descriptions in the text for Figure 2 in our paper.

Note
------
Sometimes, there will be an error report "The text progress must be initialized with a string" for the textprogressbar.m function. Just "clear" and rerun the code.

License
------
This code is covered under the MIT License.
