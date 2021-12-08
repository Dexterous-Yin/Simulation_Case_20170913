# Simulation_Case_20170913
This set of code is used to reproduce the observations for pitch-angle distribution of protons fluxes from Van Allen Probe B during the magnetic dip event on 13 September 2017. The basic algorithm of this set of code is to use the Runge-Kutta method (ode45 fuction in MATLAB) to calculate the motion of particle's guiding centers (Northrop, 1963).

Data Preparation
-----------
To run this set of code, first one should put the ephemeris datasets and HOPE measurement files from Van Allen Probe B during this event into the "/data" folder. Specially, the ephemeris datasets are "rbspb_def_MagEphem_TS04D_20170912_v1.0.0.h5" and "rbspb_def_MagEphem_TS04D_20170913_v1.0.0.h5". The HOPE-Level 3 datasets are "rbspb_rel04_ect-hope-pa-l3_20170912_v7.4.0.cdf" and "rbspb_rel04_ect-hope-pa-l3_20170913_v7.2.0.cdf".

Code Description
------
Trace_magdip_proton_rbsp_case.m: main function to simulate the proton fluxes at virtual spacecraft location (which follows the equatorial projection of the Van Allen Probe B orbit).
