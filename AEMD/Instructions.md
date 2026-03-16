1. Edit 'Thermo_AEMD_CuBHT_working.in' if necessary, but note that the next step uses the structure of this file as base to make some modifications.
2. 'run_AEMD_study_3.sh' is a script that automatically set a calculation for different supercell sizes, timesteps, detla temperatures, etc. Edit this file to make the opotunes adjustements.
3. After finish the calculations, calculate the media with 'Media_replicas_AEMD_3.py'.
4. The final fitting and the calculation of the thermal properties is done with the 'fit_aemd.py'

Notes: The chimesFF sometimes gets the atoms too close, decrease timesteps to avoid that.
