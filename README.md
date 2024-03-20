## WARNING: This repository is deprecated. Code development is continued in https://github.com/FormingWorlds/JANUS

## AEOLUS (temperature structure generator)

Generates a temperature profile using the generalised moist pseudoadiabat and a prescribed stratosphere. Calculates radiative fluxes using SOCRATES.

### Documentation
https://proteus-code.readthedocs.io

### Contributors (abbreviations & email addresses):
* TL - Tim Lichtenberg (tim.lichtenberg@rug.nl)
* MH - Mark Hammond (mark.hammond@physics.ox.ac.uk)
* RB – Ryan Boukrouche (ryan.boukrouche@astro.su.se)
* RJG – RJ Graham (arejaygraham@uchicago.edu)
* HN - Harrison Nicholls (harrison.nicholls@physics.ox.ac.uk)
* HII - Hamish Innes (hamish.innes@physics.ox.ac.uk)

### Repository structure

* `SocRadConv.py`               - Main AEOLUS Python script
* `demo_runaway_greenhouse.py`  - Demonstrate pure-steam runaway greenhouse OLR curve
* `demo_instellation.py`        - Calculate fluxes (and temperatures) for different instellations
* `README.md`                   - This file
* `INSTALL.md`                  - Installation instructions
* `AEOLUS.env`                  - Sets environment flags to run the code
* `utils/`                      - Utility python scripts
* `modules/`                    - Utility python scripts
* `output/`                     - Output folder
* `luminosity_tracks/`          - Stellar evolution data
* `plotting_tools/`             - Plotting scripts
* `rad_trans/`                  - Subfolder for radiative transfer code/s
* `spectral_files/`             - Spectral files for SOCRATES
* `tools/`                      - Useful tools

### Installation instructions
Follow environment- and AEOLUS-related steps in https://proteus-code.readthedocs.io.

### Run instructions
Only attempt to run AEOLUS after you have followed all of the instructions in INSTALL.md    
If using a fresh shell, it is necessary to perform the following steps:     
1. `source AEOLUS.env`    
2. `conda activate aeolus`    
Then you can run the code by running: `python SocRadConv.py`      
