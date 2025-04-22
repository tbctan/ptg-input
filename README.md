This repository contains data associated with the publication and is organized into the following structure. All input files are for Quantum Espresso unless otherwise noted.

**`input_files/`**
- **`opt`**
  (Contains sample pw.x input files for the following systems)
  - `ptg-agnr.pwi`: Structural optimization of A7-3 with adsorbed OH using the rev-vdW-DF2 functional.
  - `ptg-agnr-RPBE.pwi`: Self-consistent calculation using the RPBE functional using positions from the above calculation.
  - `ptg-zgnr.pwi`: Structural optimization of Z7-3 with adsorbed OH using the rev-vdW-DF2 functional.
  - `ptg-zgnr-RPBE.pwi`: Self-consistent calculation using the RPBE functional using positions from the above calculation.
- **`dos`**
(Files correspond to the A7-3 + OH system. These use the RPBE functional.)
  - `scf.pwi`: Sample input file for scf calculation
  - `nscf.pwi`: Sample input file for nscf calculation
  - `dos.pwi`: Sample input file for dos.x calculation
  - `projwfc.pwi`: Sample input file for projwfc.x calculation
  - `lobsterin`: LOBSTER input file to calculate interactions between Pt and O
- **`thermo`**
  - `thermo.py`: Python script executed to obtain the entropy and zero point energy using the thermochemistry module of the Atomic Simulation Environment software package,implemented in QE 7.2. 


## Questions
- My scf and nscf inputs for DOS calculations use PAW-PBE to make the Lobster calculations possible. Bands were also increased to 170~200 depending on the required number. (GBRV USPP does not enable COHP via lobster.) How do we address this?
- Output files to include
- All input files, or just sample template?

<!-- Kurt template
# Repository Contents

This repository contains data associated with the corresponding publication. The repository is organized into the following subdirectories:

1. **`spectroscopy`** - Includes force data and the results of the spectroscopy analysis. This directory contains outputs such as the mode/frequency locations and animations of the vibrational modes (in `.traj` and `.gif` formats). Please note that the atomic masses, which distinguish between HB and DB cases, are embedded in the `structure.traj` file.
   
2. **`input_files`** - Contains example input files for Quantum Espresso calculations, including system details and parameters. These input files were generated using the `ASE` package.

The analysis modules used for this project are part of the `VibIR-Parallel-Compute` package, which is available at the following link:

[https://morikawa-lab-osakau.gitlab.io/vibir-parallel-compute/intro.html](https://morikawa-lab-osakau.gitlab.io/vibir-parallel-compute/intro.html) -->
