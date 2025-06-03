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
  - `lobsterin`: Sample LOBSTER input file to calculate interactions between Pt and O
- **`thermo`**
  - `thermo.py`: Python script executed to obtain the entropy and zero point energy using the thermochemistry module of the Atomic Simulation Environment software package,implemented in QE 7.2. 

**`structure_files/`**
- `adsorbate_start_A5-2`: Contains starting adsorbate configurations as modeled on A5-2.
- **`AGNR`** and **`ZGNR`**: Contains output structures of vdW relaxed structures in `.xsf` format, including final adsorbate structures calculated. 
