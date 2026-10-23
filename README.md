# ASSC
ASSC stands for Automated Spin Spiral Calculations and generates input files that describe spin-propagation vectors along a given path. It then runs Slurm jobs, collects spin-spiral energies, and uses them to calculate desired magnetic properties, such as the exchange parameter, magnon dispersion, and the system's Curie temperature. The automated workflow was designed to run within VASP.

## Installation and Setup 
Use the requirements.txt file to set up all necessary Python libraries.
```
activate <env>
pip install -r pip_requirements.txt
```
Transfer the two Python files, one bash script and five base files to the working directory. The working directory design:
```
bcc-Fe
├── Input_generator.py
├── Exchange_calculator.py
├── executor.sh
├── example_INCAR
├── example_scr
├── example_cat
├── mag_data
├── const_data
└── propagations.csv    
```
**`Input_generator.py`**: This Python script generates the INCAR input file for VASP, the slurm_job to execute, and the cat file to merge energies.
**`Exchange_calculator.py`**: This Python script calculates Heisenberg exchange parameters from merged energies for each interacting magnetic basis in the system. 
**`executor.sh`**: This bash script executes each Python script in order and runs the Slurm script in HCP. After execution, it opens an interactive command line (ICL) where the user enters the necessary information (More info in the next section).
**`example_INCAR`**: This file includes input information about the system that is prepared by the user. *Input_generator.py* script adds the magnetic moment (MAGMOM-tag), its constraints (M_CONSTR-tag), and the spin-spiral propagation vector (QSPIRAL-tag) into the INCAR generation. 
