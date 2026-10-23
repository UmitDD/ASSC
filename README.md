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
**`Input_generator.py`**: This Python script generates the INCAR input file for VASP, the slurm_job to execute, and the cat file to merge energies.\
**`Exchange_calculator.py`**: This Python script calculates Heisenberg exchange parameters from merged energies for each interacting magnetic basis in the system. \
**`executor.sh`**: This bash script executes each Python script in order and runs the Slurm script in HCP. After execution, it opens an interactive command line (ICL) where the user enters the necessary information (More info in the next section).\
**`example_INCAR`**: This file includes input information about the system that the user prepares. *Input_generator.py* script adds the magnetic moment (MAGMOM-tag), its constraints (M_CONSTR-tag), and the spin-spiral propagation vector (QSPIRAL-tag) into the INCAR generation. \
**`example_scr`**: This file is the base file to use for generating slurm_jobs.\
**`example_cat`**: This file is the base file to use for generating energy-merging cat bash scripts. \
**`mag_data`**: This file includes the converged magnetic moment of each atom in the system. \
**`const_data`**: This file includes the constrained vectors of the converged magnetic moment of each atom in the system. \
**`propagations.csv`**: This file includes the spin spiral propagation path that the user selected.

## Interactive Command Line (ICL)
After *executor.sh* runs, the ICL is ready to accept commands to set up input files. The ICL:
```
[hpc a1_run]$ ./executor.sh
==========Now you will provide information that generates INCAR files for Spin-Spiral Calculations=======
Enter Number of magnetic basis interactions in the unit cell: 3
Enter position index of magnetic atom in unit-cell for each basis interaction : [[0,2],[0,2]]
Enter file name including propagation vectors: propagations.csv
Enter lattice constant of unit-cell: 5.99
Enter name of Bravais lattice: fcc
Enter atomic types (e.g atom1->[0], atom1,atom2->[0,1]): [1,1,2,3]
Enter unfolding mesh of IBZ points (e.g [10,10,10]): [10,10,10]
Enter Scaling factor of energy, meV or mRy: meV
Enter file name of output INCAR (e.g incar_ss): incar_ss
==========Enter these two attributes if there is more than one magnetic basis in the system. Otherwise, you might give 0=======
Enter SS energy as a function of phase 90 between magnetic basis: -24.41208464
Enter positional differences of two magnetic basis (0.25,0.25,0.25)
```
To show the user what the enter command lines, the Ni2MnGa Heusler alloy structure will be the demonstration subject. 

The primitive unit cell of Ni2MnGa has two magnetic basis atoms, Ni and Mn. 
