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
**`propagations.csv`**: This file includes the spin spiral propagation path that the user selected. Column arrangement: [qx,qy,qz]

## Interactive Command Line (ICL)
After *executor.sh* runs, the ICL is ready to accept commands to set up input files. The ICL:
```
[hpc a1_run]$ ./executor.sh
==========Now you will provide information that generates INCAR files for Spin-Spiral Calculations=======
Enter the Number of magnetic basis interactions in the unit cell: 3
Enter position index of magnetic atom in unit-cell for each basis interaction : [[0,2],[0,2]]
Enter file name including propagation vectors: propagations.csv
Enter lattice constant of unit-cell: 5.99
Enter name of Bravais lattice: fcc
Enter atomic types: [1,1,2,3]
Enter unfolding mesh of IBZ points: [10,10,10]
Enter Scaling factor of energy, meV or mRy: meV
Enter file name of output INCAR (e.g incar_ss): incar_ss
==========Enter these two attributes if there is more than one magnetic basis in the system. Otherwise, you might give 0=======
Enter SS energy as a function of 90 degrees magnetic moment phase between magnetic basis: -24.41208464
Enter positional differences of two magnetic basis: (0.25,0.25,0.25)
```
To show the user what the enter command lines do, we will use the Ni2MnGa Heusler alloy structure as the demonstration subject. 

The primitive unit cell of Ni2MnGa has two magnetic basis atoms, Ni and Mn. Therefore, there are 3 magnetic interaction bases: Ni-Ni, Mn-Mn, and Mn-Ni:
```
Enter the number of magnetic basis interactions in the unit cell: 3
```
The magnetic atom positions (Ni, Mn) in the Ni2MnGa POSCAR are 0 (Ni) and 2 (Mn). One Ni is sufficient to represent the magnetic basis of Ni. The first magnetic interaction is between Ni-Ni, where the interaction is represented by index 0. The second interaction, Mn-Mn, is represented by index 2. Therefore, intra-lattice Mn-Mn and Ni-Ni interactions are represented by indices 0 and 2, forming the sublist [0,2]. Inter-lattice interaction (Mn-Ni) is represented by the sublist [0,2]. Together, the sublists merge into a list to represent magnetic basis interactions, where there are two magnetic sublattices in the system:
```
Enter position index of magnetic atom in unit-cell for each basis interaction : [[0,2],[0,2]]
```
The *unfolding mesh* list should match all propagation points before being folded along the irreducible Brillouin zone. For example, 1000 propagation points fold into 47 special propagation points due to fcc symmetry. The ICL asks the user how many unfoldings will apply along the x, y and z directions. In the case of 1000 propagation points of the fcc Bravais lattice system:
```
Enter unfolding mesh of IBZ points: [10,10,10]
```

