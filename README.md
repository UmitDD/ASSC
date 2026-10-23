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
**`sample`**
**`Exchange_calculator.py`**: 
**`Exchange_calculator.py`**`: 
