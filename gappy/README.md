# GAPPY
GAPPY is a python package based on CALYPSO-GAP for calculation total energy, atomic force, cell stress.
## Structure of GAPPY
GAPPY  
|-- ASE   
|-- example    
|-- libgap   
|-- README.md    
|-- setup.py   
  

## Usages
### Installing package libgap
```shell
python setup.py build --fcompiler=intelem
python setup.py install
```
Intel fortran compiler is recommended, It has been extensively tested and work well
### Installing ASE interface
Supposing ASE was installed in dir ${ASEPATH}, usually it looks like that
```shell
${somepath}/lib/python3.7/site-packages/ase-3.18.1-py3.7.egg/ase
```
so
```shell
cp ./ASE/gap_calc.py ${ASEPATH}/calculators
cp ./ASE/logger.py ${ASEPATH}/md
```
### Using ASE with GAPPY   
There are two examples in ./example,  
(1) BC, this is a sample example to introduce how to use GAPPY
```python
form libgap.GAP import Calculator
gap = Calculator(rcut = 6.0)
lgrad = True
gap.gap_read()
ene, force, stress, _ = gap.gap_calc(species, lat, pos, lgrad)
```
rcut: (float, units: angstrom) defining the size of sphere region to calculate wACSF  
species: (python list) the chemical symbol of each atoms  
lat: (numpy array float[3,3], units: angstrom) the cell matrix   
pos: (numpy array float[natoms,3], units: angstrom) the cartesian position of each atoms   
ene: (float, units: eV) the total energy of structure  
force: (numpy array float[natoms,3], units: eV/angstrom) the atomic force components of each atoms along X, Y, Z direction  
stress:(numpy array float[6], units: GPa) the cell stress as VASP OUTCAR format (XX YY ZZ XY YZ  ZX)  
(2) ASE-GAPPY, this example introduce how to use GAPPY with ASE  
run_opt.py a python script to run fix cell optimization  
run_md.py a python script to run molecular simulation



