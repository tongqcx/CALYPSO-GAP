# GAPPY
GAPPY is a Python2/3-compatible toolkit for calculation total energy, atomic force, cell stress using Gaussian Approximation Potentialsi(GAP).
## Structure of GAPPY
GAPPY
|-- ASE   
|-- example    
|-- libgap   
|-- README.md    
|-- setup.py   
  

## Usages
* INSTALL PACKAGE libgap
```shell
python setup.py install --fcompiler=intelem
```
Intel fortran compiler is recommended, It has been extensively tested and work well
* USING libgap   
There are two methods in libgap
    * GAP.py 
This method could calculate the energy atomic force and cell stress for given structure, using as
``` python
from libgap import GAP  
gap = GAP(rcut=6.0)  
energy, force, stress, _ = gap.gap_calc(species, lat, pos, lgrad)  
```
__There must have neural.in and gap_parameters files in WorkDir__

    * WACSF.py
ACSF is a well-developmented methods for representation of atomic environment, 
it was proposed by Behler and Parrinello firstly. For using ACSF in multi-species system,
each atom is assigned a weight to distinguish their contribution in ACSF, that is called wACSF.
The aim of this method is to convert cartesian coordinate to wACSF  
which will be used in machin learning.
```python
from libgap import WACSF
wacsf = WACSF.Wacsf(nf = 33, rcut=6.0, lgrad = False)
Acsf = wacsf.car2wacsf(Lattice, Elements, Position)
```
__rcut__ define the size of region to calculate ACSF  
__nf__ define the number of symmetry function  
__Lattice__ is the cell matrix shape (3,3)  
__Position__ is the atomic position in cartesian coordinate shape (natoms,3),
natoms is the number of atoms in cell, the unit of Lattice and Position is angstrom  
__Acsf__ is the wacsf matrix shape (natoms, nf)  
__There must have neural.in and gap_parameters files in WorkDir__


## 2019.01.04
