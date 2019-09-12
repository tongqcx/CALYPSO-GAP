# weighted Atom Centered Symmetry Functions(wACSF)
ACSF is a well-developmented methods for representation of atomic environment, 
it was proposed by Behler and Parrinello firstly. For using ACSF in multi-species system,
each atom is assigned a weight to distinguish their contribution in ACSF, that is called wACSF.
The aim of this package is to convert cartesian coordinate to wACSF  
which will be used in machin learning.
## Usages
* INSTALL PACKAGE
```shell
python setup.py install
```
will install wACSF package in your machine
* USING wACSF
```python
from libwacsf.wacsf import WACSF
a = WACSF(Rcut=6.0, Nfeature= 33)
Acsf = a.car2wacsf(Lattice, Position)
```
__Rcut__ define the size of region to calculate ACSF  
__Nfearure__ define the number of symmetry function  
__Lattice__ is the cell matrix shape (3,3)  
__Position__ is the atomic position in cartesian coordinate shape (natoms,3), natoms is the number of atoms in cell, the unit of __Lattice__ and
__Position__ is angstrom  
__Acsf__ is the wacsf matrix shape (natoms, nfeature)  
