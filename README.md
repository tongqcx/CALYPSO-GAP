# CALYPSO-GAP
CALYPSO-GAP is a package of machine learning potential (MLP) to calculate the properties of multi-species materials.  
There are two key methods in machine learning potentials, one is descriptors and the other is machine learning model.  
In the CALYPSO-GAP, the weighted atom-centered symmetry functions (wACSFs) is used as descriptors. 
ACSF were proposed firstly by Parrinello and Behler to describe the atomic environment by a series of
radial symmetry functions and angular symmetry functions. To overcome the scaling problem of ACSF, the weight of each atom in ACSF is used 
to distinguish the contribution from different atoms. This
idea is inspired by the work of [Nongnuch Artrith](https://journals.aps.org/prb/supplemental/10.1103/PhysRevB.96.014112),  
We use Gaussian process regression(GPR) to map the atomic feature vector to atomic energy, total energy usually
is calculated by sum of atomic energy, and atomic force and cell stress is the derivative of total energy,
for more information about GPR, please read the paper  
[GaussianApproximation Potentials: The Accuracy of Quantum Mechanics, without the Electrons](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.136403)  
[Gaussian Approximation Potentials: A Brief Tutorial introduction](https://onlinelibrary.wiley.com/doi/full/10.1002/qua.24927)  

# Using CALYPSO-GAP
## Fitting CALYPSO-GAP potentials,
CALYPSO-GAP module is integrated in binary calypos.x, it is available free in charge for non-commerical use by individuals and academic or research institutions. Please register in
[CALYPSO](http://www.calypso.cn/)

## Using CALYPSO-GAP for materials simulation,
These are two ways to use CALYPSO-GAP potentials,  
(1) Using CALYPSO-GAP for structure prediction, which is the original intention of developing CALYPSO-GAP.  
(2) Combining the python package [GAPPY](https://github.com/tongqcx/CALYPSO-GAP/tree/master/gappy) and Atomic Simulation Environment([ASE](https://wiki.fysik.dtu.dk/ase/)) to perform calculation, such as phonon calculation, molecular simulations and geometry optimization.


