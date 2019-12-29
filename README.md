# CALYPSO-GAP
CALYPSO-GAP is a machine learning potentials(MLPs) to calculate the properties of multi-species materials.  
There are two key methods in machine learning potentials, one is descriptors another is machine learning model.  
In the CALYPSO-GAP, the weighted atom centered symmetry functions(wACSF) is used as descriptors, 
wACSF based on the original ACSF suggested by Parrinello and Behler, 
original ACSF describe the atomic environment by a series of radial symmetry functions(RSF) and
angular symmetry functions(ASF). To overcome the scaling problem of ACSF, a weight is assign to each atoms in ACSF
to distinguish the contribution from different atoms. This
idea come from the work by [Nongnuch Artrith](https://journals.aps.org/prb/supplemental/10.1103/PhysRevB.96.014112),  
We use Gaussian process regression(GPR) to map the atomic environment vector to materials propertirs such as total energy
atomic force and cell stress, for more information about GPR, please read the paper  
[GaussianApproximation Potentials: The Accuracy of Quantum Mechanics, without the Electrons](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.136403)  
[Gaussian Approximation Potentials: A Brief Tutorial introduction](https://onlinelibrary.wiley.com/doi/full/10.1002/qua.24927)  

# Using CALYPSO-GAP
## Fitting CALYPSO-GAP potentials,
CALYPSO-GAP module had been intergred in binary calypos.x, it is available free in charge for non-commerical use by individuals and academic or research institutions. Please register in
[CALYPSO](https://www.calypso.cn/getting-calypso)

## Using CALYPSO-GAP for materials simulation,
These are two ways to using CALYPSO-GAP potentials,  
(1) Using CALYPSO-GAP for structure prediction, this is the original intention for developing CALYPSO-GAP.  
(2) Combining the python package [GAPPY](https://github.com/tongqcx/CALYPSO-GAP/tree/master/gappy) and Atomic Simulation Environment([ASE](https://wiki.fysik.dtu.dk/ase/)) to performing calculation, such as phonon calculation, molecular simulations or geometry optimization.


