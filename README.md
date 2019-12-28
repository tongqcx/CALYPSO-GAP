# CALYPSO-GAP
CALYPSO-GAP is a machine learning potentials(MLPs) to calculate the properties of multi-species materials under high pressure.  
There are two key methods in CALYPSO-GAP
* Weighted atom centered symmetry function(wACSF)
* Gaussian process regression(GPR)  
wACSF based on the original ACSF suggested by Parrinello and Behler, original ACSF describe the atomic environment by a series of radial symmetry functions(RSF) and
angular symmetry functions(ASF). To overcome the scaling problem of ACSF, a weight is assign to each atoms in ACSF to distinguish the contribution from different atoms. This
idea come from the work by [Nongnuch Artrith](https://journals.aps.org/prb/supplemental/10.1103/PhysRevB.96.014112)
    
For information of Gaussian process regression, please read the paper  
[GaussianApproximation Potentials: The Accuracy of Quantum Mechanics, without the Electrons](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.136403)  
[Gaussian Approximation Potentials: A Brief Tutorial introduction](https://onlinelibrary.wiley.com/doi/full/10.1002/qua.24927)  

## Using GAP

