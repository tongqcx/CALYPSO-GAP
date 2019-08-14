# Gaussian Approximation Potentials(GAP)
Gaussian approximation potentials is an effective machine learning potentials(MLP) to calculate the properties of multi-species crystal under high pressure.  
There are two key methods in GAP
* Weighted atom centered symmetry function(wACSF)
* Gaussian process regression(GPR)  
wACSF based on the original ACSF suggested by Parrinello and Behler, original ACSF describe the atomic environment by a series of radial symmetry functions(RSF) and
angular symmetry functions(ASF). For multi-species system AxBy there have three types of interaction in RSF
A-A
A-B  
B-B  
and six types of interaction in  ASF
A-A-A  
A-A-B 
A-B-B 
B-A-A 
B-A-B 
B-B-B 
In general, for a N species system, the number of interaction is N(N+1)/2. To overcome this scaling problem, the weighted descriptors have been 
the poineered work 

For more information of GAP, please read the paper  
[GaussianApproximation Potentials: The Accuracy of Quantum Mechanics, without the Electrons](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.136403)  
[Gaussian Approximation Potentials: A Brief Tutorial introduction](https://onlinelibrary.wiley.com/doi/full/10.1002/qua.24927)  

## Using GAP
* Compling GAP
The lapack library is need for compling GAP, So be sure the lapack library or Intel MKL library have installed.  
Entering Src dir, type make  
* Input files  
  * **control** (containing the control parameter)  
  * **config** (containing training set)  
  * **test** (containing testing set)  
  * **neural.in** (the parameter file of wACSF)  
## How to cite GAP

## Others
Thie project begin at 2019.07.08 for an effect frame to construct gaussian approximation potentials (GAP)
In this schame, there are three parts, 2-body many-body and electrostatic energy.
