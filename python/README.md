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
The example of using libgap to calculate the properties of materials could be found
in ./example dir
* About libgap.so   
When you installed libgap package in your machine, you could use libgap.so's methods 
individually. There is some methods in __libgap.so__    
__fcar2wacsf__ could convert cartesian coordinate to weight atom centered symmetry function(wACSF)
wACSF could be used in Machine Learning as descriptors  

```python   
xx,dxdy,strs = fcar2wacsf(nf,lat,elements,pos,rcut,lgrad,na=len(elements))   
```

```python
write_array_2dim(a,name,n=shape(a,0),m=shape(a,1))   
```

fgap_read and fgap_calc will be used at the same time, GAP.py in dir libgap
give the interface to use this two methods   
```python
nsparsex,des_len,theta,mm,invcmm,coeff = fgap_read()  
```
```python
ene,force,stress,variance = fgap_calc(species,lat,pos,theta,mm,qmm,coeff,rcut,lgrad,na=len(species),nsparsex=shape(mm,0),des_len=len(theta))   
```


## 2019.01.04
