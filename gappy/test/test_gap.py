import numpy as np
import sys
from gap_calc  import Calculator

#sparse_cov_mat_len,des_len,theta,mm,invcmm,coeff = gap_read()

f = open('struc.dat','r')
line = f.readlines()
Na = int(line[0].split()[0])
lat = []
species = []
pos = []
for i in range(1,4):
    lat.append(line[i].split())
lat = np.array(lat,float)
for i in range(5, 5 + Na):
    species.append(line[i].split()[0])
    pos.append(line[i].split()[1:4])
pos = np.array(pos,float)
#----------------------------------------
GAP = Calculator.Gap_calculator()
ene, force, stress, variance = GAP.gap_calc(species, lat, pos)
#----------------------------------------

print ene ,ene/Na ,type(force)
print force
print stress
print variance
