"""
2019.11.03
Extracting volume, energy, pressure from config
Qunchao Tong <tqc@calypso.cn>
"""
import os, sys
import numpy as np

assert os.path.exists('./config'), "{} does not exist!".format('config')
f = open('config','r')
line = []
Volume = []
Energy = []
Pressure = []
nstruc = int(f.readline().split()[0])
print(nstruc)
for istruc in range(nstruc):
    line = f.readline()
    if len(line) == 0:
        break
    na = int(line.split()[0])
    lat = []
    for i in range(3):
        lat.append(f.readline().split())
    strs = list(map(float, f.readline().split()))
    pres = (strs[0] + strs[3] + strs[5])/30
    for i in range(na):
        f.readline()
    ene = float(f.readline().split()[0])/na
    lat = np.array(lat, float)
    volume = np.abs(np.linalg.det(lat))/na
    Volume.append(volume)
    Energy.append(ene)
    Pressure.append(pres)
f.close()

f = open('V-E-P.txt','w')
f.write('%20s, %20s, %20s\n' % ('Volume (A^3/atom)', 'Energy (eV/atom)', 'Pressure (GPa)'))
for i in range(len(Volume)):
    f.write('%20.5f, %20.5f, %20.5f\n' % (Volume[i], Energy[i], Pressure[i]))

