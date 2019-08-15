#!/usr/bin/env python
#__author__ == "TQC"
import sys
import numpy as np
import os
import glob
import xml.etree.ElementTree as ET
import random
import re
#from get_bond_min import get_bond_min

a2b = 1.0/0.529177
xx = ['X','Y','Z']
'''
extract structure from OUTCAR
Author:	tongqunchao
2017/05/23
use re module to get lattice
'''
def extract(finput):
	f = open(finput,'r')
	nstru = -1
	atoms = []
	lat = []
	#litem = []
	pos = []
	force = []
	energy = []
	stress = []
        ethalpy = []
	stress_temp = []
	stress_sort = []
	typp = []
        element = []
        elements = []
	temp = []
	templat = []
        na = 0
        ns = 0
        nv = -1
        volume = []

	while True:
		line = f.readline()
		if len(line) == 0:
			break
                if 'VRHFIN' in line:
                        eletemp = line.split()[1].strip()
                        ele = filter(str.isalpha, eletemp)
                        element.append(ele)
		if 'ions per type' in line:
			atoms = line.split('=')[1].split()
			atoms = map(int,atoms)
                        #print element
			#print  atoms
                        for i in range(len(element)):
                            for j in range(atoms[i]):
                                elements.append(element[i])
                        #print elements
        #                sys.exit(0)
			liner(atoms,typp)
			na = np.sum(atoms)
			ns = len(atoms)
	#		print ns,na
		if 'volume of cell :' in line:
			nv += 1
			if  nv == 0:
				continue
			volume.append(float(line.split()[4].strip())/na)
		if 'direct lattice vectors' in line:
			nstru = nstru + 1
			if nstru==0:
				continue
			for i in range(3):
				temp = f.readline()
				templat = re.sub(r'-',r' -',temp).split()
				lat.append(map(float,templat[0:3]))
		if 'POSITION' in line:
			f.readline()
			for i in range(na):
				litem = f.readline().split()                     
				pos.append(map(float,litem[0:3]))
				force.append(map(float,litem[3:]))
		#if 'Total' in line:
		if 'in kB' in line:
			try:
				stress_temp = map(float,line.split()[2:8])
			except:
				continue
			stress_sort = [0.0,0.0,0.0,0.0,0.0,0.0]
			stress_sort[0] = stress_temp[0]
			stress_sort[1] = stress_temp[3]
			stress_sort[2] = stress_temp[5]
			stress_sort[3] = stress_temp[1]
			stress_sort[4] = stress_temp[4]
			stress_sort[5] = stress_temp[2]
			stress.append(stress_sort)
		if 'energy  without entropy' in line:
		#if 'enthalpy is  TOTEN    =' in line:
			ene = float(line.split()[3].strip())
			energy.append(ene)
		if 'enthalpy is  TOTEN    =' in line:
                        eth = float(line.split()[4].strip())
                        ethalpy.append(eth)
	return (volume, energy, na, ns, lat, pos, force, stress, nstru, typp, elements, ethalpy)

def liner(a,b):
    for i in range(len(a)):
        for j in range(a[i]):
            b.append(i+1)
    
def writewyc(volume, energy, nat, ntype, lat, pos, force, stress, nstruct, typp, elements, xml, enthalpy):
    global nns
    f = open('tempfile', 'a')
    if (len(energy) > 10050):	
        ll = int(len(energy))
    else:
        ll = len(energy)

    try:
        ene_last = energy[0]
    except:
        return
    for i in range(0,nstruct):
        if i >= 1 and np.abs(energy[i] - ene_last)/nat < 0.05 :
            #ene_last = energy[i]
            print i,'canceled',np.abs(energy[i] - ene_last)/nat
            continue

#        bond_dis = get_bond_min(lat[3*i:3*i+3],pos[i*nat:i*nat + nat])
#        if bond_dis < 0.5:
#            print bond_dis,'Bond_dis'
#            continue

        atomic_force = force[i*nat:(i+1)*nat]        	
	atomic_force = np.array(atomic_force)
        if ( atomic_force.max() > 80.00 or atomic_force.min() < -80.00):
            print 'force '
            continue

        temp_force = atomic_force
        force_error = np.sqrt(np.sum(temp_force * temp_force)/(3*nat))
        try:
            internal_stress = np.sum(stress[i][0] + stress[i][3] + stress[i][5])/30
            if internal_stress > 1000:
                print 'ERROR ==============>',force_error, internal_stress
                continue
        except:
            print 'FILE FOEMAT ERROR***********************>'
            continue

# >>>>>>>>>> Write structures
        f.write(str(nat) + '  ' + str(ntype) + '  ' + str(nns) + '  ' + xml + '   ' + str(volume[i]/nat) + '\n')
        nns += 1
        for j in range(3):
            f.write('%15.9f %15.9f %15.9f\n' % tuple(lat[3*i+j][0:3]))

        f.write('%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n' % tuple(stress[i]))

        for j in range(nat):
            f.write(elements[j] + '  ')
            #f.write('1'+ '  ')
            f.write(' %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n' % ( tuple(pos[i*nat+j])+tuple(force[i*nat+j])))
        try:
            f.write('%12.6f %s %12.6f %12.6f\n' % (float(energy[i]),"ENE",float(energy[i])/nat, float(enthalpy[i])/nat))
        except:
            f.write('%12.6f %s %12.6f \n' % (float(energy[i]),"ENE",float(energy[i])/nat))

    f.close()

	
def readxml():
    ns = 0
    xmlfile = glob.glob('OUTCAR_*')
    xmls = xmlfile
    if len(xmls) == 0:
        print "OUTCAR or vasprun.xml not exits"
        sys.exit(0)
    f = open('tempfile', 'w')
    f1 = open('config', 'w')
    for xml in xmls:
        print xml
        (volume, energy, nat, ntype, lat, pos, force, stress, nstruct, typp, elements, enthalpy) = extract(xml)
        writewyc(volume, energy, nat, ntype, lat, pos, force, stress, nstruct, typp, elements, xml, enthalpy)
        ns += nstruct
    f1.write(str(ns) + '\n')
    f1.close()
    os.system('cat tempfile >> config')
    os.system('rm tempfile')

if __name__ == '__main__':
    nns = 1
    readxml()

