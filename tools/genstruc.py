import os, sys
import numpy as np
'''
Author: Qunchao Tong <tqc@calypso.cn>
Date: 2019.10.29
The aim of this script is to generate structures with variable components and variable density

'''

def read_input():
    f = open('input.dat','r')
    xfile = f.readlines()
    f.close()
    yfile = xfile.copy()
    lcomps = False
    lcaly = False
    for line in yfile:
        if  line[0] == '#':
            continue

        if 'NumberOfAtoms' in line:
            xfile.remove(line)

        if 'NumberOfFormula' in line:
            xfile.remove(line)

        if 'NameOfAtoms' in line:
            name = line
            xfile.remove(line)

        if 'Split' in line:
            xfile.remove(line)

        if 'Volume' in line:
            xfile.remove(line)
            vol = line.split('=')[1].split()
            if len(vol) == 1:
                print('\nPlease fill the range of volume per atom\n')
                sys.exit(0)

        if 'PopSize' in line:
            xfile.remove(line)
            nstruct = int(line.split('=')[1])

        if 'VSC' in line:
            xfile.remove(line)

        if 'NumberOfSpecies' in line:
            ns = int(line.split('=')[1])
            
        if 'NumberOfComponents' in line:
            ncomps = int(line.split('=')[1])
            lcomps = True
            xfile.remove(line)

        if 'CalypsoPath' in line:
            CALYPSO_PATH = line.split('=')[1]
            lcaly = True
            xfile.remove(line)

    if not lcaly:
        print('\nPlease fill CalypsoPath in input.dat\n')
        sys.exit(0)

    if not lcomps:
        print('\nPlease fill NumberOfComponents in input.dat\n')
        sys.exit(0)

    try:
        ctrl_index = xfile.index('@CtrlRange\n')
    except:
        print('\nPlease fill @CtrlRange in input.dat\n')
        sys.exit(0)

    CM = []
    for i in range(ns):
        CM.append(xfile[ctrl_index + 1 + i].split())
    return xfile, CM, vol, nstruct, ncomps, name, CALYPSO_PATH

def write_input(xfile, popsize, numberofatoms, volume, name):
    f = open('input2.dat','w')
    f.write('%9s %d\n' % ('PopSize =', popsize))
    f.write(name)
    f.write('%15s' % 'NumberOfAtoms =')
    for iu in numberofatoms:
        f.write(' %d' % iu)
    f.write('\n')
    f.write('%8s %8.3f\n' % ('Volume =', volume))
    f.write('NumberOfFormula = 1  1\n')
    f.write('Split = T\n')
    for iu in xfile:
        f.write(iu)
    f.close()

xfile, CM, vol, nstruct, ncomps, name, CALY_PATH= read_input()          
popsize = int(nstruct/ncomps)
CM = np.array(CM, int)
vol = np.array(vol, float)
nspecies, _ = CM.shape

f = open('genstruc.log','w')
os.system('mkdir structure')
index = 0
for icomp in range(ncomps):
    numberofatoms = []
    for i in range(nspecies):
        numberofatoms.append(np.random.randint(low=CM[i,0], high=CM[i,1]))
    natoms = sum(numberofatoms)
    voll = vol[0] * natoms
    volu = vol[1] * natoms
    x = np.random.random()
    volume = voll + (volu - voll) * x
    print('**************')
    print('Paras', numberofatoms, volume)
    print('**************')
    f.write('%d [' % icomp)
    for ina in numberofatoms:
        f.write('%d ' % ina)
    f.write('] %8.3f %8.3f\n' % (volume, volume/natoms))
    write_input(xfile, popsize, numberofatoms, volume, name)
    os.system('mkdir %s; cp input2.dat %s/input.dat; cd %s; %s' % (str(icomp), str(icomp), str(icomp), CALY_PATH))
    for istruc in range(popsize):
        index += 1
        os.system('cp %s/POSCAR_%s structure/POSCAR_%s' % (str(icomp), str(istruc + 1), str(index)))
    os.system('rm -rf  ' +str(icomp))
os.system('rm input2.dat')
f.close()



