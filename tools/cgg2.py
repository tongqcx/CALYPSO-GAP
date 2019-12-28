#!/usr/bin/python
__author__ = 'tqc'
'''
try:
except(ErrorType)
2019.10.10
2019.11.03 add module of structure generation
2019.11.07 checking the cellparas
'''

import os
import time
import sys
import numpy as np
from numpy import pi, sin, cos, arccos, sqrt, dot
from numpy.linalg import norm
import argparse

class Autocalypso(object):
    def __init__(self, lrelax = False, lwrite = False, lgen = False):
       
        self.lrelax = lrelax
        self.lwrite = lwrite
        self.lgen = lgen
        calypath, PN, paraspath, optpath, NSW, pstress, ftol, gtol = self.check_input()
        self.CalyPsoPath = calypath 
        self.GulpPath = optpath + '> log &'
        self.submit = self.GulpPath
        self.MaxTime = 1800
        self.PopSize = 100
        self.StrNum = self.PopSize
        self.MaxStep = 20
        self.PickUp = False
        self.NameOfAtoms = []
        self.NumberOfAtoms = []
        self.PN = PN
        self.paraspath = paraspath
        self.NSW = NSW
        self.pstress = pstress
        self.ftol = ftol
        self.gtol = gtol
        self.f = open('split_calypso.log','w')
        self.input = None

    def readinput(self):
        f = open('input.dat','r')
        self.input = f.readlines()

        for line in self.input:
            if 'PopSize' in line:
                self.StrNum = int(line.split('=')[1])
            if 'MaxStep' in line:
                self.GenNum = int(line.split('=')[1])
            if 'MaxTime' in line:
                self.MaxTime = int(line.split('=')[1])
            if 'NameOfAtoms' in line:
                self.NameOfAtoms = line.split('=')[1].split()
            if 'NumberOfAtoms' in line:
                self.NumberOfAtoms = map(int, line.split('=')[1].split())
            if 'PickUp' in line:
                pickup = line.split('=')[1].strip()
                if pickup == 'T':
                    self.PickUp = True
                else:
                    self.PickUp = False
        #print self.NameOfAtoms, self.NumberOfAtoms

    def lpickup(self):
        if not self.PickUp:
            os.system('rm step')

    def submit_vasp(self):
        #jobid = 3
        #self.Write_input(jobid)
        os.system(self.CalyPsoPath)
        self.control_vasp()

    def autorun(self):
        if self.lwrite:
            self.write_gulp('./', 'POSCAR')
            sys.exit(0)

        self.readinput()

        if self.lgen:
            self.generate_structure()
            sys.exit(0)

        if self.lrelax:
            self.generate_structure()
            self.control_vasp()
            sys.exit(0)

        self.lpickup()

        for i in range(self.GenNum):
            self.submit_vasp()
            #sys.exit(0)

    def control_vasp(self):
        for i in range(self.StrNum):
            if not os.path.exists(r'./%s' % str(i+1)):
                #jobid = 1
                #self.Write_input(jobid)
                os.system('rm -rf %s' %str(i+1))
                os.system('mkdir %s' % str(i+1))
                os.system('ln -s  %s  %s/gap_parameters' % (self.paraspath,str(i+1)))

            self.write_gulp('./', 'POSCAR_%s' % str(i + 1))
            os.system('cp gulpinput    %s' % str(i + 1))
        nnn = int(self.StrNum/self.PN)
        for i in range(nnn):
            #for i in range(self.StrNum):
            for j in range(self.PN * i, self.PN * (i+1) ):
                os.system('cd %s;%s cd ..' % (str(j + 1),self.submit))
            jobtime = 0
            while True:
                if jobtime > self.MaxTime:
                    self.kill()
                    break
                    #print jobtime,self.MaxTime
                nump = int(os.popen('''ps aux | grep ngap.x | grep -c -v grep''').read().split()[0])
                if nump == 0:
                    break
                jobtime = jobtime + 2
                time.sleep(2)

        for i in range(self.StrNum):
            if os.path.exists(r'./%s/CONTCAR' % str(i+1)):
                os.system('cp %s/CONTCAR  CONTCAR_%s' % (str(i+1),str(i+1)))
            else:
                print('Optimization failure ./%s/CONTCAR not exit' % str(i+1))

            if os.path.exists(r'./%s/OUTCAR' % str(i+1)):
                os.system('cp %s/OUTCAR  OUTCAR_%s' % (str(i+1),str(i+1)))
            else:
                print('Optimization failure ./%s/OUTCAR not exit' % str(i+1))

    def write_gulp(self, dirname, finput):
        pl = 0.0
        pu = self.pstress
        (na, elements, lat, pos) = self.Read_Dposcar(finput)
        element = []
        element.append(elements[0])
        for u in elements:
            if u not in element:
                element.append(u)
          #print element
        ns = len(element)
        f = open(dirname + '/gulpinput' ,'w')
        f.write('opti conj conp nnpcal ocel\n')
        f.write('vectors\n')
        for i in range(3):
            f.write('%10.5f %10.5f %10.5f\n' % tuple(lat[i]))
        f.write('cartesian\n')
        for i in range(sum(na)):
            f.write(elements[i] + '  ' + 'core' + ' ')
            f.write('%10.5f %10.5f %10.5f\n' % tuple(pos[i]))
        f.write('space group\n')
        f.write('1 \n')
        if self.lrelax:
            x = np.random.random()
            pstress = pl + (pu - pl) * x
            f.write('pressure   %15.3f\n' % pstress)
        else:
            f.write('pressure    ' + str(self.pstress) +   '\n')
        f.write('species\n')
        for  u in element:
            f.write(u + ' core 0.0\n')
        f.write('sw2\n')
        for i in range(ns):
            for j in range(i,ns):
                f.write(element[i] + ' core ' + element[j] + ' core '\
                + '15.284982 2.0951 11.6031922831 0.0 3.77118' + '\n')
        f.write('dump every gulpopt\n')
        f.write('\n')
        f.write('maxcycle   '  + str(self.NSW) + ' \n')
        f.write('ftol   '  + str(self.ftol) + ' \n')
        f.write('gtol   '  + str(self.gtol) + ' \n')
        f.close()

    def kill(self):
        killjob = map(int,os.popen("ps aux | grep  ngap.x |grep -v grep | awk '{print $2}'").read().split())
        print(killjob)
        for id in killjob:
            os.system('kill -9 %d' % id)

    def Read_Dposcar(self, fa):
        lat = []
        pos = []
        line1 = []
        elements = []
        f1 = open(fa,'r')
        line1 = f1.readlines()
        f1.close()
        lDir = False
        try:
            na = list(map(int, line1[5].split()))
        except:
            del line1[5]
            na = list(map(int,line1[5].split()))

        name = os.popen('grep NameOfAtoms input.dat').read().split()[2:]
        try:
            NS = len(na)
        except:
            NS = len(list(na))

        if 'dir' in line1[6].split()[0].lower():
            lDir = True

        for a in range(3):
            lat.append(line1[2+a].split())
        for a in range(sum(na)):
            pos.append(line1[7+a].split()[0:3])
        lat = np.array(lat,float)
        pos = np.array(pos,float)
        if lDir:
            cpos = np.dot(pos,lat)
        else:
            cpos = pos

        for i in range(NS):
            for j in range(na[i]):
                elements.append(name[i])
        return na, elements, lat, cpos
    
    def check_input(self):
        calypsopath = './'
        paraspath = './'
        pstress = './'
        optpath = './'
        ftol = 0.0005
        gtol = 0.005
        pn = 12
        nsw = 200

        if self.lgen:
            return calypsopath, pn, paraspath, optpath, nsw, pstress, ftol, gtol

        if self.lwrite:
            tags = ['Pstress']
        else:
            tags = ['CALYPSO_PATH', 'GAPP_PATH', 'Optimizer_PATH', 'Pstress']
             
        line = []
        f = open('input.dat','r')
        ufile= f.readlines()
        f.close()

        for line in ufile:
            if line[0] == '#':
                continue
            if not self.lwrite:
                if 'CALYPSO_PATH' in line:
                    calypsopath = line.split('=')[1].strip()
                    if calypsopath == '':
                        print('\nPlease fill the CALYPSO_PATH in input.dat and run again!\n')
                        sys.exit(0)
                    tags.remove('CALYPSO_PATH')

                if 'GAPP_PATH' in line:
                    paraspath = line.split('=')[1].strip()
                    if paraspath == '':
                        print('\nPlease fill the GAPP_PATH in input.dat and run again!\n')
                        sys.exit(0)
                    tags.remove('GAPP_PATH')

                if 'Optimizer_PATH' in line:
                    optpath = line.split('=')[1].strip()
                    if optpath == '':
                        print('\nPlease fill the Optimizer_PATH in input.dat and run again!\n')
                        sys.exit(0)
                    tags.remove('Optimizer_PATH')

            if 'Pstress' in line:
                pstress = line.split('=')[1]
                try:
                    pstress = float(pstress)
                except:
                    print('\nPlease fill the Pstress in input.dat and run again!\n')
                    sys.exit(0)
                tags.remove('Pstress')

            if 'Ftol' in line:
                try:
                    ftol = float(line.split('=')[1])
                except:
                    continue

            if 'Gtol' in line:
                try:
                    gtol = float(line.split('=')[1])
                except:
                    continue

            if 'NumberOfParallel' in line:
                try:
                    pn = int(line.split('=')[1])
                except:
                    continue

            if 'Maxcycle' in line:
                try:
                    nsw = int(line.split('=')[1])
                except:
                    continue

        if len(tags) != 0:
            f = open('input.dat','w')
            for line in ufile:
                f.write(line)
            for tag  in tags:
                f.write('%12s = \n' % tag)
                print('\nPlease fill the [ %s ] in input.dat and run again! \n' % tag)
            f.close()
            sys.exit(0)

        return calypsopath, pn, paraspath, optpath, nsw, pstress, ftol, gtol

    def read_input(self):
        f = open('input.dat','r')
        xfile = f.readlines()
        f.close()
        try:
            '''
            For Python3
            '''
            yfile = xfile.copy()
        except:
            '''
            For Python2
            '''
            yfile = xfile[:]

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
    
            if 'CALYPSO_PATH' in line:
                CALYPSO_PATH = line.split('=')[1]
                lcaly = True
                xfile.remove(line)
    
        if not lcaly:
            print('\nPlease fill CALYPSO_PATH in input.dat\n')
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
    
    def write_input(self, xfile, popsize, numberofatoms, volume, name):
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
    
    def generate_structure(self):
        xfile, CM, vol, nstruct, ncomps, name, CALY_PATH= self.read_input()          
        popsize = int(nstruct/ncomps)
        CM = np.array(CM, int)
        vol = np.array(vol, float)
        nspecies, _ = CM.shape
        
        f = open('genstruc.log','w')
        #os.system('mkdir structure')
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
            self.write_input(xfile, popsize, numberofatoms, volume, name)
            os.system('mkdir %s; cp input2.dat %s/input.dat; cd %s; %s' % (str(icomp), str(icomp), str(icomp), CALY_PATH))
            for istruc in range(popsize):
                index += 1
                #if self.check_cell('%s/POSCAR_%s' % (str(icomp), str(istruc + 1))):
                os.system('cp %s/POSCAR_%s POSCAR_%s' % (str(icomp), str(istruc + 1), str(index)))
                #else:
                #print('alpha or beta or gamma are smaller than 30')
            os.system('rm -rf  ' +str(icomp))
        os.system('rm input2.dat')
        f.close()

    def cell_to_cellpar(self, cell, radians=False):
        """Returns the cell parameters [a, b, c, alpha, beta, gamma].
    
        Angles are in degrees unless radian=True is used.
        """
        lengths = [np.linalg.norm(v) for v in cell]
        angles = []
        for i in range(3):
            j = i - 1
            k = i - 2
            ll = lengths[j] * lengths[k]
            if ll > 1e-16:
                x = np.dot(cell[j], cell[k]) / ll
                angle = 180.0 / pi * arccos(x)
            else:
                angle = 90.0
            angles.append(angle)
        if radians:
            angles = [angle * pi / 180 for angle in angles]
        return np.array(lengths + angles)

    def check_cell(self, poscar):
        na, elements, lat, cpos = self.Read_Dposcar(poscar)
        cellpara = self.cell_to_cellpar(lat, False)
        if cellpara[3] < 30.0 or cellpara[4] < 30.0 or cellpara[5] < 30.0:
            return False    
        else:
            return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--relax', '-r', action='store_true', help='Generating random structures and relaxation by GAPs')
    parser.add_argument('--write', '-w', action='store_true', help='Converting POSCAR to gulp input file')
    parser.add_argument('--gen', '-g', action='store_true', help='Generating random structures')
    opt = parser.parse_args()

    print(opt.relax)
    print(opt.write)
    print(opt.gen)
    cgg2 = Autocalypso(lrelax = opt.relax, lwrite = opt.write, lgen = opt.gen)
    cgg2.autorun()
