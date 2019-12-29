from ase import Atoms, Atom
from ase.io import write , read
from ase.optimize import QuasiNewton,LBFGS
from ase.calculators.gap_calc import GAP
import sys, os
import numpy as np

def get_ele_dict(elements):
    element = []
    ele = {}
    element.append(elements[0])
    for x  in elements:
        if x not in element :
            element.append(x)
    for x in element:
        ele[x] = elements.count(x)
    return element, ele

def write_contcar(element, ele, lat, pos):
    f = open('CONTCAR','w')
    f.write('ASE-GAP-Optimization-Fixcell\n')
    f.write('1.0\n')
    for i in range(3):
        f.write('%15.10f %15.10f %15.10f\n' % tuple(lat[i]))
    for x in element:
        f.write(x + '  ')
    f.write("\n")
    for x in element:
        f.write(str(ele[x]) + '  ')
    f.write("\n")
    f.write('Direct\n')
    na = sum(ele.values())
    dpos = np.dot(pos,np.linalg.inv(lat))
    for i in range(na):
        f.write('%15.10f %15.10f %15.10f\n' % tuple(dpos[i]))

def write_outcar(element, ele, volume, lat, pos, ene, force, stress, pstress):
    f = open('OUTCAR','w')
    for x in element:
        f.write('VRHFIN =' + str(x) + '\n')
    f.write('ions per type =')
    for x in element:
        f.write('%5d' % ele[x])
    f.write('\nvolume of cell :\n')
    f.write('Direction     XX             YY             ZZ             XY             YZ             ZX\n')
    f.write('in kB')
    f.write('%15.6f' % stress[0]) 
    f.write('%15.6f' % stress[1]) 
    f.write('%15.6f' % stress[2]) 
    f.write('%15.6f' % stress[3]) 
    f.write('%15.6f' % stress[4]) 
    f.write('%15.6f' % stress[5]) 
    f.write('\n')
    ext_pressure = np.sum(stress[0] + stress[3] + stress[5])/3.0 - pstress
    f.write('external pressure = %20.6f kB    Pullay stress = %20.6f  kB\n' % \
           (ext_pressure, pstress))
    f.write('volume of cell : %20.6f\n' % volume)
    f.write('direct lattice vectors\n')
    for i in range(3):
        f.write('%20.6f %20.6f %20.6f\n' % tuple(lat[i]))
    f.write('POSITION                                       TOTAL-FORCE(eV/Angst)\n')
    f.write('-------------------------------------------------------------------\n')
    na = sum(ele.values())
    for i in range(na):
        f.write('%20.6f %20.6f %20.6f' % tuple(pos[i]))
        f.write('%20.6f %20.6f %20.6f\n' % tuple(force[i]))
    f.write('-------------------------------------------------------------------\n')
    f.write('energy  without entropy= %20.6f %20.6f\n' % (ene, ene/na))
    enthalpy = ene + pstress * volume / 1602.17733
    f.write('enthalpy is  TOTEN    = %20.6f %20.6f\n' % (enthalpy, enthalpy/na))

def read_incar():
    nsw = 1
    pstress = 0.0
    tebeg = 300.0 
    potim = 1.0
    assert os.path.exists('./INCAR'), 'INCAR does not exist!'
    f = open('INCAR','r')
    xfile = f.readlines()
    f.close()
    for line in xfile:
        if line[0] == '#':
            continue
        if 'NSW' in line:
            nsw = int(line.split('=')[1])
        if 'PSTRESS' in line:
            pstress = float(line.split('=')[1])
        if 'TEBEG' in line:
            tebeg = float(line.split('=')[1])
        if 'POTIM' in line:
            potim = float(line.split('=')[1])
        if 'EDIFFG' in line:
            ediffg = float(line.split('=')[1])
    ediffg = ediffg * -1.0
    return nsw, pstress, tebeg, potim, ediffg


if __name__ == '__main__':

    '''
    stress format: XX XY XZ YY YZ ZZ 
    units: kB
    '''
    calc = GAP(rcut = 6.0)
    nsw, pstress, _, _,  ediffg = read_incar()
    my_atoms = read('POSCAR', format = 'vasp')
    my_atoms.set_calculator(calc)
    dyn = LBFGS(my_atoms, trajectory='ase.traj')
    dyn.run(fmax = ediffg, steps = nsw)
    
    atoms_lat = my_atoms.cell
    atoms_pos = my_atoms.positions
    atoms_force = my_atoms.get_forces()
    atoms_stress = my_atoms.get_stress()  
    atoms_num = my_atoms.get_atomic_numbers()
    atoms_symbols = my_atoms.get_chemical_symbols()
    atoms_formula = my_atoms.get_chemical_formula()
    atoms_ene = my_atoms.get_potential_energy()
    atoms_vol = my_atoms.get_volume()
    element, ele = get_ele_dict(atoms_symbols)
    
    write_contcar(element, ele, atoms_lat, atoms_pos)
    write_outcar(element, ele, atoms_vol, atoms_lat, atoms_pos, \
                 atoms_ene, atoms_force, atoms_stress * 10.0, pstress)
