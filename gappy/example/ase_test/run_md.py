from ase import Atoms, Atom
from ase.io import write , read
from ase.calculators.gap_calc import GAP
from ase.md.nvtberendsen import NVTBerendsen
from ase.md import MDLogger
from ase.units import  fs
from run_opt import read_incar, write_contcar, get_ele_dict
import os, sys

calc = GAP(rcut = 6.0)

nsw, pstress, tebeg, potim, _ = read_incar()
my_atoms = read('POSCAR', format = 'vasp')
my_atoms.set_calculator(calc)
"""
for NVT 
"""
#print(tebeg,'tebeg')
dyn = NVTBerendsen(my_atoms, timestep=potim * fs, temperature=tebeg, taut=0.5*1000*fs, trajectory='ase.traj')
dyn.attach(MDLogger(dyn, my_atoms, 'OUTCAR', header=True, stress=True, pstress = pstress, peratom=False, mode="w"), interval=1)
dyn.run(steps=nsw)

atoms_lat = my_atoms.cell
atoms_pos = my_atoms.positions
atoms_symbols = my_atoms.get_chemical_symbols()
element, ele = get_ele_dict(atoms_symbols)
write_contcar(element, ele, atoms_lat, atoms_pos)
