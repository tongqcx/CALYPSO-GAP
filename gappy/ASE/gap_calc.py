"""
=========================
Author: Qunchao Tong
Email: tqc@calypso.cn
Date: 2019.01.04
Modified: 2019.10.25
Description: Calculating total energy and grad information by GAP
"""
from __future__ import division

import numpy as np

from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError
import libgap.GAP  as  my_gap


class GAP(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']

    nolabel = True

    #def __init__(self, **kwargs):
    def __init__(self, rcut = 6.0):
        Calculator.__init__(self)
        self.rcut = rcut

    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        natoms = len(self.atoms)

        energy = 0.0
        forces = np.zeros((natoms, 3))
        stress = np.zeros((6))
        variance = 0.0

        gap_calc = my_gap.Calculator(rcut = self.rcut)
        species = self.atoms.get_atomic_numbers()
        lat = self.atoms.cell
        pos = self.atoms.positions
        lgrad = True
        energy, forces, stress, variance = gap_calc.gap_calc(species, lat, pos, lgrad)

        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces
        self.results['stress'] = stress
        self.results['variance'] = variance
