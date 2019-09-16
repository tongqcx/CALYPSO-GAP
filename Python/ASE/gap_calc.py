"""
=========================
Author: Qunchao Tong
Email: tqc@calypso.cn
Date: 2019.01.04
Description: Calculating total energy by GAP
"""
from __future__ import division

import numpy as np

from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError
from libgap  import  GAP


class GAP(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']

    nolabel = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        natoms = len(self.atoms)

        energy = 0.0
        forces = np.zeros((natoms, 3))
        stress = np.zeros((6))
        variance = 0.0

        gap_calc = GAP.Calculator()
        species = self.atoms.get_atomic_numbers()
        lat = self.atoms.cell
        pos = self.atoms.positions
        energy, forces, stress, variance = gap_calc.gap_calc(species, lat, pos)

        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces
        self.results['stress'] = stress
        self.results['variance'] = variance
