from ase.io import write , read
from libgap.GAP import Calculator
from libgap.WACSF import Wacsf

gap = Calculator(rcut = 6.0)
wacsf = Wacsf(nf = 66, rcut = 6.0, lgrad = False)
diamond = read('sps_all.xyz')
lgrad = True
pos =  diamond.positions
lat =  diamond.cell
species = diamond.get_atomic_numbers()
gap.gap_read()
ene, force, stress, _ = gap.gap_calc(species, lat, pos, lgrad)
acsf = wacsf.car2wacsf(lat, species, pos)
print acsf.shape, type(acsf)
print ene
print force
print stress
