from ase.io import write , read
from libgap.GAP import Calculator

gap = Calculator(rcut = 6.0)
srbc = read('sps_all.xyz')
lgrad = True
pos =  srbc.positions
lat =  srbc.cell
species = srbc.get_atomic_numbers()
gap.gap_read()
ene, force, stress, _ = gap.gap_calc(species, lat, pos, lgrad)
print ene
print force
print stress
