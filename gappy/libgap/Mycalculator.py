import numpy as np
import sys
try:
    from libgap  import fgap_read, fgap_calc
except:
    print 'ERROR: please install lib gap_calc'
    sys.exit(0)

class Gap_calculator(object):
    def __init__(self):

        self.gap_read()

    def get_elenum(self,x):
        ChemicalSymbols = [ 'X',  'H',  'He', 'Li', 'Be','B',  'C',  'N',  'O',  'F',
                            'Ne', 'Na', 'Mg', 'Al', 'Si','P',  'S',  'Cl', 'Ar', 'K',
                            'Ca', 'Sc', 'Ti', 'V',  'Cr','Mn', 'Fe', 'Co', 'Ni', 'Cu',
                            'Zn', 'Ga', 'Ge', 'As', 'Se','Br', 'Kr', 'Rb', 'Sr', 'Y',
                            'Zr', 'Nb', 'Mo', 'Tc', 'Ru','Rh', 'Pd', 'Ag', 'Cd', 'In',
                            'Sn', 'Sb', 'Te', 'I',  'Xe','Cs', 'Ba', 'La', 'Ce', 'Pr',
                            'Nd', 'Pm', 'Sm', 'Eu', 'Gd','Tb', 'Dy', 'Ho', 'Er', 'Tm',
                            'Yb', 'Lu', 'Hf', 'Ta', 'W','Re', 'Os', 'Ir', 'Pt', 'Au',
                            'Hg', 'Tl', 'Pb', 'Bi', 'Po','At', 'Rn', 'Fr', 'Ra', 'Ac',
                            'Th', 'Pa', 'U',  'Np', 'Pu','Am', 'Cm', 'Bk', 'Cf', 'Es',
                            'Fm', 'Md', 'No', 'Lr']
        atomicNum = {}
        for anum, symbol in enumerate(ChemicalSymbols):
            atomicNum[symbol] = anum

        x_num = []
        for element in x:
            x_num.append(atomicNum[element])
        return x_num

    def gap_read(self): 
        self.nsparseX,\
        self.des_len,\
        self.theta,\
        self.mm,\
        self.invcmm,\
        self.coeff = fgap_read()
    def gap_calc(self,species, lat, pos):
        try:
            Species = self.get_elenum(species)
        except:
            Species = species
        ene, force, stress, variance = fgap_calc(Species,\
                                                 lat,\
                                                 pos,\
                                                 self.theta[0:self.des_len],\
                                                 self.mm[0:self.nsparseX,0:self.des_len],\
                                                 self.invcmm[0:self.nsparseX,0:self.nsparseX],\
                                                 self.coeff[0:self.nsparseX])
        return ene, force, stress, variance
