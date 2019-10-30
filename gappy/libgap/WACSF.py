import numpy as np
import sys
try:
    '''
    For Python 3+
    '''
    from libgap.libgap  import fgap_read, fgap_calc
except:
    '''
    For Python 2+
    '''
    from libgap  import fgap_read, fgap_calc

class  Wacsf(object):

    def __init__(self, nf = None, rcut = 6.0, lgrad = True):

        self.nf = nf
        self.rcut = rcut
        self.lgrad = lgrad

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

    def car2wacsf(self, lat, elements, pos):
        """
        xx, dxdy, strs = fcar2wacsf(nf, lat, elements, pos, rcut, lgrad, na=shape(pos,0))
        """
        try:
            Elements = self.get_elenum(elements)
        except:
            Elements = elements
        xx, dxdy, strs = fcar2wacsf(self.nf, lat, Elements, pos, self.rcut, self.lgrad)
        return xx.T
