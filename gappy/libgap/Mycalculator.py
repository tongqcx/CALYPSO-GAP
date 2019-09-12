import numpy as np
import sys
try:
    from libgap  import fgap_read, fgap_calc
except:
    print 'ERROR: please install lib gap_calc'
    sys.exit(0)

class Gap_calculator(object):
    def __init__(self):

        #self.species = Species
        #self.lat = Lat
        #self.pos = Pos
        self.gap_read()

    def get_elenum(self,x):
        ele_dist = {'H':1,\
                    'Li':3,\
                    'B':5,\
                    'C':6,\
                    'O':8,\
                    'Mg':12,\
                    'Al':13,\
                    'Si':14,\
                    'P':15,\
                    'S':16,\
                    'Ca':20,\
                    'Ni':27,\
                    'Cs':55,\
                    }
        x_num = []
        for element in x:
            x_num.append(ele_dist[element])
        return x_num
    def gap_read(self): 
        self.sparse_cov_mat_len,\
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
                                                 self.mm[0:self.sparse_cov_mat_len,0:self.des_len],\
                                                 self.invcmm[0:self.sparse_cov_mat_len,0:self.sparse_cov_mat_len],\
                                                 self.coeff[0:self.sparse_cov_mat_len])
        return ene, force, stress, variance
