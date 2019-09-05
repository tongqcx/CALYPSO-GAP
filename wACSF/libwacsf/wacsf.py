try:        
    from  wacsf import  fcar2wacsf 
except:  
    print 'Please install libwacsf'
    sys.exit(0)

class  WACSF(object):


    def __init__(self, rcut = 6.0, nfeature = None):

        self.rcut = rcut
        self.nfeature = nfeature

    def car2wacsf(self, lat, pos)
    """
    xx = fcar2wacsf(nf,lat,pos,rcut,na=shape(pos,0))
    """
        return fcar2wacsf(self.nfeature, lat, pos, self.rcut)
