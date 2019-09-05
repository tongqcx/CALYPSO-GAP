from  libwacsf import  fcar2wacsf
class  WACSF(object):

    def __init__(self, rcut = 6.0, nfeature = None):

        self.rcut = rcut
        self.nfeature = nfeature

    def car2wacsf(self, lat, pos):
        """
        xx = fcar2wacsf(nf,lat,pos,rcut,na=shape(pos,0))
        """
        acsf = fcar2wacsf(self.nfeature, lat, pos, self.rcut)
        return acsf.T
