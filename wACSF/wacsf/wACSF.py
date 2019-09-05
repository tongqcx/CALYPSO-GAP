
try:        
    from  wacsf import  fcar2wacsf
except:  
    print 'Please install libwacsf'
    sys.exit(0)

class  WACSF(object):

    def __init__(self):
        pass

    def cartesian2wacsf(self, elemetns, lat, pos)
        return fcar2wacsf(elemetns, lat, pos)
