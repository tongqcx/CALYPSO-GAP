"""Logging for molecular dynamics."""

import weakref
import sys
import ase.units as units
# ase.parallel imported in __init__

class MDLogger:
    """Class for logging molecular dynamics simulations.

    Parameters:
    dyn:           The dynamics.  Only a weak reference is kept.

    atoms:         The atoms.

    logfile:       File name or open file, "-" meaning standard output.

    stress=False:  Include stress in log.

    peratom=False: Write energies per atom.

    mode="a":      How the file is opened if logfile is a filename.
    
    This part were added by tqc, 2019.10.27
    pstress:       The target pressure.
                   Type: float

    lgap=False:    Simulation by Gaussian Approximation Potentials (GAPs)?.
                   Type: logical
    """
    def __init__(self, dyn, atoms, logfile, header=True, stress=False, pstress = 0.0,
                 peratom=False, mode="a", lgap=False):
        import ase.parallel
        if ase.parallel.rank > 0:
            logfile="/dev/null"  # Only log on master
        if hasattr(dyn, "get_time"):
            self.dyn = weakref.proxy(dyn)
        else:
            self.dyn = None
        self.atoms = atoms
        self.natoms = atoms.get_number_of_atoms()  # GLOBAL number of atoms
        if logfile == "-":
            self.logfile = sys.stdout
            self.ownlogfile = False
        elif hasattr(logfile, "write"):
            self.logfile = logfile
            self.ownlogfile = False
        else:
            self.logfile = open(logfile, mode, 1)
            self.ownlogfile = True
        self.stress = stress
        self.peratom = peratom
        self.pstress = pstress
        self.lgap = lgap
        if self.dyn is not None:
            self.hdr = "%-9s " % ("Time[ps]",)
            self.fmt = "%-10.4f "
        else:
            self.hdr = ""
            self.fmt = ""
        if self.peratom:
            self.hdr += "%12s %12s %12s  %6s" % ("Etot/N[eV]", "Epot/N[eV]",
                                                 "Ekin/N[eV]", "T[K]")
            self.fmt += "%12.4f %12.4f %12.4f  %6.1f"
        elif self.lgap:
            self.hdr += "%20s %20s %20s %15s  %9s" % ("Enth[eV/atom]","Etot[eV/atom]", "Epot[eV/atom]",
                                                 "Ekin[eV]", "T[K]")
            self.fmt += "%18.4f %18.4f %18.4f %18.4f %12.1f"
        else:
            self.hdr += "%12s %12s %12s  %6s" % ("Etot[eV]", "Epot[eV]",
                                                 "Ekin[eV]", "T[K]")
            # Choose a sensible number of decimals
            if self.natoms <= 100:
                digits = 4
            elif self.natoms <= 1000:
                digits = 3
            elif self.natoms <= 10000:
                digits = 2
            else:
                digits = 1
            self.fmt += 3*("%%12.%df " % (digits,)) + " %6.1f"
        if self.stress:
            self.hdr += "      ---------------------- stress [GPa] -----------------------"
            self.fmt += 6*" %10.3f"
        self.fmt += "\n"
        if header:
            self.logfile.write(self.hdr+"\n")
            
    def __del__(self):
        self.close()

    def close(self):
        if self.ownlogfile:
            self.logfile.close()

    def __call__(self):
        epot = self.atoms.get_potential_energy()
        ekin = self.atoms.get_kinetic_energy()
        temp = ekin / (1.5 * units.kB * self.natoms)

        volume = self.atoms.get_volume()
        enth = epot + volume * self.pstress / 1602.17733

        if self.peratom:
            epot /= self.natoms
            ekin /= self.natoms
        if self.dyn is not None:
            t = self.dyn.get_time() / (1000*units.fs)
            dat = (t,)
        else:
            dat = ()
        if self.lgap:
            dat += (enth/self.natoms, (epot+ekin)/self.natoms, epot/self.natoms, ekin, temp)
        else:
            dat += (epot+ekin, epot, ekin, temp)
        if self.stress:
            if self.lgap:
                dat += tuple(self.atoms.get_stress())
            else:
                dat += tuple(self.atoms.get_stress() / units.GPa)
                
        self.logfile.write(self.fmt % dat)
        self.logfile.flush()
        
