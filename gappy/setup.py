import sys
from numpy.distutils.core import Extension, setup


__author__ = "Qunchao Tong"
__copyright__ = "Unknown"
__credits__ = ["XXX"]
__license__ = "XXX"
__version__ = "XXX"
__maintainer__ = "Qunchao Tonf"
__email__ = "Qunchao Tong"
__status__ = "Beta"
__description__ = "Calculating total energy, atomic force, cell stress by G(aussian)A(pproximation)P(optential)"
__url__ = "XXX"


FORTRAN = "f90"

# GNU (default)
COMPILER_FLAGS = ["-O3", "-fopenmp", "-m64", "-march=native", "-fPIC",
                    "-Wno-maybe-uninitialized", "-Wno-unused-function", "-Wno-cpp"]
LINKER_FLAGS = ["-lgomp"]
MATH_LINKER_FLAGS = ["-lblas", "-llapack"]


# For clang without OpenMP: (i.e. most Apple/mac system)
if sys.platform == "darwin" and all(["gnu" not in arg for arg in sys.argv]):
    COMPILER_FLAGS = ["-O3", "-m64", "-march=native", "-fPIC"]
    LINKER_FLAGS = []
    MATH_LINKER_FLAGS = ["-lblas", "-llapack"]


# Intel
if any(["intelem" in arg for arg in sys.argv]):
    COMPILER_FLAGS = ["-xHost", "-O3" , "-fopenmp"]
    LINKER_FLAGS = ["-liomp5", " -lpthread", "-lm", "-ldl"]
    MATH_LINKER_FLAGS = ["-L${MKLROOT}/lib/intel64", "-lmkl_rt"]


mytest_module = Extension(name = 'libgap',
                          sources = [
                                './libgap/wacsf.f90',
                                './libgap/gap_init.f90',
                                './libgap/gap_calc.f90',
                            ],
                          extra_f90_compile_args = COMPILER_FLAGS,
                          extra_f77_compile_args = COMPILER_FLAGS,
                          extra_compile_args = COMPILER_FLAGS ,
                          extra_link_args = LINKER_FLAGS + MATH_LINKER_FLAGS,
                          language = FORTRAN,
                          f2py_options=['--quiet'])

# use README.md as long description
def readme():
    with open('README.md') as f:
        return f.read()

def setup_pepytools():

    setup(

        name="libgap",
        packages=['libgap'],

        # metadata
        version=__version__,
        author=__author__,
        author_email=__email__,
        platforms = 'Any',
        description = __description__,
        long_description = readme(),
        keywords = ['Machine Learning', 'Quantum Chemistry'],
        classifiers = [],
        url = __url__,

        # set up package contents

        ext_package = 'libgap',
        ext_modules = [
              mytest_module,
        ],
)

if __name__ == '__main__':

    setup_pepytools()
