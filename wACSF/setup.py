import sys
from numpy.distutils.core import Extension, setup


__author__ = "Qunchao Tong"
__copyright__ = ""
__credits__ = ["XXX"]
__license__ = "XXX"
__version__ = "1.00"
__maintainer__ = "Qunchao Tong"
__email__ = "tqc@calypso.cn"
__status__ = "Testing"
__description__ = ""
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
    COMPILER_FLAGS = ["-xHost", "-O3", "-axAVX", "-qopenmp"]
    LINKER_FLAGS = ["-liomp5", " -lpthread", "-lm", "-ldl"]
    MATH_LINKER_FLAGS = ["-L${MKLROOT}/lib/intel64", "-lmkl_rt"]


#mytest_module = Extension(name = 'libwacsf',
#                          sources = [
#                                './wacsf/constants.f90',
#                                './wacsf/math.f90',
#                                './wacsf/structure.f90',
#                                './wacsf/wacsf.f90',
#                            ],
#                          extra_f90_compile_args = COMPILER_FLAGS,
#                          extra_f77_compile_args = COMPILER_FLAGS,
#                          extra_compile_args = COMPILER_FLAGS ,
#                          #extra_link_args = LINKER_FLAGS + MATH_LINKER_FLAGS,
#                          language = FORTRAN,
#                          f2py_options=['--quiet'])

mytest_module = Extension(name = 'libwacsf',
                          sources = [
                                './libwacsf/wacsf.f90',
                            ],
                          extra_f90_compile_args = COMPILER_FLAGS,
                          extra_f77_compile_args = COMPILER_FLAGS,
                          extra_compile_args = COMPILER_FLAGS ,
                          #extra_link_args = LINKER_FLAGS + MATH_LINKER_FLAGS,
                          language = FORTRAN,
                          f2py_options=['--quiet'])

# use README.md as long description
def readme():
    with open('README.md') as f:
        return f.read()

def setup_pepytools():

    setup(

        name="libwacsf",
        packages=['libwacsf'],

        # metadata
        version=__version__,
        author=__author__,
        author_email=__email__,
        platforms = 'Any',
        description = __description__,
        long_description = readme(),
        keywords = ['Generative Adversarial Networks(GANs)', 'Quantum Chemistry'],
        classifiers = [],
        url = __url__,

        # set up package contents

        ext_package = 'libwacsf',
        ext_modules = [
              mytest_module,
        ],
)

if __name__ == '__main__':

    setup_pepytools()
