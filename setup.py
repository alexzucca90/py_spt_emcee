# setup for the diversity algorithms. Testing Phase
from numpy.distutils.core import setup, Extension

wrapper = Extension('python_spt.pyspt',
                    sources=['fortran/SPT_wrapper.f90', 'fortran/PySPT.f90'],
                    #extra_f90_compile_args=["-fPIC",  "-fopenmp", "-O3"],
                    extra_f90_compile_args=["-fPIC"],
                    libraries=['lapack'],
                    extra_link_args=['-lgomp'])
setup(
      ext_modules=[wrapper]
      )
