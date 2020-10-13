PySPT 
======
This is a code to perform MCMC on SPT likelihood with PMF.


First of all, install the requirements with 
```bash
pip install -r requirements.txt
```

Run
```bash
python test_emcee.py
```

and check that the package emcee works.

Then, install the SPT python wrapper with :
```bash
f2py -llapack -c fortran/SPT_wrapper.f90 fortran/PySPT.f90 -m python_spt.pyspt
```

Finally, modify the parameters in `emc_spt.py` and run the code,
```bash
python emc_spt.py
```


