# this compares the outpur from the f90 code and the python code of PDFchem.
#%%
import numpy as np

f = np.loadtxt('output_f90.dat'); 
p = np.loadtxt('output_py.dat'); 
for i in range(2501):
    print((p[i] - f[i])) # percentage difference

# %%
