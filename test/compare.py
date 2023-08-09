# this compares the outpur from the f90 code and the python code of PDFchem.
#%%
import numpy as np

f0 = np.loadtxt('output_f90_debug1.dat'); 
f1 = np.loadtxt('output_f90_cut.dat'); 
p0 = np.loadtxt('output_py_debug1.dat'); 
p1 = np.loadtxt('output_m_cut.dat'); 
p2 = np.loadtxt('output_py.dat'); 
for i in range(100):
    # print(f0[i])
    # print(f1[i])
    # print(p0[i])
    # print(p1[i])
    # print(p2[i])
    print((f0[i] - p0[i])/f0[i]) # percentage difference

# %%
