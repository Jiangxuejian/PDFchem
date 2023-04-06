#You can use the h5py library to create an HDF5 file and store your data in it. Here's an example code snippet that demonstrates how to do this:

#python
#Copy code

'''
the dataset can be read as:

hdf5_file['zcr00_uv00/params'][:]
hdf5_file['zcr00_uv00/pdr'][:]
'''
import h5py
import numpy as np

# Create a list of all the file names
models = ['Z0p1', 'Z0p5', 'Z1p0', 'Z2p0']

file_names = []
for i in range(41):
    for j in range(41):
        zcr_uv_params = f"zcr{i:02d}_uv{j:02d}.params"
        zcr_uv_pdr_fin = f"zcr{i:02d}_uv{j:02d}.pdr.fin"
        zcr_uv_spop_fin = f"zcr{i:02d}_uv{j:02d}.spop.fin"
        file_names.append((zcr_uv_params, zcr_uv_pdr_fin, zcr_uv_spop_fin))

for model in models:
    # Create an HDF5 file
    hdf5_file = h5py.File(model+'.hdf5', 'w')

    # Loop over the file names and store the data in the HDF5 file
    for file_name in file_names:
        # Load the data from the files
        params_data = np.loadtxt(model + '/' + file_name[0])
        pdr_data = np.loadtxt(model + '/' + file_name[1])
        spop_data = np.loadtxt(model + '/' + file_name[2])

        # Create a group in the HDF5 file for this file
        group_name = f"{file_name[0].split('.')[0]}"
        group = hdf5_file.create_group(group_name)

        # Store the data in the group as datasets
        group.create_dataset('params', data=params_data, compression="gzip", compression_opts=9)
        group.create_dataset('pdr', data=pdr_data, compression="gzip", compression_opts=9)
        group.create_dataset('spop', data=spop_data, compression="gzip", compression_opts=9)

    # Close the HDF5 file
    hdf5_file.close()

'''
In the above code, we first create a list of all the file names by looping over the possible values of i and j. Each file name is a tuple containing the name of the .params file and the name of the .pdr.fin file.

We then create an HDF5 file using the h5py.File() function with mode w for write access.

Next, we loop over the file names and for each file we load the data from the .params and .pdr.fin files using the numpy.loadtxt() function. We then create a group in the HDF5 file with the same name as the .params file (without the extension) using the h5py.Group.create_group() method.

Finally, we store the data in the group as datasets using the h5py.Group.create_dataset() method, with names params and pdr for the .params
'''