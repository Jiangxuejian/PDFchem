#%%
import numpy as np
import math
import os, bisect
import datetime
#from decimal import Decimal
import h5py
from multiprocessing import Pool, Manager

np.seterr(divide='ignore', invalid='ignore')

pc2cm = 3.0856776e18
mhp = 1.6726219e-24
freq0 = np.zeros(17)
g = np.zeros((17,2))
A = np.zeros(17)
mh = np.zeros(17)
spec = [''] * 17
phi = np.zeros(17)
sigma = np.zeros(17)
tau_incr = np.zeros(17)
frac = np.zeros(17)
t_r = np.zeros(17)
t_a_i = np.zeros(17)
t_a_ip1 = np.zeros(17)
dtau = np.zeros(17)
Ntr = np.zeros(17)

# for making prefix
ipref = 0
num = [''] * 2
pref = [''] * 2501
Nspec = 33  
vturb = 1e6 # microturbulent velocity in cm/s
            # NOTE: all given PDR calculations used 1~km/s. 

# Constants for calculations
c = 2.9979246e10 # cm/s
kb = 1.380650e-16   # erg / K or g cm^2 / K s^2
hp = 6.6260696e-27  # erg s or g cm^2 / s
pi = math.pi        # 3.1415927
mhp = 1.6726218e-24 # g

g = np.zeros((Nspec, 2))
spec = [""] * Nspec

freq0 = np.array([1900.5369e9, 492.16065e9, 1301.50262e9, 809.34197e9, 4744.77749e9, 6804.84658e9, 2060.06909e9,
115.2712018e9, 230.538e9, 345.7959899e9, 461.040768e9, 576.2679305e9, 691.4730763e9, 806.6518060e9,
921.7997000e9, 1036.9123930e9, 1151.9854520e9])

g = np.array([[2.0, 4.0], [1.0, 3.0], [1.0, 5.0], [3.0, 5.0], [5.0, 3.0], [5.0, 1.0], [3.0, 1.0], [1.0, 3.0],
[3.0, 5.0], [5.0, 7.0], [7.0, 9.0], [9.0, 11.0], [11.0, 13.0], [13.0, 15.0], [15.0, 17.0], [17.0, 19.0], [19.0, 21.0]])

spec = ["CII 158um", "CI (1-0)", "CI (2-0)", "CI (2-1)", "OI 1-0 ", "OI 2-0 ", "OI 2-1 ", "CO (1-0)",
"CO (2-1)", "CO (3-2)", "CO (4-3)", "CO (5-4)", "CO (6-5)", "CO (7-6)", "CO (8-7)", "CO (9-8)", "CO (10-9)"]    

# Einstein A coefficients and mass of hydrogen atoms
A = np.array([2.321e-06, 7.880e-08, 1.810e-14, 2.650e-07, 8.910e-05, 1.340e-10, 1.750e-05, 7.203e-08, 6.910e-07,
2.497e-06, 6.126e-06, 1.221e-05, 2.137e-05, 3.422e-05, 5.134e-05, 7.330e-05, 1.006e-04])
mh = np.array([12.*mhp, 12.*mhp, 12.*mhp, 12.*mhp, 16.*mhp, 16.*mhp, 16.*mhp, 28.*mhp, 28.*mhp, 28.*mhp, 28.*mhp,
28.*mhp, 28.*mhp, 28.*mhp, 28.*mhp, 28.*mhp, 28.*mhp])

def screen_time():
    return datetime.datetime.now().strftime("\033[0;32m%H:%M:%S %Y-%m-%d > \33[0m")

def readfile(prefix, hdf5_file):
    model = h5py.File(hdf5_file, 'r')
    #directory = directory[:4]
    fuv, cosmicrays, Z = model[f'{prefix}/params'][:]
    # print(f'{hdf5_file} read in #2. Z = {Z}')

    pdr = model[f'{prefix}/pdr']
    itot = np.shape(pdr)[0] # 490

    #id =  np.int16(pdr[:,0])
    x =    pdr[:,1]
    tgas = pdr[:,3]
    tdust= pdr[:,4] 
    #etype = np.int16(pdr[:,5])
    nh =   pdr[:,6]
    abun = pdr[:,8:]

    spop = model[f'{prefix}/spop'][:].T       # 490, 28
    pop = np.zeros((17, 2, itot))
    tau = np.zeros((17, itot))

    dummy1, dummy2,  pop[0, 0], pop[0, 1], dummy3, dummy4, dummy5, pop[1, 0], pop[1, 1], pop[2, 1],\
    dummy6, dummy7, pop[4, 0], pop[4, 1], pop[5, 1], dummy8, dummy9, \
    pop[7, 0], pop[7, 1], pop[8, 1], pop[9, 1], pop[10, 1], \
    pop[11, 1], pop[12, 1], pop[13, 1], pop[14, 1], pop[15, 1], pop[16, 1] = spop[:]

    # record pairs
    pop[2, 0, :] = pop[1, 0, :]; pop[3, 1, :] = pop[2, 1, :]; pop[3, 0, :] = pop[1, 1, :]
    pop[5, 0, :] = pop[4, 0, :]; pop[6, 1, :] = pop[5, 1, :]; pop[6, 0, :] = pop[4, 1, :]
    pop[8, 0, :] = pop[7, 1, :]; pop[9, 0, :] = pop[8, 1, :]; pop[10, 0, :] = pop[9, 1, :]
    pop[11, 0, :] = pop[10, 1,:]; pop[12, 0, :] = pop[11, 1,:]; pop[13, 0, :] = pop[12, 1,:]
    pop[14, 0, :] = pop[13, 1,:]; pop[15, 0, :] = pop[14, 1,:]; pop[16, 0, :] = pop[15, 1,:]
    
    model.close()
    return fuv,cosmicrays,Z, tgas, tdust, abun, x, nh, itot, pop, tau

def iteration(ipref, pref, hdf5_file, avtot, avin, pdfin, output_strarray):
    global tau_incr
    factor1 = A*c**2/8./pi/freq0**2 # to speed up the calculation
    factor2 = c**2/2/kb/freq0**2
    prefix = pref[ipref]

    # calculate column densities of species
    N = np.zeros(Nspec + 1)
    Ntgas, Ntot_nopdf, Nrho = 0, 0, 0
    fuv,cosmicrays,Z, tgas, tdust, abun, x, nh, itot, pop, tau = readfile(prefix, hdf5_file) # call readfile function
    print(f'{hdf5_file} read in.')
    freq1 = np.tile(freq0[:, np.newaxis], (1, itot-1))
    mh1 = np.tile(mh[:, np.newaxis], (1, itot-1))
    tgas1 = np.tile(tgas[np.newaxis,:], (17, 1))
    g1 = np.tile(np.expand_dims(g.T, axis=2), (1,1, itot-1))
    Tex = np.zeros((17, itot))
    Bnu = np.zeros((17, itot))

    # Excitation temperatures and Background Radiation for radiative transfer
    mask = (pop[:,1,:-1] == 0) | (np.abs(g1[1,:,:]*pop[:,0,:-1]/pop[:,1,:-1]/g1[0,:,:]-1) < 1e-2)
    Tex = Tex[:,:-1]
    Tex[mask] = 0
    Bnu1 = Bnu[:,:-1]
    Bnu1[mask] = 0
    pop0 = pop[:,0,:-1]
    pop1 = pop[:,1,:-1]
    Tex[~mask] = (hp*freq1[~mask]/kb)/np.log(g1[1,~mask]*pop0[~mask]/pop1[~mask]/g1[0,~mask])
    Bnu1[~mask] = (2.*hp*freq1[~mask]**3/c**2)/(np.exp(hp*freq1[~mask]/kb/Tex[~mask])-1) # Black body emission

    sigma = (freq1/c)*np.sqrt(kb*tgas1[:,:-1]/mh1+vturb**2/2.) # sigma value used in phi
    phi = 1/sigma/math.sqrt(2*pi) # phi value used in optical depth
    frac = 0.5*((pop[:,0,:-1]+pop[:,0,1:])*g1[1]/g1[0]-(pop[:,1,:-1]+pop[:,1,1:]))
    step = np.abs(x[:-1] - x[1:]) * pc2cm

    Avobs = 0.06*(nh[:-1]**0.69)        # Av,obs -- nH relation
    # print('\n')
    for i in range(itot-1):
        k = bisect.bisect_left(avin, Avobs[i])  
        if k >= avtot: k = avtot-1
        factor3 = step[i]*pdfin[k]
        if Avobs[i] <= avin[k]:
            N[0] += 0.5*(nh[i]+nh[i+1])*factor3 # total column density
            N[1:Nspec+1] += 0.5*(nh[i]*abun[i,:] + nh[i+1]*abun[i+1,:])*factor3 # species
            Ntgas += 0.5*(nh[i]*tgas[i] + nh[i+1]*tgas[i+1])*factor3 # for <Tgas>
            Nrho += 0.5*(nh[i]**2*abun[i,30] + nh[i+1]**2*abun[i+1,30])*factor3
        # print(i, k, nh[i], Avobs[i], avin[k], pdfin[k], N[0], Ntgas, Nrho)
        tau_incr += phi[:,i]*( factor1 )*frac[:,i]*step[i] # optical depth calculation
        tau[:,i] = tau_incr # record value
    # solve radiative transfer equation
    tau_cii, tau_ci, tau_co = 0, 0, 0

    Ntr, Ntot, t_r, Ncol = np.zeros((17,)), 0, np.zeros((17,)), 0
    step = np.abs(x[:-2] - x[1:-1]) * pc2cm
    for i in range(itot-2):
        dtau = tau[:, i+1] - tau[:, i]
        mask = dtau > 1e10 
        t_r[mask] = Bnu[mask, i]
        mask = (dtau > 1e-6) & (dtau <= 1e10)
        t_a_i = Bnu[mask, i]*((1-np.exp(-dtau[mask]))/dtau[mask]-np.exp(-dtau[mask]))
        t_a_ip1 = Bnu[mask, i+1]*(1-(1-np.exp(-dtau[mask]))/dtau[mask])
        t_r[mask] = t_r[mask]*np.exp(-dtau[mask])+t_a_i+t_a_ip1
        mask = dtau <= 1e-6
        t_r[mask] = t_r[mask]*(1-dtau[mask])+(Bnu[mask, i]+Bnu[mask, i+1])*dtau[mask]/2

        k = bisect.bisect_left(avin, Avobs[i])  
        if k >= avtot: k = avtot-1
        if Avobs[i] <= avin[k]:
            Ntot = Ntot + pdfin[k]
            Ntr = Ntr + (t_r*factor2)*pdfin[k]
            Ncol = Ncol + (0.5*(nh[i]+nh[i+1])*step[i])
            tau_cii = tau_cii + dtau[0]*pdfin[k]
            tau_ci  = tau_ci  + dtau[1]*pdfin[k]
            tau_co  = tau_co  + dtau[7]*pdfin[k]

    sublist = output_strarray[ipref] + [f'{fuv:11.2e}{cosmicrays:11.2e}{Z:11.2e}{Ntgas/N[0]:11.2e}']
    for i in range(1, Nspec+1):
        sublist = sublist + [f'{N[i]/N[0]:11.2e}']
    sublist = sublist + [f'{Ntr[0]:11.2e}{Ntr[1]:11.2e}{Ntr[3]:11.2e}']
    for i in range(7,17):
        sublist = sublist + [f'{Ntr[i]:11.2e}']
    sublist = sublist + ['\n']
    output_strarray[ipref] = sublist

# Av,obs - PDF function
def pdf(x, s, m):
    return (1.0 / (s * math.sqrt(2.0 * pi))) * math.exp(-(math.log(x) - m) ** 2 / (2.0 * s ** 2))


# Make Av-PDF
def makepdf(av_bar, s):
    m = math.log(av_bar) - s ** 2 / 2.0
    avtot = 500
    avin = [0] * (avtot + 1)
    pdfin = [0] * (avtot + 1)
    with open('avpdf.dat', 'w') as f:
        for i in range(avtot + 1):
            lav = -2.0 + 4.0 * float(i) / float(avtot + 1)
            avin[i] = 10.0 ** lav
            pdfin[i] = pdf(float(avin[i]), s, m)
            if pdfin[i] < 1e-10: pdfin[i] = 1e-10
            f.write(str(avin[i]) + ' ' + str(pdfin[i]) + '\n')
    return avin,pdfin,avtot

# Read inputted Av-PDF
def readpdf(avpdf_file):
    f = np.loadtxt(avpdf_file)
    avtot = len(f)
    avin = f[:,0]
    pdfin = f[:,1]
    pdfin[pdfin<1e-10] = 1e-10
    return avin,pdfin,avtot

def NH_fits2pdf(fits_file, output_path):
    from astropy.io import fits
    from astropy.stats import histogram
    bins=100
    with fits.open(fits_file) as hdul:
        NH = hdul[0].data
    lg_av = np.log10(NH) + np.log10(6.3e-22)
    lg_av_p, bin_edge = histogram(lg_av, bins=bins, density=True)
    bin = bin_edge[:-1] + np.diff(bin_edge)
    # return 10**bin, lg_av_p 
    np.savetxt(f'{output_path}/avpdf_input.dat', np.array((10**bin, lg_av_p)).T)
    return f'{output_path}/avpdf_input.dat'

def main(avpdf_file='', output_file='output.dat'):
    # Make prefix for all inputs
    def makeprefix():
        num = [str(i).zfill(2) for i in range(61)]
        pref = []
        for i in range(41):
            for j in range(61):
                pref.append("zcr" + num[i] + "_uv" + num[j])
        return pref

    # Input parameters
    av_bar, s, metallicity = np.loadtxt('./pdfchem.params')
    if metallicity == 0.1:   directory = 'Z0p1/'
    elif metallicity == 0.5: directory = 'Z0p5/'
    elif metallicity == 1.0: directory = 'Z1p0/'
    elif metallicity == 2.0: directory = 'Z2p0/'

    pref = makeprefix()
    hdf5_file = f'models/{directory[:4]}.hdf5'
    # Check if user has an input AV-PDF file. 
    # If not, call makepdf to create a simulated PDF.
    if os.path.isfile(avpdf_file):
        avin, pdfin,avtot = readpdf(avpdf_file)
    else:
        avin, pdfin,avtot = makepdf(av_bar, s)

    num_iter = 2501

    # print(f"\r{screen_time()} Looping through {num_iter} Simulations")

    # Create a shared list for output_strarray
    # if parallel:
    with Manager() as manager:
        output_strarray = manager.list([[''] for _ in range(num_iter)])
        with Pool() as pool:
            pool.starmap(iteration, [(ipref, pref, hdf5_file, avtot, avin, pdfin, output_strarray) for ipref in range(num_iter)])
        output_strarray = list(output_strarray)
    # elif parallel == False:
    #     output_strarray = list([[''] for _ in range(num_iter)])
    #     for ipref in range(num_iter):
    #         print(f"{ipref}")
    #         iteration(ipref, pref, hdf5_file, avtot, avin, pdfin, output_strarray)
    #     output_strarray = list(output_strarray)
        # print(output_strarray[0])

    outfile = open(output_file, 'w')
    for sublist in output_strarray:
        outfile.write(''.join(map(str, sublist)))
    outfile.close()
    
    # print(f'{screen_time()} Finished!')

if __name__ == '__main__':
    main()
