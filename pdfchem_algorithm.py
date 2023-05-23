def main():

    import numpy as np
    import math
    import os, bisect
    import datetime
    from decimal import Decimal
    
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
    file_exists = False

    # def constants():
    # Constants for calculations
    c = 2.9979246e10 # cm/s
    kb = 1.380650e-16 # erg / K or g cm^2 / K s^2
    hp = 6.6260696e-27 # erg s or g cm^2 / s
    pi = 3.1415927
    mhp = 1.6726218e-24 # g
    Nspec = 33  

    vturb = 1e6 # microturbulent velocity in cm/s
                # NOTE: all given PDR calculations used 1~km/s. 

    freq0 = np.zeros(Nspec)
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


    # Call constants
    # constants()

    def readfile(directory, prefix, Nspec):
        fileparams = directory.strip() + prefix.strip() + '.params' # parameters
        with open(fileparams, 'r') as f:
            fuv, cosmicrays, Z = map(float, f.readline().split())

        filein = directory.strip() + prefix.strip() + '.pdr.fin' # abundances
        filepdr = directory.strip() + prefix.strip()
        f = open(filein, 'r')
        itot = 0
        for line in f:
            itot += 1
        f.seek(0)
        tgas = np.zeros(itot)
        tdust = np.zeros(itot)
        # abun = [[0.0 for _ in range(itot)] for _ in range(Nspec)]
        abun = np.zeros((itot,Nspec))
        N = np.zeros(Nspec + 1)
        x = np.zeros(itot)
        nh = np.zeros(itot)
        for i in range(itot):
            line = f.readline().split()
            id = int(line[0])
            x[i] = float(line[1])
            tgas[i] = float(line[3])
            tdust[i] = float(line[4])
            etype = int(line[5])
            nh[i] = float(line[6])
            abun[i,:] = np.array(line[8:], dtype=float)
            # for j in range(1, Nspec):
            #     abun[j][i] = float(f.readline().split()[1])
        f.close()

        filepop = directory.strip() + prefix.strip() + '.spop.fin' # level populations
        pop = np.zeros((17, 2, itot))
        tau = np.zeros((17, itot))
        with open(filepop, 'r') as file:
            for i in range(itot):
                line = file.readline().split()
                dummy1, dummy2, pop[0, 0, i], pop[0, 1, i], dummy3, dummy4, dummy5, \
                pop[1, 0, i], pop[1, 1, i], pop[2, 1, i], dummy6, dummy7, \
                pop[4, 0, i], pop[4, 1, i], pop[5, 1, i], dummy8, dummy9, \
                pop[7, 0, i], pop[7, 1, i], pop[8, 1, i], pop[9, 1, i], pop[10, 1, i], \
                pop[11, 1, i], pop[12, 1, i], pop[13, 1, i], pop[14, 1, i], pop[15, 1, i], pop[16, 1, i] = line
                # record pairs
                pop[2, 0, i] = pop[1, 0, i]; pop[3, 1, i] = pop[2, 1, i]; pop[3, 0, i] = pop[1, 1, i]
                pop[5, 0, i] = pop[4, 0, i]; pop[6, 1, i] = pop[5, 1, i]; pop[6, 0, i] = pop[4, 1, i]
                pop[8, 0, i] = pop[7, 1, i]; pop[9, 0, i] = pop[8, 1, i]; pop[10, 0, i] = pop[9, 1, i]
                pop[11, 0, i] = pop[10, 1, i]; pop[12, 0, i] = pop[11, 1, i]; pop[13, 0, i] = pop[12, 1, i]
                pop[14, 0, i] = pop[13, 1, i]; pop[15, 0, i] = pop[14, 1, i]; pop[16, 0, i] = pop[15, 1, i]

        return fuv,cosmicrays,Z, tgas, tdust, abun, N, x, nh, itot, pop, tau


    # Av,obs - PDF function
    def pdf(x, s, m):
        return (1.0 / (s * math.sqrt(2.0 * math.pi))) * math.exp(-(math.log(x) - m) ** 2 / (2.0 * s ** 2))

    # Make Av-PDF
    def makepdf():
        m = math.log(av_bar) - s ** 2 / 2.0
        avtot = 500
        avin = [0] * (avtot + 1)
        pdfin = [0] * (avtot + 1)
        with open('avpdf.dat', 'w') as f:
            for i in range(avtot + 1):
                lav = -2.0 + 4.0 * float(i) / float(avtot + 1)
                avin[i] = 10.0 ** lav
                pdfin[i] = pdf(float(avin[i]), s, m)
                if pdfin[i] < 1e-10:
                    pdfin[i] = 1e-10
                f.write(str(avin[i]) + ' ' + str(pdfin[i]) + '\n')
        return avin,pdfin,avtot

    # Read inputted Av-PDF
    def readpdf():
        IERR = 0
        with open('avpdf_input.dat', 'r') as f:
            avtot = 0
            for line in f:
                avtot += 1
            f.seek(0)
            avin = [0] * avtot
            pdfin = [0] * avtot
            i = 0
            for line in f:
                values = line.split()
                avin[i], pdfin[i] = float(values[0]), float(values[1])
                if IERR != 0:
                    continue
                i += 1
            avin = [Decimal('10') ** Decimal(x) for x in avin]      # Issue: some values are too large.
            pdfin = [10.0 ** x for x in pdfin]
        return avin,pdfin,avtot

    # Make prefix for all inputs
    def makeprefix():
        num = [str(i).zfill(2) for i in range(61)]
        pref = []
        for i in range(41):
            for j in range(61):
                pref.append("zcr" + num[i] + "_uv" + num[j])
        return pref


    # Input parameters
    with open('./pdfchem.params', 'r') as f:
        av_bar, s, metallicity = [float(x) for x in f.readline().split()]

    if metallicity == 0.1:
        directory = 'Z0p1/'
    elif metallicity == 0.5:
        directory = 'Z0p5/'
    elif metallicity == 1.0:
        directory = 'Z1p0/'
    elif metallicity == 2.0:
        directory = 'Z2p0/'

    # Check if user has an input AV-PDF file. If not, call makepdf to create a simulated PDF.
    if os.path.isfile('avpdf_input.dat'):
        avin, pdfin,avtot = readpdf()
    else:
        avin, pdfin,avtot = makepdf()

    pref = makeprefix()

    outfile = open('output_py.dat', 'w')
    #for ipref in range(0, 1): # No. PDR simulations
    for ipref in range(0, 2501): # No. PDR simulations
        screen_time = datetime.datetime.now().strftime("\033[0;32m%H:%M:%S %Y-%m-%d > \33[0m")
        print(f"\r{screen_time} Looping through {ipref+1}/{2501} Simulations", end='', flush=True)
        prefix = pref[ipref]

        # calculate column densities of species
        N, Ntgas, Ntot_nopdf, Nrho = 0, 0, 0, 0
        fuv,cosmicrays,Z, tgas, tdust, abun, N, x, nh, itot, pop, tau = readfile(directory, prefix, Nspec) # call readfile function
        Tex = np.zeros((17, itot))
        Bnu = np.zeros((17, itot))
        for i in range(itot-1):
            Avobs = 0.06*(nh[i]**0.69)
            # for k in range(avtot):
            #     if avin[k] > Avobs: break
            k = bisect.bisect_left(avin, Avobs) - 1 # ???
            #if k >= avtot: k = avtot-1
            # print(f"{k}, {avin[k]}, {Avobs}\r")

            step = abs(x[i]-x[i+1])*pc2cm
            N[0] += 0.5*(nh[i]+nh[i+1])*step*pdfin[k] # total column density
            N[1:Nspec+1] += 0.5*(nh[i]*abun[i,:] + nh[i+1]*abun[i+1,:])*step*pdfin[k] # species
            Ntgas += 0.5*(nh[i]*tgas[i] + nh[i+1]*tgas[i+1])*step*pdfin[k] # for <Tgas>
            Nrho += 0.5*(nh[i]**2*abun[i,30] + nh[i+1]**2*abun[i+1,30])*step*pdfin[k]

            # Excitation temperatures and Background Radiation for radiative transfer

            for j in range(17):
                #if pop[j,1,i] == 0 or abs(g[j,1]*pop[j,0,i]/pop[j,1,i]/g[j,0]-1) < 1e-2:
                #    Tex[j,i] = 0
                #    Bnu[j,i] = 0
                if pop[j,1,i] != 0 or abs(g[j,1]*pop[j,0,i]/pop[j,1,i]/g[j,0]-1) >= 1e-2:
                #else:
                    try:
                        Tex[j,i] = (hp*freq0[j]/kb)/math.log(g[j,1]*pop[j,0,i]/pop[j,1,i]/g[j,0]) # excitation temperature
                    except:
                        Tex[j,i] = (hp*freq0[j]/kb)/np.log(g[j,1]*pop[j,0,i]/pop[j,1,i]/g[j,0]) # excitation temperature
                    Bnu[j,i] = (2.*hp*freq0[j]**3/c**2)/(math.exp(hp*freq0[j]/kb/Tex[j,i])-1) # Black body emission
            sigma = (freq0/c)*np.sqrt(kb*tgas[i]/mh+vturb**2/2.) # sigma value used in phi
            phi = 1/sigma/math.sqrt(2*pi) # phi value used in optical depth
            frac = 0.5*((pop[:,0,i]+pop[:,0,i+1])*g[:,1]/g[:,0]-(pop[:,1,i]+pop[:,1,i+1]))
            tau_incr += phi*(A*c**2/8./math.pi/freq0**2)*frac*step # optical depth calculation
            tau[:,i] = tau_incr # record value
        # solve radiative transfer equation
        Ntr, Ntot, t_r, Ncol = np.zeros(17), 0, np.zeros(17), 0
        tau_cii, tau_ci, tau_co = 0, 0, 0
        for i in range(itot-2):
            #Avobs = 0.06*(nh[i]**0.69) # Av,obs -- nH relation
            #for k in range(avtot):
            #    if (avin[k] > Avobs): break
            dtau = tau[:, i+1] - tau[:, i]
            for j in range(17): # frequencies explored (see subroutine constants)
                if (dtau[j] > 1e10):
                    t_r[j] = Bnu[j, i]
                elif (dtau[j] > 1e-6):
                    t_a_i = Bnu[j, i]*((1-math.exp(-dtau[j]))/dtau[j]-math.exp(-dtau[j]))
                    t_a_ip1 = Bnu[j, i+1]*(1-(1-math.exp(-dtau[j]))/dtau[j])
                    t_r[j] = t_r[j]*math.exp(-dtau[j])+t_a_i+t_a_ip1
                else:
                    t_r[j] = t_r[j]*(1-dtau[j])+(Bnu[j, i]+Bnu[j, i+1])*dtau[j]/2
            step = abs(x[i]-x[i+1])*pc2cm
            Ntot = Ntot + pdfin[k]
            Ntr = Ntr + (t_r*c**2/2/kb/freq0**2)*pdfin[k]
            Ncol = Ncol + (0.5*(nh[i]+nh[i+1])*step)
            tau_cii = tau_cii + dtau[0]*pdfin[k]
            tau_ci = tau_ci + dtau[1]*pdfin[k]
            tau_co = tau_co + dtau[7]*pdfin[k]
        #print("\n",  Ntr)

        # write output.dat
        outfile.write('{:11.2e}{:11.2e}{:11.2e}{:11.2e}'.format(fuv, cosmicrays, Z, Ntgas/N[0]))  #!parameters of PDR simulation
        for i in range(1,Nspec+1):
            outfile.write('{:11.2e}'.format(N[i]/N[0]))
        outfile.write('{:11.2e}{:11.2e}{:11.2e}'.format(Ntr[0], Ntr[1], Ntr[3]))
        for i in range(7, 17):
            outfile.write('{:11.2e}'.format(Ntr[i]))
        outfile.write('\n')
    outfile.close()

    print('\nFinished!')

if __name__ == '__main__':
    main()