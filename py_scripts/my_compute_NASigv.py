import numpy as np
from scipy.interpolate import interp1d

import argparse

# need mean thermal velocity in cm/s
# need cross section in cm^2 instead of barns

#Kb = 8.617*10**-5 #eV/Kelvin
Kb = 1.380649e-23 # J/K
N_A = 6.022e23 # avogadros constant

def Mean_thermal_Velocity(Temp, mu):
    # Temp - > Kelvin
    # mu -> reduced mass of the target nucleus and neutron.
    
    return np.sqrt(2.*Kb*Temp/mu)

def Maxwellian_Mean(Temp):
    # Temp - > Kelvin
    
    return 2./(np.sqrt(np.pi)*(Kb*Temp)**2)

def Maxwellian(x,Temp):
    #x -> eV
    #y -> barns   
    #Temp -> Kelvin
    
    return x*np.exp(-x/(Kb*Temp))


def parse_args():
    parser = argparse.ArgumentParser(description="This script reads En,sigma from input.dat and writes in output.out a Maxwellian-averaged-cross-section (MACS) given temperature. The maxwellian is also print out.")    
    parser.add_argument("-input",type = str, help = 'input file with incident energy [MeV] and cross section [barns]. Default = 99Tc_sig_cap_ENDF.txt', default = './99Tc_sig_cap_ENDF.txt')


    parser.add_argument('-A',type = int, help="used to get reduced mass of neutron + target, mu = m_n*m_A/(m_n+m_A). Default = 0.9985", default = -1)
    parser.add_argument('-n_points', type=int, help="Number of points for discretization. Default = 100000", default = 100000)

    # parser.add_argument('-T',type = float, default=348149007.7753278, help="temperature in kelvin. default = 348149007.7753278")


    #group = parser.add_mutually_exclusive_group()
    #group.add_argument('-T9s', type=str, help="Temperature in T9. Default = None", default = None)
    # group.add_argument('-KT', type=float, help="KT in eV. Default = 300",default = None)


    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    T9_list = [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]

    for j in T9_list:
        xfile = args.input

        # T given in T9, so convert temp to K
        temp = j * 10**9

        # get this in kg
        A = ""
        for i in range(len(xfile)):
            if xfile[i].isdigit():
                A += xfile[i]
        A = int(A)

        mu =  ((A * 1) / (A + 1)) * (1 / N_A) * (1 / 1000)

        DATA=np.loadtxt(xfile)

        Ens = DATA[:, 0]  # Incident energy in MeV
        Sigs = DATA[:, 1]  # Cross section in mbarns BECAUSE i am using TALYS

        # Convert MeV to J for the Maxwellian calculation
        Ens_J = Ens * (1e6) * (1.60218e-19)

        # convert mb to m^2 to maintain scripts function
        Sigs_m_squared = Sigs * (1 / 1000) * (1e-28)

        interp_Sigs = interp1d(Ens_J, Sigs_m_squared, kind='linear')

        # number of points for discretization
        n_points = args.n_points
        energies = []
        maxwellian_values = []

        min_energy = Ens_J.min()
        max_energy = Ens_J.max()
        de = (max_energy - min_energy)/n_points

        integral = 0.
        norm = 0.
        for i in range(0, n_points+1):
            energy = min_energy + i * de
            sig_eval = interp_Sigs(energy)

            maxwellian_value = Maxwellian(energy, temp)
            energies.append(energy)
            maxwellian_values.append(maxwellian_value)

            integral += sig_eval*maxwellian_value
            norm += maxwellian_value

        integral = integral*de
        #norm = norm*de
        norm=1
        


        # Calculate the Maxwellian-averaged cross section
        MACS =  (integral / norm) * Maxwellian_Mean(temp)
        v_T = Mean_thermal_Velocity(temp,mu)
        sigV = MACS*v_T
        
        # DEBUGGING
        # print(v_T, Maxwellian_Mean(temp), integral, norm)

        # fw = open('output.out','w')
        # print("temp[kelvin],MACS[mb],v_T,<sigv>",file = fw)
        # print(f"{temp},{MACS*1000.},{v_T},{sigV}",file = fw)
        # fw.close()

        # print("temp = ",temp, " kelvin")
        # print("MACS = ", MACS*1000., "mb")
        # print("v_T = ", v_T)
        # print("<sig v> = ", sigV)
        # print("N_A<sig v> = ", N_A*sigV)

        # want this in cm^3/(mol*s)
        # remember sigV is (m/s)*(m^2)
        print(
            str(temp/10**9) + "   " + str(N_A*sigV*(100**3))
        )

        # fw = open('Maxwellian.out','w')
        # print("# temp = %s kelvin"%temp, file = fw)
        # print("# E[eV] \t ", file = fw)
        # for i in range(len(energies)):
        #     print(" %s \t %s"%(energies[i],maxwellian_values[i]), file = fw)
        # fw.close()