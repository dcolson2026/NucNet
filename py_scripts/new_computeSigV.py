import numpy as np
from scipy.interpolate import interp1d

import argparse



Kb = 8.617*10**-5 #eV/Kelvin

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


    parser.add_argument('-mu',type = float, default=0.9985, help="reduced mass of neutron + target, mu = m_n*m_A/(m_n+m_A). Default = 0.9985")
    parser.add_argument('-n_points', type=int, default=100000, help="Number of points for discretization. Default = 100000")

    # parser.add_argument('-T',type = float, default=348149007.7753278, help="temperature in kelvin. default = 348149007.7753278")


    group = parser.add_mutually_exclusive_group()
    group.add_argument('-T', type=float, help="Temperature in Kelvin. Default = 348149007.7753278", default=False)
    group.add_argument('-KT', type=float, help="KT in eV.",default=False)

    print(parser.parse_args())
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    xfile = args.input

    if args.KT is not None:
        temp = args.KT / Kb
    else:
        temp = args.T

    mu =  args.mu

    DATA=np.loadtxt(xfile)

    Ens = DATA[:, 0]  # Incident energy in MeV
    Sigs = DATA[:, 1]  # Cross section in barns

    # Convert MeV to eV for the Maxwellian calculation
    Ens_eV = Ens * 1e6

    interp_Sigs = interp1d(Ens_eV, Sigs, kind='linear')

    # number of points for discretization
    n_points = args.n_points
    energies = []
    maxwellian_values = []

    min_energy = Ens_eV.min()
    max_energy = Ens_eV.max()
    de = (max_energy -min_energy)/n_points

    itegral = 0.
    norm = 0.
    for i in range(0, n_points+1):
        energy = min_energy + i * de
        sig_eval = interp_Sigs(energy)

        maxwellian_value = Maxwellian(energy, temp)
        energies.append(energy)
        maxwellian_values.append(maxwellian_value)

        itegral = itegral + sig_eval*maxwellian_value
        norm = norm + maxwellian_value

    itegral = itegral*de
    norm = norm*de
    # 

    # Calculate the Maxwellian-averaged cross section
    MACS =  itegral / norm #* Maxwellian_Mean(temp)
    v_T = Mean_thermal_Velocity(temp,mu)
    sigV = MACS*Mean_thermal_Velocity(temp,mu)

                                                            
    fw = open('output.out','w')
    print("temp[kelvin],MACS[mb],v_T,<sigv>",file = fw)
    print(f"{temp},{MACS*1000.},{v_T},{sigV}",file = fw)
    fw.close()

    print("temp = ",temp, " kelvin")
    print("MCAS = ", MACS*1000., "mb")
    print("v_T = ", v_T)
    print("<sig v> = ", sigV)

    fw = open('Maxwellian.out','w')
    print("# temp = %s kelvin"%temp, file = fw)
    print("# E[eV] \t ", file = fw)

    for i in range(len(energies)):
        print(" %s \t %s"%(energies[i],maxwellian_values[i]), file = fw)

    fw.close()