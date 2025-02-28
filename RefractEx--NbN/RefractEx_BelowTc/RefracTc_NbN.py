import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.io

def epsr_integrand(E, omega, hbar, delta, kB, T, flag):
    if (kB*T) == 0:
        t1 = 0     
    else:
        exp_arg = np.clip((E+hbar*omega)/(kB*T), -700, 700)
        t1 = 1/(1 + np.exp(exp_arg))
    t2 = E**2 + delta**2 + hbar*omega*E
    if flag == 0:
        t3 = np.sqrt((delta**2 - E**2) * ((E + hbar*omega)**2 - delta**2))
    else:
        t3 = np.sqrt((E**2 - delta**2) * ((E + hbar*omega)**2 - delta**2))
    return (1-2*t1)*(t2/t3)

def epsi_integrand(E, omega, hbar, delta, kB, T):
    if (kB*T) == 0:
        t1 = 0   
        t2 = 0
    else:
        exp_arg1 = np.clip(E/(kB*T), -700, 700)
        t1 = 1/(1 + np.exp(exp_arg1))
        exp_arg2 = np.clip((E+hbar*omega)/(kB*T), -700, 700)
        t2 = 1/(1 + np.exp(exp_arg2))
    t3 = E**2 + delta**2 + hbar*omega*E
    t4 = np.sqrt((E**2 - delta**2) * ((E + hbar*omega)**2 - delta**2))
    return (t1-t2)*(t3/t4)

def epsr_integration(omegas, delta, hbar, eps0, sigmaN, kB, T):
    integral_values = []
    flag = 0
    for omega in omegas:
        if hbar*omega<2*delta:
            integral_value, _ = quad(epsr_integrand, delta-hbar*omega, delta, args=(omega, hbar, delta, kB, T, flag))
        else:
            integral_value, _ = quad(epsr_integrand, -delta, delta, args=(omega, hbar, delta, kB, T, flag))
        prefactor = sigmaN/(hbar*eps0*omega**2)
        integral_values.append(1-prefactor*integral_value)
    return integral_values

def epsi_integration(omegas, delta, hbar, eps0, sigmaN, kB, T):
    integral_values = []
    flag = 1
    for omega in omegas:
        integral_value1, _ = quad(epsi_integrand, delta, np.inf, args=(omega, hbar, delta, kB, T))
        integral_value2, _ = quad(epsr_integrand, delta-hbar*omega, -delta, args=(omega, hbar, delta, kB, T, flag))
        prefactor = sigmaN/(hbar*eps0*omega**2)
        integral_values.append(prefactor*(2*integral_value1+integral_value2))
    return integral_values

#%% NbN
eV = 1.602e-19 
hbar = 1.0546e-34
eps0 = 8.854e-12
sigmaN = 1/(2640.3e-9)                 # normal state conductivity
kB = 1.380649e-23                      # constant
T = 4                                  # NbN: Tc 8.9 K, T < Tc
Tc = 8.9
delta0 = 0.5*3.52*kB*Tc
delta = delta0*1.74*np.sqrt(1-(T/Tc))
Ereduce = np.linspace(0.1,1050, 9000)
omegas = Ereduce*2*delta/hbar
wv = 2*np.pi*3e8/omegas*1e6            # unit in um
omega_delta = 2*delta/hbar
# omega_plasma = 57.5e12
wv_delta = 2*np.pi*3e8/omega_delta*1e6            # unit in um
# wv_plasma = 2*np.pi*3e8/omega_plasma*1e6          # unit in um
epsr_result = epsr_integration(omegas, delta, hbar, eps0, sigmaN, kB, T)
epsi_result = epsi_integration(omegas, delta, hbar, eps0, sigmaN, kB, T)
nk_result = np.sqrt(np.abs(epsr_result)+1j*np.abs(epsi_result))

# Printing
# ind = [index for index, value in enumerate(epsr_result) if value > 0]
# print('Bandgap wv: '+str(wv_delta)+' um')
# # print('Plasma wv: '+str(wv_plasma)+' um')
# print('Zero crossing of epsr at Ereduced : '+str(Ereduce[ind[0]]))
# omega_zero = Ereduce[ind[0]]*2*delta/hbar
# wv_zero = 2*np.pi*3e8/omega_zero*1e6
# print('Zero crossing of epsr at wv : '+str(wv_zero)+' um')

#%% Plotting
plt.figure(1)
plt.plot(wv, np.real(nk_result), color='blue', label='BCS theory, @ 4K')
plt.ylabel('Refractive index, n')
plt.xlim(0, 25)
plt.ylim(1, 20)
plt.xlabel('wv [um]')
plt.grid(True)

plt.figure(2)
plt.plot(wv, np.imag(nk_result), 'red', label='BCS theory, @ 4K')
plt.ylabel('Extinction coefficient, k')
plt.xlim(0, 25)
plt.ylim(1, 20)
plt.xlabel('wv [um]')
plt.grid(True)
scipy.io.savemat('NbN_nk_4K.mat', {'wv': wv, 'nk_result': nk_result})

#%% Matlab compare
x = np.array([6.906980816237389, 6.716385517669770, 1.599998412858263, 
              17.847903843066636, 1.935406412562817, 19.999981207736077])
wv = np.linspace(0.1, 30, 5000)
omega = 1.2398 / wv
Drude = x[1]**2 / (omega**2 + 1j * omega * x[2])
Lorentz = x[3]**2 / (x[4]**2 - omega**2 - 1j * omega * x[5])
ncal = np.sqrt(x[0] - Drude + Lorentz)

plt.figure(1)
plt.plot(wv, np.real(ncal), color='blue', linestyle='dashed', label='FTIR, @ 290K')
plt.legend()

plt.figure(2)
plt.plot(wv, np.imag(ncal), color='red', linestyle='dashed', label='FTIR, @ 290K')
plt.legend()

