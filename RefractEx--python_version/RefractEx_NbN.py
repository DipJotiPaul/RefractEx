import tmm, numpy as np, scipy.io, matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize

#%% Define the Drude-Lorentz Model Function
def drude_lorentz(wavelengths, params):    
    epsilon_inf, x2, x3, x4, x5, x6 = params    
    freq = 1.2398 / wavelengths                    # freq --> oscillator energy, wavelength in um 
    drude = x2**2 / (freq**2 + 1j * freq * x3)
    lorentz = x4**2 / (x5**2 - freq**2 - 1j * freq * x6)
    epsilon = epsilon_inf - drude + lorentz
    return np.sqrt(epsilon)

def tmm_function(params, wavelengths, nk_CaF2, thickd, th_0, spec_flag):
    nk_NbN = drude_lorentz(wavelengths, params)         # Calculate NbN refractive indices
    thickd[0] = thickd[0]                               # Define the list of layer thicknesses
    d_list = [np.inf] + thickd + [np.inf]
    th_0 = np.deg2rad(th_0)                             # Define incident angle (in radians)
    rtsim_s = []                      # Initialize arrays to store reflection and transmission values
    rtsim_p = []
    for k1 in range(len(wavelengths)):
        n_list = [1, nk_NbN[k1], nk_CaF2[k1], 1]        # Semi-infinite layers at the start and end
        # Calculate the reflection and transmission coefficients for s/p-polarization
        # rt_s = tmm.coh_tmm('s', n_list, d_list, th_0, wavelengths[k1]*1e-6)     # Ensure wavelength is in nm
        # rt_p = tmm.coh_tmm('p', n_list, d_list, th_0, wavelengths[k1]*1e-6) 
        c_list = ['i', 'c', 'i', 'i'] 
        rt_s = tmm.inc_tmm('s', n_list, d_list, c_list, th_0, wavelengths[k1]*1e-6)     # Ensure wavelength is in nm
        rt_p = tmm.inc_tmm('p', n_list, d_list, c_list, th_0, wavelengths[k1]*1e-6) 
        rtsim_s.append(np.abs(rt_s['R']) if spec_flag else np.abs(rt_s['T']))
        rtsim_p.append(np.abs(rt_p['R']) if spec_flag else np.abs(rt_p['T']))
    # Average reflections and transmissions for s- and p-polarizations to obtain unpolarized results
    rtsim = 0.5*(np.array(rtsim_s) + np.array(rtsim_p))
    return rtsim

def objective_function(params, wavelengths, thickd, incidence_angle, rtexp, nk_CaF2, spec_flag):
    rmse_list = []
    for angle, flag, rt_exp in zip(incidence_angle, spec_flag, rtexp):
        rtsim = tmm_function(params, wavelengths, nk_CaF2, thickd, angle, flag)
        rmse_wv = np.sqrt(np.mean((rt_exp - rtsim*1e2) ** 2))
        rmse_list.append(rmse_wv)
    print("RMSE: ",sum(rmse_list))
    return sum(rmse_list)

def initial_function(params, wavelengths, real_exp, imag_exp):
    nk_model = drude_lorentz(wavelengths, params)
    rmse_real = np.sum((real_exp - np.real(nk_model))**2)
    rmse_imag = np.sum((imag_exp - np.imag(nk_model))**2)
    rmse_total = np.sqrt((rmse_real + rmse_imag) / (2 * len(wavelengths)))
    print("RMSE:", rmse_total)
    return rmse_total
    
#%% Load FTIR data
mat_data = scipy.io.loadmat('FTIR_NbN_Spectra.mat')
wv = np.squeeze(mat_data['wv'])                     # Wavelength array in um
wavelengths = np.linspace(wv[0],wv[-1],300)         # wv[-1]
incidence_angle = [0,12,30,45,60]
spec_flag = [0, 1, 1, 1, 1]
thickd = [9.93e-9, 500e-6]                          # Thickness of NbN and CaF2 layers
CaF2 = mat_data['nk_CaF2']
nk_CaF2 = np.interp(wavelengths,CaF2[0].real,CaF2[2].real) + 1j*np.interp(wavelengths,CaF2[0].real,CaF2[2].imag)

NbN_T = np.interp(wavelengths,wv,np.squeeze(mat_data['NbN_Tran']))
NbN_R12 = np.interp(wavelengths,wv,np.squeeze(mat_data['NbN_R12']))
NbN_R30 = np.interp(wavelengths,wv,np.squeeze(mat_data['NbN_R30']))
# NbN_R35 = np.interp(wavelengths,wv,np.squeeze(mat_data['NbN_R35']))
# NbN_R40 = np.interp(wavelengths,wv,np.squeeze(mat_data['NbN_R40']))
NbN_R45 = np.interp(wavelengths,wv,np.squeeze(mat_data['NbN_R45']))
# NbN_R50 = np.interp(wavelengths,wv,np.squeeze(mat_data['NbN_R50']))
# NbN_R55 = np.interp(wavelengths,wv,np.squeeze(mat_data['NbN_R55']))
NbN_R60 = np.interp(wavelengths,wv,np.squeeze(mat_data['NbN_R60']))

rtexp = np.stack((NbN_T, NbN_R12, NbN_R30, NbN_R45, NbN_R60))

#%% Perform Nonlinear Optimization with Termination Criterion
x0_initial = [6.91, 6.72, 1.6, 17.85, 1.94, 19.99] 
bounds = [(1, 50), (0.1, 50), (0.1, 50), (0.1, 50), (0.1, 50), (0.1, 50)]

result = optimize.minimize(
    lambda params: objective_function(params, wavelengths, thickd, incidence_angle, rtexp, nk_CaF2, spec_flag),
    x0_initial, method='L-BFGS-B', bounds=bounds, tol=1e-10)

# Calculate RMSE using optimized parameters
optimized_params = result.x
RMSE = objective_function(optimized_params, wavelengths, thickd, incidence_angle, rtexp, nk_CaF2, spec_flag)
print(f"Optimized Parameters:\nepsilon_inf: {optimized_params[0]}\nx1: {optimized_params[1]}\nx2: {optimized_params[2]}\nx3: {optimized_params[3]}\nx4: {optimized_params[4]}\nx5: {optimized_params[5]}\nFinal RMSE: {RMSE}")

# Calculate NbN refractive indices
nk_NbN = drude_lorentz(wavelengths, optimized_params)
rtsim0 = tmm_function(optimized_params, wavelengths, nk_CaF2, thickd, incidence_angle[0], spec_flag[0])*1e2
rtsim1 = tmm_function(optimized_params, wavelengths, nk_CaF2, thickd, incidence_angle[1], spec_flag[1])*1e2
rtsim2 = tmm_function(optimized_params, wavelengths, nk_CaF2, thickd, incidence_angle[2], spec_flag[2])*1e2
rtsim3 = tmm_function(optimized_params, wavelengths, nk_CaF2, thickd, incidence_angle[3], spec_flag[3])*1e2
rtsim4 = tmm_function(optimized_params, wavelengths, nk_CaF2, thickd, incidence_angle[4], spec_flag[4])*1e2

# Comparing FTIR data with the simulated data
plt.figure(1)
plt.plot(wavelengths, rtsim0, 'b-', label='Fitted')
plt.plot(wavelengths, rtexp[0], 'r-', label='Exp')
plt.xlabel('wv(um)');   plt.ylabel('%T');   plt.legend();   plt.show()

plt.figure(2)
plt.plot(wavelengths, rtsim1, 'b-', label='Fitted')
plt.plot(wavelengths, rtexp[1], 'r-', label='Exp')
plt.xlabel('wv(um)');   plt.ylabel('%R12');   plt.legend();   plt.show()

plt.figure(3)
plt.plot(wavelengths, rtsim2, 'b-', label='Fitted')
plt.plot(wavelengths, rtexp[2], 'r-', label='Exp')
plt.xlabel('wv(um)');   plt.ylabel('%R30');   plt.legend();   plt.show()

plt.figure(4)
plt.plot(wavelengths, rtsim3, 'b-', label='Fitted')
plt.plot(wavelengths, rtexp[3], 'r-', label='Exp')
plt.xlabel('wv(um)');   plt.ylabel('%R45');   plt.legend();   plt.show()

plt.figure(5)
plt.plot(wavelengths, rtsim4, 'b-', label='Fitted')
plt.plot(wavelengths, rtexp[4], 'r-', label='Exp')
plt.xlabel('wv(um)');   plt.ylabel('%R60');   plt.legend();   plt.show()

plt.figure(6)
plt.plot(wavelengths, nk_NbN.real, 'b-', label='Fitted (Real)')
plt.plot(wavelengths, nk_NbN.imag, 'r-', label='Fitted (Imaginary)')
plt.xlabel('wv(um)');   plt.ylabel('%nk');   plt.legend();   plt.show()
