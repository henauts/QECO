import warnings
warnings.filterwarnings('ignore')
import numpy as np
import time
from bacteria_numba import solve_model
import sys
import os
import matplotlib.pyplot as plt

if os.path.exists('Results') == False:
    os.makedirs('Results')

print("\nWelcome to the ultra cool bacterial dynamics simulator \n")

r = 0.69 # h⁻¹
k = 0.1
chi = 315/3600 # μm²/h
gamma = r/k
lambd = 1.25 # mM/OD600*h
<<<<<<< HEAD
q = 2 # μM/h
=======
q = 1 # μM/h # M is a unit of concentration representing the number of millimols dissolved in 1 liter of solvent
>>>>>>> e0e52f762bd69cfbf8f31bdf00a97f79ea24c4ba
beta = 1

D_b = 50.2/3600 # μm²/h
D_s = 800/3600 # μm²/h
#D_s = 1e-1
t_c = 100
t_max = 1000000#*(t_c/128)
dt_size = 128
x_max = 10
x_L = 5
dx = 0.1
n = 100
S_plus = 30 # μM
S_minus = 1 # μM
S_max = 1 # μM

# Initial condition
S = np.zeros(n)
rho = np.random.uniform(0.08, 0.08, n)

print('\nReady to run simulations. Making a test before starting \n')

test = solve_model(5, rho, S, 0, 5, n, D_s, D_b, chi, r, k, lambd, 1, 5, 0.6, beta, S_plus, S_minus, S_max, 32)

print("\nTest concluded. Proceeding with simulation... \n")

# Initial condition
S = np.zeros(n)
rho = np.random.uniform(0.1, 0.1, n)

start_time = time.time()

rhos, Ss, tot_rho, tot_S, idx, dx, dt, x = solve_model(t_max, rho, S, 0,
                                                        x_max, n, D_s, D_b, chi,
                                                        r, k, lambd, t_c, x_L,
                                                        q, beta, S_plus, S_minus, S_max, dt_size)

end_time = time.time()

print(f'Time taken to run simulation: {(end_time-start_time):.2f} seconds \n')

plt.plot(tot_rho)
plt.axvline(t_c/dt, ls = '--')
# plt.text(0,0, tot_rho[int(np.round(t_c/dt, 0))])
# plt.text(0,0.1, int(np.round(t_c/dt, 0)))
# plt.text(0,0.2, t_c/dt)
# plt.xscale('log')
plt.show()

# print(len(rhos))
contour = plt.contourf(rhos, levels = 50, vmax = k, vmin = 0)
plt.colorbar(contour)
# plt.yscale('log')
plt.ylim(1e0, len(rhos))
plt.show()

print('\nFinal interaction = ', idx)
print('\nΔt = ', dt)
print('\nΔx = ', dx)

# np.savetxt(f'Results/Densities_bacterial_model_v2_chi={chi}_t_c={t_c}.txt', rhos)
# np.savetxt(f'Results/Concentrations_bacterial_model_v2_chi={chi}_t_c={t_c}.txt', Ss)
# np.savetxt(f'Results/Total_pop_bacterial_model_v2_chemotatic_ones_chi={chi}_alpha={alpha}.txt', tot_rho1)
# np.savetxt(f'Results/Total_pop_bacterial_model_v2_nonchemotatic_ones_chi={chi}_alpha={alpha}.txt', tot_rho2)
# np.savetxt(f'Results/delta_t_bacterial_model_v2_chi={chi}_alpha={alpha}.txt', np.array([dt]))

print('\nDONE! :)')
