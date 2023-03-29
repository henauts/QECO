import warnings
warnings.filterwarnings('ignore')
import numpy as np
import time
from bacteria_numba_natural_selection import solve_model
import sys
import os

if os.path.exists('Results') == False:
    os.makedirs('Results')

print("\nWelcome to the ultra cool bacterial dynamics simulator \n")

r = 0.05
k = 0.8
chi = float(sys.argv[1])
gamma = r/k
lambd = 0.03
q = 0.6
beta = 0.5
alpha = 1

D_b = 1e-2
D_s = 5e-2
#D_s = 1e-1
t_c = 200
t_f = int(sys.argv[2])
t_max = 1000000#*(t_c/128)
dt_size = 128
x_max = 10
x_L = 5
dx = 0.1
n = 100

# Initial condition
S = np.zeros(n)
rho1 = np.random.uniform(0.05, 0.1, n)
rho2 = rho1.copy()

print('\nReady to run simulations. Making a test before starting \n')

test = solve_model(5, rho1, rho2, S, 0, 5, n, D_s, D_b, chi, r, k, lambd, 1, 2, 5, 0.6, beta, 32, alpha)

print("\nTest concluded. Proceeding with simulation... \n")

# Initial condition
S = np.zeros(n)
rho1 = np.random.uniform(0.05, 0.1, n)
rho2 = rho1.copy()

start_time = time.time()

rhos1, rhos2, Ss, tot_rho1, tot_rho2, tot_S, idx, dx, dt, x = solve_model(t_max, rho1, rho2, S, 0,
                                                                           x_max, n, D_s, D_b, chi,
                                                                           r, k, lambd, t_c, t_f, x_L,
                                                                           q, beta, dt_size, alpha)

end_time = time.time()

print(f'Time taken to run simulation: {(end_time-start_time):.2f} seconds \n')

print('\nFinal interaction = ', idx)
print('\nΔt = ', dt)
print('\nΔx = ', dx)

# np.savetxt(f'Results/Densities_bacterial_model_v2_chi={chi}_t_c={t_c}.txt', rhos)
# np.savetxt(f'Results/Concentrations_bacterial_model_v2_chi={chi}_t_c={t_c}.txt', Ss)
np.savetxt(f'Results/Total_pop_bacterial_model_v3_chemotatic_ones_chi={chi}_t_f={t_f}.txt', tot_rho1)
np.savetxt(f'Results/Total_pop_bacterial_model_v3_nonchemotatic_ones_chi={chi}_t_f={t_f}.txt', tot_rho2)
np.savetxt(f'Results/delta_t_bacterial_model_v3_chi={chi}_t_f={t_f}.txt', np.array([dt]))

print('\nDONE! :)')
