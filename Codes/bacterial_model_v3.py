import warnings
warnings.filterwarnings('ignore')
import numpy as np
import time
# from bacteria_numba_natural_selection_v2 import solve_model
from bacteria_numba import solve_model
import sys
import os
# import matplotlib.pyplot as plt

if os.path.exists('Results') == False:
    os.makedirs('Results')

print("\nWelcome to the ultra cool bacterial dynamics simulator \n")

r = 0.69
k = 0.1
chi = 315/3600
gamma = r/k
lambd = 1.25
q = float(sys.argv[1])
beta = 1
alpha = 1

D_b = 50.2/3600 # μm²/h
D_s = 800/3600 # μm²/h
#D_s = 1e-1
t_c = 100
t_f = 400
t_max = 1000000#*(t_c/128)
dt_size = 128
x_max = 10
x_L = float(sys.argv[2])
dx = 0.1
n = 100
S_plus = 30 # μM
S_minus = 1 # μM
S_max = 1 # μM

# Initial condition
S = np.zeros(n)
# rho1 = np.random.uniform(0.04, 0.04, n)
# rho2 = rho1.copy()
rho = np.random.uniform(0.04, 0.04, n)

print('\nReady to run simulations. Making a test before starting \n')

test = solve_model(5, rho, S, 0, 10, n, D_s, D_b, chi, r, k, lambd, 1, 5, 0.6,
                   beta, S_plus, S_minus, S_max, 32)

print("\nTest concluded. Proceeding with simulation... \n")

# Initial condition
S = np.zeros(n)
# rho1 = np.random.uniform(0.04, 0.04, n)
# rho2 = rho1.copy()
rho = np.random.uniform(0.04, 0.04, n)

start_time = time.time()

# rhos1, rhos2, Ss, tot_rho1, tot_rho2, tot_S, idx, dx, dt, x = solve_model(t_max, rho1, rho2, S, 0,
#                                                                            x_max, n, D_s, D_b, chi,
#                                                                            r, k, lambd, t_c, t_f, x_L,
#                                                                            q, beta, dt_size, alpha,
#                                                                            S_plus, S_minus, S_max)
rhos, Ss, tot_rho, tot_S, idx, dx, dt, x = solve_model(t_max, rho, S, 0,
                                                        x_max, n, D_s, D_b, chi,
                                                        r, k, lambd, t_c, x_L,
                                                        q, beta, S_plus, S_minus, S_max, dt_size)

end_time = time.time()

print(f'Time taken to run simulation: {(end_time-start_time):.2f} seconds \n')

print('\nFinal interaction = ', idx)
print('\nΔt = ', dt)
print('\nΔx = ', dx)

# plt.plot(tot_rho1)
# plt.plot(tot_rho2)
# # plt.plot(tot_S)
# plt.show()

# plt.plot(tot_rho)
# plt.plot(tot_S)
# plt.show()

# np.savetxt(f'Results/Densities_bacterial_model_v2_q={q}_x_L={x_L}.txt', rhos)
# np.savetxt(f'Results/Concentrations_bacterial_model_v2_q={q}_x_L={x_L}.txt', Ss)
# np.savetxt(f'Results/Total_pop_bacterial_model_v3_chemotatic_ones_chi={chi}_t_f={t_f}.txt', tot_rho1)
# np.savetxt(f'Results/Total_pop_bacterial_model_v3_nonchemotatic_ones_chi={chi}_t_f={t_f}.txt', tot_rho2)
np.savetxt(f'Results/Total_pop_bacterial_model_v2_q={q:.1f}_x_L={x_L:.1f}.txt', tot_rho)
np.savetxt(f'Results/delta_t_bacterial_model_v2_q={q:.1f}_x_L={x_L:.1f}.txt', np.array([dt]))

print('\nDONE! :)')
