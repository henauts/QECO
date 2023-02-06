import warnings
warnings.filterwarnings('ignore')
import numpy as np
# import matplotlib.pyplot as plt
import time
from bacteria_numba import solve_model
import sys
import os

if os.path.exists('Results') == False:
    os.makedirs('Results')

print("\nWelcome to the ultra cool bacterial dynamics simulator \n")

r = 0.05
k = 0.8
chi = 0.05
gamma = r/k
lambd = 0.03
q = 0.6
beta = 0.5

D_b = 1e-2
D_s = float(sys.argv[1])
#D_s = 1e-1
t_c = 200
t_max = 1000000
dt_size = 64
x_max = int(sys.argv[2])
x_L = int(x_max/2)
dx = 0.1
n = int(x_max/dx)

# Initial condition
S = np.zeros(n)
rho = np.random.uniform(0.05, 0.1, n)

print('\nReady to run simulations. Making a test before starting \n')

test = solve_model(5, rho, S, 0, 5, n, D_s, D_b, chi, r, k, lambd, 1, 5, 0.6, beta, 32)

print("\nTest concluded. Proceeding with simulation... \n")

# Initial condition
S = np.zeros(n)
rho = np.random.uniform(0.05, 0.1, n)

start_time = time.time()

rhos, Ss, tot_rho, tot_S, idx, dx, dt, x = solve_model(t_max, rho, S, 0, x_max, n, D_s, D_b, chi, r, k, lambd, t_c, x_L, q, beta, dt_size)

end_time = time.time()

print(f'Time taken to run simulation: {(end_time-start_time):.2f} seconds \n')

# X, Y = np.meshgrid(x, [j*dt for j in range(idx)])
# heat = plt.contourf(X, Y, rhos[1:], cmap = 'inferno', levels = 30)
# plt.colorbar(heat)
# plt.show()

# print(tot_rho)
# print(rhos)
# plt.plot(tot_rho)
# plt.show()

print('\nFinal interaction = ', idx)
print('\nΔt = ', dt)
print('\nΔx = ', dx)

# np.savetxt(f'Results/Densities_bacterial_model_v2_chi={chi}_t_c={t_c}.txt', rhos)
# np.savetxt(f'Results/Concentrations_bacterial_model_v2_chi={chi}_t_c={t_c}.txt', Ss)
np.savetxt(f'Results/Total_pop_bacterial_model_v2_L={x_max}_D_s={D_s}.txt', tot_rho)
np.savetxt(f'Results/delta_t_bacterial_model_v2_L={x_max}_D_s={D_s}.txt', np.array([dt]))

print('\nDONE! :)')
