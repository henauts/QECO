import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib.pyplot as plt
import time
from bacteria_numba import solve_model
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
q = float(sys.argv[2])
beta = 0.5

D_b = 1e-2
D_s = 5e-2
#D_s = 1e-1
t_c = 200
t_max = 2000000
dt_size = 200
x_max = 10
x_L = 5
dx = 0.1
n = 100

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

rhos = np.array(rhos)

print(f'Negative bacteria? {(rhos < 0).any()}\n')
if (rhos < 0).any():
    print(f'Terminating...')
    quit()

X, Y = np.meshgrid(x, [j*dt for j in range(idx)])
heat = plt.contourf(X, Y, rhos[1:], cmap = 'inferno', levels = 30)
plt.colorbar(heat)
plt.show()

plt.plot(tot_rho)
plt.show()

print('\nFinal interaction = ', idx)
print('\nFinal population = ', tot_rho[-1])
print('\nΔt = ', dt)
print('\nΔx = ', dx)

# np.savetxt(f'Results/Test_Densities_bacterial_model_v2_chi={chi}_q={beta}.txt', rhos)
# np.savetxt(f'Results/Test_Concentrations_bacterial_model_v2_chi={chi}_q={beta}.txt', Ss)
# np.savetxt(f'Results/Test_Total_pop_bacterial_model_v2_chi={chi}_q={beta}.txt', tot_rho)

print('\nDONE! :)')
