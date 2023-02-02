import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import pickle
from bacteria_numba import solve_model

print("\nWelcome to the ultra cool bacterial dynamics simulator \n")

n = 100
center = 0
sd = 0.4
r = 0.05
k = 0.8
chi = 0.05
gamma = r/k
lambd = 0.03
q = 0.6
beta = 0.5

D_b = 1e-2
D_s = 5e-2
t_max = 1000000
t_c = 200
dt_size = 128

# Initial condition
S = np.zeros(n)
rho = np.random.uniform(0.05, 0.1, n)

x_max = np.arange(1, 25, 0.5)
q_s = np.arange(0.05, 1, 0.05)

tot_rho_final = np.zeros((len(x_max), len(q_s)))
n = 100
final_profile_rho = {}
final_profile_S = {}
tot_rho_profile = {}

test = solve_model(5, rho, S, 0, 10, n, D_s, D_b, chi, r, k, lambd, 1, 5, 0.6, beta, 68)

print('\nReady to run simulations. Varying L and q \n')

for i in range(len(x_max)):
    for j in range(len(q_s)):
        dx = 0.1
        
        n = int(x_max[i]/dx)

        print(f"Run {i+1} of {len(x_max)} - Size: {x_max[i]} - Total steps: {n}", end = "\r")

        S = np.zeros(n)

        dt = dx**2 / (64*D_b)

        # Create the space
        x = np.arange(0, x_max[i], dx)

        # S[50] = 5
        rho = np.random.uniform(0.05, 0.1, n)
        rhos, Ss, tot_rho, tot_S, idx, dx, dt, x = solve_model(t_max, rho, S, 0, x_max[i], n, D_s, D_b, chi, r, k, lambd, t_c, x_max[i]/2, q_s[j], beta, dt_size)
        rhos = np.array(rhos)
        Ss = np.array(Ss)
        if (rhos < 0).any() == False:
            tot_rho_final[i,j] = tot_rho[-1]
            final_profile_rho[f"{x_max[i]:.2f} - {q_s[j]:.2f}"] = rhos[-1]
            final_profile_S[f"{x_max[i]:.2f} - {q_s[j]:.2f}"] = Ss[-1]
            tot_rho_profile[f"{x_max[i]:.2f} - {q_s[j]:.2f}"] = tot_rho
        else:
            print(f"There was an error. x_max = {x_max[i]:1f}")
            tot_rho_final[i,j] = np.nan
            final_profile_rho[f"{x_max[i]:.2f} - {q_s[j]:.2f}"] = np.nan
            final_profile_S[f"{x_max[i]:.2f} - {q_s[j]:.2f}"] = np.nan
            tot_rho_profile[f"{x_max[i]:.2f} - {q_s[j]:.2f}"] = np.nan

print("\nSimulations done! \n")

try:
    pd.DataFrame(final_profile_rho).to_csv('final_rho_profile_varying_L-q.csv')
    pd.DataFrame(final_profile_S).to_csv('final_S_profile_varying_L-q.csv')
    pd.DataFrame(tot_rho_profile).to_csv('temporal_N_profile_varying_L-q.csv')
except:
    print('Could not save files as DataFrames, saving using Pickle instead')
    with open('final_rho_profile_varying_L-q.pkl', "wb") as fp:
        pickle.dump(final_profile_rho, fp)
    with open('final_S_profile_varying_L-q.pkl', "wb") as fp:
        pickle.dump(final_profile_S, fp)
    with open('temporal_N_profile_varying_L-q.pkl', "wb") as fp:
        pickle.dump(tot_rho_profile, fp)

np.savetxt('final_population_varying_L-q.txt', tot_rho_final)

print('DONE! :)')
