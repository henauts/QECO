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
t_max = 1000000
t_c = 200
dt_size = 128

# Initial condition
S = np.zeros(n)
rho = np.random.uniform(0.05, 0.1, n)

q_s = np.arange(0.05, 1, 0.05)
D_ss = np.linspace(5e-3, 1e-1, 15)

tot_rho_final = np.zeros((len(q_s), len(D_ss)))
n = 100
final_profile_rho = {}
final_profile_S = {}
tot_rho_profile = {}

test = solve_model(5, rho, S, 0, 10, n, 5e-2, D_b, chi, r, k, lambd, 1, 5, 0.6, beta, 32)

print('\nReady to run simulations. Varying q and D_s \n')

for i in range(len(q_s)):
    for j in range(len(D_ss)):
        dx = 0.1
        
        n = 100

        print(f"Run {i+1} of {len(q_s)} - q: {q_s[i]:.2f} - Total steps: {n}", end = "\r")

        S = np.zeros(n)

        dt = dx**2 / (2*dt_size*D_b)

        # Create the space
        x = np.arange(0, 10, dx)

        # S[50] = 5
        rho = np.random.uniform(0.05, 0.1, n)
        rhos, Ss, tot_rho, tot_S, idx, dx, dt, x = solve_model(t_max, rho, S, 0, 10, n, D_ss[j], D_b, chi, r, k, lambd, t_c, 5, q_s[i], beta, dt_size)
        rhos = np.array(rhos)
        Ss = np.array(Ss)
        if (rhos < 0).any() == False:
            tot_rho_final[i,j] = tot_rho[-1]
            final_profile_rho[f"{q_s[i]:.2f} - {D_ss[j]:.2f}"] = rhos[-1]
            final_profile_S[f"{q_s[i]:.2f} - {D_ss[j]:.2f}"] = Ss[-1]
            tot_rho_profile[f"{q_s[i]:.2f} - {D_ss[j]:.2f}"] = tot_rho
        else:
            print(f"There was an error. q_s = {q_s[i]:.1f} - D_s = {D_ss[j]:.3f}")
            tot_rho_final[i,j] = np.nan
            final_profile_rho[f"{q_s[i]:.2f} - {D_ss[j]:.2f}"] = np.nan
            final_profile_S[f"{q_s[i]:.2f} - {D_ss[j]:.2f}"] = np.nan
            tot_rho_profile[f"{q_s[i]:.2f} - {D_ss[j]:.2f}"] = np.nan

print("\nSimulations done! \n")

try:
    pd.DataFrame(final_profile_rho).to_csv('final_rho_profile_varying_q_s-D_s.csv')
    pd.DataFrame(final_profile_S).to_csv('final_S_profile_varying_q_s-D_s.csv')
    pd.DataFrame(tot_rho_profile).to_csv('temporal_N_profile_varying_q_s-D_s.csv')
except:
    print('Could not save files as DataFrames, saving using Pickle instead')
    with open('final_rho_profile_varying_q_s-D_s.pkl', "wb") as fp:
        pickle.dump(final_profile_rho, fp)
    with open('final_S_profile_varying_q_s-D_s.pkl', "wb") as fp:
        pickle.dump(final_profile_S, fp)
    with open('temporal_N_profile_varying_q_s-D_s.pkl', "wb") as fp:
        pickle.dump(tot_rho_profile, fp)

np.savetxt('final_population_varying_q_s-D_s.txt', tot_rho_final)

print('\nDONE! :)')
