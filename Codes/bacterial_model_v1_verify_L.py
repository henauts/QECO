import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm

print("Welcome to the ultra cool bacterial dynamics simulator \n")

# Solving the differential equation
def solve_model(t, rho, S, x_min, x_max, t_max, n, D_s, D_b, chi, r, k, lambd, t_c, x_l, q, beta):
    # Defining the step in space and time
    dx = (x_max - x_min)/n
    dt = dx**2 / (64*D_b)
    # Defining space
    x = np.linspace(x_min, x_max, n)
    # Array of derivatives
    drhodt = np.empty(n)
    dSdt = np.empty(n)
    # Defining time
    t = np.arange(0, t_max, dt)
    
    # Array of solutions for density and concentration
    rhos = [rho]
    Ss = [S]
    tot_rho = [np.sum(rho*dx)]
    tot_S = [np.sum(S*dx)]
    
    #Dirac delta function
    dirac = np.zeros(n)
    dirac[int(x_l/dx)] = 1
    
    gamma = r/k
    
    # Loop in time
    for i in range(len(t)):
        # Loop in space
        if t[i] < t_c:
            q_eff = 0
        else:
            q_eff = q
        S[n-1] = S[n-2]
        S[0] = S[1]
        rho[0] = D_b*rho[1]/(chi*(S[1] - S[0]) + D_b)
        rho[n-1] = D_b*rho[n-2]/(chi*(S[n-1] - S[n-2]) + D_b)
        for j in range(1,n-1):
            # Substance diffusion and degradetion
            dSdt[j] = D_s*((S[j+1] - S[j])/dx**2 - (S[j]-S[j-1])/dx**2) - lambd*S[j]*rho[j] + q_eff*dirac[j]
            # Chemotaxis
            chem = chi*(((rho[j+1] - rho[j-1])*(S[j+1] - S[j-1])/(4*dx**2)) + rho[j]*((S[j+1] - 2*S[j] + S[j-1])/dx**2))
            # Growth
            growth = r*rho[j]
            # Competition
            competition = gamma*rho[j]**2
            # Death by consumption
            death = lambd*beta*rho[j]*S[j]
            # Bacterial equation
            drhodt[j] = D_b*((rho[j+1] - rho[j])/dx**2 - (rho[j]-rho[j-1])/dx**2) + chem + growth - competition - death
        rho = rho + drhodt*dt
        S = S + dSdt*dt

        rhos.append(rho)
        Ss.append(S)
        tot_rho.append(np.sum(rho*dx))
        tot_S.append(np.sum(S*dx))
    return rhos, Ss, tot_rho, tot_S

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
t_final = 2200
t_c = 200

# Initial condition
S = np.zeros(n)
rho = np.random.uniform(0.05, 0.1, n)

x_max = np.arange(1, 20, 0.5)

tot_rho_final = []
n = 100
final_profile_rho = {}
final_profile_S = {}
tot_rho_profile = {}

print('\nReady to run simulations. Varying L \n')

for i in range(len(x_max)):
    dx = 0.1
    
    n = int(x_max[i]/dx)

    print(f"Run {i+1} of {len(x_max)} - Size: {x_max[i]} - Total steps: {n}", end = "\r")

    S = np.zeros(n)

    dt = dx**2 / (64*D_b)

    # Create the space
    x = np.arange(0, x_max[i], dx)

    # Create time array
    t = np.arange(0, t_final, dt)

    # S[50] = 5
    rho = np.random.uniform(0.05, 0.1, n)
    rhos, Ss, tot_rho, tot_S = solve_model(t, rho, S, 0, x_max[i], t_final, n, D_s, D_b, chi, r, k, lambd, t_c, x_max[i]/2, q, beta)
    rhos = np.array(rhos)
    Ss = np.array(Ss)
    if (rhos < 0).any() == False:
        tot_rho_final.append(tot_rho[-1])
        final_profile_rho[x_max[i]] = rhos[-1]
        final_profile_S[x_max[i]] = Ss[-1]
    else:
        print(f"There was an error. x_max = {x_max[i]:1f}")
        tot_rho_final.append(np.nan)
        final_profile_rho[x_max[i]] = np.nan
        final_profile_S[x_max[i]] = np.nan

print("Simulations done! \n")

try:
    pd.DataFrame(final_profile_rho).to_csv('final_rho_profile_varying_L.csv')
    pd.DataFrame(final_profile_S).to_csv('final_S_profile_varying_L.csv')
    pd.DataFrame(tot_rho_profile).to_csv('temporal_N_profile_varying_L.csv')
except:
    print('Could not save files as DataFrames, saving using Pickle instead')
    with open('final_rho_profile_varying_L.pkl', "wb") as fp:
        pickle.dump(final_profile_rho, fp)
    with open('final_S_profile_varying_L.pkl', "wb") as fp:
        pickle.dump(final_profile_S, fp)
    with open('temporal_N_profile_varying_L.pkl', "wb") as fp:
        pickle.dump(tot_rho_profile, fp)

np.savetxt('final_population_varying_L.txt', tot_rho_final)

print('DONE! :)')