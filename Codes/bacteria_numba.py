import warnings
warnings.filterwarnings('ignore')
import numpy as np
from numba import njit

# Solving the differential equation
@njit
def solve_model(t_max, rho, S, x_min, x_max, n, D_s, D_b, chi, r, k, lambd, t_c, x_l, q, beta, dt_size):
    # Defining the step in space and time
    dx = (x_max - x_min)/n
    dt = dx**2 / (2*dt_size*D_b)
    # print(dt)
    # print(t_c/dt)
    # Defining space
    x = np.linspace(x_min, x_max, n)
    # Array of derivatives
    drhodt = np.zeros(n)
    dSdt = np.zeros(n)
    # Defining time
    #t = np.arange(0, t_max, dt)
    
    # Array of solutions for density and concentration
    rhos = [rho]
    Ss = [S]
    tot_rho = [np.sum(rho*dx), np.sum(rho*dx)]
    tot_S = [np.sum(S*dx)]
    
    #Dirac delta function
    dirac = np.zeros(n)
    dirac[int(x_l/dx)] = 1
    
    gamma = r/k
    
    # Loop in time
    i = 0
    while i < t_max:
        #print(rho)
        # if i % 20000 == 0:
        #      print("Still computing... step:", i)
        if np.abs(tot_rho[-1] - tot_rho[-2]) < 1e-8 and dt*i > 1.1*t_c:
            break
        i += 1
    #for i in range(len(t)):
        # Loop in space
        if dt*i < t_c:
            q_eff = 0
        else:
            q_eff = q
        S[n-1] = S[n-2]
        S[0] = S[1]
        rho[0] = rho[1]
        rho[n-1] = rho[n-2]
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
            drhodt[j] = D_b*((rho[j+1] - 2*rho[j] + rho[j-1])/dx**2) + chem + growth - competition - death
        #print(drhodt)
        rho = rho + drhodt*dt
        S = S + dSdt*dt

        #print(rhos)
        #print(rho)
        rhos.append(rho)
        Ss.append(S)
        tot_rho.append(np.sum(rho*dx, axis = 0))
        tot_S.append(np.sum(S*dx, axis = 0))
    return rhos, Ss, tot_rho, tot_S, i, dx, dt, x