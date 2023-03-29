import warnings
warnings.filterwarnings('ignore')
import numpy as np
from numba import njit

# Solving the differential equation
@njit
def solve_model(t_max, rho1, rho2, S, x_min, x_max, n, D_s, D_b, chi, r, k, lambd, t_c, t_f, x_l, q, beta, dt_size, alpha):
    # Defining the step in space and time
    dx = (x_max - x_min)/n
    dt = dx**2 / (2*dt_size*D_b)
    # print(dt)
    # print(t_c/dt)
    # Defining space
    x = np.linspace(x_min, x_max, n)
    # Array of derivatives
    drho1dt = np.zeros(n)
    drho2dt = np.zeros(n)
    dSdt = np.zeros(n)
    # Defining time
    #t = np.arange(0, t_max, dt)
    
    # Array of solutions for density and concentration
    rhos1 = [rho1]
    rhos2 = [rho2]
    Ss = [S]
    tot_rho1 = [np.sum(rho1*dx), np.sum(rho1*dx)]
    tot_rho2 = [np.sum(rho2*dx), np.sum(rho2*dx)]
    tot_S = [np.sum(S*dx)]
    
    #Dirac delta function
    dirac = np.zeros(n)
    dirac[int(x_l/dx)] = 1
    
    gamma = r/k
    
    # Loop in time
    i = 0
    while i < t_max:
        # print(rho1)
        if i % 50000 == 0:
             print("Still computing... step:", i)
        if np.abs(tot_rho1[-1] - tot_rho1[-2]) < 1e-9 and dt*i > 2*t_c and np.abs(tot_rho2[-1] - tot_rho2[-2]) < 1e-9:
            break
        i += 1
    #for i in range(len(t)):
        # Loop in space
        if dt*i < t_c or dt*i > t_f:
            q_eff = 0
        else:
            q_eff = q
        S[n-1] = S[n-2]
        S[0] = S[1]

        rho1[0] = rho1[1]
        rho1[n-1] = rho1[n-2]

        rho2[0] = rho2[1]
        rho2[n-1] = rho2[n-2]
        for j in range(1,n-1):
            # Substance diffusion and degradetion
            dSdt[j] = D_s*((S[j+1] - S[j])/dx**2 - (S[j]-S[j-1])/dx**2) - lambd*S[j]*rho1[j] + q_eff*dirac[j] - lambd*S[j]*rho2[j]
            # Chemotaxis
            chem = chi*(((rho1[j+1] - rho1[j-1])*(S[j+1] - S[j-1])/(4*dx**2)) + rho1[j]*((S[j+1] - 2*S[j] + S[j-1])/dx**2))
            # Growth
            growth1 = r*rho1[j]
            growth2 = r*rho2[j]
            # Competition
            competition1 = gamma*rho1[j]**2
            competition2 = gamma*rho2[j]**2

            resource_competition = r*alpha*rho1[j]*rho2[j]/k
            # Death by consumption
            death1 = lambd*beta*rho1[j]*S[j]
            death2 = lambd*beta*rho2[j]*S[j]
            # Bacterial equation
            drho1dt[j] = D_b*((rho1[j+1] - 2*rho1[j] + rho1[j-1])/dx**2) + chem + growth1 - competition1 - death1 - resource_competition
            drho2dt[j] = D_b*((rho2[j+1] - 2*rho2[j] + rho2[j-1])/dx**2) + growth2 - competition2 - death2 - resource_competition
        #print(drho1dt)
        rho1 = rho1 + drho1dt*dt
        rho2 = rho2 + drho2dt*dt
        S = S + dSdt*dt

        #print(rho1s)
        #print(rho1)
        rhos1.append(rho1)
        rhos2.append(rho2)
        Ss.append(S)
        tot_rho1.append(np.sum(rho1*dx, axis = 0))
        tot_rho2.append(np.sum(rho2*dx, axis = 0))
        tot_S.append(np.sum(S*dx, axis = 0))
    return rhos1, rhos2, Ss, tot_rho1, tot_rho2, tot_S, i, dx, dt, x