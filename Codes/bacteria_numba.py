import warnings
warnings.filterwarnings('ignore')
import numpy as np
from numba import njit, prange

@njit
def solve_model_single_iteration(i, t_max, rho, S, x_min, x_max, n, D_s, D_b, chi, r, k, lambd, t_c, x_l, q, beta, S_plus, S_minus, S_max, dt_size):
    dx = (x_max - x_min) / n
    dt = dx**2 / (2 * dt_size * D_b)
    x = np.linspace(x_min, x_max, n)

    dirac = np.zeros(n)
    dirac[int(x_l / dx)] = 1

    if dt * i < t_c:
        q_eff = 0
    else:
        q_eff = q

    gamma = r/k

    S[n - 1] = S[n - 2]
    S[0] = S[1]
    rho[0] = rho[1]
    rho[n - 1] = rho[n - 2]

    dSdt = np.zeros(n)
    drhodt = np.zeros(n)

    # Vectorized computation of substance diffusion and degradation
    dSdt[1:n-1] = D_s * ((S[2:n] - S[1:n-1]) / dx**2 - (S[1:n-1] - S[:n-2]) / dx**2) - lambd * (S[1:n-1] / (S[1:n-1] + S_max)) * rho[1:n-1] + q_eff * dirac[1:n-1]

    # Vectorized computation of logarithm gradients for chemotaxis
    f_front = np.log((1 + S[2:n] / S_minus) / (1 + S[2:n] / S_plus))
    f_middle = np.log((1 + S[1:n-1] / S_minus) / (1 + S[1:n-1] / S_plus))
    f_back = np.log((1 + S[:n-2] / S_minus) / (1 + S[:n-2] / S_plus))

    # Vectorized computation of chemotaxis
    chem = chi * (((rho[2:n] - rho[:n-2]) * (f_front - f_back) / (4 * dx**2)) + rho[1:n-1] * ((f_front - 2 * f_middle + f_back) / dx**2))

    # Vectorized computation of growth, competition, and death by consumption
    growth = r * rho[1:n-1]
    competition = gamma * rho[1:n-1]**2
    death = lambd * beta * rho[1:n-1] * (S[1:n-1] / (S[1:n-1] + S_max))

    # Vectorized computation of bacterial equation
    drhodt[1:n-1] = D_b * ((rho[2:n] - 2 * rho[1:n-1] + rho[:n-2]) / dx**2) + chem + growth - competition - death
    
    # Update bacterial density and substance concentration
    rho[1:n-1] = rho[1:n-1] + drhodt[1:n-1] * dt
    S[1:n-1] = S[1:n-1] + dSdt[1:n-1] * dt

    return rho, S

@njit(parallel=True)
def solve_model_optimized(t_max, rho, S, x_min, x_max, n, D_s, D_b, chi, r, k, lambd, t_c, x_l, q, beta, S_plus, S_minus, S_max, dt_size, automatic_stop=True, print_progress=True):
    # Defining the step in space and time
    dx = (x_max - x_min)/n
    dt = dx**2 / (2*dt_size*D_b)

    x = np.linspace(x_min, x_max, n)

    # Array of solutions for density and concentration
    rhos = [rho]
    Ss = [S]
    tot_rho = [np.sum(rho*dx), np.sum(rho*dx)]
    tot_S = [np.sum(S*dx)]

    for i in prange(1, t_max):
        if print_progress and i % 1000000 == 0:
            print("Still computing... step:", i)

        if automatic_stop and np.abs(tot_rho[-1] - tot_rho[-2]) < 1e-9 and dt * i > 2 * t_c:
            break

        rho, S = solve_model_single_iteration(i, t_max, rho, S, x_min, x_max, n, D_s, D_b, chi, r, k, lambd, t_c, x_l, q, beta, S_plus, S_minus, S_max, dt_size)

        if i % 1000 == 0:
            tot_rho.append(np.sum(rho * dx))
            tot_S.append(np.sum(S * dx))

    return rhos, Ss, tot_rho, tot_S, i, dx, dt, x

# Solving the differential equation
@njit
def solve_model(t_max, rho, S, x_min, x_max, n, D_s, D_b, chi, r, k, lambd, t_c, x_l, q, beta, S_plus, S_minus, S_max, dt_size, save_every = 100, automatic_stop = True, print_progress = True):
    # Defining the step in space and time
    dx = (x_max - x_min)/n
    print(dx)
    dt = dx**2 / (2*dt_size*D_b)
    print(dt)
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
    death_by_competition = []
    death_by_substance = []
    
    #Dirac delta function
    dirac = np.zeros(n)
    # print(x_l)
    dirac[int(x_l/dx)] = 1
    
    gamma = r/k
    
    # Loop in time
    i = 0
    while i < t_max:
        # print(rho)
        if print_progress:
            if i % 1000000 == 0:
                print("Still computing... step:", i)
        if automatic_stop:
            if np.abs(tot_rho[-1] - tot_rho[-2]) < 1e-9 and dt*i > 2*t_c:
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


        # # Vectorized computation of substance diffusion and degradation
        # dSdt[1:n-1] = D_s * ((S[2:n] - S[1:n-1]) / dx**2 - (S[1:n-1] - S[:n-2]) / dx**2) - lambd * (S[1:n-1] / (S[1:n-1] + S_max)) * rho[1:n-1] + q_eff * dirac[1:n-1]

        # # Vectorized computation of logarithm gradients for chemotaxis
        # f_front = np.log((1 + S[2:n] / S_minus) / (1 + S[2:n] / S_plus))
        # f_middle = np.log((1 + S[1:n-1] / S_minus) / (1 + S[1:n-1] / S_plus))
        # f_back = np.log((1 + S[:n-2] / S_minus) / (1 + S[:n-2] / S_plus))

        # # Vectorized computation of chemotaxis
        # chem = chi * (((rho[2:n] - rho[:n-2]) * (f_front - f_back) / (4 * dx**2)) + rho[1:n-1] * ((f_front - 2 * f_middle + f_back) / dx**2))

        # # Vectorized computation of growth, competition, and death by consumption
        # growth = r * rho[1:n-1]
        # competition = gamma * rho[1:n-1]**2
        # death = lambd * beta * rho[1:n-1] * (S[1:n-1] / (S[1:n-1] + S_max))

        # # Vectorized computation of bacterial equation
        # drhodt[1:n-1] = D_b * ((rho[2:n] - 2 * rho[1:n-1] + rho[:n-2]) / dx**2) + chem + growth - competition - death

        # # Update bacterial density and substance concentration
        # rho[1:n-1] = rho[1:n-1] + drhodt[1:n-1] * dt
        # S[1:n-1] = S[1:n-1] + dSdt[1:n-1] * dt




        for j in range(1,n-1):
            compet = 0
            intake = 0
            # Substance diffusion and degradetion
            dSdt[j] = D_s*((S[j+1] - S[j])/dx**2 - (S[j]-S[j-1])/dx**2) - lambd*(S[j]/(S[j] + S_max))*rho[j] + q_eff*dirac[j]
            # Chemotaxis
            f_front = np.log((1 + S[j+1]/S_minus)/(1 + S[j+1]/S_plus)) # Logarithm gradient
            f_middle = np.log((1 + S[j]/S_minus)/(1 + S[j]/S_plus)) # Logarithm gradient
            f_back = np.log((1 + S[j-1]/S_minus)/(1 + S[j-1]/S_plus)) # Logarithm gradient
            chem = chi*(((rho[j+1] - rho[j-1])*(f_front - f_back)/(4*dx**2)) + rho[j]*((f_front - 2*f_middle + f_back)/dx**2))
            # Growth
            growth = r*rho[j]
            # Competition
            competition = gamma*rho[j]**2
            compet += competition
            # Death by consumption
            death = lambd*beta*rho[j]*(S[j]/(S[j] + S_max))
            intake += death
            # Bacterial equation
            drhodt[j] = D_b*((rho[j+1] - 2*rho[j] + rho[j-1])/dx**2) + chem + growth - competition - death
        #print(drhodt)
        rho = rho + drhodt*dt
        S = S + dSdt*dt

        #print(rhos)
        #print(rho)
        if i % save_every == 0:
            death_by_competition.append(compet)
            death_by_substance.append(intake)
            rhos.append(rho)
            Ss.append(S)
            tot_rho.append(np.sum(rho*dx, axis = 0))
            tot_S.append(np.sum(S*dx, axis = 0))
    return rhos, Ss, tot_rho, tot_S, death_by_competition, death_by_substance, i, dx, dt, x