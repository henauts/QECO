import subprocess
import time
import psutil
import numpy as np
import os
import matplotlib.pyplot as plt

print('\nRunning loop for checking chi and t_f!\n')

final_populations1 = np.zeros((50,40))
final_populations2 = np.zeros((50,40))
final_times = np.zeros((50,40))
# final_populations = np.zeros(200)
# final_times = np.zeros(200)
if os.path.exists('Results') == False:
    os.makedirs('Results')
np.savetxt('Results/final_populations_chi_and_t_f_substance_middle_chemotatic.txt', final_populations1)
np.savetxt('Results/final_populations_chi_and_t_f_substance_middle_nonchemotatic.txt', final_populations2)
np.savetxt('Results/final_times_chi_and_t_f_substance_middle.txt', final_times)
# np.savetxt('Results/final_populations_chi_substance_middle_q=1.txt', final_populations)

# start_time = time.time()
# for i in range(200):
#     chi = float(0.005*i)
#     if i <= 40:
#         end_time = time.time()
#         print(f'\nRun i = {i+1} of 200... excution time: {(end_time - start_time)/60:.2f} min... RAM usage: {psutil.virtual_memory()[3]/1000000000:.2f} GB ({psutil.virtual_memory()[2]:.2f}%)')
#         bashCommand = f"python bacterial_model_v2.py {chi}"
#         process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
#         output, error = process.communicate()
#         pop = np.loadtxt(f'Results/Total_pop_bacterial_model_v2_chi={chi}_q=1.txt')
#         final_pop = np.loadtxt('Results/final_populations_chi_substance_middle_q=1.txt')
#         final_pop[i] = pop[-1]
#         dt = np.loadtxt(f'Results/delta_t_bacterial_model_v2_chi={chi}_q=1.txt')
#         final_time = np.loadtxt('Results/final_times_chi_substance_middle_q=1.txt')
#         final_time[i] = len(pop)*dt
#         np.set_printoptions(precision=3)
#         print('\n', final_pop[:40])
#         print('\n', final_time[:40])
#         np.savetxt('Results/final_populations_chi_substance_middle_q=1.txt', final_pop)
#         np.savetxt('Results/final_times_chi_substance_middle_q=1.txt', final_time)
#         os.remove(f'Results/Total_pop_bacterial_model_v2_chi={chi}_q=1.txt')
#         os.remove(f'Results/delta_t_bacterial_model_v2_chi={chi}_q=1.txt')
#         print('\nFile updated, proceeding to next interation...')
#     else:
#         end_time = time.time()
#         print(f'Run i = {i+1} of 200... excution time: {(end_time - start_time)/60:.2f} min... RAM usage: {psutil.virtual_memory()[3]/1000000000:.2f} GB ({psutil.virtual_memory()[2]:.2f}%)', end = "\r")
#         bashCommand = f"python bacterial_model_v2.py {chi}"
#         process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
#         output, error = process.communicate()
#         pop = np.loadtxt(f'Results/Total_pop_bacterial_model_v2_chi={chi}_q=1.txt')
#         final_pop = np.loadtxt('Results/final_populations_chi_substance_middle_q=1.txt')
#         final_pop[i] = pop[-1]
#         dt = np.loadtxt(f'Results/delta_t_bacterial_model_v2_chi={chi}_q=1.txt')
#         final_time = np.loadtxt('Results/final_times_chi_substance_middle_q=1.txt')
#         final_time[i] = len(pop)*dt
#         np.savetxt('Results/final_populations_chi_substance_middle_q=1.txt', final_pop)
#         np.savetxt('Results/final_times_chi_substance_middle_q=1.txt', final_time)
#         os.remove(f'Results/Total_pop_bacterial_model_v2_chi={chi}_q=1.txt')
#         os.remove(f'Results/delta_t_bacterial_model_v2_chi={chi}_q=1.txt')

idx = 0
start_time = time.time()
for i in range(50):
    chi = float(0.005*i)
    for j in range(40):
        t_f = 220 + 20*j
        idx += 1
        end_time = time.time()
        if idx <= 40:
            print(f'\nRun in i = {i+1} of 50, j = {j+1} of 40... excution time: {(end_time - start_time)/60:.2f} min... RAM usage: {psutil.virtual_memory()[3]/1000000000:.2f} GB ({psutil.virtual_memory()[2]:.2f}%)')
        else:
            print(f'Run in i = {i+1} of 50, j = {j+1} of 40... excution time: {(end_time - start_time)/60:.2f} min... RAM usage: {psutil.virtual_memory()[3]/1000000000:.2f} GB ({psutil.virtual_memory()[2]:.2f}%)', end = '\r')
        bashCommand = f"python bacterial_model_v3.py {chi} {t_f}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        pop1 = np.loadtxt(f'Results/Total_pop_bacterial_model_v3_chemotatic_ones_chi={chi}_t_f={t_f}.txt')
        pop2 = np.loadtxt(f'Results/Total_pop_bacterial_model_v3_nonchemotatic_ones_chi={chi}_t_f={t_f}.txt')
        final_pop1 = np.loadtxt('Results/final_populations_chi_and_t_f_substance_middle_chemotatic.txt')
        final_pop2 = np.loadtxt('Results/final_populations_chi_and_t_f_substance_middle_nonchemotatic.txt')
        final_pop1[i,j] = pop1[-1]
        final_pop2[i,j] = pop2[-1]
        dt = np.loadtxt(f'Results/delta_t_bacterial_model_v3_chi={chi}_t_f={t_f}.txt')
        final_time = np.loadtxt('Results/final_times_chi_and_t_f_substance_middle.txt')
        final_time[i,j] = len(pop1)*dt
        if idx <= 80:
            np.set_printoptions(precision=3)
            print('\n', np.array(final_pop1[:2]) - np.array(final_pop2[:2]))
            print('\n', final_time[:2])
        np.savetxt('Results/final_times_chi_and_t_f_substance_middle.txt', final_time)
        np.savetxt('Results/final_populations_chi_and_t_f_substance_middle_chemotatic.txt', final_pop1)
        np.savetxt('Results/final_populations_chi_and_t_f_substance_middle_nonchemotatic.txt', final_pop2)
        os.remove(f'Results/Total_pop_bacterial_model_v3_chemotatic_ones_chi={chi}_t_f={t_f}.txt')
        os.remove(f'Results/Total_pop_bacterial_model_v3_nonchemotatic_ones_chi={chi}_t_f={t_f}.txt')
        os.remove(f'Results/delta_t_bacterial_model_v3_chi={chi}_t_f={t_f}.txt')
        #print('\nFile updated, proceeding to next interation...\n')

print('\nLoop finished :) trying preliminary plot of results\n')


data1 = np.loadtxt('Results/final_populations_chi_and_t_f_substance_middle_chemotatic.txt')
data2 = np.loadtxt('Results/final_populations_chi_and_t_f_substance_middle_nonchemotatic.txt')
chi = np.array([float(0.005*i) for i in range(50)])
t_f = np.array([float(200 + 20*j) for j in range(40)])

try:
    X, Y = np.meshgrid(t_f, chi)
    heat = plt.contourf(X, Y, data1 - data2, cmap = 'Spectral', levels = 40, vmax = 8, vmin = -8)
    plt.colorbar(heat)
    plt.xlabel('Substance Injection Stop')
    plt.ylabel('Chemotatic Response Strengh')
    plt.show()
except TypeError:
    X, Y = np.meshgrid(chi, t_f)
    heat = plt.contourf(X, Y, data1 - data2, cmap = 'Spectral', levels = 40, vmax = 8, vmin = -8)
    plt.colorbar(heat)
    plt.ylabel('Substance Injection Stop')
    plt.xlabel('Chemotatic Response Strengh')
    plt.show()
except:
    print("\nSorry, couldn't make the plot\n")

print('\nDone!')