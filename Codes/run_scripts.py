import subprocess
import time
import psutil
import numpy as np
import os
import matplotlib.pyplot as plt

print('\n\n\n\n\n\n\n')

print('Welcome to the ultra cool simulator. The code is about to run, to stop it at any moment press ctrl+C.')

print('\nRunning loop for checking L and D_s!\n')

D_ss = np.logspace(-2, 1, 250)
Ls = np.arange(2, 41, 2)
final_populations1 = np.zeros((len(D_ss),len(Ls)))
# final_populations2 = np.zeros((51,31))
final_times = np.zeros((len(D_ss),len(Ls)))
# final_populations = np.zeros(200)
# final_times = np.zeros(200)
if os.path.exists('Results') == False:
    os.makedirs('Results')
np.savetxt('Results/final_populations_D_s_and_L.txt', final_populations1)
# np.savetxt('Results/final_populations_chi_and_t_f_substance_middle_nonchemotatic.txt', final_populations2)
# np.savetxt('Results/populations_at_t_f_chi_and_t_f_substance_middle_chemotatic.txt', population1_tf)
# np.savetxt('Results/populations_at_t_f_chi_and_t_f_substance_middle_nonchemotatic.txt', population2_tf)
np.savetxt('Results/final_times_D_s_and_L.txt', final_times)
# np.savetxt('Results/final_populations_chi_substance_middle_D_s=1.txt', final_populations)

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

print('\nInitiating tests to estimate time...\n')
times_tests = []
for i in range(2):
    for j in range(2):
        start_test_time = time.time()
        bashCommand = f"python bacterial_model_v3.py {D_ss[(2*i+1)]} {Ls[2*(j+1)]}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        end_test_time = time.time()
        exec_test_time = (end_test_time - start_test_time)
        times_tests.append(exec_test_time)
print(f'\nOne execution takes {np.mean(times_tests):.2f} +/- {np.std(times_tests):.2f} seconds. {int(final_populations1.shape[0]*final_populations1.shape[1])} executions will take an estimated time of {int(final_populations1.shape[0]*final_populations1.shape[1])*np.mean(times_tests)/3600:.2f} +/- {int(final_populations1.shape[0]*final_populations1.shape[1])*np.std(times_tests)/3600:.2f} hours.\n')

print('\nStarting loop, relax and chill ;D\n')
idx = 0
start_time = time.time()
for i in range(len(D_ss)):
    D_s = D_ss[i]
    for j in range(len(Ls)):
        L = Ls[j]
        idx += 1
        end_time = time.time()
        # if idx <= 50:
        #     print(f'\nRun in i = {i+1} of 250, j = {j+1} of 30... Parameter values - chi: {chi}, t_f: {t_f}... excution time: {(end_time - start_time)/60:.2f} min... RAM usage: {psutil.virtual_memory()[3]/1000000000:.2f} GB ({psutil.virtual_memory()[2]:.2f}%)')
        # else:
        print(f'Run in i = {i+1} of {len(D_ss)}, j = {j+1} of {len(Ls)}... Parameter values - D_s: {D_s:.4f}, L: {L}... excution time: {(end_time - start_time)/3600:.2f} hours... RAM usage: {psutil.virtual_memory()[3]/1000000000:.2f} GB ({psutil.virtual_memory()[2]:.2f}%)', end = '\r')
        bashCommand = f"python bacterial_model_v3.py {D_s} {L}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        pop1 = np.loadtxt(f'Results/Total_pop_bacterial_model_v2_D_s={D_s:.4f}_L={L}.txt')
        # pop2 = np.loadtxt(f'Results/Total_pop_bacterial_model_v3_nonchemotatic_ones_q={q}_L={L}.txt')
        final_pop1 = np.loadtxt('Results/final_populations_D_s_and_L.txt')
        # final_pop2 = np.loadtxt('Results/final_populations_q_and_L_substance_middle_nonchemotatic.txt')
        # pop1_tf = np.loadtxt('Results/populations_at_L_q_and_L_substance_middle_chemotatic.txt')
        # pop2_tf = np.loadtxt('Results/populations_at_L_q_and_L_substance_middle_nonchemotatic.txt')
        final_pop1[i,j] = pop1[-1]
        # final_pop2[i,j] = pop2[-1]
        dt = np.loadtxt(f'Results/delta_t_bacterial_model_v2_D_s={D_s:.4f}_L={L}.txt')
        # pop1_tf[i,j] = pop1[int(np.round(L/dt))]
        # pop2_tf[i,j] = pop2[int(np.round(L/dt))]
        final_time = np.loadtxt('Results/final_times_D_s_and_L.txt')
        final_time[i,j] = len(pop1)*dt
        # if idx <= 100:
        #     np.set_printoptions(precision=3)
        #     print('\n', np.array(final_pop1[:2]) - np.array(final_pop2[:2]))
        #     print('\n', final_time[:2])
            # print(f'\n Values - chi: {chi}, L: {L}')
        np.savetxt('Results/final_times_D_s_and_L.txt', final_time)
        np.savetxt('Results/final_populations_D_s_and_L.txt', final_pop1)
        # np.savetxt('Results/final_populations_D_s_and_L_substance_middle_nonchemotatic.txt', final_pop2)
        # np.savetxt('Results/populations_at_L_D_s_and_L_substance_middle_chemotatic.txt', pop1_tf)
        # np.savetxt('Results/populations_at_L_D_s_and_L_substance_middle_nonchemotatic.txt', pop2_tf)
        os.remove(f'Results/Total_pop_bacterial_model_v2_D_s={D_s:.4f}_L={L}.txt')
        # os.remove(f'Results/Total_pop_bacterial_model_v3_nonchemotatic_ones_D_s={D_s}_L={L}.txt')
        os.remove(f'Results/delta_t_bacterial_model_v2_D_s={D_s:.4f}_L={L}.txt')
        #print('\nFile updated, proceeding to next interation...\n')

print('\nLoop finished :) trying preliminary plot of results\n')


data1 = np.loadtxt('Results/final_populations_D_s_and_L.txt')
# data2 = np.loadtxt('Results/final_populations_D_s_and_L_substance_middle_nonchemotatic.txt')

try:
    X, Y = np.meshgrid(Ls, D_ss)
    heat = plt.contourf(X, Y, data1, cmap = 'inferno', levels = 40, vmax = 1, vmin = 0)
    plt.colorbar(heat)
    plt.xlabel('Substance Injection Stop')
    plt.ylabel('Chemotatic Response Strengh')
    plt.show()
except TypeError:
    X, Y = np.meshgrid(D_ss, Ls)
    heat = plt.contourf(X, Y, data1, cmap = 'inferno', levels = 40, vmax = 1, vmin = 0)
    plt.colorbar(heat)
    plt.ylabel('Substance Injection Stop')
    plt.xlabel('Chemotatic Response Strengh')
    plt.show()
except:
    print("\nSorry, couldn't make the plot\n")

print('\nDone!')