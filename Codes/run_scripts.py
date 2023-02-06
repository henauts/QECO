import subprocess
import time
import psutil
import numpy as np
import os

print('\nRunning loop for checking L and D_s!\n')

final_populations = np.zeros((13,20))
final_times = np.zeros((13,20))
np.savetxt('Results/final_populations_L_and_D_s_substance_middle.txt', final_populations)
np.savetxt('Results/final_times_L_and_D_s_substance_middle.txt', final_times)

idx = 0
start_time = time.time()
for i in range(13):
    L = int(4 + 2*(i+1))
    for j in range(20):
        D_s = float(0.005*(j+1))
        idx += 1
        if idx <= 20:
            end_time = time.time()
            print(f'\nRun in i = {i+1} of 13, j = {j+1} of 20... excution time: {(end_time - start_time)/60:.2f} min... RAM usage: {psutil.virtual_memory()[3]/1000000000:.2f} GB ({psutil.virtual_memory()[2]:.2f}%)')
            bashCommand = f"python bacterial_model_v2.py {D_s} {L}"
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            pop = np.loadtxt(f'Results/Total_pop_bacterial_model_v2_L={L}_D_s={D_s}.txt')
            final_pop = np.loadtxt('Results/final_populations_L_and_D_s_substance_middle.txt')
            final_pop[i,j] = pop[-1]
            dt = np.loadtxt(f'Results/delta_t_bacterial_model_v2_L={L}_D_s={D_s}.txt')
            final_time = np.loadtxt('Results/final_times_L_and_D_s_substance_middle.txt')
            final_time[i,j] = len(pop)*dt
            np.set_printoptions(precision=3)
            print('\n', final_pop[:1])
            print('\n', final_time[:1])
            np.savetxt('Results/final_times_L_and_D_s_substance_middle.txt', final_time)
            np.savetxt('Results/final_populations_L_and_D_s_substance_middle.txt', final_pop)
            os.remove(f'Results/Total_pop_bacterial_model_v2_L={L}_D_s={D_s}.txt')
            print('\nFile updated, proceeding to next interation...')
        else:
            end_time = time.time()
            print(f'Run in i = {i+1} of 13, j = {j+1} of 20... excution time: {(end_time - start_time)/60:.2f} min... RAM usage: {psutil.virtual_memory()[3]/1000000000:.2f} GB ({psutil.virtual_memory()[2]:.2f}%)', end = "\r")
            bashCommand = f"python bacterial_model_v2.py {D_s} {L}"
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            pop = np.loadtxt(f'Results/Total_pop_bacterial_model_v2_L={L}_D_s={D_s}.txt')
            final_pop = np.loadtxt('Results/final_populations_L_and_D_s_substance_middle.txt')
            final_pop[i,j] = pop[-1]
            dt = np.loadtxt(f'Results/delta_t_bacterial_model_v2_L={L}_D_s={D_s}.txt')
            final_time = np.loadtxt('Results/final_times_L_and_D_s_substance_middle.txt')
            final_time[i,j] = len(pop)*dt[0]
            np.savetxt('Results/final_populations_L_and_D_s_substance_middle.txt', final_pop)
            os.remove(f'Results/Total_pop_bacterial_model_v2_L={L}_D_s={D_s}.txt')

print('\nDone!')