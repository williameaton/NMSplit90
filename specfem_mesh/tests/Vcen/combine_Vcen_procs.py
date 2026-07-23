# Combines processors from different Vcen matrix outputs for SEM run 
import numpy as np 
import matplotlib.pyplot as plt 


nprocs = 6

N = 0
L = 2

tl1 = 2*L + 1 

name = f"{N}S{L}"

Vcen = np.zeros((tl1*2, tl1))


for iproc in range(nprocs): 
    Vcen += np.loadtxt(f"../rot_mat/VcenProcs/Vcen_{name}_{name}_proc{iproc}.txt")


# Kernel: 
V0 = np.loadtxt(f"vcen_woodhouse")[:tl1,:].diagonal()



Vcen = Vcen[:tl1,:].diagonal()

print(Vcen)
print()
print(V0)
print()
print(Vcen/V0)
