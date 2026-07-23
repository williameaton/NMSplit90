# Combines the processor wise Vcens 
import numpy as np 


n = 0
l = 2

vcen = np.zeros((2*(2*l+1), 2*l+1))


for i in range(20): 
    vcen += np.loadtxt(f"VcenProcs/Vcen_{n}S{l}_{n}S{l}_proc{i}.txt")


SEM = vcen[:2*l+1, :2*l+1].diagonal()


analytical = np.loadtxt(f"../Vcen/vcen_woodhouse_{n}S{l}.txt")[:2*l+1, :2*l+1].diagonal()

print(SEM)
print(SEM/analytical)
