import numpy as np 
import matplotlib.pyplot as plt 

Nprocs = 36

# Load for a given mode and check the b parameter: 

thisN = 0
thisT = 'S'
thisL = 2
thisTL1 = 2*thisL + 1

name = f"{thisN}{thisT}{thisL}"


Wmat = np.zeros(thisTL1)

for iproc in range(Nprocs): 


    Wmat += np.loadtxt(f"WmatProcs/Wmat_{name}_proc{iproc}.txt")[:thisTL1, :].diagonal()

print(Wmat)
print(Wmat[thisL] - Wmat[thisl-1])
