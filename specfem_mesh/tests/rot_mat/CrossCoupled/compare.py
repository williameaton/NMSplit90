import numpy as np 
import matplotlib.pyplot as plt 

n1 = 2
t1 = 'S'

n2 = 3
t2 = 'T'
l2 = 2
name2 = f"{n2}{t2}{l2}"
tl2 = 2*l2 + 1


fig, ax = plt.subplots(2,2)


for ibench in range(2): 
    if ibench == 0: 
        l1 = 3
    else: 
        l1 = 4

    tl1 = 2*l1 + 1

    name1 = f"{n1}{t1}{l1}"
    

    # load Radial eigens: 
    # Note it is called ICOnly but its actually whole Earth: 
    # It will be imaginary only and diagonal but not on the true diag 
    re = np.loadtxt(f"IConly_{name1}_{name2}.txt")[tl1:, :]

    # Load the 36 procs from the SEM: 
    sem = np.zeros_like(re)
    for iproc in range(36):
        sem += np.loadtxt(f"CCwhole/CC_Wmat_{name1}_{name2}_proc{iproc}.txt")[tl1:, :]


    print(np.min(sem), np.max(sem))

    ax[0,ibench].imshow(re)
    ax[1,ibench].imshow(sem)


fig.savefig("cross_benchmark.pdf", format='pdf')