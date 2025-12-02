import numpy as np 


def test_tmat_ic(): 
    n = 0
    t = 'S'
    l = 2
    tag = f"{n}{t}{l}"

    ddir = './Tmat/'

    radial = np.loadtxt(f"{ddir}/radial_{tag}_{tag}.txt")
    sem    = np.loadtxt(f"{ddir}/Tmat_{tag}_{tag}.txt")


    (rows, cols) = np.shape(radial)
    for row in range(rows):
        for col in range(cols):

            # Sometimes true answer is zero: 
            if radial[row,col] == 0: 
                err = abs(sem[row,col])
                assert (err < 1e-10 )
            else: 
                err = 100*abs(sem[row,col]-radial[row,col])/radial[row,col] # percent
                assert (err < 2e-2 )
            print(err)


