import numpy as np 
import matplotlib.pyplot as plt 
SCALE_T = 930.12663175873513

ylims = {'2S3':  [1.234, 1.244], '3S2': [1.092, 1.112],  '6S3': [2.810, 2.823],
         '8S5':  [4.151, 4.167], '9S3': [3.541, 3.556],  '11S4': [4.755, 4.776], 
         '11S5': [5.064, 5.082], '13S2': [4.837, 4.856], '13S3': [5.184, 5.204], 
         '14S4': [5.530, 5.557], '15S3': [6.020, 6.045], '16S6': [7.130, 7.170], 
         '18S4': [7.225, 7.250], '20S5': [8.440, 8.485], '21S6': [8.830, 8.860], 
         '23S5': [9.274, 9.305], '25S2': [9.010, 9.040], '27S2': [9.855, 9.885]
         }

modeNs = [2, 3, 6, 8, 9, 11, 11, 13, 13, 14, 15, 16, 18, 20, 21, 23, 25, 27]
modeLs = [3, 2, 3, 5, 3, 4, 5, 2, 3, 4, 3, 6, 4, 5, 6, 5, 2, 2]
nmodes = len(modeLs)

fig, axes = plt.subplots(6, 3, figsize=(6, 12))
fig.set_tight_layout(True)
rowctr = 0
colctr = 0
for imode in range(nmodes): 

    
    ax = axes[rowctr, colctr]
    ax.tick_params(labelsize=6)  # Adjust tick labels size
    ax.xaxis.label.set_size(6)  # Adjust x-axis label size
    ax.yaxis.label.set_size(6)  # Adjust y-axis label size




    n = modeNs[imode]
    t = 'S'
    l = modeLs[imode]
    tl1 = 2*l +1
    name = f"{n}{t}{l}"

    # Load frequency (stored in Hz)
    freq = np.loadtxt(f"../output/freqs/{name}.txt", skiprows=1) 
    omega0       = freq * 2*np.pi
    omega0nondim = omega0 * SCALE_T

    # Load the matrices (real parts):
    Hmat = np.loadtxt(f"matrices/RotEll_{name}_{name}.txt")[:tl1, :]

    # a b c parameters that define the quartic of m 
    abc = np.loadtxt(f"abcparams/abc{name}_{name}.txt")

    m     = np.arange(-l, l+1)
    mcont = np.linspace(-l, l, 100)

    dfcont = omega0*(abc[0] + abc[1]*mcont + abc[2]*mcont*mcont)/(2*np.pi)


    # Non dimensionalised Hmat 
    df = np.diag(Hmat)

    ax.axhline(freq*1000)
    ax.plot(mcont, (dfcont + freq)*1000, ':r')
    ax.plot(m, (df + freq)*1000, 'x')

    ax.set_ylim(ylims[name])
    ax.set_title(name, fontsize=6)

    

    colctr += 1
    if colctr==3: 
        colctr = 0
        rowctr += 1

fig.savefig(f'Tromp93_rot_ell.pdf', format='pdf')