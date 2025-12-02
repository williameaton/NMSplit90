import numpy as np 
import matplotlib.pyplot as plt 
#SCALE_T = 930.12663175873513

ylims = {'2S3':  [1.234, 1.244], '3S2': [1.092, 1.112],  '6S3': [2.810, 2.823],
         '8S5':  [4.151, 4.167], '9S3': [3.541, 3.556],  '11S4': [4.755, 4.776], 
         '11S5': [5.064, 5.082], '13S2': [4.837, 4.856], '13S3': [5.184, 5.204], 
         '14S4': [5.530, 5.557], '15S3': [6.020, 6.045], '16S6': [7.130, 7.170], 
         '18S4': [7.225, 7.250], '20S5': [8.440, 8.485], '21S6': [8.830, 8.860], 
         '23S5': [9.274, 9.305], '25S2': [9.010, 9.040], '27S2': [9.855, 9.885]
         }

modeNs = [16]
modeLs = [5]
nmodes = len(modeLs)




fig, ax = plt.subplots(figsize=(10, 10))
fig.set_tight_layout(True)

#ax.tick_params(labelsize=6)  # Adjust tick labels size
#ax.xaxis.label.set_size(6)   # Adjust x-axis label size
#ax.yaxis.label.set_size(6)  # Adjust y-axis label size

imode = 0
n = modeNs[imode]
t = 'S'
l = modeLs[imode]
tl1 = 2*l +1
name = f"{n}{t}{l}"

# Load frequency (stored in Hz)
freq = np.loadtxt(f"../output/freqs/{name}.txt", skiprows=1) 
omega0       = freq * 2*np.pi
#omega0nondim = omega0 * SCALE_T

# Load the matrices (real parts):
# Note these are dimensionalised and in units of angular freq
Vcen = np.loadtxt(f"matrices/Vcen_{name}_{name}.txt")[:tl1, :]
Vell = np.loadtxt(f"matrices/Vell_{name}_{name}.txt")[:tl1, :]
Tell = np.loadtxt(f"matrices/Tell_{name}_{name}.txt")[:tl1, :]
Wmat = np.loadtxt(f"matrices/Wmat_{name}_{name}.txt")[:tl1, :]

# No longer need nondim omega0 since matrices are dimensionalised
#14.84
Hmat2 = Wmat + (Vell+Vcen - omega0*omega0*Tell)/(2.0*omega0)

Hmat = np.loadtxt(f"matrices/RotEll_{name}_{name}.txt")[:tl1, :]



# a b c parameters that define the quartic of m 
abc = np.loadtxt(f"abcparams/abc{name}_{name}.txt")

m     = np.arange(-l, l+1)
mcont = np.linspace(-l, l, 100)

dfcont = omega0*(abc[0] + abc[1]*mcont + abc[2]*mcont*mcont)/(2*np.pi)


# Hmat converted to Hz 
df  = np.diag(Hmat)/(2*np.pi)
df2 = np.diag(Hmat2)/(2*np.pi)

ax.axhline(freq*1000)
ax.plot(mcont, (dfcont + freq)*1000, '--r')
ax.plot(m, (df + freq)*1000, 'x')
ax.plot(m, (df2 + freq)*1000, 'o')

#ax.set_ylim(ylims[name])
ax.set_title(name, fontsize=6)



fig.savefig(f'{name}.pdf', format='pdf')