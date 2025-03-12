import numpy as np 
import matplotlib.pyplot as plt 
SCALE_T = 930.12663175873513

ylims = {'3S2': [1.092, 1.112]}


n = 3
t = 'S'
l = 2
tl1 = 2*l +1
name = f"{n}{t}{l}"

# Load frequency (stored in Hz)
freq = np.loadtxt(f"../output/freqs/{name}.txt", skiprows=1) 
omega0       = freq * 2*np.pi
omega0nondim = omega0 * SCALE_T

# Load the matrices (real parts):
Wmat = np.loadtxt(f"matrices/Wmat_{name}_{name}.txt")[:tl1, :]
Vcen = np.loadtxt(f"matrices/Vcen_{name}_{name}.txt")[:tl1, :]
Vell = np.loadtxt(f"matrices/Vell_{name}_{name}.txt")[:tl1, :]
Tell = np.loadtxt(f"matrices/Tell_{name}_{name}.txt")[:tl1, :]

print(Wmat.min(), Wmat.max())
print(Vcen.min(), Vcen.max())
print(Vell.min(), Vell.max())
print(Tell.min(), Tell.max())

# Non dimensionalised Hmat 
Hmat = Wmat # + (Vell + Vcen - omega0*omega0*Tell)/(2*omega0)

Hmatdim = Hmat/SCALE_T      # in radians?

HmatdimHz = Hmatdim/(2*np.pi)

df = np.diag(HmatdimHz)

fig, ax = plt.subplots()

ax.axhline(freq*1000)
ax.plot(np.arange(-l, l+1), (df + freq)*1000)

ax.set_ylim(ylims[name])


fig.savefig(f'{name}.pdf', format='pdf')