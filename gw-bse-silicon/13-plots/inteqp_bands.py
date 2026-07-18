import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# It reads the bandstructure.dat file from the inteqp.real.x code
data = np.loadtxt('../12-inteqp/bandstructure.dat')
bands = data[:,1]
kpts = data[:,2:5]
emf = data[:,5]
eqp = data[:,6]
emf -= np.amax(emf[bands==4])
eqp -= np.amax(eqp[bands==4])

def get_x(ks):
    # X axis is \int_0^k |dk|
    dk_vec = np.diff(ks, axis=0)
    dk_len = np.linalg.norm(dk_vec, axis=1)
    return np.insert(np.cumsum(dk_len), 0, 0.)

plt.figure(figsize=(6,4))

bands_uniq = np.unique(bands).astype(np.int32)
for ib in bands_uniq:
    cond = bands==ib
    x = get_x(kpts[cond])
    lqp, = plt.plot(x, eqp[cond], '-', color='orange', lw=1.5)
    lmf, = plt.plot(x, emf[cond], '--', color='blue', lw=1.5)

ind = [0, 45, 95, 145, 173, 199]
plt.xticks(x[ind], ['W', 'L', r'$\Gamma$', 'X', 'W', 'K'], fontsize=14)
for idx in ind:
    plt.axvline(x[idx], color='k', linestyle='--', linewidth=0.7)

plt.axhline(0, color='grey', linewidth=1)
plt.xlim(x[0], x[-1])
plt.ylim(-12.5, 7.5)
plt.legend([lmf, lqp], ['LDA', 'GW'], loc='lower right', fontsize=12)
plt.ylabel('Energy (eV)', fontsize=14)
plt.tight_layout()
plt.savefig('si_bands.png', dpi=300)
