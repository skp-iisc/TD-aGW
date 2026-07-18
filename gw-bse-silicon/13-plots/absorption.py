import numpy as np
import matplotlib.pyplot as plt

# Read the data (lines starting with '#' are ignored)
data1 = np.loadtxt("../10-absorption/absorption_eh.dat", comments="#")
data2 = np.loadtxt("../10-absorption/absorption_noeh.dat", comments="#")

# Extract columns
omega_1 = data1[:, 0]
eps2_1  = data1[:, 1]
eps1_1  = data1[:, 2]
jdos_1  = data1[:, 3]

omega_2 = data2[:, 0]
eps2_2  = data2[:, 1]
eps1_2  = data2[:, 2]
jdos_2  = data2[:, 3]

# Plot: Imaginary Dielectric Function (epsilon2) vs omega

plt.figure(figsize=(6, 4))
# plt.axhline(0.0, color='k', linewidth=1)
# plt.axvline(0.0, color='k', linewidth=1)
plt.plot(omega_1, eps2_1, linewidth=1.5, label='with e-h interaction')
plt.plot(omega_2, eps2_2, linewidth=1.5, label='without e-h interaction')
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel(r'$\epsilon_2(\omega)$', fontsize=14)
plt.tick_params(direction='in', top=True, right=True)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("si_gw_bse_eps2.png", dpi=300)

# Plot: Real Dielectric Function (epsilon1) vs omega

plt.figure(figsize=(6, 4))
# plt.axhline(0.0, color='k', linewidth=1)
# plt.axvline(0.0, color='k', linewidth=1)
plt.plot(omega_1, eps1_1, linewidth=1.5, label='with e-h interaction')
plt.plot(omega_2, eps1_2, linewidth=1.5, label='without e-h interaction')
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel(r'$\epsilon_1(\omega)$', fontsize=14)
plt.tick_params(direction='in', top=True, right=True)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("si_gw_bse_eps1.png", dpi=300)

# Plot: Joint Density of States (JDOS) vs omega

plt.figure(figsize=(6, 4))
# plt.axhline(0.0, color='k', linewidth=1)
# plt.axvline(0.0, color='k', linewidth=1)
plt.plot(omega_1, jdos_1, linewidth=1.5, label='with e-h interaction')
plt.plot(omega_2, jdos_2, linewidth=1.5, label='without e-h interaction')
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel('JDOS', fontsize=14)
plt.tick_params(direction='in', top=True, right=True)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("si_gw_bse_jdos.png", dpi=300)
