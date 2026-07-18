import numpy as np
import matplotlib.pyplot as plt

# Read the data (lines starting with '#' are ignored)
data = np.loadtxt("../10-absorption/absorption_eh.dat", comments="#")
# data = np.loadtxt("../10-absorption/absorption_noeh.dat", comments="#")

# Extract columns
omega = data[:, 0]
eps2  = data[:, 1]
eps1  = data[:, 2]
jdos  = data[:, 3]

# Plot: Imaginary Dielectric Function (epsilon2) vs omega

plt.figure(figsize=(6, 4))
# plt.axhline(0.0, color='k', linewidth=1)
# plt.axvline(0.0, color='k', linewidth=1)
plt.plot(omega, eps2, linewidth=1.5)
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel(r'$\epsilon_2(\omega)$', fontsize=14)
plt.tick_params(direction='in', top=True, right=True)
plt.tight_layout()
plt.savefig("si_eps2_eh.png", dpi=300)
# plt.savefig("si_eps2_noeh.png", dpi=300)

# Plot: Real Dielectric Function (epsilon1) vs omega

plt.figure(figsize=(6, 4))
# plt.axhline(0.0, color='k', linewidth=1)
# plt.axvline(0.0, color='k', linewidth=1)
plt.plot(omega, eps1, linewidth=1.5)
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel(r'$\epsilon_1(\omega)$', fontsize=14)
plt.tick_params(direction='in', top=True, right=True)
plt.tight_layout()
plt.savefig("si_eps1_eh.png", dpi=300)
# plt.savefig("si_eps1_noeh.png", dpi=300)

# Plot: Joint Density of States (JDOS) vs omega

plt.figure(figsize=(6, 4))
# plt.axhline(0.0, color='k', linewidth=1)
# plt.axvline(0.0, color='k', linewidth=1)
plt.plot(omega, jdos, linewidth=1.5)
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel('JDOS', fontsize=14)
plt.tick_params(direction='in', top=True, right=True)
plt.tight_layout()
plt.savefig("si_jdos_eh.png", dpi=300)
# plt.savefig("si_jdos_noeh.png", dpi=300)
