import numpy as np

# see 12-inteqp/bandstructure.dat
# VBM and CBM at Gamma point
vbmMF, vbmQP = 6.006978020, 5.566072063
cbmMF, cbmQP = 8.567229331, 8.888389412

EgQP = cbmQP - vbmQP
EgMF = cbmMF - vbmMF

print(f"MF/DFT band-gap = {EgMF} eV.")
print(f"QP/QP band-gap = {EgQP} eV.")
