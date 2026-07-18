# Time-dependent adiabatic GW

Introduction to time-dependent adiabatic GW (TD-aGW).

The goal is to calculate optical response of a material. It is obtained in three following ways.

## By GW-BSE

The calculation is done on BerkeleyGW for silicon. Codes are in the [`gw-bse-silicon/`](gw-bse-silicon) folder.

Main reference: Jack Deslippe et al. In: Computer Physics Communications 183.6 (2012), pp. 1269–1289. ISSN: 0010-4655.

## Within 3-band tight-binding model

This involves calculation of dynamical polarization for a 3-band tight-binding model. The codes are on https://github.com/skp-iisc/berry-phase-polarization/tree/main/tight-binding.

Main reference: Ivo Souza, Jorge Iniguez, and David Vanderbilt. In: Phys. Rev. B 69 (2004), p. 085106.

## Dynamical solution

The goal is to calculate optical response within TD-aGW. Here, the self-energy parts are not included. So, technically we can name it time-dependent Hartree method. The codes are in the [`dynamics-silicon/`](dynamics-silicon) folder.

Main reference: Yang-Hao Chan et al. In: Proceedings of the National Academy of Sciences 118.25 (2021), e1906938118.

All the results are saved in the [`results/`](results/) folder.

