#!/bin/bash -l

ln -sf ../../../pseudo/Si.UPF ./01-scf

for d in 0{2..8}*/; do
  mkdir -p $d/silicon.save
  ln -sf ../01-scf/Si.UPF $d/
  ln -sf ../../01-scf/silicon.save/data-file-schema.xml $d/silicon.save
  ln -sf ../../01-scf/silicon.save/charge-density.dat $d/silicon.save
done

ln -sf ../ESPRESSO/02-wfn/wfn.real ../07-epsilon/WFN
ln -sf ../ESPRESSO/03-wfnq/wfn.real ../07-epsilon/WFNq

ln -sf ../ESPRESSO/04-wfn_co/vxc.dat ../08-sigma/vxc.dat
ln -sf ../ESPRESSO/04-wfn_co/rho.real ../08-sigma/RHO
ln -sf ../ESPRESSO/04-wfn_co/wfn.real ../08-sigma/WFN_inner
ln -sf ../07-epsilon/eps0mat ../08-sigma
ln -sf ../07-epsilon/epsmat ../08-sigma
ln -sf ../07-epsilon/eps0mat.h5 ../08-sigma
ln -sf ../07-epsilon/epsmat.h5 ../08-sigma

ln -sf ../ESPRESSO/04-wfn_co/wfn.real ../09-kernel/WFN_co
ln -sf ../07-epsilon/eps0mat ../09-kernel
ln -sf ../07-epsilon/epsmat ../09-kernel
ln -sf ../07-epsilon/eps0mat.h5 ../09-kernel
ln -sf ../07-epsilon/epsmat.h5 ../09-kernel

ln -sf ../ESPRESSO/04-wfn_co/wfn.real ../10-absorption/WFN_co
ln -sf ../ESPRESSO/05-wfn_fi/wfn.real ../10-absorption/WFN_fi
ln -sf ../ESPRESSO/06-wfnq_fi/wfn.real ../10-absorption/WFNq_fi
ln -sf ../07-epsilon/eps0mat ../10-absorption
ln -sf ../07-epsilon/epsmat ../10-absorption
ln -sf ../09-kernel/bsedmat ../10-absorption
ln -sf ../09-kernel/bsexmat ../10-absorption
ln -sf ../07-epsilon/eps0mat.h5 ../10-absorption
ln -sf ../07-epsilon/epsmat.h5 ../10-absorption
ln -sf ../09-kernel/bsemat.h5 ../10-absorption

ln -sf ../ESPRESSO/04-wfn_co/wfn.real ../12-inteqp/WFN_co
ln -sf ../ESPRESSO/08-path/wfn.real ../12-inteqp/WFN_fi

