#!/bin/bash -l

cd ESPRESSO

for d in 0{0..8}*/; do
    # echo "Cleaning $d..."
    rm -f $d/Si.UPF
    rm -f $d/*.log
    rm -f $d/*out
    rm -rf $d/silicon.save
    rm -f $d/silicon.xml
    rm -f $d/*.real
    rm -f $d/*.dat
    rm -f $d/silicon_band*
    rm -f $d/*.amn
    rm -f $d/*.chk
    rm -f $d/*.eig
    rm -f $d/*.mmn
    rm -f $d/*.nnkp
done

cd ..

for d in *{7..12}*/; do
    # echo "Cleaning $d..."
    rm -f $d/*.log
    rm -f $d/*out
    rm -f $d/*.dat
    rm -f $d/RHO
    rm -f $d/WFN*
    rm -f $d/eps0mat*
    rm -f $d/epsmat*
    rm -f $d/bsemat*
    rm -f $d/bsedmat*
    rm -f $d/bsexmat*
    rm -f $d/dtmat
    rm -f $d/vmtxel
    rm -f $d/*.pdf
    rm -f $d/*.win
done
