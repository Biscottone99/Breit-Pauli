#!/bin/bash
ifx basis.f90 -o basis.e
./basis.e
ifx geometria.f90 -o prova.e
./prova.e
python3 rotate.py
ulimit -s unlimited
ifx newmodule.f90 ppp-breit-pauli.f90 -o breit.e -qmkl -qopenmp -traceback 
./breit.e

cd unitary
ifx Unitary.f90 -o unitary.e -qmkl -qopenmp
./unitary.e <<EOF
70
rho.bin
eigen.bin
1
spin-density.bin
20000
1000
EOF
echo "Script completato."
