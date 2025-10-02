#!/bin/bash
ifx basis.f90 -o basis.e
./basis.e
ifx geometria.f90 -o prova.e
./prova.e
python3 rotate.py
ulimit -s unlimited
ifx newmodule.f90 ppp-breit-pauli.f90 -o vb.e -qmkl -qopenmp -traceback 
./vb.e

#cd redfield-ap
#./compile.sh
echo "Script completato."
