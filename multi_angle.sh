#!/usr/bin/env bash
set -Eeuo pipefail

#############################################
# LISTA DEGLI INPUT (.inp) DA PROCESSARE
#############################################
INPS=(
  "15.dat"
  "30.dat"
  "45.dat"
  "60.dat"
  "75.dat"
  "90.dat"
)

#############################################
# PARAMETRI per Unitary (modificabili)
#############################################
UNITARY_N=70
UNITARY_RHO="rho.bin"           # usa ../rho.bin se il file è nella dir superiore
UNITARY_EIGEN="eigen.bin"       # usa ../eigen.bin se serve
UNITARY_FLAG=1
UNITARY_SD="spin-density.bin"   # usa ../spin-density.bin se serve
UNITARY_PARAM_A=20000
UNITARY_PARAM_B=1000

#############################################
# COMPILAZIONE (una sola volta)
#############################################
echo "Compilo i programmi"

ifx basis.f90 -o basis.e
ifx newmodule.f90 ppp-breit-pauli.f90 -o breit.e -qmkl -qopenmp -traceback

pushd unitary >/dev/null
ifx Unitary.f90 -o unitary.e -qmkl -qopenmp
popd >/dev/null
echo "Compilazione completata."

#############################################
# LOOP PRINCIPALE
#############################################
ulimit -s unlimited

for INP in "${INPS[@]}"; do
  BASENAME="$(basename "$INP")"
  STEM="${BASENAME%.dat}"

  echo "=============================="
  echo "Elaboro: $INP  (stem: $STEM)"
  echo "=============================="

  # Copia l'input corrente nel nome atteso dai programmi
  cp -- "${STEM}_geom.dat" geom.dat

  echo "[${STEM}] Eseguo basis.e..."
  ./basis.e

  echo "[${STEM}] Eseguo rotate.py..."
  python3 rotate.py

  echo "[${STEM}] Eseguo breit.e..."
  ./breit.e

  # Salva output.out con nome specifico (se esiste)
  if [[ -f output.out ]]; then
    cp output.out "${STEM}.out"
  else
    echo "[${STEM}] Attenzione: output.out non trovato dopo breit.e" >&2
  fi

  echo "[${STEM}] Eseguo unitary.e..."
  pushd unitary >/dev/null
  ./unitary.e <<EOF
${UNITARY_N}
${UNITARY_RHO}
${UNITARY_EIGEN}
${UNITARY_FLAG}
${UNITARY_SD}
${UNITARY_PARAM_A}
${UNITARY_PARAM_B}
EOF

  # Salva prop1.dat con nome specifico (copia nella dir superiore)
  if [[ -f prop1.dat ]]; then
    cp prop1.dat "../${STEM}.dat"
  else
    echo "[${STEM}] Attenzione: prop1.dat non trovato in unitary/" >&2
  fi
  popd >/dev/null

  echo "[${STEM}] Completato."
done

echo "Tutti i job sono stati completati."
