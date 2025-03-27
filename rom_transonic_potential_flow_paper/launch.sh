#!/bin/bash
set -e  # Detener el script si ocurre un error

BASE_DIR=$(pwd)  # Guarda el directorio base donde se ejecuta el script

declare -a CASES=(
    "Rom_transonic_potential_flow_Naca0012/50"
    "Rom_transonic_potential_flow_Naca0012/75"
    "Rom_transonic_potential_flow_Naca0012/100"
    "Rom_transonic_potential_flow_Onera_M6/25"
    "Rom_transonic_potential_flow_Onera_M6/50"
    "Rom_transonic_potential_flow_Onera_M6/75"
)

for case in "${CASES[@]}"; do
    cd "$BASE_DIR/$case"
    python3 rom_manager.py
done
