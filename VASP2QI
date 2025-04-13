#!/bin/bash
file="${1:-.}"
num_processes="${2:-1}"
test -f "$file/WAVECAR" || echo "Error: WAVECAR not found!" >&2
test -f "$file/WAVECAR" || exit 1
test -f "$file/KPOINTS" || echo "Error: KPOINTS not found!" >&2
test -f "$file/KPOINTS" || exit 1
test -f "$file/OUTCAR" || echo "Error: OUTCAR not found!" >&2
test -f "$file/OUTCAR" || exit 1
num_kpoints="$(grep "Following reciprocal coordinates" $file/OUTCAR -B3|grep Found|cut -d " " -f 6)"
weights="$(grep "Following reciprocal coordinates" Test_Data/OUTCAR -A "$(( $num_kpoints + 1 ))"|tail -n $num_kpoints| tr "-" " "| tr -s " "| cut -d " " -f 5)"
./VASP2QI.py "$weights" "$num_kpoints" "$file" "1E16" "$num_processes"
