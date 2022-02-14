#!/usr/bin/env bash
# Copyright 2022- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.

if [ "$#" -eq 0 ]; then
    echo "Usage:"
    echo ""
    echo "capture-topology.sh <output-root-path>"
    echo
    echo "CPU topology files will be copied from /sys/devices/system/cpu"
    echo "to <output-root-path>/sys/devices/system/cpu"
    echo ""
    exit 0
fi

declare -a files=(
    /sys/devices/system/cpu/possible
    /sys/devices/system/cpu/cpu*/topology/physical_package_id
    /sys/devices/system/cpu/cpu*/topology/core_id
)

mkdir -p "$1"

for file in "${files[@]}"; do
    cp -v --no-preserve=all --parents "$file" "$1"
done
