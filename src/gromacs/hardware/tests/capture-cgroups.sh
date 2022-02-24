#!/usr/bin/env bash
# Copyright 2022- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.

if [ "$#" -eq 0 ]; then
    echo "Usage:"
    echo ""
    echo "capture-cgroups.sh <output-root-path>"
    echo
    echo "/etc/mtab, /proc/self/cgroup and files under the cgroup mount points"
    echo "specificed in /etc/mtab will be copied to <output-root-path>"
    echo ""
    exit 0
fi

# cgroup1 - note that cgroups can be mounted on multiple mount/subgroup points
cgroup1RootPaths=$(awk '/cgroup .*(cpu |cpu,cpuacct )/ {print $2}' /etc/mtab)
cgroup1SubPaths=$(awk -F ":" '/:(cpu|cpu,cpuacct):/ {print $3}' /proc/self/cgroup)
# Add empty subgroup
cgroup1SubPaths+=("/")

# cgroup2
cgroup2RootPaths=$(awk '/cgroup2/ {print $2}' /etc/mtab)
Cgroup2SubPath=$(awk -F ":" '/::/ {print $3}' /proc/self/cgroup)
# Add empty subgroup
cgroup2SubPaths+=("/")

mkdir -p "$1"

# common files
filesToCopy="/etc/mtab /proc/self/cgroup /proc/self/stat "

# Create the directory for all cgroups1 root/sub combinations we found,
# and copy cgroup.procs and cpu.cfs_period_us & cpu.cfs_quota_us if found in these.
for cgroup1RootPath in "${cgroup1RootPaths[@]}"; do
    for cgroup1SubPath in "${cgroup1SubPaths[@]}"; do
        path="${cgroup1RootPath}/${cgroup1SubPath}"
        [[ -d "${path}" ]] && mkdir -v -p "$1/${path}"
        [[ -f "${path}/cgroup.procs" ]] && filesToCopy+="${path}/cgroup.procs "
        [[ -f "${path}/cpu.cfs_period_us" ]] && filesToCopy+="${path}/cpu.cfs_period_us "
        [[ -f "${path}/cpu.cfs_quota_us" ]] && filesToCopy+="${path}/cpu.cfs_quota_us "
    done
done

# Create the directory for all cgroups2 root/sub combinations we found,
# and copy cgroup.procs and cpu.max if found in these.
for cgroup2RootPath in "${cgroup2RootPaths[@]}"; do
    for cgroup2SubPath in "${cgroup2SubPaths[@]}"; do
       	path="${cgroup2RootPath}/${cgroup2SubPath}"
       	[[ -d "${path}" ]] && mkdir -v -p "$1/${path}"
        [[ -f "${path}/cgroup.procs" ]] && filesToCopy+="${path}/cgroup.procs "
        [[ -f "${path}/cpu.max" ]] && filesToCopy+="${path}/cpu.max "
    done
done

# Copy all files in a single operation, since we need the PID of the current
# process (cp) to match between /proc/self/stat and cgroup.procs
cp -v --no-preserve=all --parents $filesToCopy $1
