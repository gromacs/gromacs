#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2018, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

# The gmx_detect_gpu() macro aims to detect GPUs available in the build machine
# and provide the number, names, and compute-capabilities of these devices.
#
# The current version is limited to checking the availability of ROCM GPUs
# without compute-capability information.
#
# The current detection relies on the following checks in the order of listing:
# - output of rocm-smi (if available);
# - presence and content of of /proc/driver/nvidia/gpus/*/information (Linux)
# - output of lspci (Linux)
#
# If any of the checks succeeds in finding devices, consecutive checks will not
# be carried out. Additionally, when lspci is used and a device with unknown
# PCI ID is encountered, lspci tries to check the online PCI ID database. If
# this is not possible or the device is simply not recognized, no device names
# will be available.
#
# The following advanced variables are defined:
# - GMX_DETECT_GPU_AVAILABLE - TRUE if any GPUs were detected, otherwise FALSE
# - GMX_DETECT_GPU_COUNT     - # of GPUs detected
# - GMX_DETECT_GPU_INFO      - list of information strings of the detected GPUs
#
# NOTE: The proper solution is to detect hardware compatible with the native
# GPU acceleration. However, this requires checking the compute capability
# of the device which is not possible with the current checks and requires
# interfacing with the CUDA driver API.
#

# check whether the number of GPUs machetes the number of elements in the GPU info list
macro(check_num_gpu_info NGPU GPU_INFO)
    list(LENGTH ${GPU_INFO} _len)
    if (NOT NGPU EQUAL _len)
        list(APPEND ${GMX_DETECT_GPU_INFO} "NOTE: information about some GPU(s) missing!")
    endif()
endmacro()

macro(gmx_detect_gpu)

    if (NOT DEFINED GMX_DETECT_GPU_AVAILABLE OR
        NOT DEFINED GMX_DETECT_GPU_COUNT OR
        NOT DEFINED GMX_DETECT_GPU_INFO)

        set(GMX_DETECT_GPU_COUNT 0)
        set(GMX_DETECT_GPU_INFO  "")

        message(STATUS "Looking for ROCM GPUs present in the system")

        # rocm-smi-based detection.
        # Requires the rocm-smi tool to be installed and available in the path
        # or in one of the default search locations
        if (NOT DEFINED GMX_DETECT_GPU_COUNT_ROCM_SMI)
            # try to find the rocm-smi binary
            # TODO add location hints
            find_program(_rocm_smi "rocm-smi" PATHS /opt/rocm/bin/)
            if (_rocm_smi)
                set(GMX_DETECT_GPU_COUNT_ROCM_SMI 0)
                # execute rocm-smi -L to get a short list of GPUs available
                exec_program(${_rocm_smi_path} ARGS -i
                    OUTPUT_VARIABLE _rocm_smi_out
                    RETURN_VALUE    _rocm_smi_ret)
                # process the stdout of rocm-smi
                if (_rocm_smi_ret EQUAL 0)
                    # convert string with newlines to list of strings
                    string(REGEX REPLACE "\n" ";" _rocm_smi_out "${_rocm_smi_out}")
                    foreach(_line ${_rocm_smi_out})
                        if (_line MATCHES "^GPU\[[0-9]+\]")
                            math(EXPR GMX_DETECT_GPU_COUNT_ROCM_SMI "${GMX_DETECT_GPU_COUNT_ROCM_SMI}+1")
                        endif()
                    endforeach()

                    set(GMX_DETECT_GPU_COUNT ${GMX_DETECT_GPU_COUNT_ROCM_SMI})
                endif()
            endif()

            unset(_rocm_smi CACHE)
            unset(_rocm_smi_ret)
            unset(_rocm_smi_out)
            unset(_line)
        endif()

        if (UNIX AND NOT (APPLE OR CYGWIN))
            # lspci-based detection (does not provide GPU information).
            # Requires lspci and for GPU names to be fetched from the central
            # PCI ID db if not available locally.
            if (NOT DEFINED GMX_DETECT_GPU_COUNT_LSPCI AND GMX_DETECT_GPU_COUNT EQUAL 0)
                set(GMX_DETECT_GPU_COUNT_LSPCI 0)
                exec_program(lspci ARGS -q
                    OUTPUT_VARIABLE _lspci_out
                    RETURN_VALUE    _lspci_ret)
                # prehaps -q is not supported, try running without
                if (NOT RETURN_VALUE EQUAL 0)
                    exec_program(lspci
                        OUTPUT_VARIABLE _lspci_out
                        RETURN_VALUE    _lspci_ret)
                endif()
                if (_lspci_ret EQUAL 0)
                    # convert string with newlines to list of strings
                    STRING(REGEX REPLACE ";" "\\\\;" _lspci_out "${_lspci_out}")
                    string(REGEX REPLACE "\n" ";" _lspci_out "${_lspci_out}")
                    foreach(_line ${_lspci_out})
                        string(TOUPPER "${_line}" _line_upper)
                        if (_line_upper MATCHES ".*VGA.*AMD/ATI.*")
                            math(EXPR GMX_DETECT_GPU_COUNT_LSPCI "${GMX_DETECT_GPU_COUNT_LSPCI}+1")
                            # Try to parse out the device name which should be
                            # included in the lspci -q output between []-s
                            string(REGEX REPLACE ".*\\[(.*)\\].*$" "\\1" _gpu_name "${_line}")
                            if (NOT _gpu_name EQUAL "")
                                list(APPEND GMX_DETECT_GPU_INFO "${_gpu_name}")
                            endif()
                        endif()
                    endforeach()

                    check_num_gpu_info(${GMX_DETECT_GPU_COUNT_LSPCI} GMX_DETECT_GPU_INFO)
                    set(GMX_DETECT_GPU_COUNT ${GMX_DETECT_GPU_COUNT_LSPCI})
                endif()

                unset(_lspci_ret)
                unset(_lspci_out)
                unset(_gpu_name)
                unset(_line)
                unset(_line_upper)
            endif()
        endif()

        if (GMX_DETECT_GPU_COUNT GREATER 0)
            set(GMX_DETECT_GPU_AVAILABLE YES)
        else()
            set(GMX_DETECT_GPU_AVAILABLE NO)
        endif()
        set(GMX_DETECT_GPU_AVAILABLE ${GMX_DETECT_GPU_AVAILABLE} CACHE BOOL "Whether any ROCM GPU was detected" FORCE)

        set(GMX_DETECT_GPU_COUNT ${GMX_DETECT_GPU_COUNT}
            CACHE STRING "Number of ROCM GPUs detected")
        set(GMX_DETECT_GPU_INFO ${GMX_DETECT_GPU_INFO}
            CACHE STRING "basic information on the detected ROCM GPUs")

        set(GMX_DETECT_GPU_COUNT_ROCM_SMI ${GMX_DETECT_GPU_COUNT_ROCM_SMI}
            CACHE INTERNAL "Number of ROCM GPUs detected using rocm-smi")
        set(GMX_DETECT_GPU_COUNT_PROC ${GMX_DETECT_GPU_COUNT_PROC}
            CACHE INTERNAL "Number of ROCM GPUs detected in /proc/driver/nvidia/gpus")
        set(GMX_DETECT_GPU_COUNT_LSPCI ${GMX_DETECT_GPU_COUNT_LSPCI}
            CACHE INTERNAL "Number of ROCM GPUs detected using lspci")

        mark_as_advanced(GMX_DETECT_GPU_AVAILABLE
                         GMX_DETECT_GPU_COUNT
                         GMX_DETECT_GPU_INFO)

        if (GMX_DETECT_GPU_AVAILABLE)
            message(STATUS "Number of ROCM GPUs detected: ${GMX_DETECT_GPU_COUNT} ")
        else()
            message(STATUS "Could not detect ROCM GPUs")
        endif()

    endif()
endmacro(gmx_detect_gpu)
