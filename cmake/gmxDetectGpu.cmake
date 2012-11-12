# The gmx_detect_gpu() macro aims to detect GPUs available in the build machine.
# The current version is limited to checking the availability of NVIDIA GPUs
# and is limited in its ability to do so.
# It uses the following checks:
# - output of nvidia-smi (if available);
# - presence and content of of /proc/driver/nvidia/gpus/*/information (Linux)
# - presence of "/dev/nvidia[0-9]+" (Linux)
# - output of lspci (Linux)
#
# If any of the above checks is positive, meaning that the user does have an
# NVIDIA GPU, ????. 
#
# NOTE: The proper solution is to detect hardware compatible with the native
# GPU acceleration. However, this requires checking the compute capability
# of the device which is not possible with the current checks.
#

macro(gmx_detect_gpu)
    # nvidia-smi
    if (NOT DEFINED GMX_DETECT_GPU_NVIDIA_SMI_GPU_COUNT)
        set(GMX_DETECT_GPU_NVIDIA_SMI_GPU_COUNT 0)
        find_program(NVIDIA_SMI_PATH nvidia-smi)
        if (NVIDIA_SMI_PATH)
            exec_program(${NVIDIA_SMI_PATH}
                         ARGS -L
                         OUTPUT_VARIABLE _nvidia_smi_out
                         RETURN_VALUE    _nvidia_smi_ret)
            if (_nvidia_smi_ret EQUAL 0)
                string(REGEX REPLACE "\n" ";" _nvidia_smi_out "${_nvidia_smi_out}")
                foreach(_line ${_nvidia_smi_out})
                    if (_line MATCHES "^GPU [0-9]+:")
                        math(EXPR GMX_DETECT_GPU_NVIDIA_SMI_GPU_COUNT "${GMX_DETECT_GPU_NVIDIA_SMI_GPU_COUNT}+1") 
                    endif()
                endforeach()
            endif()
        endif()
    endif()

    if (UNIX AND NOT (APPLE OR CYGWIN))
        # /proc/driver/nvidia/gpus/*/information
        if (NOT DEFINED GMX_DETECT_GPU_PROC_GPU_COUNT)
            set(GMX_DETECT_GPU_PROC_GPU_COUNT 0)
            file(GLOB _proc_nv_gpu_info "/proc/driver/nvidia/gpus/*/information") 
            foreach (_file ${_proc_nv_gpu_info})
                math(EXPR GMX_DETECT_GPU_PROC_GPU_COUNT "${GMX_DETECT_GPU_PROC_GPU_COUNT}+1")
            endforeach()
        endif()

        # lspci
        if (NOT DEFINED GMX_DETECT_GPU_LSPCI_GPU_COUNT)
            set(GMX_DETECT_GPU_LSPCI_GPU_COUNT 0) 
            exec_program(lspci
                         OUTPUT_VARIABLE _lspci_out
                         RETURN_VALUE    _lspci_ret)
            if (_lspci_ret EQUAL 0)
                string(TOUPPER "${_lspci_out}" _lspci_out)
                STRING(REGEX REPLACE ";" "\\\\;" _lspci_out "${_lspci_out}")
                string(REGEX REPLACE "\n" ";" _lspci_out "${_lspci_out}")
                foreach(_line ${_lspci_out})
                    if (_line MATCHES ".*VGA.*NVIDIA.*" OR _line MATCHES ".*3D.*NVIDIA.*")
                        math(EXPR GMX_DETECT_GPU_LSPCI_GPU_COUNT "${GMX_DETECT_GPU_LSPCI_GPU_COUNT}+1") 
                    endif()
                endforeach()
            endif()
        endif()

        set(GMX_DETECT_GPU_DEV_GPU_COUNT 0)
    else()
        set(GMX_DETECT_GPU_PROC_GPU_COUNT 0)
        set(GMX_DETECT_GPU_DEV_GPU_COUNT 0)
        set(GMX_DETECT_GPU_LSPCI_GPU_COUNT 0)
    endif()
    
    message(">>> ${GMX_DETECT_GPU_NVIDIA_SMI_GPU_COUNT} <<<")
    message(">>> ${GMX_DETECT_GPU_PROC_GPU_COUNT} <<<")
    message(">>> ${GMX_DETECT_GPU_DEV_GPU_COUNT} <<<")
    message(">>> ${GMX_DETECT_GPU_LSPCI_GPU_COUNT} <<<")

    set(GMX_DETECT_GPU_NVIDIA_SMI_GPU_COUNT ${GMX_DETECT_GPU_NVIDIA_SMI_GPU_COUNT}
        CACHE STRING "Number of NVIDIA GPUs detected using nvidia-smi")
    set(GMX_DETECT_GPU_PROC_GPU_COUNT ${GMX_DETECT_GPU_PROC_GPU_COUNT}
        CACHE STRING "Number of NVIDIA GPUs detected in /proc/driver/nvidia/gpus")
    set(GMX_DETECT_GPU_DEV_GPU_COUNT ${GMX_DETECT_GPU_DEV_GPU_COUNT}
        CACHE STRING "Number of NVIDIA GPUs detected in /dev")
    set(GMX_DETECT_GPU_LSPCI_GPU_COUNT ${GMX_DETECT_GPU_LSPCI_GPU_COUNT}
        CACHE STRING "Number of NVIDIA GPUs detected using lspci")

    mark_as_advanced(GMX_DETECT_GPU_NVIDIA_SMI_GPU_COUNT
                     GMX_DETECT_GPU_PROC_GPU_COUNT
                     GMX_DETECT_GPU_DEV_GPU_COUNT
                     GMX_DETECT_GPU_LSPCI_GPU_COUNT)
endmacro(gmx_detect_gpu)
