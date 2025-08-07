#!/usr/bin/env bash
set -e
cd regressiontests
# Needed to run MPI enabled code in the docker images, until we set up different users
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export ASAN_OPTIONS="check_initialization_order=1:detect_invalid_pointer_pairs=1:strict_init_order=true:strict_string_checks=true:detect_stack_use_after_return=true"
export LSAN_OPTIONS="suppressions=$CI_PROJECT_DIR/admin/lsan-suppressions.txt:print_suppressions=0"

if grep -qF 'NVIDIA' <<< "$GPU_VENDOR"; then
    nvidia-smi -L && nvidia-smi || true;
    if [[ "$GMX_ENABLE_NVSHMEM" != "" ]] && [[ "$GPU_COUNT" -eq "2" ]]; then
        # if NVSHMEM uses NCCL it requires more than 4 GB of GPU RAM as min memory,
        # in CI we've GPUs like T400 which have 4 GB RAM so we disable NVSHMEM's
        # NCCL usage, as NCCL isn't required in our NVSHMEM usage.
        export NVSHMEM_DISABLE_NCCL=1
    fi
    if [[ "$GMX_GPU_PME_DECOMPOSITION" != "" ]]; then
        # cuFFTMp uses NVSHMEM, which in turn will use NCCL on single-node when available
        # But NCCL has extra restrictions (such as only 1 rank per GPU) so we disable for CI
        export NVSHMEM_DISABLE_NCCL=1
        # Default NVSHMEM_SYMMETRIC_SIZE is around 1.5GB which causes mem allocation error on CI
        # GPUs having 4GB vidmem also this size is overkill for PME decomp tests, hence we
        # restrict the default NVSHMEM heap size to 32MB.
        export NVSHMEM_SYMMETRIC_SIZE=32M
        # set CUDA Module load to be lazy in order to save on CUDA libraries memory footprint
        # as without lazy loading whole cuFFTMp library may get loaded into GPU RAM.
        export CUDA_MODULE_LOADING=LAZY
        # set nvshmem cumem allocation granularity to 2MB, default is 512MB.
        export NVSHMEM_CUMEM_GRANULARITY=2097152
    fi
fi
if grep -qF 'AMD' <<< "$GPU_VENDOR"; then
    clinfo -l || true;
    rocm-smi || true;
    if [[ "$GPU_COUNT" -eq "2" ]]; then
        # Our CI has some issues with ROCm (and CUDA) IPC; disable rocm_ipc
        export UCX_TLS="sm,self,rocm_copy"
    fi
    export HSA_ENABLE_SDMA=0  # Work around CI issues, Issue #5341
fi
if grep -qF 'INTEL' <<< "$GPU_VENDOR"; then
    sycl-ls || true;
    export SYCL_CACHE_PERSISTENT=1  # Issue #4218
fi

perl gmxtest.pl $REGRESSIONTEST_PARALLEL $REGRESSIONTEST_TOTAL_RANK_NUMBER -ntomp $REGRESSIONTEST_OMP_RANK_NUMBER -npme $REGRESSIONTEST_PME_RANK_NUMBER $REGRESSIONTEST_DOUBLE $REGRESSIONTEST_MPI_RUN_COMMAND -xml all
