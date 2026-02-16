#!/usr/bin/env bash
set -e
echo "Running on host:" $KUBERNETES_HOSTNAME
cd $BUILD_DIR
# Needed to run MPI enabled code in the docker images, until we set up different users
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export ASAN_OPTIONS="check_initialization_order=1:detect_invalid_pointer_pairs=1:strict_init_order=true:strict_string_checks=true:detect_stack_use_after_return=true"
export UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1:suppressions=$CI_PROJECT_DIR/admin/ubsan-suppressions.txt

# If $GMX_TEST_REQUIRED_NUMBER_OF_DEVICES is not set and we have GPUs, set it
if [[ -z "$GMX_TEST_REQUIRED_NUMBER_OF_DEVICES" ]] && [[ -n "$GPU_VENDOR" ]] ; then
    echo "export GMX_TEST_REQUIRED_NUMBER_OF_DEVICES=\"$GPU_COUNT\"";
    export GMX_TEST_REQUIRED_NUMBER_OF_DEVICES="$GPU_COUNT";
fi

export PARALLEL_TEST_EXECUTION=$KUBERNETES_CPU_LIMIT
if grep -qF 'NVIDIA' <<< "$GPU_VENDOR"; then
    nvidia-smi -L && nvidia-smi || true;
    if [[ "$GMX_ENABLE_NVSHMEM" != "" ]] && [[ "$GPU_COUNT" -eq "2" ]]; then
        # In CI with dual GPUs NVSHMEM cannot support more than 2 MPI processes/GPU
        # which happens in multi sim tests so we disable them
        echo "Disabling MdrunMultiSim tests as with GMX_ENABLE_NVSHMEM it does not work on dual GPU setup without MPS"
        EXTRA_FLAGS="--exclude-regex MdrunMultiSim "
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
    # Speed up device re-initialization, especially when running multiple tests in parallel
    export CUDA_DEVICE_MAX_CONNECTIONS=2  # default is 8
fi
if grep -qF 'AMD' <<< "$GPU_VENDOR"; then
    clinfo -l || true;
    rocm-smi || true;
    if [[ "$GPU_COUNT" -eq "2" ]]; then
        # Our CI has some issues with ROCm (and CUDA) IPC; disable rocm_ipc
        export UCX_TLS="sm,self,rocm_copy"
    fi
    export HSA_ENABLE_SDMA=0  # Work around CI issues, Issue #5341
    export GPU_MAX_HW_QUEUES=2  # Prevent "amdgpu: Runlist is getting oversubscribed", Issue #5341
    export PARALLEL_TEST_EXECUTION=1 # Prevent tests getting stuck
fi
if grep -qF 'INTEL' <<< "$GPU_VENDOR"; then
    sycl-ls || true;
    export SYCL_CACHE_PERSISTENT=1  # Issue #4218
fi

LABEL_REGEX=
if [[ -n "$GMX_TEST_LABELS" ]] ; then
    LABEL_REGEX="--label-regex $GMX_TEST_LABELS"
fi

TESTS_REGEX=
# Check GMX_TESTS_TO_RUN_REGEX to see whether the caller wants to run
# all tests, or just a subset.
if [[ -n "$GMX_TESTS_TO_RUN_REGEX" ]] ; then
    TESTS_REGEX="--tests-regex $GMX_TESTS_TO_RUN_REGEX"
fi

ctest -D $CTEST_RUN_MODE $LABEL_REGEX $TESTS_REGEX $EXTRA_FLAGS --parallel $PARALLEL_TEST_EXECUTION --schedule-random --output-on-failure | tee ctestLog.log || true

EXITCODE=$?

if [[ "$CTEST_RUN_MODE" == "ExperimentalMemCheck" ]] ; then
    TEST_XML_OUTPUT="DynamicAnalysis.xml"
else
    TEST_XML_OUTPUT="Test.xml"
fi

awk '/The following tests FAILED/,/^Errors while running CTest|^$/' ctestLog.log | tee ctestErrors.log
xsltproc $CI_PROJECT_DIR/scripts/CTest2JUnit.xsl Testing/`head -n 1 < Testing/TAG`/$TEST_XML_OUTPUT > JUnitTestResults.xml
if [ -s ctestErrors.log ] || [ $EXITCODE != 0 ] ; then
    echo "Error during running ctest";
    exit 1;
fi
cd .
