#!/usr/bin/env bash
set -e
CMAKE=${CMAKE:-$(which cmake)}
echo "Running on host:" $KUBERNETES_HOSTNAME
cd $BUILD_DIR
export UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1:suppressions=$CI_PROJECT_DIR/admin/ubsan-suppressions.txt
# Needed to run MPI enabled code in the docker images, until we set up different users
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export ASAN_OPTIONS="check_initialization_order=1:detect_invalid_pointer_pairs=1:strict_init_order=true:strict_string_checks=true:detect_stack_use_after_return=true"
# If $GMX_TEST_REQUIRED_NUMBER_OF_DEVICES is not set and we have GPUs, set it
if [[ -z "$GMX_TEST_REQUIRED_NUMBER_OF_DEVICES" ]] && [[ -n "$GPU_VENDOR" ]] ; then
    echo "export GMX_TEST_REQUIRED_NUMBER_OF_DEVICES=\"$GPU_COUNT\"";
    export GMX_TEST_REQUIRED_NUMBER_OF_DEVICES="$GPU_COUNT";
fi
if grep -qF 'NVIDIA' <<< "$GPU_VENDOR"; then
    nvidia-smi -L && nvidia-smi || true;
    if [ "$GMX_CI_DISABLE_CUFFTMP_DECOMPOSITION_ON_INCOMPATIBLE_DEVICES" != "" ] 
    then
        echo "DUE TO LIMITATIONS OF CUFFTMP, THIS JOB RUNS IN DIFFERENT CONFIGURATIONS DEPENDING ON THE VERSION OF GPU AVAILABLE. Now running:" 
        computeCapability=`nvidia-smi -i 0 --query-gpu=compute_cap --format=csv | tail -1 | sed 's/\.//g'`    
        if [ "$computeCapability" -lt "70" ]
        then
            echo "    without PME decomposition, since compute Capability is less than 7.0"
            unset GMX_GPU_PME_DECOMPOSITION
            export LD_PRELOAD=$CUFFTLIB #TODO remove this when cuFFTMp is fixed regarding "regular" ffts for older GPUs #3884
        else
            echo "    with PME decomposition"
        fi
    fi
    if [[ "$GMX_ENABLE_NVSHMEM" != "" ]] && [[ "$GPU_COUNT" -eq "2" ]]
    then
        # In CI with dual GPUs NVSHMEM cannot support more than 2 MPI processes/GPU
        # which happens in multi sim tests so we disable them
        echo "Disabling MdrunMultiSim tests as with GMX_ENABLE_NVSHMEM it does not work on dual GPU setup without MPS"
        EXTRA_FLAGS="--exclude-regex MdrunMultiSim "
    fi
    # Speed up device re-initialization, especially when running multiple tests in parallel
    export CUDA_DEVICE_MAX_CONNECTIONS=2  # default is 8
fi
if grep -qF 'AMD' <<< "$GPU_VENDOR"; then
    clinfo -l || true;
    rocm-smi || true;
    export GPU_MAX_HW_QUEUES=2  # Prevent "amdgpu: Runlist is getting oversubscribed", Issue #5341
fi
if grep -qF 'INTEL' <<< "$GPU_VENDOR"; then
    sycl-ls || true;
    export SYCL_CACHE_PERSISTENT=1; # Issue #4218
fi
LABEL_REGEX=
if [[ -n "$GMX_TEST_LABELS" ]] ; then
    LABEL_REGEX="--label-regex $GMX_TEST_LABELS"
fi
ctest -D $CTEST_RUN_MODE $LABEL_REGEX $EXTRA_FLAGS --parallel $KUBERNETES_CPU_LIMIT --output-on-failure | tee ctestLog.log || true

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
