#!/usr/bin/env bash
set -e
CMAKE=${CMAKE:-$(which cmake)}
cd $BUILD_DIR
export UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1:suppressions=$CI_PROJECT_DIR/admin/ubsan-suppressions.txt
# Needed to run MPI enabled code in the docker images, until we set up different users
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export ASAN_OPTIONS="check_initialization_order=1:detect_invalid_pointer_pairs=1:strict_init_order=true:strict_string_checks=true:detect_stack_use_after_return=true"
# If $GMX_TEST_REQUIRED_NUMBER_OF_DEVICES is not set and we have GPUs, set it
if [ -z $GMX_TEST_REQUIRED_NUMBER_OF_DEVICES ] && [ -n $KUBERNETES_EXTENDED_RESOURCE_NAME ] ; then
    if grep -q '/gpu$' <<< "$KUBERNETES_EXTENDED_RESOURCE_NAME"; then
        echo "export GMX_TEST_REQUIRED_NUMBER_OF_DEVICES=\"$KUBERNETES_EXTENDED_RESOURCE_LIMIT\"";
        export GMX_TEST_REQUIRED_NUMBER_OF_DEVICES="$KUBERNETES_EXTENDED_RESOURCE_LIMIT";
    fi
fi
if grep -qF 'nvidia.com/gpu' <<< "$KUBERNETES_EXTENDED_RESOURCE_NAME"; then
    nvidia-smi || true;
fi
if grep -qF 'amd.com/gpu' <<< "$KUBERNETES_EXTENDED_RESOURCE_NAME"; then
    clinfo -l || true;
fi
if grep -qF 'intel.com/gpu' <<< "$KUBERNETES_EXTENDED_RESOURCE_NAME"; then
    sycl-ls || true;
    export SYCL_CACHE_PERSISTENT=1; # Issue #4218
fi
ctest -D $CTEST_RUN_MODE --output-on-failure | tee ctestLog.log || true
awk '/The following tests FAILED/,/^Errors while running CTest|^$/' ctestLog.log | tee ctestErrors.log
xsltproc $CI_PROJECT_DIR/scripts/CTest2JUnit.xsl Testing/`head -n 1 < Testing/TAG`/*.xml > JUnitTestResults.xml
if [ -s ctestErrors.log ] ; then
    echo "Error during running ctest";
    exit 1;
fi
cd .
