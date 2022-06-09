#!/usr/bin/env bash
set -eo pipefail

CMAKE="${CMAKE:-$(which cmake)}"
cd "${BUILD_DIR}"

# If $GMX_TEST_REQUIRED_NUMBER_OF_DEVICES is not set and we have GPUs, set it
if [ -z "${GMX_TEST_REQUIRED_NUMBER_OF_DEVICES}" ] && [ -n "${KUBERNETES_EXTENDED_RESOURCE_NAME}" ] ; then
    if grep -q '/gpu$' <<< "${KUBERNETES_EXTENDED_RESOURCE_NAME}"; then
        echo "export GMX_TEST_REQUIRED_NUMBER_OF_DEVICES=\"${KUBERNETES_EXTENDED_RESOURCE_LIMIT}\"";
        export GMX_TEST_REQUIRED_NUMBER_OF_DEVICES="${KUBERNETES_EXTENDED_RESOURCE_LIMIT}";
    fi
fi
if grep -qF 'nvidia.com/gpu' <<< "$KUBERNETES_EXTENDED_RESOURCE_NAME"; then
    nvidia-smi || true
else
    echo "Can not use CUDA Compute Sanitizer without an NVIDIA GPU"
    exit 1
fi

# Path to the compute-sanitizer binary
COMPUTE_SANITIZER_BIN="$(which compute-sanitizer)"
# Common flags: non-zero exit code on error; require that CUDA is actually used in the tests; trace child processes.
COMPUTE_SANITIZER_FLAGS='--error-exitcode=1 --require-cuda-init=yes --target-processes=all'
# Compute Sanitizer slows things down, so we only run a selected subset of tests
TEST_LABELS='QuickGpuTest'
# Flag to mark that any
TOOLS_FAILED=""

# TODO: Initcheck will only be enabled for GROMACS 2023, see #4453
for TOOL in memcheck racecheck synccheck; do
    echo "Running CUDA Compute Sanitizer in ${TOOL} mode"
    ctest -T MemCheck \
      --overwrite MemoryCheckCommand="${COMPUTE_SANITIZER_BIN}" \
      --overwrite MemoryCheckCommandOptions="--tool=${TOOL} ${COMPUTE_SANITIZER_FLAGS}" \
      --overwrite MemoryCheckType=CudaSanitizer \
      --label-regex "${TEST_LABELS}" \
      | tee ctestLog.log || TOOLS_FAILED="${TOOL},${TOOLS_FAILED}"
done

xsltproc "${CI_PROJECT_DIR}/scripts/CTest2JUnit.xsl" Testing/"$(head -n 1 < Testing/TAG)"/*.xml > JUnitTestResults.xml

if [ -n "${TOOLS_FAILED}" ] ; then
    echo -e '\n'
    echo "The following Compute Sanitizer tools failed: ${TOOLS_FAILED} more details above."
    echo -e '\n\n'
    exit 1
fi
