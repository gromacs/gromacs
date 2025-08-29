#!/usr/bin/env bash
set -eo pipefail

CMAKE="${CMAKE:-$(which cmake)}"
CTEST="$(dirname $CMAKE)/ctest"
cd "${BUILD_DIR}"

# If $GMX_TEST_REQUIRED_NUMBER_OF_DEVICES is not set and we have GPUs, set it
if [ -z "${GMX_TEST_REQUIRED_NUMBER_OF_DEVICES}" ] && [ -n "${GPU_VENDOR}" ] ; then
    if grep -q 'NVIDIA' <<< "${GPU_VENDOR}"; then
        echo "export GMX_TEST_REQUIRED_NUMBER_OF_DEVICES=\"${GPU_COUNT}\"";
        export GMX_TEST_REQUIRED_NUMBER_OF_DEVICES="${GPU_COUNT}";
    fi
fi
if grep -qF 'NVIDIA' <<< "$GPU_VENDOR"; then
    nvidia-smi -L && nvidia-smi || true
else
    echo "Can not use CUDA Compute Sanitizer without an NVIDIA GPU"
    exit 1
fi

# Path to the compute-sanitizer binary
COMPUTE_SANITIZER_BIN="$(which compute-sanitizer)"
# Common flags: non-zero exit code on error; require that CUDA is actually used in the tests; trace child processes.
COMPUTE_SANITIZER_FLAGS='--error-exitcode=1 --target-processes=all'
# Compute Sanitizer slows things down, so we only run a selected subset of tests
TEST_LABELS='QuickGpuTest'
# Flag to mark that any
TOOLS_FAILED=""

for TOOL in memcheck ; do # racecheck synccheck initcheck; do
    echo "Running CUDA Compute Sanitizer in ${TOOL} mode"
    "${CTEST}" -T MemCheck \
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
