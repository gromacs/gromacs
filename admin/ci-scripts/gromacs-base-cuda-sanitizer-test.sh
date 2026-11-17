#!/usr/bin/env bash
set -eo pipefail

CMAKE="${CMAKE:-$(which cmake)}"
CTEST="$(dirname $CMAKE)/ctest"
echo "Running on host:" $KUBERNETES_HOSTNAME
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

# Let test binaries know they are running in the GROMACS CI compute-sanitizer job
export GMX_CI_COMPUTE_SANITIZER=1

# Each pass would overwrite the logs and XML files of the previous one in Testing/,
# so after each pass its output is moved to a directory per tool. CTest recreates
# the directories in the next pass. Testing/ is archived as a job artifact.
PER_TOOL_RESULTS_DIR="Testing/ComputeSanitizer"
rm -rf "${PER_TOOL_RESULTS_DIR}"
# List of the CTest result files of all passes, combined into one JUnit report below
RUNS_LIST="${PER_TOOL_RESULTS_DIR}/runs.xml"
mkdir -p "${PER_TOOL_RESULTS_DIR}"
echo "<runs>" > "${RUNS_LIST}"

for TOOL in memcheck racecheck synccheck initcheck; do
    echo "Running CUDA Compute Sanitizer in ${TOOL} mode"
    "${CTEST}" -T MemCheck \
      --overwrite MemoryCheckCommand="${COMPUTE_SANITIZER_BIN}" \
      --overwrite MemoryCheckCommandOptions="--tool=${TOOL} ${COMPUTE_SANITIZER_FLAGS}" \
      --overwrite MemoryCheckType=CudaSanitizer \
      --label-regex "${TEST_LABELS}" \
      | tee ctestLog.log || TOOLS_FAILED="${TOOL},${TOOLS_FAILED}"

    TAG_DIR="$(head -n 1 < Testing/TAG)"
    TOOL_RESULTS_DIR="${PER_TOOL_RESULTS_DIR}/${TOOL}"
    mkdir -p "${TOOL_RESULTS_DIR}"
    mv ctestLog.log Testing/Temporary Testing/"${TAG_DIR}" "${TOOL_RESULTS_DIR}/"
    # The href is relative to the list file
    if [ -f "${TOOL_RESULTS_DIR}/${TAG_DIR}/DynamicAnalysis-Test.xml" ] ; then
        echo "  <run tool=\"${TOOL}\" href=\"${TOOL}/${TAG_DIR}/DynamicAnalysis-Test.xml\"/>" >> "${RUNS_LIST}"
    fi
done
echo "</runs>" >> "${RUNS_LIST}"

# Report the results of all passes, with one test suite per tool
xsltproc "${CI_PROJECT_DIR}/admin/ci-scripts/CTestSanitizerRuns2JUnit.xsl" "${RUNS_LIST}" > JUnitTestResults.xml

if [ -n "${TOOLS_FAILED}" ] ; then
    echo -e '\n'
    echo "The following Compute Sanitizer tools failed: ${TOOLS_FAILED} more details above."
    echo -e '\n\n'
    exit 1
fi
