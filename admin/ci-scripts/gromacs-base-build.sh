#!/usr/bin/env bash
set -e
CMAKE=${CMAKE:-$(which cmake)}
echo "Running on host:" $KUBERNETES_HOSTNAME
cd $BUILD_DIR

# If there is a parameter, use that value, else default to $KUBERNETES_CPU_LIMIT
BUILD_PARALLELISM=${1:-$KUBERNETES_CPU_LIMIT}
BUILD_TEST_PARALLELISM=${2:-$KUBERNETES_CPU_LIMIT}

$CMAKE --build . -- -j$BUILD_PARALLELISM 2>&1 | tee buildLogFile.log

# Accepts string containing a space-separated list of CMake test
# targets to build, which can be "tests" (in order to make all of
# them).  "tests" is also the default if the environment variable
# GMX_TESTS_TO_BUILD is empty.
$CMAKE --build . --target ${GMX_TESTS_TO_BUILD:-tests} -- -j$BUILD_TEST_PARALLELISM 2>&1 | tee testBuildLogFile.log

# Find compiler warnings
awk '/warning/,/warning.*generated|^$/' buildLogFile.log testBuildLogFile.log \
      | grep -v "CMake" | tee buildErrors.log || true
grep "cannot be built" buildLogFile.log testBuildLogFile.log | tee -a buildErrors.log || true
grep "fatal error" buildLogFile.log testBuildLogFile.log | tee -a buildErrors.log || true
grep "error generated when compiling" buildLogFile.log testBuildLogFile.log | tee -a buildErrors.log || true
grep "error:" buildLogFile.log testBuildLogFile.log | tee -a buildErrors.log || true

# Find linking errors:
grep "^/usr/bin/ld:" buildLogFile.log testBuildLogFile.log | tee -a buildErrors.log || true

# Install GROMACS
set +e -o pipefail # Make "$?" work correctly
$CMAKE --build . --target install 2>&1 | tee installBuildLogFile.log
EXITCODE=$?
set -e +o pipefail

# Fail if there were warnings or errors reported
if [ -s buildErrors.log ] || [ $EXITCODE != 0 ] ; then echo "Found compiler warning during build"; cat buildErrors.log; exit 1; fi
# Remove object files to minimize artifact size
find . -mindepth 1 -name '*.o' ! -type l -printf '%p\n' -delete 2>&1 > remove-build-objects.log
cd ..
