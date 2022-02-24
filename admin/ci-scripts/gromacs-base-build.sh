#!/usr/bin/env bash
set -e
CMAKE=${CMAKE:-$(which cmake)}
cd $BUILD_DIR
$CMAKE --build . -- -j$KUBERNETES_CPU_LIMIT 2>&1 | tee buildLogFile.log
$CMAKE --build . --target tests -- -j$KUBERNETES_CPU_LIMIT 2>&1 | tee testBuildLogFile.log

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
$CMAKE --build . --target install 2>&1 | tee installBuildLogFile.log

EXITCODE=$?

# Fail if there were warnings or errors reported
if [ -s buildErrors.log ] || [ $EXITCODE != 0 ] ; then echo "Found compiler warning during build"; cat buildErrors.log; exit 1; fi
# Remove object files to minimize artifact size
find . -mindepth 1 -name '*.o' ! -type l -printf '%p\n' -delete 2>&1 > remove-build-objects.log
cd ..
