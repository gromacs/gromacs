#!/usr/bin/env bash
set -e
CMAKE=${CMAKE:-$(which cmake)}
cd $BUILD_DIR
$CMAKE --build . -- -j$KUBERNETES_CPU_LIMIT 2>&1 | tee buildLogFile.log
$CMAKE --build . --target tests -- -j$KUBERNETES_CPU_LIMIT 2>&1 | tee testBuildLogFile.log
awk '/warning/,/warning.*generated|^$/' buildLogFile.log testBuildLogFile.log \
      | grep -v "CMake" | tee buildErrors.log || true
grep "cannot be built" buildLogFile.log testBuildLogFile.log | tee -a buildErrors.log || true
$CMAKE --build . --target install 2>&1 | tee installBuildLogFile.log
if [ -s buildErrors.log ] ; then echo "Found compiler warning during build"; cat buildErrors.log; exit 1; fi
for file in `find . -mindepth 1 -name "*.o" ! -type l` ; do echo $file ; rm $file ; done 2>&1 > remove-build-objects.log
cd ..
