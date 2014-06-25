#!/bin/sh

# Change to root of the source tree, no matter from where the script is
# invoked.
cd `dirname $0`
cd ../../

for destdir in analysisdata selection trajectoryanalysis ; do
    cp -f src/testutils/common-referencedata.xsl \
          src/gromacs/$destdir/tests/refdata/
done

for destdir in trajectoryanalysis ; do
    cp -f src/gromacs/analysisdata/tests/refdata/analysisdata-referencedata.xsl \
          src/gromacs/$destdir/tests/refdata/
done
