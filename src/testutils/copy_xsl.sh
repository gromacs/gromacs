#!/bin/sh

# Change to root of the source tree, no matter from where the script is
# invoked.
cd `dirname $0`
cd ../../

for destdir in analysisdata selection ; do
    cp -f src/testutils/common-referencedata.xsl \
          src/gromacs/$destdir/tests/refdata/
done

for destdir in analysisdata ; do
    cp -f src/testutils/analysisdata-referencedata.xsl \
          src/gromacs/$destdir/tests/refdata/
done
