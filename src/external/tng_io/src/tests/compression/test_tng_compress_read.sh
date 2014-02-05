#!/bin/sh
numtests=78
for x in $(seq 1 $numtests); do
    ./test_tng_compress_read$x
done