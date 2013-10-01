#!/bin/sh
if [ -z "$3" ]; then
    echo $0 STARTTEST ENDTEST TNGFILEDIR
else
    STARTTEST="$1"
    ENDTEST="$2"
    TNGFILEDIR="$3"
    for testnum in $(seq $STARTTEST $ENDTEST); do
	if [ -r $TNGFILEDIR/test$testnum.tng_compress ]; then
	    grep -v "EXPECTED_FILESIZE" test$testnum.h >tmp$$.h
	    echo "#define EXPECTED_FILESIZE" $(ls -l $TNGFILEDIR/test$testnum.tng_compress |awk '{print $5}'). >>tmp$$.h
	    mv tmp$$.h test$testnum.h
	fi
    done
fi