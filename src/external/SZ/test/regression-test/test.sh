#!/bin/bash

if [[ $# < 2 ]] ; then
	echo "Usage: test.sh [config_file] [SZ root package]"
	echo "Example: test.sh sz.config /home/sdi/Development/SZ_C_Version/sz-1.4.9-beta-normalsize"
	exit
fi

java -cp lib/SZ_RegressionTest.jar test.CheckConfiguration $1 $2
