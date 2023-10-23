#!/bin/bash

if [ ! -d "SZ-travis-testdata-master" ]; then
	#wget http://www.mcs.anl.gov/~shdi/download/travis-testdata.tar.gz
	#git clone https://github.com/disheng222/SZ-travis-testdata.git
	wget https://github.com/disheng222/SZ-travis-testdata/archive/master.zip
	unzip master.zip
	rm master.zip
	cd SZ-travis-testdata-master
	tar -xzvf travis-testdata.tar.gz	
	rm travis-testdata.tar.gz
fi
