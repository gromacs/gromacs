#!/bin/bash

if [ -z "$*" ]
then
    echo "usage: `basename "$0"` name-of-target"
    exit
fi

docker login
docker build $1 --target $1 -t gromacs/gromacs:$1
docker push gromacs/gromacs:$1
