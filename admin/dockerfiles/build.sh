#!/bin/bash

if [ -z "$*" ]
then
    echo "usage: `basename "$0"` name-of-target"
    exit
fi

docker login
docker build $1 --target $1 -t gromacs/continuous-integration:$1
docker push gromacs/continuous-integration:$1
