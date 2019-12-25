#!/bin/bash

export TARGET=$0
export TARGET_VERSION=$1
export MATRIX="$TARGET-$TARGET_VERSION"
export SLUG="ci-$MATRIX"

docker login

docker pull gromacs/base || true
docker build -t gromacs/base --cache-from gromacs/base base
docker push gromacs/base

tool=clang
for tool_version in 6 7 8; do
  MATRIX="$tool-$tool_version"
  SLUG="ci-$MATRIX"
  docker build \
    -t gromacs/continuous-integration:$SLUG \
    --build-arg TOOL_VERSION=$tool_version \
    ci-$tool
done

tool=gcc
for tool_version in 5 6 7 8; do
  MATRIX="$tool-$tool_version"
  SLUG="ci-$MATRIX"
  docker build \
    -t gromacs/continuous-integration:$SLUG \
    --build-arg TOOL_VERSION=$tool_version \
    ci-$tool
done

docker build -t gromacs/continuous-integration:ci-docs-clang \
             ci-docs-clang

docker build -t gromacs/continuous-integration:ci-docs-gcc \
             ci-docs-gcc

docker push gromacs/continuous-integration
