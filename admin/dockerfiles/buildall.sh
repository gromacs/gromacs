#!/bin/bash

export TARGET=$0
export TARGET_VERSION=$1
export MATRIX="$TARGET-$TARGET_VERSION"
export SLUG="ci-$MATRIX"

docker login

tags[0]=gromacs/base:2020
docker pull ${tags[0]} || true
docker build -t ${tags[0]} --cache-from ${tags[0]} base

tool=clang
for tool_version in 6 7 8; do
  MATRIX="$tool-$tool_version"
  SLUG="ci-$MATRIX"
  tag=gromacs/$SLUG:2020
  tags[${#tags[@]}]=$tag
  docker build \
    -t $tag \
    --build-arg TOOL_VERSION=$tool_version \
    ci-$tool
done

tool=gcc
for tool_version in 5 6 7 8; do
  MATRIX="$tool-$tool_version"
  SLUG="ci-$MATRIX"
  tag=gromacs/$SLUG:2020
  tags[${#tags[@]}]=$tag
  docker build \
    -t $tag \
    --build-arg TOOL_VERSION=$tool_version \
    ci-$tool
done

tag=gromacs/ci-docs-clang:2020
tags[${#tags[@]}]=$tag
docker build -t $tag \
             ci-docs-clang

tag=gromacs/ci-docs-gcc:2020
tags[${#tags[@]}]=$tag
docker build -t $tag \
             ci-docs-gcc

for tag in ${tags[@]}; do
  docker push $tag
done
