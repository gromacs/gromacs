# Build and serve documentation for the current working files in a tidy container.
#
# Use the CI image from https://gitlab.com/gromacs/gromacs/container_registry/
# and build the GROMACS 'webpage' CMake target. Copy the web content to an
# httpd image.
#
# Docker build context should be the local GROMACS repository root. I.e.
#
#     docker build -t gromacs/docs:some-useful-tag -f docs.dockerfile ..
#
# If you do not have authorization to `docker push` to the `gromacs/` space
# on DockerHub, be sure to tag the image with a more appropriate name,
# or invite collaborators to build the image locally, themselves.
#
# Note that initial download of the base build image from registry.gitlab.com
# may take quite a while, but it will be cached for subsequent builds.
# The final web server image should be much smaller.
# The base image is fairly stable. If you need to update the `latest` tag,
# though, you can explicitly
#     docker pull registry.gitlab.com/gromacs/gromacs/ci-ubuntu-20.04-llvm-9-docs
# To use a specific base image from
# https://gitlab.com/gromacs/gromacs/container_registry/?search[]=docs,
# specify the image and/or tag with the BASE_IMAGE and BASE_TAG build args,
# respectively. E.g.
#     docker build -t gromacs/docs --build-arg BASE_TAG=release-2021 -f docs.dockerfile ..
# Alternatively, build a suitable base image locally, such as by providing
# `--doxygen --mpi` to `admin/containers/scripted_gmx_docker_builds.py`.
# (See `admin/containers/buildall.sh`)
#
# Peruse the ARG lines below for variables that can be passed to the build
# via --build-arg. In particular, it can be very helpful to direct parallel
# builds with BUILD_ARGS, which is appended to the `cmake --build` line after
# `--` to be passed along to the build system.
#    docker build -t gromacs/docs -f docs.dockerfile --build-arg BUILD_ARGS="-j4" ..
# You may have to set a smaller number of cores if your Docker build environment
# has limited memory.
#
# Launch and bind the host port 8080 to the web server in the container.
#
#     docker run --rm -p 8080:80 gromacs/docs
#
# Connect by browsing to http://localhost:8080/
#
# Web service is provided by a base `httpd` Docker image. Refer to
# https://hub.docker.com/_/httpd for additional details.
#
# Hint: For minor updates copy new content instead of rebuilding.
#
# The image can take a long time to rebuild since a full gromacs build is
# required (any additional caching would be unnecessarily complex and may
# hurt reproducibility). To update the web content,
#  1. create a container from an existing image
#  2. copy a local docs build into the container
#  3. save a snapshot of the container as a new image
#
#     BUILD=/path/to/gromacs/build
#     docker create --name docs_container gromacs/docs:my-original
#     docker cp $BUILD/docs/html/. docs_container:/usr/local/apache2/htdocs/
#     docker commit docs_container gromacs/docs:my-update
#     docker rm docs_container
#
# Note that the exact syntax for the paths to `docker cp` is a bit sensitive
# when the destination directory exists.
# Ref https://docs.docker.com/engine/reference/commandline/cp/#description
#
# If you push images updated this, it will be less expensive for others
# to pull successive updates. But try to update the same base image when only
# content is changing so that you don't add unnecessary extra layers that bloat
# the overall docker image size for first-time pullers.

ARG BASE_IMAGE=ci-ubuntu-20.04-llvm-9-docs
ARG BASE_TAG=latest
FROM registry.gitlab.com/gromacs/gromacs/$BASE_IMAGE:$BASE_TAG as build

# We let the sources stay in the intermediate `build` stage, since a little bloat here shouldn't
#impact the size of the final containers unless we copy explicitly (see below)
ENV SRC_DIR /gromacs-source
COPY . $SRC_DIR

# Note that BUILD_DIR must be duplicated for the final build below.
ENV BUILD_DIR /tmp/gromacs-build

RUN mkdir -p $BUILD_DIR
WORKDIR $BUILD_DIR

# Enable (or disable) Sphinx "todo" output.
ARG TODOS=1
# Allow the build type to be specified with `docker build --build-arg TYPE=something`
ARG TYPE=Release
# Allow arbitrary CMake args with `--build-arg CMAKE_ARGS="..."`
ARG CMAKE_ARGS=""
RUN . /root/venv/py3.9/bin/activate && \
    cmake $SRC_DIR \
        -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs \
        -DGMXAPI=ON \
        -DGMX_PYTHON_PACKAGE=ON \
        -DGMX_THREAD_MPI=ON \
        -DGMX_USE_RDTSCP=OFF \
        -DGMX_INSTALL_LEGACY_API=ON \
        -DCMAKE_BUILD_TYPE=$TYPE \
        -DSPHINX_CONFIG_OVERRIDES="-Dtodo_include_todos=$TODOS" \
        $CMAKE_ARGS

# Additional arguments to pass to the build system.
ARG BUILD_ARGS=""

RUN cmake --build . -- $BUILD_ARGS

RUN cmake --build . --target webpage -- $BUILD_ARGS


FROM httpd

ENV BUILD_DIR /tmp/gromacs-build
COPY --from=build $BUILD_DIR/docs/html/ /usr/local/apache2/htdocs/
