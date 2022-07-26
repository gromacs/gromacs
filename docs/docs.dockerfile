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
#     docker pull registry.gitlab.com/gromacs/gromacs/ci-ubuntu-20.04-llvm-7-docs
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

ARG BASE_IMAGE=ci-ubuntu-20.04-llvm-7-docs
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
RUN . /root/venv/py3.7/bin/activate && \
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
RUN cmake --build . --target webpage -- $BUILD_ARGS


FROM httpd

ENV BUILD_DIR /tmp/gromacs-build
COPY --from=build $BUILD_DIR/docs/html/ /usr/local/apache2/htdocs/
