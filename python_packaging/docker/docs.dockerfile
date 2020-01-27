# Docker build context should be the GROMACS repository root. I.e.
#
#     docker build -t gmxapi/docs -f docs.dockerfile ../..
#
# Note that this image depends on the image gmxapi/ci-mpich, which is built from
# ci.dockerfile. Build errors will occur if the current repository is too
# different from the version used to build gmxapi/ci-mpich. Either (re)build a
# local copy of gmxapi/ci-mpich or specify a particular tagged version
# gmxapi/ci-mpich:$REF by passing the build argument REF set to one of the tags
# at https://cloud.docker.com/u/gmxapi/repository/docker/gmxapi/ci-mpich
# I.e.
#
#     docker build -t gmxapi/docs -f docs.dockerfile --build-arg REF=sometag ../..
#
# Launch and bind the host port 8080 to the web server in the container.
#
#     docker run --rm -p 8080:80 gmxapi/docs
#
# Connect by browsing to http://localhost:8080/
#
# Web service is provided by a base `httpd` Docker image. Refer to
# https://hub.docker.com/_/httpd for additional details.
#

ARG REF=latest
FROM gmxapi/gromacs-dependencies-mpich:$REF as docbuild-base

RUN apt-get update && \
    apt-get -yq --no-install-suggests --no-install-recommends install \
        plantuml \
        python3 \
        python3-venv && \
    rm -rf /var/lib/apt/lists/*


FROM gmxapi/ci-mpich:$REF as cibuild

RUN . $VENV/bin/activate && \
    pip install -r /home/testing/gmxapi/requirements-docs.txt --no-cache-dir


FROM docbuild-base as docbuild
COPY --from=cibuild /usr/local/gromacs /usr/local/gromacs

RUN groupadd -r testing && useradd -m -s /bin/bash -g testing testing
USER testing
ENV HOME /home/testing
ENV VENV $HOME/venv

COPY --chown=testing:testing --from=cibuild $VENV $VENV

COPY --chown=testing:testing docs/gmxapi $HOME/documentation
COPY --chown=testing:testing python_packaging/documentation/conf.py $HOME/documentation
RUN cd $HOME && \
    . $VENV/bin/activate && \
    sphinx-build -b html documentation html


FROM httpd

COPY --from=docbuild /home/testing/html/ /usr/local/apache2/htdocs/
