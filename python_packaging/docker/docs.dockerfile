# Docker build context should be one directory up from the "docker" directory. I.e.
#
#     docker build -t gmxapi/docs -f docs.dockerfile ..
#
# Launch and bind the host port 8080 to the web server in the container.
#
#     docker run --rm -tp 8080:80 gmxapi/docs
#
# Connect by browsing to http://localhost:8080/
#

FROM gmxapi/ci-mpich as docsbuild

USER root

RUN apt-get update && \
    apt-get -yq --no-install-suggests --no-install-recommends install \
        plantuml && \
    rm -rf /var/lib/apt/lists/*

USER testing

RUN . $VENV/bin/activate && \
    pip install -r /home/testing/gmxapi/requirements-docs.txt --no-cache-dir

COPY documentation /home/testing/gmxapi/documentation
RUN cd /home/testing/gmxapi && \
    . $VENV/bin/activate && \
    sphinx-build -b html documentation html


FROM httpd

COPY --from=docsbuild /home/testing/gmxapi/html/ /usr/local/apache2/htdocs/
