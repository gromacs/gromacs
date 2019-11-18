# Provide an easy-to-reproduce environment in which to test full Python functionality.

# Run with default command and map the http port to the host.
#
#    docker run --rm -ti -p 8888:8888 gmxapi/notebook
#
# Building:
#
# Requires Docker 17.05 or higher.
#
# Use this (the "docker" directory) as the context directory.
#
#    docker build -t gmxapi/notebook -f notebook.dockerfile .
#

ARG REF=latest
FROM gmxapi/ci-mpich:$REF

RUN . $VENV/bin/activate && \
    pip install --no-cache-dir jupyter

ADD --chown=testing:testing notebook /docker_entry_points

CMD ["notebook"]
