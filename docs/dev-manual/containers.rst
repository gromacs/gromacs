==========
Containers
==========

|Gromacs| project infrastructure uses Docker containerization to
isolate automated tasks.
A number of images are maintained to provide a breadth of testing coverage.

Scripts and configuration files for building images are stored in the repository
under :file:`admin/containers/`
Images are (re)built manually by |Gromacs| project staff and pushed to
DockerHub and GitLab.
See https://hub.docker.com/u/gromacs and https://gitlab.com/gromacs/gromacs/container_registry

GitLab Container Registry
=========================

CI Pipelines use a GitLab container registry instead of pulling from Docker Hub.

Project members with role ``Developer`` or higher privilege can
`push images <https://docs.gitlab.com/ee/user/packages/container_registry/index.html#build-and-push-images-by-using-docker-commands>`__
to the container registry.

Steps:

1. Create a `personal access token <https://gitlab.com/-/profile/personal_access_tokens>`__ (`docs <https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html>`__)
   with ``write_registry`` and ``read_registry`` scopes. Save the hash!
2. Authenticate from the command line with ``docker login registry.gitlab.com -u <user name> -p <hash>``
3. ``docker push registry.gitlab.com/gromacs/gromacs/<imagename>``

Refer to :file:`buildall.sh` in the ``main`` branch for the set of images
currently built.

Within :doc:`pipeline jobs <gitlab-ci>`, jobs specify a Docker image with the *image* property.
For image naming convention, see :py:func:`utility.image_name`.
Images from the GitLab registry
are easily accessible with the same identifier as above.
For portability, CI environment variables may be preferable for parts of the image identifier.
Example::

    some_job:
      image: ${CI_REGISTRY_IMAGE}/ci-<configuration>
      ...

For more granularity,
consider equivalent expressions ``${CI_REGISTRY}/${CI_PROJECT_PATH}``
or ``${CI_REGISTRY}/${CI_PROJECT_NAMESPACE}/${CI_PROJECT_NAME}``
Ref: https://docs.gitlab.com/ee/ci/variables/predefined_variables.html

Utilities
=========

:file:`utility.py`
------------------

.. automodule:: utility
    :members:

:file:`scripted_gmx_docker_builds.py`
-------------------------------------

.. automodule:: scripted_gmx_docker_builds
