GitLab
======

|Gromacs| is transitioning to GitLab for source code management, issue tracking,
and integrated automation for testing and documentation.

The repository contains DockerFiles and GitLab Runner configuration
files to support automated testing and documentation builds.
General information on configuring GitLab CI pipelines can be found
in the official `Gitlab documentation <https://docs.gitlab.com/ee/ci/yaml/>`_.

The GitLab CI configuration entry point is the :file:`.gitlab-ci.yml` file
at the root of the source tree.
Configuration templates are found in the files in the
:file:`admin/ci-templates/` directory.

Docker images used by GitLab Runner are available on `Docker Hub <https://hub.docker.com/u/gromacs>`__.
Images are (re)built manually from DockerFiles in :file:`admin/dockerfiles`.

This documentation is incomplete, pending resolution of :issue:`3275`.

..  todo:: Expand this documentation to resolve :issue:`3275`
