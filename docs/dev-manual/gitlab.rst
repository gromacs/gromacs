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

Pipeline execution
------------------

.. todo:: Discuss the distinct characteristics of |Gromacs| CI pipelines to relevant to job configuration.

.. todo:: Comment on the number of pipelines that can be or which are likely to be running at the same time.

Variables
~~~~~~~~~

The GitLab CI framework, GitLab Runner, plugins, and our own scripts set and
use several variables (usually as a key under a *variables* parameter in
the YAML configuration).

Some default values are specified for all jobs in :file:`.gitlab-ci.yml`.
Many of the mix-in / template jobs in :file:`admin/gitlab-ci/global.gitlab-ci.yml`
provide additional or overriding definitions.
Other variables may be set when making final job definitions.

Variables may control the behvior of GitLab-CI (those beginning with ``CI_``),
GitLab Runner and supporting infrastructure, or may be used by job definitions,
or passed along to the environment of executed commands.

*variables* keys beginning with ``KUBERNETES_`` relate to the GitLab Runner
`Kubernets executor <https://docs.gitlab.com/runner/executors/kubernetes.html#the-kubernetes-executor>`__

Other important variable keys are as follows.

.. glossary::
    CMAKE_MPI_OPTIONS
        Provide CMake command line arguments to define GROMACS MPI build options.

.. todo:: Define common variables.
    ``BUILD_DIR``, ``INSTALL_DIR``, ``CACHE_FALLBACK_KEY``, ...
