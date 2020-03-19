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

Configuration files
~~~~~~~~~~~~~~~~~~~

At the root of the repository, :file:`.gitlab-ci.yml` defines the stages and
some default parameters, then includes files from :file:`admin/gitlab-ci/` to
define jobs to be executed in the pipelines.

Note that job names beginning with a period (``.``) are
`"hidden" <https://docs.gitlab.com/ee/ci/yaml/#hidden-keys-jobs>`_.
Such jobs are not directly eligible to run, but may be used as templates
via the `*extends* job property <https://docs.gitlab.com/ee/ci/yaml/#extends>`_.

Job parameters
~~~~~~~~~~~~~~

Refer to https://docs.gitlab.com/ee/ci/yaml for complete documentation on
GitLab CI job parameters, but note the following GROMACS-specific conventions.

.. glossary::

    before_script
        Used by several of our templates to prepend shell commands to
        a job *script* parameter.
        Avoid using *before-script* directly, and be cautious
        about nested *extends* overriding multiple *before_script* definitions.

    image
        Part of the tool chain configuration. Instead of setting *image*
        directly, *extend* a *.use_<toolchain>* template from
        :file:`admin/gitlab-ci/global.gitlab-ci.yml`

    variables
        Many job definitions will add or override keys in *variables*.
        Refer to `GitLab <https://docs.gitlab.com/ee/ci/yaml/#variables>`__
        for details of the merging behavior. Refer to :ref:`variables` for local usage.

In addition to the templates in the main job definition files,
common "mix-in" functionality and behavioral templates are defined in
:file:`admin/gitlab-ci/global.gitlab-ci.yml`.

.. _variables:

Variables
~~~~~~~~~

The GitLab CI framework, GitLab Runner, plugins, and our own scripts set and
use several `variables <https://docs.gitlab.com/ee/ci/variables/README.html>`__.

Default values are available from the ``.variables:default`` definition in
:file:`admin/gitlab-ci/global.gitlab-ci.yml`.
Many of the mix-in / template jobs provide additional or overriding definitions.
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
