GitLab
======

The repository contains DockerFiles and GitLab Runner configuration
files to support automated testing and documentation builds.
General information on configuring GitLab CI pipelines can be found
in the official `Gitlab documentation <https://docs.gitlab.com/ee/ci/yaml/>`_.

The GitLab CI configuration entry point is the :file:`.gitlab-ci.yml` file
at the root of the source tree.
Configuration templates are found in the files in the
:file:`admin/ci-templates/` directory.

Docker images used by GitLab Runner are available on `Docker Hub <https://hub.docker.com/u/gromacs>`__.
Images are (re)built manually using details in :file:`admin/containers`.

This documentation is incomplete, pending resolution of :issue:`3275`.

..  todo:: Expand this documentation to resolve :issue:`3275`

Pipeline execution
------------------

.. todo:: Discuss the distinct characteristics of |Gromacs| CI pipelines to relevant to job configuration.

.. todo:: Comment on the number of pipelines that can be or which are likely to be running at the same time.

.. note::

    Full automated testing is only available for merge requests originating from
    branches of the main https://gitlab.com/gromacs/gromacs repository.
    GitLab CI pipelines created for forked repositories will include fewer jobs
    in the testing pipeline. Non-trivial merge requests may need to be issued
    from a branch in the ``gromacs`` project namespace in order to receive
    sufficient testing before acceptance.

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

    cache
        There is no global default, but jobs that build software will likely
        set *cache*. To explicitly unset *cache* directives, specify a job
        parameter of ``cache: {}``.
        Refer to `GitLab docs <https://docs.gitlab.com/ee/ci/yaml/#cache>`__
        for details. In particular, note the details of cache identity according
        to `cache:key <https://docs.gitlab.com/ee/ci/yaml/#cachekey>`__

    image
        Part of the tool chain configuration. Instead of setting *image*
        directly, *extend* a *.use_<toolchain>* template from
        :file:`admin/gitlab-ci/global.gitlab-ci.yml`

    rules
    only
    except
    when
        *Job* parameters for controlling the circumstances under which jobs run.
        (Some key words may have different meanings when occurring as elements
        of other parameters, such as *archive:when*, to which this note is not
        intended to apply.)
        Instead of setting any of these directly in a job definition, try to use
        one of the pre-defined behaviors (defined as ``.rules:<something>`` in
        :file:`admin/gitlab-ci/global.gitlab-ci.yml`).
        Errors or unexpected behavior will occur if you specify more than one
        *.rules:...* template, or if you use these parameters in combination
        with a *.rules...* template.
        To reduce errors and unexpected behavior, restrict usage of these controls
        to regular job definitions (don't use in "hidden" or parent jobs).

    tags
        By `default <https://docs.gitlab.com/ee/ci/yaml/#setting-default-parameters>`__,
        jobs require the ``k8s-scilifelab`` tag, which identifies Runners in the
        |Gromacs| infrastructure. A small number of jobs in the first pipeline
        stage override the default with an empty tag list so that all GitLab
        users can run basic tests in their forked project.

    variables
        Many job definitions will add or override keys in *variables*.
        Refer to `GitLab <https://docs.gitlab.com/ee/ci/yaml/#variables>`__
        for details of the merging behavior. Refer to :ref:`variables` for local usage.

Schedules and triggers
~~~~~~~~~~~~~~~~~~~~~~

Pipeline `schedules <https://gitlab.com/help/ci/pipelines/schedules>`__ are
configured through the GitLab web interface.
Scheduled pipelines may provide different variable definitions through the
environment to jobs that run under the ``schedules``
`condition <https://gitlab.com/help/ci/pipelines/schedules#using-only-and-except>`__.

Nightly scheduled pipelines run against ``master`` and *release* branches in
the GROMACS repository.

Global templates
~~~~~~~~~~~~~~~~

In addition to the templates in the main job definition files,
common "mix-in" functionality and behavioral templates are defined in
:file:`admin/gitlab-ci/global.gitlab-ci.yml`.

Jobs beginning with ``.use-`` provide mix-in behavior, such as boilerplate for
jobs using a particular tool chain.

Jobs beginning with a `parameter <https://docs.gitlab.com/ee/ci/yaml>`__
name allow parameters to be set in a single place for common job characteristics.
If providing more than a default parameter value, the job name should be suffixed
by a meaningful descriptor and documented within
:file:`admin/gitlab-ci/global.gitlab-ci.yml`

Job names
~~~~~~~~~

Job names should

1. Indicate the purpose of the job.
2. Indicate relationships between multi-stage tasks.
3. Distinguish jobs in the same stage.
4. Distinguish job definitions throughout the configuration.

Jobs may be reassigned to different stages over time, so including the stage
name in the job name is not helpful, generally. If tags like "pre" and "post,"
or "build" and "test" are necessary to distinguish phases of, say, "webpage,"
then such tags can be buried at the end of the job name.

Stylistically, it is helpful to use delimiters like ``:`` to distinguish the
basic job name from qualifiers or details. Also consider
`grouping jobs <https://docs.gitlab.com/ee/ci/pipelines/index.html#grouping-jobs>`__

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
    CI_PROJECT_NAMESPACE
        Distinguishes pipelines created for repositories in the ``gromacs``
        GitLab project space. May be used to pre-screen jobs to determine
        whether |Gromacs| GitLab infrastructure is available to the pipeline
        before the job is created.

    COMPILER_MAJOR_VERSION
        Integer version number provided by toolchain mix-in for convenience and
        internal use.

    CMAKE_COMPILER_SCRIPT
        CMake command line options for a tool chain. A definition is provided by
        the mix-in toolchain definitions (e.g. ``.use-gcc8``) to be appended to
        :command:`cmake` calls in a job's *script*.

    CMAKE_MPI_OPTIONS
        Provide CMake command line arguments to define GROMACS MPI build options.

    GROMACS_RELEASE
        Read-only environment variable that can be checked to see if a job is
        executing in a pipeline for preparing a tagged release.
        Can be set when launching pipelines via the GitLab web interface.
        For example, see *rules* mix-ins in :file:`admin/gitlab-ci/global.gitlab-ci.yml`.

    EXTRA_INSTALLS
        List additional OS package requirements. Used in *before_script* for some
        mix-in job definitions to install additional software dependencies. If
        using such a job with *extends*, override this variable key with a
        space-delimited list of packages (default: ``""``). Consider proposing a
        patch to the base Docker images to include the dependency to reduce
        pipeline execution time.

.. todo:: Define common variables.
    ``BUILD_DIR``, ``INSTALL_DIR``, ``CACHE_FALLBACK_KEY``, ...
