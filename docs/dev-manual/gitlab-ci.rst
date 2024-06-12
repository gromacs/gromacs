
GitLab CI Pipeline Execution
============================

The repository contains DockerFiles and GitLab Runner configuration
files to support automated testing and documentation builds.
General information on configuring GitLab CI pipelines can be found
in the official `Gitlab documentation <https://docs.gitlab.com/ee/ci/yaml/>`_.

The GitLab CI configuration entry point is the :file:`.gitlab-ci.yml` file
at the root of the source tree.
Configuration templates are found in the files in the
:file:`admin/ci-templates/` directory.

Docker images used by GitLab Runner are available in our
`GitLab Container Registry <https://gitlab.com/gromacs/gromacs/container_registry>`__.
(See :ref:`containers`.)
Images are (re)built manually using details in :file:`admin/containers`.
(See :ref:`gitlab-ci tools`.)

.. todo:: (:issue:`3617`) Comment on the number of pipelines that can be or which are likely to be running at the same time.

.. note::

    Full automated testing is only available for merge requests originating from
    branches of the main https://gitlab.com/gromacs/gromacs repository.
    GitLab CI pipelines created for forked repositories will include fewer jobs
    in the testing pipeline. Non-trivial merge requests may need to be issued
    from a branch in the ``gromacs`` project namespace in order to receive
    sufficient testing before acceptance.

Configuration files
-------------------

At the root of the repository, :file:`.gitlab-ci.yml` defines the stages and
some default parameters, then includes files from :file:`admin/gitlab-ci/` to
define jobs to be executed in the pipelines.

Note that job names beginning with a period (``.``) are
`"hidden" <https://docs.gitlab.com/ee/ci/yaml/#hidden-keys-jobs>`_.
Such jobs are not directly eligible to run, but may be used as templates
via the `*extends* job property <https://docs.gitlab.com/ee/ci/yaml/#extends>`_.

Job parameters
""""""""""""""

Refer to https://docs.gitlab.com/ee/ci/yaml for complete documentation on
GitLab CI job parameters, but note the following |Gromacs|-specific conventions.

.. glossary::

    before_script
        Used by several of our templates to prepend shell commands to
        a job *script* parameter.
        Avoid using *before-script* directly, and be cautious
        about nested *extends* overriding multiple *before_script* definitions.

    job cache
        There is no global default, but jobs that build software will likely
        set *cache*. To explicitly unset *cache* directives, specify a job
        parameter of ``cache: {}``.
        Refer to `GitLab docs <https://docs.gitlab.com/ee/ci/yaml/#cache>`__
        for details. In particular, note the details of cache identity according
        to `cache:key <https://docs.gitlab.com/ee/ci/yaml/#cachekey>`__

    image
        See :ref:`containers` for more about the Docker images used for the
        CI pipelines. If a job depends on artifacts from previous jobs, be sure
        to use the same (or a compatible) image as the dependency!

    rules
    only
    except
    when
        *Job* parameters for controlling the circumstances under which jobs run.
        (Some key words may have different meanings when occurring as elements
        of other parameters, such as *archive:when*, to which this note is not
        intended to apply.)
        Rules in GitLab are special, since the first matching rule will cause
        a job to trigger, and then all remaining rules are ignored. To create
        rules to skip jobs, write rules that use the execution time "never".
        Errors or unexpected behavior will occur if you specify more than one
        *.rules:...* template, or if you use these parameters in combination
        with a *.rules...* template - it is thus NOT possible to combine
        rules through inheritance with the ``extends`` tag.
        Instead, to combine sequences of rules we recommend using a plain
        rules tag where you reference rule entries with the !reference tag,
        e.g. ``!reference [.rules:<something>, rules]``. Each such reference
        can be used as an individual rule in the list.
        To reduce errors and unexpected behavior, restrict usage of these controls
        to regular job definitions (don't use in "hidden" or parent jobs).
        Note that *rules* is not compatible with the older *only* and *except*
        parameters. We have standardized on the (newer) *rules* mechanism.

    tags
        We no longer use any special tags for general (meaning CPU-only)
        |Gromacs| CI jobs, to make sure at least the CPU jobs can still run even
        if somebody clones the repo. For testing you can still add a default tag
        at the start of the top-level :file:`.gitlab-ci.yml`, but this should
        only be used to check that a specific runner works - for production we
        handle it in GitLab instead by selecting what runners accept untagged
        jobs. By default we currently run those on the infrastructure in
        Stockholm for the |Gromacs| project, but please design all CPU jobs so
        they will work on the shared runners too.

    variables
        Many job definitions will add or override keys in *variables*.
        Refer to `GitLab <https://docs.gitlab.com/ee/ci/yaml/#variables>`__
        for details of the merging behavior. Refer to :ref:`variables` for local usage.

Schedules and triggers
""""""""""""""""""""""

Pipeline `schedules <https://gitlab.com/help/ci/pipelines/schedules>`__ are
configured through the GitLab web interface.
Scheduled pipelines may provide different variable definitions through the
environment to jobs that run under the ``schedules``
`condition <https://gitlab.com/help/ci/pipelines/schedules#using-only-and-except>`__.

Nightly scheduled pipelines run against ``main`` and *release* branches in
the |Gromacs| repository.

Some of the rules defined in :file:`rules.gitlab-ci.yml` restrict jobs
to run *only* for scheduled pipelines, or only for *specific* schedules
according to the variables defined for that schedule in the web interface.
For example, the rule element ``if-weekly-then-on-success`` causes a job
to run only if the schedule sets ``GMX_PIPELINE_SCHEDULE=weekly``.

.. admonition:: Running post-merge-acceptance pipelines
    :class: tip

    The Gitlab CI for |Gromacs| runs a set of jobs by default only after a MR has been
    accepted and the resulting commit is included in the target branch if it is ``main``
    or one of the *release* branches. Those jobs can be triggered manually using the
    ``POST_MERGE_ACCEPTANCE`` input variable documented below when executing a new pipeline
    through the Gitlab web interface.

    See also :ref:`trigger-post-merge`.

Global templates
""""""""""""""""

In addition to the templates in the main job definition files,
common "mix-in" functionality and behavioral templates are defined in
:file:`admin/gitlab-ci/global.gitlab-ci.yml`.
For readability, some parameters may be separated into their own files, named
according to the parameter (e.g. :file:`rules.gitlab-ci.yml`).

Jobs beginning with ``.use-`` provide mix-in behavior, such as boilerplate for
jobs using a particular tool chain.

Jobs beginning with a `parameter <https://docs.gitlab.com/ee/ci/yaml>`__
name allow parameters to be set in a single place for common job characteristics.
If providing more than a default parameter value, the job name should be suffixed
by a meaningful descriptor and documented within
:file:`admin/gitlab-ci/global.gitlab-ci.yml`

Job names
"""""""""

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

Updating regression tests
-------------------------

Changes in |Gromacs| that require changes in regression-tests are notoriously hard,
because a merge request that tests against the non-updated version of the
regression tests will necessarily fail, while updating regression tests while
the current change is not integrated into main, might cause other
merge request pipelines to fail.

The solution is a new regression-test branch or commit, uploaded to gitlab.
Then set that regression test branch with REGRESSIONTESTBRANCH or
the specific commit with REGRESSIONTESTCOMMIT when
running the specific pipeline that requires the regressiontest-update.
See below on how to set variables for specific pipelines.

Variables
---------

The GitLab CI framework, GitLab Runner, plugins, and our own scripts set and
use several `variables <https://docs.gitlab.com/ee/ci/variables/README.html>`__.

Default values are available from the top level ``variables`` definition in
:file:`global.gitlab-ci.yml`.
Many of the mix-in / template jobs provide additional or overriding definitions.
Other variables may be set when making final job definitions.

Variables may control the behvior of GitLab-CI (those beginning with ``CI_``),
GitLab Runner and supporting infrastructure, or may be used by job definitions,
or passed along to the environment of executed commands.

*variables* keys beginning with ``KUBERNETES_`` relate to the GitLab Runner
`Kubernets executor <https://docs.gitlab.com/runner/executors/kubernetes.html#the-kubernetes-executor>`__

Other important variable keys are as follows.

.. glossary::
    BUILD_DIR
        |Gromacs| specific directory to perform configuration, building and testing in.
        Usually job dependent, needs to be the same for all tasks of dependent jobs.

    CI_PROJECT_NAMESPACE
        Distinguishes pipelines created for repositories in the ``gromacs``
        GitLab project space. May be used to pre-screen jobs to determine
        whether |Gromacs| GitLab infrastructure is available to the pipeline
        before the job is created.

    COMPILER_MAJOR_VERSION
        Integer version number provided by toolchain mix-in for convenience and
        internal use.

    CMAKE
        ``gromacs/ci-...`` Docker images built after October 2020 have several
        versions of CMake installed. The most recent version of CMake in the
        container will be appear first in ``PATH``. To allow individual jobs to
        use specific versions of CMake, please write the job *script* sections
        using ``$CMAKE`` instead of ``cmake`` and begin the *script* section with
        a line such as ``- CMAKE=${CMAKE:-$(which cmake)}``. Specify a CMake
        version by setting the *CMAKE* variable to the full executable path for
        the CMake version you would like to use. See also :ref:`containers`.

    CMAKE_COMPILER_SCRIPT
        CMake command line options for a tool chain. A definition is provided by
        the mix-in toolchain definitions (e.g. ``.use-gcc8``) to be appended to
        :command:`cmake` calls in a job's *script*.

    CMAKE_MPI_OPTIONS
        Provide CMake command line arguments to define |Gromacs| MPI build options.

    DRY_RUN
        Read-only environment variable used to control behaviour of script uploading
        artifact files to the ftp and web servers. Set to false to actually upload
        files. This is usually done through the pipeline submission script, but can
        be done manual as well through the web interface.                                   

    GROMACS_MAJOR_VERSION
        Read-only environment variable for CI scripts to check the
        library API version to expect from the ``build`` job artifacts.
        Initially, this variable is only defined in
        :file:`admin/gitlab-ci/api-client.matrix/gromacs-main.gitlab-ci.yml`
        but could be moved to :file:`admin/gitlab-ci/global.gitlab-ci.yml` if found
        to be of general utility.

    GROMACS_RELEASE
        Read-only environment variable that can be checked to see if a job is
        executing in a pipeline for preparing a tagged release.
        Can be set when launching pipelines via the GitLab web interface.
        For example, see *rules* mix-ins in :file:`admin/gitlab-ci/global.gitlab-ci.yml`.

    REGRESSIONTESTBRANCH
        Use this branch of the regressiontests rather than main to allow for
        merge requests that require updated regression tests with valid CI tests.

    REGRESSIONTESTCOMMIT
        Use this commit to the regressiontests rather than the head on main to
        allow for merge requests that require updated regression tests with
        valid CI tests.

    POST_MERGE_ACCEPTANCE
        Read-only environment variable that indicates that only jobs scheduled to
        run after a commit has been merged into its target branch should be executed.
        Can be set to run pipelines through the web interface or as schedules.
        For use please see the *rules* mix-ins in :file:`admin/gitlab-ci/global.gitlab-ci.yml`.

    GMX_PIPELINE_SCHEDULE
        Read-only environment variable used exclusively by job rules.
        Rule elements of the form ``if-<value>-then-on-success`` check
        whether ``GMX_PIPELINE_SCHEDULE==value``. Allowed values
        are determined by the rule elements available in
        :file:`admin/gitlab-ci/rules.gitlab-ci.yml`, and include
        ``nightly`` and ``weekly`` to restrict jobs to only run
        in the corresponding schedules.

Setting variables
"""""""""""""""""

Variables for individual piplelines are set in the gitlab interface under 
``CI/CD``; ``Pipelines``. Then chose in the top right corner ``Run Piplelines``.
Under ``Run for``, the desired branch may be selected, and variables may be set
in the fields below.

Using GPUs in Gitlab-runner
"""""""""""""""""""""""""""

Previously, |Gromacs| used a hacked local version of Gitlab-runner where we had
added support for Kubernetes extended resources. However, Gitlab has unfortunately
not shown interest in merging these, and as the runner has evolved it is
difficult to keep up. In the future it might be possible to select GPUs directly
in the job configuration, but for now we use the ability to specify it in each
Gitlab-runner configuration and thus have separate runners going for CPU-only as
well as single or dual GPU devices from Nvidia, AMD, and Intel.

To enable both us and other users to also use the shared Gitlab runners, the
top-level configuration :file:`.gitlab-ci.yml` now contains a few variables where
you can select what tags to use for Gitlab-runners to get single or dual devices
from each vendor. There are also variables that allow you to set the largest
number of devices you have (on single nodes) in these runners; if any tests
cannot be run because you do not have the right hardware, we will simply skip
those tests.

.. _containers:

Containers
----------

|Gromacs| project infrastructure uses Docker containerization to
isolate automated tasks.
A number of images are maintained to provide a breadth of testing coverage.

Scripts and configuration files for building images are stored in the repository
under :file:`admin/containers/`.
Images are (re)built manually by |Gromacs| project staff and pushed to
DockerHub and GitLab.
See https://hub.docker.com/u/gromacs and https://gitlab.com/gromacs/gromacs/container_registry

GitLab Container Registry
"""""""""""""""""""""""""

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

.. _gitlab-ci tools:

Tools
-----
(in :file:`admin/`)

.. autoprogram:: make-release-build:parser
    :prog: make-release-build.py

.. _trigger-post-merge:

.. autoprogram:: trigger-post-merge:parser
    :prog: trigger-post-merge.py

admin/containers/buildall.sh
""""""""""""""""""""""""""""

Uses NVidia's
`HPC Container Maker <https://github.com/NVIDIA/hpc-container-maker/tree/master/docs>`__
to generate DockerFiles using our :py:mod:`scripted_gmx_docker_builds` module.
Refer to the contents of :file:`admin/buildall.sh` for the flags currently in use.
Run the script to see the tagged images currently being produced.

scripted_gmx_docker_builds.py
"""""""""""""""""""""""""""""
(in :file:`admin/containers/`)

.. argparse::
    :module: scripted_gmx_docker_builds
    :func: parser
    :prog: scripted_gmx_docker_builds.py
    :nodefault:

Supporting modules in :file:`admin/containers`
----------------------------------------------

:file:`scripted_gmx_docker_builds.py`
"""""""""""""""""""""""""""""""""""""

.. automodule:: scripted_gmx_docker_builds
    :members:

:file:`utility.py`
""""""""""""""""""

.. automodule:: utility
    :members:
