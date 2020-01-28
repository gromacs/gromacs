==============================
Full installation instructions
==============================

.. highlight:: bash

Installation instructions for the :py:mod:`gmxapi` Python package,
built on GROMACS.

Command line examples assume the `bash <https://www.gnu.org/software/bash/>`_ shell.

.. note:: Regarding multiple GROMACS installations:
    Many GROMACS users switch between multiple GROMACS installations on the same
    computer using an HPC module system and/or a GMXRC configuration script.
    For the equivalent sort of environment switching with the :py:mod:`gmxapi` Python package,
    we recommend installing it in a different
    `Python virtual environment <https://www.google.com/search?q=python+virtual+environment>`_
    for each GROMACS installation.
    Once built, a particular copy of the :py:mod:`gmxapi` Python package always refers to the
    same GROMACS installation.

.. contents:: Contents
    :local:
    :depth: 2

.. note::

    The following documentation contains frequent references to the :command:`pip` tool
    for installing Python packages. In some cases, an unprivileged user should
    use the ``--user`` command line flag to tell :py:mod:`pip` to install packages
    into the user site-packages directory rather than the default site-packages
    directory for the Python installation. This flag is not appropriate when
    running :command:`pip` in a virtual environment (as recommended) and is omitted in
    this documentation. If you need the ``--user`` flag, you should modify the
    example commands to look something like :command:`pip install --upgrade somepackage --user`

.. note::

    These instructions use the executable names :command:`python` and :command:`pip`
    instead of :command:`python3` or :command:`pip3`. Some Python installations require the ``3``
    suffix, but it is usually not necessary if you have already activated a Python
    virtual environment (recommended).

Overview
========

Typically, setting up the *gmxapi* Python package follows these three steps.
If this overview is sufficient for your computing environment,
you may disregard the rest of this document.

Install GROMACS
---------------

Locate your GROMACS installation, or build and install GROMACS 2020 or higher.

.. seealso:: `GROMACS installation <http://manual.gromacs.org/documentation/current/install-guide/index.html>`_

The following assumes GROMACS is installed to :file:`/path/to/gromacs`

Set up a Python virtual environment
-----------------------------------

::

    python3 -m venv $HOME/myvenv
    . $HOME/myvenv/bin/activate
    python -m ensurepip --default-pip
    pip install --upgrade pip setuptools
    pip install --upgrade cmake scikit-build

.. seealso:: :ref:`gmxapi venv`

Install the gmxapi Python package
---------------------------------

::

    . /path/to/gromacs/bin/GMXRC
    pip install gmxapi

.. seealso:: :ref:`installation`

Background
==========

*gmxapi* comes in three parts:

* GROMACS gmxapi library for C++.
* This Python package, supporting Python 3.5 and higher
* MD restraint plugins and sample gmxapi client code

GROMACS requirements
--------------------

The Python package requires a GROMACS installation.
Locate an existing GROMACS installation, or
`build and install GROMACS <http://manual.gromacs.org/documentation/current/install-guide/index.html>`_
before proceeding.

.. note::

    Note that gmxapi requires that GROMACS is configured with ``GMXAPI=ON`` and ``BUILD_SHARED_LIBS=ON``.
    These are enabled by default in most cases. If these options were overridden
    for your GROMACS installation, you will see CMake errors when trying to build
    and install the gmxapi Python package or other client software.

Then, "source" the :file:`GMXRC` file from the GROMACS installation as you normally would
before using GROMACS, or note its installation location so that you can pass it
to the build configuration.

Build system requirements
-------------------------

gmxapi can be built for Python 3.5 and higher.

You will need a C++ 14 compatible compiler and a reasonably up-to-date version
of CMake.
Full gmxapi functionality may also require an MPI compiler (e.g. :command:`mpicc`).

Important: To build a module that can be imported by Python, you need a Python
installation that includes the Python headers. Unfortunately, it is not always
obvious whether these headers are present or where to find them. The simplest
answer is to just try to build the Python package using these instructions, and
if gmxapi is unable to find the Python tools it needs, try a different Python
installation or install the additional development packages.

On a Linux system, this may require installing packages such as ``python-dev``
and/or ``python3-dev``.
If you are building Python, either from scratch or with a tool like
:command:`pyenv install` (see
`wiki entry <https://github.com/pyenv/pyenv/wiki#how-to-build-cpython-with---enable-shared>`_
),
be sure to enable installation of the Python C library with the
``--enable-shared`` flag.
Alternatively, various Python distributions provide a
sufficient build environment while only requiring installation into a user
home directory. (Some examples below.)

If you are using an HPC system with software available through modules you may
be able to just :command:`module load` a different Python installation and find one
that works.

Python environment requirements
-------------------------------

gmxapi requires Python 3.5 or higher. Check your version with
:command:`python3 --version` or :command:`python --version`.

..  note::

    The following documentation assumes you do not need to use a trailing '3' to
    access a Python 3 interpreter on your system.
    The default Python interpreter on your system may use :command:`python3` and :command:`pip3`
    instead of :command:`python` and :command:`pip`. You can check the version with
    :command:`python3 --version` or :command:`python --version` and :command:`pip --version`.

To build and install, you also need the packages :py:mod:`cmake`,
:py:mod:`setuptools`, :py:mod:`networkx`, and :py:mod:`scikit-build`.

For full functionality, you should also have :py:mod:`mpi4py` and :py:mod:`numpy`.
These requirements and version numbers are listed in :file:`requirements.txt`.

The easiest way to make sure you have the requirements installed, first update
:py:mod:`pip`, then use the :file:`requirements.txt` file provided with the repository.
File paths in this section are relative to the root directory of your local copy
of the GROMACS source.

Confirm that :py:mod:`pip` is available, install :py:mod:`pip` if it is missing, or get
instructions on how to install :py:mod:`pip`::

    python -m ensurepip --default-pip

Install or upgrade required components::

    python -m pip install --upgrade pip
    pip install --upgrade setuptools

"requirements" files in GROMACS source tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are building from source code in a local copy of the GROMACS source
repository, some helpful files allow you to preinstall the Python requirements
before installing the :py:mod:`gmxapi` package.

    pip install -r python_packaging/src/requirements.txt

If building documentation or running tests,
:command:`pip install -r python_packaging/requirements-docs.txt` or
:command:`pip install -r python_packaging/requirements-test.txt`,
respectively, or see below.

Documentation build requirements
--------------------------------

See :ref:`gmxapi_package_documentation`

.. _testing_requirements:

Testing requirements
--------------------

Testing is performed with `pytest <https://docs.pytest.org/en/latest/>`_.
Tests also require :py:mod:`numpy`.
You can probably install both with :command:`pip`::

    pip install pytest numpy

To test the full functionality also requires an MPI parallel environment.
You will need the :py:mod:`mpi4py` Python package and an MPI launcher
(such as :command:`mpiexec`, :command:`mpirun`, a launcher provided by your HPC queuing system,
or whatever is provided by your favorite MPI package for your operating system).

.. _mpi_requirements:

MPI requirements
----------------

For the ensemble simulations features, you will need an MPI installation.
On an HPC system, this means you will probably have to use :command:`module load`
to load a compatible set of MPI tools and compilers.
Check your HPC documentation or try :command:`module avail` to look for an
``openmpi``, ``mpich``, or ``mvapich`` module and matching compiler module.
This may be as simple as::

    module load gcc
    module load mpicc

Note that the compilers loaded might not be the first compilers discovered
automatically by the build tools we will use below,
so you may have to specify compilers on the command line for consistency.
It may be necessary to require that GROMACS, gmxapi,
and the sample code are built with the same compiler(s).

Note that strange errors have been known to occur when :py:mod:`mpi4py` is built with
different a different tool set than has been used to build Python and gmxapi.
If the default compilers on your system are not sufficient for GROMACS or gmxapi,
you may need to build, e.g., OpenMPI or MPICH, and/or build :py:mod:`mpi4py` with a
specific MPI compiler wrapper. This can complicate building in environments such
as Conda_.

Set the MPICC environment variable to the MPI compiler wrapper and forcibly
reinstall :py:mod:`mpi4py`::

    export MPICC=`which mpicc`
    pip install --no-cache-dir --upgrade --no-binary \":all:\" --force-reinstall mpi4py

If you have a different MPI C compiler wrapper, substitute it for :command:`mpicc` above.

.. _installation:

Installing the Python package
=============================

We recommend using Python's `pip <https://pip.pypa.io/en/stable/>`_
package installer to automatically download, build, and install the latest
version of the gmxapi package into a Python
`virtual environment <https://docs.python.org/3/tutorial/venv.html>`_,
though it is also possible to install without a virtual environment.
If installing without a virtual environment as an un-privileged user,
you may need to set the CMake variable ``GMXAPI_USER_INSTALL``
(``-DGMXAPI_USER_INSTALL=ON`` on the :command:`cmake` command line)
and / or use the ``--user`` option with :command:`pip install`.

Recommended installation
------------------------

The instructions in this section assume that *pip* is able to download files
from the internet. Alternatively, refer to :ref:`gmxapi offline install`.

Locate or install GROMACS
^^^^^^^^^^^^^^^^^^^^^^^^^

You need a GROMACS installation that includes the gmxapi headers and library.
If GROMACS 2020 or higher is already installed,
*and* was configured with ``GMXAPI=ON`` at build time,
you can just source the GMXRC
(so that the Python package knows where to find GROMACS)
and skip to the next section.

Otherwise, install a supported version of GROMACS.
When building GROMACS from source, be sure to configure cmake with the flag
``-DGMXAPI=ON`` (default).

Set the environment variables for the GROMACS installation so that the gmxapi
headers and library can be found when building the Python package.
If you installed to a :file:`gromacs-gmxapi` directory in your home directory as
above and you use the :command:`bash` shell, do::

    source $HOME/gromacs-gmxapi/bin/GMXRC

.. _gmxapi venv:

Set up a Python virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We recommend installing the Python package in a virtual environment.
If not installing in a virtual environment, you may not be able to install
necessary prerequisites (e.g. if you are not an administrator of the system you are on).

The following instructions use the :py:mod:`venv` module.
Alternative virtual environments, such as Conda_,
should work fine, but are beyond the scope of this document.
(We welcome contributed recipes!)

Depending on your computing environment, the Python 3 interpreter may be accessed
with the command :command:`python` or :command:`python3`. Use :command:`python --version` and
:command:`python3 --version` to figure out which you need to use. The following assumes
the Python 3 interpreter is accessed with :command:`python3`.

Create a Python 3 virtual environment::

    python3 -m venv $HOME/myvenv

Activate the virtual environment. Your shell prompt will probably be updated with the name of the environment you
created to make it more obvious.

.. code-block:: none

    $ source $HOME/myvenv/bin/activate
    (myvenv)$

..  note::

    After activating the *venv*, :command:`python` and :command:`pip` are sufficient.
    (The '3' suffix will no longer be necessary and will be omitted in the rest
    of this document.)

Activating the virtual environment may change your shell prompt to indicate the
environment is active. The prompt is omitted from the remaining examples, but
the remaining examples assume the virtual environment is still active.
(Don't do it now, but you can deactivate the environment by running :command:`deactivate`.)

Install dependencies
^^^^^^^^^^^^^^^^^^^^

It is always a good idea to update :py:mod:`pip` and :py:mod:`setuptools` before installing
new Python packages::

    pip install --upgrade pip setuptools

The gmxapi installer requires a few additional packages. It is best to make sure
they are installed and up to date before proceeding.

::

    pip install --upgrade cmake scikit-build

For MPI, we use :py:mod:`mpi4py`.
Make sure it is using the same MPI installation that we are building
GROMACS against and building with compatible compilers.

::

    python -m pip install --upgrade pip setuptools
    MPICC=`which mpicc` pip install --upgrade mpi4py

.. seealso:: :ref:`mpi_requirements`

Install the latest version of gmxapi
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fetch and install the latest version of gmxapi from the Python Packaging Index::

    pip install gmxapi

If :command:`pip` does not find your GROMACS installation, use one of the following
environment variables to provide a hint.

gmxapi_DIR
~~~~~~~~~~

If you have a single GROMACS installation at :file:`/path/to/gromacs`, it is usually
sufficient to provide this location to :command:`pip` through the :envvar:`gmxapi_DIR`
environment variable.

Example::

    gmxapi_DIR=/path/to/gromacs pip install gmxapi

GMXTOOLCHAINDIR
~~~~~~~~~~~~~~~

If you have multiple builds of GROMACS distinguished by suffixes
(e.g. *_d*, *_mpi*, etcetera), or if you need to provide extra hints to :command:`pip`
about the software tools that were used to build GROMACS, you can specify a
directory in which the installer can find a CMake "tool chain".

In the following example, ``${SUFFIX}`` is the suffix that distinguishes the
particular build of GROMACS you want to target (refer to GROMACS installation
instructions for more information.) ``${SUFFIX}`` may simply be empty, or ``''``.

::

    GMXTOOLCHAINDIR=/path/to/gromacs/share/cmake/gromacs${SUFFIX} pip install gmxapi

Install from source
-------------------

You can also install the :py:mod:`gmxapi` Python package from within a local copy of
the GROMACS source repository. Assuming you have already obtained the GROMACS
source code and you are in the root directory of the source tree, you will find
the :py:mod`gmxapi` Python package sources in the :file:`python_packaging/src` directory.

::

    cd python_packaging/src
    pip install -r requirements.txt
    pip install .

.. _gmxapi offline install:

Offline install
---------------

If the required dependencies are already installed, you can do a quick installation
without internet access, either from the source directory or from a source archive.

For example, the last line of the previous example could be replaced with::

    pip install --no-cache-dir --no-deps --no-index --no-build-isolation .

Refer to :py:mod:`pip` documentation for descriptions of these options.

If you have built or downloaded a source distribution archive, you can provide
the archive file to :command:`pip` instead of the ``.`` argument::

    pip install gmxapi-0.1.0.tar.gz

In this example, the archive file name is as was downloaded from
`PyPI <https://pypi.org/project/gmxapi/#history>`_ or as built locally,
according to the following instructions.

Building a source archive
-------------------------

A source archive for the gmxapi python package can be built from the GROMACS
source repository using Python ``setuptools`` and ``scikit-build``.

Example::

    pip install --upgrade setuptools scikit-build
    cd python_packaging/src
    python setup.py sdist

This command will create a ``dist`` directory containing a source distribution
archive file. The file name has the form *gmxapi-<version>.<suffix>*, where
*<version>* is the version from the ``setup.py`` file, and *<suffix>* is
determined by the local environment or by additional arguments to ``setup.py``.

.. seealso::

    Python documentation for
    `creating a source distribution
    <https://docs.python.org/3/distutils/sourcedist.html#creating-a-source-distribution>`_

Package maintainers may update the online respository by uploading a freshly
built ``sdist`` with ``python -m twine upload dist/*``

.. _gmxapi_package_documentation:

Accessing gmxapi documentation
==============================

Documentation for the Python classes and functions in the gmx module can
be accessed in the usual ways, using ``pydoc`` from the command line or
``help()`` in an interactive Python session.

The complete documentation (which you are currently reading)
can be browsed `online <http://manual.gromacs.org/current/gmxapi/>`__
or built from a copy of the GROMACS source repository.

Documentation is built from a combination of Python module documentation and
static content, and requires a local copy of the GROMACS source repository.

Build with GROMACS
------------------

To build the full gmxapi documentation with GROMACS, configure GROMACS with
``-DGMX_PYTHON_PACKAGE=ON`` and build the GROMACS documentation normally.
This will first build the *gmxapi* Python package and install it to a temporary
location in the build tree. Sphinx can then import the package to automatically
extract Python docstrings.

Sometimes the build environment can choose a different Python interpreter than
the one you intended.
You can set the ``PYTHON_EXECUTABLE`` CMake variable to explicitly choose the
Python interpreter for your chosen installation.
For example: ``-DPYTHON_EXECUTABLE=\`which python\```

Docker web server
-----------------

Alternatively, build the ``docs`` Docker image from ``python_packaging/docker/docs.dockerfile``
or pull a prebuilt image from DockerHub. Refer to the dockerfile or to
https://hub.docker.com/r/gmxapi/docs for more information.

.. todo::

    Document sample_restraint package. Reference issue
    `3027 <https://redmine.gromacs.org/issues/3027>`_

.. _gmxapi install troubleshooting:

Troubleshooting
===============

Couldn't find the ``gmxapi`` support library?
If you don't want to "source" your ``GMXRC`` file, you
can tell the package where to find a gmxapi compatible GROMACS installation with
``gmxapi_DIR``. E.g. ``gmxapi_DIR=/path/to/gromacs pip install .``

Before updating the ``gmxapi`` package it is generally a good idea to remove the
previous installation and to start with a fresh build directory. You should be
able to just ``pip uninstall gmxapi``.

Do you see something like the following?

.. code-block:: none

   CMake Error at gmx/core/CMakeLists.txt:45 (find_package):
      Could not find a package configuration file provided by "gmxapi" with any
      of the following names:

        gmxapiConfig.cmake
        gmxapi-config.cmake

      Add the installation prefix of "gmxapi" to CMAKE_PREFIX_PATH or set
      "gmxapi_DIR" to a directory containing one of the above files.  If "gmxapi"
      provides a separate development package or SDK, be sure it has been
      installed.

This could be because

* GROMACS is not already installed
* GROMACS was built without the CMake variable ``GMXAPI=ON``
* or if ``gmxapi_DIR`` (or ``GROMACS_DIR``) is not a path containing directories
  like ``bin`` and ``share``.

If you are not a system administrator you are encouraged to install in a Python
virtual environment, created with virtualenv or Conda_.
Otherwise, you will need to specify the ``--user`` flag to ``pip``.

Two of the easiest problems to run into are incompatible compilers and
incompatible Python. Try to make sure that you use the same C and C++
compilers for GROMACS, for the Python package, and for the sample
plugin. These compilers should also correspond to the :command:`mpicc` compiler
wrapper used to compile :py:mod:`mpi4py`. In order to build the Python
package, you will need the Python headers or development installation,
which might not already be installed on the machine you are using. (If
not, then you will get an error about missing :file:`Python.h` at some
point.) If you have multiple Python installations (or modules available
on an HPC system), you could try one of the other Python installations,
or you or a system administrator could install an appropriate Python dev
package. Alternatively, you might try installing your own Anaconda or
MiniConda in your home directory.

If an attempted installation fails with CMake errors about missing
“gmxapi”, make sure that Gromacs is installed and can be found during
installation. For instance,

::

    gmxapi_DIR=/Users/eric/gromacs python setup.py install --verbose

Pip and related Python package management tools can be a little too
flexible and ambiguous sometimes. If things get really messed up, try
explicitly uninstalling the :py:mod:`gmxapi` module and its dependencies, then do
it again and repeat until :command:`pip` can no longer find any version of any
of the packages.

::

    pip uninstall gmxapi
    pip uninstall cmake
    # ...

Successfully running the test suite is not essential to having a working
:py:mod:`gmxapi` package. We are working to make the testing more robust, but
right now the test suite is a bit delicate and may not work right, even
though you have a successfully built the :py:mod:`gmxapi` package. If you want to
troubleshoot, though, the main problems seem to be that automatic
installation of required python packages may not work (requiring manual
installations, such as with :command:`pip install somepackage`) and ambiguities
between python versions. 

If you are working in a development branch of the repository, note that
the upstream branch may be reset to ``master`` after a new release is
tagged. In general, but particularly on the ``devel`` branch, when you
do a :command:`git pull`, you should use the ``--rebase`` flag.

If you fetch this repository and then see a git status like this::

    $ git status
    On branch devel
    Your branch and 'origin/devel' have diverged,
    and have 31 and 29 different commits each, respectively.

then :py:mod:`gmxapi` has probably entered a new development cycle. You can
do :command:`git pull --rebase` to update to the latest development branch.

If you do a :command:`git pull` while in ``devel`` and get a bunch of unexpected
merge conflicts, do :command:`git merge --abort; git pull --rebase` and you should
be back on track.

If you are developing code for gmxapi, this should be an indication to
rebase your feature branches for the new development cycle.

.. _Conda: https://docs.conda.io/en/latest/
