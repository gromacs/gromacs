==============================
Full installation instructions
==============================

.. highlight:: bash

Installation instructions for the :py:mod:`gmxapi` Python package,
built on |Gromacs|.

Command line examples assume the `bash <https://www.gnu.org/software/bash/>`_ shell.

.. admonition:: Regarding multiple |Gromacs| installations
    :class: note

    Many |Gromacs| users switch between multiple |Gromacs| installations on the same
    computer using an HPC module system and/or a :ref:`GMXRC <getting access to |Gromacs|>` configuration script.
    For the equivalent sort of environment switching with the :py:mod:`gmxapi` Python package,
    we recommend installing it in a different
    `Python virtual environment <https://www.google.com/search?q=python+virtual+environment>`_
    for each |Gromacs| installation.
    Once built, a particular copy of the :py:mod:`gmxapi` Python package always refers to the
    same |Gromacs| installation.

.. admonition:: Unprivileged :command:`pip install`
    :class: tip

    The following documentation contains frequent references to the pip_ tool
    for installing Python packages. In some cases, an unprivileged user should
    use the ``--user`` command line flag to tell pip_ to install packages
    into the user site-packages directory rather than the default site-packages
    directory for the Python installation. This flag is not appropriate when
    running :command:`pip` in a virtual environment (as recommended) and is omitted in
    this documentation. If you need the ``--user`` flag, you should modify the
    example commands to look something like :command:`pip install --upgrade somepackage --user`

.. admonition:: Python 3 executable names
    :class: info

    These instructions use the executable names :command:`python` and :command:`pip`
    instead of :command:`python3` or :command:`pip3`. Some Python installations require the ``3``
    suffix, but it is usually not necessary if you have already activated a Python
    virtual environment (recommended).

Overview
========

Typically, setting up the *gmxapi* Python package follows these three steps.
If this overview is sufficient for your computing environment,
you may disregard the rest of this document.

Install |Gromacs|
-----------------

Locate your |Gromacs| installation, or build and install.
|Gromacs| 2022 or higher is recommended.

.. seealso:: `GROMACS installation <http://manual.gromacs.org/documentation/current/install-guide/index.html>`_

The following assumes |Gromacs| is installed to :file:`/path/to/gromacs`

Set up a Python virtual environment
-----------------------------------

.. seealso:: :ref:`gmxapi venv`

.. note:: :py:mod:`mpi4py` may require additional arguments (compiler hints).
    See :ref:`mpi_requirements`

::

    python3 -m venv $HOME/myvenv
    . $HOME/myvenv/bin/activate
    python -m ensurepip --default-pip
    pip install --upgrade pip setuptools wheel
    pip install mpi4py

Install the gmxapi Python package
---------------------------------

::

    . /path/to/gromacs/bin/GMXRC
    pip install --no-cache-dir gmxapi

.. seealso:: :ref:`installation`

Background
==========

*gmxapi* comes in three parts:

* |Gromacs| gmxapi library for C++.
* This Python package, supporting Python 3.7 and higher
* MD restraint plugins and sample gmxapi client code

|Gromacs| requirements
----------------------

The Python package requires a |Gromacs| installation.
Locate an existing |Gromacs| installation, or
`build and install GROMACS <http://manual.gromacs.org/documentation/current/install-guide/index.html>`_
before proceeding.

.. note::

    Note that gmxapi requires that |Gromacs| is configured with ``GMXAPI=ON`` and ``BUILD_SHARED_LIBS=ON``.
    These are enabled by default in most cases. If these options were overridden
    for your |Gromacs| installation, you will see CMake errors when trying to build
    and install the gmxapi Python package or other client software.

If your installation has a :file:`GMXRC` file, "source" the file
:ref:`as you normally would <getting access to |Gromacs|>` before using |Gromacs|.
Otherwise, note the installation location so that you can provide it when
building the gmxapi package.

Build system requirements
-------------------------

gmxapi can be built for Python 3.7 and higher.

You will need a C++ 17 compatible compiler and a reasonably up-to-date version
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

gmxapi requires Python 3.7 or higher. Check your version with
:command:`python3 --version` or :command:`python --version`.

..  note::

    The following documentation assumes you do not need to use a trailing '3' to
    access a Python 3 interpreter on your system.
    The default Python interpreter on your system may use :command:`python3` and :command:`pip3`
    instead of :command:`python` and :command:`pip`. You can check the version with
    :command:`python3 --version` or :command:`python --version` and :command:`pip --version`.

To build and install, you need the Python packages for
cmake_, networkx_, and setuptools_
(all available from `PyPI with pip <https://pip.pypa.io/en/stable/>`_).

For full functionality, you should also have mpi4py_ and numpy_.
These requirements and version numbers are listed in :file:`requirements.txt`.

The easiest way to make sure you have the requirements installed, first update
pip_, then use the :file:`requirements.txt` file provided with the repository.
File paths in this section are relative to the root directory of your local copy
of the |Gromacs| source.

Confirm that pip_ is available, install pip_ if it is missing, or get
instructions on how to install pip_::

    python -m ensurepip --default-pip

Install or upgrade required components::

    python -m pip install --upgrade pip
    pip install --upgrade setuptools wheel

"requirements" files in |Gromacs| source tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are building from source code in a local copy of the |Gromacs| source
repository, a :file:`requirements.txt` allows you to preinstall the Python
requirements before installing the :py:mod:`gmxapi` package.

    pip install -r python_packaging/gmxapi/requirements.txt

Documentation build requirements
--------------------------------

See :ref:`gmxapi_package_documentation`

.. _testing requirements:

Testing requirements
--------------------

Note that the test suite is only available in the |Gromacs| source tree.
(It is not part of the installed package.)
Acquire the |Gromacs| sources with :command:`git` or by downloading an archive, as documented elsewhere.

Testing is performed with `pytest <https://docs.pytest.org/en/latest/>`_.

:file:`python_packaging/gmxapi/requirements.txt` lists additional requirements for testing.
With pip_::

    pip install -r python_packaging/gmxapi/requirements.txt

To test the full functionality also requires an MPI parallel environment.
You will need the mpi4py_ Python package and an MPI launcher
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

If you are using a |Gromacs| installation that is already available through
``module load``, try to find a Python installation with the ``mpi4py`` package
that is also available through ``module load``. The *module* system will
generally enforce toolchain compatibility between the loaded modules. If you
``module load`` mpi4py or a Python installation with mpi4py, you will probably
want to use this version of the package in your venv. (See :ref:`gmxapi venv`)
If you ``module load`` an MPI-enabled |Gromacs| installation, ``gmxapi`` will
try to check ``mpi4py`` for compatibility.

Note that the compilers loaded might not be the first compilers discovered
automatically by the build tools we will use below,
so you may have to specify compilers on the command line for consistency.
It may be necessary to require that |Gromacs|, gmxapi,
and the sample code are built with the same compiler(s).

Note that strange errors have been known to occur when mpi4py_ is built with
a different tool set than has been used to build Python and gmxapi.
If the default compilers on your system are not sufficient for |Gromacs| or gmxapi,
you may need to build, e.g., OpenMPI or MPICH, and/or
`build mpi4py <https://mpi4py.readthedocs.io/en/stable/install.html>`__ with a
specific MPI compiler wrapper. This can complicate building in environments such
as Conda_. You should be able to confirm that your MPI compiler wrapper is consistent
with your |Gromacs| tool chain by comparing the output of :command:`mpicc --version`
with the compiler information reported by :command:`gmx --version`.

Set the ``MPICC`` environment variable to the MPI compiler wrapper and forcibly
reinstall mpi4py_::

    export MPICC=`which mpicc`
    pip install --no-cache-dir --upgrade --no-binary ":all:" --force-reinstall mpi4py

If you have a different MPI C compiler wrapper, substitute it for :command:`mpicc` above.

While ``gmxapi`` is configuring its build system during installation, it will
try to confirm the compatibility of the ``mpi4py`` toolchain with that of the
|Gromacs| installation. If they appear incompatible, you should see a ``CMake``
message that includes a guess at what you might try using for ``MPICC``.
(If using ``pip``, consider using the ``--verbose`` option for more build output.)

.. _installation:

Installing the Python package
=============================

We recommend using Python's `pip <https://pip.pypa.io/en/stable/>`_
package installer to automatically download, build, and install the latest
version of the gmxapi package into a Python
`virtual environment <https://docs.python.org/3/tutorial/venv.html>`_,
though it is also possible to install without a virtual environment.
If installing without a virtual environment as an un-privileged user,
you may need to use the ``--user`` option with :command:`pip install`.

Recommended installation
------------------------

The instructions in this section assume that *pip* is able to download files
from the internet. Alternatively, refer to :ref:`gmxapi offline install`.

Locate or install |Gromacs|
^^^^^^^^^^^^^^^^^^^^^^^^^^^

You need a |Gromacs| installation that includes the gmxapi headers and library.

.. warning:: gmxapi does not recognize multiple |Gromacs| installations to the same ``CMAKE_INSTALL_PREFIX``.

    The Python package uses files installed to ``.../share/cmake/gmxapi/`` to configure its C++
    component. These configuration files are overwritten when installing |Gromacs| to the same
    `CMAKE_INSTALL_PREFIX <https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html>`__.
    Overlapping |Gromacs| installations may occur when |Gromacs| is installed for multiple
    configurations of MPI support and floating point precision.
    (See :issue:`4334` and related issues.)

If |Gromacs| 2020 or higher is already installed,
*and* was configured with ``GMXAPI=ON`` at build time (the default),
you can just source the :ref:`GMXRC <getting access to |Gromacs|>`
(so that the Python package knows where to find |Gromacs|)
and skip to the next section.

Otherwise, install a supported version of |Gromacs|.
When building |Gromacs| from source, be sure to configure cmake with the flag
``-DGMXAPI=ON`` (default).

Set the environment variables for the |Gromacs| installation so that the gmxapi
headers and library can be found when building the Python package.
If you installed to a :file:`gromacs-gmxapi` directory in your home directory as
above and you use the :command:`bash` shell, do::

    source $HOME/gromacs-gmxapi/bin/GMXRC

If you are using a |Gromacs| installation that does not provide ``GMXRC``, see
`gmxapi cmake hints`_ and additional CMake hints below.

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

.. _system-site-packages:

.. admonition:: --system-site-packages
    :class: tip

    It can be tricky to properly or optimally build MPI enabled software in
    computing clusters, and administrators often provide prebuilt packages like
    :py:mod:`mpi4py`. If your computing environment has multiple Python installations,
    try to choose one that already includes ``mpi4py``. When you are using a
    Python installation that provides ``mpi4py``, generally, you should be sure
    to use the existing ``mpi4py`` installation in your new virtual environment
    by creating the ``venv`` with the ``--system-site-packages`` option.

    In personal computing environments (laptops and workstations), it is common to
    have multiple Python installations, and it can be hard to keep packages in the
    different installations from conflicting with each other. Unless you know that
    you want to inherit the ``mpi4py`` package from the system installation, it is
    generally cleaner *not* to inherit the system site-packages.

Create a Python 3 virtual environment::

    python3 -m venv $HOME/myvenv

*or* (see note)::

    python3 -m venv --system-site-packages $HOME/myvenv

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

It is always a good idea to update pip_, setuptools_, and wheel_ before installing
new Python packages::

    pip install --upgrade pip setuptools wheel

The gmxapi installer requires a few additional packages. It is best to make sure
they are installed and up to date before proceeding.

::

    pip install --upgrade cmake pybind11

We use mpi4py_ for some features and to ensure compatible MPI bindings
throughout your Python environment.
**If you did not inherit mpi4py from system site-packages**
(see :ref:`above <system-site-packages>`),
make sure to
`install it <https://mpi4py.readthedocs.io/en/stable/install.html>`__
using the same MPI installation that we are building
|Gromacs| against, and build with compatible compilers.

::

    MPICC=`which mpicc` pip install --no-cache-dir --upgrade mpi4py

.. seealso:: :ref:`mpi_requirements`

Install the latest version of gmxapi
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fetch and install the latest official version of gmxapi from the Python Packaging Index.
Avoid locally cached previously-built packages that may be incompatible with
your current environment or |Gromacs| installation::

    # Get the latest official release.
    pip install --no-cache-dir gmxapi

or::

    pip download gmxapi
    pip install gmxapi-<version>.tar.gz

substituting the name of the downloaded source distribution archive.

.. admonition:: Avoid cached "wheel" packages.
    :class: warning

    ``pip`` downloads a source distribution archive for gmxapi, then builds a
    "wheel" package for your |Gromacs| installation.
    This "wheel" normally gets cached, and will be used by any later attempt to
    ``pip install gmxapi`` instead of rebuilding. This is not what you want,
    if you upgrade |Gromacs| or if you want to install the Python package for a
    different |Gromacs| configuration (e.g. double-precision or different MPI option.)

    You can use ``--no-cache-dir`` to force rebuild of the package and its
    build dependencies. This may be slow, however, and you may want to use
    cached dependencies. You can
    `avoid wheel cache <https://pip.pypa.io/en/stable/topics/caching/#avoiding-caching>`__
    for just one target package by installing from the filesystem
    instead of directly from PyPI.

    See also :issue:`4335`

The `PyPI repository <https://pypi.org/project/gmxapi/#history>`_
may include pre-release versions,
but :command:`pip` will ignore them unless you use the ``--pre`` flag::

    # Get the latest version, including pre-release versions.
    pip install --no-cache-dir --pre gmxapi

If :command:`pip` does not find your |Gromacs| installation, use one of the following
environment variables to provide a hint.

The installer will also look for a ``CMAKE_ARGS`` environment variable. If found,
The ``$CMAKE_ARGS`` string will be split into additional arguments that will be
provided to CMake when building the *gmxapi* package.

.. _gmxapi cmake hints:

gmxapi_ROOT
~~~~~~~~~~~

If you have a single |Gromacs| installation at :file:`/path/to/gromacs`, it is usually
sufficient to provide this location to :command:`pip` through the :envvar:`gmxapi_ROOT`
environment variable.

Example::

    gmxapi_ROOT=/path/to/gromacs pip install --no-cache-dir gmxapi

Note that this is equivalent to providing the CMake variable definition::

    CMAKE_ARGS="-Dgmxapi_ROOT=/path/to/gromacs" pip install --no-cache-dir gmxapi

|Gromacs| CMake hints
~~~~~~~~~~~~~~~~~~~~~

If you have multiple builds of |Gromacs| distinguished by suffixes
(e.g. *_d*, *_mpi*, etcetera), or if you need to provide extra hints to :command:`pip`
about the software tools that were used to build |Gromacs|, you can specify a
CMake "hints" file by including a ``-C <initial-cache>`` option with your ``CMAKE_ARGS``.
(For more information, read about the ``-C``
`command line option <https://cmake.org/cmake/help/latest/manual/cmake.1.html#options>`__
for CMake.)

In the following example, ``${UNIQUE_PREFIX}`` is the path to the directory that holds the
|Gromacs| ``bin``, ``lib``, ``share`` directories, *etc*.
It is *unique* because |Gromacs| provides CMake support for only one build configuration at a time
through ``.../share/cmake/gmxapi/``, even if there are multiple library configurations installed to
the same location. See :issue:`4334`.

``${SUFFIX}`` is the suffix that distinguishes the
particular build of |Gromacs| you want to target (refer to |Gromacs| installation
instructions for more information.) ``${SUFFIX}`` may simply be empty, or ``''``.

You can export ``CMAKE_ARGS`` in your environment, or just provide it at the beginning
of the ``pip install`` command line::

    CMAKE_ARGS="-Dgmxapi_ROOT=${UNIQUE_PREFIX} -C ${UNIQUE_PREFIX}/share/cmake/gromacs${SUFFIX}/gromacs-hints${SUFFIX}.cmake" \
        pip install --no-cache-dir gmxapi

Install from source
-------------------

You can also install the :py:mod:`gmxapi` Python package from within a local copy of
the |Gromacs| source repository. Assuming you have already obtained the |Gromacs|
source code and you are in the root directory of the source tree, you will find
the :py:mod:`gmxapi` Python package sources in the :file:`python_packaging/gmxapi` directory.

::

    cd python_packaging/gmxapi
    pip install -r requirements.txt
    pip install .

.. _gmxapi offline install:

Offline install
---------------

.. admonition:: Recommended, first:
    :class: tip

    :command:`pip install --upgrade build pip setuptools wheel`

You can use :command:`python -m build --skip-dependency-check` to build a binary
distribution archive (from the source distribution) for just the *gmxapi* package,
but then you will have to manually satisfy (separate) dependencies in both the
build and installation environments.

While you have internet access, you need to get access to the *gmxapi* source
distribution and the package dependencies.
You will also want the ``wheel`` and ``build`` packages in environments where
the package(s) will be built.
Only ``pip`` is necessary once a gmxapi ``wheel`` is built.

The following instructions are paraphrased from
https://pip.pypa.io/en/stable/user_guide/#installing-from-local-packages

To build with internet access and then install without::

    # Remove any locally cached (previously built) wheels.
    pip cache remove gmxapi

    # Download gmxapi and dependencies from pypi.
    pip wheel --wheel-dir DIR gmxapi
    # or, using package source from the GROMACS repository
    cd python_packaging/gmxapi
    pip wheel --wheel-dir DIR .

    # Later, install.
    pip install --no-index --find-links=DIR DIR/gmxapi*whl

To download packages and dependencies for later build and installation::

    # if in the GROMACS source repository
    cd python_packaging/gmxapi
    # or download and expand the archive
    pip download --destination-directory DIR gmxapi
    tar xf DIR/gmxapi*
    cd gmxapi*

    # Pre-fetch dependencies to DIR
    pip download --destination-directory DIR .

    # Build and install from the source directory.
    pip install --no-index --find-links=DIR .

Building a source archive
-------------------------

A source archive for the gmxapi python package can be built from the |Gromacs|
source repository using the Python
`build <https://pypa-build.readthedocs.io/en/latest/>`__ module.

Example::

    pip install --upgrade setuptools build
    cd python_packaging/gmxapi
    python -m build --sdist

This command will create a ``dist`` directory containing a source distribution
archive file. The file name has the form
:file:`gmxapi-{version}.{suffix}`, where
*version* is the version from the package metadata, and *suffix* is an
archive file extension determined by the local environment and the current
packaging specifications.

The version information is derived from :py:data:`gmxapi.__version__`
defined by the :py:mod:`gmxapi.version` module.
Pending refinement under :issue:`3851`,
the gmxapi version information is hard coded in the :file:`version.py`.
Make sure you have an up-to-date version of the sources and that the version
information is appropriate before distributing a new release.

.. seealso::

    Python documentation for
    `creating a source distribution
    <https://docs.python.org/3/distutils/sourcedist.html#creating-a-source-distribution>`_

Package maintainers may update the
`online repository <https://pypi.org/project/gmxapi/>`__
by uploading a freshly built ``sdist`` with
``python -m twine upload dist/gmxapi-{version}.{suffix}``.
To update the repository at the PyPI test server, use
``python -m twine upload --repository testpypi dist/gmxapi-{version}.{suffix}``.

.. _gmxapi_package_documentation:

Accessing gmxapi documentation
==============================

Documentation for the Python classes and functions in the gmx module can
be accessed in the usual ways, using ``pydoc`` from the command line or
``help()`` in an interactive Python session.

The complete documentation (which you are currently reading)
can be browsed `online <http://manual.gromacs.org/current/gmxapi/>`__
or built from a copy of the |Gromacs| source repository.

Documentation is built from a combination of Python module documentation and
static content, and requires a local copy of the |Gromacs| source repository.

Build with |Gromacs|
--------------------

To build the full gmxapi documentation with |Gromacs|, configure |Gromacs| with
``-DGMX_PYTHON_PACKAGE=ON`` and build the |Gromacs| documentation normally.
This will first build the *gmxapi* Python package and install it to a temporary
location in the build tree. Sphinx can then import the package to automatically
extract Python docstrings.

Note that this is an entirely CMake-driven installation and Python dependencies
will not be installed automatically. You can update your Python environment
(before configuring with CMake) using the :file:`requirements.txt` files provided
in the :file:`python_packaging/` directory of the repository. Example::

    pip install -r python_packaging/gmxapi/requirements.txt

Sometimes the build environment can choose a different Python interpreter than
the one you intended.
You can set the ``Python3_ROOT_DIR`` or ``CMAKE_PREFIX_PATH`` CMake variable to
explicitly choose the Python installation or *venv* directory.
See also
`CMake FindPython3 <https://cmake.org/cmake/help/latest/module/FindPython3.html>`__.

If you use pyenv or pyenv-virtualenv to dynamically manage your Python version,
you can help identify a particular version with ``pyenv version-name`` and the
directory with ``pyenv prefix {version}``. For example::

    -DPython3_ROOT_DIR=$(pyenv prefix $(pyenv version-name))

.. todo::

    Document sample_restraint package. Reference :issue:`3027`

Testing
=======

Note `testing requirements`_ above.

After installing the :py:mod:`gmxapi` Python package,
you can run the Python test suite from the |Gromacs| source tree.
Example::

    # Assuming you are in the root directory of the repository:
    pytest python_packaging/gmxapi/test/

Refer to :file:`python_packaging/README.md` for more detailed information.

.. _gmxapi install troubleshooting:

Troubleshooting
===============

ImportError at run time with dynamic linking error
--------------------------------------------------

Symptom: Python fails with a weird ``ImportError`` citing something like ``dlopen``::

    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ImportError: dlopen(/.../gmxapi/_gmxapi.so, 0x0002): Symbol not found:
    __ZN12gmxapicompat11readTprFileERKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
      Referenced from: /.../gmxapi/_gmxapi.so
      Expected in: /path/to/gromacs/lib/libgmxapi_mpi_d.0.3.1.dylib

Inconsistencies in the build and run time environments can cause dynamic linking problems at run time.
This could occur if you reinstall |Gromacs| built with a different compiler,
or if ``pip`` or ``CMake`` somehow get tricked into using the wrong compiler tool chain.

Refer to the `gmxapi cmake hints`_ for notes about compiler toolchains.
Rebuild and reinstall the gmxapi Python package with ``--no-cache-dir``
and provide the ``gromacs-hints.cmake`` file for the |Gromacs| installation
you intend to use.

AttributeError: module 'enum' has no attribute 'IntFlag'
--------------------------------------------------------

If you had older versions of some of the dependencies installed,
you might have picked up a transitive dependency on the ``enum34`` package.
Try::

    pip uninstall -y enum34

and see if that fixes the problem. If not, try a fresh virtual environment
(see above) to help narrow down the problem before you
`open an issue <https://gitlab.com/gromacs/gromacs/-/issues/>`_.

Errors regarding pybind11
-------------------------

An error may occur in ``setup.py`` with output that contains something like the following::

      ModuleNotFoundError: No module named 'pybind11'
      Building wheel for gmxapi (pyproject.toml): finished with status 'error'
      ERROR: Failed building wheel for gmxapi
    Failed to build gmxapi
    ERROR: Could not build wheels for gmxapi, which is required to install pyproject.toml-based projects

The important information here is that ``pybind11`` was not found.

Build dependencies aren't always automatically installed.
Even if you are using ``pip``, you may have disabled automatic dependency fulfillment with an option like ``--no-build-isolation`` or ``--no-deps``.

In any case, the problem should be resolved by explicitly installing the ``pybind11``
Python package before attempting to build ``gmxapi``::

    pip install --upgrade pybind11

Couldn't find the ``gmxapi`` support library?
---------------------------------------------

If you don't want to "source" your :ref:`GMXRC <getting access to |Gromacs|>` file, you
can tell the package where to find a gmxapi compatible |Gromacs| installation with
``gmxapi_ROOT``. E.g. ``gmxapi_ROOT=/path/to/gromacs pip install .``

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
      "gmxapi_ROOT" to a directory containing one of the above files.  If "gmxapi"
      provides a separate development package or SDK, be sure it has been
      installed.

This could be because

* |Gromacs| is not already installed
* |Gromacs| was built without the CMake variable ``GMXAPI=ON``
* or if ``gmxapi_ROOT`` (or ``GROMACS_DIR``) is not a path containing directories
  like ``bin`` and ``share``.

If you are not a system administrator you are encouraged to install in a Python
virtual environment, created with virtualenv or Conda_.
Otherwise, you will need to specify the ``--user`` flag to ``pip``.

Two of the easiest problems to run into are incompatible compilers and
incompatible Python. Try to make sure that you use the same C and C++
compilers for |Gromacs|, for the Python package, and for the sample
plugin. These compilers should also correspond to the :command:`mpicc` compiler
wrapper used to
`compile mpi4py <https://mpi4py.readthedocs.io/en/stable/install.html>`__.
In order to build the Python
package, you will need the Python headers or development installation,
which might not already be installed on the machine you are using. (If
not, then you will get an error about missing :file:`Python.h` at some
point.) If you have multiple Python installations (or modules available
on an HPC system), you could try one of the other Python installations,
or you or a system administrator could install an appropriate Python dev
package. Alternatively, you might try installing your own Anaconda or
MiniConda in your home directory.

If an attempted installation fails with CMake errors about missing
“gmxapi”, make sure that |Gromacs| is installed and can be found during
installation. For instance,

::

    gmxapi_ROOT=/Users/eric/gromacs pip install --verbose gmxapi

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
the upstream branch may be reset to ``main`` after a new release is
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

.. _cmake: https://pypi.org/project/cmake/

.. _Conda: https://docs.conda.io/en/latest/

.. _mpi4py: https://pypi.org/project/mpi4py/

.. _networkx: https://pypi.org/project/networkx/

.. _numpy: https://www.numpy.org/

.. _pip: https://pip.pypa.io/en/stable/

.. _scikit-build: https://pypi.org/project/scikit-build/

.. _setuptools: https://pypi.org/project/setuptools/

.. _wheel: https://pypi.org/project/wheel/
