Sample MD restraint plugin
==========================

This [repository](https://github.com/kassonlab/sample_restraint)
provides a complete and working implementation of a few GROMACS
restraint potentials. It is intended as both a tutorial and as a
template for implementing new custom restraint potentials.

Restraint potentials are implemented as \"plugins\" to GROMACS.
GROMACS must be [configured with
GMXAPI=ON](http://manual.gromacs.org/current/install-guide/index.html#gmxapi-external-api)

The plugin potentials are loaded and configured via Python and are
compatible with the [gmxapi](https://github.com/kassonlab/gmxapi) Python
package for MD simulation workflows.

For a quick start, consider pulling a recent Docker image that has
already been configured for gmxapi and this plug-in.
**todo:** check and update (ref: [GitHub issue 230](https://github.com/kassonlab/gmxapi/issues/230))

Reference:

Irrgang, M. E., Hays, J. M., & Kasson, P. M. gmxapi: a high-level
interface for advanced control and extension of molecular dynamics
simulations. *Bioinformatics* 2018. DOI:
[10.1093/bioinformatics/bty484](https://doi.org/10.1093/bioinformatics/bty484)

Repository Contents
-------------------

This repository uses CMake to build and install a Python C++ extension
package.

-   `CMakeLists.txt`, `cmake/FindGROMACS.cmake`, and
    `src/CMakeLists.txt` provide necessary CMake infrastructure. You
    should not need to edit these.
-   `src/cpp` contains a header and `cpp` file for each restraint
    potential built with this module. When adding new potentials, you
    will update `CMakeLists.txt` to create build targets. Use the
    existing potentials as examples.
-   `src/pythonmodule/` contains `CMakeLists.txt`, `export_plugin.h`,
    and `export_plugin.cpp`. When you have written a new potential, you
    can add it to `CMakeLists.txt` and `export_plugin.cpp`. This is the
    code that produces the C++ extension for Python.
    `EnsemblePotential` applies a restrained ensemble potential and
    uses additional facilities provided by gmxapi.
-   <strike>`src/pybind11` is just a copy of the Python bindings framework from
    the Pybind project (ref <https://github.com/pybind/pybind11> ). It
    is used to wrap the C++ restraint code and give it a Python
    interface.</strike> Note: pybind is currently retrieved while configuring
    with CMake. Ref redmine issues [3027](https://redmine.gromacs.org/issues/3027)
    and [3033](https://redmine.gromacs.org/issues/3033)
-   `tests/` contains C++ and Python tests for the provided code. Update
    `CMakeLists.txt` to add your own, based on these examples. C++ unit
    tests use [googletest](https://github.com/google/googletest). Python
    tests use the [pytest](https://docs.pytest.org/en/latest/). Refer to
    those respective projects for more about how they make test-writing
    easier. Note: googletest is currently downloaded while configuring with
    CMake. Ref [3033](https://redmine.gromacs.org/issues/3033)
-   `examples` contains a sample SLURM job script and
    `restrained-ensemble.py` gmxapi script that have been used to do
    restrained ensemble simulations. `example.py` and `example.ipynb`
    explore a toy alanine dipeptide system. `strip_notebook.py` is a
    helper script to remove extra output and state data from an iPython
    notebook before checking updates back into the repository.
-   `Dockerfile` is a recipe to build a Docker image from the root of
    the repository. **todo:** Check and update.
    ref: GitHub issue [230](https://github.com/kassonlab/gmxapi/issues/230)

Docker quick-start
------------------

**todo: check and update** ref: [GitHub issue 230](https://github.com/kassonlab/gmxapi/issues/230)

Pull the docker image and launch a container with port 8888 on the host
mapped to port 8888 in the container. :

    docker run --rm -ti -p 8888:8888 gmxapi/sample_restraint:devel

Note that the `--rm` option tells docker not to save any changes you
make after launching the container. You can, however, download any
changes you make to the notebook through the web interface. Refer to the
[Docker documentation](https://docs.docker.com) for more options on
managing containers.

You should then see something like the following, but with a different
`token` for the URL. Open the URL in a browser on the same (host)
machine to access the notebook server. Browse to `sample_restraint` and
`examples` and then launch the `example` notebook for an interactive
walk-through. Example output:

    Execute the command: jupyter notebook
    [I 15:26:07.683 NotebookApp] Writing notebook server cookie secret to /home/jovyan/.local/share/jupyter/runtime/notebook_cookie_secret
    [W 15:26:08.184 NotebookApp] WARNING: The notebook server is listening on all IP addresses and not using encryption. This is not recommended.
    [I 15:26:08.223 NotebookApp] JupyterLab alpha preview extension loaded from /opt/conda/lib/python3.6/site-packages/jupyterlab
    [I 15:26:08.230 NotebookApp] Serving notebooks from local directory: /home/jovyan
    [I 15:26:08.230 NotebookApp] 0 active kernels
    [I 15:26:08.230 NotebookApp] The Jupyter Notebook is running at:
    [I 15:26:08.230 NotebookApp] http://[all ip addresses on your system]:8888/?token=948d611453ea3f03ad406dc375bfc186c4315fa68c50e23d
    [I 15:26:08.230 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
    [C 15:26:08.231 NotebookApp]

        Copy/paste this URL into your browser when you connect for the first time,
        to login with a token:
            http://localhost:8888/?token=948d611453ea3f03ad406dc375bfc186c4315fa68c50e23d

The basics
----------

This repository provides a potentially useful plugin, but also serves as
documentation by example and as a template for developing GROMACS
extension code in the gmxapi framework.

### Build and install

To download, build, and install, you may need to first install `wget`,
`git`, and/or `cmake`.

The plugin requires libgmxapi to build. See
[gmxapi](http://manual.gromacs.org/current/install-guide/index.html#gmxapi-external-api).
Download an official release from http://www.gromacs.org or the latest gmxapi
development branch from https://github.com/kassonlab/gmxapi/

We use CMake to configure and build a C++ library and a Python module
for interacting with it.

After installing GROMACS, either source the GMXRC file provided with the GROMACS
installation or set `gmxapi_DIR` to the GROMACS installation path.

The GROMACS installation provides some additional CMake infrastructure to help us build compatible client software.
To help set the correct compilers, specify the CMake toolchains file with,
*e.g.*, `-DCMAKE_TOOLCHAIN_FILE=/usr/local/gromacs/share/cmake/gromacs/gromacs-toolchain.cmake` (for GROMACS installed
 to `/usr/local/gromacs`).
**todo:** Link to GROMACS docs for the toolchains file.

We recommend installing and using this code in a Python virtual
environment. (See the documentation for your `gmxapi` distribution or
<http://gmxapi.readthedocs.io/en/latest/install.html> ) Accordingly, if
you choose to install the plugin rather than just to use it out of
its build directory, consider whether you want to have to set your
`PYTHONPATH` environment variable or where you can install it that
Python will find it. You can explicitly set the installation location by
setting `-DGMXPLUGIN_INSTALL_PATH=/path/to/install/directory` or you can
let CMake determine an appropriate location automatically for your
Python interpreter. If you have administrative privileges (such as when
running on a desktop computer) or if you are using a Python virtual
environment (recommended), you don\'t need to specify anything
additional. If you are an unprivileged user (such as on a shared
machine) and are not in a Python virtual environment, set
-DGMXPLUGIN\_USER\_INSTALL=ON to install into the \"user\" Python
packages directory in your home directory. (Equivalent to the `--user`
option to `pip`)

If you have multiple Python installations or just want to be
unambiguous, provide CMake with the Python interpreter you wish to use
(the same as you are using for `gmxapi`) with
`-DPYTHON_EXECUTABLE=/path/to/python3`.

From the root directory of the GROMACS source, the sample_restraint source code is in
`python_packaging/sample_restraint`

    cd python_packaging/sample_restraint
    mkdir build
    cd build
    # Get the GROMACS environment settings.
    source $HOME/gromacs/bin/GMXRC
    # Configure the build environment with CMake
    cmake ..
    # or
    # cmake .. -DGMXPLUGIN_INSTALL_PATH=/path/to/install/directory
    # or
    # cmake .. -DGMXPLUGIN_USER_INSTALL=ON -DPYTHON_EXECUTABLE=`which python3`
    # Build myplugin.
    make
    # run C++ tests
    make test
    # optionally, install
    make install

If you choose not to install the plugin module, you can tell Python
where to find it by setting your PYTHONPATH environment variable. For
instance, while still in the build directory:

    export PYTHONPATH=`pwd`/src/pythonmodule

The Python `gmxapi` package is required for testing.
See the [README.md](../README.md) 
file in the parent directory.

### Running

The `examples` directory contains some sample scripts for running
`gmxapi` workflows using the restraint potential samples in this
repository. You may also find [tests/test_binding.py](tests/test_binding.py) informative.

For a basic walk-through with a toy system, launch a Jupyter notebook
server and navigate to `examples/example.py`

**todo** These scripts have not been checked since migrating to the GROMACS source repository.

### What\'s going on

This sample project builds several C++ object files, which are used to build a
Python module named `myplugin`.

When setting up a workflow, a Python script provides gmxapi with
parameters and a factory function for a plugin restraint potential. This
Python interface is defined in `src/pythonmodule/export_plugin.cpp`.
When a Session is launched, a C++ object that performs restraint force
calculations is created and given to the GROMACS library. During each MD
step, part of the MD force evaluation includes a call to the
calculations performed by the restraint. For the pair restraints
demonstrated here, GROMACS provides relative coordinates of two atomic
sites to the calculation code in the plugin. If multiple restrained
pairs are needed, multiple restraints are attached to the simulation.
Coordination across an ensemble of simulations is possible using
resources provided by the Session.

Fundamentally, a new restraint potential is implemented by creating a
class that provides a `calculate()` method and using wrappers to give it
interfaces to GROMACS and to Python. C++ wrappers allow the basic class
implementing the potential to be presented to the GROMACS library in a
way that can be used to evaluate forces during a simulation. Other C++
template code wraps the potential in a portable way so that it can be
passed to GROMACS through a Python interface and to receive parameters
from the Python interpreter. Pybind11 syntax in `export_plugin.cpp`
provides the code to actually expose the plugin as a class in a Python
module that is compatible with the `gmx` package provided in the
`gmxapi` project.

By version 0.1.0, additional wrappers and boilerplate code will be
migrated out of the files that define the `calculate()` methods. Until
then, some amount of copy-and-paste or editing is necessary to implement
a new potential. Refer to `src/cpp/harmonicpotential.h` and to
`src/cpp/harmonicpotential.cpp` for a documented example of a simple
pair restraint. A more complex example is found in the
`ensemblepotential` files. The code in `src/cpp` is sufficient to
produce testable object code, but the Python module is exported in
`src/pythonmodule/export_plugin.cpp`. If you add additional source files
for a new potential, you will need to update `src/cpp/CMakeLists.txt` as
well.

Python tests
------------

For the Python-level testing, you will need `pytest` and `gmxapi`. We
recommend setting up a Python virtual environment as described in the gmxapi installation instructions.

You will also need a functioning MPI installation and the `mpi4py`
package.

Python tests can be run from the root directory of the repository after
building. Assuming you built in a subdirecory of the repository named
`build` (as above):

    PYTHONPATH=build/src/pythonmodule/ python -m pytest tests

This command causes the directory named `tests` to be explored for
Python files with names like `test_*.py` or `*_test.py`. Matching files
will be imported and any functions with similarly obvious names will be
run and errors reported. In particular, `assert` statements will be
evaluated to perform individual tests. See also
<https://docs.pytest.org/en/latest/goodpractices.html#test-discovery>

The tests assume that the package is already installed or is available
on the default Python path (such as by setting the `PYTHONPATH`
environment variable). If you just run `pytest` with no arguments, it
will discover and try to run tests from elsewhere in the repository that
were not intended, and they will fail.

To run the full set of tests for the ensemble workflow features, first
make sure that you have an MPI-capable environment and `mpi4py`
installed. Refer to <http://mpi4py.readthedocs.io/en/stable/> and
<https://github.com/kassonlab/gmxapi> for more information.

The ensemble tests assume that 2 ranks are available. After installing
the plugin, run (for example):

    mpiexec -n 2 python -m mpi4py -m pytest

**todo** check and update the following. (ref: [GitHub issue 230](https://github.com/kassonlab/gmxapi/issues/230))

If you do not have MPI set up for your system, you could build a docker
image using the Dockerfile in this repository.

    docker build -t samplerestraint . Dockerfile
    docker run --cpus 2 --rm -ti samplerestraint bash -c \
        "cd /home/jovyan/sample_restraint/tests && 
        mpiexec -n 2 python -m mpi4py -m pytest"

To test with a pre-built image from our docker hub repository, do

    docker run --cpus 2 --rm -ti gmxapi/sample_restraint bash -c \
            "cd /home/jovyan/sample_restraint/tests && 
            mpiexec -n 2 python -m mpi4py -m pytest"
