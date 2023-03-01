Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Default values for temperature and pressure coupling intervals are now 10
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
With the default mdp input value of -1 for nsttcouple and nstpcouple, grompp would
set these values to nstlist. Now these are set to 10 and thus independent of nstlist
(note that grompp may choose smaller values when needed for accurate integration).

Uniform and manual CMake GPU-support configuration
""""""""""""""""""""""""""""""""""""""""""""""""""
The GPU accelerations setup has been changed to be uniform for CUDA and OpenCL. Either
option is now enabled by setting GMX_GPU to CUDA or OpenCL in the CMake configuration.
To simplify the CMake code, we have also moved away from automated option selection
based on the build host. In particular, this means that CUDA will not be enabled unless
the GMX_GPU option is explicitly enabled, and CMake will no longer perform the extra
steps of trying to detect hardware and propose to install CUDA if hardware is available.
Apart from the simplification, this should also make it easier to handle multiple
different accelerator APIs targeting e.g. NVIDIA hardware.

Configuration-time trivalue options changed from autodetection to boolean on/off
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
To simplify the CMake configuration and avoid having multiple settings that
change outside of the users direct control we have removed the support for
automatically setting booleans. GMX_BUILD_HELP and GMX_HWLOC are now
disabled by default, while GMX_LOAD_PLUGINS is enabled by default.

gmxapi C++ interface
""""""""""""""""""""

``gmxapi::Context`` is now created with ``gmxapi::createContext()``, which allows
the client to provide an MPI communicator for the library to use instead of its default
(e.g MPI_COMM_WORLD). MPI-enabled clients may use the :file:`gmxapi/mpi/gmxapi_mpi.h`
template header and the ``assignResource()`` helper to generate the argument to
``createContext``.

Unification of several CUDA and OpenCL environment variables
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The environment variables that had exactly the same meaning in OpenCL and CUDA were unified:

* GMX_CUDA_NB_ANA_EWALD and GMX_OCL_NB_ANA_EWALD into GMX_GPU_NB_ANA_EWALD
* GMX_CUDA_NB_TAB_EWALD and GMX_OCL_NB_TAB_EWALD into GMX_GPU_NB_TAB_EWALD
* GMX_CUDA_NB_EWALD_TWINCUT and GMX_OCL_NB_EWALD_TWINCUT into GMX_GPU_NB_EWALD_TWINCUT

Dysfunctional parts of the QMMM interface has been removed
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Currently, |Gromacs| supports QM/MM officially only via MiMiC; a new CP2K QM/MM interface is being
developed within BioExcel. All other QM/MM
support has been untested and likely dysfunctional for years and has now been removed from .mdp
input and output, resulting in smaller .mdp output files from grompp.
