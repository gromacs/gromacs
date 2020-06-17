Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

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

