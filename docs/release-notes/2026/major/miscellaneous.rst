Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!


Status of MDModules that require an external library is now logged
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Support status of CP2K QM/MM, Colvars, Plumed & Neural Network Potential modules
along with the version of libraries they use are now reported
in ``md.log`` and ``gmx --version``.

Modular simulator reports why it cannot be used
"""""""""""""""""""""""""""""""""""""""""""""""

When a Velocity-verlet integrator simulation cannot use the modular
simulator, the reasons are logged. If the modular simulator was
required by the user, a fatal error is given along with the reasons
why it cannot be run.

:issue:`5339`

Invalid ``-nstlist`` value is now an error instead of a warning
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

An invalid (e.g., too large) ``-nstlist`` value now triggers a fatal error instead
of falling back to the fixed value from the TPR.

:issue:`5365`

FFTW3 now builds with ```CMAKE_C_COMPILER```
""""""""""""""""""""""""""""""""""""""""""""

Previously, the cmake flag ```-DGMX_BUILD_OWN_FFTW=ON -DGMX_FFT_LIBRARY=fftw3```
leads to FFTW3 being compiled with the system default compiler. The ```CMAKE_C_COMPILER```
is now passed to the FFTW3 configure command.

:issue:`5344`
