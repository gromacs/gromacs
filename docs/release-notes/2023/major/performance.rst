Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Increased default T- and P-coupling intervals
"""""""""""""""""""""""""""""""""""""""""""""

The default maximum values temperature and pressure coupling intervals
have been increased from 10 to 100 steps. These values are used when
the default value of -1 is specified in the mdp file and a lower value
is used when required for accurate integration. The improves the performance
of both GPU runs and parallel runs.


PME decomposition support with CUDA backend
""""""""""""""""""""""""""""""""""""""""""""

PME decomposition support has been added to the CUDA backend. With PME offloaded to the GPU, the number of PME ranks can
now be configured with ``-npme`` option (previously limited to 1). The implementation requires building |Gromacs|
with CUDA-aware MPI and with :ref:`NVIDIA's cuFFTMp library <cufftmp installation>`. GPU-based PME decomposition support still lacks substantial testing,
hence is included in the current release as an experimental feature and should be used with caution (with results compared to 
those from equivalent runs using a single PME GPU). This feature can be enabled using the ``GMX_GPU_PME_DECOMPOSITION`` environment 
variable. The |Gromacs| development team welcomes any feedback to help mature this feature.

:issue:`3884`
