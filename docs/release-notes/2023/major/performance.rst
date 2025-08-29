Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Update will run on GPU by default
"""""""""""""""""""""""""""""""""

The mdrun ``-update auto`` will by default map to GPU if supported.
This gives a significant performance improvement with a single MPI rank.

Increased default T- and P-coupling intervals
"""""""""""""""""""""""""""""""""""""""""""""

The default maximum values temperature and pressure coupling intervals
have been increased from 10 to 100 steps. These values are used when
the default value of -1 is specified in the mdp file and a lower value
is used when required for accurate integration. The improves the performance
of both GPU runs and parallel runs.

The global communication frequency is independent of nstlist
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The global communication frequency no longer depends on nstlist.
This can improve performance in simulations using GPUs in particular.

PME decomposition support with CUDA and SYCL backends
"""""""""""""""""""""""""""""""""""""""""""""""""""""

PME decomposition support has been added to CUDA and SYCL
backends. With PME offloaded to the GPU, the number of PME ranks can
now be configured with ``-npme`` option (previously limited to 1). The
implementation requires building |Gromacs| with GPU-aware MPI and
either with the :ref:`cuFFTMp library <cufftmp installation>` in a
CUDA build configuration, or with :ref:`heFFTe <heffte installation>`
in either a CUDA or SYCL build configuration.

GPU-based PME decomposition support still lacks substantial testing,
hence is included in the current release as an experimental feature
and should be used with caution (with results compared to those from
equivalent runs using a single PME GPU). This feature can be enabled
using the ``GMX_GPU_PME_DECOMPOSITION`` environment variable. The
|Gromacs| development team welcomes any feedback to help mature this
feature.

:issue:`3884`
:issue:`4090`

CUDA Graphs for GPU-resident Steps
""""""""""""""""""""""""""""""""""

New CUDA functionality has been introduced, allowing GPU activities
to be launched as a single CUDA graph on each step rather than multiple
activities scheduled to multiple CUDA streams. It only works for those
cases which already support GPU-resident steps (where all force and
update calculations are GPU-accelerated). This offers performance
advantages, especially for small cases, through reduction in both CPU
and GPU side scheduling overheads. The feature can optionally be
activated via the ``GMX_CUDA_GRAPH`` environment variable. 

:issue:`4277`

VkFFT support
"""""""""""""

For AMD GPUs, VkFFT has been integrated to provide performance improvements.
Using this library is supported in all non-decomposed PME simulations (either
single rank or single separate PME rank) and can be enabled with ``-DGMX_GPU_FFT_LIBRARY=VKFFT``
when using AdaptiveCpp_ (hipSYCL).

:issue:`4052`
