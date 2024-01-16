Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Improved performance for inhomogeneous systems
""""""""""""""""""""""""""""""""""""""""""""""

The performance of systems with a lot of empty space is improved
by optimizing the pair search grid for the effective atom density.

:issue:`4805`

Flexible hydrogen mass repartitioning using grompp
""""""""""""""""""""""""""""""""""""""""""""""""""

Instead of using pdb2gmx, which modifies the topology, a flexible
scheme for hydrogen mass repartitioning is now available in ``grompp``
through the ``mass-repartition-factor`` ``mdp`` option. This provides
easy access to a performance improvement of close to a factor two.

:issue:`4866`

Small performance regression to achieve more accurate pressure
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To lower the effect on Lennard-Jones pair interaction on the pressure,
the Verlet buffer has been increased for most simulations using default mdp settings.
This can lead to a few percent performance loss, in particular when using GPUs.
The effect will be strongest for systems with no or weak electrostatics,
which includes most coarse-grained systems.

:issue:`4861`
       
Reduced grompp and mdrun setup time for systems with many atom types
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Verlet buffer calculation could take many minutes for systems with thousands
of different Verlet buffer atom types (different atom type and charge).
Such times have now been reduced to seconds.

:issue:`4892`

With wall potentials, bonded interactions can now be run on GPUs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

HeFFTe multi-GPU FFT plan options are now configurable
""""""""""""""""""""""""""""""""""""""""""""""""""""""

New environment variables ``GMX_HEFFTE_RESHAPE_ALGORITHM``,
``GMX_HEFFTE_USE_GPU_AWARE``, ``GMX_HEFFTE_USE_PENCILS``, and
``GMX_HEFFTE_USE_REORDER`` permit the HeFFTe plan options to be
configured at run time. The performance obtained can vary with the
quality of implementation of e.g. the GPU-aware MPI library, as well
as the layout and number of the GPUs participating in the 3D-FFT.
Users can now find and use the best settings for their case. See
the HeFFTe documentation for more details.
