Highlights
^^^^^^^^^^

|Gromacs| 2024.0 was released on January 30th, 2024. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, with or
without GPUs, all enabled and automated by default. In addition,
several new features are available for running simulations. We are extremely
interested in your feedback on how well the new release works on your
simulations and hardware. The new features are:

* The `Colvars <https://colvars.github.io>`_ library can now be used natively
  from |Gromacs|. This simplifies the use of advanced enhanced sampling simulations.

* Reduced artifacts from Lennard-Jones pair interactions on the pressure by
  a configurable increase of the Verlet buffer. Can lead to a slight performance
  loss, especially for coarse-grained systems.

* Corrected several aspects of the deform option. Now simulations with box deformation
  behave correctly under high shear or when a solid or membrane fractures. This also
  means that the deform option is now suitable for computing viscosities.

* New option for hydrogen mass repartitioning in ``grompp`` enables easy access
  to performance improvements.

* Improvements to AWH, such as better control of the histogram growth factor as
  well as enabling automatic scaling of the target distribution based on the AWH
  friction metric.

* Configurable HeFFTe multi-GPU FFT options lets users fine-tune the settings for
  specific use-cases.

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!
