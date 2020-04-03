Highlights
^^^^^^^^^^

|Gromacs| 2020 was released on January 1, 2020. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, with or
without GPUs, all enabled and automated by default. In addition,
several new features are available for running simulations. We are extremely
interested in your feedback on how well the new release works on your
simulations and hardware. The new features are:

* Density-guided simulations allow "fitting" atoms into three-dimensional
  density maps. 
* Inclusion of gmxapi 0.1, an API and user interface for managing
  complex simulations, data flow, and pluggable molecular dynamics extension code.
* New modular simulator that can be built from individual objects describing different
  calculations happening at each simulation step.
* Parrinello-Rahman pressure coupling is now also available for the md-vv integrator.
* Running almost the entire simulation step on a single CUDA compatible GPU for supported
  types of simulations, including coordinate update and constraint calculation.


.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!
