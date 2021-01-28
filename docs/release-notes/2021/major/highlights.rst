Highlights
^^^^^^^^^^

|Gromacs| 2021 was released on January 28th, 2021. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, with or
without GPUs, all enabled and automated by default. In addition,
several new features are available for running simulations. We are extremely
interested in your feedback on how well the new release works on your
simulations and hardware. The new features are:

* Support for multiple time stepping, allowing for simple near doubling of simulation speed and is intended to replace the virtual site treatment
* Ability to use stochastic cell rescaling barostat for equilibration and production simulations
* Preliminary support for using SYCL as accelerator framework
* Support for performing free energy perturbation with AWH
* Support PME offloading to GPU for free energy simulations
* Support for ARM SVE and Fujitsu A64FX (contribution by Research Organization for Information Science and Technology (RIST))
* New nonbonded interaction API with NB-LIB (in collaboration with PRACE)
* New |Gromacs| logo!

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!
