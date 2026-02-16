Highlights
^^^^^^^^^^

|Gromacs| 2026 was released on January 19th, 2026. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, with or
without GPUs, all enabled and automated by default. In addition,
several new features are available for running simulations. We are extremely
interested in your feedback on how well the new release works on your
simulations and hardware. The new features are:

* Two `protein force fields from AMBER <https://ambermd.org/AmberModels_proteins.php>`_,
  ff14SB and ff19SB, have been ported to |Gromacs| as AMBER14SB and AMBER19SB.
  The ports also include the new OPC and OPC3 water models as well as several others.

* Expanded support for running simulations with Neural Network Potential models,
  now including link atom treatment for NNP/MM, pairlist input, and electrostatic
  embedding models.

* Experimental support for `H5MD`_ as a trajectory output format for ``mdrun``.

* Full support for using HIP as the GPU backend for AMD devices.

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!
