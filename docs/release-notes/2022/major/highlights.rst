Highlights
^^^^^^^^^^

|Gromacs| 2022 was released on INSERT DATE HERE. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, with or
without GPUs, all enabled and automated by default. In addition,
several new features are available for running simulations. We are extremely
interested in your feedback on how well the new release works on your
simulations and hardware. The new features are:

* Free-energy kernels are accelerated using SIMD, which make free-energy
  calculations up to three times as fast when using GPUs
* A new formulation of the soft-cored non-bonded interactions for free-energy calculations allows for a finer control of the alchemical transformation pathways
* New transformation pull coordinate allows arbibrary mathematical transformations of one of more other pull coordinates
* Cool quotes music play list


.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!
