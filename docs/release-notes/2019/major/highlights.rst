Highlights
^^^^^^^^^^

|Gromacs| 2019 was released on December 31st, 2018. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, with or
without GPUs, all enabled and automated by default. We are extremely
interested in your feedback on how well this worked on your
simulations and hardware. They are:

* Simulations now automatically run using update groups of atoms whose
  coordinate updates have only intra-group dependencies. These can
  include both constraints and virtual sites. This improves performance
  by eliminating overheads during the update, at no cost.
* Intel integrated GPUs are now supported with OpenCL for offloading
  non-bonded interactions.
* PME long-ranged interactions can now also run on a single AMD GPU
  using OpenCL, which means many fewer CPU cores are needed for good
  performance with such hardware.
