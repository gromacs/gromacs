Highlights
^^^^^^^^^^

|Gromacs| 2019 was released on INSERT DATE HERE. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, including
allowing GPU offloading with more accelerators and of more calculation
types. Also, domain decomposition can now be performed using the concept
of update groups instead of single atoms.

We are extremely interested in your feedback on how well this worked on
your simulations and hardware. They are:

* Running simulations with OpenCL on Nvidia and Intel hardware.
* SIMD accelerations of bonded interactions.
* PME on GPU with OpenCL.
* Update groups to reduce global communication.

There are some new features available also:

* Updates to gmx cluster code.
* Initial support for GMX Python API.
* Documentation completely available as HTML (including reference manual!)
* |Gromacs| will now be build more reproducible.
* COM can be used as reference for PBC in pull simulations.
