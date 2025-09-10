.. _plumed:

Using PLUMED
------------

PLUMED functionality is enabled by using ``gmx mdrun -plumed plumed.dat``.

The interface will look for an environment variable ``PLUMED_KERNEL`` that should 
contain the path and the name of a shared object that contains the PLUMED kernel,
and usually is called ``libPlumedKernel.so``.

If the library is not present an error message will inform the user to export the ``PLUMED_KERNEL`` variable.

Usually the PLUMED kernel is stored in ``$plumed_prefix/lib/libPlumedKernel.so``, 
so it should be enough to ``export PLUMED_KERNEL=$plumed_prefix/lib/libPlumedKernel.so``, 
where ``$plumed_prefix`` is the PLUMED installation prefix.

Configuration files for input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLUMED input files are plain text files describing the analysis process and how PLUMED should influence the |Gromacs| simulation.
The full documentation is available `here <https://www.plumed.org/doc>`_

Limitations
^^^^^^^^^^^

The current implementation of the PLUMED interface does not support the following features:

* interaction with |Gromacs| energy
* PLUMED is incompatible with thread-MPI (PLUMED exits with an error when more than 1 thread-MPI rank is used)
* The various replica exchange flavors are not yet implemented
