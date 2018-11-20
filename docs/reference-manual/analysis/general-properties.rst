General properties
------------------

| :ref:`gmx energy <gmx energy>`, :ref:`gmx traj <gmx traj>`
| To analyze some or all *energies* and other properties, such as *total
  pressure*, *pressure tensor*, *density*, *box-volume* and *box-sizes*,
  use the program :ref:`gmx energy <gmx energy>`. A choice can be made from a
  list a set of energies, like potential, kinetic or total energy, or
  individual contributions, like Lennard-Jones or dihedral energies.

The *center-of-mass velocity*, defined as

.. math:: {\bf v}_{com} = {1 \over M} \sum_{i=1}^N m_i {\bf v}_i
          :label: eqncomvelocity

with :math:`M = \sum_{i=1}^N m_i` the total mass of the system, can be
monitored in time by the program :ref:`gmx traj <gmx traj>` ``-com -ov``. It is however
recommended to remove the center-of-mass velocity every step (see
chapter :ref:`algorithms`)!
