.. _pbc:

Periodic boundary conditions
----------------------------

.. _fig-pbc:

.. figure:: plots/pbctric.*
   :width: 9.00000cm

   Periodic boundary conditions in two dimensions.

The classical way to minimize edge effects in a finite system is to
apply *periodic boundary conditions*. The atoms of the system to be
simulated are put into a space-filling box, which is surrounded by
translated copies of itself (:numref:`Fig. %s <fig-pbc>`). Thus
there are no boundaries of the system; the artifact caused by unwanted
boundaries in an isolated cluster is now replaced by the artifact of
periodic conditions. If the system is crystalline, such boundary
conditions are desired (although motions are naturally restricted to
periodic motions with wavelengths fitting into the box). If one wishes
to simulate non-periodic systems, such as liquids or solutions, the
periodicity by itself causes errors. The errors can be evaluated by
comparing various system sizes; they are expected to be less severe than
the errors resulting from an unnatural boundary with vacuum.

There are several possible shapes for space-filling unit cells. Some,
like the *rhombic dodecahedron* and the *truncated octahedron* :ref:`20 <refAdams79>` are closer to being a sphere than a cube is, and
are therefore better suited to the study of an approximately spherical
macromolecule in solution, since fewer solvent molecules are required to
fill the box given a minimum distance between macromolecular images. At
the same time, rhombic dodecahedra and truncated octahedra are special
cases of *triclinic* unit cells; the most general space-filling unit
cells that comprise all possible space-filling shapes \ :ref:`21 <refBekker95>`.
For this reason, |Gromacs| is based on the triclinic unit cell.

|Gromacs| uses periodic boundary conditions, combined with the 
*minimum image convention*: only one – the nearest – image of each particle is
considered for short-range non-bonded interaction terms. For long-range
electrostatic interactions this is not always accurate enough, and
|Gromacs| therefore also incorporates lattice sum methods such as Ewald
Sum, PME and PPPM.

|Gromacs| supports triclinic boxes of any shape. The simulation box (unit
cell) is defined by the 3 box vectors :math:`{\bf a}`,\ :math:`{\bf b}`
and :math:`{\bf c}`. The box vectors must satisfy the following
conditions:

.. math:: a_y = a_z = b_z = 0
          :label: eqnboxrot

.. math:: a_x>0,~~~~b_y>0,~~~~c_z>0
          :label: eqnboxshift

.. math:: |b_x| \leq \frac{1}{2} \, a_x,~~~~
          |c_x| \leq \frac{1}{2} \, a_x,~~~~
          |c_y| \leq \frac{1}{2} \, b_y
          :label: eqnboxshift2

Equations :eq:`%s <eqnboxrot>` can always be satisfied by
rotating the box. Inequalities (:eq:`%s <eqnboxshift>`) and
(:eq:`%s <eqnboxshift2>`) can always be satisfied by adding
and subtracting box vectors.

Even when simulating using a triclinic box, |Gromacs| always keeps the
particles in a brick-shaped volume for efficiency, as illustrated in
:numref:`Fig. %s <fig-pbc>` for a 2-dimensional system. Therefore,
from the output trajectory it might seem that the simulation was done in
a rectangular box. The program :ref:`trjconv <gmx trjconv>` can be used to
convert the trajectory to a different unit-cell representation.

It is also possible to simulate without periodic boundary conditions,
but it is usually more efficient to simulate an isolated cluster of
molecules in a large periodic box, since fast grid searching can only be
used in a periodic system.

.. _fig-boxshapes:

.. figure:: plots/rhododec.*
        :width: 5.00000cm

        A rhombic dodecahedron (arbitrary orientation).


.. figure:: plots/truncoct.*
        :width: 5.00000cm

        A truncated octahedron (arbitrary orientation).

Some useful box types
~~~~~~~~~~~~~~~~~~~~~

.. |mathd| replace:: :math:`d`
.. |mathd3| replace:: :math:`d^{3}`
.. |mathd23| replace:: :math:`\frac{1}{2}\sqrt{2}~d^{3}`
.. |mathd70| replace:: :math:`0.707~d^{3}`
.. |mathd43| replace:: :math:`\frac{4}{9}\sqrt{3}~d^{3}`
.. |mathd77| replace:: :math:`0.770~d^{3}`
.. |math12d| replace:: :math:`\frac{1}{2}~d`
.. |math13d| replace:: :math:`\frac{1}{3}~d`
.. |math13dn| replace:: :math:`-\frac{1}{3}~d`
.. |math12s2| replace:: :math:`\frac{1}{2}\sqrt{2}~d`
.. |math12s3| replace:: :math:`\frac{1}{2}\sqrt{3}~d`
.. |math16s3| replace:: :math:`\frac{1}{6}\sqrt{3}~d`
.. |math13s6| replace:: :math:`\frac{1}{3}\sqrt{6}~d`
.. |math23s2| replace:: :math:`\frac{2}{3}\sqrt{2}~d`
.. |math13s2| replace:: :math:`\frac{1}{3}\sqrt{2}~d`
.. |angbc| replace:: :math:`\angle` **bc** 
.. |angac| replace:: :math:`\angle` **ac** 
.. |angab| replace:: :math:`\angle` **ab** 
.. |90deg| replace:: :math:`90^\circ`
.. |60deg| replace:: :math:`60^\circ`
.. |70deg| replace:: :math:`70.53^\circ`
.. |109deg| replace:: :math:`109.47^\circ`

.. _table-boxtypes:

.. table:: Overview over different box types
    :align: center
    :widths: auto

    +--------------+-----------+-----------+-----------------------------------+------------------------------+
    | box type     | image     | box       | box vectors                       | box vector angles            | 
    |              |           |           +---------+------------+------------+---------+----------+---------+
    |              | distance  | volume    | **a**   | **b**      | **c**      | |angbc| | |angac|  | |angab| |
    +==============+===========+===========+=========+============+============+=========+==========+=========+
    |              |           |           | |mathd| |   0        |   0        |         |          |         |
    |              |           |           +---------+------------+------------+         |          |         |
    | cubic        | |mathd|   | |mathd3|  |   0     | |mathd|    |   0        | |90deg| | |90deg|  | |90deg| |
    |              |           |           +---------+------------+------------+         |          |         |
    |              |           |           |   0     |   0        | |mathd|    |         |          |         |
    +--------------+-----------+-----------+---------+------------+------------+---------+----------+---------+
    | rhombic      |           | |mathd23| | |mathd| | 0          | |math12d|  |         |          |         |
    |              |           |           +---------+------------+------------+         |          |         |
    | dodecahedron | |mathd|   | |mathd70| | 0       | |mathd|    | |math12d|  | |60deg| | |60deg|  | |90deg| |
    |              |           |           +---------+------------+------------+         |          |         |
    | (xy-square)  |           |           | 0       | 0          | |math12s2| |         |          |         |
    +--------------+-----------+-----------+---------+------------+------------+---------+----------+---------+
    | rhombic      |           | |mathd23| | |mathd| | |math12d|  | |math12d|  |         |          |         |
    |              |           |           +---------+------------+------------+         |          |         |
    | dodecahedron | |mathd|   | |mathd70| | 0       | |math12s3| | |math16s3| | |60deg| | |60deg|  | |60deg| |
    |              |           |           +---------+------------+------------+         |          |         |
    | (xy-         |           |           | 0       | 0          | |math13s6| |         |          |         |
    | hexagon)     |           |           |         |            |            |         |          |         |
    +--------------+-----------+-----------+---------+------------+------------+---------+----------+---------+
    | truncated    |           | |mathd43| | |mathd| | |math13d|  | |math13dn| |         |          |         |
    |              |           |           +---------+------------+------------+         |          |         |
    | octahedron   | |mathd|   | |mathd77| | 0       | |math23s2| | |math13s2| | |70deg| | |109deg| | |70deg| |
    |              |           |           +---------+------------+------------+         |          |         |
    |              |           |           | 0       | 0          | |math13s6| |         |          |         |
    +--------------+-----------+-----------+---------+------------+------------+---------+----------+---------+

The three most useful box types for simulations of solvated systems are
described in :numref:`Table %s <table-boxtypes>`. The rhombic
dodecahedron (:numref:`Fig. %s <fig-boxshapes>`) is the smallest and
most regular space-filling unit cell. Each of the 12 image cells is at
the same distance. The volume is 71% of the volume of a cube having the
same image distance. This saves about 29% of CPU-time when simulating a
spherical or flexible molecule in solvent. There are two different
orientations of a rhombic dodecahedron that satisfy equations
:eq:`%s <eqnboxrot>`, :eq:`%s <eqnboxshift>` and
:eq:`%s <eqnboxshift2>`. The program :ref:`editconf <gmx editconf>`
produces the orientation which has a square intersection with the
xy-plane. This orientation was chosen because the first two box vectors
coincide with the x and y-axis, which is easier to comprehend. The other
orientation can be useful for simulations of membrane proteins. In this
case the cross-section with the xy-plane is a hexagon, which has an area
which is 14% smaller than the area of a square with the same image
distance. The height of the box (:math:`c_z`) should be changed to
obtain an optimal spacing. This box shape not only saves CPU time, it
also results in a more uniform arrangement of the proteins.

Cut-off restrictions
~~~~~~~~~~~~~~~~~~~~

The minimum image convention implies that the cut-off radius used to
truncate non-bonded interactions may not exceed half the shortest box
vector:

.. math:: R_c < {\frac{1}{2}}\min(\|{\bf a}\|,\|{\bf b}\|,\|{\bf c}\|),
          :label: eqnphysicalrc

because otherwise more than one image would be within the cut-off
distance of the force. When a macromolecule, such as a protein, is
studied in solution, this restriction alone is not sufficient: in
principle, a single solvent molecule should not be able to ‘see’ both
sides of the macromolecule. This means that the length of each box
vector must exceed the length of the macromolecule in the direction of
that edge *plus* two times the cut-off radius :math:`R_c`. It is,
however, common to compromise in this respect, and make the solvent
layer somewhat smaller in order to reduce the computational cost. For
efficiency reasons the cut-off with triclinic boxes is more restricted.
For grid search the extra restriction is weak:

.. math:: R_c < \min(a_x,b_y,c_z)
         :label: eqngridrc
   

For simple search the extra restriction is stronger:

.. math:: R_c < {\frac{1}{2}}\min(a_x,b_y,c_z)
          :label: eqnsimplerc

Each unit cell (cubic, rectangular or triclinic) is surrounded by 26
translated images. A particular image can therefore always be identified
by an index pointing to one of 27 *translation vectors* and constructed
by applying a translation with the indexed vector (see :ref:`forces`).
Restriction :eq:`%s <eqngridrc>` ensures that only 26 images need to be
considered.
