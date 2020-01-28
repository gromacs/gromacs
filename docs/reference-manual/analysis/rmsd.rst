.. _rmsd:

Root mean square deviations in structure
----------------------------------------

| :ref:`gmx rms <gmx rms>`, :ref:`gmx rmsdist <gmx rmsdist>`
| The *root mean square deviation* (:math:`RMSD`) of certain atoms in a
  molecule with respect to a reference structure can be calculated with
  the program :ref:`gmx rms <gmx rms>` by least-square fitting the structure to the
  reference structure (:math:`t_2 = 0`) and subsequently calculating the
  :math:`RMSD` (:eq:`eqn. %s <eqnrmsd>`).

  .. math:: RMSD(t_1,t_2) ~=~ \left[\frac{1}{M} \sum_{i=1}^N m_i \|{\bf r}_i(t_1)-{\bf r}_i(t_2)\|^2 \right]^{\frac{1}{2}}
            :label: eqnrmsd

| where :math:`M = \sum_{i=1}^N m_i` and :math:`{\bf r}_i(t)` is the
  position of atom :math:`i` at time :math:`t`. **Note** that fitting
  does not have to use the same atoms as the calculation of the
  :math:`RMSD`; *e.g.* a protein is usually fitted on the backbone atoms
  (N, C\ :math:`_{\alpha}`, C), but the :math:`RMSD` can be computed of the
  backbone or of the whole protein.

Instead of comparing the structures to the initial structure at time
:math:`t=0` (so for example a crystal structure), one can also calculate
:eq:`eqn. %s <eqnrmsd>` with a structure at time :math:`t_2=t_1-\tau`. This
gives some insight in the mobility as a function of :math:`\tau`. A
matrix can also be made with the :math:`RMSD` as a function of
:math:`t_1` and :math:`t_2`, which gives a nice graphical interpretation
of a trajectory. If there are transitions in a trajectory, they will
clearly show up in such a matrix.

Alternatively the :math:`RMSD` can be computed using a fit-free method
with the program :ref:`gmx rmsdist <gmx rmsdist>`:

.. math:: RMSD(t) ~=~     \left[\frac{1}{N^2}\sum_{i=1}^N \sum_{j=1}^N    \|{\bf r}_{ij}(t)-{\bf r}_{ij}(0)\|^2\right]^{\frac{1}{2}}
          :label: eqnrmsdff

where the *distance* **r**\ :math:`_{ij}` between atoms at time
:math:`t` is compared with the distance between the same atoms at time
:math:`0`.
