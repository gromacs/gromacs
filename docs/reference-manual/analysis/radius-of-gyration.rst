.. _rg:

Radius of gyration and distances
--------------------------------

| :ref:`gmx gyrate <gmx gyrate>`, :ref:`gmx distance <gmx distance>`, 
  :ref:`gmx mindist <gmx mindist>`, :ref:`gmx mdmat <gmx mdmat>`,
  :ref:`gmx pairdist <gmx pairdist>`, :ref:`gmx xpm2ps <gmx xpm2ps>`
| To have a rough measure for the compactness of a structure, you can
  calculate the *radius of gyration* with the program
  :ref:`gmx gyrate <gmx gyrate>` as follows:

  .. math:: R_g ~=~ \left({\frac{\sum_i \|{\bf r}_i\|^2 m_i}{\sum_i m_i}}\right)^{{\frac{1}{2}}}
            :label: eqnrg

| where :math:`m_i` is the mass of atom :math:`i` and :math:`{\bf r}_i`
  the position of atom :math:`i` with respect to the center of mass of
  the molecule. It is especially useful to characterize polymer
  solutions and proteins. The program will also provide the radius of
  gyration around the coordinate axis (or, optionally, principal axes)
  by only summing the radii components orthogonal to each axis, for
  instance

  .. math:: R_{g,x} ~=~ \left({\frac{\sum_i \left( r_{i,y}^2 + r_{i,z}^2 \right) m_i}{\sum_i m_i}}\right)^{{\frac{1}{2}}}
            :label: eqnrgaxis

Sometimes it is interesting to plot the *distance* between two atoms, or
the *minimum* distance between two groups of atoms (*e.g.*: protein
side-chains in a salt bridge). To calculate these distances between
certain groups there are several possibilities:

*   The *distance between the geometrical centers* of two groups can be
    calculated with the program :ref:`gmx distance <gmx distance>`, as explained in
    sec. :ref:`bad`.

*   The *minimum distance* between two groups of atoms during time can
    be calculated with the program :ref:`gmx mindist <gmx mindist>`. It also calculates the
    *number of contacts* between these groups within a certain radius
    :math:`r_{max}`.

*   :ref:`gmx pairdist <gmx pairdist>` is a selection-enabled version of :ref:`gmx mindist <gmx mindist>`.

*   To monitor the *minimum distances between amino acid residues*
    within a (protein) molecule, you can use the program :ref:`gmx mdmat <gmx mdmat>`. This
    minimum distance between two residues A\ :math:`_i` and
    A\ :math:`_j` is defined as the smallest distance between any pair
    of atoms (i :math:`\in` A\ :math:`_i`, j :math:`\in` A\ :math:`_j`).
    The output is a symmetrical matrix of smallest distances between all
    residues. To visualize this matrix, you can use a program such as
    ``xv``. If you want to view the axes and legend or if you want to print
    the matrix, you can convert it with :ref:`xpm2ps <gmx xpm2ps>` into a Postscript
    :numref:`Fig. %s <fig-distm>`. 

.. _fig-distm:

.. figure:: plots/distm.*
       :width: 6.50000cm

       A minimum distance matrix for a
       peptide \ :ref:`168 <refSpoel96b>`.

*   Plotting these matrices for different time-frames, one can analyze
    changes in the structure, and *e.g.* forming of salt bridges.


