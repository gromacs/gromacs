Dihedral principal component analysis
-------------------------------------

| :ref:`gmx angle <gmx angle>`, :ref:`gmx covar <gmx covar>`, 
  :ref:`gmx anaeig <gmx anaeig>`
| Principal component analysis can be performed in dihedral
  space \ :ref:`172 <refMu2005a>` using |Gromacs|. You start by defining the
  dihedral angles of interest in an index file, either using
  :ref:`gmx mk_angndx <gmx mk_angndx>` or otherwise. Then you use the
  :ref:`gmx angle <gmx angle>` program with the ``-or`` flag to
  produce a new :ref:`trr` file containing the cosine and sine
  of each dihedral angle in two coordinates, respectively. That is, in
  the :ref:`trr` file you will have a series of numbers
  corresponding to: cos(\ :math:`\phi_1`), sin(\ :math:`\phi_1`),
  cos(\ :math:`\phi_2`), sin(\ :math:`\phi_2`), ...,
  cos(\ :math:`\phi_n`), sin(\ :math:`\phi_n`), and the array is padded
  with zeros, if necessary. Then you can use this :ref:`trr`
  file as input for the :ref:`gmx covar <gmx covar>` program and perform
  principal component analysis as usual. For this to work you will need
  to generate a reference file (:ref:`tpr`,
  :ref:`gro`, :ref:`pdb` etc.) containing the same
  number of “atoms” as the new :ref:`trr` file, that is for
  :math:`n` dihedrals you need 2\ :math:`n`/3 atoms (rounded up if not
  an integer number). You should use the ``-nofit`` option
  for :ref:`gmx covar <gmx covar>` since the coordinates in the dummy
  reference file do not correspond in any way to the information in the
  :ref:`trr` file. Analysis of the results is done using
  :ref:`gmx anaeig <gmx anaeig>`.

Hydrogen bonds
--------------

| :ref:`gmx hbond <gmx hbond>`
| The program :ref:`gmx hbond <gmx hbond>`
  analyzes the *hydrogen bonds* (H-bonds) between all possible donors D
  and acceptors A. To determine if an H-bond exists, a geometrical
  criterion is used, see also :numref:`Fig. %s <fig-hbond>`:

  .. math::

     \begin{array}{rclcl}
     r       & \leq  & r_{HB}        & = & 0.35~\mbox{nm}    \\
     \alpha  & \leq  & \alpha_{HB}   & = & 30^o              \\
     \end{array}

.. _fig-hbond:

.. figure:: plots/hbond.*
   :width: 2.50000cm

   Geometrical Hydrogen bond criterion.

The value of :math:`r_{HB} = 0.35 \mathrm{nm}` corresponds to the first minimum
of the RDF of SPC water (see also :numref:`Fig. %s <fig-hbondinsert>`).

The program :ref:`gmx hbond <gmx hbond>` analyzes all hydrogen bonds
existing between two groups of atoms (which must be either identical or
non-overlapping) or in specified donor-hydrogen-acceptor triplets, in
the following ways:

.. _fig-hbondinsert:

.. figure:: plots/hbond-insert.*
    :width: 3.50000cm

    Insertion of water into an H-bond. (1) Normal H-bond between two
    residues. (2) H-bonding bridge via a water molecule.

-  Donor-Acceptor distance (:math:`r`) distribution of all H-bonds

-  Hydrogen-Donor-Acceptor angle (:math:`\alpha`) distribution of all
   H-bonds

-  The total number of H-bonds in each time frame

-  The number of H-bonds in time between residues, divided into groups
   :math:`n`-:math:`n`\ +\ :math:`i` where :math:`n` and
   :math:`n`\ +\ :math:`i` stand for residue numbers and :math:`i` goes
   from 0 to 6. The group for :math:`i=6` also includes all H-bonds for
   :math:`i>6`. These groups include the
   :math:`n`-:math:`n`\ +\ :math:`3`, :math:`n`-:math:`n`\ +\ :math:`4`
   and :math:`n`-:math:`n`\ +\ :math:`5` H-bonds, which provide a
   measure for the formation of :math:`\alpha`-helices or
   :math:`\beta`-turns or strands.

-  The lifetime of the H-bonds is calculated from the average over all
   autocorrelation functions of the existence functions (either 0 or 1)
   of all H-bonds:

   .. math:: C(\tau) ~=~ \langle s_i(t)~s_i (t + \tau) \rangle
             :label: eqnhbcorr

-  with :math:`s_i(t) = \{0,1\}` for H-bond :math:`i` at time
   :math:`t`. The integral of :math:`C(\tau)` gives a rough estimate of
   the average H-bond lifetime :math:`\tau_{HB}`:

   .. math::  \tau_{HB} ~=~ \int_{0}^{\infty} C(\tau) d\tau
              :label: eqnhblife

-  Both the integral and the complete autocorrelation function
   :math:`C(\tau)` will be output, so that more sophisticated analysis
   (*e.g.* using multi-exponential fits) can be used to get better
   estimates for :math:`\tau_{HB}`. A more complete analysis is given in
   ref. \ :ref:`173 <refSpoel2006b>`; one of the more fancy option is the Luzar
   and Chandler analysis of hydrogen bond kinetics \ :ref:`174 <refLuzar96b>`, :ref:`175 <refLuzar2000a>`.

-  An H-bond existence map can be generated of dimensions
   *# H-bonds*\ :math:`\times`\ *# frames*. The ordering is identical to
   the index file (see below), but reversed, meaning that the last
   triplet in the index file corresponds to the first row of the
   existence map.

-  Index groups are output containing the analyzed groups, all
   donor-hydrogen atom pairs and acceptor atoms in these groups,
   donor-hydrogen-acceptor triplets involved in hydrogen bonds between
   the analyzed groups and all solvent atoms involved in insertion.

