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
