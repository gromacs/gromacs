New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Virtual site with single constructing atom
""""""""""""""""""""""""""""""""""""""""""

Added a virtual site that is constructed on top if its single constructing
atom. This can be useful for free-energy calculations.

Lower energy drift due to SETTLE
""""""""""""""""""""""""""""""""

|Gromacs| already applied an improvement to the center of mass calculation in
SETTLE to reduce energy drift in single precision. Now the center of mass
calculation is completely avoided, which significantly reduces the energy
drift when large coordinate values are present. This allows for accurate
simulations of systems with SETTLE up to 1000 nm in size (but note that
constraining with LINCS and SHAKE still introduces significant drift,
which limits the system size to 100 to 200 nm).
