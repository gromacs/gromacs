New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

A feature-limited version of the PLUMED interface is available
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A basic version of the `PLUMED <https://www.plumed.org/>`_ interface is now bundled by default in |Gromacs|.
With a non-Windows installation, it is possible to invoke PLUMED directly from the command line using the 
``-plumed`` option of the ``gmx mdrun`` command, followed by the path to a PLUMED input file.
This can be done **without the need to apply a patch** as in previous |Gromacs| versions.
Importantly, this interface is not feature complete, see :ref:`the section in the manual <plumed>` for the details.
