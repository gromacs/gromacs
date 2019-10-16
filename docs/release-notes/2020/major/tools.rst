Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on redmine, without the
   a space between the colon and number!

Fixed bug in gmx order -calcdist
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The reference position for the distance calculation was calculated
wrongly.

Improved grompp usability by rejecting more invalid .mdp lines
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Lines like

    ref-t 298
    = 0.1
    =

are now all rejected with a descriptive message, which will help
prevent some kinds of errors in constructing .mdp inputs. Note that an
.mdp parameter name with a missing value is still accepted, and leads
to the default behavior for that parameter.

Added convert-trj
"""""""""""""""""""""""""""""""""""""""
A new tool :ref:`convert-trj <gmx convert-trj>` has been added to allow
users to interchange trajectory formats without having to use legacy :ref:`gmx trjconv`.
Supported actions include the generation of slimmed down output trajectories, as well
as the replacement of particle information in individual frames with data from a structure file.
The new tool allows the usage of command line selections, meaning it is no longer
necessary to write :ref:`index <ndx>` files to select certain atoms.
It is part of the drive to split up the :ref:`trjconv <gmx trjconv>` tool
into smaller parts.

Added extract-cluster
"""""""""""""""""""""""""""""""""""""""

Added a dedicated tool to extract trajectory frames corresponding to different clusters obtained
from :ref:`gmx cluster`. The new :ref:`extract-cluster <gmx extract-cluster>` tool
generates new trajectories that contain only those frames that correspond to the correct cluster.
The corresponding option **-sub** in :ref:`gmx trjconv` has been removed.

Changed behaviour of genion
"""""""""""""""""""""""""""

Functionality of genion was altered to prevent swapping ions for solvent closer
than -rmin from any other non-solvent atom.
This improvement prevents situations where an ion could be placed at the core
of a protein, which would potentially render the folded protein less stable or
may require long equilibration times.
