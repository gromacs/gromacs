New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

|Gromacs| build is now more reproducible
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The build system no longer embeds information about who built the
binary and where.  We used to include this information to help
troubleshoot problems and ensure checkpoint continuations are exact
where possible, but this does not seem necessary. This makes the build
closer to ``reproducible by default`` which is useful for projects
that offer distributions of reproducible software, including
|Gromacs|.

Update gmx cluster to write correct PDB files and index files with cluster frames
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

:ref:`PDB <pdb>` files from gmx cluster were missing the CRYST header for box information, making
it more difficult than needed to use them with our |Gromacs| tools. Also, the :ref:`index <ndx>`
files needed for :ref:`gmx trjconv` to split up trajectories into frames corresponding
to the clusters were not written. This adds support for writing this :ref:`index <ndx>` file
as well as proper :ref:`PDB <pdb>` files.

Allow using COM of previous step as PBC reference
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Added an option (``pull-pbc-ref-from-prev-step-com``), when pulling, to use
the COM of the group of the previous step, to calculate PBC jumps instead of a
reference atom, which can sometimes move a lot during the simulation.
With this option the PBC reference atom is only used at initialization.
This can be of use when using large pull groups or groups with potentially
large relative movement of atoms.
