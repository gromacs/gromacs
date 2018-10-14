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

Transitional external API headers and library
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Library access to |Gromacs| is transitioning to new infrastructure.
gmxapi 0.0.7 provides abstractions for execution environment and simulation work,
as well as development tools for extending MD simulation code without patching
the |Gromacs| source.
Client code may be built against a |Gromacs| installation.
MD plugin code may apply externally calculated forces (see restraint module) or
issue simulation stop signals through session resources available at run time
to registered plugins.
For more project information and use cases,
refer to the tracked :issue:`2585` and to
DOI `10.1093/bioinformatics/bty484 <https://doi.org/10.1093/bioinformatics/bty484>`_.
For a few examples of building on and extending |Gromacs|, refer to the
`Python package <https://github.com/kassonlab/gmxapi>`_ and sample
`restraint plugin <https://github.com/kassonlab/sample_restraint>`_ repository.

Restraint module for gmxapi MD extension code
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Provides functionality that was previously accessed by modifying the "pull" code in the
|Gromacs| source.
Client software may be built against an unmodified |Gromacs| installation.
Separately compiled MD extensions can be registered with the new Restraint
functionality at run time using simulation client code built with the new ``gmxapi`` tools.
(See above.)

Enable output of average pull forces and positions
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Normally the pull module writes instantaneous output of positions and forces, however
now it can write the average of these values over the period since the last output.
This works correctly even if a checkpoint restart intervened. This is enabled using the
new options ``pull-fout-average`` and ``pull-xout-average``.
