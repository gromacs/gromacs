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
Client code may be built against a |Gromacs| installation (CMake support through
`find_package(gmxapi)` after sourcing `GMXRC`).
MD plugin code may apply externally calculated forces (see restraint module) or
issue simulation stop signals through session resources available at run time
to registered plugins.
In supported environments,
CMake build option `-DGMXAPI=ON` (default) causes headers to be
installed in a `gmxapi` directory next to the previously available `gromacs`
headers directory.
Build targets `gmxapi_cppdocs` and `gmxapi_cppdocs_dev` produce documentation in
`docs/api-user` and `docs/api-dev`, respectively.
For more project information and use cases,
refer to the tracked :issue:`#2585`,
associated GitHub `gmxapi`_ projects,
or DOI `10.1093/bioinformatics/bty484 <gmxapiDOI>`_

Restraint module
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Provides functionality that was previously accessed by modifying the "pull" code in the
|Gromacs| source.
Client software may be built against an unmodified |Gromacs| installation.
Separately compiled MD extensions can be registered with the new Restraint
functionality at run time using simulation client code built with the new `gmxapi` tools.
For more project information and use cases,
refer to the tracked :issue:`#2585`,
associated GitHub `gmxapi`_ projects,
or DOI `10.1093/bioinformatics/bty484 <gmxapiDOI>`_

.. _gmxapi: https://github.com/kassonlab/gmxapi
.. _gmxapiDOI: https://doi.org/10.1093/bioinformatics/bty484
