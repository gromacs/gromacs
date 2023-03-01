Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Supported replacing solvent in ``gmx insert-molecules``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Make it possible to specify the solvent (or other set of atoms) with
``-replace`` (as a selection) for ``gmx insert-molecules``, and make the tool
replace residues from this set with the inserted molecules, instead of
not inserting there. It is assumed that the solvent consists of
single-residue molecules, since molecule information would require a tpr
input, which might not be commonly available when preparing the system.

Default random seeds have changed for some analysis tools
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
See individual tools documentation for their functionality. In some
cases, the magic value to obtain a generated seed has changed (or is
now documented.)

Made ``gmx solvate`` and ``gmx insert-molecules`` work better with PDB inputs
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
When both ``-f`` and ``-o`` were .pdb files, the pdbinfo struct got
out-of-sync when the atoms were added/removed.

:issue:`1887`

Tools in the new analysis framework can read trajectory files with subsets
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Make tools written for the new C++ analysis framework support analyzing
trajectories that contain an arbitrary subset of atoms.

:issue:`1861`

Made moleculetype name case sensitive
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This is useful in case you have more than 36 chains in your system
with chain IDs set. PDB allows using both uppercase letters, lowercase
letters and numbers for chain identifiers. Now we can use the maximum
of 62 chains.

Added number density normalization option for ``gmx rdf``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Add an option to ``gmx rdf`` that allows selecting a radial number density
as the normalization for the output (in addition to current raw
neighbor counts and the actual RDF).

Simplified ``gmx genconf`` by removing ``-block``, ``-sort`` and ``-shuffle``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Option ``-block`` isn't useful since particle decomposition was removed.
Options ``-sort`` and ``-shuffle`` were undocumented and don't seem very
useful - these days they would be somebody's simple python script.

Used macros for units and conversions in ``gmx wham``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Also :issue:`1841`

Improved ``gmx sasa`` error message
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Print more information when an output group is not part of the group
selected for calculation, which should help the user diagnosing the issue.

Made ``gmx vanhove`` work without PBC
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Fix ``gmx hbond`` group overlap check
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
``gmx hbond`` does not support partially overlapping analysis groups.
The check in the code was broken and never caught this, resulting
incorrect output that might OK at first sight.
Also corrected bitmasks = enums that (intentionally?) seemed to give
correct results by not using non power of 2 enum index entries.

Made ``gmx dos`` work again.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Due to an error in the index handling ``gmx dos`` always stopped with a fatal
error.

:issue:`1996`

Add checks for too much memory in ``gmx nmeig``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
``gmx nmeig`` could request storage for eigenvector output and matrices
for more than ``INT_MAX`` elements, but nearly all loop variables are int.
Now a fatal error is produced in this case. This also avoids the
confusing error message when too much memory is requested; the allocation
routine will get the correct size, but gmx_fatal prints it as a smaller
integer.
Added support for ``-first`` > 1 with sparse matrices.

