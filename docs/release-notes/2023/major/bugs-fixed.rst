Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Density-guided simulation normalization
"""""""""""""""""""""""""""""""""""""""

With the ``.mdp`` option ``density-guided-simulation-normalize-densities = yes``
, the reference density and the simulated density values were previously divided
by the sum of their values.

This lead to surprising behavior for reference densities with lots of negative
voxel values: the density started to "repel" the protein structure
instead of attracting it, if the total sum of voxel values was smaller
than zero. The negative normalization constant lead to a sign change in voxel
values.

To avoid this behavior, the reference density is now normalized so that the
sum of *positive* values is unity, ensuring that the normalization constant is
always positive.

Apart from avoiding the unexpected behavior, we expect that this also leads 
to smaller absolute differences between reference density and simulated density,
with some small benefits for numerical stability.

This change affects all simulations where voxel values are negative
(usually this excludes synthetic data) and that are run with
``density-guided-simulation-normalize-densities = yes``, but only has a larger
effect for: first, similarity  measure ``inner-product`` as an effective
force-constant scaling and, second, for all similarity measures where the sum
of all voxel values was negative.   
