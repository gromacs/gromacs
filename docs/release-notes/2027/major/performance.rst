Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!


Improved performance of systems that are inhomogeneous along z
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The non-bonded kernels got very inefficient with systems that are have
a very inhomogeneous atom distribution along the z-dimension. Now
all three dimensions are considered when estimating atom density.

:issue:`5622`


Reduced host-side blocking in the GPU update and constraints setup
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The mass and temperature-coupling group data that the GPU update and
constraints upload on every neighbour-search step is now held in pinned
host memory and copied asynchronously, rather than with blocking copies.
This cuts the host-side cost of the GPU update and constraints setup by
roughly a factor of 2.6, with the largest benefit for runs with frequent
pair-list updates.
