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

