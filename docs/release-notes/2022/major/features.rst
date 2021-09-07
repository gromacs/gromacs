New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!


Transformation pull coordinate for mathematical transformations of pull coordinates
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A new pull coordinate type named transformation has been added. This enables mathematical
transformation of previously defined pull coordinates using a user supplied formula
in a string. This allows for example non-linear transformation of a distance, e.g.
a contact coordinate or (non-)linear combinations of multiple pull coordinates.
This is a powerful tool for defining complex reaction coordinates and it can be combined
with the Accelerated Weight Histogram Method to enhance sampling.

Replica-exchange molecular dynamics simulations with GPU update
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Replica-exchange molecular dynamics now works with GPU update.
