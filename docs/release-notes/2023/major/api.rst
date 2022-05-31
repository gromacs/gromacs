Changes to the API
^^^^^^^^^^^^^^^^^^

Legacy aggregating headers have been removed.
"""""""""""""""""""""""""""""""""""""""""""""

Previously, some of the legacy API headers existed _only_ to aggregate ``#include``
lines for other installed headers.
No guidance was provided regarding which header to include for a given feature.
These redundant headers have been removed.
Client software relying on ``#include "gromacs/module.h"`` will need to be updated
with more specific ``#include "gromacs/module/feature.h"`` directives.

:issue:`4487`

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

