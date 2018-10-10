Bugs fixed
^^^^^^^^^^

Fix type bug in trilinic DD code
""""""""""""""""""""""""""""""""""""""""""""""""""

Fix bug with unusual off-diagonal elements communicating too few atoms.

Ensure domains are large enough for atom motion
""""""""""""""""""""""""""""""""""""""""""""""""""

Domain decomposition now makes sure that domains will always be large
enough so that atoms will not move across additional domains.

:issue:`2614`

Fix chainsep behaviour of pdb2gmx
""""""""""""""""""""""""""""""""""""""""""""""""""

:issue:`2577`
