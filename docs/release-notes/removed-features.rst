Removed features
^^^^^^^^^^^^^^^^

Removed hybrid GPU+CPU nonbonded mode
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This mode was not very useful, since it ran the non-local non-bonded
interactions on the CPU. The fraction of non-local interaction is set
by the domain decomposition, so this is not flexible.
Also this mode was not being tested.

QM/MM: removed optimization and transition-state search
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
These functionalities used to only work with old versions of Orca,
had very limited use and will possibly not work any longer now.

Updated application clock handling on Pascal+ GPUs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Starting with Pascal (CC >= 6.0) it is no longer possible to change
application clocks without root privileges. Application
clocks are still reported for Pascal+, but there is no longer
suggestions about changing them.

Removed continuation from :ref:`gmx convert-tpr`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Removed the obsolete option of :ref:`gmx convert-tpr` to write a tpr
file for continuation using a trajectory and energy file. This is
superseded by checkpointing.
