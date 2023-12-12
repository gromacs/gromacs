.. _colvars:

Collective Variable simulations with the Colvars module
-------------------------------------------------------

The Colvars module enables on-the-fly computation of low-dimension quantities (collective
variables or colvars) in simulations, the application of external forces (biases) to these
colvars for restraining or enhanced sampling purposes, and the computation of free energy
profiles and other properties.
The Colvars library and module are described in ref. \ :ref:`195 <refFiorin13>`
as well as in other references that are reported in the log file when the corresponding features
are used.

Using Colvars
^^^^^^^^^^^^^

Colvars simulations are enabled by the following :ref:`mdp` file options: ``colvars-active``,
``colvars-configfile``, and  ``colvars-seed``.

Setting :mdp-value:`colvars-active = true` enables Colvars, using a configuration that can be
defined by specifying a Colvars configuration file using :mdp-value:`colvars-configfile`.

See :ref:`this section of the documentation <mdp-colvars>` for detailed usage of these options.


Configuration files for input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Colvars configuration files are plain text files describing specific collective variables and
biasing and analysis algorithms to be applied onto them.
Full documentation is available `here <https://colvars.github.io/gromacs-2024/colvars-refman-gromacs.html>`_.

Additionally, the Colvars Dashboard extension within VMD can be used to prepare a Colvars
configuration file, leveraging input templates for many features; VMD version 1.9.4 or later is
strongly recommended.


Colvars output files
^^^^^^^^^^^^^^^^^^^^

When Colvars is active, additional output files are written during a
|Gromacs| simulation.  Their file names share the same prefix as the
:ref:`edr` file, which is ``ener`` by default, followed by a suffix specific
to their content (e.g. ``ener.colvars.traj`` for the trajectory of the
collective variables). These files are useful for analysis purposes, but
they are not required for continuing a simulation, because the relevant data
from Colvars is included in the binary checkpoint file.


Colvars checkpointing
^^^^^^^^^^^^^^^^^^^^^

The state of the Colvars library is written to the checkpoint file and read automatically upon restarting.
