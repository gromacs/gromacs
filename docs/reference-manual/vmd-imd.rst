
.. this really does not belong here in my opinion

Using VMD plug-ins for trajectory file I/O
------------------------------------------

|Gromacs|
tools are able to use the plug-ins found in an existing installation of
`VMD`_ in order to read and write
trajectory files in formats that are not native to |Gromacs|. You will be
able to supply an AMBER DCD-format trajectory filename directly to
|Gromacs| tools, for example.

This requires a VMD installation not older than version 1.8, that your
system provides the dlopen function so that programs can determine at
run time what plug-ins exist, and that you build shared libraries when
building |Gromacs|. CMake will find the vmd executable in your path, and
from it, or the environment variable ``VMDDIR`` at
configuration or run time, locate the plug-ins. Alternatively, the
``VMD_PLUGIN_PATH`` can be used at run time to specify a
path where these plug-ins can be found. Note that these plug-ins are in
a binary format, and that format must match the architecture of the
machine attempting to use them.

Interactive Molecular Dynamics
------------------------------

|Gromacs| supports the interactive molecular dynamics (IMD) protocol as
implemented by `VMD`_ to control
a running simulation in NAMD. IMD allows to monitor a running |Gromacs|
simulation from a VMD client. In addition, the user can interact with
the simulation by pulling on atoms, residues or fragments with a mouse
or a force-feedback device. Additional information about the |Gromacs|
implementation and an exemplary |Gromacs| IMD system can be found `on this
homepage <http://www.mpibpc.mpg.de/grubmueller/interactivemd>`__.

Simulation input preparation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The |Gromacs| implementation allows transmission and interaction with a
part of the running simulation only, e.g. in cases where no water
molecules should be transmitted or pulled. The group is specified via
the :ref:`mdp` option ``IMD-group``. When
``IMD-group`` is empty, the IMD protocol is disabled and
cannot be enabled via the switches in :ref:`mdrun <gmx mdrun>`. To interact
with the entire system, ``IMD-group`` can be set to
``System``. When using :ref:`grompp <gmx grompp>`, a
:ref:`gro` file to be used as VMD input is written out
(``-imd`` switch of :ref:`grompp <gmx grompp>`).

Starting the simulation
^^^^^^^^^^^^^^^^^^^^^^^

Communication between VMD and |Gromacs| is achieved via TCP sockets and
thus enables controlling an :ref:`mdrun <gmx mdrun>` running locally or on
a remote cluster. The port for the connection can be specified with the
``-imdport`` switch of :ref:`mdrun <gmx mdrun>`, 8888 is the
default. If a port number of 0 or smaller is provided, |Gromacs|
automatically assigns a free port to use with IMD.

Every :math:`N` steps, the :ref:`mdrun <gmx mdrun>` client receives the
applied forces from VMD and sends the new positions to the client. VMD
permits increasing or decreasing the communication frequency
interactively. By default, the simulation starts and runs even if no IMD
client is connected. This behavior is changed by the
``-imdwait`` switch of :ref:`mdrun <gmx mdrun>`. After startup
and whenever the client has disconnected, the integration stops until
reconnection of the client. When the ``-imdterm`` switch is
used, the simulation can be terminated by pressing the stop button in
VMD. This is disabled by default. Finally, to allow interacting with the
simulation (i.e. pulling from VMD) the ``-imdpull`` switch
has to be used. Therefore, a simulation can only be monitored but not
influenced from the VMD client when none of ``-imdwait``,
``-imdterm`` or ``-imdpull`` are set. However,
since the IMD protocol requires no authentication, it is not advisable
to run simulations on a host directly reachable from an insecure
environment. Secure shell forwarding of TCP can be used to connect to
running simulations not directly reachable from the interacting host.
Note that the IMD command line switches of :ref:`mdrun <gmx mdrun>` are
hidden by default and show up in the help text only with
:ref:`gmx mdrun` ``-h -hidden``.

Connecting from VMD
^^^^^^^^^^^^^^^^^^^

In VMD, first the structure corresponding to the IMD group has to be
loaded (*File* :math:`\rightarrow` *New Molecule*). Then the IMD
connection window has to be used (*Extensions* :math:`\rightarrow`
*Simulation* :math:`\rightarrow` *IMD Connect (NAMD)*). In the IMD
connection window, hostname and port have to be specified and followed
by pressing *Connect*. *Detach Sim* allows disconnecting without
terminating the simulation, while *Stop Sim* ends the simulation on the
next neighbor searching step (if allowed by ``-imdterm``).

The timestep transfer rate allows adjusting the communication frequency
between simulation and IMD client. Setting the keep rate loads every
:math:`N^\mathrm{th}` frame into VMD instead of discarding them when a
new one is received. The displayed energies are in SI units in contrast
to energies displayed from NAMD simulations.s
