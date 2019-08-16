.. _electric fields:

Electric fields
---------------

A pulsed and oscillating electric field can be applied according to:

.. math:: E(t) = E_0 \exp\left[-\frac{(t-t_0)^2}{2\sigma^2}\right]\cos\left[\omega (t-t_0)\right]
          :label: eq-efield

where :math:`E_0` is the field strength, the angular frequency
:math:`\omega = 2\pi c/\lambda`, :math:`t_0` is the time
at of the peak in the field strength and :math:`\sigma` is the width of
the pulse. Special cases occur when :math:`\sigma` = 0 (non-pulsed
field) and for :math:`\omega` is 0 (static field). See
:mdp:`electric-field-x` for more details.

This simulated laser-pulse was applied to simulations of melting
ice \ :ref:`146 <refCaleman2008a>`. A pulsed electric field may look ike
:numref:`Fig. %s <fig-field>`. In the supporting information of that paper the impact
of an applied electric field on a system under periodic boundary
conditions is analyzed. It is described that the effective electric
field under PBC is larger than the applied field, by a factor depending
on the size of the box and the dielectric properties of molecules in the
box. For a system with static dielectric properties this factor can be
corrected for. But for a system where the dielectric varies over time,
for example a membrane protein with a pore that opens and closes during
the simulation, this way of applying an electric field is not useful.
In such cases one can use the computational electrophysiology protocol
described in the next section (sec. :ref:`compel`).

.. _fig-field:

.. figure:: plots/field.*
   :width: 8.00000cm

   A simulated laser pulse in |Gromacs|.

Electric fields are applied when the following options are specified in
the :ref:`grompp <gmx grompp>` :ref:`mdp` file. You specify, in order, :math:`E_0`,
:math:`\omega`, :math:`t_0` and :math:`\sigma`:

::

    electric-field-x = 0.04 0       0     0

yields a static field with :math:`E_0` = 0.04 V/nm in the X-direction.
In contrast,

::

    electric-field-x = 2.0  150     5     0

yields an oscillating electric field with :math:`E_0` = 2 V/nm,
:math:`\omega` = 150/ps and :math:`t_0` = 5 ps. Finally

::

    electric-field-x = 2.0  150     5     1

yields an pulsed-oscillating electric field with :math:`E_0` = 2 V/nm,
:math:`\omega` = 150/ps and :math:`t_0` = 5 ps and :math:`\sigma` = 1
ps. Read more in ref. \ :ref:`146 <refCaleman2008a>`. Note that the input file
format is changed from the undocumented older version. A figure like
:numref:`Fig. %s <fig-field>` may be produced by passing the
``-field`` option to :ref:`gmx mdrun`.
