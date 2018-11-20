Polarization
------------

Polarization can be treated by |Gromacs| by attaching shell (Drude)
particles to atoms and/or virtual sites. The energy of the shell
particle is then minimized at each time step in order to remain on the
Born-Oppenheimer surface.

Simple polarization
~~~~~~~~~~~~~~~~~~~

This is implemented as a harmonic potential with equilibrium distance 0.
The input given in the topology file is the polarizability
:math:`\alpha` (in |Gromacs| units) as follows:

::

    [ polarization ]
    ; Atom i  j  type  alpha
    1         2  1     0.001

in this case the polarizability volume is 0.001 nm\ :math:`^3` (or 1
Å\ :math:`^3`). In order to compute the harmonic force constant
:math:`k_{cs}` (where :math:`cs` stands for core-shell), the following
is used \ :ref:`45 <refMaaren2001a>`:

.. math:: k_{cs} ~=~ \frac{q_s^2}{\alpha}
          :label: eqnsimplepol

where :math:`q_s` is the charge on the shell particle.

Anharmonic polarization
~~~~~~~~~~~~~~~~~~~~~~~

For the development of the Drude force field by Roux and
McKerell \ :ref:`93 <refLopes2013a>` it was found that some particles can
overpolarize and this was fixed by introducing a higher order term in
the polarization energy:

.. math:: \begin{aligned}
          V_{pol} ~=& \frac{k_{cs}}{2} r_{cs}^2 & r_{cs} \le \delta \\
                      =& \frac{k_{cs}}{2} r_{cs}^2 + k_{hyp} (r_{cs}-\delta)^4 & r_{cs} > \delta\end{aligned}
          :label: eqnanharmpol

where :math:`\delta` is a user-defined constant that is set to 0.02 nm
for anions in the Drude force field \ :ref:`94 <refHYu2010>`. Since this
original introduction it has also been used in other atom
types \ :ref:`93 <refLopes2013a>`.

::

    [ polarization ]
    ;Atom i j    type   alpha (nm^3)    delta  khyp
    1       2       2       0.001786     0.02  16.736e8

The above force constant :math:`k_{hyp}` corresponds to
4\ :math:`\cdot`\ 10\ :math:`^8` kcal/mol/nm\ :math:`^4`, hence the
strange number.

Water polarization
~~~~~~~~~~~~~~~~~~

A special potential for water that allows anisotropic polarization of a
single shell particle \ :ref:`45 <refMaaren2001a>`.

Thole polarization
~~~~~~~~~~~~~~~~~~

Based on early work by Thole :ref:`95 <refThole81>`, Roux and coworkers
have implemented potentials for molecules like
ethanol \ :ref:`96 <refLamoureux2003a>`\ :ref:`98 <refNoskov2005a>`.
Within such molecules, there are intra-molecular interactions between
shell particles, however these must be screened because full Coulomb
would be too strong. The potential between two shell particles :math:`i`
and :math:`j` is:

.. math:: V_{thole} ~=~ \frac{q_i q_j}{r_{ij}}\left[1-\left(1+\frac{{\bar{r}_{ij}}}{2}\right){\rm exp}^{-{\bar{r}_{ij}}}\right]
          :label: eqntholepol

**Note** that there is a sign error in Equation 1 of Noskov
*et al.*  :ref:`98 <refNoskov2005a>`:

.. math:: {\bar{r}_{ij}}~=~ a\frac{r_{ij}}{(\alpha_i \alpha_j)^{1/6}}
          :label: eqntholsignerror

where :math:`a` is a magic (dimensionless) constant, usually chosen to
be 2.6 \ :ref:`98 <refNoskov2005a>`; :math:`\alpha_i` and
:math:`\alpha_j` are the polarizabilities of the respective shell
particles.
