.. _analysis:

Analysis
========

In this chapter different ways of analyzing your trajectory are
described. The names of the corresponding analysis programs are given.
Specific information on the in- and output of these programs can be
found in the tool documentation :ref:`here <gmx-cmdline>`.
The output files are often
produced as finished Grace/Xmgr graphs.

First, in sec. :ref:`usinggroups`, the group concept in analysis is
explained. :ref:`selections` explains a newer concept of dynamic
selections, which is currently supported by a few tools. Then, the
different analysis tools are presented.

.. _usinggroups:

Using Groups
------------

| In chapter :ref:`algorithms`, it was explained how *groups of atoms* can
  be used in mdrun (see sec. :ref:`groupconcept`). In most analysis
  programs, groups of atoms must also be chosen. Most programs can
  generate several default index groups, but groups can always be read
  from an index file. Let’s consider the example of a simulation of a
  binary mixture of components A and B. When we want to calculate the
  radial distribution function (RDF) :math:`g_{AB}(r)` of A with respect
  to B, we have to calculate:

  .. math:: 4\pi r^2 g_{AB}(r)      ~=~     V~\sum_{i \in A}^{N_A} \sum_{j \in B}^{N_B} P(r)

| where :math:`V` is the volume and :math:`P(r)` is the probability of
  finding a B atom at distance :math:`r` from an A atom.

By having the user define the *atom numbers* for groups A and B in a
simple file, we can calculate this :math:`g_{AB}` in the most general
way, without having to make any assumptions in the RDF program about the
type of particles.

Groups can therefore consist of a series of *atom numbers*, but in some
cases also of *molecule numbers*. It is also possible to specify a
series of angles by *triples* of *atom numbers*, dihedrals by
*quadruples* of *atom numbers* and bonds or vectors (in a molecule) by
*pairs* of *atom numbers*. When appropriate the type of index file will
be specified for the following analysis programs. To help creating such
:ref:`index file <ndx>` ``index.ndx``), there are a couple of programs to generate
them, using either your input configuration or the topology. To generate
an index file consisting of a series of *atom numbers* (as in the
example of :math:`g_{AB}`), use :ref:`gmx make_ndx`
or :ref:`gmx select`. To generate an index file with
angles or dihedrals, use :ref:`gmx mk_angndx`. Of course you can also
make them by hand. The general format is presented here:

::

    [ Oxygen ]
       1       4       7

    [ Hydrogen ]
       2       3       5       6
       8       9

First, the group name is written between square brackets. The following
atom numbers may be spread out over as many lines as you like. The atom
numbering starts at 1.

Each tool that can use groups will offer the available alternatives for
the user to choose. That choice can be made with the number of the
group, or its name. In fact, the first few letters of the group name
will suffice if that will distinguish the group from all others. There
are ways to use Unix shell features to choose group names on the command
line, rather than interactively. Consult our `webpage`_ for suggestions.

.. _defaultgroups:

Default Groups
~~~~~~~~~~~~~~

When no index file is supplied to analysis tools or
:ref:`grompp <gmx grompp>`, a number of default
groups are generated to choose from:

``System``
    | all atoms in the system

``Protein``
    | all protein atoms

``Protein-H``
    | protein atoms excluding hydrogens

``C-alpha``
    | C\ :math:`_{\alpha}` atoms

``Backbone``
    | protein backbone atoms; N, C\ :math:`_{\alpha}` and C

``MainChain``
    | protein main chain atoms: N, C\ :math:`_{\alpha}`, C and O,
      including oxygens in C-terminus

``MainChain+Cb``
    | protein main chain atoms including C\ :math:`_{\beta}`

``MainChain+H``
    | protein main chain atoms including backbone amide hydrogens and
      hydrogens on the N-terminus

``SideChain``
    | protein side chain atoms; that is all atoms except N,
      C\ :math:`_{\alpha}`, C, O, backbone amide hydrogens, oxygens in
      C-terminus and hydrogens on the N-terminus

``SideChain-H``
    | protein side chain atoms excluding all hydrogens

``Prot-Masses``
    | protein atoms excluding dummy masses (as used in virtual site
      constructions of NH\ :math:`_3` groups and tryptophan
      side-chains), see also sec. :ref:`vsitetop`; this group is only
      included when it differs from the ``Protein`` group

``Non-Protein``
    | all non-protein atoms

``DNA``
    | all DNA atoms

``RNA``
    | all RNA atoms

``Water``
    | water molecules (names like ``SOL``, ``WAT``, ``HOH``, etc.) See
      ``residuetypes.dat`` for a full listing

``non-Water``
    | anything not covered by the ``Water`` group

``Ion``
    | any name matching an Ion entry in
      ``residuetypes.dat``

``Water_and_Ions``
    | combination of the ``Water`` and ``Ions``
      groups

``molecule_name``
    | for all residues/molecules which are not recognized as protein,
      DNA, or RNA; one group per residue/molecule name is generated

``Other``
    | all atoms which are neither protein, DNA, nor RNA.

Empty groups will not be generated. Most of the groups only contain
protein atoms. An atom is considered a protein atom if its residue name
is listed in the
``residuetypes.dat``
file and is listed as a “Protein” entry. The process for determinding
DNA, RNA, etc. is analogous. If you need to modify these
classifications, then you can copy the file from the library directory
into your working directory and edit the local copy.

.. _selections:

Selections
~~~~~~~~~~

| :ref:`gmx select <gmx select>`
| Currently, a few analysis tools support an extended concept of
  *(dynamic) selections*. There are three
  main differences to traditional index groups:

-  The selections are specified as text instead of reading fixed atom
   indices from a file, using a syntax similar to VMD. The text can be
   entered interactively, provided on the command line, or from a file.

-  The selections are not restricted to atoms, but can also specify that
   the analysis is to be performed on, e.g., center-of-mass positions of
   a group of atoms. Some tools may not support selections that do not
   evaluate to single atoms, e.g., if they require information that is
   available only for single atoms, like atom names or types.

-  The selections can be dynamic, i.e., evaluate to different atoms for
   different trajectory frames. This allows analyzing only a subset of
   the system that satisfies some geometric criteria.

As an example of a simple selection, ``resname ABC`` and
``within 2 of resname DEF`` selects all atoms in residues named ABC that are
within 2nm of any atom in a residue named DEF.

Tools that accept selections can also use traditional index files
similarly to older tools: it is possible to give an :ref:`ndx`
file to the tool, and directly select a group from the index file as a
selection, either by group number or by group name. The index groups can
also be used as a part of a more complicated selection.

To get started, you can run :ref:`gmx select <gmx select>` with a single
structure, and use the interactive prompt to try out different
selections. The tool provides, among others, output options
``-on`` and ``-ofpdb`` to write out the selected
atoms to an index file and to a :ref:`pdb` file, respectively.
This does not allow testing selections that evaluate to center-of-mass
positions, but other selections can be tested and the result examined.

The detailed syntax and the individual keywords that can be used in
selections can be accessed by typing ``help`` in the
interactive prompt of any selection-enabled tool, as well as with
:ref:`gmx help <gmx help>` selections. The help is divided into subtopics
that can be accessed with, e.g., ``help syntax``/
:ref:`gmx help <gmx help>` ``selections syntax``. Some individual selection
keywords have extended help as well, which can be accessed with, e.g.,
``help keywords`` within.

The interactive prompt does not currently provide much editing
capabilities. If you need them, you can run the program under
``rlwrap``.

For tools that do not yet support the selection syntax, you can use
:ref:`gmx select <gmx select>` -on to generate static index groups to pass
to the tool. However, this only allows for a small subset (only the
first bullet from the above list) of the flexibility that fully
selection-aware tools offer.

It is also possible to write your own analysis tools to take advantage
of the flexibility of these selections: see the
``template.cpp`` file in the
``share/gromacs/template`` directory of your installation
for an example.

Looking at your trajectory
--------------------------

.. _fig-ngmxdump:

.. figure:: plots/ngmxdump.*
   :width: 8.00000cm

   The window of :ref:`gmx view <gmx view>` showing a box of water.

| :ref:`gmx view <gmx view>`
| Before analyzing your trajectory it is often informative to look at
  your trajectory first. |Gromacs| comes with a simple trajectory viewer
  :ref:`gmx view <gmx view>`; the advantage
  with this one is that it does not require OpenGL, which usually isn’t
  present on *e.g.* supercomputers. It is also possible to generate a
  hard-copy in Encapsulated Postscript format (see
  :numref:`Fig. %s <fig-ngmxdump>`). If you want a faster and more
  fancy viewer there are several programs that can read the |Gromacs|
  trajectory formats – have a look at our `webpage`_ for updated links.

General properties
------------------

| :ref:`gmx energy <gmx energy>`, :ref:`gmx traj <gmx traj>`
| To analyze some or all *energies* and other properties, such as *total
  pressure*, *pressure tensor*, *density*, *box-volume* and *box-sizes*,
  use the program :ref:`gmx energy <gmx energy>`. A choice can be made from a
  list a set of energies, like potential, kinetic or total energy, or
  individual contributions, like Lennard-Jones or dihedral energies.

The *center-of-mass velocity*, defined as

.. math:: {\bf v}_{com} = {1 \over M} \sum_{i=1}^N m_i {\bf v}_i

with :math:`M = \sum_{i=1}^N m_i` the total mass of the system, can be
monitored in time by the program :ref:`gmx traj <gmx traj>` ``-com -ov``. It is however
recommended to remove the center-of-mass velocity every step (see
chapter :ref:`algorithms`)!

Radial distribution functions
-----------------------------

| :ref:`gmx rdf <gmx rdf>`
| The *radial distribution function* (RDF) or pair correlation function
  :math:`g_{AB}(r)` between particles of type :math:`A` and :math:`B` is
  defined in the following way:

.. math::

   \begin{array}{rcl}
   g_{AB}(r)&=&    {\displaystyle \frac{\langle \rho_B(r) \rangle}{\langle\rho_B\rangle_{local}}}         \\
            &=&    {\displaystyle \frac{1}{\langle\rho_B\rangle_{local}}}{\displaystyle \frac{1}{N_A}}
                   \sum_{i \in A}^{N_A} \sum_{j \in B}^{N_B} 
                   {\displaystyle \frac{\delta( r_{ij} - r )}{4 \pi r^2}}         \\
   \end{array}

with :math:`\langle\rho_B(r)\rangle` the particle density of type
:math:`B` at a distance :math:`r` around particles :math:`A`, and
:math:`\langle\rho_B\rangle_{local}` the particle density of type
:math:`B` averaged over all spheres around particles :math:`A` with
radius :math:`r_{max}` (see :numref:`Fig. %s <fig-rdfex>` C).

.. _fig-rdfex:

.. figure:: plots/rdf.*
    :width: 7.00000cm

    Definition of slices in :ref:`gmx rdf <gmx rdf>`: A. :math:`g_{AB}(r)`.
    B. :math:`g_{AB}(r,\theta)`. The slices are colored gray. C.
    Normalization :math:`\langle\rho_B\rangle_{local}`. D. Normalization
    :math:`\langle\rho_B\rangle_{local,\:\theta }`. Normalization volumes
    are colored gray.

Usually the value of :math:`r_{max}` is half of the box length. The
averaging is also performed in time. In practice the analysis program
:ref:`gmx rdf <gmx rdf>` divides the system
into spherical slices (from :math:`r` to :math:`r+dr`, see
:numref:`Fig. %s <fig-rdfex>` A) and makes a histogram in stead of
the :math:`\delta`-function. An example of the RDF of oxygen-oxygen in
SPC water \ :ref::ref:`80 <refBerendsen81>` is given in :numref:`Fig. %s <fig-rdf>`

.. _fig-rdf:

.. figure:: plots/rdfO-O.*
    :width: 8.00000cm

    :math:`g_{OO}(r)` for Oxygen-Oxygen of SPC-water.

With :ref:`gmx rdf <gmx rdf>` it is also possible to calculate an angle
dependent rdf :math:`g_{AB}(r,\theta)`, where the angle :math:`\theta`
is defined with respect to a certain laboratory axis :math:`{\bf e}`,
see :numref:`Fig. %s <fig-rdfex>` B.

.. math::

   \begin{aligned}
   g_{AB}(r,\theta) &=& {1 \over \langle\rho_B\rangle_{local,\:\theta }} {1 \over N_A} \sum_{i \in A}^{N_A} \sum_{j \in B}^{N_B} {\delta( r_{ij} - r ) \delta(\theta_{ij} -\theta) \over 2 \pi r^2 sin(\theta)}\\
   cos(\theta_{ij}) &=& {{\bf r}_{ij} \cdot {\bf e} \over \|r_{ij}\| \;\| e\| }\end{aligned}

This :math:`g_{AB}(r,\theta)` is useful for analyzing anisotropic
systems. **Note** that in this case the normalization
:math:`\langle\rho_B\rangle_{local,\:\theta}` is the average density in
all angle slices from :math:`\theta` to :math:`\theta + d\theta` up to
:math:`r_{max}`, so angle dependent, see :numref:`Fig. %s <fig-rdfex>` D.

Correlation functions
---------------------

Theory of correlation functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The theory of correlation functions is well established \ :ref:`108 <refAllen87>`.
We describe here the implementation of the various
correlation function flavors in the |Gromacs| code. The definition of the
autocorrelation function (ACF) :math:`C_f(t)` for a property
:math:`f(t)` is:

.. math:: C_f(t)  ~=~     \left\langle f(\xi) f(\xi+t)\right\rangle_{\xi}
          :label: eqncorr

where the notation on the right hand side indicates averaging over
:math:`\xi`, *i.e.* over time origins. It is also possible to compute
cross-correlation function from two properties :math:`f(t)` and
:math:`g(t)`:

.. math:: C_{fg}(t) ~=~   \left\langle f(\xi) g(\xi+t)\right\rangle_{\xi}

however, in |Gromacs| there is no standard mechanism to do this
(**note:** you can use the ``xmgr`` program to compute cross correlations).
The integral of the correlation function over time is the correlation
time :math:`\tau_f`:

.. math:: \tau_f  ~=~     \int_0^{\infty} C_f(t) {\rm d} t
          :label: eqncorrtime

In practice, correlation functions are calculated based on data points
with discrete time intervals :math:`\Delta`\ t, so that the ACF from an
MD simulation is:

.. math:: C_f(j\Delta t)  ~=~     \frac{1}{N-j}\sum_{i=0}^{N-1-j} f(i\Delta t) f((i+j)\Delta t)
          :label: eqncorrmd

where :math:`N` is the number of available time frames for the
calculation. The resulting ACF is obviously only available at time
points with the same interval :math:`\Delta`\ t. Since, for many
applications, it is necessary to know the short time behavior of the ACF
(*e.g.* the first 10 ps) this often means that we have to save the data
with intervals much shorter than the time scale of interest. Another
implication of :eq:`eqn. %s <eqncorrmd>` is that in principle we can not compute
all points of the ACF with the same accuracy, since we have :math:`N-1`
data points for :math:`C_f(\Delta t)` but only 1 for
:math:`C_f((N-1)\Delta t)`. However, if we decide to compute only an ACF
of length :math:`M\Delta t`, where :math:`M \leq N/2` we can compute all
points with the same statistical accuracy:

.. math:: C_f(j\Delta t)  ~=~ \frac{1}{M}\sum_{i=0}^{N-1-M} f(i\Delta t)f((i+j)\Delta t)

Here of course :math:`j < M`. :math:`M` is sometimes referred to as the
time lag of the correlation function. When we decide to do this, we
intentionally do not use all the available points for very short time
intervals (:math:`j << M`), but it makes it easier to interpret the
results. Another aspect that may not be neglected when computing ACFs
from simulation is that usually the time origins :math:`\xi`
(:eq:`eqn. %s <eqncorr>`) are not statistically independent, which may introduce
a bias in the results. This can be tested using a block-averaging
procedure, where only time origins with a spacing at least the length of
the time lag are included, *e.g.* using :math:`k` time origins with
spacing of :math:`M\Delta t` (where :math:`kM \leq N`):

.. math:: C_f(j\Delta t)  ~=~ \frac{1}{k}\sum_{i=0}^{k-1} f(iM\Delta t)f((iM+j)\Delta t)

However, one needs very long simulations to get good accuracy this way,
because there are many fewer points that contribute to the ACF.

Using FFT for computation of the ACF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The computational cost for calculating an ACF according to
:eq:`eqn. %s <eqncorrmd>` is proportional to :math:`N^2`, which is considerable.
However, this can be improved by using fast Fourier transforms to do the
convolution \ :ref:`108 <refAllen87>`.

Special forms of the ACF
~~~~~~~~~~~~~~~~~~~~~~~~

There are some important varieties on the ACF, *e.g.* the ACF of a
vector :math:`\mathbf{p}`:

.. math:: C_{\mathbf{p}}(t) ~=~       \int_0^{\infty} P_n(\cos\angle\left(\mathbf{p}(\xi),\mathbf{p}(\xi+t)\right) {\rm d} \xi
          :label: eqncorrleg

where :math:`P_n(x)` is the :math:`n^{th}` order Legendre
polynomial. [1]_ Such correlation times can actually be obtained
experimentally using *e.g.* NMR or other relaxation experiments. |Gromacs|
can compute correlations using the 1\ :math:`^{st}` and 2\ :math:`^{nd}`
order Legendre polynomial (:eq:`eqn. %s <eqncorrleg>`). This can also be used
for rotational autocorrelation (:ref:`gmx rotacf`) and dipole autocorrelation
(:ref:`gmx dipoles <gmx dipoles>`).

In order to study torsion angle dynamics, we define a dihedral
autocorrelation function as \ :ref:`159 <refSpoel97a>`:

.. math:: C(t)    ~=~     \left\langle \cos(\theta(\tau)-\theta(\tau+t))\right\rangle_{\tau}
          :label: eqncoenk

**Note** that this is not a product of two functions as is generally
used for correlation functions, but it may be rewritten as the sum of
two products:

.. math:: C(t)    ~=~     \left\langle\cos(\theta(\tau))\cos(\theta(\tau+t))\,+\,\sin(\theta(\tau))\sin(\theta(\tau+t))\right\rangle_{\tau}
          :label: eqncot

Some Applications
~~~~~~~~~~~~~~~~~

The program :ref:`gmx velacc <gmx velacc>`
calculates the *velocity autocorrelation function*.

.. math:: C_{\mathbf{v}} (\tau) ~=~ \langle {\mathbf{v}}_i(\tau) \cdot {\mathbf{v}}_i(0) \rangle_{i \in A}

The self diffusion coefficient can be calculated using the Green-Kubo
relation \ :ref:`108 <refAllen87>`:

.. math:: D_A ~=~ {1\over 3} \int_0^{\infty} \langle {\bf v}_i(t) \cdot {\bf v}_i(0) \rangle_{i \in A} \; dt

which is just the integral of the velocity autocorrelation function.
There is a widely-held belief that the velocity ACF converges faster
than the mean square displacement (sec. :ref:`msd`), which can also be
used for the computation of diffusion constants. However, Allen &
Tildesley \ :ref:`108 <refAllen87>` warn us that the long-time
contribution to the velocity ACF can not be ignored, so care must be
taken.

Another important quantity is the dipole correlation time. The *dipole
correlation function* for particles of type :math:`A` is calculated as
follows by :ref:`gmx dipoles <gmx dipoles>`:

.. math::

   C_{\mu} (\tau) ~=~
   \langle {\bf \mu}_i(\tau) \cdot {\bf \mu}_i(0) \rangle_{i \in A}

with :math:`{\bf \mu}_i = \sum_{j \in i} {\bf r}_j q_j`. The dipole
correlation time can be computed using :eq:`eqn. %s <eqncorrtime>`. 
For some applications
see (**???**).

The viscosity of a liquid can be related to the correlation time of the
Pressure tensor
:math:`\mathbf{P}` :ref:`160 <refPSmith93c>`,
:ref:`161 <refBalasubramanian96>`. :ref:`gmx energy` can compute the viscosity,
but this is not very accurate \ :ref:`149 <refHess2002a>`, and actually
the values do not converge.

Curve fitting in |Gromacs|
--------------------------

Sum of exponential functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes it is useful to fit a curve to an analytical function, for
example in the case of autocorrelation functions with noisy tails.
|Gromacs| is not a general purpose curve-fitting tool however and
therefore |Gromacs| only supports a limited number of functions.
:numref:`Table %s <table-fitfn>`  lists the available options with the corresponding
command-line options. The underlying routines for fitting use the
Levenberg-Marquardt algorithm as implemented in the lmfit package \ :ref:`162 <reflmfit>`
(a bare-bones version of which is included in |Gromacs| in which an
option for error-weighted fitting was implemented).

.. |exp|  replace:: :math:`e^{-t/{a_0}}`                                                       
.. |aexp| replace:: :math:`a_1e^{-t/{a_0}}`                                                    
.. |exp2| replace:: :math:`a_1e^{-t/{a_0}}+(1-a_1)e^{-t/{a_2}}`                                
.. |exp5| replace:: :math:`a_1e^{-t/{a_0}}+a_3e^{-t/{a_2}}+a_4`                                
.. |exp7| replace:: :math:`a_1e^{-t/{a_0}}+a_3e^{-t/{a_2}}+a_5e^{-t/{a_4}}+a_6`                
.. |exp9| replace:: :math:`a_1e^{-t/{a_0}}+a_3e^{-t/{a_2}}+a_5e^{-t/{a_4}}+a_7e^{-t/{a_6}}+a_8`
.. |nexp2| replace:: :math:`a_2\ge a_0\ge 0`               
.. |nexp5| replace:: :math:`a_2\ge a_0\ge 0`               
.. |nexp7| replace:: :math:`a_4\ge a_2\ge a_0 \ge0`        
.. |nexp9| replace:: :math:`a_6\ge a_4\ge a_2\ge a_0\ge 0` 

.. _table-fitfn:

.. table:: Overview of fitting functions supported in (most) analysis tools 
    that compute autocorrelation functions. The **Note** column describes 
    properties of the output parameters.
    :align: center
    :widths: auto

    +-------------+------------------------------+---------------------+
    | Command     | Functional form :math:`f(t)` | Note                |
    | line option |                              |                     |
    +=============+==============================+=====================+
    | exp         | |exp|                        |                     |
    +-------------+------------------------------+---------------------+
    | aexp        | |aexp|                       |                     |
    +-------------+------------------------------+---------------------+
    | exp_exp     | |exp2|                       | |nexp2|             |
    +-------------+------------------------------+---------------------+
    | exp5        | |exp5|                       | |nexp5|             |
    +-------------+------------------------------+---------------------+
    | exp7        | |exp7|                       | |nexp7|             |
    +-------------+------------------------------+---------------------+
    | exp9        | |exp9|                       | |nexp9|             |
    +-------------+------------------------------+---------------------+


Error estimation
~~~~~~~~~~~~~~~~

Under the hood |Gromacs| implements some more fitting functions, namely a
function to estimate the error in time-correlated data due to Hess \ :ref:`149 <refHess2002a>`:

.. math::

   \varepsilon^2(t) =
   \alpha\tau_1\left(1+\frac{\tau_1}{t}\left(e^{-t/\tau_1}-1\right)\right)
         + (1-\alpha)\tau_2\left(1+\frac{\tau_2}{t}\left(e^{-t/\tau_2}-1\right)\right)

where :math:`\tau_1` and :math:`\tau_2` are time constants (with
:math:`\tau_2 \ge \tau_1`) and :math:`\alpha` usually is close to 1 (in
the fitting procedure it is enforced that :math:`0\leq\alpha\leq 1`).
This is used in :ref:`gmx analyze <gmx analyze>` for error estimation using

.. math:: \lim_{t\rightarrow\infty}\varepsilon(t) = \sigma\sqrt{\frac{2(\alpha\tau_1+(1-\alpha)\tau_2)}{T}}

where :math:`\sigma` is the standard deviation of the data set and
:math:`T` is the total simulation time \ :ref:`149 <refHess2002a>`.

Interphase boundary demarcation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to determine the position and width of an interface,
Steen-Sæthre *et al.* fitted a density profile to the following function

.. math::

   f(x) ~=~ \frac{a_0+a_1}{2} - \frac{a_0-a_1}{2}{\rm
     erf}\left(\frac{x-a_2}{a_3^2}\right)

where :math:`a_0` and :math:`a_1` are densities of different phases,
:math:`x` is the coordinate normal to the interface, :math:`a_2` is the
position of the interface and :math:`a_3` is the width of the
interface \ :ref:`163 <refSteen-Saethre2014a>`. This is implemented
in :ref:`gmx densorder <gmx densorder>`.

Transverse current autocorrelation function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to establish the transverse current autocorrelation function
(useful for computing viscosity  \ :ref:`164 <refPalmer1994a>`) the following function is
fitted:

.. math::

   f(x) ~=~ e^{-\nu}\left({\rm cosh}(\omega\nu)+\frac{{\rm
       sinh}(\omega\nu)}{\omega}\right)

with :math:`\nu = x/(2a_0)` and :math:`\omega = \sqrt{1-a_1}`. This is
implemented in :ref:`gmx tcaf <gmx tcaf>`.

Viscosity estimation from pressure autocorrelation function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The viscosity is a notoriously difficult property to extract from
simulations \ :ref:`149 <refHess2002a>`, :ref:`165 <refWensink2003a>`. It is *in principle*
possible to determine it by integrating the pressure autocorrelation
function \ :ref:`160 <refPSmith93c>`, however this is often hampered by
the noisy tail of the ACF. A workaround to this is fitting the ACF to
the following function \ :ref:`166 <refGuo2002b>`:

.. math::

   f(t)/f(0) = (1-C) {\rm cos}(\omega t) e^{-(t/\tau_f)^{\beta_f}} + C
   e^{-(t/\tau_s)^{\beta_s}}

where :math:`\omega` is the frequency of rapid pressure oscillations
(mainly due to bonded forces in molecular simulations), :math:`\tau_f`
and :math:`\beta_f` are the time constant and exponent of fast
relaxation in a stretched-exponential approximation, :math:`\tau_s` and
:math:`\beta_s` are constants for slow relaxation and :math:`C` is the
pre-factor that determines the weight between fast and slow relaxation.
After a fit, the integral of the function :math:`f(t)` is used to
compute the viscosity:

.. math:: \eta = \frac{V}{k_B T}\int_0^{\infty} f(t) dt

This equation has been applied to computing the bulk and shear
viscosity using different elements from the pressure tensor \ :ref:`167 <refFanourgakis2012a>`.

.. _msd:

Mean Square Displacement
------------------------

| :ref:`gmx msd <gmx msd>`
| To determine the self diffusion
  coefficient :math:`D_A` of
  particles of type :math:`A`, one can use the Einstein
  relation :ref:`108 <refAllen87>`:

  .. math::

     \lim_{t \rightarrow \infty} \langle
     \|{\bf r}_i(t) - {\bf r}_i(0)\|^2 \rangle_{i \in A} ~=~ 6 D_A t

| This *mean square displacement* and :math:`D_A` are calculated by the
  program :ref:`gmx msd <gmx msd>`. Normally
  an index file containing atom numbers is used and the MSD is averaged
  over these atoms. For molecules consisting of more than one atom,
  :math:`{\bf r}_i` can be taken as the center of mass positions of the
  molecules. In that case, you should use an index file with molecule
  numbers. The results will be nearly identical to averaging over atoms,
  however. The :ref:`gmx msd <gmx msd>` program can also be used for
  calculating diffusion in one or two dimensions. This is useful for
  studying lateral diffusion on interfaces.

An example of the mean square displacement of SPC water is given in
:numref:`Fig. %s <fig-msdwater>`.

.. _fig-msdwater:

.. figure:: plots/msdwater.*
    :width: 8.00000cm

    Mean Square Displacement of SPC-water.

.. _bad:

Bonds/distances, angles and dihedrals
-------------------------------------

| :ref:`gmx distance <gmx distance>`, :ref:`gmx angle <gmx angle>`, 
  :ref:`gmx gangle <gmx gangle>`
| To monitor specific *bonds* in your modules, or more generally
  distances between points, the program 
  :ref:`gmx distance <gmx distance>` can calculate distances as a
  function of time, as well as the distribution of the distance. With a
  traditional index file, the groups should consist of pairs of atom
  numbers, for example:

::

    [ bonds_1 ]
     1     2
     3     4
     9    10

    [ bonds_2 ]
    12    13

Selections are also supported, with first two positions defining the
first distance, second pair of positions defining the second distance
and so on. You can calculate the distances between CA and CB atoms in
all your residues (assuming that every residue either has both atoms, or
neither) using a selection such as:

::

    name CA CB

The selections also allow more generic distances to be computed. For
example, to compute the distances between centers of mass of two
residues, you can use:

::

    com of resname AAA plus com of resname BBB

The program :ref:`gmx angle <gmx angle>`
calculates the distribution of *angles* and *dihedrals* in time. It also
gives the average angle or dihedral. The index file consists of triplets
or quadruples of atom numbers:

::

    [ angles ]
     1     2     3
     2     3     4
     3     4     5

    [ dihedrals ]
     1     2     3     4
     2     3     5     5

For the dihedral angles you can use either the “biochemical convention”
(:math:`\phi = 0 \equiv cis`) or “polymer convention”
(:math:`\phi = 0 \equiv trans`), see
:numref:`Fig. %s <fig-dihdef>`.

.. _fig-dihdef:

.. figure:: plots/dih-def.*
    :width: 5.00000cm

    Dihedral conventions: A. “Biochemical convention”. B. “Polymer
    convention”.

The program :ref:`gmx gangle <gmx gangle>`
provides a selection-enabled version to compute angles. This tool can
also compute angles and dihedrals, but does not support all the options
of :ref:`gmx angle <gmx angle>`, such as autocorrelation or other time
series analyses. In addition, it supports angles between two vectors, a
vector and a plane, two planes (defined by 2 or 3 points, respectively),
a vector/plane and the :math:`z` axis, or a vector/plane and the normal
of a sphere (determined by a single position). Also the angle between a
vector/plane compared to its position in the first frame is supported.
For planes, :ref:`gmx gangle <gmx gangle>`
uses the normal vector perpendicular to the plane. See
:numref:`Fig. %s <fig-sgangle>` A, B, C) for the definitions.

.. _fig-sgangle:

.. figure:: plots/sgangle.*
    :width: 3.50000cm

    Angle options of :ref:`gmx gangle <gmx gangle>`: A. Angle between two
    vectors. B. Angle between two planes. C. Angle between a vector and the
    :math:`z` axis. D. Angle between a vector and the normal of a sphere.
    Also other combinations are supported: planes and vectors can be used
    interchangeably.

.. _rg:

Radius of gyration and distances
--------------------------------

| :ref:`gmx gyrate <gmx gyrate>`, :ref:`gmx distance <gmx distance>`, 
  :ref:`gmx mindist <gmx mindist>`, :ref:`gmx mdmat <gmx mdmat>`,
  :ref:`gmx pairdist <gmx pairdist>`, :ref:`gmx xpm2ps <gmx xpm2ps>`
| To have a rough measure for the compactness of a structure, you can
  calculate the *radius of gyration* with the program
  :ref:`gmx gyrate <gmx gyrate>` as follows:

  .. math:: R_g ~=~ \left({\frac{\sum_i \|{\bf r}_i\|^2 m_i}{\sum_i m_i}}\right)^{{\frac{1}{2}}}
            :label: eqnrg

| where :math:`m_i` is the mass of atom :math:`i` and :math:`{\bf r}_i`
  the position of atom :math:`i` with respect to the center of mass of
  the molecule. It is especially useful to characterize polymer
  solutions and proteins. The program will also provide the radius of
  gyration around the coordinate axis (or, optionally, principal axes)
  by only summing the radii components orthogonal to each axis, for
  instance

  .. math:: R_{g,x} ~=~ \left({\frac{\sum_i \left( r_{i,y}^2 + r_{i,z}^2 \right) m_i}{\sum_i m_i}}\right)^{{\frac{1}{2}}}
            :label: eqnrgaxis

Sometimes it is interesting to plot the *distance* between two atoms, or
the *minimum* distance between two groups of atoms (*e.g.*: protein
side-chains in a salt bridge). To calculate these distances between
certain groups there are several possibilities:

*   The *distance between the geometrical centers* of two groups can be
    calculated with the program :ref:`gmx distance <gmx distance>`, as explained in
    sec. :ref:`bad`.

*   The *minimum distance* between two groups of atoms during time can
    be calculated with the program :ref:`gmx mindist <gmx mindist>`. It also calculates the
    *number of contacts* between these groups within a certain radius
    :math:`r_{max}`.

*   :ref:`gmx pairdist <gmx pairdist>` is a selection-enabled version of :ref:`gmx mindist <gmx mindist>`.

*   To monitor the *minimum distances between amino acid residues*
    within a (protein) molecule, you can use the program :ref:`gmx mdmat <gmx mdmat>`. This
    minimum distance between two residues A\ :math:`_i` and
    A\ :math:`_j` is defined as the smallest distance between any pair
    of atoms (i :math:`\in` A\ :math:`_i`, j :math:`\in` A\ :math:`_j`).
    The output is a symmetrical matrix of smallest distances between all
    residues. To visualize this matrix, you can use a program such as
    ``xv``. If you want to view the axes and legend or if you want to print
    the matrix, you can convert it with :ref:`xpm2ps <gmx xpm2ps>` into a Postscript
    :numref:`Fig. %s <fig-distm>`. 

.. _fig-distm:

.. figure:: plots/distm.*
       :width: 6.50000cm

       A minimum distance matrix for a
       peptide \ :ref:`168 <refSpoel96b>`.

*   Plotting these matrices for different time-frames, one can analyze
    changes in the structure, and *e.g.* forming of salt bridges.

.. _rmsd:

Root mean square deviations in structure
----------------------------------------

| :ref:`gmx rms <gmx rms>`, :ref:`gmx rmsdist <gmx rmsdist>`
| The *root mean square deviation* (:math:`RMSD`) of certain atoms in a
  molecule with respect to a reference structure can be calculated with
  the program :ref:`gmx rms <gmx rms>` by least-square fitting the structure to the
  reference structure (:math:`t_2 = 0`) and subsequently calculating the
  :math:`RMSD` (:eq:`eqn. %s <eqnrmsd>`).

  .. math:: RMSD(t_1,t_2) ~=~ \left[\frac{1}{M} \sum_{i=1}^N m_i \|{\bf r}_i(t_1)-{\bf r}_i(t_2)\|^2 \right]^{\frac{1}{2}}
            :label: eqnrmsd

| where :math:`M = \sum_{i=1}^N m_i` and :math:`{\bf r}_i(t)` is the
  position of atom :math:`i` at time :math:`t`. **Note** that fitting
  does not have to use the same atoms as the calculation of the
  :math:`RMSD`; *e.g.* a protein is usually fitted on the backbone atoms
  (N,C:math:`_{\alpha}`,C), but the :math:`RMSD` can be computed of the
  backbone or of the whole protein.

Instead of comparing the structures to the initial structure at time
:math:`t=0` (so for example a crystal structure), one can also calculate
:eq:`eqn. %s <eqnrmsd>` with a structure at time :math:`t_2=t_1-\tau`. This
gives some insight in the mobility as a function of :math:`\tau`. A
matrix can also be made with the :math:`RMSD` as a function of
:math:`t_1` and :math:`t_2`, which gives a nice graphical interpretation
of a trajectory. If there are transitions in a trajectory, they will
clearly show up in such a matrix.

Alternatively the :math:`RMSD` can be computed using a fit-free method
with the program :ref:`gmx rmsdist <gmx rmsdist>`:

.. math:: RMSD(t) ~=~     \left[\frac{1}{N^2}\sum_{i=1}^N \sum_{j=1}^N    \|{\bf r}_{ij}(t)-{\bf r}_{ij}(0)\|^2\right]^{\frac{1}{2}}
          :label: eqnrmsdff

where the *distance* **r**\ :math:`_{ij}` between atoms at time
:math:`t` is compared with the distance between the same atoms at time
:math:`0`.

.. _covanal:

Covariance analysis
-------------------

Covariance analysis, also called principal component analysis or
essential dynamics :ref:`169 <refAmadei93>`\ , can find
correlated motions. It uses the covariance matrix :math:`C` of the
atomic coordinates:

.. math::

   C_{ij} = \left \langle 
   M_{ii}^{\frac{1}{2}} (x_i - \langle x_i \rangle)
   M_{jj}^{\frac{1}{2}}  (x_j - \langle x_j \rangle)
   \right \rangle

where :math:`M` is a diagonal matrix containing the masses of the atoms
(mass-weighted analysis) or the unit matrix (non-mass weighted
analysis). :math:`C` is a symmetric :math:`3N \times 3N` matrix, which
can be diagonalized with an orthonormal transformation matrix :math:`R`:

.. math::

   R^T C R = \mbox{diag}(\lambda_1,\lambda_2,\ldots,\lambda_{3N})
   ~~~~\mbox{where}~~\lambda_1 \geq \lambda_2 \geq \ldots \geq \lambda_{3N}

The columns of :math:`R` are the eigenvectors, also called principal or
essential modes. :math:`R` defines a transformation to a new coordinate
system. The trajectory can be projected on the principal modes to give
the principal components :math:`p_i(t)`:

.. math:: {\bf p}(t) = R^T M^{\frac{1}{2}} ({\bf x}(t) - \langle {\bf x} \rangle)

The eigenvalue :math:`\lambda_i` is the mean square fluctuation of
principal component :math:`i`. The first few principal modes often
describe collective, global motions in the system. The trajectory can be
filtered along one (or more) principal modes. For one principal mode
:math:`i` this goes as follows:

.. math::

   {\bf x}^f(t) =
   \langle {\bf x} \rangle + M^{-\frac{1}{2}} R_{ * i} \, p_i(t)

When the analysis is performed on a macromolecule, one often wants to
remove the overall rotation and translation to look at the internal
motion only. This can be achieved by least square fitting to a reference
structure. Care has to be taken that the reference structure is
representative for the ensemble, since the choice of reference structure
influences the covariance matrix.

One should always check if the principal modes are well defined. If the
first principal component resembles a half cosine and the second
resembles a full cosine, you might be filtering noise (see below). A
good way to check the relevance of the first few principal modes is to
calculate the overlap of the sampling between the first and second half
of the simulation. **Note** that this can only be done when the same
reference structure is used for the two halves.

A good measure for the overlap has been defined in \ :ref:`170 <refHess2002b>`. The
elements of the covariance matrix are proportional to the square of the
displacement, so we need to take the square root of the matrix to
examine the extent of sampling. The square root can be calculated from
the eigenvalues :math:`\lambda_i` and the eigenvectors, which are the
columns of the rotation matrix :math:`R`. For a symmetric and
diagonally-dominant matrix :math:`A` of size :math:`3N \times 3N` the
square root can be calculated as:

.. math::

   A^\frac{1}{2} = 
   R \, \mbox{diag}(\lambda_1^\frac{1}{2},\lambda_2^\frac{1}{2},\ldots,\lambda_{3N}^\frac{1}{2}) \, R^T

It can be verified easily that the product of this matrix with itself
gives :math:`A`. Now we can define a difference :math:`d` between
covariance matrices :math:`A` and :math:`B` as follows:

.. math::

   \begin{aligned}
   d(A,B) & = & \sqrt{\mbox{tr}\left(\left(A^\frac{1}{2} - B^\frac{1}{2}\right)^2\right)
   }
   \\ & = &
   \sqrt{\mbox{tr}\left(A + B - 2 A^\frac{1}{2} B^\frac{1}{2}\right)}
   \\ & = &
   \left( \sum_{i=1}^N \left( \lambda_i^A + \lambda_i^B \right)
   - 2 \sum_{i=1}^N \sum_{j=1}^N \sqrt{\lambda_i^A \lambda_j^B}
   \left(R_i^A \cdot R_j^B\right)^2 \right)^\frac{1}{2}\end{aligned}

where tr is the trace of a matrix. We can now define the overlap
:math:`s` as:

.. math:: s(A,B) = 1 - \frac{d(A,B)}{\sqrt{\mbox{tr}A + \mbox{tr} B}}

The overlap is 1 if and only if matrices :math:`A` and :math:`B` are
identical. It is 0 when the sampled subspaces are completely orthogonal.

A commonly-used measure is the subspace overlap of the first few
eigenvectors of covariance matrices. The overlap of the subspace spanned
by :math:`m` orthonormal vectors :math:`{\bf w}_1,\ldots,{\bf w}_m` with
a reference subspace spanned by :math:`n` orthonormal vectors
:math:`{\bf v}_1,\ldots,{\bf v}_n` can be quantified as follows:

.. math::

   \mbox{overlap}({\bf v},{\bf w}) =
   \frac{1}{n} \sum_{i=1}^n \sum_{j=1}^m ({\bf v}_i \cdot {\bf w}_j)^2

The overlap will increase with increasing :math:`m` and will be 1 when
set :math:`{\bf v}` is a subspace of set :math:`{\bf w}`. The
disadvantage of this method is that it does not take the eigenvalues
into account. All eigenvectors are weighted equally, and when degenerate
subspaces are present (equal eigenvalues), the calculated overlap will
be too low.

Another useful check is the cosine content. It has been proven that the
the principal components of random diffusion are cosines with the number
of periods equal to half the principal component
index \ :ref:`170 <refHess2002b>`, :ref:`171 <refHess2000>`.
The eigenvalues are proportional to the index to the power
:math:`-2`. The cosine content is defined as:

.. math::

   \frac{2}{T}
   \left( \int_0^T \cos\left(\frac{i \pi t}{T}\right) \, p_i(t) \mbox{d} t \right)^2
   \left( \int_0^T p_i^2(t) \mbox{d} t \right)^{-1}

When the cosine content of the first few principal components is close
to 1, the largest fluctuations are not connected with the potential, but
with random diffusion.

The covariance matrix is built and diagonalized by
:ref:`gmx covar <gmx covar>`. The principal components and
overlap (and many more things) can be plotted and analyzed with
:ref:`gmx anaeig <gmx anaeig>`. The cosine
content can be calculated with
:ref:`gmx analyze <gmx analyze>`.

Dihedral principal component analysis
-------------------------------------

| :ref:`gmx angle <gmx angle>`, :ref:`gmx covar <gmx covar>`, 
  :ref:`gmx anaeig <gmx anaeig>`
| Principal component analysis can be performed in dihedral
  space \ :ref:`172 <refMu2005a>` using |Gromacs|. You start by defining the
  dihedral angles of interest in an index file, either using
  :ref:`gmx mk_angndx <gmx mk_angndx>` or otherwise. Then you use the
  :ref:`gmx angle <gmx angle>` program with the ``-or`` flag to
  produce a new :ref:`trr` file containing the cosine and sine
  of each dihedral angle in two coordinates, respectively. That is, in
  the :ref:`trr` file you will have a series of numbers
  corresponding to: cos(\ :math:`\phi_1`), sin(\ :math:`\phi_1`),
  cos(\ :math:`\phi_2`), sin(\ :math:`\phi_2`), ...,
  cos(\ :math:`\phi_n`), sin(\ :math:`\phi_n`), and the array is padded
  with zeros, if necessary. Then you can use this :ref:`trr`
  file as input for the :ref:`gmx covar <gmx covar>` program and perform
  principal component analysis as usual. For this to work you will need
  to generate a reference file (:ref:`tpr`,
  :ref:`gro`, :ref:`pdb` etc.) containing the same
  number of “atoms” as the new :ref:`trr` file, that is for
  :math:`n` dihedrals you need 2\ :math:`n`/3 atoms (rounded up if not
  an integer number). You should use the ``-nofit`` option
  for :ref:`gmx covar <gmx covar>` since the coordinates in the dummy
  reference file do not correspond in any way to the information in the
  :ref:`trr` file. Analysis of the results is done using
  :ref:`gmx anaeig <gmx anaeig>`.

Hydrogen bonds
--------------

| :ref:`gmx hbond <gmx hbond>`
| The program :ref:`gmx hbond <gmx hbond>`
  analyzes the *hydrogen bonds* (H-bonds) between all possible donors D
  and acceptors A. To determine if an H-bond exists, a geometrical
  criterion is used, see also :numref:`Fig. %s <fig-hbond>`:

  .. math::

     \begin{array}{rclcl}
     r       & \leq  & r_{HB}        & = & 0.35~\mbox{nm}    \\
     \alpha  & \leq  & \alpha_{HB}   & = & 30^o              \\
     \end{array}

.. _fig-hbond:

.. figure:: plots/hbond.*
   :width: 2.50000cm

   Geometrical Hydrogen bond criterion.

The value of :math:`r_{HB} = 0.35 \mathrm{nm}` corresponds to the first minimum
of the RDF of SPC water (see also :numref:`Fig. %s <fig-hbondinsert>`).

The program :ref:`gmx hbond <gmx hbond>` analyzes all hydrogen bonds
existing between two groups of atoms (which must be either identical or
non-overlapping) or in specified donor-hydrogen-acceptor triplets, in
the following ways:

.. _fig-hbondinsert:

.. figure:: plots/hbond-insert.*
    :width: 3.50000cm

    Insertion of water into an H-bond. (1) Normal H-bond between two
    residues. (2) H-bonding bridge via a water molecule.

-  Donor-Acceptor distance (:math:`r`) distribution of all H-bonds

-  Hydrogen-Donor-Acceptor angle (:math:`\alpha`) distribution of all
   H-bonds

-  The total number of H-bonds in each time frame

-  The number of H-bonds in time between residues, divided into groups
   :math:`n`-:math:`n`\ +\ :math:`i` where :math:`n` and
   :math:`n`\ +\ :math:`i` stand for residue numbers and :math:`i` goes
   from 0 to 6. The group for :math:`i=6` also includes all H-bonds for
   :math:`i>6`. These groups include the
   :math:`n`-:math:`n`\ +\ :math:`3`, :math:`n`-:math:`n`\ +\ :math:`4`
   and :math:`n`-:math:`n`\ +\ :math:`5` H-bonds, which provide a
   measure for the formation of :math:`\alpha`-helices or
   :math:`\beta`-turns or strands.

-  The lifetime of the H-bonds is calculated from the average over all
   autocorrelation functions of the existence functions (either 0 or 1)
   of all H-bonds:

   .. math:: C(\tau) ~=~ \langle s_i(t)~s_i (t + \tau) \rangle
             :label: eqnhbcorr

-  with :math:`s_i(t) = \{0,1\}` for H-bond :math:`i` at time
   :math:`t`. The integral of :math:`C(\tau)` gives a rough estimate of
   the average H-bond lifetime :math:`\tau_{HB}`:

   .. math::  \tau_{HB} ~=~ \int_{0}^{\infty} C(\tau) d\tau
              :label: eqnhblife

-  Both the integral and the complete autocorrelation function
   :math:`C(\tau)` will be output, so that more sophisticated analysis
   (*e.g.* using multi-exponential fits) can be used to get better
   estimates for :math:`\tau_{HB}`. A more complete analysis is given in
   ref. \ :ref:`173 <refSpoel2006b>`; one of the more fancy option is the Luzar
   and Chandler analysis of hydrogen bond kinetics \ :ref:`174 <refLuzar96b>`, :ref:`175 <refLuzar2000a>`.

-  An H-bond existence map can be generated of dimensions
   *# H-bonds*\ :math:`\times`\ *# frames*. The ordering is identical to
   the index file (see below), but reversed, meaning that the last
   triplet in the index file corresponds to the first row of the
   existence map.

-  Index groups are output containing the analyzed groups, all
   donor-hydrogen atom pairs and acceptor atoms in these groups,
   donor-hydrogen-acceptor triplets involved in hydrogen bonds between
   the analyzed groups and all solvent atoms involved in insertion.

Protein-related items
---------------------

| :ref:`gmx do_dssp <gmx do_dssp>`, :ref:`gmx rama <gmx rama>`,
  :ref:`gmx wheel <gmx wheel>`
| To analyze structural changes of a protein, you can calculate the
  radius of gyration or the minimum residue distances over time (see
  sec. :ref:`rg`), or calculate the RMSD (sec. :ref:`rmsd`).

You can also look at the changing of *secondary structure elements*
during your run. For this, you can use the program 
:ref:`gmx do_dssp <gmx do_dssp>`, which is an interface for the
commercial program ``DSSP``  :ref:`176 <refKabsch83>`. For
further information, see the ``DSSP`` manual. A typical
output plot of :ref:`gmx do_dssp <gmx do_dssp>` is given in
:numref:`Fig. %s <fig-dssp>`.

.. _fig-dssp: 

.. figure:: plots/dssp.*
   :width: 12.00000cm

   Analysis of the secondary structure elements of a peptide in time.

One other important analysis of proteins is the so-called *Ramachandran
plot*. This is the projection of the structure on the two dihedral
angles :math:`\phi` and :math:`\psi` of the protein backbone, see
:numref:`Fig. %s <fig-phipsi>`: 

.. _fig-phipsi:

.. figure:: plots/phipsi.*
   :width: 5.00000cm

   Definition of the dihedral angles :math:`\phi` and :math:`\psi` of
   the protein backbone.

To evaluate this Ramachandran plot you can use the program
:ref:`gmx rama <gmx rama>`. A typical output
is given in :numref:`Fig. %s <fig-rama>`.

.. _fig-rama:

.. figure:: plots/rama.* 
    :width: 5.00000cm

    Ramachandran plot of a small protein.

When studying :math:`\alpha`-helices it is useful to have a *helical
wheel* projection of your peptide, to see whether a peptide is
amphipathic. This can be done using the :ref:`gmx wheel <gmx wheel>`
program. Two examples are plotted in
:numref:`Fig. %s <fig-hprwheel>`.

.. _fig-hprwheel:

.. figure:: plots/hpr-wheel.*
   :width: 5.00000cm

   Helical wheel projection of the N-terminal helix of HPr.

Interface-related items
-----------------------

| :ref:`gmx order <gmx order>`, :ref:`gmx density <gmx density>`, 
  :ref:`gmx potential <gmx potential>`, :ref:`gmx traj <gmx traj>`
| When simulating molecules with long carbon tails, it can be
  interesting to calculate their average orientation. There are several
  flavors of order parameters, most of which are related. The program
  :ref:`gmx order <gmx order>` can calculate
  order parameters using the equation:

.. math:: S_{z} = \frac{3}{2}\langle {\cos^2{\theta_z}} \rangle - \frac{1}{2}
          :label: eqnSgr

where :math:`\theta_z` is the angle between the :math:`z`-axis of the
simulation box and the molecular axis under consideration. The latter is
defined as the vector from C\ :math:`_{n-1}` to C\ :math:`_{n+1}`. The
parameters :math:`S_x` and :math:`S_y` are defined in the same way. The
brackets imply averaging over time and molecules. Order parameters can
vary between 1 (full order along the interface normal) and :math:`-1/2`
(full order perpendicular to the normal), with a value of zero in the
case of isotropic orientation.

The program can do two things for you. It can calculate the order
parameter for each CH\ :math:`_2` segment separately, for any of three
axes, or it can divide the box in slices and calculate the average value
of the order parameter per segment in one slice. The first method gives
an idea of the ordering of a molecule from head to tail, the second
method gives an idea of the ordering as function of the box length.

The electrostatic potential (:math:`\psi`) across the interface can be
computed from a trajectory by evaluating the double integral of the
charge density (:math:`\rho(z)`):

.. math:: \psi(z) - \psi(-\infty) = - \int_{-\infty}^z dz' \int_{-\infty}^{z'} \rho(z'')dz''/ \epsilon_0 
          :label: eqnelpotgr

where the position :math:`z=-\infty` is far enough in the bulk phase
such that the field is zero. With this method, it is possible to “split”
the total potential into separate contributions from lipid and water
molecules. The program :ref:`gmx potential <gmx potential>` divides the box in slices and sums
all charges of the atoms in each slice. It then integrates this charge
density to give the electric field, which is in turn integrated to give
the potential. Charge density, electric field, and potential are written
to xvgr input files.

The program :ref:`gmx traj <gmx traj>` is a very simple analysis program. All it does is
print the coordinates, velocities, or forces of selected atoms. It can
also calculate the center of mass of one or more molecules and print the
coordinates of the center of mass to three files. By itself, this is
probably not a very useful analysis, but having the coordinates of
selected molecules or atoms can be very handy for further analysis, not
only in interfacial systems.

The program :ref:`gmx density <gmx density>`
calculates the mass density of groups and gives a plot of the density
against a box axis. This is useful for looking at the distribution of
groups or atoms across the interface.

.. raw:: latex

    \clearpage


.. [1]
   :math:`P_0(x) = 1`, :math:`P_1(x) = x`, :math:`P_2(x) = (3x^2-1)/2`
