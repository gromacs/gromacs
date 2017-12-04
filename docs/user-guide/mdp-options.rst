.. README
   See the "run control" section for a working example of the
   syntax to use when making .mdp entries, with and without detailed
   documentation for values those entries might take. Everything can
   be cross-referenced, see the examples there. TODO Make more
   cross-references.

Molecular dynamics parameters (.mdp options)
============================================

.. _mdp-general:

General information
-------------------

Default values are given in parentheses, or listed first among
choices. The first option in the list is always the default
option. Units are given in square brackets. The difference between a
dash and an underscore is ignored.

A :ref:`sample mdp file <mdp>` is available. This should be
appropriate to start a normal simulation. Edit it to suit your
specific needs and desires.


Preprocessing
^^^^^^^^^^^^^

.. mdp:: include

   directories to include in your topology. Format:
   ``-I/home/john/mylib -I../otherlib``

.. mdp:: define

   defines to pass to the preprocessor, default is no defines. You can
   use any defines to control options in your customized topology
   files. Options that act on existing :ref:`top` file mechanisms
   include

      ``-DFLEXIBLE`` will use flexible water instead of rigid water
      into your topology, this can be useful for normal mode analysis.

      ``-DPOSRES`` will trigger the inclusion of ``posre.itp`` into
      your topology, used for implementing position restraints.


Run control
^^^^^^^^^^^

.. mdp:: integrator

   (Despite the name, this list includes algorithms that are not
   actually integrators over time. :mdp-value:`integrator=steep` and
   all entries following it are in this category)

   .. mdp-value:: md

      A leap-frog algorithm for integrating Newton's equations of motion.

   .. mdp-value:: md-vv

      A velocity Verlet algorithm for integrating Newton's equations
      of motion.  For constant NVE simulations started from
      corresponding points in the same trajectory, the trajectories
      are analytically, but not binary, identical to the
      :mdp-value:`integrator=md` leap-frog integrator. The the kinetic
      energy, which is determined from the whole step velocities and
      is therefore slightly too high. The advantage of this integrator
      is more accurate, reversible Nose-Hoover and Parrinello-Rahman
      coupling integration based on Trotter expansion, as well as
      (slightly too small) full step velocity output. This all comes
      at the cost off extra computation, especially with constraints
      and extra communication in parallel. Note that for nearly all
      production simulations the :mdp-value:`integrator=md` integrator
      is accurate enough.

   .. mdp-value:: md-vv-avek

      A velocity Verlet algorithm identical to
      :mdp-value:`integrator=md-vv`, except that the kinetic energy is
      determined as the average of the two half step kinetic energies
      as in the :mdp-value:`integrator=md` integrator, and this thus
      more accurate.  With Nose-Hoover and/or Parrinello-Rahman
      coupling this comes with a slight increase in computational
      cost.

   .. mdp-value:: sd

      An accurate and efficient leap-frog stochastic dynamics
      integrator. With constraints, coordinates needs to be
      constrained twice per integration step. Depending on the
      computational cost of the force calculation, this can take a
      significant part of the simulation time. The temperature for one
      or more groups of atoms (:mdp:`tc-grps`) is set with
      :mdp:`ref-t`, the inverse friction constant for each group is
      set with :mdp:`tau-t`.  The parameter :mdp:`tcoupl` is
      ignored. The random generator is initialized with
      :mdp:`ld-seed`. When used as a thermostat, an appropriate value
      for :mdp:`tau-t` is 2 ps, since this results in a friction that
      is lower than the internal friction of water, while it is high
      enough to remove excess heat NOTE: temperature deviations decay
      twice as fast as with a Berendsen thermostat with the same
      :mdp:`tau-t`.

   .. mdp-value:: bd

      An Euler integrator for Brownian or position Langevin dynamics,
      the velocity is the force divided by a friction coefficient
      (:mdp:`bd-fric`) plus random thermal noise (:mdp:`ref-t`). When
      :mdp:`bd-fric` is 0, the friction coefficient for each particle
      is calculated as mass/ :mdp:`tau-t`, as for the integrator
      :mdp-value:`integrator=sd`. The random generator is initialized
      with :mdp:`ld-seed`.

   .. mdp-value:: steep

      A steepest descent algorithm for energy minimization. The
      maximum step size is :mdp:`emstep`, the tolerance is
      :mdp:`emtol`.

   .. mdp-value:: cg

      A conjugate gradient algorithm for energy minimization, the
      tolerance is :mdp:`emtol`. CG is more efficient when a steepest
      descent step is done every once in a while, this is determined
      by :mdp:`nstcgsteep`. For a minimization prior to a normal mode
      analysis, which requires a very high accuracy, |Gromacs| should be
      compiled in double precision.

   .. mdp-value:: l-bfgs

      A quasi-Newtonian algorithm for energy minimization according to
      the low-memory Broyden-Fletcher-Goldfarb-Shanno approach. In
      practice this seems to converge faster than Conjugate Gradients,
      but due to the correction steps necessary it is not (yet)
      parallelized.

   .. mdp-value:: nm

      Normal mode analysis is performed on the structure in the :ref:`tpr`
      file.  |Gromacs| should be compiled in double precision.

   .. mdp-value:: tpi

      Test particle insertion. The last molecule in the topology is
      the test particle. A trajectory must be provided to ``mdrun
      -rerun``. This trajectory should not contain the molecule to be
      inserted. Insertions are performed :mdp:`nsteps` times in each
      frame at random locations and with random orientiations of the
      molecule. When :mdp:`nstlist` is larger than one,
      :mdp:`nstlist` insertions are performed in a sphere with radius
      :mdp:`rtpi` around a the same random location using the same
      neighborlist. Since neighborlist construction is expensive,
      one can perform several extra insertions with the same list
      almost for free. The random seed is set with
      :mdp:`ld-seed`. The temperature for the Boltzmann weighting is
      set with :mdp:`ref-t`, this should match the temperature of the
      simulation of the original trajectory. Dispersion correction is
      implemented correctly for TPI. All relevant quantities are
      written to the file specified with ``mdrun -tpi``. The
      distribution of insertion energies is written to the file
      specified with ``mdrun -tpid``. No trajectory or energy file is
      written. Parallel TPI gives identical results to single-node
      TPI. For charged molecules, using PME with a fine grid is most
      accurate and also efficient, since the potential in the system
      only needs to be calculated once per frame.

   .. mdp-value:: tpic

      Test particle insertion into a predefined cavity location. The
      procedure is the same as for :mdp-value:`integrator=tpi`, except
      that one coordinate extra is read from the trajectory, which is
      used as the insertion location. The molecule to be inserted
      should be centered at 0,0,0. |Gromacs| does not do this for you,
      since for different situations a different way of centering
      might be optimal. Also :mdp:`rtpi` sets the radius for the
      sphere around this location. Neighbor searching is done only
      once per frame, :mdp:`nstlist` is not used. Parallel
      :mdp-value:`integrator=tpic` gives identical results to
      single-rank :mdp-value:`integrator=tpic`.

.. mdp:: tinit

        (0) \[ps\]
        starting time for your run (only makes sense for time-based
        integrators)

.. mdp:: dt

        (0.001) \[ps\]
        time step for integration (only makes sense for time-based
        integrators)

.. mdp:: nsteps

        (0)
        maximum number of steps to integrate or minimize, -1 is no
        maximum

.. mdp:: init-step

        (0)
        The starting step. The time at an step i in a run is
        calculated as: t = :mdp:`tinit` + :mdp:`dt` *
        (:mdp:`init-step` + i). The free-energy lambda is calculated
        as: lambda = :mdp:`init-lambda` + :mdp:`delta-lambda` *
        (:mdp:`init-step` + i). Also non-equilibrium MD parameters can
        depend on the step number. Thus for exact restarts or redoing
        part of a run it might be necessary to set :mdp:`init-step` to
        the step number of the restart frame. :ref:`gmx convert-tpr`
        does this automatically.

.. mdp:: simulation-part

         (0)
         A simulation can consist of multiple parts, each of which has
         a part number. This option specifies what that number will
         be, which helps keep track of parts that are logically the
         same simulation. This option is generally useful to set only
         when coping with a crashed simulation where files were lost.

.. mdp:: comm-mode

   .. mdp-value:: Linear

      Remove center of mass translational velocity

   .. mdp-value:: Angular

      Remove center of mass translational and rotational velocity around
      the center of mass

   .. mdp-value:: Linear-acceleration-correction

      Remove center of mass translational velocity. Correct the center of
      mass position assuming linear acceleration over :mdp:`nstcomm` steps.
      This is useful for cases where an acceleration is expected on the
      center of mass which is nearly constant over mdp:`nstcomm` steps.
      This can occur for example when pulling on a group using an absolute
      reference.

   .. mdp-value:: None

      No restriction on the center of mass motion

.. mdp:: nstcomm

   (100) \[steps\]
   frequency for center of mass motion removal

.. mdp:: comm-grps

   group(s) for center of mass motion removal, default is the whole
   system


Langevin dynamics
^^^^^^^^^^^^^^^^^

.. mdp:: bd-fric

   (0) \[amu ps-1\]
   Brownian dynamics friction coefficient. When :mdp:`bd-fric` is 0,
   the friction coefficient for each particle is calculated as mass/
   :mdp:`tau-t`.

.. mdp:: ld-seed

   (-1) \[integer\]
   used to initialize random generator for thermal noise for
   stochastic and Brownian dynamics. When :mdp:`ld-seed` is set to -1,
   a pseudo random seed is used. When running BD or SD on multiple
   processors, each processor uses a seed equal to :mdp:`ld-seed` plus
   the processor number.


Energy minimization
^^^^^^^^^^^^^^^^^^^

.. mdp:: emtol

   (10.0) \[kJ mol-1 nm-1\]
   the minimization is converged when the maximum force is smaller
   than this value

.. mdp:: emstep

   (0.01) \[nm\]
   initial step-size

.. mdp:: nstcgsteep

   (1000) \[steps\]
   frequency of performing 1 steepest descent step while doing
   conjugate gradient energy minimization.

.. mdp:: nbfgscorr

   (10)
   Number of correction steps to use for L-BFGS minimization. A higher
   number is (at least theoretically) more accurate, but slower.


Shell Molecular Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^

When shells or flexible constraints are present in the system the
positions of the shells and the lengths of the flexible constraints
are optimized at every time step until either the RMS force on the
shells and constraints is less than :mdp:`emtol`, or a maximum number
of iterations :mdp:`niter` has been reached. Minimization is converged
when the maximum force is smaller than :mdp:`emtol`. For shell MD this
value should be 1.0 at most.

.. mdp:: niter

   (20)
   maximum number of iterations for optimizing the shell positions and
   the flexible constraints.

.. mdp:: fcstep

   (0) \[ps^2\]
   the step size for optimizing the flexible constraints. Should be
   chosen as mu/(d2V/dq2) where mu is the reduced mass of two
   particles in a flexible constraint and d2V/dq2 is the second
   derivative of the potential in the constraint direction. Hopefully
   this number does not differ too much between the flexible
   constraints, as the number of iterations and thus the runtime is
   very sensitive to fcstep. Try several values!


Test particle insertion
^^^^^^^^^^^^^^^^^^^^^^^

.. mdp:: rtpi

   (0.05) \[nm\]
   the test particle insertion radius, see integrators
   :mdp-value:`integrator=tpi` and :mdp-value:`integrator=tpic`


Output control
^^^^^^^^^^^^^^

.. mdp:: nstxout

   (0) \[steps\]
   number of steps that elapse between writing coordinates to output
   trajectory file, the last coordinates are always written

.. mdp:: nstvout

   (0) \[steps\]
   number of steps that elapse between writing velocities to output
   trajectory, the last velocities are always written

.. mdp:: nstfout

   (0) \[steps\]
   number of steps that elapse between writing forces to output
   trajectory.

.. mdp:: nstlog

   (1000) \[steps\]
   number of steps that elapse between writing energies to the log
   file, the last energies are always written

.. mdp:: nstcalcenergy

   (100)
   number of steps that elapse between calculating the energies, 0 is
   never. This option is only relevant with dynamics. This option affects the
   performance in parallel simulations, because calculating energies
   requires global communication between all processes which can
   become a bottleneck at high parallelization.

.. mdp:: nstenergy

   (1000) \[steps\]
   number of steps that else between writing energies to energy file,
   the last energies are always written, should be a multiple of
   :mdp:`nstcalcenergy`. Note that the exact sums and fluctuations
   over all MD steps modulo :mdp:`nstcalcenergy` are stored in the
   energy file, so :ref:`gmx energy` can report exact energy averages
   and fluctuations also when :mdp:`nstenergy` > 1

.. mdp:: nstxout-compressed

   (0) \[steps\]
   number of steps that elapse between writing position coordinates
   using lossy compression

.. mdp:: compressed-x-precision

   (1000) \[real\]
   precision with which to write to the compressed trajectory file

.. mdp:: compressed-x-grps

   group(s) to write to the compressed trajectory file, by default the
   whole system is written (if :mdp:`nstxout-compressed` > 0)

.. mdp:: energygrps

   group(s) for which to write to write short-ranged non-bonded
   potential energies to the energy file (not supported on GPUs)


Neighbor searching
^^^^^^^^^^^^^^^^^^

.. mdp:: cutoff-scheme

   .. mdp-value:: Verlet

      Generate a pair list with buffering. The buffer size is
      automatically set based on :mdp:`verlet-buffer-tolerance`,
      unless this is set to -1, in which case :mdp:`rlist` will be
      used. This option has an explicit, exact cut-off at :mdp:`rvdw`
      equal to :mdp:`rcoulomb`, unless PME or Ewald is used, in which
      case :mdp:`rcoulomb` > :mdp:`rvdw` is allowed. Currently only
      cut-off, reaction-field, PME or Ewald electrostatics and plain
      LJ are supported. Some :ref:`gmx mdrun` functionality is not yet
      supported with the :mdp-value:`cutoff-scheme=Verlet` scheme, but :ref:`gmx grompp`
      checks for this. Native GPU acceleration is only supported with
      :mdp-value:`cutoff-scheme=Verlet`. With GPU-accelerated PME or with separate PME
      ranks, :ref:`gmx mdrun` will automatically tune the CPU/GPU load
      balance by scaling :mdp:`rcoulomb` and the grid spacing. This
      can be turned off with ``mdrun -notunepme``. :mdp-value:`cutoff-scheme=Verlet` is
      faster than :mdp-value:`cutoff-scheme=group` when there is no water, or if
      :mdp-value:`cutoff-scheme=group` would use a pair-list buffer to conserve energy.

   .. mdp-value:: group

      Generate a pair list for groups of atoms. These groups
      correspond to the charge groups in the topology. This was the
      only cut-off treatment scheme before version 4.6, and is
      **deprecated in |gmx-version|**. There is no explicit buffering of
      the pair list. This enables efficient force calculations for
      water, but energy is only conserved when a buffer is explicitly
      added.

.. mdp:: nstlist

   \(10) \[steps\]

   .. mdp-value:: >0

      Frequency to update the neighbor list. When this is 0, the
      neighbor list is made only once. With energy minimization the
      neighborlist will be updated for every energy evaluation when
      :mdp:`nstlist` is greater than 0. With :mdp-value:`cutoff-scheme=Verlet` and
      :mdp:`verlet-buffer-tolerance` set, :mdp:`nstlist` is actually
      a minimum value and :ref:`gmx mdrun` might increase it, unless
      it is set to 1. With parallel simulations and/or non-bonded
      force calculation on the GPU, a value of 20 or 40 often gives
      the best performance. With :mdp-value:`cutoff-scheme=group` and non-exact
      cut-off's, :mdp:`nstlist` will affect the accuracy of your
      simulation and it can not be chosen freely.

   .. mdp-value:: 0

      The neighbor list is only constructed once and never
      updated. This is mainly useful for vacuum simulations in which
      all particles see each other.

   .. mdp-value:: <0

      Unused.

.. mdp:: ns-type

   .. mdp-value:: grid

      Make a grid in the box and only check atoms in neighboring grid
      cells when constructing a new neighbor list every
      :mdp:`nstlist` steps. In large systems grid search is much
      faster than simple search.

   .. mdp-value:: simple

      Check every atom in the box when constructing a new neighbor
      list every :mdp:`nstlist` steps (only with :mdp-value:`cutoff-scheme=group`
      cut-off scheme).

.. mdp:: pbc

   .. mdp-value:: xyz

      Use periodic boundary conditions in all directions.

   .. mdp-value:: no

      Use no periodic boundary conditions, ignore the box. To simulate
      without cut-offs, set all cut-offs and :mdp:`nstlist` to 0. For
      best performance without cut-offs on a single MPI rank, set
      :mdp:`nstlist` to zero and :mdp:`ns-type` =simple.

   .. mdp-value:: xy

      Use periodic boundary conditions in x and y directions
      only. This works only with :mdp:`ns-type` =grid and can be used
      in combination with walls_. Without walls or with only one wall
      the system size is infinite in the z direction. Therefore
      pressure coupling or Ewald summation methods can not be
      used. These disadvantages do not apply when two walls are used.

.. mdp:: periodic-molecules

   .. mdp-value:: no

      molecules are finite, fast molecular PBC can be used

   .. mdp-value:: yes

      for systems with molecules that couple to themselves through the
      periodic boundary conditions, this requires a slower PBC
      algorithm and molecules are not made whole in the output

.. mdp:: verlet-buffer-tolerance

   (0.005) \[kJ/mol/ps\]

   Useful only with the :mdp-value:`cutoff-scheme=Verlet` :mdp:`cutoff-scheme`. This sets
   the maximum allowed error for pair interactions per particle caused
   by the Verlet buffer, which indirectly sets :mdp:`rlist`. As both
   :mdp:`nstlist` and the Verlet buffer size are fixed (for
   performance reasons), particle pairs not in the pair list can
   occasionally get within the cut-off distance during
   :mdp:`nstlist` -1 steps. This causes very small jumps in the
   energy. In a constant-temperature ensemble, these very small energy
   jumps can be estimated for a given cut-off and :mdp:`rlist`. The
   estimate assumes a homogeneous particle distribution, hence the
   errors might be slightly underestimated for multi-phase
   systems. (See the `reference manual`_ for details). For longer
   pair-list life-time (:mdp:`nstlist` -1) * :mdp:`dt` the buffer is
   overestimated, because the interactions between particles are
   ignored. Combined with cancellation of errors, the actual drift of
   the total energy is usually one to two orders of magnitude
   smaller. Note that the generated buffer size takes into account
   that the |Gromacs| pair-list setup leads to a reduction in the
   drift by a factor 10, compared to a simple particle-pair based
   list. Without dynamics (energy minimization etc.), the buffer is 5%
   of the cut-off. For NVE simulations the initial temperature is
   used, unless this is zero, in which case a buffer of 10% is
   used. For NVE simulations the tolerance usually needs to be lowered
   to achieve proper energy conservation on the nanosecond time
   scale. To override the automated buffer setting, use
   :mdp:`verlet-buffer-tolerance` =-1 and set :mdp:`rlist` manually.

.. mdp:: rlist

   (1) \[nm\]
   Cut-off distance for the short-range neighbor list. With the
   :mdp-value:`cutoff-scheme=Verlet` :mdp:`cutoff-scheme`, this is by default set by the
   :mdp:`verlet-buffer-tolerance` option and the value of
   :mdp:`rlist` is ignored.


Electrostatics
^^^^^^^^^^^^^^

.. mdp:: coulombtype

   .. mdp-value:: Cut-off

      Plain cut-off with neighborlist radius :mdp:`rlist` and
      Coulomb cut-off :mdp:`rcoulomb`, where :mdp:`rlist` >=
      :mdp:`rcoulomb`.

   .. mdp-value:: Ewald

      Classical Ewald sum electrostatics. The real-space cut-off
      :mdp:`rcoulomb` should be equal to :mdp:`rlist`. Use *e.g.*
      :mdp:`rlist` =0.9, :mdp:`rcoulomb` =0.9. The highest magnitude
      of wave vectors used in reciprocal space is controlled by
      :mdp:`fourierspacing`. The relative accuracy of
      direct/reciprocal space is controlled by :mdp:`ewald-rtol`.

      NOTE: Ewald scales as O(N^3/2) and is thus extremely slow for
      large systems. It is included mainly for reference - in most
      cases PME will perform much better.

   .. mdp-value:: PME

      Fast smooth Particle-Mesh Ewald (SPME) electrostatics. Direct
      space is similar to the Ewald sum, while the reciprocal part is
      performed with FFTs. Grid dimensions are controlled with
      :mdp:`fourierspacing` and the interpolation order with
      :mdp:`pme-order`. With a grid spacing of 0.1 nm and cubic
      interpolation the electrostatic forces have an accuracy of
      2-3*10^-4. Since the error from the vdw-cutoff is larger than
      this you might try 0.15 nm. When running in parallel the
      interpolation parallelizes better than the FFT, so try
      decreasing grid dimensions while increasing interpolation.

   .. mdp-value:: P3M-AD

      Particle-Particle Particle-Mesh algorithm with analytical
      derivative for for long range electrostatic interactions. The
      method and code is identical to SPME, except that the influence
      function is optimized for the grid. This gives a slight increase
      in accuracy.

   .. mdp-value:: Reaction-Field

      Reaction field electrostatics with Coulomb cut-off
      :mdp:`rcoulomb`, where :mdp:`rlist` >= :mdp:`rvdw`. The
      dielectric constant beyond the cut-off is
      :mdp:`epsilon-rf`. The dielectric constant can be set to
      infinity by setting :mdp:`epsilon-rf` =0.

   .. mdp-value:: Generalized-Reaction-Field

      Generalized reaction field with Coulomb cut-off
      :mdp:`rcoulomb`, where :mdp:`rlist` >= :mdp:`rcoulomb`. The
      dielectric constant beyond the cut-off is
      :mdp:`epsilon-rf`. The ionic strength is computed from the
      number of charged (*i.e.* with non zero charge) charge
      groups. The temperature for the GRF potential is set with
      :mdp:`ref-t`.

   .. mdp-value:: Reaction-Field-zero

      In |Gromacs|, normal reaction-field electrostatics with
      :mdp:`cutoff-scheme` = :mdp-value:`cutoff-scheme=group` leads to bad energy
      conservation. :mdp-value:`coulombtype=Reaction-Field-zero` solves this by making
      the potential zero beyond the cut-off. It can only be used with
      an infinite dielectric constant (:mdp:`epsilon-rf` =0), because
      only for that value the force vanishes at the
      cut-off. :mdp:`rlist` should be 0.1 to 0.3 nm larger than
      :mdp:`rcoulomb` to accommodate for the size of charge groups
      and diffusion between neighbor list updates. This, and the fact
      that table lookups are used instead of analytical functions make
      :mdp-value:`coulombtype=Reaction-Field-zero` computationally more expensive than
      normal reaction-field.

   .. mdp-value:: Shift

      Analogous to :mdp-value:`vdwtype=Shift` for :mdp:`vdwtype`. You
      might want to use :mdp-value:`coulombtype=Reaction-Field-zero` instead, which has
      a similar potential shape, but has a physical interpretation and
      has better energies due to the exclusion correction terms.

   .. mdp-value:: Encad-Shift

      The Coulomb potential is decreased over the whole range, using
      the definition from the Encad simulation package.

   .. mdp-value:: Switch

      Analogous to :mdp-value:`vdwtype=Switch` for
      :mdp:`vdwtype`. Switching the Coulomb potential can lead to
      serious artifacts, advice: use :mdp-value:`coulombtype=Reaction-Field-zero`
      instead.

   .. mdp-value:: User

      :ref:`gmx mdrun` will now expect to find a file ``table.xvg``
      with user-defined potential functions for repulsion, dispersion
      and Coulomb. When pair interactions are present, :ref:`gmx
      mdrun` also expects to find a file ``tablep.xvg`` for the pair
      interactions. When the same interactions should be used for
      non-bonded and pair interactions the user can specify the same
      file name for both table files. These files should contain 7
      columns: the ``x`` value, ``f(x)``, ``-f'(x)``, ``g(x)``,
      ``-g'(x)``, ``h(x)``, ``-h'(x)``, where ``f(x)`` is the Coulomb
      function, ``g(x)`` the dispersion function and ``h(x)`` the
      repulsion function. When :mdp:`vdwtype` is not set to User the
      values for ``g``, ``-g'``, ``h`` and ``-h'`` are ignored. For
      the non-bonded interactions ``x`` values should run from 0 to
      the largest cut-off distance + :mdp:`table-extension` and
      should be uniformly spaced. For the pair interactions the table
      length in the file will be used. The optimal spacing, which is
      used for non-user tables, is ``0.002 nm`` when you run in mixed
      precision or ``0.0005 nm`` when you run in double precision. The
      function value at ``x=0`` is not important. More information is
      in the printed manual.

   .. mdp-value:: PME-Switch

      A combination of PME and a switch function for the direct-space
      part (see above). :mdp:`rcoulomb` is allowed to be smaller than
      :mdp:`rlist`. This is mainly useful constant energy simulations
      (note that using PME with :mdp:`cutoff-scheme` = :mdp-value:`cutoff-scheme=Verlet`
      will be more efficient).

   .. mdp-value:: PME-User

      A combination of PME and user tables (see
      above). :mdp:`rcoulomb` is allowed to be smaller than
      :mdp:`rlist`. The PME mesh contribution is subtracted from the
      user table by :ref:`gmx mdrun`. Because of this subtraction the
      user tables should contain about 10 decimal places.

   .. mdp-value:: PME-User-Switch

      A combination of PME-User and a switching function (see
      above). The switching function is applied to final
      particle-particle interaction, *i.e.* both to the user supplied
      function and the PME Mesh correction part.

.. mdp:: coulomb-modifier

   .. mdp-value:: Potential-shift-Verlet

      Selects Potential-shift with the Verlet cutoff-scheme, as it is
      (nearly) free; selects None with the group cutoff-scheme.

   .. mdp-value:: Potential-shift

      Shift the Coulomb potential by a constant such that it is zero
      at the cut-off. This makes the potential the integral of the
      force. Note that this does not affect the forces or the
      sampling.

   .. mdp-value:: None

      Use an unmodified Coulomb potential. With the group scheme this
      means no exact cut-off is used, energies and forces are
      calculated for all pairs in the neighborlist.

.. mdp:: rcoulomb-switch

   (0) \[nm\]
   where to start switching the Coulomb potential, only relevant
   when force or potential switching is used

.. mdp:: rcoulomb

   (1) \[nm\]
   distance for the Coulomb cut-off

.. mdp:: epsilon-r

   (1)
   The relative dielectric constant. A value of 0 means infinity.

.. mdp:: epsilon-rf

   (0)
   The relative dielectric constant of the reaction field. This
   is only used with reaction-field electrostatics. A value of 0
   means infinity.


Van der Waals
^^^^^^^^^^^^^

.. mdp:: vdwtype

   .. mdp-value:: Cut-off

      Twin range cut-offs with neighbor list cut-off :mdp:`rlist` and
      VdW cut-off :mdp:`rvdw`, where :mdp:`rvdw` >= :mdp:`rlist`.

   .. mdp-value:: PME

      Fast smooth Particle-mesh Ewald (SPME) for VdW interactions. The
      grid dimensions are controlled with :mdp:`fourierspacing` in
      the same way as for electrostatics, and the interpolation order
      is controlled with :mdp:`pme-order`. The relative accuracy of
      direct/reciprocal space is controlled by :mdp:`ewald-rtol-lj`,
      and the specific combination rules that are to be used by the
      reciprocal routine are set using :mdp:`lj-pme-comb-rule`.

   .. mdp-value:: Shift

      This functionality is deprecated and replaced by
      :mdp:`vdw-modifier` = Force-switch. The LJ (not Buckingham)
      potential is decreased over the whole range and the forces decay
      smoothly to zero between :mdp:`rvdw-switch` and
      :mdp:`rvdw`. The neighbor search cut-off :mdp:`rlist` should
      be 0.1 to 0.3 nm larger than :mdp:`rvdw` to accommodate for the
      size of charge groups and diffusion between neighbor list
      updates.

   .. mdp-value:: Switch

      This functionality is deprecated and replaced by
      :mdp:`vdw-modifier` = Potential-switch. The LJ (not Buckingham)
      potential is normal out to :mdp:`rvdw-switch`, after which it
      is switched off to reach zero at :mdp:`rvdw`. Both the
      potential and force functions are continuously smooth, but be
      aware that all switch functions will give rise to a bulge
      (increase) in the force (since we are switching the
      potential). The neighbor search cut-off :mdp:`rlist` should be
      0.1 to 0.3 nm larger than :mdp:`rvdw` to accommodate for the
      size of charge groups and diffusion between neighbor list
      updates.

   .. mdp-value:: Encad-Shift

      The LJ (not Buckingham) potential is decreased over the whole
      range, using the definition from the Encad simulation package.

   .. mdp-value:: User

      See user for :mdp:`coulombtype`. The function value at zero is
      not important. When you want to use LJ correction, make sure
      that :mdp:`rvdw` corresponds to the cut-off in the user-defined
      function. When :mdp:`coulombtype` is not set to User the values
      for the ``f`` and ``-f'`` columns are ignored.

.. mdp:: vdw-modifier

   .. mdp-value:: Potential-shift-Verlet

      Selects Potential-shift with the Verlet cutoff-scheme, as it is
      (nearly) free; selects None with the group cutoff-scheme.

   .. mdp-value:: Potential-shift

      Shift the Van der Waals potential by a constant such that it is
      zero at the cut-off. This makes the potential the integral of
      the force. Note that this does not affect the forces or the
      sampling.

   .. mdp-value:: None

      Use an unmodified Van der Waals potential. With the group scheme
      this means no exact cut-off is used, energies and forces are
      calculated for all pairs in the neighborlist.

   .. mdp-value:: Force-switch

      Smoothly switches the forces to zero between :mdp:`rvdw-switch`
      and :mdp:`rvdw`. This shifts the potential shift over the whole
      range and switches it to zero at the cut-off. Note that this is
      more expensive to calculate than a plain cut-off and it is not
      required for energy conservation, since Potential-shift
      conserves energy just as well.

   .. mdp-value:: Potential-switch

      Smoothly switches the potential to zero between
      :mdp:`rvdw-switch` and :mdp:`rvdw`. Note that this introduces
      articifically large forces in the switching region and is much
      more expensive to calculate. This option should only be used if
      the force field you are using requires this.

.. mdp:: rvdw-switch

   (0) \[nm\]

   where to start switching the LJ force and possibly the potential,
   only relevant when force or potential switching is used

.. mdp:: rvdw

   (1) \[nm\]
   distance for the LJ or Buckingham cut-off

.. mdp:: DispCorr

   .. mdp-value:: no

      don't apply any correction

   .. mdp-value:: EnerPres

      apply long range dispersion corrections for Energy and Pressure

   .. mdp-value:: Ener

      apply long range dispersion corrections for Energy only


Tables
^^^^^^

.. mdp:: table-extension

   (1) \[nm\]
   Extension of the non-bonded potential lookup tables beyond the
   largest cut-off distance. The value should be large enough to
   account for charge group sizes and the diffusion between
   neighbor-list updates. Without user defined potential the same
   table length is used for the lookup tables for the 1-4
   interactions, which are always tabulated irrespective of the use of
   tables for the non-bonded interactions. The value of
   :mdp:`table-extension` in no way affects the values of
   :mdp:`rlist`, :mdp:`rcoulomb`, or :mdp:`rvdw`.

.. mdp:: energygrp-table

   When user tables are used for electrostatics and/or VdW, here one
   can give pairs of energy groups for which seperate user tables
   should be used. The two energy groups will be appended to the table
   file name, in order of their definition in :mdp:`energygrps`,
   seperated by underscores. For example, if ``energygrps = Na Cl
   Sol`` and ``energygrp-table = Na Na Na Cl``, :ref:`gmx mdrun` will
   read ``table_Na_Na.xvg`` and ``table_Na_Cl.xvg`` in addition to the
   normal ``table.xvg`` which will be used for all other energy group
   pairs.


Ewald
^^^^^

.. mdp:: fourierspacing

   (0.12) \[nm\]
   For ordinary Ewald, the ratio of the box dimensions and the spacing
   determines a lower bound for the number of wave vectors to use in
   each (signed) direction. For PME and P3M, that ratio determines a
   lower bound for the number of Fourier-space grid points that will
   be used along that axis. In all cases, the number for each
   direction can be overridden by entering a non-zero value for that
   :mdp:`fourier-nx` direction. For optimizing the relative load of
   the particle-particle interactions and the mesh part of PME, it is
   useful to know that the accuracy of the electrostatics remains
   nearly constant when the Coulomb cut-off and the PME grid spacing
   are scaled by the same factor.

.. mdp:: fourier-nx
.. mdp:: fourier-ny
.. mdp:: fourier-nz

   (0)
   Highest magnitude of wave vectors in reciprocal space when using Ewald.
   Grid size when using PME or P3M. These values override
   :mdp:`fourierspacing` per direction. The best choice is powers of
   2, 3, 5 and 7. Avoid large primes.

.. mdp:: pme-order

   (4)
   Interpolation order for PME. 4 equals cubic interpolation. You
   might try 6/8/10 when running in parallel and simultaneously
   decrease grid dimension.

.. mdp:: ewald-rtol

   (1e-5)
   The relative strength of the Ewald-shifted direct potential at
   :mdp:`rcoulomb` is given by :mdp:`ewald-rtol`. Decreasing this
   will give a more accurate direct sum, but then you need more wave
   vectors for the reciprocal sum.

.. mdp:: ewald-rtol-lj

   (1e-3)
   When doing PME for VdW-interactions, :mdp:`ewald-rtol-lj` is used
   to control the relative strength of the dispersion potential at
   :mdp:`rvdw` in the same way as :mdp:`ewald-rtol` controls the
   electrostatic potential.

.. mdp:: lj-pme-comb-rule

   (Geometric)
   The combination rules used to combine VdW-parameters in the
   reciprocal part of LJ-PME. Geometric rules are much faster than
   Lorentz-Berthelot and usually the recommended choice, even when the
   rest of the force field uses the Lorentz-Berthelot rules.

   .. mdp-value:: Geometric

      Apply geometric combination rules

   .. mdp-value:: Lorentz-Berthelot

      Apply Lorentz-Berthelot combination rules

.. mdp:: ewald-geometry

   .. mdp-value:: 3d

      The Ewald sum is performed in all three dimensions.

   .. mdp-value:: 3dc

      The reciprocal sum is still performed in 3D, but a force and
      potential correction applied in the `z` dimension to produce a
      pseudo-2D summation. If your system has a slab geometry in the
      `x-y` plane you can try to increase the `z`-dimension of the box
      (a box height of 3 times the slab height is usually ok) and use
      this option.

.. mdp:: epsilon-surface

   (0)
   This controls the dipole correction to the Ewald summation in
   3D. The default value of zero means it is turned off. Turn it on by
   setting it to the value of the relative permittivity of the
   imaginary surface around your infinite system. Be careful - you
   shouldn't use this if you have free mobile charges in your
   system. This value does not affect the slab 3DC variant of the long
   range corrections.


Temperature coupling
^^^^^^^^^^^^^^^^^^^^

.. mdp:: tcoupl

   .. mdp-value:: no

      No temperature coupling.

   .. mdp-value:: berendsen

      Temperature coupling with a Berendsen-thermostat to a bath with
      temperature :mdp:`ref-t`, with time constant
      :mdp:`tau-t`. Several groups can be coupled separately, these
      are specified in the :mdp:`tc-grps` field separated by spaces.

   .. mdp-value:: nose-hoover

      Temperature coupling using a Nose-Hoover extended ensemble. The
      reference temperature and coupling groups are selected as above,
      but in this case :mdp:`tau-t` controls the period of the
      temperature fluctuations at equilibrium, which is slightly
      different from a relaxation time. For NVT simulations the
      conserved energy quantity is written to energy and log file.

   .. mdp-value:: andersen

      Temperature coupling by randomizing a fraction of the particles
      at each timestep. Reference temperature and coupling groups are
      selected as above. :mdp:`tau-t` is the average time between
      randomization of each molecule. Inhibits particle dynamics
      somewhat, but little or no ergodicity issues. Currently only
      implemented with velocity Verlet, and not implemented with
      constraints.

   .. mdp-value:: andersen-massive

      Temperature coupling by randomizing all particles at infrequent
      timesteps. Reference temperature and coupling groups are
      selected as above. :mdp:`tau-t` is the time between
      randomization of all molecules. Inhibits particle dynamics
      somewhat, but little or no ergodicity issues. Currently only
      implemented with velocity Verlet.

   .. mdp-value:: v-rescale

      Temperature coupling using velocity rescaling with a stochastic
      term (JCP 126, 014101). This thermostat is similar to Berendsen
      coupling, with the same scaling using :mdp:`tau-t`, but the
      stochastic term ensures that a proper canonical ensemble is
      generated. The random seed is set with :mdp:`ld-seed`. This
      thermostat works correctly even for :mdp:`tau-t` =0. For NVT
      simulations the conserved energy quantity is written to the
      energy and log file.

.. mdp:: nsttcouple

   (-1)
   The frequency for coupling the temperature. The default value of -1
   sets :mdp:`nsttcouple` equal to :mdp:`nstlist`, unless
   :mdp:`nstlist` <=0, then a value of 10 is used. For velocity
   Verlet integrators :mdp:`nsttcouple` is set to 1.

.. mdp:: nh-chain-length

   (10)
   The number of chained Nose-Hoover thermostats for velocity Verlet
   integrators, the leap-frog :mdp-value:`integrator=md` integrator
   only supports 1. Data for the NH chain variables is not printed
   to the :ref:`edr` file by default, but can be turned on with the
   :mdp:`print-nose-hoover-chains` option.

.. mdp:: print-nose-hoover-chain-variables

   .. mdp-value:: no

      Do not store Nose-Hoover chain variables in the energy file.

   .. mdp-value:: yes

      Store all positions and velocities of the Nose-Hoover chain
      in the energy file.

.. mdp:: tc-grps

   groups to couple to separate temperature baths

.. mdp:: tau-t

   \[ps\]
   time constant for coupling (one for each group in
   :mdp:`tc-grps`), -1 means no temperature coupling

.. mdp:: ref-t

   \[K\]
   reference temperature for coupling (one for each group in
   :mdp:`tc-grps`)


Pressure coupling
^^^^^^^^^^^^^^^^^

.. mdp:: pcoupl

   .. mdp-value:: no

      No pressure coupling. This means a fixed box size.

   .. mdp-value:: Berendsen

      Exponential relaxation pressure coupling with time constant
      :mdp:`tau-p`. The box is scaled every timestep. It has been
      argued that this does not yield a correct thermodynamic
      ensemble, but it is the most efficient way to scale a box at the
      beginning of a run.

   .. mdp-value:: Parrinello-Rahman

      Extended-ensemble pressure coupling where the box vectors are
      subject to an equation of motion. The equation of motion for the
      atoms is coupled to this. No instantaneous scaling takes
      place. As for Nose-Hoover temperature coupling the time constant
      :mdp:`tau-p` is the period of pressure fluctuations at
      equilibrium. This is probably a better method when you want to
      apply pressure scaling during data collection, but beware that
      you can get very large oscillations if you are starting from a
      different pressure. For simulations where the exact fluctation
      of the NPT ensemble are important, or if the pressure coupling
      time is very short it may not be appropriate, as the previous
      time step pressure is used in some steps of the |Gromacs|
      implementation for the current time step pressure.

   .. mdp-value:: MTTK

      Martyna-Tuckerman-Tobias-Klein implementation, only useable with
      :mdp-value:`integrator=md-vv` or :mdp-value:`integrator=md-vv-avek`, very similar to
      Parrinello-Rahman. As for Nose-Hoover temperature coupling the
      time constant :mdp:`tau-p` is the period of pressure
      fluctuations at equilibrium. This is probably a better method
      when you want to apply pressure scaling during data collection,
      but beware that you can get very large oscillations if you are
      starting from a different pressure. Currently (as of version
      5.1), it only supports isotropic scaling, and only works without
      constraints.

.. mdp:: pcoupltype

   Specifies the kind of isotropy of the pressure coupling used. Each
   kind takes one or more values for :mdp:`compressibility` and
   :mdp:`ref-p`. Only a single value is permitted for :mdp:`tau-p`.

   .. mdp-value:: isotropic

      Isotropic pressure coupling with time constant
      :mdp:`tau-p`. One value each for :mdp:`compressibility` and
      :mdp:`ref-p` is required.

   .. mdp-value:: semiisotropic

      Pressure coupling which is isotropic in the ``x`` and ``y``
      direction, but different in the ``z`` direction. This can be
      useful for membrane simulations. Two values each for
      :mdp:`compressibility` and :mdp:`ref-p` are required, for
      ``x/y`` and ``z`` directions respectively.

   .. mdp-value:: anisotropic

      Same as before, but 6 values are needed for ``xx``, ``yy``, ``zz``,
      ``xy/yx``, ``xz/zx`` and ``yz/zy`` components,
      respectively. When the off-diagonal compressibilities are set to
      zero, a rectangular box will stay rectangular. Beware that
      anisotropic scaling can lead to extreme deformation of the
      simulation box.

   .. mdp-value:: surface-tension

      Surface tension coupling for surfaces parallel to the
      xy-plane. Uses normal pressure coupling for the `z`-direction,
      while the surface tension is coupled to the `x/y` dimensions of
      the box. The first :mdp:`ref-p` value is the reference surface
      tension times the number of surfaces ``bar nm``, the second
      value is the reference `z`-pressure ``bar``. The two
      :mdp:`compressibility` values are the compressibility in the
      `x/y` and `z` direction respectively. The value for the
      `z`-compressibility should be reasonably accurate since it
      influences the convergence of the surface-tension, it can also
      be set to zero to have a box with constant height.

.. mdp:: nstpcouple

   (-1)
   The frequency for coupling the pressure. The default value of -1
   sets :mdp:`nstpcouple` equal to :mdp:`nstlist`, unless
   :mdp:`nstlist` <=0, then a value of 10 is used. For velocity
   Verlet integrators :mdp:`nstpcouple` is set to 1.

.. mdp:: tau-p

   (1) \[ps\]
   The time constant for pressure coupling (one value for all
   directions).

.. mdp:: compressibility

   \[bar^-1\]
   The compressibility (NOTE: this is now really in bar^-1) For water at 1
   atm and 300 K the compressibility is 4.5e-5 bar^-1. The number of
   required values is implied by :mdp:`pcoupltype`.

.. mdp:: ref-p

   \[bar\]
   The reference pressure for coupling. The number of required values
   is implied by :mdp:`pcoupltype`.

.. mdp:: refcoord-scaling

   .. mdp-value:: no

      The reference coordinates for position restraints are not
      modified. Note that with this option the virial and pressure
      will depend on the absolute positions of the reference
      coordinates.

   .. mdp-value:: all

      The reference coordinates are scaled with the scaling matrix of
      the pressure coupling.

   .. mdp-value:: com

      Scale the center of mass of the reference coordinates with the
      scaling matrix of the pressure coupling. The vectors of each
      reference coordinate to the center of mass are not scaled. Only
      one COM is used, even when there are multiple molecules with
      position restraints. For calculating the COM of the reference
      coordinates in the starting configuration, periodic boundary
      conditions are not taken into account.


Simulated annealing
^^^^^^^^^^^^^^^^^^^

Simulated annealing is controlled separately for each temperature
group in |Gromacs|. The reference temperature is a piecewise linear
function, but you can use an arbitrary number of points for each
group, and choose either a single sequence or a periodic behaviour for
each group. The actual annealing is performed by dynamically changing
the reference temperature used in the thermostat algorithm selected,
so remember that the system will usually not instantaneously reach the
reference temperature!

.. mdp:: annealing

   Type of annealing for each temperature group

   .. mdp-value:: no

       No simulated annealing - just couple to reference temperature value.

   .. mdp-value:: single

       A single sequence of annealing points. If your simulation is
       longer than the time of the last point, the temperature will be
       coupled to this constant value after the annealing sequence has
       reached the last time point.

   .. mdp-value:: periodic

       The annealing will start over at the first reference point once
       the last reference time is reached. This is repeated until the
       simulation ends.

.. mdp:: annealing-npoints

   A list with the number of annealing reference/control points used
   for each temperature group. Use 0 for groups that are not
   annealed. The number of entries should equal the number of
   temperature groups.

.. mdp:: annealing-time

   List of times at the annealing reference/control points for each
   group. If you are using periodic annealing, the times will be used
   modulo the last value, *i.e.* if the values are 0, 5, 10, and 15,
   the coupling will restart at the 0ps value after 15ps, 30ps, 45ps,
   etc. The number of entries should equal the sum of the numbers
   given in :mdp:`annealing-npoints`.

.. mdp:: annealing-temp

   List of temperatures at the annealing reference/control points for
   each group. The number of entries should equal the sum of the
   numbers given in :mdp:`annealing-npoints`.

Confused? OK, let's use an example. Assume you have two temperature
groups, set the group selections to ``annealing = single periodic``,
the number of points of each group to ``annealing-npoints = 3 4``, the
times to ``annealing-time = 0 3 6 0 2 4 6`` and finally temperatures
to ``annealing-temp = 298 280 270 298 320 320 298``. The first group
will be coupled to 298K at 0ps, but the reference temperature will
drop linearly to reach 280K at 3ps, and then linearly between 280K and
270K from 3ps to 6ps. After this is stays constant, at 270K. The
second group is coupled to 298K at 0ps, it increases linearly to 320K
at 2ps, where it stays constant until 4ps. Between 4ps and 6ps it
decreases to 298K, and then it starts over with the same pattern
again, *i.e.* rising linearly from 298K to 320K between 6ps and
8ps. Check the summary printed by :ref:`gmx grompp` if you are unsure!


Velocity generation
^^^^^^^^^^^^^^^^^^^

.. mdp:: gen-vel

   .. mdp-value:: no

        Do not generate velocities. The velocities are set to zero
        when there are no velocities in the input structure file.

   .. mdp-value:: yes

        Generate velocities in :ref:`gmx grompp` according to a
        Maxwell distribution at temperature :mdp:`gen-temp`, with
        random seed :mdp:`gen-seed`. This is only meaningful with
        integrator :mdp-value:`integrator=md`.

.. mdp:: gen-temp

   (300) \[K\]
   temperature for Maxwell distribution

.. mdp:: gen-seed

   (-1) \[integer\]
   used to initialize random generator for random velocities,
   when :mdp:`gen-seed` is set to -1, a pseudo random seed is
   used.


Bonds
^^^^^

.. mdp:: constraints

   .. mdp-value:: none

      No constraints except for those defined explicitly in the
      topology, *i.e.* bonds are represented by a harmonic (or other)
      potential or a Morse potential (depending on the setting of
      :mdp:`morse`) and angles by a harmonic (or other) potential.

   .. mdp-value:: h-bonds

      Convert the bonds with H-atoms to constraints.

   .. mdp-value:: all-bonds

      Convert all bonds to constraints.

   .. mdp-value:: h-angles

      Convert all bonds and additionally the angles that involve
      H-atoms to bond-constraints.

   .. mdp-value:: all-angles

      Convert all bonds and angles to bond-constraints.

.. mdp:: constraint-algorithm

   .. mdp-value:: LINCS

      LINear Constraint Solver. With domain decomposition the parallel
      version P-LINCS is used. The accuracy in set with
      :mdp:`lincs-order`, which sets the number of matrices in the
      expansion for the matrix inversion. After the matrix inversion
      correction the algorithm does an iterative correction to
      compensate for lengthening due to rotation. The number of such
      iterations can be controlled with :mdp:`lincs-iter`. The root
      mean square relative constraint deviation is printed to the log
      file every :mdp:`nstlog` steps. If a bond rotates more than
      :mdp:`lincs-warnangle` in one step, a warning will be printed
      both to the log file and to ``stderr``. LINCS should not be used
      with coupled angle constraints.

   .. mdp-value:: SHAKE

      SHAKE is slightly slower and less stable than LINCS, but does
      work with angle constraints. The relative tolerance is set with
      :mdp:`shake-tol`, 0.0001 is a good value for "normal" MD. SHAKE
      does not support constraints between atoms on different nodes,
      thus it can not be used with domain decompositon when inter
      charge-group constraints are present. SHAKE can not be used with
      energy minimization.

.. mdp:: continuation

   This option was formerly known as unconstrained-start.

   .. mdp-value:: no

      apply constraints to the start configuration and reset shells

   .. mdp-value:: yes

      do not apply constraints to the start configuration and do not
      reset shells, useful for exact coninuation and reruns

.. mdp:: shake-tol

   (0.0001)
   relative tolerance for SHAKE

.. mdp:: lincs-order

   (4)
   Highest order in the expansion of the constraint coupling
   matrix. When constraints form triangles, an additional expansion of
   the same order is applied on top of the normal expansion only for
   the couplings within such triangles. For "normal" MD simulations an
   order of 4 usually suffices, 6 is needed for large time-steps with
   virtual sites or BD. For accurate energy minimization an order of 8
   or more might be required. With domain decomposition, the cell size
   is limited by the distance spanned by :mdp:`lincs-order` +1
   constraints. When one wants to scale further than this limit, one
   can decrease :mdp:`lincs-order` and increase :mdp:`lincs-iter`,
   since the accuracy does not deteriorate when (1+ :mdp:`lincs-iter`
   )* :mdp:`lincs-order` remains constant.

.. mdp:: lincs-iter

   (1)
   Number of iterations to correct for rotational lengthening in
   LINCS. For normal runs a single step is sufficient, but for NVE
   runs where you want to conserve energy accurately or for accurate
   energy minimization you might want to increase it to 2.

.. mdp:: lincs-warnangle

   (30) \[deg\]
   maximum angle that a bond can rotate before LINCS will complain

.. mdp:: morse

   .. mdp-value:: no

      bonds are represented by a harmonic potential

   .. mdp-value:: yes

      bonds are represented by a Morse potential


Energy group exclusions
^^^^^^^^^^^^^^^^^^^^^^^

.. mdp:: energygrp-excl

   Pairs of energy groups for which all non-bonded interactions are
   excluded. An example: if you have two energy groups ``Protein`` and
   ``SOL``, specifying ``energygrp-excl = Protein Protein SOL SOL``
   would give only the non-bonded interactions between the protein and
   the solvent. This is especially useful for speeding up energy
   calculations with ``mdrun -rerun`` and for excluding interactions
   within frozen groups.


Walls
^^^^^

.. mdp:: nwall

   (0)
   When set to 1 there is a wall at ``z=0``, when set to 2 there is
   also a wall at ``z=z-box``. Walls can only be used with :mdp:`pbc`
   ``=xy``. When set to 2 pressure coupling and Ewald summation can be
   used (it is usually best to use semiisotropic pressure coupling
   with the ``x/y`` compressibility set to 0, as otherwise the surface
   area will change). Walls interact wit the rest of the system
   through an optional :mdp:`wall-atomtype`. Energy groups ``wall0``
   and ``wall1`` (for :mdp:`nwall` =2) are added automatically to
   monitor the interaction of energy groups with each wall. The center
   of mass motion removal will be turned off in the ``z``-direction.

.. mdp:: wall-atomtype

   the atom type name in the force field for each wall. By (for
   example) defining a special wall atom type in the topology with its
   own combination rules, this allows for independent tuning of the
   interaction of each atomtype with the walls.

.. mdp:: wall-type

   .. mdp-value:: 9-3

      LJ integrated over the volume behind the wall: 9-3 potential

   .. mdp-value:: 10-4

      LJ integrated over the wall surface: 10-4 potential

   .. mdp-value:: 12-6

      direct LJ potential with the ``z`` distance from the wall

.. mdp:: table

   user defined potentials indexed with the ``z`` distance from the
   wall, the tables are read analogously to the
   :mdp:`energygrp-table` option, where the first name is for a
   "normal" energy group and the second name is ``wall0`` or
   ``wall1``, only the dispersion and repulsion columns are used

.. mdp:: wall-r-linpot

   (-1) \[nm\]
   Below this distance from the wall the potential is continued
   linearly and thus the force is constant. Setting this option to a
   postive value is especially useful for equilibration when some
   atoms are beyond a wall. When the value is <=0 (<0 for
   :mdp:`wall-type` =table), a fatal error is generated when atoms
   are beyond a wall.

.. mdp:: wall-density

   \[nm^-3/nm^-2\]
   the number density of the atoms for each wall for wall types 9-3
   and 10-4

.. mdp:: wall-ewald-zfac

   (3)
   The scaling factor for the third box vector for Ewald summation
   only, the minimum is 2. Ewald summation can only be used with
   :mdp:`nwall` =2, where one should use :mdp:`ewald-geometry`
   ``=3dc``. The empty layer in the box serves to decrease the
   unphysical Coulomb interaction between periodic images.


COM pulling
^^^^^^^^^^^

Note that where pulling coordinate are applicable, there can be more
than one (set with :mdp:`pull-ncoords`) and multiple related :ref:`mdp`
variables will exist accordingly. Documentation references to things
like :mdp:`pull-coord1-vec` should be understood to apply to to the
applicable pulling coordinate.

.. mdp:: pull

   .. mdp-value:: no

      No center of mass pulling. All the following pull options will
      be ignored (and if present in the :ref:`mdp` file, they unfortunately
      generate warnings)

   .. mdp-value:: yes

       Center of mass pulling will be applied on 1 or more groups using
       1 or more pull coordinates.

.. mdp:: pull-cylinder-r

   (1.5) \[nm\]
   the radius of the cylinder for
   :mdp:`pull-coord1-geometry` = :mdp-value:`pull-coord1-geometry=cylinder`

.. mdp:: pull-constr-tol

   (1e-6)
   the relative constraint tolerance for constraint pulling

.. mdp:: pull-print-com

   .. mdp-value:: no

      do not print the COM for any group

   .. mdp-value:: yes

      print the COM of all groups for all pull coordinates

.. mdp:: pull-print-ref-value

   .. mdp-value:: no

      do not print the reference value for each pull coordinate

   .. mdp-value:: yes

      print the reference value for each pull coordinate

.. mdp:: pull-print-components

   .. mdp-value:: no

      only print the distance for each pull coordinate
   
   .. mdp-value:: yes

      print the distance and Cartesian components selected in
      :mdp:`pull-coord1-dim`

.. mdp:: pull-nstxout

   (50)
   frequency for writing out the COMs of all the pull group (0 is
   never)

.. mdp:: pull-nstfout

   (50)
   frequency for writing out the force of all the pulled group
   (0 is never)


.. mdp:: pull-ngroups

   (1)
   The number of pull groups, not including the absolute reference
   group, when used. Pull groups can be reused in multiple pull
   coordinates. Below only the pull options for group 1 are given,
   further groups simply increase the group index number.

.. mdp:: pull-ncoords

   (1)
   The number of pull coordinates. Below only the pull options for
   coordinate 1 are given, further coordinates simply increase the
   coordinate index number.

.. mdp:: pull-group1-name

   The name of the pull group, is looked up in the index file or in
   the default groups to obtain the atoms involved.

.. mdp:: pull-group1-weights

   Optional relative weights which are multiplied with the masses of
   the atoms to give the total weight for the COM. The number should
   be 0, meaning all 1, or the number of atoms in the pull group.

.. mdp:: pull-group1-pbcatom

   (0)
   The reference atom for the treatment of periodic boundary
   conditions inside the group (this has no effect on the treatment of
   the pbc between groups). This option is only important when the
   diameter of the pull group is larger than half the shortest box
   vector. For determining the COM, all atoms in the group are put at
   their periodic image which is closest to
   :mdp:`pull-group1-pbcatom`. A value of 0 means that the middle
   atom (number wise) is used. This parameter is not used with
   :mdp:`pull-coord1-geometry` cylinder. A value of -1 turns on cosine
   weighting, which is useful for a group of molecules in a periodic
   system, *e.g.* a water slab (see Engin et al. J. Chem. Phys. B
   2010).

.. mdp:: pull-coord1-type

   .. mdp-value:: umbrella

      Center of mass pulling using an umbrella potential between the
      reference group and one or more groups.

   .. mdp-value:: constraint

      Center of mass pulling using a constraint between the reference
      group and one or more groups. The setup is identical to the
      option umbrella, except for the fact that a rigid constraint is
      applied instead of a harmonic potential.

   .. mdp-value:: constant-force

      Center of mass pulling using a linear potential and therefore a
      constant force. For this option there is no reference position
      and therefore the parameters :mdp:`pull-coord1-init` and
      :mdp:`pull-coord1-rate` are not used.

   .. mdp-value:: flat-bottom

      At distances above :mdp:`pull-coord1-init` a harmonic potential
      is applied, otherwise no potential is applied.

   .. mdp-value:: flat-bottom-high

      At distances below :mdp:`pull-coord1-init` a harmonic potential
      is applied, otherwise no potential is applied.

   .. mdp-value:: external-potential

      An external potential that needs to be provided by another
      module.

.. mdp:: pull-coord1-potential-provider

      The name of the external module that provides the potential for
      the case where :mdp:`pull-coord1-type` is external-potential.

.. mdp:: pull-coord1-geometry

   .. mdp-value:: distance

      Pull along the vector connecting the two groups. Components can
      be selected with :mdp:`pull-coord1-dim`.

   .. mdp-value:: direction

      Pull in the direction of :mdp:`pull-coord1-vec`.

   .. mdp-value:: direction-periodic

      As :mdp-value:`pull-coord1-geometry=direction`, but allows the distance to be larger
      than half the box size. With this geometry the box should not be
      dynamic (*e.g.* no pressure scaling) in the pull dimensions and
      the pull force is not added to virial.

   .. mdp-value:: direction-relative

      As :mdp-value:`pull-coord1-geometry=direction`, but the pull vector is the vector
      that points from the COM of a third to the COM of a fourth pull
      group. This means that 4 groups need to be supplied in
      :mdp:`pull-coord1-groups`. Note that the pull force will give
      rise to a torque on the pull vector, which is turn leads to
      forces perpendicular to the pull vector on the two groups
      defining the vector. If you want a pull group to move between
      the two groups defining the vector, simply use the union of
      these two groups as the reference group.

   .. mdp-value:: cylinder

      Designed for pulling with respect to a layer where the reference
      COM is given by a local cylindrical part of the reference group.
      The pulling is in the direction of :mdp:`pull-coord1-vec`. From
      the first of the two groups in :mdp:`pull-coord1-groups` a
      cylinder is selected around the axis going through the COM of
      the second group with direction :mdp:`pull-coord1-vec` with
      radius :mdp:`pull-cylinder-r`. Weights of the atoms decrease
      continously to zero as the radial distance goes from 0 to
      :mdp:`pull-cylinder-r` (mass weighting is also used). The radial
      dependence gives rise to radial forces on both pull groups.
      Note that the radius should be smaller than half the box size.
      For tilted cylinders they should be even smaller than half the
      box size since the distance of an atom in the reference group
      from the COM of the pull group has both a radial and an axial
      component. This geometry is not supported with constraint
      pulling.

   .. mdp-value:: angle

      Pull along an angle defined by four groups. The angle is
      defined as the angle between two vectors: the vector connecting
      the COM of the first group to the COM of the second group and
      the vector connecting the COM of the third group to the COM of
      the fourth group.

   .. mdp-value:: angle-axis

      As :mdp-value:`pull-coord1-geometry=angle` but the second vector is given by :mdp:`pull-coord1-vec`.
      Thus, only the two groups that define the first vector need to be given.

   .. mdp-value:: dihedral

      Pull along a dihedral angle defined by six groups. These pairwise
      define three vectors: the vector connecting the COM of group 1
      to the COM of group 2, the COM of group 3 to the COM of group 4,
      and the COM of group 5 to the COM group 6. The dihedral angle is
      then defined as the angle between two planes: the plane spanned by the
      the two first vectors and the plane spanned the two last vectors.


.. mdp:: pull-coord1-groups

   The group indices on which this pull coordinate will operate.
   The number of group indices required is geometry dependent.
   The first index can be 0, in which case an
   absolute reference of :mdp:`pull-coord1-origin` is used. With an
   absolute reference the system is no longer translation invariant
   and one should think about what to do with the center of mass
   motion.

.. mdp:: pull-coord1-dim

   (Y Y Y)
   Selects the dimensions that this pull coordinate acts on and that
   are printed to the output files when
   :mdp:`pull-print-components` = :mdp-value:`pull-coord1-start=yes`. With
   :mdp:`pull-coord1-geometry` = :mdp-value:`pull-coord1-geometry=distance`, only Cartesian
   components set to Y contribute to the distance. Thus setting this
   to Y Y N results in a distance in the x/y plane. With other
   geometries all dimensions with non-zero entries in
   :mdp:`pull-coord1-vec` should be set to Y, the values for other
   dimensions only affect the output.

.. mdp:: pull-coord1-origin

   (0.0 0.0 0.0)
   The pull reference position for use with an absolute reference.

.. mdp:: pull-coord1-vec

   (0.0 0.0 0.0)
   The pull direction. :ref:`gmx grompp` normalizes the vector.

.. mdp:: pull-coord1-start

   .. mdp-value:: no

      do not modify :mdp:`pull-coord1-init`

   .. mdp-value:: yes

      add the COM distance of the starting conformation to
      :mdp:`pull-coord1-init`

.. mdp:: pull-coord1-init

   (0.0) \[nm\] / \[deg\]
   The reference distance at t=0.

.. mdp:: pull-coord1-rate

   (0) \[nm/ps\] / \[deg/ps\]
   The rate of change of the reference position.

.. mdp:: pull-coord1-k

   (0) \[kJ mol-1 nm-2\] / \[kJ mol-1 nm-1\] / \[kJ mol-1 rad-2\] / \[kJ mol-1 rad-1\]
   The force constant. For umbrella pulling this is the harmonic force
   constant in kJ mol-1 nm-2 (or kJ mol-1 rad-2 for angles). For constant force pulling this is the
   force constant of the linear potential, and thus the negative (!)
   of the constant force in kJ mol-1 nm-1 (or kJ mol-1 rad-1 for angles).
   Note that for angles the force constant is expressed in terms of radians
   (while :mdp:`pull-coord1-init` and :mdp:`pull-coord1-rate` are expressed in degrees).

.. mdp:: pull-coord1-kB

   (pull-k1) \[kJ mol-1 nm-2\] / \[kJ mol-1 nm-1\] / \[kJ mol-1 rad-2\] / \[kJ mol-1 rad-1\]
   As :mdp:`pull-coord1-k`, but for state B. This is only used when
   :mdp:`free-energy` is turned on. The force constant is then (1 -
   lambda) * :mdp:`pull-coord1-k` + lambda * :mdp:`pull-coord1-kB`.

AWH adaptive biasing
^^^^^^^^^^^^^^^^^^^^

.. mdp:: awh

   .. mdp-value:: no

      No biasing.

   .. mdp-value:: yes

      Adaptively bias a reaction coordinate using the AWH method and estimate
      the corresponding PMF. The PMF and other AWH data are written to energy
      file at an interval set by :mdp:`awh-nstout` and can be extracted with
      the ``gmx awh`` tool. The AWH coordinate can be
      multidimensional and is defined by mapping each dimension to a pull coordinate index.
      This is only allowed if :mdp-value:`pull-coord1-type=external-potential` and
      :mdp:`pull-coord1-potential-provider` = ``awh`` for the concerned pull coordinate
      indices.

.. mdp:: awh-potential

   .. mdp-value:: convolved

      The applied biasing potential is the convolution of the bias function and a
      set of harmonic umbrella potentials (see :mdp-value:`awh-potential=umbrella` below). This results
      in a smooth potential function and force. The resolution of the potential is set
      by the force constant of each umbrella, see :mdp:`awh1-dim1-force-constant`.

   .. mdp-value:: umbrella

      The potential bias is applied by controlling the position of an harmonic potential
      using Monte-Carlo sampling.  The force constant is set with
      :mdp:`awh1-dim1-force-constant`. The umbrella location
      is sampled using Monte-Carlo every :mdp:`awh-nstsample` steps.
      There are no advantages to using an umbrella.
      This option is mainly for comparison and testing purposes.

.. mdp:: awh-share-multisim

   .. mdp-value:: no

      AWH will not share biases across simulations started with
      :ref:`gmx mdrun` option ``-multidir``. The biases will be independent.

   .. mdp-value:: yes

      With :ref:`gmx mdrun` and option ``-multidir`` the bias and PMF estimates
      for biases with :mdp:`awh1-share-group` >0 will be shared across simulations
      with the biases with the same :mdp:`awh1-share-group` value.
      The simulations should have the same AWH settings for sharing to make sense.
      :ref:`gmx mdrun` will check whether the simulations are technically
      compatible for sharing, but the user should check that bias sharing
      physically makes sense.

.. mdp:: awh-seed

   (-1) Random seed for Monte-Carlo sampling the umbrella position,
   where -1 indicates to generate a seed. Only used with
   :mdp-value:`awh-potential=umbrella`.

.. mdp:: awh-nstout

   (100000)
   Number of steps between printing AWH data to the energy file, should be
   a multiple of :mdp:`nstenergy`.

.. mdp:: awh-nstsample

   (10)
   Number of steps between sampling of the coordinate value. This sampling
   is the basis for updating the bias and estimating the PMF and other AWH observables.

.. mdp:: awh-nsamples-update

   (10)
   The number of coordinate samples used for each AWH update.
   The update interval in steps is :mdp:`awh-nstsample` times this value.

.. mdp:: awh-nbias

   (1)
   The number of biases, each acting on its own coordinate.
   The following options should be specified
   for each bias although below only the options for bias number 1 is shown. Options for
   other bias indices are  obtained by replacing '1' by the bias index.

.. mdp:: awh1-error-init

   (10.0) \[kJ mol-1\]
   Estimated initial average error of the PMF for this bias. This value together with the
   given diffusion constant(s) :mdp:`awh1-dim1-diffusion` determine the initial biasing rate.
   The error is obviously not known *a priori*. Only a rough estimate of :mdp:`awh1-error-init`
   is needed however.
   As a  general guideline, leave :mdp:`awh1-error-init` to its default value when starting a new
   simulation. On the other hand, when there is *a priori* knowledge of the PMF (e.g. when
   an initial PMF estimate is provided, see the :mdp:`awh1-user-data` option)
   then :mdp:`awh1-error-init` should reflect that knowledge.

.. mdp:: awh1-growth

   .. mdp-value:: exp-linear

   Each bias keeps a reference weight histogram for the coordinate samples.
   Its size sets the magnitude of the bias function and free energy estimate updates
   (few samples corresponds to large updates and vice versa).
   Thus, its growth rate sets the maximum convergence rate.
   By default, there is an initial stage in which the histogram grows close to exponentially (but slower than the sampling rate).
   In the final stage that follows, the growth rate is linear and equal to the sampling rate (set by :mdp:`awh-nstsample`).
   The initial stage is typically necessary for efficient convergence when starting a new simulation where
   high free energy barriers have not yet been flattened by the bias.

   .. mdp-value:: linear

   As :mdp-value:`awh1-growth=exp-linear` but skip the initial stage. This may be useful if there is *a priori*
   knowledge (see :mdp:`awh1-error-init`) which eliminates the need for an initial stage. This is also
   the setting compatible with :mdp-value:`awh1-target=local-boltzmann`.

.. mdp:: awh1-equilibrate-histogram

   .. mdp-value:: no

      Do not equilibrate histogram.

   .. mdp-value:: yes

      Before entering the initial stage (see :mdp-value:`awh1-growth=exp-linear`), make sure the
      histogram of sampled weights is following the target distribution closely enough (specifically,
      at least 80% of the target region needs to have a local relative error of less than 20%). This
      option would typically only be used when :mdp:`awh1-share-group` > 0
      and the initial configurations poorly represent the target
      distribution.

.. mdp:: awh1-target

   .. mdp-value:: constant

      The bias is tuned towards a constant (uniform) coordinate distribution
      in the defined sampling interval (defined by  \[:mdp:`awh1-dim1-start`, :mdp:`awh1-dim1-end`\]).

   .. mdp-value:: cutoff

      Similar to :mdp-value:`awh1-target=constant`, but the target
      distribution is proportional to 1/(1 + exp(F - :mdp-value:`awh1-target=cutoff`)),
      where F is the free energy relative to the estimated global minimum.
      This provides a smooth switch of a flat target distribution in
      regions with free energy lower than the cut-off to a Boltzmann
      distribution in regions with free energy higher than the cut-off.

   .. mdp-value:: boltzmann

      The target distribution is a Boltzmann distribtution with a scaled beta (inverse temperature)
      factor given by :mdp:`awh1-target-beta-scaling`. *E.g.*, a value of 0.1
      would give the same coordinate distribution as sampling with a simulation temperature
      scaled by 10.

   .. mdp-value:: local-boltzmann

      Same target distribution and use of :mdp:`awh1-target-beta-scaling`
      but the convergence towards the target distribution is inherently local *i.e.*, the rate of
      change of the bias only depends on the local sampling. This local convergence property is
      only compatible with :mdp-value:`awh1-growth=linear`, since for
      :mdp-value:`awh1-growth=exp-linear` histograms are globally rescaled in the initial stage.

.. mdp:: awh1-target-beta-scaling

   [0] \[\]
   For :mdp-value:`awh1-target=boltzmann` and :mdp-value:`awh1-target=local-boltzmann`
   it is the unitless beta scaling factor taking values in (0,1).

.. mdp:: awh1-target-cutoff

   [0] \[kJ mol-1\]
   For :mdp-value:`awh1-target=cutoff` this is the cutoff, should be > 0.

.. mdp:: awh1-user-data

   .. mdp-value:: no

      Initialize the PMF and target distribution with default values.

   .. mdp-value:: yes

      Initialize the PMF and target distribution with user provided data. For :mdp:`awh-nbias` = 1,
      :ref:`gmx mdrun` will expect a file ``awhinit.xvg`` to be present in the run directory.
      For multiple biases, :ref:`gmx mdrun` expects files ``awhinit1.xvg``, ``awhinit2.xvg``, etc.
      The file name can be changed with the ``-awh`` option.
      The first :mdp:`awh1-ndim` columns of
      each input file should contain the coordinate values, such that each row defines a point in
      coordinate space. Column :mdp:`awh1-ndim` + 1 should contain the PMF value for each point.
      The target distribution column can either follow the PMF (column  :mdp:`awh1-ndim` + 2) or
      be in the same column as written by :ref:`gmx awh`.

.. mdp:: awh1-share-group

   .. mdp-value:: 0

      Do not share the bias.

   .. mdp-value:: positive

      Share the bias and PMF estimates within and/or between simulations.
      Within a simulation, the bias will be shared between biases that have the
      same :mdp:`awh1-share-group` index (note that the current code does not support this).
      With :mdp-value:`awh-share-multisim=yes` and
      :ref:`gmx mdrun` option ``-multidir`` the bias will also be shared across simulations.
      Sharing may increase convergence initially, although the starting configurations
      can be critical, especially when sharing between many biases.
      Currently, positive group values should start at 1 and increase
      by 1 for each subsequent bias that is shared.

.. mdp:: awh1-ndim

   (1) \[integer\]
   Number of dimensions of the coordinate, each dimension maps to 1 pull coordinate.
   The following options should be specified for each such dimension. Below only
   the options for dimension number 1 is shown. Options for other dimension indices are
   obtained by replacing '1' by the dimension index.

.. mdp:: awh1-dim1-coord-provider

   .. mdp-value:: pull

      The module providing the reaction coordinate for this dimension.
      Currently AWH can only act on pull coordinates.

.. mdp:: awh1-dim1-coord-index

   (1)
   Index of the pull coordinate defining this coordinate dimension.

.. mdp:: awh1-dim1-force-constant

   (0) \[kJ/mol/nm^2\] or \[kJ/mol/rad^2\]
   Force constant for the (convolved) umbrella potential(s) along this
   coordinate dimension.

.. mdp:: awh1-dim1-start

   (0.0) \[nm\]/\[rad\]
   Start value of the sampling interval along this dimension. The range of allowed
   values depends on the relevant pull geometry (see :mdp:`pull-coord1-geometry`).
   For periodic geometries :mdp:`awh1-dim1-start` greater than :mdp:`awh1-dim1-end`
   is allowed. The interval will then wrap around from +period/2 to -period/2.

.. mdp:: awh1-dim1-end

   (0.0) \[nm\]/\[rad\]
   End value defining the sampling interval together with :mdp:`awh1-dim1-start`.

.. mdp:: awh1-dim1-period

   (0.0) \[nm\]/\[rad\]
   The period of this reaction coordinate, use 0 when the coordinate is not periodic.

.. mdp:: awh1-dim1-diffusion

   (1e-5) \[nm^2/ps\]/\[rad^2/ps\]
   Estimated diffusion constant for this coordinate dimension determining the initial
   biasing rate. This needs only be a rough estimate and should not critically
   affect the results unless it is set to something very low, leading to slow convergence,
   or very high, forcing the system far from equilibrium. Not setting this value
   explicitly generates a warning.

.. mdp:: awh1-dim1-cover-diameter

   (0.0)) \[nm\]/\[rad\]
   Diameter that needs to be sampled by a single simulation around a coordinate value
   before the point is considered covered in the initial stage (see :mdp-value:`awh1-growth=exp-linear`).
   A value > 0  ensures that for each covering there is a continuous transition of this diameter
   across each coordinate value.
   This is trivially true for independent simulations but not for for multiple bias-sharing simulations
   (:mdp:`awh1-share-group`>0).
   For a diameter = 0, covering occurs as soon as the simulations have sampled the whole interval, which
   for many sharing simulations does not guarantee transitions across free energy barriers.
   On the other hand, when the diameter >= the sampling interval length, covering occurs when a single simulation
   has independently sampled the whole interval.

Enforced rotation
^^^^^^^^^^^^^^^^^

These :ref:`mdp` parameters can be used enforce the rotation of a group of atoms,
e.g. a protein subunit. The `reference manual`_ describes in detail 13 different potentials
that can be used to achieve such a rotation.

.. mdp:: rotation

   .. mdp-value:: no

      No enforced rotation will be applied. All enforced rotation options will
      be ignored (and if present in the :ref:`mdp` file, they unfortunately
      generate warnings).

   .. mdp-value:: yes

      Apply the rotation potential specified by :mdp:`rot-type0` to the group of atoms given
      under the :mdp:`rot-group0` option.

.. mdp:: rot-ngroups

   (1)
   Number of rotation groups.

.. mdp:: rot-group0

   Name of rotation group 0 in the index file.

.. mdp:: rot-type0

   (iso)
   Type of rotation potential that is applied to rotation group 0. Can be of of the following:
   ``iso``, ``iso-pf``, ``pm``, ``pm-pf``, ``rm``, ``rm-pf``, ``rm2``, ``rm2-pf``,
   ``flex``, ``flex-t``, ``flex2``, or ``flex2-t``.

.. mdp:: rot-massw0

   (no)
   Use mass weighted rotation group positions.

.. mdp:: rot-vec0

   (1.0 0.0 0.0)
   Rotation vector, will get normalized.

.. mdp:: rot-pivot0

   (0.0 0.0 0.0)
   Pivot point (nm) for the potentials ``iso``, ``pm``, ``rm``, and ``rm2``.

.. mdp:: rot-rate0

   (0)
   Reference rotation rate (degree/ps) of group 0.

.. mdp:: rot-k0

   (0)
   Force constant (kJ/(mol*nm^2)) for group 0.

.. mdp:: rot-slab-dist0

   (1.5)
   Slab distance (nm), if a flexible axis rotation type was chosen.

.. mdp:: rot-min-gauss0

   (0.001)
   Minimum value (cutoff) of Gaussian function for the force to be evaluated
   (for the flexible axis potentials).

.. mdp:: rot-eps0

   (0.0001)
   Value of additive constant epsilon' (nm^2) for ``rm2*`` and ``flex2*`` potentials.

.. mdp:: rot-fit-method0

   (rmsd)
   Fitting method when determining the actual angle of a rotation group
   (can be one of ``rmsd``, ``norm``, or ``potential``).

.. mdp:: rot-potfit-nsteps0

   (21)
   For fit type ``potential``, the number of angular positions around the reference angle for which the
   rotation potential is evaluated.

.. mdp:: rot-potfit-step0

   (0.25)
   For fit type ``potential``, the distance in degrees between two angular positions.

.. mdp:: rot-nstrout

   (100)
   Output frequency (in steps) for the angle of the rotation group, as well as for the torque
   and the rotation potential energy.

.. mdp:: rot-nstsout

   (1000)
   Output frequency for per-slab data of the flexible axis potentials, i.e. angles, torques and slab centers.


NMR refinement
^^^^^^^^^^^^^^

.. mdp:: disre

   .. mdp-value:: no

      ignore distance restraint information in topology file

   .. mdp-value:: simple

      simple (per-molecule) distance restraints.

   .. mdp-value:: ensemble

      distance restraints over an ensemble of molecules in one
      simulation box. Normally, one would perform ensemble averaging
      over multiple subsystems, each in a separate box, using ``mdrun
      -multi``. Supply ``topol0.tpr``, ``topol1.tpr``, ... with
      different coordinates and/or velocities. The environment
      variable ``GMX_DISRE_ENSEMBLE_SIZE`` sets the number of systems
      within each ensemble (usually equal to the ``mdrun -multi``
      value).

.. mdp:: disre-weighting

   .. mdp-value:: equal

      divide the restraint force equally over all atom pairs in the
      restraint

   .. mdp-value:: conservative

      the forces are the derivative of the restraint potential, this
      results in an weighting of the atom pairs to the reciprocal
      seventh power of the displacement. The forces are conservative
      when :mdp:`disre-tau` is zero.

.. mdp:: disre-mixed

   .. mdp-value:: no

      the violation used in the calculation of the restraint force is
      the time-averaged violation

   .. mdp-value:: yes

      the violation used in the calculation of the restraint force is
      the square root of the product of the time-averaged violation
      and the instantaneous violation

.. mdp:: disre-fc

   (1000) \[kJ mol-1 nm-2\]
   force constant for distance restraints, which is multiplied by a
   (possibly) different factor for each restraint given in the `fac`
   column of the interaction in the topology file.

.. mdp:: disre-tau

   (0) \[ps\]
   time constant for distance restraints running average. A value of
   zero turns off time averaging.

.. mdp:: nstdisreout

   (100) \[steps\]
   period between steps when the running time-averaged and
   instantaneous distances of all atom pairs involved in restraints
   are written to the energy file (can make the energy file very
   large)

.. mdp:: orire

   .. mdp-value:: no

      ignore orientation restraint information in topology file

   .. mdp-value:: yes

      use orientation restraints, ensemble averaging can be performed
      with `mdrun -multi`

.. mdp:: orire-fc

   (0) \[kJ mol\]
   force constant for orientation restraints, which is multiplied by a
   (possibly) different weight factor for each restraint, can be set
   to zero to obtain the orientations from a free simulation

.. mdp:: orire-tau

   (0) \[ps\]
   time constant for orientation restraints running average. A value
   of zero turns off time averaging.

.. mdp:: orire-fitgrp

   fit group for orientation restraining. This group of atoms is used
   to determine the rotation **R** of the system with respect to the
   reference orientation. The reference orientation is the starting
   conformation of the first subsystem. For a protein, backbone is a
   reasonable choice

.. mdp:: nstorireout

   (100) \[steps\]
   period between steps when the running time-averaged and
   instantaneous orientations for all restraints, and the molecular
   order tensor are written to the energy file (can make the energy
   file very large)


Free energy calculations
^^^^^^^^^^^^^^^^^^^^^^^^

.. mdp:: free-energy

   .. mdp-value:: no

      Only use topology A.

   .. mdp-value:: yes

      Interpolate between topology A (lambda=0) to topology B
      (lambda=1) and write the derivative of the Hamiltonian with
      respect to lambda (as specified with :mdp:`dhdl-derivatives`),
      or the Hamiltonian differences with respect to other lambda
      values (as specified with foreign lambda) to the energy file
      and/or to ``dhdl.xvg``, where they can be processed by, for
      example :ref:`gmx bar`. The potentials, bond-lengths and angles
      are interpolated linearly as described in the manual. When
      :mdp:`sc-alpha` is larger than zero, soft-core potentials are
      used for the LJ and Coulomb interactions.

.. mdp:: expanded

   Turns on expanded ensemble simulation, where the alchemical state
   becomes a dynamic variable, allowing jumping between different
   Hamiltonians. See the expanded ensemble options for controlling how
   expanded ensemble simulations are performed. The different
   Hamiltonians used in expanded ensemble simulations are defined by
   the other free energy options.

.. mdp:: init-lambda

   (-1)
   starting value for lambda (float). Generally, this should only be
   used with slow growth (*i.e.* nonzero :mdp:`delta-lambda`). In
   other cases, :mdp:`init-lambda-state` should be specified
   instead. Must be greater than or equal to 0.

.. mdp:: delta-lambda

   (0)
   increment per time step for lambda

.. mdp:: init-lambda-state

   (-1)
   starting value for the lambda state (integer). Specifies which
   columm of the lambda vector (:mdp:`coul-lambdas`,
   :mdp:`vdw-lambdas`, :mdp:`bonded-lambdas`,
   :mdp:`restraint-lambdas`, :mdp:`mass-lambdas`,
   :mdp:`temperature-lambdas`, :mdp:`fep-lambdas`) should be
   used. This is a zero-based index: :mdp:`init-lambda-state` 0 means
   the first column, and so on.

.. mdp:: fep-lambdas

   \[array\]
   Zero, one or more lambda values for which Delta H values will be
   determined and written to dhdl.xvg every :mdp:`nstdhdl`
   steps. Values must be between 0 and 1. Free energy differences
   between different lambda values can then be determined with
   :ref:`gmx bar`. :mdp:`fep-lambdas` is different from the
   other -lambdas keywords because all components of the lambda vector
   that are not specified will use :mdp:`fep-lambdas` (including
   :mdp:`restraint-lambdas` and therefore the pull code restraints).

.. mdp:: coul-lambdas

   \[array\]
   Zero, one or more lambda values for which Delta H values will be
   determined and written to dhdl.xvg every :mdp:`nstdhdl`
   steps. Values must be between 0 and 1. Only the electrostatic
   interactions are controlled with this component of the lambda
   vector (and only if the lambda=0 and lambda=1 states have differing
   electrostatic interactions).

.. mdp:: vdw-lambdas

   \[array\]
   Zero, one or more lambda values for which Delta H values will be
   determined and written to dhdl.xvg every :mdp:`nstdhdl`
   steps. Values must be between 0 and 1. Only the van der Waals
   interactions are controlled with this component of the lambda
   vector.

.. mdp:: bonded-lambdas

   \[array\]
   Zero, one or more lambda values for which Delta H values will be
   determined and written to dhdl.xvg every :mdp:`nstdhdl`
   steps. Values must be between 0 and 1. Only the bonded interactions
   are controlled with this component of the lambda vector.

.. mdp:: restraint-lambdas

   \[array\]
   Zero, one or more lambda values for which Delta H values will be
   determined and written to dhdl.xvg every :mdp:`nstdhdl`
   steps. Values must be between 0 and 1. Only the restraint
   interactions: dihedral restraints, and the pull code restraints are
   controlled with this component of the lambda vector.

.. mdp:: mass-lambdas

   \[array\]
   Zero, one or more lambda values for which Delta H values will be
   determined and written to dhdl.xvg every :mdp:`nstdhdl`
   steps. Values must be between 0 and 1. Only the particle masses are
   controlled with this component of the lambda vector.

.. mdp:: temperature-lambdas

   \[array\]
   Zero, one or more lambda values for which Delta H values will be
   determined and written to dhdl.xvg every :mdp:`nstdhdl`
   steps. Values must be between 0 and 1. Only the temperatures
   controlled with this component of the lambda vector. Note that
   these lambdas should not be used for replica exchange, only for
   simulated tempering.

.. mdp:: calc-lambda-neighbors

   (1)
   Controls the number of lambda values for which Delta H values will
   be calculated and written out, if :mdp:`init-lambda-state` has
   been set. A positive value will limit the number of lambda points
   calculated to only the nth neighbors of :mdp:`init-lambda-state`:
   for example, if :mdp:`init-lambda-state` is 5 and this parameter
   has a value of 2, energies for lambda points 3-7 will be calculated
   and writen out. A value of -1 means all lambda points will be
   written out. For normal BAR such as with :ref:`gmx bar`, a value of
   1 is sufficient, while for MBAR -1 should be used.

.. mdp:: sc-alpha

   (0)
   the soft-core alpha parameter, a value of 0 results in linear
   interpolation of the LJ and Coulomb interactions

.. mdp:: sc-r-power

   (6)
   the power of the radial term in the soft-core equation. Possible
   values are 6 and 48. 6 is more standard, and is the default. When
   48 is used, then sc-alpha should generally be much lower (between
   0.001 and 0.003).

.. mdp:: sc-coul

   (no)
   Whether to apply the soft-core free energy interaction
   transformation to the Columbic interaction of a molecule. Default
   is no, as it is generally more efficient to turn off the Coulomic
   interactions linearly before turning off the van der Waals
   interactions. Note that it is only taken into account when lambda
   states are used, not with :mdp:`couple-lambda0` /
   :mdp:`couple-lambda1`, and you can still turn off soft-core
   interactions by setting :mdp:`sc-alpha` to 0.

.. mdp:: sc-power

   (0)
   the power for lambda in the soft-core function, only the values 1
   and 2 are supported

.. mdp:: sc-sigma

   (0.3) \[nm\]
   the soft-core sigma for particles which have a C6 or C12 parameter
   equal to zero or a sigma smaller than :mdp:`sc-sigma`

.. mdp:: couple-moltype

   Here one can supply a molecule type (as defined in the topology)
   for calculating solvation or coupling free energies. There is a
   special option ``system`` that couples all molecule types in the
   system. This can be useful for equilibrating a system starting from
   (nearly) random coordinates. :mdp:`free-energy` has to be turned
   on. The Van der Waals interactions and/or charges in this molecule
   type can be turned on or off between lambda=0 and lambda=1,
   depending on the settings of :mdp:`couple-lambda0` and
   :mdp:`couple-lambda1`. If you want to decouple one of several
   copies of a molecule, you need to copy and rename the molecule
   definition in the topology.

.. mdp:: couple-lambda0

   .. mdp-value:: vdw-q

      all interactions are on at lambda=0

   .. mdp-value:: vdw

      the charges are zero (no Coulomb interactions) at lambda=0

   .. mdp-value:: q

      the Van der Waals interactions are turned at lambda=0; soft-core
      interactions will be required to avoid singularities

   .. mdp-value:: none

      the Van der Waals interactions are turned off and the charges
      are zero at lambda=0; soft-core interactions will be required to
      avoid singularities.

.. mdp:: couple-lambda1

   analogous to :mdp:`couple-lambda1`, but for lambda=1

.. mdp:: couple-intramol

   .. mdp-value:: no

      All intra-molecular non-bonded interactions for moleculetype
      :mdp:`couple-moltype` are replaced by exclusions and explicit
      pair interactions. In this manner the decoupled state of the
      molecule corresponds to the proper vacuum state without
      periodicity effects.

   .. mdp-value:: yes

      The intra-molecular Van der Waals and Coulomb interactions are
      also turned on/off. This can be useful for partitioning
      free-energies of relatively large molecules, where the
      intra-molecular non-bonded interactions might lead to
      kinetically trapped vacuum conformations. The 1-4 pair
      interactions are not turned off.

.. mdp:: nstdhdl

   (100)
   the frequency for writing dH/dlambda and possibly Delta H to
   dhdl.xvg, 0 means no ouput, should be a multiple of
   :mdp:`nstcalcenergy`.

.. mdp:: dhdl-derivatives

   (yes)

   If yes (the default), the derivatives of the Hamiltonian with
   respect to lambda at each :mdp:`nstdhdl` step are written
   out. These values are needed for interpolation of linear energy
   differences with :ref:`gmx bar` (although the same can also be
   achieved with the right foreign lambda setting, that may not be as
   flexible), or with thermodynamic integration

.. mdp:: dhdl-print-energy

   (no)

   Include either the total or the potential energy in the dhdl
   file. Options are 'no', 'potential', or 'total'. This information
   is needed for later free energy analysis if the states of interest
   are at different temperatures. If all states are at the same
   temperature, this information is not needed. 'potential' is useful
   in case one is using ``mdrun -rerun`` to generate the ``dhdl.xvg``
   file. When rerunning from an existing trajectory, the kinetic
   energy will often not be correct, and thus one must compute the
   residual free energy from the potential alone, with the kinetic
   energy component computed analytically.

.. mdp:: separate-dhdl-file

   .. mdp-value:: yes

      The free energy values that are calculated (as specified with
      the foreign lambda and :mdp:`dhdl-derivatives` settings) are
      written out to a separate file, with the default name
      ``dhdl.xvg``. This file can be used directly with :ref:`gmx
      bar`.

   .. mdp-value:: no

      The free energy values are written out to the energy output file
      (``ener.edr``, in accumulated blocks at every :mdp:`nstenergy`
      steps), where they can be extracted with :ref:`gmx energy` or
      used directly with :ref:`gmx bar`.

.. mdp:: dh-hist-size

   (0)
   If nonzero, specifies the size of the histogram into which the
   Delta H values (specified with foreign lambda) and the derivative
   dH/dl values are binned, and written to ener.edr. This can be used
   to save disk space while calculating free energy differences. One
   histogram gets written for each foreign lambda and two for the
   dH/dl, at every :mdp:`nstenergy` step. Be aware that incorrect
   histogram settings (too small size or too wide bins) can introduce
   errors. Do not use histograms unless you're certain you need it.

.. mdp:: dh-hist-spacing

   (0.1)
   Specifies the bin width of the histograms, in energy units. Used in
   conjunction with :mdp:`dh-hist-size`. This size limits the
   accuracy with which free energies can be calculated. Do not use
   histograms unless you're certain you need it.


Expanded Ensemble calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. mdp:: nstexpanded

   The number of integration steps beween attempted moves changing the
   system Hamiltonian in expanded ensemble simulations. Must be a
   multiple of :mdp:`nstcalcenergy`, but can be greater or less than
   :mdp:`nstdhdl`.

.. mdp:: lmc-stats

   .. mdp-value:: no

      No Monte Carlo in state space is performed.

   .. mdp-value:: metropolis-transition

      Uses the Metropolis weights to update the expanded ensemble
      weight of each state. Min{1,exp(-(beta_new u_new - beta_old
      u_old)}

   .. mdp-value:: barker-transition

      Uses the Barker transition critera to update the expanded
      ensemble weight of each state i, defined by exp(-beta_new
      u_new)/(exp(-beta_new u_new)+exp(-beta_old u_old))

   .. mdp-value:: wang-landau

      Uses the Wang-Landau algorithm (in state space, not energy
      space) to update the expanded ensemble weights.

   .. mdp-value:: min-variance

      Uses the minimum variance updating method of Escobedo et al. to
      update the expanded ensemble weights. Weights will not be the
      free energies, but will rather emphasize states that need more
      sampling to give even uncertainty.

.. mdp:: lmc-mc-move

   .. mdp-value:: no

      No Monte Carlo in state space is performed.

   .. mdp-value:: metropolis-transition

      Randomly chooses a new state up or down, then uses the
      Metropolis critera to decide whether to accept or reject:
      Min{1,exp(-(beta_new u_new - beta_old u_old)}

   .. mdp-value:: barker-transition

      Randomly chooses a new state up or down, then uses the Barker
      transition critera to decide whether to accept or reject:
      exp(-beta_new u_new)/(exp(-beta_new u_new)+exp(-beta_old u_old))

   .. mdp-value:: gibbs

       Uses the conditional weights of the state given the coordinate
       (exp(-beta_i u_i) / sum_k exp(beta_i u_i) to decide which state
       to move to.

   .. mdp-value:: metropolized-gibbs

       Uses the conditional weights of the state given the coordinate
       (exp(-beta_i u_i) / sum_k exp(beta_i u_i) to decide which state
       to move to, EXCLUDING the current state, then uses a rejection
       step to ensure detailed balance. Always more efficient that
       Gibbs, though only marginally so in many situations, such as
       when only the nearest neighbors have decent phase space
       overlap.

.. mdp:: lmc-seed

   (-1)
   random seed to use for Monte Carlo moves in state space. When
   :mdp:`lmc-seed` is set to -1, a pseudo random seed is us

.. mdp:: mc-temperature

   Temperature used for acceptance/rejection for Monte Carlo moves. If
   not specified, the temperature of the simulation specified in the
   first group of :mdp:`ref-t` is used.

.. mdp:: wl-ratio

   (0.8)
   The cutoff for the histogram of state occupancies to be reset, and
   the free energy incrementor to be changed from delta to delta *
   :mdp:`wl-scale`. If we define the Nratio = (number of samples at
   each histogram) / (average number of samples at each
   histogram). :mdp:`wl-ratio` of 0.8 means that means that the
   histogram is only considered flat if all Nratio > 0.8 AND
   simultaneously all 1/Nratio > 0.8.

.. mdp:: wl-scale

   (0.8)
   Each time the histogram is considered flat, then the current value
   of the Wang-Landau incrementor for the free energies is multiplied
   by :mdp:`wl-scale`. Value must be between 0 and 1.

.. mdp:: init-wl-delta

   (1.0)
   The initial value of the Wang-Landau incrementor in kT. Some value
   near 1 kT is usually most efficient, though sometimes a value of
   2-3 in units of kT works better if the free energy differences are
   large.

.. mdp:: wl-oneovert

   (no)
   Set Wang-Landau incrementor to scale with 1/(simulation time) in
   the large sample limit. There is significant evidence that the
   standard Wang-Landau algorithms in state space presented here
   result in free energies getting 'burned in' to incorrect values
   that depend on the initial state. when :mdp:`wl-oneovert` is true,
   then when the incrementor becomes less than 1/N, where N is the
   mumber of samples collected (and thus proportional to the data
   collection time, hence '1 over t'), then the Wang-Lambda
   incrementor is set to 1/N, decreasing every step. Once this occurs,
   :mdp:`wl-ratio` is ignored, but the weights will still stop
   updating when the equilibration criteria set in
   :mdp:`lmc-weights-equil` is achieved.

.. mdp:: lmc-repeats

   (1)
   Controls the number of times that each Monte Carlo swap type is
   performed each iteration. In the limit of large numbers of Monte
   Carlo repeats, then all methods converge to Gibbs sampling. The
   value will generally not need to be different from 1.

.. mdp:: lmc-gibbsdelta

   (-1)
   Limit Gibbs sampling to selected numbers of neighboring states. For
   Gibbs sampling, it is sometimes inefficient to perform Gibbs
   sampling over all of the states that are defined. A positive value
   of :mdp:`lmc-gibbsdelta` means that only states plus or minus
   :mdp:`lmc-gibbsdelta` are considered in exchanges up and down. A
   value of -1 means that all states are considered. For less than 100
   states, it is probably not that expensive to include all states.

.. mdp:: lmc-forced-nstart

   (0)
   Force initial state space sampling to generate weights. In order to
   come up with reasonable initial weights, this setting allows the
   simulation to drive from the initial to the final lambda state,
   with :mdp:`lmc-forced-nstart` steps at each state before moving on
   to the next lambda state. If :mdp:`lmc-forced-nstart` is
   sufficiently long (thousands of steps, perhaps), then the weights
   will be close to correct. However, in most cases, it is probably
   better to simply run the standard weight equilibration algorithms.

.. mdp:: nst-transition-matrix

   (-1)
   Frequency of outputting the expanded ensemble transition matrix. A
   negative number means it will only be printed at the end of the
   simulation.

.. mdp:: symmetrized-transition-matrix

   (no)
   Whether to symmetrize the empirical transition matrix. In the
   infinite limit the matrix will be symmetric, but will diverge with
   statistical noise for short timescales. Forced symmetrization, by
   using the matrix T_sym = 1/2 (T + transpose(T)), removes problems
   like the existence of (small magnitude) negative eigenvalues.

.. mdp:: mininum-var-min

   (100)
   The min-variance strategy (option of :mdp:`lmc-stats` is only
   valid for larger number of samples, and can get stuck if too few
   samples are used at each state. :mdp:`mininum-var-min` is the
   minimum number of samples that each state that are allowed before
   the min-variance strategy is activated if selected.

.. mdp:: init-lambda-weights

   The initial weights (free energies) used for the expanded ensemble
   states. Default is a vector of zero weights. format is similar to
   the lambda vector settings in :mdp:`fep-lambdas`, except the
   weights can be any floating point number. Units are kT. Its length
   must match the lambda vector lengths.

.. mdp:: lmc-weights-equil

   .. mdp-value:: no

      Expanded ensemble weights continue to be updated throughout the
      simulation.

   .. mdp-value:: yes

      The input expanded ensemble weights are treated as equilibrated,
      and are not updated throughout the simulation.

   .. mdp-value:: wl-delta

      Expanded ensemble weight updating is stopped when the
      Wang-Landau incrementor falls below this value.

   .. mdp-value:: number-all-lambda

      Expanded ensemble weight updating is stopped when the number of
      samples at all of the lambda states is greater than this value.

   .. mdp-value:: number-steps

      Expanded ensemble weight updating is stopped when the number of
      steps is greater than the level specified by this value.

   .. mdp-value:: number-samples

      Expanded ensemble weight updating is stopped when the number of
      total samples across all lambda states is greater than the level
      specified by this value.

   .. mdp-value:: count-ratio

      Expanded ensemble weight updating is stopped when the ratio of
      samples at the least sampled lambda state and most sampled
      lambda state greater than this value.

.. mdp:: simulated-tempering

   (no)
   Turn simulated tempering on or off. Simulated tempering is
   implemented as expanded ensemble sampling with different
   temperatures instead of different Hamiltonians.

.. mdp:: sim-temp-low

   (300) \[K\]
   Low temperature for simulated tempering.

.. mdp:: sim-temp-high

   (300) \[K\]
   High temperature for simulated tempering.

.. mdp:: simulated-tempering-scaling

   Controls the way that the temperatures at intermediate lambdas are
   calculated from the :mdp:`temperature-lambdas` part of the lambda
   vector.

   .. mdp-value:: linear

      Linearly interpolates the temperatures using the values of
      :mdp:`temperature-lambdas`, *i.e.* if :mdp:`sim-temp-low`
      =300, :mdp:`sim-temp-high` =400, then lambda=0.5 correspond to
      a temperature of 350. A nonlinear set of temperatures can always
      be implemented with uneven spacing in lambda.

   .. mdp-value:: geometric

      Interpolates temperatures geometrically between
      :mdp:`sim-temp-low` and :mdp:`sim-temp-high`. The i:th state
      has temperature :mdp:`sim-temp-low` * (:mdp:`sim-temp-high` /
      :mdp:`sim-temp-low`) raised to the power of
      (i/(ntemps-1)). This should give roughly equal exchange for
      constant heat capacity, though of course things simulations that
      involve protein folding have very high heat capacity peaks.

   .. mdp-value:: exponential

      Interpolates temperatures exponentially between
      :mdp:`sim-temp-low` and :mdp:`sim-temp-high`. The i:th state
      has temperature :mdp:`sim-temp-low` + (:mdp:`sim-temp-high` -
      :mdp:`sim-temp-low`)*((exp(:mdp:`temperature-lambdas`
      (i))-1)/(exp(1.0)-i)).


Non-equilibrium MD
^^^^^^^^^^^^^^^^^^

.. mdp:: acc-grps

   groups for constant acceleration (*e.g.* ``Protein Sol``) all atoms
   in groups Protein and Sol will experience constant acceleration as
   specified in the :mdp:`accelerate` line

.. mdp:: accelerate

   (0) \[nm ps^-2\]
   acceleration for :mdp:`acc-grps`; x, y and z for each group
   (*e.g.* ``0.1 0.0 0.0 -0.1 0.0 0.0`` means that first group has
   constant acceleration of 0.1 nm ps-2 in X direction, second group
   the opposite).

.. mdp:: freezegrps

   Groups that are to be frozen (*i.e.* their X, Y, and/or Z position
   will not be updated; *e.g.* ``Lipid SOL``). :mdp:`freezedim`
   specifies for which dimension the freezing applies. To avoid
   spurious contibrutions to the virial and pressure due to large
   forces between completely frozen atoms you need to use energy group
   exclusions, this also saves computing time. Note that coordinates
   of frozen atoms are not scaled by pressure-coupling algorithms.

.. mdp:: freezedim

   dimensions for which groups in :mdp:`freezegrps` should be frozen,
   specify `Y` or `N` for X, Y and Z and for each group (*e.g.* ``Y Y
   N N N N`` means that particles in the first group can move only in
   Z direction. The particles in the second group can move in any
   direction).

.. mdp:: cos-acceleration

   (0) \[nm ps^-2\]
   the amplitude of the acceleration profile for calculating the
   viscosity. The acceleration is in the X-direction and the magnitude
   is :mdp:`cos-acceleration` cos(2 pi z/boxheight). Two terms are
   added to the energy file: the amplitude of the velocity profile and
   1/viscosity.

.. mdp:: deform

   (0 0 0 0 0 0) \[nm ps-1\]
   The velocities of deformation for the box elements: a(x) b(y) c(z)
   b(x) c(x) c(y). Each step the box elements for which :mdp:`deform`
   is non-zero are calculated as: box(ts)+(t-ts)*deform, off-diagonal
   elements are corrected for periodicity. The coordinates are
   transformed accordingly. Frozen degrees of freedom are (purposely)
   also transformed. The time ts is set to t at the first step and at
   steps at which x and v are written to trajectory to ensure exact
   restarts. Deformation can be used together with semiisotropic or
   anisotropic pressure coupling when the appropriate
   compressibilities are set to zero. The diagonal elements can be
   used to strain a solid. The off-diagonal elements can be used to
   shear a solid or a liquid.


Electric fields
.. mdp:: electric-field-x ; electric-field-y ; electric-field-z

   Here you can specify an electric field that optionally can be
   alternating and pulsed. The general expression for the field
   has the form of a gaussian laser pulse:

   E(t) = E0 exp ( -(t-t0)^2/(2 sigma^2) ) cos(omega (t-t0))

   For example, the four parameters for direction x are set in the
   three fields of :mdp:`electric-field-x` (and similar for y and z) 
   like

   electric-field-x  = E0 omega t0 sigma

   In the special case that sigma = 0, the exponential term is omitted
   and only the cosine term is used. If also omega = 0 a static
   electric field is applied.

   More details in Carl Caleman and David van der Spoel: Picosecond
   Melting of Ice by an Infrared Laser Pulse - A Simulation Study
   Angew. Chem. Intl. Ed. 47 pp. 14 17-1420 (2008)



Mixed quantum/classical molecular dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. MDP:: QMMM

   .. mdp-value:: no

      No QM/MM.

   .. mdp-value:: yes

      Do a QM/MM simulation. Several groups can be described at
      different QM levels separately. These are specified in the
      :mdp:`QMMM-grps` field separated by spaces. The level of *ab
      initio* theory at which the groups are described is specified by
      :mdp:`QMmethod` and :mdp:`QMbasis` Fields. Describing the
      groups at different levels of theory is only possible with the
      ONIOM QM/MM scheme, specified by :mdp:`QMMMscheme`.

.. mdp:: QMMM-grps

   groups to be descibed at the QM level

.. mdp:: QMMMscheme

   .. mdp-value:: normal

      normal QM/MM. There can only be one :mdp:`QMMM-grps` that is
      modelled at the :mdp:`QMmethod` and :mdp:`QMbasis` level of
      *ab initio* theory. The rest of the system is described at the
      MM level. The QM and MM subsystems interact as follows: MM point
      charges are included in the QM one-electron hamiltonian and all
      Lennard-Jones interactions are described at the MM level.

   .. mdp-value:: ONIOM

      The interaction between the subsystem is described using the
      ONIOM method by Morokuma and co-workers. There can be more than
      one :mdp:`QMMM-grps` each modeled at a different level of QM
      theory (:mdp:`QMmethod` and :mdp:`QMbasis`).

.. mdp:: QMmethod

   (RHF)
   Method used to compute the energy and gradients on the QM
   atoms. Available methods are AM1, PM3, RHF, UHF, DFT, B3LYP, MP2,
   CASSCF, and MMVB. For CASSCF, the number of electrons and orbitals
   included in the active space is specified by :mdp:`CASelectrons`
   and :mdp:`CASorbitals`.

.. mdp:: QMbasis

   (STO-3G)
   Basis set used to expand the electronic wavefuntion. Only Gaussian
   basis sets are currently available, *i.e.* ``STO-3G, 3-21G, 3-21G*,
   3-21+G*, 6-21G, 6-31G, 6-31G*, 6-31+G*,`` and ``6-311G``.

.. mdp:: QMcharge

   (0) \[integer\]
   The total charge in `e` of the :mdp:`QMMM-grps`. In case there are
   more than one :mdp:`QMMM-grps`, the total charge of each ONIOM
   layer needs to be specified separately.

.. mdp:: QMmult

   (1) \[integer\]
   The multiplicity of the :mdp:`QMMM-grps`. In case there are more
   than one :mdp:`QMMM-grps`, the multiplicity of each ONIOM layer
   needs to be specified separately.

.. mdp:: CASorbitals

   (0) \[integer\]
   The number of orbitals to be included in the active space when
   doing a CASSCF computation.

.. mdp:: CASelectrons

   (0) \[integer\]
   The number of electrons to be included in the active space when
   doing a CASSCF computation.

.. MDP:: SH

   .. mdp-value:: no

      No surface hopping. The system is always in the electronic
      ground-state.

   .. mdp-value:: yes

      Do a QM/MM MD simulation on the excited state-potential energy
      surface and enforce a *diabatic* hop to the ground-state when
      the system hits the conical intersection hyperline in the course
      the simulation. This option only works in combination with the
      CASSCF method.


Implicit solvent
^^^^^^^^^^^^^^^^

.. mdp:: implicit-solvent

   .. mdp-value:: no

      No implicit solvent

   .. mdp-value:: GBSA

      Do a simulation with implicit solvent using the Generalized Born
      formalism. Three different methods for calculating the Born
      radii are available, Still, HCT and OBC. These are specified
      with the :mdp:`gb-algorithm` field. The non-polar solvation is
      specified with the :mdp:`sa-algorithm` field.

.. mdp:: gb-algorithm

   .. mdp-value:: Still

      Use the Still method to calculate the Born radii

   .. mdp-value:: HCT

      Use the Hawkins-Cramer-Truhlar method to calculate the Born
      radii

   .. mdp-value:: OBC

      Use the Onufriev-Bashford-Case method to calculate the Born
      radii

.. mdp:: nstgbradii

   (1) \[steps\]
   Frequency to (re)-calculate the Born radii. For most practial
   purposes, setting a value larger than 1 violates energy
   conservation and leads to unstable trajectories.

.. mdp:: rgbradii

   (1.0) \[nm\]
   Cut-off for the calculation of the Born radii. Currently must be
   equal to rlist

.. mdp:: gb-epsilon-solvent

   (80)
   Dielectric constant for the implicit solvent

.. mdp:: gb-saltconc

   (0) \[M\]
   Salt concentration for implicit solvent models, currently not used

.. mdp:: gb-obc-alpha
.. mdp:: gb-obc-beta
.. mdp:: gb-obc-gamma

   Scale factors for the OBC model. Default values of 1, 0.78 and 4.85
   respectively are for OBC(II). Values for OBC(I) are 0.8, 0 and 2.91
   respectively

.. mdp:: gb-dielectric-offset

   (0.009) \[nm\]
   Distance for the di-electric offset when calculating the Born
   radii. This is the offset between the center of each atom the
   center of the polarization energy for the corresponding atom

.. mdp:: sa-algorithm

   .. mdp-value:: Ace-approximation

      Use an Ace-type approximation

   .. mdp-value:: None

      No non-polar solvation calculation done. For GBSA only the polar
      part gets calculated

.. mdp:: sa-surface-tension

   \[kJ mol-1 nm-2\]
   Default value for surface tension with SA algorithms. The default
   value is -1; Note that if this default value is not changed it will
   be overridden by :ref:`gmx grompp` using values that are specific
   for the choice of radii algorithm (0.0049 kcal/mol/Angstrom^2 for
   Still, 0.0054 kcal/mol/Angstrom2 for HCT/OBC) Setting it to 0 will
   while using an sa-algorithm other than None means no non-polar
   calculations are done.


Computational Electrophysiology
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use these options to switch on and control ion/water position exchanges in "Computational
Electrophysiology" simulation setups. (See the `reference manual`_ for details).

.. mdp:: swapcoords

   .. mdp-value:: no

      Do not enable ion/water position exchanges.

   .. mdp-value:: X ; Y ; Z

      Allow for ion/water position exchanges along the chosen direction.
      In a typical setup with the membranes parallel to the x-y plane,
      ion/water pairs need to be exchanged in Z direction to sustain the
      requested ion concentrations in the compartments.

.. mdp:: swap-frequency

   (1) The swap attempt frequency, i.e. every how many time steps the ion counts
   per compartment are determined and exchanges made if necessary.
   Normally it is not necessary to check at every time step.
   For typical Computational Electrophysiology setups, a value of about 100 is
   sufficient and yields a negligible performance impact.

.. mdp:: split-group0

   Name of the index group of the membrane-embedded part of channel #0.
   The center of mass of these atoms defines one of the compartment boundaries
   and should be chosen such that it is near the center of the membrane.

.. mdp:: split-group1

   Channel #1 defines the position of the other compartment boundary.

.. mdp:: massw-split0

   (no) Defines whether or not mass-weighting is used to calculate the split group center.

   .. mdp-value:: no

      Use the geometrical center.

   .. mdp-value:: yes

      Use the center of mass.

.. mdp:: massw-split1

   (no) As above, but for split-group #1.

.. mdp:: solvent-group

   Name of the index group of solvent molecules.

.. mdp:: coupl-steps

   (\10) Average the number of ions per compartment over these many swap attempt steps.
   This can be used to prevent that ions near a compartment boundary
   (diffusing through a channel, e.g.) lead to unwanted back and forth swaps.

.. mdp:: iontypes

   (1) The number of different ion types to be controlled. These are during the
   simulation exchanged with solvent molecules to reach the desired reference numbers.

.. mdp:: iontype0-name

   Name of the first ion type.

.. mdp:: iontype0-in-A

   (-1) Requested (=reference) number of ions of type 0 in compartment A.
   The default value of -1 means: use the number of ions as found in time step 0
   as reference value.

.. mdp:: iontype0-in-B

   (-1) Reference number of ions of type 0 for compartment B.

.. mdp:: bulk-offsetA

   (0.0) Offset of the first swap layer from the compartment A midplane.
   By default (i.e. bulk offset = 0.0), ion/water exchanges happen between layers
   at maximum distance (= bulk concentration) to the split group layers. However,
   an offset b (-1.0 < b < +1.0) can be specified to offset the bulk layer from the middle at 0.0
   towards one of the compartment-partitioning layers (at +/- 1.0).

.. mdp:: bulk-offsetB

   (0.0) Offset of the other swap layer from the compartment B midplane.


.. mdp:: threshold

   (\1) Only swap ions if threshold difference to requested count is reached.

.. mdp:: cyl0-r

   (2.0) \[nm\] Radius of the split cylinder #0.
   Two split cylinders (mimicking the channel pores) can optionally be defined
   relative to the center of the split group. With the help of these cylinders
   it can be counted which ions have passed which channel. The split cylinder
   definition has no impact on whether or not ion/water swaps are done.

.. mdp:: cyl0-up

   (1.0) \[nm\] Upper extension of the split cylinder #0.

.. mdp:: cyl0-down

   (1.0) \[nm\] Lower extension of the split cylinder #0.

.. mdp:: cyl1-r

   (2.0) \[nm\] Radius of the split cylinder #1.

.. mdp:: cyl1-up

   (1.0) \[nm\] Upper extension of the split cylinder #1.

.. mdp:: cyl1-down

   (1.0) \[nm\] Lower extension of the split cylinder #1.


User defined thingies
^^^^^^^^^^^^^^^^^^^^^

.. mdp:: user1-grps
.. mdp:: user2-grps
.. mdp:: userint1 (0)
.. mdp:: userint2 (0)
.. mdp:: userint3 (0)
.. mdp:: userint4 (0)
.. mdp:: userreal1 (0)
.. mdp:: userreal2 (0)
.. mdp:: userreal3 (0)
.. mdp:: userreal4 (0)

   These you can use if you modify code. You can pass integers and
   reals and groups to your subroutine. Check the inputrec definition
   in ``src/gromacs/mdtypes/inputrec.h``

Removed features
^^^^^^^^^^^^^^^^

This feature has been removed from |Gromacs|, but so that old
:ref:`mdp` and :ref:`tpr` files cannot be mistakenly misused, we still
parse this option. :ref:`gmx grompp` and :ref:`gmx mdrun` will issue a
fatal error if this is set.

.. mdp:: adress

   (no)

.. _reference manual: gmx-manual-parent-dir_
