Molecular dynamics parameters (.mdp options)
============================================

General information
-------------------
Default values are given in parentheses, or listed first among choices. The first option in the list is always the default option. Units are given in square brackets The difference between a dash and an underscore is ignored.

A `sample mdp file`_ is available. This should be appropriate to start a normal simulation. Edit it to suit your specific needs and desires.


Preprocessing
^^^^^^^^^^^^^

.. glossary::

    include
        directories to include in your topology. Format: ``-I/home/john/mylib -I../otherlib``

    define
        defines to pass to the preprocessor, default is no defines. You can use any defines to control options in your customized topology files. Options that act on existing top_ file mechanisms include

            ``-DFLEXIBLE`` will use flexible water instead of rigid water into your topology, this can be useful for normal mode analysis.

            ``-DPOSRES`` will trigger the inclusion of ``posre.itp`` into your topology, used for implementing position restraints.


Run control
^^^^^^^^^^^

.. glossary::

    integrator
        (Despite the name, this list includes algorithms that are not actually integrators over time. :term:`steep` and all entries following it are in this category)

        .. glossary::

            md
                A leap-frog algorithm for integrating Newton's equations of motion.

            md-vv
                A velocity Verlet algorithm for integrating Newton's equations of motion. For constant NVE simulations started from corresponding points in the same trajectory, the trajectories are analytically, but not binary, identical to the :term:`md` leap-frog integrator. The the kinetic energy, which is determined from the whole step velocities and is therefore slightly too high. The advantage of this integrator is more accurate, reversible Nose-Hoover and Parrinello-Rahman coupling integration based on Trotter expansion, as well as (slightly too small) full step velocity output. This all comes at the cost off extra computation, especially with constraints and extra communication in parallel. Note that for nearly all production simulations the :term:`md` integrator is accurate enough.

            md-vv-avek
                A velocity Verlet algorithm identical to :term:`md-vv`, except that the kinetic energy is determined as the average of the two half step kinetic energies as in the :term:`md` integrator, and this thus more accurate. With Nose-Hoover and/or Parrinello-Rahman coupling this comes with a slight increase in computational cost.

            sd
                An accurate and efficient leap-frog stochastic dynamics integrator. With constraints, coordinates needs to be constrained twice per integration step. Depending on the computational cost of the force calculation, this can take a significant part of the simulation time. The temperature for one or more groups of atoms (:term:`tc-grps`) is set with :term:`ref-t`, the inverse friction constant for each group is set with :term:`tau-t`. The parameter :term:`tcoupl` is ignored. The random generator is initialized with :term:`ld-seed`. When used as a thermostat, an appropriate value for :term:`tau-t` is 2 ps, since this results in a friction that is lower than the internal friction of water, while it is high enough to remove excess heat NOTE: temperature deviations decay twice as fast as with a Berendsen thermostat with the same :term:`tau-t`.

            sd2
                This used to be the default sd integrator, but is now deprecated. Four Gaussian random numbers are required per coordinate per step. With constraints, the temperature will be slightly too high.

            bd
                An Euler integrator for Brownian or position Langevin dynamics, the velocity is the force divided by a friction coefficient (:term:`bd-fric`) plus random thermal noise (:term:`ref-t`). When :term:`bd-fric` is 0, the friction coefficient for each particle is calculated as mass/ :term:`tau-t`, as for the integrator :term:`sd`. The random generator is initialized with :term:`ld-seed`.

            steep
                 A steepest descent algorithm for energy minimization. The maximum step size is :term:`emstep`, the tolerance is :term:`emtol`.

            cg
                 A conjugate gradient algorithm for energy minimization, the tolerance is :term:`emtol`. CG is more efficient when a steepest descent step is done every once in a while, this is determined by :term:`nstcgsteep`. For a minimization prior to a normal mode analysis, which requires a very high accuracy, GROMACS should be compiled in double precision.

            l-bfgs
                 A quasi-Newtonian algorithm for energy minimization according to the low-memory Broyden-Fletcher-Goldfarb-Shanno approach. In practice this seems to converge faster than Conjugate Gradients, but due to the correction steps necessary it is not (yet) parallelized.

            nm
                 Normal mode analysis is performed on the structure in the tpr_ file. GROMACS should be compiled in double precision.

            tpi
                 Test particle insertion. The last molecule in the topology is the test particle. A trajectory must be provided to ``mdrun -rerun``. This trajectory should not contain the molecule to be inserted. Insertions are performed :term:`nsteps` times in each frame at random locations and with random orientiations of the molecule. When :term:`nstlist` is larger than one, :term:`nstlist` insertions are performed in a sphere with radius :term:`rtpi` around a the same random location using the same neighborlist (and the same long-range energy when :term:`rvdw` or :term:`rcoulomb` > :term:`rlist`, which is only allowed for single-atom molecules). Since neighborlist construction is expensive, one can perform several extra insertions with the same list almost for free. The random seed is set with :term:`ld-seed`. The temperature for the Boltzmann weighting is set with :term:`ref-t`, this should match the temperature of the simulation of the original trajectory. Dispersion correction is implemented correctly for TPI. All relevant quantities are written to the file specified with ``mdrun -tpi``. The distribution of insertion energies is written to the file specified with ``mdrun -tpid``. No trajectory or energy file is written. Parallel TPI gives identical results to single-node TPI. For charged molecules, using PME with a fine grid is most accurate and also efficient, since the potential in the system only needs to be calculated once per frame.

            tpic
                Test particle insertion into a predefined cavity location. The procedure is the same as for :term:`tpi`, except that one coordinate extra is read from the trajectory, which is used as the insertion location. The molecule to be inserted should be centered at 0,0,0. Gromacs does not do this for you, since for different situations a different way of centering might be optimal. Also :term:`rtpi` sets the radius for the sphere around this location. Neighbor searching is done only once per frame, :term:`nstlist` is not used. Parallel :term:`tpic` gives identical results to single-rank :term:`tpic`.

    tinit
        (0) \[ps\]
        starting time for your run (only makes sense for time-based integrators)

    dt
        (0.001) \[ps\]
        time step for integration (only makes sense for time-based integrators)

    nsteps
        (0)
        maximum number of steps to integrate or minimize, -1 is no maximum

    init-step
        (0)
        The starting step. The time at an step i in a run is calculated as: t = :term:`tinit` + :term:`dt` * (:term:`init-step` + i). The free-energy lambda is calculated as: lambda = :term:`init-lambda` + :term:`delta-lambda` * (:term:`init-step` + i). Also non-equilibrium MD parameters can depend on the step number. Thus for exact restarts or redoing part of a run it might be necessary to set :term:`init-step` to the step number of the restart frame. `gmx convert-tpr`_ does this automatically.

    comm-mode
        Linear
            Remove center of mass translation

        Angular
            Remove center of mass translation and rotation around the center of mass

        None
            No restriction on the center of mass motion

    nstcomm
        (100) \[steps\]
        frequency for center of mass motion removal

    comm-grps
        group(s) for center of mass motion removal, default is the whole system


Langevin dynamics
^^^^^^^^^^^^^^^^^

.. glossary::

    bd-fric
        (0) \[amu ps-1\]
        Brownian dynamics friction coefficient. When :term:`bd-fric` is 0, the friction coefficient for each particle is calculated as mass/ :term:`tau-t`.

    ld-seed
        (-1) \[integer\]
        used to initialize random generator for thermal noise for stochastic and Brownian dynamics. When :term:`ld-seed` is set to -1, a pseudo random seed is used. When running BD or SD on multiple processors, each processor uses a seed equal to :term:`ld-seed` plus the processor number.


Energy minimization
^^^^^^^^^^^^^^^^^^^

.. glossary::

    emtol
        (10.0) \[kJ mol-1 nm-1\]
        the minimization is converged when the maximum force is smaller than this value

    emstep
        (0.01) \[nm\]
        initial step-size

    nstcgsteep
        (1000) \[steps\]
        frequency of performing 1 steepest descent step while doing conjugate gradient energy minimization.

    nbfgscorr
        (10)
        Number of correction steps to use for L-BFGS minimization. A higher number is (at least theoretically) more accurate, but slower.


Shell Molecular Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^

When shells or flexible constraints are present in the system the positions of the shells and the lengths of the flexible constraints are optimized at every
time step until either the RMS force on the shells and constraints is less than emtol, or a maximum number of iterations :term:`niter` has been reached

.. glossary::

    emtol
        (10.0) \[kJ mol-1 nm-1\]
        the minimization is converged when the maximum force is smaller than this value. For shell MD this value should be 1.0 at most, but since the variable is used for energy minimization as well the default is 10.0.

    niter
        (20)
        maximum number of iterations for optimizing the shell positions and the flexible constraints.

    fcstep
        (0) \[ps^2\]
        the step size for optimizing the flexible constraints. Should be chosen as mu/(d2V/dq2) where mu is the reduced mass of two particles in a flexible constraint and d2V/dq2 is the second derivative of the potential in the constraint direction. Hopefully this number does not differ too much between the flexible constraints, as the number of iterations and thus the runtime is very sensitive to fcstep. Try several values!


Test particle insertion
^^^^^^^^^^^^^^^^^^^^^^^

.. glossary::
    rtpi
        (0.05) \[nm\]
        the test particle insertion radius, see integrators :term:`tpi` and :term:`tpic`


Output control
^^^^^^^^^^^^^^

.. glossary::

    nstxout
        (0) \[steps\]
        number of steps that elapse between writing coordinates to output trajectory file, the last coordinates are always written

    nstvout
        (0) \[steps\]
        number of steps that elapse between writing velocities to output trajectory, the last velocities are always written

    nstfout
        (0) \[steps\]
        number of steps that elapse between writing forces to output trajectory.

    nstlog
        (1000) \[steps\]
        number of steps that elapse between writing energies to the log file, the last energies are always written

    nstcalcenergy
        (100)
        number of steps that elapse between calculating the energies, 0 is never. This option is only relevant with dynamics. With a twin-range cut-off setup :term:`nstcalcenergy` should be equal to or a multiple of :term:`nstlist`. This option affects the performance in parallel simulations, because calculating energies requires global communication between all processes which can become a bottleneck at high parallelization.

    nstenergy
        (1000) \[steps\]
        number of steps that else between writing energies to energy file, the last energies are always written, should be a multiple of :term:`nstcalcenergy`. Note that the exact sums and fluctuations over all MD steps modulo :term:`nstcalcenergy` are stored in the energy file, so `gmx energy`_ can report exact energy averages and fluctuations also when :term:`nstenergy` > 1

    nstxout-compressed
        (0) \[steps\]
        number of steps that elapse between writing position coordinates using lossy compression

    compressed-x-precision
        (1000) \[real\]
        precision with which to write to the compressed trajectory file

    compressed-x-grps
        group(s) to write to the compressed trajectory file, by default the whole system is written (if :term:`nstxout-compressed` > 0)

    energygrps
        group(s) to write to energy file


Neighbor searching
^^^^^^^^^^^^^^^^^^

.. glossary::

    cutoff-scheme

        .. glossary::

            Verlet
                Generate a pair list with buffering. The buffer size is automatically set based on :term:`verlet-buffer-tolerance`, unless this is set to -1, in which case :term:`rlist` will be used. This option has an explicit, exact cut-off at :term:`rvdw` equal to :term:`rcoulomb`. Currently only cut-off, reaction-field, PME electrostatics and plain LJ are supported. Some mdrun_ functionality is not yet supported with the :term:`Verlet` scheme, but `gmx grompp`_ checks for this. Native GPU acceleration is only supported with :term:`Verlet`. With GPU-accelerated PME or with separate PME ranks, mdrun_ will automatically tune the CPU/GPU load balance by scaling :term:`rcoulomb` and the grid spacing. This can be turned off with ``mdrun -notunepme``. :term:`Verlet` is faster than :term:`group` when there is no water, or if :term:`group` would use a pair-list buffer to conserve energy.

            group
                Generate a pair list for groups of atoms. These groups correspond to the charge groups in the topology. This was the only cut-off treatment scheme before version 4.6, and is **deprecated in 5.0**. There is no explicit buffering of the pair list. This enables efficient force calculations for water, but energy is only conserved when a buffer is explicitly added.

    nstlist
        \(10) \[steps\]

        >0
            Frequency to update the neighbor list (and the long-range forces, when using twin-range cut-offs). When this is 0, the neighbor list is made only once. With energy minimization the neighborlist will be updated for every energy evaluation when :term:`nstlist` is greater than 0. With :term:`Verlet` and :term:`verlet-buffer-tolerance` set, :term:`nstlist` is actually a minimum value and mdrun_ might increase it, unless it is set to 1. With parallel simulations and/or non-bonded force calculation on the GPU, a value of 20 or 40 often gives the best performance. With :term:`group` and non-exact cut-off's, :term:`nstlist` will affect the accuracy of your simulation and it can not be chosen freely.

        0
            The neighbor list is only constructed once and never updated. This is mainly useful for vacuum simulations in which all particles see each other.

        <0
            Unused.

    nstcalclr
        (-1) \[steps\]
        Controls the period between calculations of long-range forces when using the group cut-off scheme.

        1
            Calculate the long-range forces every single step. This is useful to have separate neighbor lists with buffers for electrostatics and Van der Waals interactions, and in particular it makes it possible to have the Van der Waals cutoff longer than electrostatics (useful *e.g.* with PME). However, there is no point in having identical long-range cutoffs for both interaction forms and update them every step - then it will be slightly faster to put everything in the short-range list.

        >1
            Calculate the long-range forces every :term:`nstcalclr` steps and use a multiple-time-step integrator to combine forces. This can now be done more frequently than :term:`nstlist` since the lists are stored, and it might be a good idea *e.g.* for Van der Waals interactions that vary slower than electrostatics.

        \-1
            Calculate long-range forces on steps where neighbor searching is performed. While this is the default value, you might want to consider updating the long-range forces more frequently.

        Note that twin-range force evaluation might be enabled automatically by PP-PME load balancing. This is done in order to maintain the chosen Van der Waals interaction radius even if the load balancing is changing the electrostatics cutoff. If the mdp_ file already specifies twin-range interactions (*e.g.* to evaluate Lennard-Jones interactions with a longer cutoff than the PME electrostatics every 2-3 steps), the load balancing will have also a small effect on Lennard-Jones, since the short-range cutoff (inside which forces are evaluated every step) is changed.

    ns-type
        grid
            Make a grid in the box and only check atoms in neighboring grid cells when constructing a new neighbor list every :term:`nstlist` steps. In large systems grid search is much faster than simple search.

        simple
            Check every atom in the box when constructing a new neighbor list every :term:`nstlist` steps (only with :term:`group` cut-off scheme).

    pbc
        xyz
            Use periodic boundary conditions in all directions.

        no
            Use no periodic boundary conditions, ignore the box. To simulate without cut-offs, set all cut-offs and :term:`nstlist` to 0. For best performance without cut-offs on a single MPI rank, set :term:`nstlist` to zero and :term:`ns-type` =simple.

        xy
            Use periodic boundary conditions in x and y directions only. This works only with :term:`ns-type` =grid and can be used in combination with walls_. Without walls or with only one wall the system size is infinite in the z direction. Therefore pressure coupling or Ewald summation methods can not be used. These disadvantages do not apply when two walls are used.

    periodic-molecules
        no
            molecules are finite, fast molecular PBC can be used

        yes
            for systems with molecules that couple to themselves through the periodic boundary conditions, this requires a slower PBC algorithm and molecules are not made whole in the output

    verlet-buffer-tolerance
        (0.005) \[kJ/mol/ps\]
        Useful only with the :term:`Verlet` :term:`cutoff-scheme`. This sets the maximum allowed error for pair interactions per particle caused by the Verlet buffer, which indirectly sets :term:`rlist`. As both :term:`nstlist` and the Verlet buffer size are fixed (for performance reasons), particle pairs not in the pair list can occasionally get within the cut-off distance during :term:`nstlist` -1 steps. This causes very small jumps in the energy. In a constant-temperature ensemble, these very small energy jumps can be estimated for a given cut-off and :term:`rlist`. The estimate assumes a homogeneous particle distribution, hence the errors might be slightly underestimated for multi-phase systems. For longer pair-list life-time (:term:`nstlist` -1) * :term:`dt` the buffer is overestimated, because the interactions between particles are ignored. Combined with cancellation of errors, the actual drift of the total energy is usually one to two orders of magnitude smaller. Note that the generated buffer size takes into account that the GROMACS pair-list setup leads to a reduction in the drift by a factor 10, compared to a simple particle-pair based list. Without dynamics (energy minimization etc.), the buffer is 5% of the cut-off. For NVE simulations the initial temperature is used, unless this is zero, in which case a buffer of 10% is used. For NVE simulations the tolerance usually needs to be lowered to achieve proper energy conservation on the nanosecond time scale. To override the automated buffer setting, use :term:`verlet-buffer-tolerance` =-1 and set :term:`rlist` manually.

    rlist
        (1) \[nm\]
        Cut-off distance for the short-range neighbor list. With the :term:`Verlet` :term:`cutoff-scheme`, this is by default set by the :term:`verlet-buffer-tolerance` option and the value of :term:`rlist` is ignored.

    rlistlong
        (-1) \[nm\]
        Cut-off distance for the long-range neighbor list. This parameter is only relevant for a twin-range cut-off setup with switched potentials. In that case a buffer region is required to account for the size of charge groups. In all other cases this parameter is automatically set to the longest cut-off distance.


Electrostatics
^^^^^^^^^^^^^^

.. glossary::

    coulombtype

        .. glossary::

            Cut-off
                Twin range cut-offs with neighborlist cut-off :term:`rlist` and Coulomb cut-off :term:`rcoulomb`, where :term:`rcoulomb` >= :term:`rlist`.

            Ewald
                Classical Ewald sum electrostatics. The real-space cut-off :term:`rcoulomb` should be equal to :term:`rlist`. Use *e.g.* :term:`rlist` =0.9, :term:`rcoulomb` =0.9. The highest magnitude of wave vectors used in reciprocal space is controlled by :term:`fourierspacing`. The relative accuracy of direct/reciprocal space is controlled by :term:`ewald-rtol`.
                NOTE: Ewald scales as O(N^3/2) and is thus extremely slow for large systems. It is included mainly for reference - in most cases PME will perform much better.

            PME
                Fast smooth Particle-Mesh Ewald (SPME) electrostatics. Direct space is similar to the Ewald sum, while the reciprocal part is performed with FFTs. Grid dimensions are controlled with :term:`fourierspacing` and the interpolation order with :term:`pme-order`. With a grid spacing of 0.1 nm and cubic interpolation the electrostatic forces have an accuracy of 2-3*10^-4. Since the error from the vdw-cutoff is larger than this you might try 0.15 nm. When running in parallel the interpolation parallelizes better than the FFT, so try decreasing grid dimensions while increasing interpolation.

            P3M-AD
                Particle-Particle Particle-Mesh algorithm with analytical derivative for for long range electrostatic interactions. The method and code is identical to SPME, except that the influence function is optimized for the grid. This gives a slight increase in accuracy.

            Reaction-Field
                Reaction field electrostatics with Coulomb cut-off :term:`rcoulomb`, where :term:`rcoulomb` >= :term:`rlist`. The dielectric constant beyond the cut-off is :term:`epsilon-rf`. The dielectric constant can be set to infinity by setting :term:`epsilon-rf` =0.

            Generalized-Reaction-Field
                Generalized reaction field with Coulomb cut-off :term:`rcoulomb`, where :term:`rcoulomb` >= :term:`rlist`. The dielectric constant beyond the cut-off is :term:`epsilon-rf`. The ionic strength is computed from the number of charged (*i.e.* with non zero charge) charge groups. The temperature for the GRF potential is set with :term:`ref-t`.

            Reaction-Field-zero
                In GROMACS, normal reaction-field electrostatics with :term:`cutoff-scheme` = :term:`group` leads to bad energy conservation. :term:`Reaction-Field-zero` solves this by making the potential zero beyond the cut-off. It can only be used with an infinite dielectric constant (:term:`epsilon-rf` =0), because only for that value the force vanishes at the cut-off. :term:`rlist` should be 0.1 to 0.3 nm larger than :term:`rcoulomb` to accommodate for the size of charge groups and diffusion between neighbor list updates. This, and the fact that table lookups are used instead of analytical functions make :term:`Reaction-Field-zero` computationally more expensive than normal reaction-field.

            Reaction-Field-nec
                The same as :term:`Reaction-Field`, but implemented as in GROMACS versions before 3.3. No reaction-field correction is applied to excluded atom pairs and self pairs. The 1-4 interactions are calculated using a reaction-field. The missing correction due to the excluded pairs that do not have a 1-4 interaction is up to a few percent of the total electrostatic energy and causes a minor difference in the forces and the pressure.

            Shift
                Analogous to Shift for :term:`vdwtype`. You might want to use :term:`Reaction-Field-zero` instead, which has a similar potential shape, but has a physical interpretation and has better energies due to the exclusion correction terms.

            Encad-Shift
                The Coulomb potential is decreased over the whole range, using the definition from the Encad simulation package.

            Switch
                Analogous to Switch for :term:`vdwtype`. Switching the Coulomb potential can lead to serious artifacts, advice: use :term:`Reaction-Field-zero` instead.

            User
                mdrun_ will now expect to find a file ``table.xvg`` with user-defined potential functions for repulsion, dispersion and Coulomb. When pair interactions are present, mdrun_ also expects to find a file ``tablep.xvg`` for the pair interactions. When the same interactions should be used for non-bonded and pair interactions the user can specify the same file name for both table files. These files should contain 7 columns: the ``x`` value, ``f(x)``, ``-f'(x)``, ``g(x)``, ``-g'(x)``, ``h(x)``, ``-h'(x)``, where ``f(x)`` is the Coulomb function, ``g(x)`` the dispersion function and ``h(x)`` the repulsion function. When :term:`vdwtype` is not set to User the values for ``g``, ``-g'``, ``h`` and ``-h'`` are ignored. For the non-bonded interactions ``x`` values should run from 0 to the largest cut-off distance + :term:`table-extension` and should be uniformly spaced. For the pair interactions the table length in the file will be used. The optimal spacing, which is used for non-user tables, is ``0.002 nm`` when you run in mixed precision or ``0.0005 nm`` when you run in double precision. The function value at ``x=0`` is not important. More information is in the printed manual.

            PME-Switch
                A combination of PME and a switch function for the direct-space part (see above). :term:`rcoulomb` is allowed to be smaller than :term:`rlist`. This is mainly useful constant energy simulations (note that using PME with :term:`cutoff-scheme` = :term:`Verlet` will be more efficient).

            PME-User
                A combination of PME and user tables (see above). :term:`rcoulomb` is allowed to be smaller than :term:`rlist`. The PME mesh contribution is subtracted from the user table by mdrun_. Because of this subtraction the user tables should contain about 10 decimal places.

            PME-User-Switch
                A combination of PME-User and a switching function (see above). The switching function is applied to final particle-particle interaction, *i.e.* both to the user supplied function and the PME Mesh correction part.

    coulomb-modifier
        Potential-shift-Verlet
            Selects Potential-shift with the Verlet cutoff-scheme, as it is (nearly) free; selects None with the group cutoff-scheme.

        Potential-shift
            Shift the Coulomb potential by a constant such that it is zero at the cut-off. This makes the potential the integral of the force. Note that this does not affect the forces or the sampling.

        None
            Use an unmodified Coulomb potential. With the group scheme this means no exact cut-off is used, energies and forces are calculated for all pairs in the neighborlist.

    rcoulomb-switch
        (0) \[nm\]
        where to start switching the Coulomb potential, only relevant when force or potential switching is used

    rcoulomb
        (1) \[nm\]
        distance for the Coulomb cut-off

    epsilon-r
        (1)
        The relative dielectric constant. A value of 0 means infinity.

    epsilon-rf
        (0)
        The relative dielectric constant of the reaction field. This is only used with reaction-field electrostatics. A value of 0 means infinity.


Van der Waals
^^^^^^^^^^^^^

.. glossary::

    vdwtype
        Cut-off
            Twin range cut-offs with neighbor list cut-off :term:`rlist` and VdW cut-off :term:`rvdw`, where :term:`rvdw` >= :term:`rlist`.

        PME
            Fast smooth Particle-mesh Ewald (SPME) for VdW interactions. The grid dimensions are controlled with :term:`fourierspacing` in the same way as for electrostatics, and the interpolation order is controlled with :term:`pme-order`. The relative accuracy of direct/reciprocal space is controlled by :term:`ewald-rtol-lj`, and the specific combination rules that are to be used by the reciprocal routine are set using :term:`lj-pme-comb-rule`.

        Shift
            This functionality is deprecated and replaced by :term:`vdw-modifier` = Force-switch. The LJ (not Buckingham) potential is decreased over the whole range and the forces decay smoothly to zero between :term:`rvdw-switch` and :term:`rvdw`. The neighbor search cut-off :term:`rlist` should be 0.1 to 0.3 nm larger than :term:`rvdw` to accommodate for the size of charge groups and diffusion between neighbor list updates.

        Switch
            This functionality is deprecated and replaced by :term:`vdw-modifier` = Potential-switch. The LJ (not Buckingham) potential is normal out to :term:`rvdw-switch`, after which it is switched off to reach zero at :term:`rvdw`. Both the potential and force functions are continuously smooth, but be aware that all switch functions will give rise to a bulge (increase) in the force (since we are switching the potential). The neighbor search cut-off :term:`rlist` should be 0.1 to 0.3 nm larger than :term:`rvdw` to accommodate for the size of charge groups and diffusion between neighbor list updates.

        Encad-Shift
            The LJ (not Buckingham) potential is decreased over the whole range, using the definition from the Encad simulation package.

        User
            See user for :term:`coulombtype`. The function value at zero is not important. When you want to use LJ correction, make sure that :term:`rvdw` corresponds to the cut-off in the user-defined function. When :term:`coulombtype` is not set to User the values for the ``f`` and ``-f'`` columns are ignored.

    vdw-modifier
        Potential-shift-Verlet
            Selects Potential-shift with the Verlet cutoff-scheme, as it is (nearly) free; selects None with the group cutoff-scheme.

        Potential-shift
            Shift the Van der Waals potential by a constant such that it is zero at the cut-off. This makes the potential the integral of the force. Note that this does not affect the forces or the sampling.

        None
            Use an unmodified Van der Waals potential. With the group scheme this means no exact cut-off is used, energies and forces are calculated for all pairs in the neighborlist.

        Force-switch
            Smoothly switches the forces to zero between :term:`rvdw-switch` and :term:`rvdw`. This shifts the potential shift over the whole range and switches it to zero at the cut-off. Note that this is more expensive to calculate than a plain cut-off and it is not required for energy conservation, since Potential-shift conserves energy just as well.

        Potential-switch
            Smoothly switches the potential to zero between :term:`rvdw-switch` and :term:`rvdw`. Note that this introduces articifically large forces in the switching region and is much more expensive to calculate. This option should only be used if the force field you are using requires this.

    rvdw-switch
        (0) \[nm\]
        where to start switching the LJ force and possibly the potential, only relevant when force or potential switching is used

    rvdw
        (1) \[nm\]
        distance for the LJ or Buckingham cut-off

    DispCorr
        no
            don't apply any correction

        EnerPres
            apply long range dispersion corrections for Energy and Pressure

        Ener
            apply long range dispersion corrections for Energy only


Tables
^^^^^^

.. glossary::

    table-extension
        (1) \[nm\]
        Extension of the non-bonded potential lookup tables beyond the largest cut-off distance. The value should be large enough to account for charge group sizes and the diffusion between neighbor-list updates. Without user defined potential the same table length is used for the lookup tables for the 1-4 interactions, which are always tabulated irrespective of the use of tables for the non-bonded interactions. The value of :term:`table-extension` in no way affects the values of :term:`rlist`, :term:`rcoulomb`, or :term:`rvdw`.

    energygrp-table
        When user tables are used for electrostatics and/or VdW, here one can give pairs of energy groups for which seperate user tables should be used. The two energy groups will be appended to the table file name, in order of their definition in :term:`energygrps`, seperated by underscores. For example, if ``energygrps = Na Cl Sol`` and ``energygrp-table = Na Na Na Cl``, mdrun_ will read ``table_Na_Na.xvg`` and ``table_Na_Cl.xvg`` in addition to the normal ``table.xvg`` which will be used for all other energy group pairs.


Ewald
^^^^^

.. glossary::

    fourierspacing
        (0.12) \[nm\]
        For ordinary Ewald, the ratio of the box dimensions and the spacing determines a lower bound for the number of wave vectors to use in each (signed) direction. For PME and P3M, that ratio determines a lower bound for the number of Fourier-space grid points that will be used along that axis. In all cases, the number for each direction can be overridden by entering a non-zero value for that :term:`fourier-nx` direction. For optimizing the relative load of the particle-particle interactions and the mesh part of PME, it is useful to know that the accuracy of the electrostatics remains nearly constant when the Coulomb cut-off and the PME grid spacing are scaled by the same factor.

    fourier-nx
    fourier-ny
    fourier-nz
        (0)
        Highest magnitude of wave vectors in reciprocal space when using Ewald.
        Grid size when using PME or P3M. These values override :term:`fourierspacing` per direction. The best choice is powers of 2, 3, 5 and 7. Avoid large primes.

    pme-order
        (4)
        Interpolation order for PME. 4 equals cubic interpolation. You might try 6/8/10 when running in parallel and simultaneously decrease grid dimension.

    ewald-rtol
        (1e-5)
        The relative strength of the Ewald-shifted direct potential at :term:`rcoulomb` is given by :term:`ewald-rtol`. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

    ewald-rtol-lj
        (1e-3)
        When doing PME for VdW-interactions, :term:`ewald-rtol-lj` is used to control the relative strength of the dispersion potential at :term:`rvdw` in the same way as :term:`ewald-rtol` controls the electrostatic potential.

    lj-pme-comb-rule
        (Geometric)
        The combination rules used to combine VdW-parameters in the reciprocal part of LJ-PME. Geometric rules are much faster than Lorentz-Berthelot and usually the recommended choice, even when the rest of the force field uses the Lorentz-Berthelot rules.

        Geometric
            Apply geometric combination rules

        Lorentz-Berthelot
            Apply Lorentz-Berthelot combination rules

    ewald-geometry
        3d
            The Ewald sum is performed in all three dimensions.

        3dc
            The reciprocal sum is still performed in 3D, but a force and potential correction applied in the `z` dimension to produce a pseudo-2D summation. If your system has a slab geometry in the `x-y` plane you can try to increase the `z`-dimension of the box (a box height of 3 times the slab height is usually ok) and use this option.

    epsilon-surface
        (0)
        This controls the dipole correction to the Ewald summation in 3D. The default value of zero means it is turned off. Turn it on by setting it to the value of the relative permittivity of the imaginary surface around your infinite system. Be careful - you shouldn't use this if you have free mobile charges in your system. This value does not affect the slab 3DC variant of the long range corrections.


Temperature coupling
^^^^^^^^^^^^^^^^^^^^

.. glossary::

    tcoupl
        no
            No temperature coupling.

        berendsen
            Temperature coupling with a Berendsen-thermostat to a bath with temperature :term:`ref-t`, with time constant :term:`tau-t`. Several groups can be coupled separately, these are specified in the :term:`tc-grps` field separated by spaces.

        nose-hoover
            Temperature coupling using a Nose-Hoover extended ensemble. The reference temperature and coupling groups are selected as above, but in this case :term:`tau-t` controls the period of the temperature fluctuations at equilibrium, which is slightly different from a relaxation time. For NVT simulations the conserved energy quantity is written to energy and log file.

        andersen
            Temperature coupling by randomizing a fraction of the particles at each timestep. Reference temperature and coupling groups are selected as above. :term:`tau-t` is the average time between randomization of each molecule. Inhibits particle dynamics somewhat, but little or no ergodicity issues. Currently only implemented with velocity Verlet, and not implemented with constraints.

        andersen-massive
            Temperature coupling by randomizing all particles at infrequent timesteps. Reference temperature and coupling groups are selected as above. :term:`tau-t` is the time between randomization of all molecules. Inhibits particle dynamics somewhat, but little or no ergodicity issues. Currently only implemented with velocity Verlet.

        v-rescale
            Temperature coupling using velocity rescaling with a stochastic term (JCP 126, 014101). This thermostat is similar to Berendsen coupling, with the same scaling using :term:`tau-t`, but the stochastic term ensures that a proper canonical ensemble is generated. The random seed is set with :term:`ld-seed`. This thermostat works correctly even for :term:`tau-t` =0. For NVT simulations the conserved energy quantity is written to the energy and log file.

    nsttcouple
        (-1)
        The frequency for coupling the temperature. The default value of -1 sets :term:`nsttcouple` equal to :term:`nstlist`, unless :term:`nstlist` <=0, then a value of 10 is used. For velocity Verlet integrators :term:`nsttcouple` is set to 1.

    nh-chain-length
        (10)
        The number of chained Nose-Hoover thermostats for velocity Verlet integrators, the leap-frog :term:`md` integrator only supports 1. Data for the NH chain variables is not printed to the edr_ file, but can be using the ``GMX_NOSEHOOVER_CHAINS`` environment variable

    tc-grps
        groups to couple to separate temperature baths

    tau-t
        \[ps\]
        time constant for coupling (one for each group in :term:`tc-grps`), -1 means no temperature coupling

    ref-t
        \[K\]
        reference temperature for coupling (one for each group in :term:`tc-grps`)


Pressure coupling
^^^^^^^^^^^^^^^^^

.. glossary::

    pcoupl
        no
            No pressure coupling. This means a fixed box size.

        berendsen
            Exponential relaxation pressure coupling with time constant :term:`tau-p`. The box is scaled every timestep. It has been argued that this does not yield a correct thermodynamic ensemble, but it is the most efficient way to scale a box at the beginning of a run.

        Parrinello-Rahman
            Extended-ensemble pressure coupling where the box vectors are subject to an equation of motion. The equation of motion for the atoms is coupled to this. No instantaneous scaling takes place. As for Nose-Hoover temperature coupling the time constant :term:`tau-p` is the period of pressure fluctuations at equilibrium. This is probably a better method when you want to apply pressure scaling during data collection, but beware that you can get very large oscillations if you are starting from a different pressure. For simulations where the exact fluctation of the NPT ensemble are important, or if the pressure coupling time is very short it may not be appropriate, as the previous time step pressure is used in some steps of the GROMACS implementation for the current time step pressure.

        MTTK
            Martyna-Tuckerman-Tobias-Klein implementation, only useable with :term:`md-vv` or :term:`md-vv-avek`, very similar to Parrinello-Rahman. As for Nose-Hoover temperature coupling the time constant :term:`tau-p` is the period of pressure fluctuations at equilibrium. This is probably a better method when you want to apply pressure scaling during data collection, but beware that you can get very large oscillations if you are starting from a different pressure. Currently (as of version 5.1), it only supports isotropic scaling, and only works without constraints.

    pcoupltype
        isotropic
            Isotropic pressure coupling with time constant :term:`tau-p`. The compressibility and reference pressure are set with :term:`compressibility` and :term:`ref-p`, one value is needed.

        semiisotropic
            Pressure coupling which is isotropic in the ``x`` and ``y`` direction, but different in the ``z`` direction. This can be useful for membrane simulations. 2 values are needed for ``x/y`` and ``z`` directions respectively.

        anisotropic
            Idem, but 6 values are needed for ``xx``, ``yy``, ``zz``, ``xy/yx``, ``xz/zx`` and ``yz/zy`` components, respectively. When the off-diagonal compressibilities are set to zero, a rectangular box will stay rectangular. Beware that anisotropic scaling can lead to extreme deformation of the simulation box.

        surface-tension
            Surface tension coupling for surfaces parallel to the xy-plane. Uses normal pressure coupling for the `z`-direction, while the surface tension is coupled to the `x/y` dimensions of the box. The first :term:`ref-p` value is the reference surface tension times the number of surfaces ``bar nm``, the second value is the reference `z`-pressure ``bar``. The two :term:`compressibility` values are the compressibility in the `x/y` and `z` direction respectively. The value for the `z`-compressibility should be reasonably accurate since it influences the convergence of the surface-tension, it can also be set to zero to have a box with constant height.

    nstpcouple
        (-1)
        The frequency for coupling the pressure. The default value of -1 sets :term:`nstpcouple` equal to :term:`nstlist`, unless :term:`nstlist` <=0, then a value of 10 is used. For velocity Verlet integrators :term:`nstpcouple` is set to 1.

    tau-p
        (1) \[ps\]
        time constant for coupling

    compressibility
        \[bar^-1\]
        compressibility (NOTE: this is now really in bar-1) For water at 1 atm and 300 K the compressibility is 4.5e-5 bar^-1.

    ref-p
        \[bar\]
        reference pressure for coupling

    refcoord-scaling
        no
            The reference coordinates for position restraints are not modified. Note that with this option the virial and pressure will depend on the absolute positions of the reference coordinates.

        all
            The reference coordinates are scaled with the scaling matrix of the pressure coupling.

        com
            Scale the center of mass of the reference coordinates with the scaling matrix of the pressure coupling. The vectors of each reference coordinate to the center of mass are not scaled. Only one COM is used, even when there are multiple molecules with position restraints. For calculating the COM of the reference coordinates in the starting configuration, periodic boundary conditions are not taken into account.


Simulated annealing
^^^^^^^^^^^^^^^^^^^

Simulated annealing is controlled separately for each temperature group in GROMACS. The reference temperature is a piecewise linear function, but you can use an arbitrary number of points for each group, and choose either a single sequence or a periodic behaviour for each group. The actual annealing is performed by dynamically changing the reference temperature used in the thermostat algorithm selected, so remember that the system will usually not instantaneously reach the reference temperature!

.. glossary::

    annealing
        Type of annealing for each temperature group

        no
             No simulated annealing - just couple to reference temperature value.

        single
             A single sequence of annealing points. If your simulation is longer than the time of the last point, the temperature will be coupled to this constant value after the annealing sequence has reached the last time point.

        periodic
             The annealing will start over at the first reference point once the last reference time is reached. This is repeated until the simulation ends.

    annealing-npoints
         A list with the number of annealing reference/control points used for each temperature group. Use 0 for groups that are not annealed. The number of entries should equal the number of temperature groups.

    annealing-time
        List of times at the annealing reference/control points for each group. If you are using periodic annealing, the times will be used modulo the last value, *i.e.* if the values are 0, 5, 10, and 15, the coupling will restart at the 0ps value after 15ps, 30ps, 45ps, etc. The number of entries should equal the sum of the numbers given in :term:`annealing-npoints`.

    annealing-temp
        List of temperatures at the annealing reference/control points for each group. The number of entries should equal the sum of the numbers given in :term:`annealing-npoints`.

Confused? OK, let's use an example. Assume you have two temperature groups, set the group selections to ``annealing = single periodic``, the number of points of each group to ``annealing-npoints = 3 4``, the times to ``annealing-time = 0 3 6 0 2 4 6`` and finally temperatures to ``annealing-temp = 298 280 270 298 320 320 298``. The first group will be coupled to 298K at 0ps, but the reference temperature will drop linearly to reach 280K at 3ps, and then linearly between 280K and 270K from 3ps to 6ps. After this is stays constant, at 270K. The second group is coupled to 298K at 0ps, it increases linearly to 320K at 2ps, where it stays constant until 4ps. Between 4ps and 6ps it decreases to 298K, and then it starts over with the same pattern again, *i.e.* rising linearly from 298K to 320K between 6ps and 8ps. Check the summary printed by `gmx grompp`_ if you are unsure!


Velocity generation
^^^^^^^^^^^^^^^^^^^

.. glossary::

    gen-vel
        no
        Do not generate velocities. The velocities are set to zero when there are no velocities in the input structure file.

        yes
        Generate velocities in `gmx grompp`_ according to a Maxwell distribution at temperature :term:`gen-temp`, with random seed :term:`gen-seed`. This is only meaningful with integrator :term:`md`.

    gen-temp
        (300) \[K\]
        temperature for Maxwell distribution

    gen-seed
        (-1) \[integer\]
        used to initialize random generator for random velocities, when :term:`gen-seed` is set to -1, a pseudo random seed is used.


Bonds
^^^^^

.. glossary::

    constraints
        none
            No constraints except for those defined explicitly in the topology, *i.e.* bonds are represented by a harmonic (or other) potential or a Morse potential (depending on the setting of :term:`morse`) and angles by a harmonic (or other) potential.

        h-bonds
            Convert the bonds with H-atoms to constraints.

        all-bonds
            Convert all bonds to constraints.

        h-angles
            Convert all bonds and additionally the angles that involve H-atoms to bond-constraints.

        all-angles
            Convert all bonds and angles to bond-constraints.

    constraint-algorithm
        LINCS
            LINear Constraint Solver. With domain decomposition the parallel version P-LINCS is used. The accuracy in set with :term:`lincs-order`, which sets the number of matrices in the expansion for the matrix inversion. After the matrix inversion correction the algorithm does an iterative correction to compensate for lengthening due to rotation. The number of such iterations can be controlled with :term:`lincs-iter`. The root mean square relative constraint deviation is printed to the log file every :term:`nstlog` steps. If a bond rotates more than :term:`lincs-warnangle` in one step, a warning will be printed both to the log file and to ``stderr``. LINCS should not be used with coupled angle constraints.

        SHAKE
            SHAKE is slightly slower and less stable than LINCS, but does work with angle constraints. The relative tolerance is set with :term:`shake-tol`, 0.0001 is a good value for "normal" MD. SHAKE does not support constraints between atoms on different nodes, thus it can not be used with domain decompositon when inter charge-group constraints are present. SHAKE can not be used with energy minimization.

    continuation
        This option was formerly known as unconstrained-start.

        no
            apply constraints to the start configuration and reset shells

        yes
            do not apply constraints to the start configuration and do not reset shells, useful for exact coninuation and reruns

    shake-tol
        (0.0001)
        relative tolerance for SHAKE

    lincs-order
        (4)
        Highest order in the expansion of the constraint coupling matrix. When constraints form triangles, an additional expansion of the same order is applied on top of the normal expansion only for the couplings within such triangles. For "normal" MD simulations an order of 4 usually suffices, 6 is needed for large time-steps with virtual sites or BD. For accurate energy minimization an order of 8 or more might be required. With domain decomposition, the cell size is limited by the distance spanned by :term:`lincs-order` +1 constraints. When one wants to scale further than this limit, one can decrease :term:`lincs-order` and increase :term:`lincs-iter`, since the accuracy does not deteriorate when (1+ :term:`lincs-iter` )* :term:`lincs-order` remains constant.

    lincs-iter
        (1)
        Number of iterations to correct for rotational lengthening in LINCS. For normal runs a single step is sufficient, but for NVE runs where you want to conserve energy accurately or for accurate energy minimization you might want to increase it to 2.

    lincs-warnangle
        (30) \[degrees\]
        maximum angle that a bond can rotate before LINCS will complain

    morse
        no
            bonds are represented by a harmonic potential

        yes
            bonds are represented by a Morse potential


Energy group exclusions
^^^^^^^^^^^^^^^^^^^^^^^

.. glossary::

    energygrp-excl:
        Pairs of energy groups for which all non-bonded interactions are excluded. An example: if you have two energy groups ``Protein`` and ``SOL``, specifying ``energygrp-excl = Protein Protein  SOL SOL`` would give only the non-bonded interactions between the protein and the solvent. This is especially useful for speeding up energy calculations with ``mdrun -rerun`` and for excluding interactions within frozen groups.


Walls
^^^^^

.. glossary::

    nwall
        (0)
        When set to 1 there is a wall at ``z=0``, when set to 2 there is also a wall at ``z=z-box``. Walls can only be used with :term:`pbc` ``=xy``. When set to 2 pressure coupling and Ewald summation can be used (it is usually best to use semiisotropic pressure coupling with the ``x/y`` compressibility set to 0, as otherwise the surface area will change). Walls interact wit the rest of the system through an optional :term:`wall-atomtype`. Energy groups ``wall0`` and ``wall1`` (for :term:`nwall` =2) are added automatically to monitor the interaction of energy groups with each wall. The center of mass motion removal will be turned off in the ``z``-direction.

    wall-atomtype
        the atom type name in the force field for each wall. By (for example) defining a special wall atom type in the topology with its own combination rules, this allows for independent tuning of the interaction of each atomtype with the walls.

    wall-type
        9-3
            LJ integrated over the volume behind the wall: 9-3 potential

        10-4
            LJ integrated over the wall surface: 10-4 potential

        12-6
            direct LJ potential with the ``z`` distance from the wall

    table
        user defined potentials indexed with the ``z`` distance from the wall, the tables are read analogously to the :term:`energygrp-table` option, where the first name is for a "normal" energy group and the second name is ``wall0`` or ``wall1``, only the dispersion and repulsion columns are used

    wall-r-linpot
        (-1) \[nm\]
        Below this distance from the wall the potential is continued linearly and thus the force is constant. Setting this option to a postive value is especially useful for equilibration when some atoms are beyond a wall. When the value is <=0 (<0 for :term:`wall-type` =table), a fatal error is generated when atoms are beyond a wall.

    wall-density
        \[nm^-3/nm^-2\]
        the number density of the atoms for each wall for wall types 9-3 and 10-4

    wall-ewald-zfac
        (3)
        The scaling factor for the third box vector for Ewald summation only, the minimum is 2. Ewald summation can only be used with :term:`nwall` =2, where one should use :term:`ewald-geometry` ``=3dc``. The empty layer in the box serves to decrease the unphysical Coulomb interaction between periodic images.


COM pulling
^^^^^^^^^^^
Note that where pulling coordinate are applicable, there can be more than one (set with :term:`pull-ncoords`) and multiple related mdp_ variables will exist accordingly. Documentation references to things like :term:`pull-coord1-vec` should be understood to apply to to the applicable pulling coordinate.

.. glossary::

    pull
        no
            No center of mass pulling. All the following pull options will be ignored (and if present in the mdp_ file, they unfortunately generate warnings)

        umbrella
            Center of mass pulling using an umbrella potential between the reference group and one or more groups.

        constraint
            Center of mass pulling using a constraint between the reference group and one or more groups. The setup is identical to the option umbrella, except for the fact that a rigid constraint is applied instead of a harmonic potential.

        constant-force
            Center of mass pulling using a linear potential and therefore a constant force. For this option there is no reference position and therefore the parameters :term:`pull-coord1-init` and :term:`pull-coord1-rate` are not used.

    pull-geometry
        distance
            Pull along the vector connecting the two groups. Components can be selected with :term:`pull-dim`.

        direction
            Pull in the direction of :term:`pull-coord1-vec`.

        direction-periodic
            As direction, but allows the distance to be larger than half the box size. With this geometry the box should not be dynamic (*e.g.* no pressure scaling) in the pull dimensions and the pull force is not added to virial.

        cylinder
            Designed for pulling with respect to a layer where the reference COM is given by a local cylindrical part of the reference group. The pulling is in the direction of :term:`pull-coord1-vec`. From the reference group a cylinder is selected around the axis going through the pull group with direction :term:`pull-coord1-vec` using two radii. The radius :term:`pull-r1` gives the radius within which all the relative weights are one, between :term:`pull-r1` and :term:`pull-r0` the weights are switched to zero. Mass weighting is also used. Note that the radii should be smaller than half the box size. For tilted cylinders they should be even smaller than half the box size since the distance of an atom in the reference group from the COM of the pull group has both a radial and an axial component.

    pull-dim
        (Y Y Y)
        the distance components to be used with :term:`pull-geometry` distance, and also sets which components are printed to the output files

    pull-r1
        (1) \[nm\]
        the inner radius of the cylinder for :term:`pull-geometry` cylinder

    pull-r0
        (1) \[nm\]
        the outer radius of the cylinder for :term:`pull-geometry` cylinder

    pull-constr-tol
        (1e-6)
        the relative constraint tolerance for constraint pulling

    pull-start
        no
            do not modify :term:`pull-coord1-init`

        yes
            add the COM distance of the starting conformation to :term:`pull-coord1-init`

    pull-print-reference
        no
            do not print the COM of the first group in each pull coordinate

        yes
            print the COM of the first group in each pull coordinate

    pull-nstxout
        (10)
        frequency for writing out the COMs of all the pull group

    pull-nstfout
        (1)
        frequency for writing out the force of all the pulled group

    pull-ngroups
        (1)
        The number of pull groups, not including the absolute reference group, when used. Pull groups can be reused in multiple pull coordinates. Below only the pull options for group 1 are given, further groups simply increase the group index number.

    pull-ncoords
        (1)
        The number of pull coordinates. Below only the pull options for coordinate 1 are given, further coordinates simply increase the coordinate index number.

    pull-group1-name
        The name of the pull group, is looked up in the index file or in the default groups to obtain the atoms involved.

    pull-group1-weights
        Optional relative weights which are multiplied with the masses of the atoms to give the total weight for the COM. The number should be 0, meaning all 1, or the number of atoms in the pull group.

    pull-group1-pbcatom
        (0)
        The reference atom for the treatment of periodic boundary conditions inside the group (this has no effect on the treatment of the pbc between groups). This option is only important when the diameter of the pull group is larger than half the shortest box vector. For determining the COM, all atoms in the group are put at their periodic image which is closest to :term:`pull-group1-pbcatom`. A value of 0 means that the middle atom (number wise) is used. This parameter is not used with :term:`pull-geometry` cylinder. A value of -1 turns on cosine weighting, which is useful for a group of molecules in a periodic system, *e.g.* a water slab (see Engin et al. J. Chem. Phys. B 2010).

    pull-coord1-groups
        The two groups indices should be given on which this pull coordinate will operate. The first index can be 0, in which case an absolute reference of :term:`pull-coord1-origin` is used. With an absolute reference the system is no longer translation invariant and one should think about what to do with the center of mass motion.

    pull-coord1-origin
        (0.0 0.0 0.0)
        The pull reference position for use with an absolute reference.

    pull-coord1-vec
        (0.0 0.0 0.0)
        The pull direction. `gmx grompp`_ normalizes the vector.

    pull-coord1-init
        (0.0) \[nm\]
        The reference distance at t=0.

    pull-coord1-rate
        (0) \[nm/ps\]
        The rate of change of the reference position.

    pull-coord1-k
        (0) \[kJ mol-1 nm-2\] / \[kJ mol-1 nm-1\]
        The force constant. For umbrella pulling this is the harmonic force constant in kJ mol-1 nm-2. For constant force pulling this is the force constant of the linear potential, and thus the negative (!) of the constant force in kJ mol-1 nm-1.

    pull-coord1-kB
        (pull-k1) \[kJ mol-1 nm-2\] / \[kJ mol-1 nm-1\]
        As :term:`pull-coord1-k`, but for state B. This is only used when :term:`free-energy` is turned on. The force constant is then (1 - lambda) * :term:`pull-coord1-k` + lambda * :term:`pull-coord1-kB`.


NMR refinement
^^^^^^^^^^^^^^

.. glossary::

    disre
        no
            ignore distance restraint information in topology file

        simple
            simple (per-molecule) distance restraints.

        ensemble
            distance restraints over an ensemble of molecules in one simulation box. Normally, one would perform ensemble averaging over multiple subsystems, each in a separate box, using ``mdrun -multi``. Supply ``topol0.tpr``, ``topol1.tpr`, ... with different coordinates and/or velocities. The environment variable ``GMX_DISRE_ENSEMBLE_SIZE`` sets the number of systems within each ensemble (usually equal to the ``mdrun -multi`` value).

    disre-weighting
        equal
            divide the restraint force equally over all atom pairs in the restraint

        conservative
            the forces are the derivative of the restraint potential, this results in an weighting of the atom pairs to the reciprocal seventh power of the displacement. The forces are conservative when :term:`disre-tau` is zero.

    disre-mixed
        no
            the violation used in the calculation of the restraint force is the time-averaged violation

        yes
            the violation used in the calculation of the restraint force is the square root of the product of the time-averaged violation and the instantaneous violation

    disre-fc
        (1000) \[kJ mol-1 nm-2\]
        force constant for distance restraints, which is multiplied by a (possibly) different factor for each restraint given in the `fac` column of the interaction in the topology file.

    disre-tau
        (0) \[ps\]
        time constant for distance restraints running average. A value of zero turns off time averaging.

    nstdisreout
        (100) \[steps\]
        period between steps when the running time-averaged and instantaneous distances of all atom pairs involved in restraints are written to the energy file (can make the energy file very large)

    orire
        no
            ignore orientation restraint information in topology file

        yes
            use orientation restraints, ensemble averaging can be performed with `mdrun -multi`

    orire-fc
        (0) \[kJ mol\]
        force constant for orientation restraints, which is multiplied by a (possibly) different weight factor for each restraint, can be set to zero to obtain the orientations from a free simulation

    orire-tau
        (0) \[ps\]
        time constant for orientation restraints running average. A value of zero turns off time averaging.

    orire-fitgrp
        fit group for orientation restraining. This group of atoms is used to determine the rotation **R** of the system with respect to the reference orientation. The reference orientation is the starting conformation of the first subsystem. For a protein, backbone is a reasonable choice

    nstorireout
        (100) \[steps\]
        period between steps when the running time-averaged and instantaneous orientations for all restraints, and the molecular order tensor are written to the energy file (can make the energy file very large)


Free energy calculations
^^^^^^^^^^^^^^^^^^^^^^^^

.. glossary::

    free-energy
        no
            Only use topology A.

        yes
            Interpolate between topology A (lambda=0) to topology B (lambda=1) and write the derivative of the Hamiltonian with respect to lambda (as specified with :term:`dhdl-derivatives`), or the Hamiltonian differences with respect to other lambda values (as specified with foreign lambda) to the energy file and/or to ``dhdl.xvg``, where they can be processed by, for example `gmx bar`_. The potentials, bond-lengths and angles are interpolated linearly as described in the manual. When :term:`sc-alpha` is larger than zero, soft-core potentials are used for the LJ and Coulomb interactions.

    expanded
        Turns on expanded ensemble simulation, where the alchemical state becomes a dynamic variable, allowing jumping between different Hamiltonians. See the expanded ensemble options for controlling how expanded ensemble simulations are performed. The different Hamiltonians used in expanded ensemble simulations are defined by the other free energy options.

    init-lambda
        (-1)
        starting value for lambda (float). Generally, this should only be used with slow growth (*i.e.* nonzero :term:`delta-lambda`). In other cases, :term:`init-lambda-state` should be specified instead. Must be greater than or equal to 0.

    delta-lambda
        (0)
        increment per time step for lambda

    init-lambda-state
        (-1)
        starting value for the lambda state (integer). Specifies which columm of the lambda vector (:term:`coul-lambdas`, :term:`vdw-lambdas`, :term:`bonded-lambdas`, :term:`restraint-lambdas`, :term:`mass-lambdas`, :term:`temperature-lambdas`, :term:`fep-lambdas`) should be used. This is a zero-based index: :term:`init-lambda-state` 0 means the first column, and so on.

    fep-lambdas
        \[array\]
        Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every :term:`nstdhdl` steps. Values must be between 0 and 1. Free energy differences between different lambda values can then be determined with `gmx bar`_. :term:`fep-lambdas` is different from the other -lambdas keywords because all components of the lambda vector that are not specified will use :term:`fep-lambdas` (including :term:`restraint-lambdas` and therefore the pull code restraints).

    coul-lambdas
        \[array\]
        Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every :term:`nstdhdl` steps. Values must be between 0 and 1. Only the electrostatic interactions are controlled with this component of the lambda vector (and only if the lambda=0 and lambda=1 states have differing electrostatic interactions).

    vdw-lambdas
        \[array\]
        Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every :term:`nstdhdl` steps. Values must be between 0 and 1. Only the van der Waals interactions are controlled with this component of the lambda vector.

    bonded-lambdas
        \[array\]
        Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every :term:`nstdhdl` steps. Values must be between 0 and 1. Only the bonded interactions are controlled with this component of the lambda vector.

    restraint-lambdas
        \[array\]
        Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every :term:`nstdhdl` steps. Values must be between 0 and 1. Only the restraint interactions: dihedral restraints, and the pull code restraints are controlled with this component of the lambda vector.

    mass-lambdas
        \[array\]
        Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every :term:`nstdhdl` steps. Values must be between 0 and 1. Only the particle masses are controlled with this component of the lambda vector.

    temperature-lambdas
        \[array\]
        Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every :term:`nstdhdl` steps. Values must be between 0 and 1. Only the temperatures controlled with this component of the lambda vector. Note that these lambdas should not be used for replica exchange, only for simulated tempering.

    calc-lambda-neighbors
        (1)
        Controls the number of lambda values for which Delta H values will be calculated and written out, if :term:`init-lambda-state` has been set. A positive value will limit the number of lambda points calculated to only the nth neighbors of :term:`init-lambda-state`: for example, if :term:`init-lambda-state` is 5 and this parameter has a value of 2, energies for lambda points 3-7 will be calculated and writen out. A value of -1 means all lambda points will be written out. For normal BAR such as with `gmx bar`_, a value of 1 is sufficient, while for MBAR -1 should be used.

    sc-alpha
        (0)
        the soft-core alpha parameter, a value of 0 results in linear interpolation of the LJ and Coulomb interactions

    sc-r-power
        (6)
        the power of the radial term in the soft-core equation. Possible values are 6 and 48. 6 is more standard, and is the default. When 48 is used, then sc-alpha should generally be much lower (between 0.001 and 0.003).

    sc-coul
        (no)
        Whether to apply the soft core free energy interaction transformation to the Columbic interaction of a molecule. Default is no, as it is generally more efficient to turn off the Coulomic interactions linearly before turning off the van der Waals interactions.

    sc-power
        (0)
        the power for lambda in the soft-core function, only the values 1 and 2 are supported

    sc-sigma
        (0.3) \[nm\]
        the soft-core sigma for particles which have a C6 or C12 parameter equal to zero or a sigma smaller than :term:`sc-sigma`

    couple-moltype
        Here one can supply a molecule type (as defined in the topology) for calculating solvation or coupling free energies. There is a special option ``system`` that couples all molecule types in the system. This can be useful for equilibrating a system starting from (nearly) random coordinates. :term:`free-energy` has to be turned on. The Van der Waals interactions and/or charges in this molecule type can be turned on or off between lambda=0 and lambda=1, depending on the settings of :term:`couple-lambda0` and :term:`couple-lambda1`. If you want to decouple one of several copies of a molecule, you need to copy and rename the molecule definition in the topology.

    couple-lambda0
        vdw-q
            all interactions are on at lambda=0

        vdw
            the charges are zero (no Coulomb interactions) at lambda=0

        q
            the Van der Waals interactions are turned at lambda=0; soft-core interactions will be required to avoid singularities

        none
            the Van der Waals interactions are turned off and the charges are zero at lambda=0; soft-core interactions will be required to avoid singularities.

    couple-lambda1
        analogous to :term:`couple-lambda1`, but for lambda=1

    couple-intramol
        no
            All intra-molecular non-bonded interactions for moleculetype :term:`couple-moltype` are replaced by exclusions and explicit pair interactions. In this manner the decoupled state of the molecule corresponds to the proper vacuum state without periodicity effects.

        yes
            The intra-molecular Van der Waals and Coulomb interactions are also turned on/off. This can be useful for partitioning free-energies of relatively large molecules, where the intra-molecular non-bonded interactions might lead to kinetically trapped vacuum conformations. The 1-4 pair interactions are not turned off.

    nstdhdl
        (100)
        the frequency for writing dH/dlambda and possibly Delta H to dhdl.xvg, 0 means no ouput, should be a multiple of :term:`nstcalcenergy`.

    dhdl-derivatives
        (yes)
        If yes (the default), the derivatives of the Hamiltonian with respect to lambda at each :term:`nstdhdl` step are written out. These values are needed for interpolation of linear energy differences with `gmx bar`_ (although the same can also be achieved with the right foreign lambda setting, that may not be as flexible), or with thermodynamic integration

    dhdl-print-energy
        (no)
        Include either the total or the potential energy in the dhdl file. Options are 'no', 'potential', or 'total'. This information is needed for later free energy analysis if the states of interest are at different temperatures. If all states are at the same temperature, this information is not needed. 'potential' is useful in case one is using ``mdrun -rerun`` to generate the ``dhdl.xvg`` file. When rerunning from an existing trajectory, the kinetic energy will often not be correct, and thus one must compute the residual free energy from the potential alone, with the kinetic energy component computed analytically.

    separate-dhdl-file
        yes
            The free energy values that are calculated (as specified with the foreign lambda and :term:`dhdl-derivatives` settings) are written out to a separate file, with the default name ``dhdl.xvg``. This file can be used directly with `gmx bar`_.

        no
            The free energy values are written out to the energy output file (``ener.edr``, in accumulated blocks at every :term:`nstenergy` steps), where they can be extracted with `gmx energy`_ or used directly with `gmx bar`_.

    dh-hist-size
        (0)
        If nonzero, specifies the size of the histogram into which the Delta H values (specified with foreign lambda) and the derivative dH/dl values are binned, and written to ener.edr. This can be used to save disk space while calculating free energy differences. One histogram gets written for each foreign lambda and two for the dH/dl, at every :term:`nstenergy` step. Be aware that incorrect histogram settings (too small size or too wide bins) can introduce errors. Do not use histograms unless you're certain you need it.

    dh-hist-spacing
        (0.1)
        Specifies the bin width of the histograms, in energy units. Used in conjunction with :term:`dh-hist-size`. This size limits the accuracy with which free energies can be calculated. Do not use histograms unless you're certain you need it.


Expanded Ensemble calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. glossary::

    nstexpanded
        The number of integration steps beween attempted moves changing the system Hamiltonian in expanded ensemble simulations. Must be a multiple of :term:`nstcalcenergy`, but can be greater or less than :term:`nstdhdl`.

    lmc-stats
        no
            No Monte Carlo in state space is performed.

        metropolis-transition
            Uses the Metropolis weights to update the expanded ensemble weight of each state. Min{1,exp(-(beta_new u_new - beta_old u_old)}

        barker-transition
            Uses the Barker transition critera to update the expanded ensemble weight of each state i, defined by exp(-beta_new u_new)/(exp(-beta_new u_new)+exp(-beta_old u_old))

        wang-landau
            Uses the Wang-Landau algorithm (in state space, not energy space) to update the expanded ensemble weights.

        min-variance
            Uses the minimum variance updating method of Escobedo et al. to update the expanded ensemble weights. Weights will not be the free energies, but will rather emphasize states that need more sampling to give even uncertainty.

    lmc-mc-move
        no
            No Monte Carlo in state space is performed.

        metropolis-transition
            Randomly chooses a new state up or down, then uses the Metropolis critera to decide whether to accept or reject: Min{1,exp(-(beta_new u_new - beta_old u_old)}

        barker-transition
            Randomly chooses a new state up or down, then uses the Barker transition critera to decide whether to accept or reject: exp(-beta_new u_new)/(exp(-beta_new u_new)+exp(-beta_old u_old))

        gibbs
             Uses the conditional weights of the state given the coordinate (exp(-beta_i u_i) / sum_k exp(beta_i u_i) to decide which state to move to.

        metropolized-gibbs
             Uses the conditional weights of the state given the coordinate (exp(-beta_i u_i) / sum_k exp(beta_i u_i) to decide which state to move to, EXCLUDING the current state, then uses a rejection step to ensure detailed balance. Always more efficient that Gibbs, though only marginally so in many situations, such as when only the nearest neighbors have decent phase space overlap.

    lmc-seed
        (-1)
        random seed to use for Monte Carlo moves in state space. When :term:`lmc-seed` is set to -1, a pseudo random seed is us

    mc-temperature
        Temperature used for acceptance/rejection for Monte Carlo moves. If not specified, the temperature of the simulation specified in the first group of :term:`ref-t` is used.

    wl-ratio
        (0.8)
        The cutoff for the histogram of state occupancies to be reset, and the free energy incrementor to be changed from delta to delta * :term:`wl-scale`. If we define the Nratio = (number of samples at each histogram) / (average number of samples at each histogram). :term:`wl-ratio` of 0.8 means that means that the histogram is only considered flat if all Nratio > 0.8 AND simultaneously all 1/Nratio > 0.8.

    wl-scale
        (0.8)
        Each time the histogram is considered flat, then the current value of the Wang-Landau incrementor for the free energies is multiplied by :term:`wl-scale`. Value must be between 0 and 1.

    init-wl-delta
        (1.0)
        The initial value of the Wang-Landau incrementor in kT. Some value near 1 kT is usually most efficient, though sometimes a value of 2-3 in units of kT works better if the free energy differences are large.

    wl-oneovert
        (no)
        Set Wang-Landau incrementor to scale with 1/(simulation time) in the large sample limit. There is significant evidence that the standard Wang-Landau algorithms in state space presented here result in free energies getting 'burned in' to incorrect values that depend on the initial state. when :term:`wl-oneovert` is true, then when the incrementor becomes less than 1/N, where N is the mumber of samples collected (and thus proportional to the data collection time, hence '1 over t'), then the Wang-Lambda incrementor is set to 1/N, decreasing every step. Once this occurs, :term:`wl-ratio` is ignored, but the weights will still stop updating when the equilibration criteria set in :term:`lmc-weights-equil` is achieved.

    lmc-repeats
        (1)
        Controls the number of times that each Monte Carlo swap type is performed each iteration. In the limit of large numbers of Monte Carlo repeats, then all methods converge to Gibbs sampling. The value will generally not need to be different from 1.

    lmc-gibbsdelta
        (-1)
        Limit Gibbs sampling to selected numbers of neighboring states. For Gibbs sampling, it is sometimes inefficient to perform Gibbs sampling over all of the states that are defined. A positive value of :term:`lmc-gibbsdelta` means that only states plus or minus :term:`lmc-gibbsdelta` are considered in exchanges up and down. A value of -1 means that all states are considered. For less than 100 states, it is probably not that expensive to include all states.

    lmc-forced-nstart
        (0)
        Force initial state space sampling to generate weights. In order to come up with reasonable initial weights, this setting allows the simulation to drive from the initial to the final lambda state, with :term:`lmc-forced-nstart` steps at each state before moving on to the next lambda state. If :term:`lmc-forced-nstart` is sufficiently long (thousands of steps, perhaps), then the weights will be close to correct. However, in most cases, it is probably better to simply run the standard weight equilibration algorithms.

    nst-transition-matrix
        (-1)
        Frequency of outputting the expanded ensemble transition matrix. A negative number means it will only be printed at the end of the simulation.

    symmetrized-transition-matrix
        (no)
        Whether to symmetrize the empirical transition matrix. In the infinite limit the matrix will be symmetric, but will diverge with statistical noise for short timescales. Forced symmetrization, by using the matrix T_sym = 1/2 (T + transpose(T)), removes problems like the existence of (small magnitude) negative eigenvalues.

    mininum-var-min
        (100)
        The min-variance strategy (option of :term:`lmc-stats` is only valid for larger number of samples, and can get stuck if too few samples are used at each state. :term:`mininum-var-min` is the minimum number of samples that each state that are allowed before the min-variance strategy is activated if selected.

    init-lambda-weights:
        The initial weights (free energies) used for the expanded ensemble states. Default is a vector of zero weights. format is similar to the lambda vector settings in :term:`fep-lambdas`, except the weights can be any floating point number. Units are kT. Its length must match the lambda vector lengths.

    lmc-weights-equil
        no
            Expanded ensemble weights continue to be updated throughout the simulation.

        yes
            The input expanded ensemble weights are treated as equilibrated, and are not updated throughout the simulation.

        wl-delta
            Expanded ensemble weight updating is stopped when the Wang-Landau incrementor falls below this value.

        number-all-lambda
            Expanded ensemble weight updating is stopped when the number of samples at all of the lambda states is greater than this value.

        number-steps
            Expanded ensemble weight updating is stopped when the number of steps is greater than the level specified by this value.

        number-samples
            Expanded ensemble weight updating is stopped when the number of total samples across all lambda states is greater than the level specified by this value.

        count-ratio
            Expanded ensemble weight updating is stopped when the ratio of samples at the least sampled lambda state and most sampled lambda state greater than this value.

    simulated-tempering
        (no)
        Turn simulated tempering on or off. Simulated tempering is implemented as expanded ensemble sampling with different temperatures instead of different Hamiltonians.

    sim-temp-low
        (300) \[K\]
        Low temperature for simulated tempering.

    sim-temp-high
        (300) \[K\]
        High temperature for simulated tempering.

    simulated-tempering-scaling
        Controls the way that the temperatures at intermediate lambdas are calculated from the :term:`temperature-lambdas` part of the lambda vector.

        linear
            Linearly interpolates the temperatures using the values of :term:`temperature-lambdas`, *i.e.* if :term:`sim-temp-low` =300, :term:`sim-temp-high` =400, then lambda=0.5 correspond to a temperature of 350. A nonlinear set of temperatures can always be implemented with uneven spacing in lambda.

        geometric
            Interpolates temperatures geometrically between :term:`sim-temp-low` and :term:`sim-temp-high`. The i:th state has temperature :term:`sim-temp-low` * (:term:`sim-temp-high` / :term:`sim-temp-low`) raised to the power of (i/(ntemps-1)). This should give roughly equal exchange for constant heat capacity, though of course things simulations that involve protein folding have very high heat capacity peaks.

        exponential
            Interpolates temperatures exponentially between :term:`sim-temp-low` and :term:`sim-temp-high`. The i:th state has temperature :term:`sim-temp-low` + (:term:`sim-temp-high` - :term:`sim-temp-low`)*((exp(:term:`temperature-lambdas` (i))-1)/(exp(1.0)-i)).


Non-equilibrium MD
^^^^^^^^^^^^^^^^^^

.. glossary::

    acc-grps
        groups for constant acceleration (*e.g.* ``Protein Sol``) all atoms in groups Protein and Sol will experience constant acceleration as specified in the :term:`accelerate` line

    accelerate
        (0) \[nm ps^-2\]
        acceleration for :term:`acc-grps`; x, y and z for each group (*e.g.* ``0.1 0.0 0.0 -0.1 0.0 0.0`` means that first group has constant acceleration of 0.1 nm ps-2 in X direction, second group the opposite).

    freezegrps
        Groups that are to be frozen (*i.e.* their X, Y, and/or Z position will not be updated; *e.g.* ``Lipid SOL``). :term:`freezedim` specifies for which dimension the freezing applies. To avoid spurious contibrutions to the virial and pressure due to large forces between completely frozen atoms you need to use energy group exclusions, this also saves computing time. Note that coordinates of frozen atoms are not scaled by pressure-coupling algorithms.

    freezedim
        dimensions for which groups in :term:`freezegrps` should be frozen, specify `Y` or `N` for X, Y and Z and for each group (*e.g.* ``Y Y N N N N`` means that particles in the first group can move only in Z direction. The particles in the second group can move in any direction).

    cos-acceleration
        (0) \[nm ps^-2\]
        the amplitude of the acceleration profile for calculating the viscosity. The acceleration is in the X-direction and the magnitude is :term:`cos-acceleration` cos(2 pi z/boxheight). Two terms are added to the energy file: the amplitude of the velocity profile and 1/viscosity.

    deform
        (0 0 0 0 0 0) \[nm ps-1\]
        The velocities of deformation for the box elements: a(x) b(y) c(z) b(x) c(x) c(y). Each step the box elements for which :term:`deform` is non-zero are calculated as: box(ts)+(t-ts)*deform, off-diagonal elements are corrected for periodicity. The coordinates are transformed accordingly. Frozen degrees of freedom are (purposely) also transformed. The time ts is set to t at the first step and at steps at which x and v are written to trajectory to ensure exact restarts. Deformation can be used together with semiisotropic or anisotropic pressure coupling when the appropriate compressibilities are set to zero. The diagonal elements can be used to strain a solid. The off-diagonal elements can be used to shear a solid or a liquid.


Electric fields
^^^^^^^^^^^^^^^

.. glossary::

    E-x ; E-y ; E-z
        If you want to use an electric field in a direction, enter 3 numbers after the appropriate E-direction, the first number: the number of cosines, only 1 is implemented (with frequency 0) so enter 1, the second number: the strength of the electric field in V nm^-1, the third number: the phase of the cosine, you can enter any number here since a cosine of frequency zero has no phase.

    E-xt; E-yt; E-zt:
        not implemented yet


Mixed quantum/classical molecular dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. glossary::

    QMMM
        no
            No QM/MM.

        yes
            Do a QM/MM simulation. Several groups can be described at different QM levels separately. These are specified in the :term:`QMMM-grps` field separated by spaces. The level of *ab initio* theory at which the groups are described is specified by :term:`QMmethod` and :term:`QMbasis` Fields. Describing the groups at different levels of theory is only possible with the ONIOM QM/MM scheme, specified by :term:`QMMMscheme`.

    QMMM-grps
        groups to be descibed at the QM level

    QMMMscheme
        normal
            normal QM/MM. There can only be one :term:`QMMM-grps` that is modelled at the :term:`QMmethod` and :term:`QMbasis` level of *ab initio* theory. The rest of the system is described at the MM level. The QM and MM subsystems interact as follows: MM point charges are included in the QM one-electron hamiltonian and all Lennard-Jones interactions are described at the MM level.

        ONIOM
            The interaction between the subsystem is described using the ONIOM method by Morokuma and co-workers. There can be more than one :term:`QMMM-grps` each modeled at a different level of QM theory (:term:`QMmethod` and :term:`QMbasis`).

    QMmethod
        (RHF)
        Method used to compute the energy and gradients on the QM atoms. Available methods are AM1, PM3, RHF, UHF, DFT, B3LYP, MP2, CASSCF, and MMVB. For CASSCF, the number of electrons and orbitals included in the active space is specified by :term:`CASelectrons` and :term:`CASorbitals`.

    QMbasis
        (STO-3G)
        Basis set used to expand the electronic wavefuntion. Only Gaussian basis sets are currently available, *i.e.* ``STO-3G, 3-21G, 3-21G*, 3-21+G*, 6-21G, 6-31G, 6-31G*, 6-31+G*,`` and ``6-311G``.

    QMcharge
        (0) \[integer\]
        The total charge in `e` of the :term:`QMMM-grps`. In case there are more than one :term:`QMMM-grps`, the total charge of each ONIOM layer needs to be specified separately.

    QMmult
        (1) \[integer\]
        The multiplicity of the :term:`QMMM-grps`. In case there are more than one :term:`QMMM-grps`, the multiplicity of each ONIOM layer needs to be specified separately.

    CASorbitals
        (0) \[integer\]
        The number of orbitals to be included in the active space when doing a CASSCF computation.

    CASelectrons
        (0) \[integer\]
        The number of electrons to be included in the active space when doing a CASSCF computation.

    SH
        no
            No surface hopping. The system is always in the electronic ground-state.

        yes
            Do a QM/MM MD simulation on the excited state-potential energy surface and enforce a *diabatic* hop to the ground-state when the system hits the conical intersection hyperline in the course the simulation. This option only works in combination with the CASSCF method.


Implicit solvent
^^^^^^^^^^^^^^^^

.. glossary::

    implicit-solvent
        no
            No implicit solvent

        GBSA
            Do a simulation with implicit solvent using the Generalized Born formalism. Three different methods for calculating the Born radii are available, Still, HCT and OBC. These are specified with the :term:`gb-algorithm` field. The non-polar solvation is specified with the :term:`sa-algorithm` field.

    gb-algorithm
        Still
            Use the Still method to calculate the Born radii

        HCT
            Use the Hawkins-Cramer-Truhlar method to calculate the Born radii

        OBC
            Use the Onufriev-Bashford-Case method to calculate the Born radii

    nstgbradii
        (1) \[steps\]
        Frequency to (re)-calculate the Born radii. For most practial purposes, setting a value larger than 1 violates energy conservation and leads to unstable trajectories.

    rgbradii
        (1.0) \[nm\]
        Cut-off for the calculation of the Born radii. Currently must be equal to rlist

    gb-epsilon-solvent
        (80)
        Dielectric constant for the implicit solvent

    gb-saltconc
        (0) \[M\]
        Salt concentration for implicit solvent models, currently not used

    gb-obc-alpha
    gb-obc-beta
    gb-obc-gamma
        Scale factors for the OBC model. Default values of 1, 0.78 and 4.85 respectively are for OBC(II). Values for OBC(I) are 0.8, 0 and 2.91 respectively

    gb-dielectric-offset
        (0.009) \[nm\]
        Distance for the di-electric offset when calculating the Born radii. This is the offset between the center of each atom the center of the polarization energy for the corresponding atom

    sa-algorithm
        Ace-approximation
            Use an Ace-type approximation

        None
            No non-polar solvation calculation done. For GBSA only the polar part gets calculated

    sa-surface-tension
        \[kJ mol-1 nm-2\]
        Default value for surface tension with SA algorithms. The default value is -1; Note that if this default value is not changed it will be overridden by `gmx grompp`_ using values that are specific for the choice of radii algorithm (0.0049 kcal/mol/Angstrom^2 for Still, 0.0054 kcal/mol/Angstrom2 for HCT/OBC) Setting it to 0 will while using an sa-algorithm other than None means no non-polar calculations are done.


Adaptive Resolution Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. glossary::

    adress
        (no)
        Decide whether the AdResS feature is turned on.

    adress-type
        Off
            Do an AdResS simulation with weight equal 1, which is equivalent to an explicit (normal) MD simulation. The difference to disabled AdResS is that the AdResS variables are still read-in and hence are defined.

        Constant
            Do an AdResS simulation with a constant weight, :term:`adress-const-wf` defines the value of the weight

        XSplit
            Do an AdResS simulation with simulation box split in x-direction, so basically the weight is only a function of the x coordinate and all distances are measured using the x coordinate only.

        Sphere
            Do an AdResS simulation with spherical explicit zone.

    adress-const-wf
        (1)
        Provides the weight for a constant weight simulation (:term:`adress-type` =Constant)

    adress-ex-width
        (0)
        Width of the explicit zone, measured from :term:`adress-reference-coords`.

    adress-hy-width
        (0)
        Width of the hybrid zone.

    adress-reference-coords
        (0,0,0)
        Position of the center of the explicit zone. Periodic boundary conditions apply for measuring the distance from it.

    adress-cg-grp-names
        The names of the coarse-grained energy groups. All other energy groups are considered explicit and their interactions will be automatically excluded with the coarse-grained groups.

    adress-site
        The mapping point from which the weight is calculated.

        COM
           The weight is calculated from the center of mass of each charge group.

        COG
           The weight is calculated from the center of geometry of each charge group.

        Atom
           The weight is calculated from the position of 1st atom of each charge group.

        AtomPerAtom
           The weight is calculated from the position of each individual atom.

    adress-interface-correction
        Off
            Do not a apply any interface correction.

        thermoforce
            Apply thermodynamic force interface correction. The table can be specified using the ``-tabletf`` option of mdrun_. The table should contain the potential and force (acting on molecules) as function of the distance from :term:`adress-reference-coords`.

    adress-tf-grp-names
        The names of the energy groups to which the thermoforce is applied if enabled in :term:`adress-interface-correction`. If no group is given the default table is applied.

    adress-ex-forcecap
        (0)
        Cap the force in the hybrid region, useful for big molecules. 0 disables force capping.


User defined thingies
^^^^^^^^^^^^^^^^^^^^^

.. glossary::

    user1-grps; user2-grps; userint1 (0); userint2 (0); userint3 (0); userint4 (0); userreal1 (0); userreal2 (0); userreal3 (0); userreal4 (0)
        These you can use if you modify code. You can pass integers and reals to your subroutine. Check the inputrec definition in ``src/gromacs/legacyheaders/types/inputrec.h``
