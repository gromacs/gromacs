# Molecular dynamics parameters (.mdp options) {#mdp-options}

## General information

Default values are given in parentheses. The first option in the list is
always the default option. Units are given in square brackets The difference
between a dash and an underscore is ignored.

A [sample `.mdp` file](online/mdp.html) is available. This should be appropriate to
start a normal simulation. Edit it to suit your specific needs and desires.

* * *

### Preprocessing

* <a name="include">`include`</a>  
    directories to include in your topology. Format: `-I/home/john/mylib -I../otherlib`

* <a name="define">`define`</a>  
    defines to pass to the preprocessor, default is no defines. You can use any defines to control options in your customized topology files. Options that are already available by default are:  
    `-DFLEXIBLE` will use flexible water instead of rigid water into your topology, this can be useful for normal mode analysis.  
    `-DPOSRES` will trigger the inclusion of `posre.itp` into your topology, used for position restraints.

* * *

### Run control

* <a name="integrator">`integrator`</a>  
    (Despite the name, this list includes algorithms that are not actually integrators. [`steep`](#steep) and all entries following it are in this category)

    + <a name ="md">`md`</a>  
    A leap-frog algorithm for integrating Newton's equations of motion.

    + <a name="md-vv">`md-vv`</a>  
    A velocity Verlet algorithm for integrating Newton's equations of motion. For constant NVE simulations started from corresponding points in the same trajectory, the trajectories are analytically, but not binary, identical to the [`md`](#md) leap-frog integrator. The the kinetic energy, which is determined from the whole step velocities and is therefore slightly too high. The advantage of this integrator is more accurate, reversible Nose-Hoover and Parrinello-Rahman coupling integration based on Trotter expansion, as well as (slightly too small) full step velocity output. This all comes at the cost off extra computation, especially with constraints and extra communication in parallel. Note that for nearly all production simulations the [`md`](#md) integrator is accurate enough. 

    + <a name="md-vv-avek">`md-vv-avek`</a>  
    A velocity Verlet algorithm identical to [`md-vv`](#md-vv), except that the kinetic energy is determined as the average of the two half step kinetic energies as in the [`md`](#md) integrator, and this thus more accurate. With Nose-Hoover and/or Parrinello-Rahman coupling this comes with a slight increase in computational cost.

    + <a name="sd">`sd`</a>  
     An accurate and efficient leap-frog stochastic dynamics integrator. With constraints, coordinates needs to be constrained twice per integration step. Depending on the computational cost of the force calculation, this can take a significant part of the simulation time. The temperature for one or more groups of atoms ([`tc-grps`](#tc-grps)) is set with [`ref-t`](#ref-t), the inverse friction constant for each group is set with [`tau-t`](#tau-t). The parameter [`tcoupl`](#tcoupl) is ignored. The random generator is initialized with [`ld-seed`](#ld-seed). When used as a thermostat, an appropriate value for [`tau-t`](#tau-t) is 2 ps, since this results in a friction that is lower than the internal friction of water, while it is high enough to remove excess heat NOTE: temperature deviations decay twice as fast as with a Berendsen thermostat with the same [`tau-t`](#tau-t).

     + <a name="sd2">`sd2`</a>  
     This used to be the default sd integrator, but is now deprecated. Four Gaussian random numbers are required per coordinate per step. With constraints, the temperature will be slightly too high.

     + <a name="bd">`bd`</a>  
    An Euler integrator for Brownian or position Langevin dynamics, the velocity is the force divided by a friction coefficient ([`bd-fric`](#bd-fric) [amu ps-1]) plus random thermal noise ([`ref-t`](#ref-t)). When [`bd-fric`](#bd-fric) is zero, the friction coefficient for each particle is calculated as `mass/tau-t`, as for the integrator [`sd`](#sd). The random generator is initialized with [`ld-seed`](#ld-seed).

    + <a name="steep">`steep`</a>  
    A steepest descent algorithm for energy minimization. The maximum step size is [`emstep`](#emstep), the tolerance is [`emtol`](#emtol).

    + <a name="cg">`cg`</a>  
    A conjugate gradient algorithm for energy minimization, the tolerance is [`emtol`](#emtol). CG is more efficient when a steepest descent step is done every once in a while, this is determined by [`nstcgsteep`](#nstcgsteep). For a minimization prior to a normal mode analysis, which requires a very high accuracy, GROMACS should be compiled in double precision.

    + <a name="l-bfgs">`l-bfgs`</a>  
    A quasi-Newtonian algorithm for energy minimization according to the low-memory Broyden-Fletcher-Goldfarb-Shanno approach. In practice this seems to converge faster than Conjugate Gradients, but due to the correction steps necessary it is not (yet) parallelized. 

    + <a name="nm">`nm`</a>  
    Normal mode analysis is performed on the structure in the [.tpr] file. GROMACS should be compiled in double precision.

    + <a name="tpi">`tpi`</a>  
     Test particle insertion. The last molecule in the topology is the test particle. A trajectory must be provided to [mdrun] by doing a [rerun](#rerun). This trajectory should not contain the molecule to be inserted. Insertions are performed [`nsteps`](#nsteps) times in each frame at random locations and with random orientiations of the molecule. When [`nstlist`](#nstlist) is larger than one, [`nstlist`](#nstlist) insertions are performed in a sphere with radius [`rtpi`](#rtpi) around a the same random location using the same neighborlist (and the same long-range energy when [`rvdw`](#rvdw) or [`rcoulomb`](#rcoulomb) is greater than [`rlist`](#rlist), which is only allowed for single-atom molecules). Since neighborlist construction is expensive, one can perform several extra insertions with the same list almost for free. The random seed is set with [`ld-seed`](#ld-seed). The temperature for the Boltzmann weighting is set with [`ref-t`](#ref-t), this should match the temperature of the simulation of the original trajectory. Dispersion correction is implemented correctly for TPI. All relevant quantities are written to the file specified with `mdrun -tpi`. The distribution of insertion energies is written to the file specified with `mdrun -tpid`. No trajectory or energy file is written. Parallel TPI gives identical results to single-node TPI. For charged molecules, using PME with a fine grid is most accurate and also efficient, since the potential in the system only needs to be calculated once per frame.

     + <a name="tpic">`tpic`</a>  
     Test particle insertion into a predefined cavity location. The procedure is the same as for [`tpi`](#tpi), except that one coordinate extra is read from the trajectory, which is used as the insertion location. The molecule to be inserted should be centered at 0,0,0. Gromacs does not do this for you, since for different situations a different way of centering might be optimal. Also [`rtpi`](#rtpi) sets the radius for the sphere around this location. Neighbor searching is done only once per frame, [`nstlist`](#nstlist) is not used. Parallel [`tpic`](#tpic) gives identical results to single-rank [`tpic`](#tpic).

* <a name="tinit">`tinit`</a> (0) [ps]  
    starting time for your run (only makes sense for integrators [`md`](#md), [`sd`](#sd) and [`bd`](#bd))

* <a name="dt">`dt`</a> (0.001) [ps]  
    time step for integration (only makes sense for integrators [`md`](#md), [`sd`](#sd) and [`bd`](#bd))

* <a name="nsteps">`nsteps`</a> (0)  
    maximum number of steps to integrate or minimize, -1 is no maximum

* <a name="init-step">`init-step`</a> (0)  
    The starting step. The time at an step i in a run is calculated as: `t = tinit + dt * (init-step + i)`. The free-energy lambda is calculated as: `lambda = init-lambda + delta-lambda * (init-step + i)`. Also non-equilibrium MD parameters can depend on the step number. Thus for exact restarts or redoing part of a run it might be necessary to set [`init-step`](#init-step) to the step number of the restart frame. [gmx convert-tpr] does this automatically.

* <a name="comm-mode">`comm-mode`</a>

    + `Linear`  
    Remove center of mass translation

    + `Angular`  
    Remove center of mass translation and rotation around the center of mass 

    + `None`  
    No restriction on the center of mass motion 

* <a name="nstcomm">`nstcomm`</a> (100) [steps]  
    frequency for center of mass motion removal

* <a name="comm-grps">`comm-grps`</a>
    group(s) for center of mass motion removal, default is the whole system

* * *

### Langevin dynamics

* <a name="bd-fric">`bd-fric`</a> (0) [amu ps-1]  
    Brownian dynamics friction coefficient. When [`bd-fric`](#bd-fric) =0, the friction coefficient for each particle is calculated as `mass/tau-t`.

* <a name="ld-seed">`ld-seed`</a> (-1) [integer]  
    used to initialize random generator for thermal noise for stochastic and Brownian dynamics. When [`ld-seed`](#ld-seed) is set to -1, a pseudo random seed is used. When running BD or SD on multiple processors, each processor uses a seed equal to [`ld-seed`](#ld-seed) plus the processor number.

* * *

### Energy minimization

* <a name="emtol">`emtol`</a> (10.0) [kJ mol-1 nm-1]  
    the minimization is converged when the maximum force is smaller than this value

* <a name="emstep">`emstep`</a> (0.01) [nm]  
    initial step-size

* <a name="nstcgsteep">`nstcgsteep`</a> (1000) [steps]  
    frequency of performing 1 steepest descent step while doing conjugate gradient energy minimization.

* <a name="nbfgscorr">`nbfgscorr`</a> (10)  
    Number of correction steps to use for L-BFGS minimization. A higher number is (at least theoretically) more accurate, but slower.

* * *

### Shell Molecular Dynamics

When shells or flexible constraints are present in the system the positions of
the shells and the lengths of the flexible constraints are optimized at every
time step until either the RMS force on the shells and constraints is less
than emtol, or a maximum number of iterations [`niter`](#niter) has been reached

* <a name="emtol">`emtol`</a> (10.0) [kJ mol-1 nm-1]  
    the minimization is converged when the maximum force is smaller than this value. For shell MD this value should be 1.0 at most, but since the variable is used for energy minimization as well the default is 10.0.

* <a name="niter">`niter`</a> (20)  
    maximum number of iterations for optimizing the shell positions and the flexible constraints.

* <a name="fcstep">`fcstep`</a> (0) [ps2]  
    the step size for optimizing the flexible constraints. Should be chosen as mu/(d2V/dq2) where mu is the reduced mass of two particles in a flexible constraint and d2V/dq2 is the second derivative of the potential in the constraint direction. Hopefully this number does not differ too much between the flexible constraints, as the number of iterations and thus the runtime is very sensitive to [`fcstep`](#fcstep). Try several values!
  

* * *

### Test particle insertion

* <a name="rtpi">`rtpi`</a> (0.05) [nm]
    the test particle insertion radius see integrators [`tpi`](#tpi) and [`tpic`](#tpic)

* * *

### Output control

* <a name="nstxout">`nstxout`</a> (0) [steps]  
    number of steps that elapse between writing coordinates to output trajectory file, the last coordinates are always written

* <a name="nstvout">`nstvout`</a> (0) [steps]  
    number of steps that elapse between writing velocities to output trajectory, the last velocities are always written

* <a name="nstfout">`nstfout`</a> (0) [steps]  
    number of steps that elapse between writing forces to output trajectory.

* <a name="nstlog">`nstlog`</a> (1000) [steps]  
    number of steps that elapse between writing energies to the log file, the last energies are always written

* <a name="nstcalcenergy">`nstcalcenergy`</a> (100)  
    number of steps that elapse between calculating the energies, 0 is never. This option is only relevant with dynamics. With a twin-range cut-off setup [`nstcalcenergy`](#nstcalcenergy) should be equal to or a multiple of [`nstlist`](#nstlist). This option affects the performance in parallel simulations, because calculating energies requires global communication between all processes which can become a bottleneck at high parallelization. 

* <a name="nstenergy">`nstenergy`</a> (1000) [steps]  
    number of steps that else between writing energies to energy file, the last energies are always written, should be a multiple of [`nstcalcenergy`](#nstcalcenergy). Note that the exact sums and fluctuations over all MD steps modulo [`nstcalcenergy`](#nstcalcenergy) are stored in the energy file, so [gmx energy] can report exact energy averages and fluctuations also when [`nstenergy`](#nstenergy) is greater than 1.

* <a name="nstxout-compressed">`nstxout-compressed`</a> (0) [steps]  
    number of steps that elapse between writing position coordinates using lossy compression

* <a name="compressed-x-precision">`compressed-x-precision`</a> (1000) [real]  
    precision with which to write to the compressed trajectory file

* <a name="compressed-x-grps">`compressed-x-grps`</a>  
    group(s) to write to the compressed trajectory file, by default the whole system is written (if [`nstxout-compressed`](#nstxout-compressed) is positive)

* <a name="energygrps">`energygrps`</a>  
    group(s) to write to energy file
  

* * *

### Neighbor searching

* <a name="cutoff-scheme">`cutoff-scheme`</a>  

    + <a name="Verlet">`Verlet`</a>  
    Generate a pair list with buffering. The buffer size is automatically set based on [`verlet-buffer-tolerance`](#verlet-buffer-tolerance), unless this is set to -1, in which case [`rlist`](#rlist) will be used. This option has an explicit, exact cut-off when [`rvdw`](#rvdw) equals [`rcoulomb`](#rcoulomb). Currently only cut-off, reaction-field, PME electrostatics and plain and PME LJ are supported. Some [mdrun] functionality is not yet supported with the [`Verlet`](#Verlet) scheme, but [gmx grompp] checks for this. Native GPU acceleration is only supported with [`Verlet`](#Verlet). With GPU-accelerated PME or with separate PME ranks, [mdrun] will automatically tune the CPU/GPU load balance by scaling [`rcoulomb`](#rcoulomb) and the grid spacing. This can be turned off with `mdrun -notunepme`. [`Verlet`](#Verlet) is faster than [`group`](#group) when there is no water, or if [`group`](#group) would use a pair-list buffer to conserve energy.

    + <a name="group">`group`</a>  
    Generate a pair list for groups of atoms. These groups correspond to the charge groups in the topology. This was the only cut-off treatment scheme before version 4.6, and is **deprecated in 5.0**. There is no explicit buffering of the pair list. This enables efficient force calculations for water, but energy is only conserved when a buffer is explicitly added.

* <a name="nstlist">`nstlist`</a> (10) [steps]  

    + `>0`  
    Frequency to update the neighbor list (and the long-range forces, when using twin-range cut-offs). When this is 0, the neighbor list is made only once. With energy minimization the neighborlist will be updated for every energy evaluation when [`nstlist`](#nstlist) is positive. With [`cutoff-scheme=Verlet`](#cutoff-scheme=Verlet) and [`verlet-buffer-tolerance`](#verlet-buffer-tolerance) set, [`nstlist`](#nstlist) is actually a minimum value and [mdrun] might increase it, unless it is set to 1. With parallel simulations and/or non-bonded force calculation on the GPU, a value of 20 or 40 often gives the best performance. With [`cutoff-scheme=Group`](#cutoff-scheme=Group) and non-exact cut-off's, [`nstlist`](#nstlist) will affect the accuracy of your simulation and it can not be chosen freely. 

    + `0`  
    The neighbor list is only constructed once and never updated. This is mainly useful for vacuum simulations in which all particles see each other.

    + `-1`  
    Automated update frequency, only supported with the [`group`](#group) [`cutoff-scheme`](#cutoff-scheme). This can only be used with switched, shifted or user potentials where the cut-off can be smaller than [`rlist`](#rlist). One then has a buffer of size [`rlist`](#rlist) minus the longest cut-off. The neighbor list is only updated when one or more particles have moved further than half the buffer size from the center of geometry of their charge group as determined at the previous neighbor search. Coordinate scaling due to pressure coupling or the [`deform`](#deform) option is taken into account. This option guarantees that their are no cut-off artifacts, but for larger systems this can come at a high computational cost, since the neighbor list update frequency will be determined by just one or two particles moving slightly beyond the half buffer length (which does not necessarily imply that the neighbor list is invalid), while 99.99% of the particles are fine. 

* <a name="nstcalclr">`nstcalclr`</a> (-1) [steps]  
     Controls the period between calculations of long-range forces when using the group cut-off scheme. 

    + `1`  
    Calculate the long-range forces every single step. This is useful to have separate neighbor lists with buffers for electrostatics and Van der Waals interactions, and in particular it makes it possible to have the Van der Waals cutoff longer than electrostatics (useful _e.g._ with PME). However, there is no point in having identical long-range cutoffs for both interaction forms and update them every step - then it will be slightly faster to put everything in the short-range list.

    + `>1`  
    Calculate the long-range forces every [`nstcalclr`](#nstcalclr) steps and use a multiple-time-step integrator to combine forces. This can now be done more frequently than [`nstlist`](#nstlist) since the lists are stored, and it might be a good idea _e.g._ for Van der Waals interactions that vary slower than electrostatics.

    + `-1`  
    Calculate long-range forces on steps where neighbor searching is performed. While this is the default value, you might want to consider updating the long-range forces more frequently.
Note that twin-range force evaluation might be enabled automatically by PP-PME
load balancing. This is done in order to maintain the chosen Van der Waals
interaction radius even if the load balancing is changing the electrostatics
cutoff. If the [.mdp] file already specifies twin-range interactions (_e.g._ to
evaluate Lennard-Jones interactions with a longer cutoff than the PME
electrostatics every 2-3 steps), the load balancing will have also a small
effect on Lennard-Jones, since the short-range cutoff (inside which forces are
evaluated every step) is changed.

* <a name="ns-type">`ns-type`</a>  
    + `grid`  
    Make a grid in the box and only check atoms in neighboring grid cells when constructing a new neighbor list every [`nstlist`](#nstlist) steps. In large systems grid search is much faster than simple search.

    + `simple`  
    Check every atom in the box when constructing a new neighbor list every [`nstlist`](#nstlist) steps (only with [`group`](#group) [`cutoff-scheme`](#cutoff-scheme)).

* <a name="pbc">`pbc`</a>  

    + `xyz`  
    Use periodic boundary conditions in all directions.

    + `no`  
    Use no periodic boundary conditions, ignore the box. To simulate without cut-offs, set all cut-offs and [`nstlist`](#nstlist) to 0. For best performance without cut-offs on a single MPI rank, set [`nstlist`](#nstlist) to zero and `ns-type=simple`.

    + `xy`  
    Use periodic boundary conditions in x and y directions only. This works only with `ns-type=grid` and can be used in combination with [`walls`](#walls). Without walls or with only one wall the system size is infinite in the z direction. Therefore pressure coupling or Ewald summation methods can not be used. These disadvantages do not apply when two walls are used.

* <a name="periodic-molecules">`periodic-molecules`</a>  
    + `no`  
    molecules are finite, fast molecular PBC can be used

    + `yes`  
    for systems with molecules that couple to themselves through the periodic boundary conditions, this requires a slower PBC algorithm and molecules are not made whole in the output

* <a name="verlet-buffer-tolerance">`verlet-buffer-tolerance`</a> (0.005) [kJ/mol/ps]  
    Useful only with [`Verlet`](#Verlet) [`cutoff-scheme`](#cutoff-scheme). This sets the maximum allowed error for pair interactions per particle caused by the Verlet buffer, which indirectly sets [`rlist`](#rlist). As both [`nstlist`](#nstlist) and the Verlet buffer size are fixed (for performance reasons), particle pairs not in the pair list can occasionally get within the cut-off distance during the `nstlist-1` steps between list updates. This causes very small jumps in the energy. In a constant-temperature ensemble, these very small energy jumps can be estimated for a given cut-off and [`rlist`](#rlist). The estimate assumes a homogeneous particle distribution, hence the errors might be slightly underestimated for multi-phase systems. For longer pair-list life-time `(nstlist-1)*dt`, the buffer is overestimated, because the interactions between particles are ignored. Combined with cancellation of errors, the actual drift of the total energy is usually one to two orders of magnitude smaller. Note that the generated buffer size takes into account that the GROMACS pair-list setup leads to a reduction in the drift by a factor 10, compared to a simple particle-pair based list. Without dynamics (energy minimization etc.), the buffer is 5% of the cut-off. For NVE simulations the initial temperature is used, unless this is zero, in which case a buffer of 10% is used. For NVE simulations the tolerance usually needs to be lowered to achieve proper energy conservation on the nanosecond time scale. To override the automated buffer setting, use `verlet-buffer-tolerance = -1` and set [`rlist`](#rlist) manually.

* <a name="rlist">`rlist`</a> (1) [nm]  
    Cut-off distance for the short-range neighbor list. With the [`Verlet`](#Verlet) [`cutoff-scheme`](#cutoff-scheme), this is by default set by the [`verlet-buffer-tolerance`](#verlet-buffer-tolerance) option and the value of [`rlist`](#rlist) is ignored.

* <a name="rlistlong">`rlistlong`</a> (-1) [nm]  
    Cut-off distance for the long-range neighbor list. This parameter is only relevant for a twin-range cut-off setup with switched potentials. In that case a buffer region is required to account for the size of charge groups. In all other cases this parameter is automatically set to the longest cut-off distance.
  
* * *

### Electrostatics

* <a name="coulombtype">`coulombtype`</a>  
    + `Cut-off`  
    Twin range cut-offs with neighborlist cut-off [`rlist`](#rlist) and Coulomb cut-off [`rcoulomb`](#rcoulomb), where `rcoulomb >= rlist`.

    + `Ewald`  
    Classical Ewald sum electrostatics. The real-space cut-off [`rcoulomb`](#rcoulomb) should be equal to [`rlist`](#rlist). Use _e.g._ `rlist =0.9, rcoulomb =0.9`. The highest magnitude of wave vectors used in reciprocal space is controlled by [`fourierspacing`](#fourierspacing). The relative accuracy of direct/reciprocal space is controlled by [`ewald-rtol`](#ewald-rtol).  
NOTE: Ewald scales as O(N3/2) and is thus extremely slow for large systems. It
is included mainly for reference - in most cases PME will perform much better.

    + `PME`  
    Fast smooth Particle-Mesh Ewald (SPME) electrostatics. Direct space is similar to the Ewald sum, while the reciprocal part is performed with FFTs. Grid dimensions are controlled with [`fourierspacing`](#fourierspacing) and the interpolation order with [`pme-order`](#pme-order). With a grid spacing of 0.1 nm and cubic interpolation the electrostatic forces have an accuracy of 2-3*10-4. Since the error from the vdw-cutoff is larger than this you might try 0.15 nm. When running in parallel the interpolation parallelizes better than the FFT, so try decreasing grid dimensions while increasing interpolation.

    + `P3M-AD`  
    Particle-Particle Particle-Mesh algorithm with analytical derivative for for long range electrostatic interactions. The method and code is identical to SPME, except that the influence function is optimized for the grid. This gives a slight increase in accuracy.

    + `Reaction-Field` electrostatics  
    Reaction field with Coulomb cut-off [`rcoulomb`](#rcoulomb), where `rcoulomb >= rlist`. The dielectric constant beyond the cut-off is [`epsilon-rf`](#epsilon-rf). The dielectric constant can be set to infinity by setting [`epsilon-rf`](#epsilon-rf) to zero.

    + `Generalized-Reaction-Field`  
    Generalized reaction field with Coulomb cut-off [`rcoulomb`](#rcoulomb), where `rcoulomb >= rlist`. The dielectric constant beyond the cut-off is [`epsilon-rf`](#epsilon-rf). The ionic strength is computed from the number of charged (_i.e._ with non zero charge) charge groups. The temperature for the GRF potential is set with [`ref-t`](#ref-t).

    + `Reaction-Field-zero`  
    In GROMACS, normal reaction-field electrostatics with the [`group`](#group) [`cutoff-scheme`](#cutoff-scheme) leads to bad energy conservation. [`Reaction-Field-zero`](#Reaction-Field-zero) solves this by making the potential zero beyond the cut-off. It can only be used with an infinite dielectric constant (`epsilon-rf = 0`), because only for that value the force vanishes at the cut-off. [`rlist`](#rlist) should be 0.1 to 0.3 nm larger than [`rcoulomb`](#rcoulomb) to accommodate for the size of charge groups and diffusion between neighbor list updates. This, and the fact that table lookups are used instead of analytical functions make [`Reaction-Field-zero`](#Reaction-Field-zero) computationally more expensive than normal reaction-field.

    + `Reaction-Field-nec`  
    The same as [`Reaction-Field`](#Reaction-Field), but implemented as in GROMACS versions before 3.3. No reaction-field correction is applied to excluded atom pairs and self pairs. The 1-4 interactions are calculated using a reaction-field. The missing correction due to the excluded pairs that do not have a 1-4 interaction is up to a few percent of the total electrostatic energy and causes a minor difference in the forces and the pressure.

    + `Shift`  
    Analogous to `Shift` for [`vdwtype`](#vdwtype). You might want to use [`Reaction-Field-zero`](#Reaction-Field-zero) instead, which has a similar potential shape, but has a physical interpretation and has better energies due to the exclusion correction terms. 

    + `Encad-Shift`  
    The Coulomb potential is decreased over the whole range, using the definition from the Encad simulation package.

    + `Switch`  
    Analogous to `Switch` for [`vdwtype`](#vdwtype). Switching the Coulomb potential can lead to serious artifacts, advice: use [`Reaction-Field-zero`](#Reaction-Field-zero) instead.

    + `User`  
    [mdrun] will now expect to find a file `table.xvg` with user-defined potential functions for repulsion, dispersion and Coulomb. When pair interactions are present, [mdrun] also expects to find a file `tablep.xvg` for the pair interactions. When the same interactions should be used for non-bonded and pair interactions the user can specify the same file name for both table files. These files should contain 7 columns: the `x` value, `f(x)`, `-f'(x)`, `g(x)`, `-g'(x)`, `h(x)`, `-h'(x)`, where `f(x)` is the Coulomb function, `g(x)` the dispersion function and `h(x)` the repulsion function. When [`vdwtype`](#vdwtype) is not set to `User` the values for `g`, `-g'`, `h` and `-h'` are ignored. For the non-bonded interactions `x` values should run from 0 to the largest cut-off distance + [`table-extension`](#table-extension) and should be uniformly spaced. For the pair interactions the table length in the file will be used. The optimal spacing, which is used for non-user tables, is `0.002 nm` when you run in mixed precision or `0.0005 nm` when you run in double precision. The function value at `x=0` is not important. More information is in the printed manual.

    + `PME-Switch`  
    A combination of PME and a switch function for the direct-space part (see above). [`rcoulomb`](#rcoulomb) is allowed to be smaller than [`rlist`](#rlist). This is mainly useful constant energy simulations (note that using [`PME`](#PME) with the [`Verlet`](#Verlet) [`cutoff-scheme`](#cutoff-scheme) will be more efficient). 

    + `PME-User`  
    A combination of PME and user tables (see above). [`rcoulomb`](#rcoulomb) is allowed to be smaller than [`rlist`](#rlist). The PME mesh contribution is subtracted from the user table by [mdrun]. Because of this subtraction the user tables should contain about 10 decimal places.

    + `PME-User-Switch`  
    A combination of PME-User and a switching function (see above). The switching function is applied to final particle-particle interaction, _i.e._ both to the user supplied function and the PME Mesh correction part.

* <a name="coulomb-modifier">`coulomb-modifier`</a>  

    + `Potential-shift-Verlet`  
    Selects [`Potential-shift`](#Potential-shift) with the Verlet cutoff-scheme, as it is (nearly) free; selects `None` with the group cutoff-scheme.

    + `Potential-shift`  
    Shift the Coulomb potential by a constant such that it is zero at the cut-off. This makes the potential the integral of the force. Note that this does not affect the forces or the sampling.

    + `None`  
    Use an unmodified Coulomb potential. With the group scheme this means no exact cut-off is used, energies and forces are calculated for all pairs in the neighborlist.

* <a name="rcoulomb-switch">`rcoulomb-switch`</a> (0) [nm]  
    where to start switching the Coulomb potential, only relevant when force or potential switching is used

* <a name="rcoulomb">`rcoulomb`</a> (1) [nm]  
    distance for the Coulomb cut-off

* <a name="epsilon-r">`epsilon-r`</a> (1)  
    The relative dielectric constant. A value of 0 means infinity.

* <a name="epsilon-rf">`epsilon-rf`</a> (0)  
    The relative dielectric constant of the reaction field. This is only used with reaction-field electrostatics. A value of 0 means infinity.

* * *

### VdW

* <a name="vdwtype">`vdwtype`</a>  

    + `Cut-off`    
    Twin range cut-offs with neighbor list cut-off [`rlist`](#rlist) and VdW cut-off [`rvdw`](#rvdw), where `rvdw >= rlist`.

    + `PME`    
    Fast smooth Particle-mesh Ewald (SPME) for VdW interactions. The grid dimensions are controlled with [`fourierspacing`](#fourierspacing) in the same way as for electrostatics, and the interpolation order is controlled with [`pme-order`](#pme-order). The relative accuracy of direct/reciprocal space is controlled by [`ewald-rtol-lj`](#ewald-rtol-lj), and the specific combination rules that are to be used by the reciprocal routine are set using [`lj-pme-comb-rule`](#lj-pme-comb-rule).

    + `Shift`    
    This functionality is deprecated and replaced by the `Force-switch` setting for [`vdw-modifier`](#vdw-modifier). The LJ (not Buckingham) potential is decreased over the whole range and the forces decay smoothly to zero between [`rvdw-switch`](#rvdw-switch) and [`rvdw`](#rvdw). The neighbor search cut-off [`rlist`](#rlist) should be 0.1 to 0.3 nm larger than [`rvdw`](#rvdw) to accommodate for the size of charge groups and diffusion between neighbor list updates.

    + `Switch`    
    This functionality is deprecated and replaced by `vdw-modifier = Potential-switch`. The LJ (not Buckingham) potential is normal out to [`rvdw-switch`](#rvdw-switch), after which it is switched off to reach zero at [`rvdw`](#rvdw). Both the potential and force functions are continuously smooth, but be aware that all switch functions will give rise to a bulge (increase) in the force (since we are switching the potential). The neighbor search cut-off [`rlist`](#rlist) should be 0.1 to 0.3 nm larger than [`rvdw`](#rvdw) to accommodate for the size of charge groups and diffusion between neighbor list updates.

    + `Encad-Shift`    
    The LJ (not Buckingham) potential is decreased over the whole range, using the definition from the Encad simulation package.

    + `User`    
    See `User` for [`coulombtype`](#coulombtype). The function value at zero is not important. When you want to use LJ correction, make sure that [`rvdw`](#rvdw) corresponds to the cut-off in the user-defined function. When [`coulombtype`](#coulombtype) is not set to `User` the values for `f` and `-f'` are ignored.

* <a name="vdw-modifier">`vdw-modifier`</a>  

    + `Potential-shift-Verlet`    
    Selects `Potential-shift` with the Verlet cutoff-scheme, as it is (nearly) free; selects `None` with the group cutoff-scheme.

    + `Potential-shift`    
    Shift the Van der Waals potential by a constant such that it is zero at the cut-off. This makes the potential the integral of the force. Note that this does not affect the forces or the sampling.

    + `None`    
    Use an unmodified Van der Waals potential. With the group scheme this means no exact cut-off is used, energies and forces are calculated for all pairs in the neighborlist.

    + `Force-switch`    
    Smoothly switches the forces to zero between [`rvdw-switch`](#rvdw-switch) and [`rvdw`](#rvdw). This shifts the potential shift over the whole range and switches it to zero at the cut-off. Note that this is more expensive to calculate than a plain cut-off and it is not required for energy conservation, since `Potential-shift` conserves energy just as well.

    + `Potential-switch`    
    Smoothly switches the potential to zero between [`rvdw-switch`](#rvdw-switch) and [`rvdw`](#rvdw). Note that this introduces articifically large forces in the switching region and is much more expensive to calculate. This option should only be used if the force field you are using requires this.

* <a name="rvdw-switch">`rvdw-switch`</a> (0) [nm]  
    where to start switching the LJ force and possibly the potential, only relevant when force or potential switching is used

* <a name="rvdw">`rvdw`</a> (1) [nm]  
    distance for the LJ or Buckingham cut-off

* <a name="DispCorr">`DispCorr`</a>  

    + `no`    
    don't apply any correction

    + `EnerPres`    
    apply long range dispersion corrections for Energy and Pressure

    + `Ener`    
    apply long range dispersion corrections for Energy only

* * *

### Tables

* <a name="table-extension">`table-extension`</a> (1) [nm]  
    Extension of the non-bonded potential lookup tables beyond the largest cut-off distance. The value should be large enough to account for charge group sizes and the diffusion between neighbor-list updates. Without user defined potential the same table length is used for the lookup tables for the 1-4 interactions, which are always tabulated irrespective of the use of tables for the non-bonded interactions. The value of [`table-extension`](#table-extension) in no way affects the values of [`rlist`](#rlist), [`rcoulomb`](#rcoulomb), or [`rvdw`](#rvdw). 

* <a name="energygrp-table">`energygrp-table`</a>  
    When user tables are used for electrostatics and/or VdW, here one can give pairs of energy groups for which seperate user tables should be used. The two energy groups will be appended to the table file name, in order of their definition in [`energygrps`](#energygrps), seperated by underscores. For example, if `energygrps = Na Cl Sol` and `energygrp-table = Na Na Na Cl`, [mdrun] will read `table_Na_Na.xvg` and `table_Na_Cl.xvg` in addition to the normal `table.xvg` which will be used for all other energy group pairs. 

* * *

### Ewald

* <a name="fourierspacing">`fourierspacing`</a> (0.12) [nm]  
    For ordinary Ewald, the ratio of the box dimensions and the spacing determines a lower bound for the number of wave vectors to use in each (signed) direction. For PME and P3M, that ratio determines a lower bound for the number of Fourier-space grid points that will be used along that axis. In all cases, the number for each direction can be overridden by entering a non-zero value for [`fourier-nx`](#fourier-nx), etc. For optimizing the relative load of the particle-particle interactions and the mesh part of PME, it is useful to know that the accuracy of the electrostatics remains nearly constant when the Coulomb cut-off and the PME grid spacing are scaled by the same factor.

* <a name="fourier-nx">`fourier-nx`</a> (0)  
  <a name="fourier-ny">`fourier-ny`</a> (0)  
  <a name="fourier-nz">`fourier-nz`</a> (0)  
    Highest magnitude of wave vectors in reciprocal space when using Ewald.
    Grid size when using PME or P3M. These values override [`fourierspacing`](#fourierspacing) per direction. The best choice is powers of 2, 3, 5 and 7. Avoid large primes.

* <a name="pme-order">`pme-order`</a> (4)  
    Interpolation order for PME. 4 equals cubic interpolation. You might try 6/8/10 when running in parallel and simultaneously decrease grid dimension.

* <a name="ewald-rtol">`ewald-rtol`</a> (1e-5)  
    The relative strength of the Ewald-shifted direct potential at [`rcoulomb`](#rcoulomb) is given by [`ewald-rtol`](#ewald-rtol). Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

* <a name="ewald-rtol-lj">`ewald-rtol-lj`</a> (1e-3)  
    When doing PME for VdW-interactions, [`ewald-rtol-lj`](#ewald-rtol-lj) is used to control the relative strength of the dispersion potential at [`rvdw`](#rvdw) in the same way as [`ewald-rtol`](#ewald-rtol) controls the electrostatic potential.

* <a name="lj-pme-comb-rule">`lj-pme-comb-rule`</a> (Geometric)  
    The combination rules used to combine VdW-parameters in the reciprocal part of LJ-PME. Geometric rules are much faster than Lorentz-Berthelot and usually the recommended choice, even when the rest of the force field uses the Lorentz-Berthelot rules.

    + `Geometric`    
    Apply geometric combination rules

    + `Lorentz-Berthelot`    
    Apply Lorentz-Berthelot combination rules

* <a name="ewald-geometry">`ewald-geometry`</a> (3d)    
    
    + `3d`    
    The Ewald sum is performed in all three dimensions.

    + `3dc`    
    The reciprocal sum is still performed in 3D, but a force and potential correction applied in the `z` dimension to produce a pseudo-2D summation. If your system has a slab geometry in the `x-y` plane you can try to increase the `z`-dimension of the box (a box height of 3 times the slab height is usually ok) and use this option.

* <a name="epsilon-surface">`epsilon-surface`</a> (0)  
    This controls the dipole correction to the Ewald summation in 3D. The default value of zero means it is turned off. Turn it on by setting it to the value of the relative permittivity of the imaginary surface around your infinite system. Be careful - you shouldn't use this if you have free mobile charges in your system. This value does not affect the slab 3DC variant of the long range corrections.

* <a name="optimize-fft">`optimize-fft`</a>  

    + `no`    
    Don't calculate the optimal FFT plan for the grid at startup.

    + `yes`    
    Calculate the optimal FFT plan for the grid at startup. This saves a few percent for long simulations, but takes a couple of minutes at start.

* * *

### Temperature coupling

* <a name="tcoupl">`tcoupl`</a>  
    
    + `no`    
    No temperature coupling.

    + `berendsen`    
    Temperature coupling with a Berendsen-thermostat to a bath with temperature [`ref-t`](#ref-t), with time constant [`tau-t`](#tau-t). Several groups can be coupled separately, these are specified in the [`tc-grps`](#tc-grps) field separated by spaces.

    + `nose-hoover`    
    Temperature coupling using a Nose-Hoover extended ensemble. The reference temperature and coupling groups are selected as above, but in this case [`tau-t`](#tau-t) controls the period of the temperature fluctuations at equilibrium, which is slightly different from a relaxation time. For NVT simulations the conserved energy quantity is written to energy and log file.

    + `andersen`    
    Temperature coupling by randomizing a fraction of the particles at each timestep. Reference temperature and coupling groups are selected as above. [`tau-t`](#tau-t) is the average time between randomization of each molecule. Inhibits particle dynamics somewhat, but little or no ergodicity issues. Currently only implemented with velocity Verlet, and not implemented with constraints.

    + `andersen-massive`    
    Temperature coupling by randomizing all particles at infrequent timesteps. Reference temperature and coupling groups are selected as above. [`tau-t`](#tau-t) is the time between randomization of all molecules. Inhibits particle dynamics somewhat, but little or no ergodicity issues. Currently only implemented with velocity Verlet.

    + `v-rescale`    
    Temperature coupling using velocity rescaling with a stochastic term (JCP 126, 014101). This thermostat is similar to Berendsen coupling, with the same scaling using [`tau-t`](#tau-t), but the stochastic term ensures that a proper canonical ensemble is generated. The random seed is set with [`ld-seed`](#ld-seed). This thermostat works correctly even for zero [`tau-t`](#tau-t). For NVT simulations the conserved energy quantity is written to the energy and log file.

* <a name="nsttcouple">`nsttcouple`</a> (-1)  
    The frequency for coupling the temperature. The default value of -1 sets [`nsttcouple`](#nsttcouple) equal to [`nstlist`](#nstlist), unless `nstlist <= 0`, then a value of 10 is used. For velocity-Verlet integrators, [`nsttcouple`](#nsttcouple) is set to 1.

* <a name="nh-chain-length">`nh-chain-length`</a> (10)  
    the number of chained Nose-Hoover thermostats for velocity Verlet integrators, the leap-frog [`md`](#md) integrator only supports 1. Data for the NH chain variables is not printed to the .edr, but can be using the `GMX_NOSEHOOVER_CHAINS` environment variable

* <a name="tc-grps">`tc-grps`</a>  
    groups to couple to separate temperature baths

* <a name="tau-t">`tau-t`</a> [ps]  
    time constant for coupling (one for each group in [`tc-grps`](#tc-grps)), -1 means no temperature coupling

* <a name="ref-t">`ref-t`</a> [K]  
    reference temperature for coupling (one for each group in [`tc-grps`](#tc-grps))
  

* * *

### Pressure coupling

* <a name="pcoupl">`pcoupl`</a>  
    
    + `no`    
    No pressure coupling. This means a fixed box size.

    + `berendsen`    
    Exponential relaxation pressure coupling with time constant [`tau-p`](#tau-p). The box is scaled every timestep. It has been argued that this does not yield a correct thermodynamic ensemble, but it is the most efficient way to scale a box at the beginning of a run.

    + `Parrinello-Rahman`    
    Extended-ensemble pressure coupling where the box vectors are subject to an equation of motion. The equation of motion for the atoms is coupled to this. No instantaneous scaling takes place. As for Nose-Hoover temperature coupling the time constant [`tau-p`](#tau-p) is the period of pressure fluctuations at equilibrium. This is probably a better method when you want to apply pressure scaling during data collection, but beware that you can get very large oscillations if you are starting from a different pressure. For simulations where the exact fluctation of the NPT ensemble are important, or if the pressure coupling time is very short it may not be appropriate, as the previous time step pressure is used in some steps of the GROMACS implementation for the current time step pressure.

    + `MTTK`    
    Martyna-Tuckerman-Tobias-Klein implementation, only useable with [`md-vv`](#md-vv) or [`md-vv-avek`](#md-vv-avek), very similar to Parrinello-Rahman. As for Nose-Hoover temperature coupling the time constant [`tau-p`](#tau-p) is the period of pressure fluctuations at equilibrium. This is probably a better method when you want to apply pressure scaling during data collection, but beware that you can get very large oscillations if you are starting from a different pressure. Currently only supports isotropic scaling.

* <a name="pcoupltype">`pcoupltype`</a>  
    
    + `isotropic`    
    Isotropic pressure coupling with time constant [`tau-p`](#tau-p). The compressibility and reference pressure are set with [`compressibility`](#compressibility) and [`ref-p`](#ref-p), one value is needed.

    + `semiisotropic`    
    Pressure coupling which is isotropic in the `x` and `y` direction, but different in the `z` direction. This can be useful for membrane simulations. 2 values are needed for `x/y` and `z` directions respectively.

    + `anisotropic`    
    Idem, but 6 values are needed for `xx`, `yy`, `zz`, `xy/yx`, `xz/zx` and `yz/zy` components, respectively. When the off-diagonal compressibilities are set to zero, a rectangular box will stay rectangular. Beware that anisotropic scaling can lead to extreme deformation of the simulation box.

    + `surface-tension`    
    Surface tension coupling for surfaces parallel to the xy-plane. Uses normal pressure coupling for the `z`-direction, while the surface tension is coupled to the `x/y` dimensions of the box. The first [`ref-p`](#ref-p) value is the reference surface tension times the number of surfaces [bar nm], the second value is the reference `z`-pressure [bar]. The two [`compressibility`](#compressibility) values are the compressibility in the `x/y` and `z` direction respectively. The value for the `z`-compressibility should be reasonably accurate since it influences the convergence of the surface-tension, it can also be set to zero to have a box with constant height.

* <a name="nstpcouple">`nstpcouple`</a> (-1)  
    The frequency for coupling the pressure. The default value of -1 sets [`nstpcouple`](#nstpcouple) equal to [`nstlist`](#nstlist), unless `nstlist <= 0`, then a value of 10 is used. For velocity Verlet integrators [`nstpcouple`](#nstpcouple) is set to 1.

* <a name="tau-p">`tau-p`</a> (1) [ps]  
    time constant for coupling

* <a name="compressibility">`compressibility`</a> [bar-1]  
    compressibility (NOTE: this is now really in bar-1) For water at 1 atm and 300 K the compressibility is 4.5e-5 [bar-1].

* <a name="ref-p">`ref-p`</a> [bar]  
    reference pressure for coupling

* <a name="refcoord-scaling">`refcoord-scaling`</a>  
    
    + `no`    
    The reference coordinates for position restraints are not modified. Note that with this option the virial and pressure will depend on the absolute positions of the reference coordinates.

    + `all`    
    The reference coordinates are scaled with the scaling matrix of the pressure coupling.

    + `com`    
    Scale the center of mass of the reference coordinates with the scaling matrix of the pressure coupling. The vectors of each reference coordinate to the center of mass are not scaled. Only one COM is used, even when there are multiple molecules with position restraints. For calculating the COM of the reference coordinates in the starting configuration, periodic boundary conditions are not taken into account. 

* * *

### Simulated annealing

Simulated annealing is controlled separately for each temperature group in
GROMACS. The reference temperature is a piecewise linear function, but you can
use an arbitrary number of points for each group, and choose either a single
sequence or a periodic behaviour for each group. The actual annealing is
performed by dynamically changing the reference temperature used in the
thermostat algorithm selected, so remember that the system will usually not
instantaneously reach the reference temperature!

* <a name="annealing">`annealing`</a>  
    Type of annealing for each temperature group
    
    + `no`    
    No simulated annealing - just couple to reference temperature value.

    + `single`    
    A single sequence of annealing points. If your simulation is longer than the time of the last point, the temperature will be coupled to this constant value after the annealing sequence has reached the last time point.

    + `periodic`    
    The annealing will start over at the first reference point once the last reference time is reached. This is repeated until the simulation ends. 

* <a name="annealing-npoints">`annealing-npoints`</a>  
    A list with the number of annealing reference/control points used for each temperature group. Use 0 for groups that are not annealed. The number of entries should equal the number of temperature groups.

* <a name="annealing-time">`annealing-time`</a>  
    List of times at the annealing reference/control points for each group. If you are using periodic annealing, the times will be used modulo the last value, _i.e._ if the values are 0, 5, 10, and 15, the coupling will restart at the 0ps value after 15ps, 30ps, 45ps, etc. The number of entries should equal the sum of the numbers given in [`annealing-npoints`](#annealing-npoints).

* <a name="annealing-temp">`annealing-temp`</a>  
    List of temperatures at the annealing reference/control points for each group. The number of entries should equal the sum of the numbers given in [`annealing-npoints`](#annealing-npoints).
  
Confused? OK, let's use an example. Assume you have two temperature groups,
set the group selections to `annealing = single periodic`, the number of
points of each group to `annealing-npoints = 3 4`, the times to `annealing-
time = 0 3 6 0 2 4 6` and finally temperatures to `annealing-temp = 298 280
270 298 320 320 298`. The first group will be coupled to 298K at 0ps, but the
reference temperature will drop linearly to reach 280K at 3ps, and then
linearly between 280K and 270K from 3ps to 6ps. After this is stays constant,
at 270K. The second group is coupled to 298K at 0ps, it increases linearly to
320K at 2ps, where it stays constant until 4ps. Between 4ps and 6ps it
decreases to 298K, and then it starts over with the same pattern again, _i.e._
rising linearly from 298K to 320K between 6ps and 8ps. Check the summary
printed by [gmx grompp] if you are unsure!  

* * *

### Velocity generation

* <a name="gen-vel">`gen-vel`</a>  
    
    + `no`    
     Do not generate velocities. The velocities are set to zero when there are no velocities in the input structure file.

    + `yes`    
    Generate velocities in [gmx grompp] according to a Maxwell distribution at temperature [`gen-temp`](#gen-temp), with random seed [`gen-seed`](#gen-seed). This is only meaningful with integrator [`md`](#md).

* <a name="gen-temp">`gen-temp`</a> (300) [K]  
    temperature for Maxwell distribution

* <a name="gen-seed">`gen-seed`</a> (-1) [integer]  
    used to initialize random generator for random velocities, when [`gen-seed`](#gen-seed) is set to -1, a pseudo random seed is used.    

* * *

### Bonds

* <a name="constraints">`constraints`</a>  
    
    + `none`    
    No constraints except for those defined explicitly in the topology, _i.e._ bonds are represented by a harmonic (or other) potential or a Morse potential (depending on the setting of [`morse`](#morse)) and angles by a harmonic (or other) potential.

    + `h-bonds`
    Convert the bonds with H-atoms to constraints.

    + `all-bonds`    
    Convert all bonds to constraints.

    + `h-angles`    
    Convert all bonds and additionally the angles that involve H-atoms to bond-constraints.

    + `all-angles`    
    Convert all bonds and angles to bond-constraints.

* <a name="constraint-algorithm">`constraint-algorithm`</a>  

    + <a name="LINCS">`LINCS`</A>    
    LINear Constraint Solver. With domain decomposition the parallel version P-LINCS is used. The accuracy in set with [`lincs-order`](#lincs-order), which sets the number of matrices in the expansion for the matrix inversion. After the matrix inversion correction the algorithm does an iterative correction to compensate for lengthening due to rotation. The number of such iterations can be controlled with [`lincs-iter`](#lincs-iter). The root mean square relative constraint deviation is printed to the log file every [`nstlog`](#nstlog) steps. If a bond rotates more than [`lincs-warnangle`](#lincs-warnangle) in one step, a warning will be printed both to the log file and to `stderr`. LINCS should not be used with coupled angle constraints.

    + <a name="SHAKE">`SHAKE`</A>    
    SHAKE is slightly slower and less stable than LINCS, but does work with angle constraints. The relative tolerance is set with [`shake-tol`](#shake-tol), 0.0001 is a good value for "normal" MD. SHAKE does not support constraints between atoms on different nodes, thus it can not be used with domain decompositon when inter charge-group constraints are present. SHAKE can not be used with energy minimization.

* <a name="continuation">`continuation`</a>  
    This option was formerly known as `unconstrained-start`.

    + `no`    
    apply constraints to the start configuration and reset shells

    + `yes`    
    do not apply constraints to the start configuration and do not reset shells, useful for exact coninuation and reruns

* <a name="shake-tol">`shake-tol`</a> (0.0001)  
    relative tolerance for SHAKE

* <a name="lincs-order">`lincs-order`</a> (4)  
    Highest order in the expansion of the constraint coupling matrix. When constraints form triangles, an additional expansion of the same order is applied on top of the normal expansion only for the couplings within such triangles. For "normal" MD simulations an order of 4 usually suffices, 6 is needed for large time-steps with virtual sites or BD. For accurate energy minimization an order of 8 or more might be required. With domain decomposition, the cell size is limited by the distance spanned by `lincs-order + 1` constraints. When one wants to scale further than this limit, one can decrease [`lincs-order`](#lincs-order) and increase [`lincs-iter`](#lincs-iter), since the accuracy does not deteriorate when `(1 + lincs-iter) * lincs-order` remains constant.

* <a name="lincs-iter">`lincs-iter`</a> (1)  
    Number of iterations to correct for rotational lengthening in LINCS. For normal runs a single step is sufficient, but for NVE runs where you want to conserve energy accurately or for accurate energy minimization you might want to increase it to 2.

* <a name="lincs-warnangle">`lincs-warnangle`</a> (30) [degrees]  
    maximum angle that a bond can rotate before LINCS will complain

* <a name="morse">`morse`</a>  
    
    + `no`    
    bonds are represented by a harmonic potential

    + `yes`    
    bonds are represented by a Morse potential
  
* * *

### Energy group exclusions

* <a name="energygrp-excl">`energygrp-excl`</a>   
    Pairs of energy groups for which all non-bonded interactions are excluded. An example: if you have two energy groups `Protein` and `SOL`, specifying   
`energygrp-excl = Protein Protein  SOL SOL`  
would give only the non-bonded interactions between the protein and the
solvent. This is especially useful for speeding up energy calculations with
`mdrun -rerun` and for excluding interactions within frozen groups.

  

* * *

### Walls

* <a name="nwall">`nwall`</a> 0  
    When set to 1, there is a wall at `z=0`, when set to 2, there is also a wall at `z=z-box`. Walls can only be used with [`pbc`](#pbc) set to `xy`. When set to 2, pressure coupling and Ewald summation can be used (it is usually best to use semiisotropic pressure coupling with the `x/y` compressibility set to 0, as otherwise the surface area will change). Walls interact wit the rest of the system through an optional [`wall-atomtype`](#wall-atomtype). Energy groups `wall0` and `wall1` (for `nwall = 2`) are added automatically to monitor the interaction of energy groups with each wall. The center of mass motion removal will be turned off in the `z`-direction.

* <a name="wall-atomtype">`wall-atomtype`</a>  
    the atom type name in the force field for each wall. By (for example) defining a special wall atom type in the topology with its own combination rules, this allows for independent tuning of the interaction of each atomtype with the walls.

* <a name="wall-type">`wall-type`</a>  
    
    + `9-3`    
    LJ integrated over the volume behind the wall: 9-3 potential

    + `10-4`    
    LJ integrated over the wall surface: 10-4 potential

    + `12-6`    
    direct LJ potential with the `z` distance from the wall

    + `table`    
    user defined potentials indexed with the `z` distance from the wall, the tables are read analogously to the [`energygrp-table`](#energygrp-table) option, where the first name is for a "normal" energy group and the second name is `wall0` or `wall1`, only the dispersion and repulsion columns are used

* <a name="wall-r-linpot">`wall-r-linpot`</a> -1 (nm)  
    Below this distance from the wall the potential is continued linearly and thus the force is constant. Setting this option to a postive value is especially useful for equilibration when some atoms are beyond a wall. When the value is <=0 (<0 for [`wall-type=table`](#wall-type=table)), a fatal error is generated when atoms are beyond a wall. 

* <a name="wall-density">`wall-density`</a> [nm-3/nm-2]  
    the number density of the atoms for each wall for wall types `9-3` and `10-4` 

* <a name="wall-ewald-zfac">`wall-ewald-zfac`</a> 3  
    The scaling factor for the third box vector for Ewald summation only, the minimum is 2. Ewald summation can only be used with `nwall=2`, where one should use [`ewald-geometry`](#ewald-geometry) `3dc`. The empty layer in the box serves to decrease the unphysical Coulomb interaction between periodic images.
  

* * *

### COM pulling

* <a name="pull">`pull`</a>  
    
    + `no`    
    No center of mass pulling. All the following pull options will be ignored (and if present in the [.mdp] file, they unfortunately generate warnings)

    + `umbrella`    
    Center of mass pulling using an umbrella potential between the reference group and one or more groups.

    + `constraint`    
    Center of mass pulling using a constraint between the reference group and one or more groups. The setup is identical to the option `umbrella`, except for the fact that a rigid constraint is applied instead of a harmonic potential.

    + `constant-force`    
    Center of mass pulling using a linear potential and therefore a constant force. For this option there is no reference position and therefore the parameters [`pull-init`](#pull-init) and [`pull-rate`](#pull-rate) are not used.

* <a name="pull-geometry">`pull-geometry`</a>  
    
    + `distance`    
    Pull along the vector connecting the two groups. Components can be selected with [`pull-dim`](#pull-dim).

    + `direction`    
    Pull in the direction of [`pull-vec`](#pull-vec).

    + `direction-periodic`    
    As `direction`, but allows the distance to be larger than half the box size. With this geometry the box should not be dynamic (_e.g._ no pressure scaling) in the pull dimensions and the pull force is not added to virial.

    + `cylinder`    
    Designed for pulling with respect to a layer where the reference COM is given by a local cylindrical part of the reference group. The pulling is in the direction of [`pull-vec`](#pull-vec). From the reference group a cylinder is selected around the axis going through the pull group with direction [`pull-vec`](#pull-vec) using two radii. The radius [`pull-r1`](#pull-r1) gives the radius within which all the relative weights are one, between [`pull-r1`](#pull-r1) and [`pull-r0`](#pull-r0) the weights are switched to zero. Mass weighting is also used. Note that the radii should be smaller than half the box size. For tilted cylinders they should be even smaller than half the box size since the distance of an atom in the reference group from the COM of the pull group has both a radial and an axial component.

* <a name="pull-dim">`pull-dim`</a> (Y Y Y)  
    the distance components to be used with geometry `distance`, and also sets which components are printed to the output files

* <a name="pull-r1">`pull-r1`</a> (1) [nm]  
    the inner radius of the cylinder for geometry `cylinder`

* <a name="pull-r0">`pull-r0`</a> (1) [nm]  
    the outer radius of the cylinder for geometry `cylinder`

* <a name="pull-constr-tol">`pull-constr-tol`</a> (1e-6)  
    the relative constraint tolerance for constraint pulling

* <a name="pull-start">`pull-start`</a>  
    
    + `no`    
    do not modify [`pull-init`](#pull-init)

    + `yes`  
    add the COM distance of the starting conformation to [`pull-init`](#pull-init)

* <a name="pull-print-reference">`pull-print-reference`</a> (10)  
    
    + `no`    
    do not print the COM of the first group in each pull coordinate

    + `yes`    
    print the COM of the first group in each pull coordinate

* <a name="pull-nstxout">`pull-nstxout`</a> (10)  
    frequency for writing out the COMs of all the pull group

* <a name="pull-nstfout">`pull-nstfout`</a> (1)  
    frequency for writing out the force of all the pulled group

* <a name="pull-ngroups">`pull-ngroups`</a> (1)  
    The number of pull groups, not including the absolute reference group, when used. Pull groups can be reused in multiple pull coordinates. Below only the pull options for group 1 are given, further groups simply increase the group index number.

* <a name="pull-ncoords">`pull-ncoords`</a> (1)  
    The number of pull coordinates. Below only the pull options for coordinate 1 are given, further coordinates simply increase the coordinate index number.

* <a name="pull-group1-name">`pull-group1-name`</a>   
    The name of the pull group, is looked up in the index file or in the default groups to obtain the atoms involved.

* <a name="pull-group1-weights">`pull-group1-weights`</a>   
    Optional relative weights which are multiplied with the masses of the atoms to give the total weight for the COM. The number should be 0, meaning all 1, or the number of atoms in the pull group.

* <a name="pull-group1-pbcatom">`pull-group1-pbcatom`</a> (0)  
    The reference atom for the treatment of periodic boundary conditions inside the group (this has no effect on the treatment of the pbc between groups). This option is only important when the diameter of the pull group is larger than half the shortest box vector. For determining the COM, all atoms in the group are put at their periodic image which is closest to [`pull-group1-pbcatom`](#pull-group1-pbcatom). A value of 0 means that the middle atom (number wise) is used. This parameter is not used with geometry `cylinder`. A value of -1 turns on cosine weighting, which is useful for a group of molecules in a periodic system, _e.g._ a water slab (see Engin et al. J. Chem. Phys. B 2010).

* <a name="pull-coord1-groups">`pull-coord1-groups`</a>   
    The two groups indices should be given on which this pull coordinate will operate. The first index can be 0, in which case an absolute reference of [`pull-coord1-origin`](#pull-coord1-origin) is used. With an absolute reference the system is no longer translation invariant and one should think about what to do with the center of mass motion.

* <a name="pull-coord1-origin">`pull-coord1-origin`</a> (0.0 0.0 0.0)  
    The pull reference position for use with an absolute reference.

* <a name="pull-coord1-vec">`pull-coord1-vec`</a> (0.0 0.0 0.0)  
    The pull direction. [gmx grompp] normalizes the vector.

* <a name="pull-coord1-init">`pull-coord1-init`</a> (0.0) [nm]  
    The reference distance at t=0.

* <a name="pull-coord1-rate">`pull-coord1-rate`</a> (0) [nm/ps]  
    The rate of change of the reference position.

* <a name="pull-coord1-k">`pull-coord1-k`</a> (0) [kJ mol-1 nm-2] / [kJ mol-1 nm-1]  
    The force constant. For umbrella pulling this is the harmonic force constant in [kJ mol-1 nm-2]. For constant force pulling this is the force constant of the linear potential, and thus the negative (!) of the constant force in [kJ mol-1 nm-1].

* <a name="pull-coord1-kB">`pull-coord1-kB`</a> (pull-k1) [kJ mol-1 nm-2] / [kJ mol-1 nm-1]  
    As [`pull-coord1-k`](#pull-coord1-k), but for state B. This is only used when [`free-energy`](#free-energy) is turned on. The force constant is then `(1 - lambda) * pull-coord1-k + lambda * pull-coord1-kB`.

* * *

### NMR refinement

* <a name="disre">`disre`</a>  
    
    + `no`    
    ignore distance restraint information in topology file

    + `simple`    
    simple (per-molecule) distance restraints.

    + `ensemble`  
    distance restraints over an ensemble of molecules in one simulation box. Normally, one would perform ensemble averaging over multiple subsystems, each in a separate box, using a [multi-simulation](#multi-simulation). Supply `topol0.tpr`, `topol1.tpr`, ... with different coordinates and/or velocities. The environment variable `GMX_DISRE_ENSEMBLE_SIZE` sets the number of systems within each ensemble (usually equal to the `mdrun -multi` value).

* <a name="disre-weighting">`disre-weighting`</a>  
    
    + `equal`     (default)
    divide the restraint force equally over all atom pairs in the restraint

    + `conservative`    
    the forces are the derivative of the restraint potential, this results in an weighting of the atom pairs to the reciprocal seventh power of the displacement. The forces are conservative when [`disre-tau`](#disre-tau) is zero.

* <a name="disre-mixed">`disre-mixed`</a>  

    + `no`    
    the violation used in the calculation of the restraint force is the time-averaged violation 

    + `yes`    
    the violation used in the calculation of the restraint force is the square root of the product of the time-averaged violation and the instantaneous violation

* <a name="disre-fc">`disre-fc`</a> (1000) [kJ mol-1 nm-2]  
    force constant for distance restraints, which is multiplied by a (possibly) different factor for each restraint given in the `fac` column of the interaction in the topology file.

* <a name="disre-tau">`disre-tau`</a> (0) [ps]  
    time constant for distance restraints running average. A value of zero turns off time averaging.

* <a name="nstdisreout">`nstdisreout`</a> (100) [steps]  
    period between steps when the running time-averaged and instantaneous distances of all atom pairs involved in restraints are written to the energy file (can make the energy file very large)

* <a name="orire">`orire`</a>  

    + `no`    
    ignore orientation restraint information in topology file

    + `yes`    
    use orientation restraints, ensemble averaging can be performed with a [multi-simulation](#multi-simulation)

* <a name="orire-fc">`orire-fc`</a> (0) [kJ mol]  
    force constant for orientation restraints, which is multiplied by a (possibly) different weight factor for each restraint, can be set to zero to obtain the orientations from a free simulation

* <a name="orire-tau">`orire-tau`</a> (0) [ps]  
    time constant for orientation restraints running average. A value of zero turns off time averaging.

* <a name="orire-fitgrp">`orire-fitgrp`</a>   
    fit group for orientation restraining. This group of atoms is used to determine the rotation `R` of the system with respect to the reference orientation. The reference orientation is the starting conformation of the first subsystem. For a protein, backbone is a reasonable choice

* <a name="nstorireout">`nstorireout`</a> (100) [steps]  
    period between steps when the running time-averaged and instantaneous orientations for all restraints, and the molecular order tensor are written to the energy file (can make the energy file very large)
  

* * *

### Free energy calculations

* <a name="free-energy">`free-energy`</a>  
    
    + `no`    
    Only use topology A.

    + `yes`    
    Interpolate between topology A (lambda=0) to topology B (lambda=1) and write the derivative of the Hamiltonian with respect to lambda (as specified with [`dhdl-derivatives`](#dhdl-derivatives)), or the Hamiltonian differences with respect to other lambda values (as specified with [`foreign-lambda`](#foreign-lambda)) to the energy file and/or to `dhdl.xvg`, where they can be processed by, for example [gmx bar]. The potentials, bond-lengths and angles are interpolated linearly as described in the manual. When [`sc-alpha`](#sc-alpha) is larger than zero, soft-core potentials are used for the LJ and Coulomb interactions.

    + `expanded`    
     Turns on expanded ensemble simulation, where the alchemical state becomes a dynamic variable, allowing jumping between different Hamiltonians. See the expanded ensemble options for controlling how expanded ensemble simulations are performed. The different Hamiltonians used in expanded ensemble simulations are defined by the other free energy options.

* <a name="init-lambda">`init-lambda`</a> (-1)  
    starting value for lambda (float). Generally, this should only be used with slow growth (_i.e._ nonzero [`delta-lambda`](#delta-lambda)). In other cases, [`init-lambda-state`](#init-lambda-state) should be specified instead. Must be greater than or equal to 0.

* <a name="delta-lambda">`delta-lambda`</a> (0)  
    increment per time step for lambda

* <a name="init-lambda-state">`init-lambda-state`</a> (-1)  
    starting value for the lambda state (integer). Specifies which columm of the lambda vector ([`coul-lambdas`](#coul-lambdas), [`vdw-lambdas`](#vdw-lambdas), [`bonded-lambdas`](#bonded-lambdas), [`restraint-lambdas`](#restraint-lambdas), [`mass-lambdas`](#mass-lambdas), [`temperature-lambdas`](#temperature-lambdas), [`fep-lambdas`](#fep-lambdas)) should be used. This is a zero-based index: [`init-lambda-state`](#init-lambda-state) 0 means the first column, and so on.

* <a name="fep-lambdas">`fep-lambdas`</a> ()  
    Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every [`nstdhdl`](#nstdhdl) steps. Values must be between 0 and 1. Free energy differences between different lambda values can then be determined with [gmx bar]. [`fep-lambdas`](#fep-lambdas) is different from the other `-lambdas` keywords because all components of the lambda vector that are not specified will use [`fep-lambdas`](#fep-lambdas) (including restraint-lambdas and therefore the pull code restraints).

* <a name="coul-lambdas">`coul-lambdas`</a> ()  
    Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every [`nstdhdl`](#nstdhdl) steps. Values must be between 0 and 1. Only the electrostatic interactions are controlled with this component of the lambda vector (and only if the lambda=0 and lambda=1 states have differing electrostatic interactions).

* <a name="vdw-lambdas">`vdw-lambdas`</a> ()  
    Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every [`nstdhdl`](#nstdhdl) steps. Values must be between 0 and 1. Only the van der Waals interactions are controlled with this component of the lambda vector.

* <a name="bonded-lambdas">`bonded-lambdas`</a> ()  
    Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every [`nstdhdl`](#nstdhdl) steps. Values must be between 0 and 1. Only the bonded interactions are controlled with this component of the lambda vector.

* <a name="restraint-lambdas">`restraint-lambdas`</a> ()  
    Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every [`nstdhdl`](#nstdhdl) steps. Values must be between 0 and 1. Only the restraint interactions: dihedral restraints, and the pull code restraints are controlled with this component of the lambda vector. 

* <a name="mass-lambdas">`mass-lambdas`</a> ()  
    Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every [`nstdhdl`](#nstdhdl) steps. Values must be between 0 and 1. Only the particle masses are controlled with this component of the lambda vector.

* <a name="temperature-lambdas">`temperature-lambdas`</a> ()  
    Zero, one or more lambda values for which Delta H values will be determined and written to dhdl.xvg every [`nstdhdl`](#nstdhdl) steps. Values must be between 0 and 1. Only the temperatures controlled with this component of the lambda vector. Note that these lambdas should not be used for replica exchange, only for simulated tempering.

* <a name="calc-lambda-neighbors">`calc-lambda-neighbors`</a> (1)  
    Controls the number of lambda values for which Delta H values will be calculated and written out, if [`init-lambda-state`](#init-lambda-state) has been set. A positive value will limit the number of lambda points calculated to only the nth neighbors of [`init-lambda-state`](#init-lambda-state): for example, if [`init-lambda-state`](#init-lambda-state) is 5 and this parameter has a value of 2, energies for lambda points 3-7 will be calculated and writen out. A value of -1 means all lambda points will be written out. For normal BAR such as with [gmx bar], a value of 1 is sufficient, while for MBAR -1 should be used.

* <a name="sc-alpha">`sc-alpha`</a> (0)  
    the soft-core alpha parameter, a value of 0 results in linear interpolation of the LJ and Coulomb interactions

* <a name="sc-r-power">`sc-r-power`</a> (6)  
    the power of the radial term in the soft-core equation. Possible values are 6 and 48. 6 is more standard, and is the default. When 48 is used, then sc-alpha should generally be much lower (between 0.001 and 0.003).

* <a name="sc-coul">`sc-coul`</a> (no)  
    Whether to apply the soft core free energy interaction transformation to the Columbic interaction of a molecule. Default is no, as it is generally more efficient to turn off the Coulomic interactions linearly before turning off the van der Waals interactions.

* <a name="sc-power">`sc-power`</a> (0)  
    the power for lambda in the soft-core function, only the values 1 and 2 are supported

* <a name="sc-sigma">`sc-sigma`</a> (0.3) [nm]  
    the soft-core sigma for particles which have a C6 or C12 parameter equal to zero or a sigma smaller than [`sc-sigma`](#sc-sigma)

* <a name="couple-moltype">`couple-moltype`</a>  
    Here one can supply a molecule type (as defined in the topology) for calculating solvation or coupling free energies. There is a special option `system` that couples all molecule types in the system. This can be useful for equilibrating a system starting from (nearly) random coordinates. [`free-energy`](#free-energy) has to be turned on. The Van der Waals interactions and/or charges in this molecule type can be turned on or off between lambda=0 and lambda=1, depending on the settings of [`couple-lambda0`](#couple-lambda0) and [`couple-lambda1`](#couple-lambda1). If you want to decouple one of several copies of a molecule, you need to copy and rename the molecule definition in the topology.

* <a name="couple-lambda0">`couple-lambda0`</a>  

    + `vdw-q`    
    all interactions are on at lambda=0

    + `vdw`  
    the charges are zero (no Coulomb interactions) at lambda=0

    + `q`  
    the Van der Waals interactions are turned at lambda=0; soft-core interactions will be required to avoid singularities

    + `none`  
    the Van der Waals interactions are turned off and the charges are zero at lambda=0; soft-core interactions will be required to avoid singularities.

* <a name="couple-lambda1">`couple-lambda1`</a>  
     analogous to [`couple-lambda0`](#couple-lambda0), but for lambda=1

* <a name="couple-intramol">`couple-intramol`</a>  

    + `no`   
    All intra-molecular non-bonded interactions for moleculetype [`couple-moltype`](#couple-moltype) are replaced by exclusions and explicit pair interactions. In this manner the decoupled state of the molecule corresponds to the proper vacuum state without periodicity effects.

    + `yes`  
    The intra-molecular Van der Waals and Coulomb interactions are also turned on/off. This can be useful for partitioning free-energies of relatively large molecules, where the intra-molecular non-bonded interactions might lead to kinetically trapped vacuum conformations. The 1-4 pair interactions are not turned off.

* <a name="nstdhdl">`nstdhdl`</a> (100)  
    the frequency for writing dH/dlambda and possibly Delta H to dhdl.xvg, 0 means no ouput, should be a multiple of [`nstcalcenergy`](#nstcalcenergy).

* <a name="dhdl-derivatives">`dhdl-derivatives`</a> (yes)  
    If yes (the default), the derivatives of the Hamiltonian with respect to lambda at each [`nstdhdl`](#nstdhdl) step are written out. These values are needed for interpolation of linear energy differences with [gmx bar] (although the same can also be achieved with the right [`foreign-lambda`](#foreign-lambda) setting, that may not be as flexible), or with thermodynamic integration

* <a name="dhdl-print-energy">`dhdl-print-energy`</a> (no)  
     Include either the total or the potential energy in the dhdl file. Options are 'no', 'potential', or 'total'. This information is needed for later free energy analysis if the states of interest are at different temperatures. If all states are at the same temperature, this information is not needed. 'potential' is useful in case one is using a [rerun](#rerun) to generate the `dhdl.xvg` file. When rerunning from an existing trajectory, the kinetic energy will often not be correct, and thus one must compute the residual free energy from the potential alone, with the kinetic energy component computed analytically.

* <a name="separate-dhdl-file">`separate-dhdl-file`</a> (yes)  
    
    + `yes`    
    the free energy values that are calculated (as specified with the [`foreign-lambda`](#foreign-lambda) and [`dhdl-derivatives`](#dhdl-derivatives) settings) are written out to a separate file, with the default name `dhdl.xvg`. This file can be used directly with [gmx bar].

    + `no`    
    The free energy values are written out to the energy output file (`ener.edr`, in accumulated blocks at every [`nstenergy`](#nstenergy) steps), where they can be extracted with [gmx energy] or used directly with [gmx bar].

* <a name="dh-hist-size">`dh-hist-size`</a> (0)  
    If nonzero, specifies the size of the histogram into which the Delta H values (specified with [`foreign-lambda`](#foreign-lambda)) and the derivative dH/dl values are binned, and written to ener.edr. This can be used to save disk space while calculating free energy differences. One histogram gets written for each [`foreign-lambda`](#foreign-lambda) and two for the dH/dl, at every [`nstenergy`](#nstenergy) step. Be aware that incorrect histogram settings (too small size or too wide bins) can introduce errors. Do not use histograms unless you're certain you need it.

* <a name="dh-hist-spacing">`dh-hist-spacing`</a> (0.1)  
    Specifies the bin width of the histograms, in energy units. Used in conjunction with [`dh-hist-size`](#dh-hist-size). This size limits the accuracy with which free energies can be calculated. Do not use histograms unless you're certain you need it.
  
* * *

### Expanded Ensemble calculations

* <a name="nstexpanded">`nstexpanded`</a>  
    The number of integration steps beween attempted moves changing the system Hamiltonian in expanded ensemble simulations. Must be a multiple of [`nstcalcenergy`](#nstcalcenergy), but can be greater or less than [`nstdhdl`](#nstdhdl).
* <a name="lmc-stats">`lmc-stats`</a>  
    
    + `no`    
    No Monte Carlo in state space is performed.

    + `metropolis-transition`    
     Uses the Metropolis weights to update the expanded ensemble weight of each state. Min{1,exp(-(beta_new u_new - beta_old u_old)}

    + `barker-transition`    
     Uses the Barker transition critera to update the expanded ensemble weight of each state i, defined by exp(-beta_new u_new)/[exp(-beta_new u_new)+exp(-beta_old u_old)

    + `wang-landau`    
    Uses the Wang-Landau algorithm (in state space, not energy space) to update the expanded ensemble weights.

    + `min-variance`    
    Uses the minimum variance updating method of Escobedo et al. to update the expanded ensemble weights. Weights will not be the free energies, but will rather emphasize states that need more sampling to give even uncertainty.

* <a name="lmc-mc-move">`lmc-mc-move`</a>  
    
    + `no`    
    No Monte Carlo in state space is performed.

    + `metropolis-transition`    
     Randomly chooses a new state up or down, then uses the Metropolis critera to decide whether to accept or reject: Min{1,exp(-(beta_new u_new - beta_old u_old)}

    + `barker-transition`    
     Randomly chooses a new state up or down, then uses the Barker transition critera to decide whether to accept or reject: exp(-beta_new u_new)/[exp(-beta_new u_new)+exp(-beta_old u_old)]

    + `gibbs`    
     Uses the conditional weights of the state given the coordinate (exp(-beta_i u_i) / sum_k exp(beta_i u_i) to decide which state to move to.

    + `metropolized-gibbs`    
          Uses the conditional weights of the state given the coordinate (exp(-beta_i u_i) / sum_k exp(beta_i u_i) to decide which state to move to, EXCLUDING the current state, then uses a rejection step to ensure detailed balance. Always more efficient that Gibbs, though only marginally so in many situations, such as when only the nearest neighbors have decent phase space overlap.

* <a name="lmc-seed">`lmc-seed`</a> (-1)  
     random seed to use for Monte Carlo moves in state space. When [`lmc-seed`](#lmc-seed) is set to -1, a pseudo random seed is us

* <a name="mc-temperature">`mc-temperature`</a>  
     Temperature used for acceptance/rejection for Monte Carlo moves. If not specified, the temperature of the simulation specified in the first group of [`ref-t`](#ref-t) is used.

* <a name="wl-ratio">`wl-ratio`</a> (0.8)  
    The cutoff for the histogram of state occupancies to be reset, and the free energy incrementor to be changed from delta to delta * [`wl-scale`](#wl-scale). If we define the Nratio = (number of samples at each histogram) / (average number of samples at each histogram). [`wl-ratio`](#wl-ratio) of 0.8 means that means that the histogram is only considered flat if all Nratio > 0.8 AND simultaneously all 1/Nratio > 0.8.

* <a name="wl-scale">`wl-scale`</a> (0.8)  
     Each time the histogram is considered flat, then the current value of the Wang-Landau incrementor for the free energies is multiplied by [`wl-scale`](#wl-scale). Value must be between 0 and 1.

* <a name="init-wl-delta">`init-wl-delta`</a> (1.0)  
    The initial value of the Wang-Landau incrementor in kT. Some value near 1 kT is usually most efficient, though sometimes a value of 2-3 in units of kT works better if the free energy differences are large.

* <a name="wl-oneovert">`wl-oneovert`</a> (no)  
    Set Wang-Landau incrementor to scale with 1/(simulation time) in the large sample limit. There is significant evidence that the standard Wang-Landau algorithms in state space presented here result in free energies getting 'burned in' to incorrect values that depend on the initial state. when [`wl-oneovert`](#wl-oneovert) is true, then when the incrementor becomes less than 1/N, where N is the mumber of samples collected (and thus proportional to the data collection time, hence '1 over t'), then the Wang-Lambda incrementor is set to 1/N, decreasing every step. Once this occurs, [`wl-ratio`](#wl-ratio) is ignored, but the weights will still stop updating when the equilibration criteria set in [`lmc-weights-equil`](#lmc-weights-equil) is achieved.

* <a name="lmc-repeats">`lmc-repeats`</a> (1)  
    Controls the number of times that each Monte Carlo swap type is performed each iteration. In the limit of large numbers of Monte Carlo repeats, then all methods converge to Gibbs sampling. The value will generally not need to be different from 1.

* <a name="lmc-gibbsdelta">`lmc-gibbsdelta`</a> (-1)  
     Limit Gibbs sampling to selected numbers of neighboring states. For Gibbs sampling, it is sometimes inefficient to perform Gibbs sampling over all of the states that are defined. A positive value of [`lmc-gibbsdelta`](#lmc-gibbsdelta) means that only states plus or minus [`lmc-gibbsdelta`](#lmc-gibbsdelta) are considered in exchanges up and down. A value of -1 means that all states are considered. For less than 100 states, it is probably not that expensive to include all states.

* <a name="lmc-forced-nstart">`lmc-forced-nstart`</a> (0)  
     Force initial state space sampling to generate weights. In order to come up with reasonable initial weights, this setting allows the simulation to drive from the initial to the final lambda state, with [`lmc-forced-nstart`](#lmc-forced-nstart) steps at each state before moving on to the next lambda state. If [`lmc-forced-nstart`](#lmc-forced-nstart) is sufficiently long (thousands of steps, perhaps), then the weights will be close to correct. However, in most cases, it is probably better to simply run the standard weight equilibration algorithms.

* <a name="nst-transition-matrix">`nst-transition-matrix`</a> (-
    Frequency of outputting the expanded ensemble transition matrix. A negative number means it will only be printed at the end of the simulation.

* <a name="symmetrized-transition-matrix">`symmetrized-transition-matrix`</a> (no)   
    Whether to symmetrize the empirical transition matrix. In the infinite limit the matrix will be symmetric, but will diverge with statistical noise for short timescales. Forced symmetrization, by using the matrix T_sym = 1/2 (T + transpose(T)), removes problems like the existence of (small magnitude) negative eigenvalues.

* <a name="mininum-var-min">`mininum-var-min`</a> (100)  
     The `min-variance` strategy (option of [`lmc-stats`](#lmc-stats) is only valid for larger number of samples, and can get stuck if too few samples are used at each state. [`mininum-var-min`](#mininum-var-min) is the minimum number of samples that each state that are allowed before the `min-variance` strategy is activated if selected.

* <a name="init-lambda-weights">`init-lambda-weights`</a>   
    The initial weights (free energies) used for the expanded ensemble states. Default is a vector of zero weights. format is similar to the lambda vector settings in [`fep-lambdas`](#fep-lambdas), except the weights can be any floating point number. Units are kT. Its length must match the lambda vector lengths.

* <a name="lmc-weights-equil">`lmc-weights-equil`</a> (no)  

    + `no`    
    Expanded ensemble weights continue to be updated throughout the simulation.

    + `yes`    
    The input expanded ensemble weights are treated as equilibrated, and are not updated throughout the simulation.

    + `wl-delta`    
    Expanded ensemble weight updating is stopped when the Wang-Landau incrementor falls below the value specified by [`weight-equil-wl-delta`](#weight-equil-wl-delta).

    + `number-all-lambda`    
    Expanded ensemble weight updating is stopped when the number of samples at all of the lambda states is greater than the value specified by [`weight-equil-number-all-lambda`](#weight-equil-number-all-lambda).

    + `number-steps`    
    Expanded ensemble weight updating is stopped when the number of steps is greater than the level specified by [`weight-equil-number-steps`](#weight-equil-number-steps).

    + `number-samples`    
    Expanded ensemble weight updating is stopped when the number of total samples across all lambda states is greater than the level specified by [`weight-equil-number-samples`](#weight-equil-number-samples).

    + `count-ratio`    
    Expanded ensemble weight updating is stopped when the ratio of samples at the least sampled lambda state and most sampled lambda state greater than the value specified by [`weight-equil-count-ratio`](#weight-equil-count-ratio).

* <a name="simulated-tempering">`simulated-tempering`</a> (no)  
    Turn simulated tempering on or off. Simulated tempering is implemented as expanded ensemble sampling with different temperatures instead of different Hamiltonians.

* <a name="sim-temp-low">`sim-temp-low`</a> (300) [K]  
    Low temperature for simulated tempering.

* <a name="sim-temp-high">`sim-temp-high`</a> (300) [K]  
    High temperature for simulated tempering.

* <a name="simulated-tempering-scaling">`simulated-tempering-scaling`</a> (linear)  
    Controls the way that the temperatures at intermediate lambdas are calculated from the [`temperature-lambda`](#temperature-lambda) part of the lambda vector.

    + `linear`    
    Linearly interpolates the temperatures using the values of [`temperature-lambda`](#temperature-lambda), _i.e._ if `sim-temp-low=300, sim-temp-high=400`, then `lambda=0.5` correspond to a temperature of 350. A nonlinear set of temperatures can always be implemented with uneven spacing in lambda.

    + `geometric`    
     Interpolates temperatures geometrically between [`sim-temp-low`](#sim-temp-low) and [`sim-temp-high`](#sim-temp-high). The i:th state has temperature `sim-temp-low * (sim-temp-high / sim-temp-low)` raised to the power of `(i/(ntemps-1))`. This should give roughly equal exchange for constant heat capacity, though of course things simulations that involve protein folding have very high heat capacity peaks.

    + `exponential`    
     Interpolates temperatures exponentially between [`sim-temp-low`](#sim-temp-low) and [`sim-temp-high`](#sim-temp-high). The ith state has temperature `sim-temp-low + (sim-temp-high - sim-temp-low) * ((exp(temperature-lambdas[i]) - 1) / (exp(1.0) - 1))`.
  

* * *

### Non-equilibrium MD


* <a name="acc-grps">`acc-grps`</a>   
    groups for constant acceleration (_e.g._ `Protein Sol`) all atoms in groups Protein and Sol will experience constant acceleration as specified in the [`accelerate`](#accelerate) line

* <a name="accelerate">`accelerate`</a> (0) [nm ps-2]  
    acceleration for [`acc-grps`](#acc-grps); x, y and z for each group (_e.g._ `0.1 0.0 0.0 -0.1 0.0 0.0` means that first group has constant acceleration of 0.1 nm ps-2 in X direction, second group the opposite).

* <a name="freezegrps">`freezegrps`</a>   
    Groups that are to be frozen (_i.e._ their X, Y, and/or Z position will not be updated; _e.g._ `Lipid SOL`). [`freezedim`](#freezedim) specifies for which dimension the freezing applies. To avoid spurious contibrutions to the virial and pressure due to large forces between completely frozen atoms you need to use energy group exclusions, this also saves computing time. Note that coordinates of frozen atoms are not scaled by pressure-coupling algorithms.

* <a name="freezedim">`freezedim`</a>   
    dimensions for which groups in [`freezegrps`](#freezegrps) should be frozen, specify `Y` or `N` for X, Y and Z and for each group (_e.g._ `Y Y N N N N` means that particles in the first group can move only in Z direction. The particles in the second group can move in any direction).

* <a name="cos-acceleration">`cos-acceleration`</a> (0) [nm ps-2]  
    the amplitude of the acceleration profile for calculating the viscosity. The acceleration is in the X-direction and the magnitude is [`cos-acceleration`](#cos-acceleration) `cos(2 pi z/boxheight)`. Two terms are added to the energy file: the amplitude of the velocity profile and 1/viscosity.

* <a name="deform">`deform`</a> (0 0 0 0 0 0) [nm ps-1]  
    The velocities of deformation for the box elements: a(x) b(y) c(z) b(x) c(x) c(y). Each step the box elements for which [`deform`](#deform) is non-zero are calculated as: box(ts)+(t-ts)*deform, off-diagonal elements are corrected for periodicity. The coordinates are transformed accordingly. Frozen degrees of freedom are (purposely) also transformed. The time ts is set to t at the first step and at steps at which x and v are written to trajectory to ensure exact restarts. Deformation can be used together with semiisotropic or anisotropic pressure coupling when the appropriate compressibilities are set to zero. The diagonal elements can be used to strain a solid. The off-diagonal elements can be used to shear a solid or a liquid.
  

* * *

### Electric fields

* <a name="E-x">`E-x`</a>  
  <a name="E-y">`E-y`</a>  
  <a name="E-z">`E-z`</a>  
    If you want to use an electric field in a direction, enter 3 numbers after the appropriate `E-`, the first number: the number of cosines, only 1 is implemented (with frequency 0) so enter 1, the second number: the strength of the electric field in [V nm-1], the third number: the phase of the cosine, you can enter any number here since a cosine of frequency zero has no phase.

* <a name="E-xt">`E-xt`</a>  
  <a name="E-yt">`E-yt`</a>  
  <a name="E-zt">`E-zt`</a>  
    not implemented yet
  

* * *

  

### Mixed quantum/classical molecular dynamics

* <a name="QMMM">`QMMM`</A>  

    + `no`    
    No QM/MM.

    + `yes`    
    Do a QM/MM simulation. Several groups can be described at different QM levels separately. These are specified in the [`QMMM-grps`](#QMMM-grps) field separated by spaces. The level of _ab initio_ theory at which the groups are described is specified by [`QMmethod`](#QMmethod) and [`QMbasis`](#QMbasis) Fields. Describing the groups at different levels of theory is only possible with the ONIOM QM/MM scheme, specified by [`QMMMscheme`](#QMMMscheme).

* <a name="QMMM-grps">`QMMM-grps`</a>  
    groups to be descibed at the QM level

* <a name="QMMMscheme">`QMMMscheme`</a>  

    + `normal`    
    normal QM/MM. There can only be one [`QMMM-grps`](#QMMM-grps) that is modelled at the [`QMmethod`](#QMmethod) and [`QMbasis`](#QMbasis) level of _ab initio_ theory. The rest of the system is described at the MM level. The QM and MM subsystems interact as follows: MM point charges are included in the QM one-electron hamiltonian and all Lennard-Jones interactions are described at the MM level.

    + `ONIOM`    
    The interaction between the subsystem is described using the ONIOM method by Morokuma and co-workers. There can be more than one [`QMMM-grps`](#QMMM-grps) each modeled at a different level of QM theory ([`QMmethod`](#QMmethod) and [`QMbasis`](#QMbasis)). 

* <a name="QMmethod">`QMmethod`</a> (RHF)  
    Method used to compute the energy and gradients on the QM atoms. Available methods are AM1, PM3, RHF, UHF, DFT, B3LYP, MP2, CASSCF, and MMVB. For CASSCF, the number of electrons and orbitals included in the active space is specified by [`CASelectrons`](#CASelectrons) and [`CASorbitals`](#CASorbitals). 

* <a name="QMbasis">`QMbasis`</a> (STO-3G)  
    Basis set used to expand the electronic wavefuntion. Only Gaussian basis sets are currently available, _i.e._ `STO-3G, 3-21G, 3-21G*, 3-21+G*, 6-21G, 6-31G, 6-31G*, 6-31+G*,` and `6-311G`.

* <a name="QMcharge">`QMcharge`</a> (0) [integer]  
    The total charge in `e` of the [`QMMM-grps`](#QMMM-grps). In case there are more than one [`QMMM-grps`](#QMMM-grps), the total charge of each ONIOM layer needs to be specified separately.

* <a name="QMmult">`QMmult`</a> (1) [integer]  
    The multiplicity of the [`QMMM-grps`](#QMMM-grps). In case there are more than one [`QMMM-grps`](#QMMM-grps), the multiplicity of each ONIOM layer needs to be specified separately.

* <a name="CASorbitals">`CASorbitals`</a> (0) [integer]  
    The number of orbitals to be included in the active space when doing a CASSCF computation.

* <a name="CASelectrons">`CASelectrons`</a> (0) [integer]  
    The number of electrons to be included in the active space when doing a CASSCF computation.

* <a name="SH">`SH`</A>  

    + `no`    
    No surface hopping. The system is always in the electronic ground-state.

    + `yes`    
    Do a QM/MM MD simulation on the excited state-potential energy surface and enforce a _diabatic_ hop to the ground-state when the system hits the conical intersection hyperline in the course the simulation. This option only works in combination with the CASSCF method.
  

* * *

### Implicit solvent

* <a name="implicit-solvent">`implicit-solvent`</a>  

    + `no`    
    No implicit solvent

    + `GBSA`    
    Do a simulation with implicit solvent using the Generalized Born formalism. Three different methods for calculating the Born radii are available, Still, HCT and OBC. These are specified with the [`gb-algorithm`](#gb-algorithm) field. The non-polar solvation is specified with the [`sa-algorithm`](#sa-algorithm) field.

* <a name="gb-algorithm">`gb-algorithm`</a>  

    + `Still`    
    Use the Still method to calculate the Born radii

    + `HCT`    
    Use the Hawkins-Cramer-Truhlar method to calculate the Born radii

    + `OBC`    
    Use the Onufriev-Bashford-Case method to calculate the Born radii

* <a name="nstgbradii">`nstgbradii`</a> (1) [steps]  
    Frequency to (re)-calculate the Born radii. For most practial purposes, setting a value larger than 1 violates energy conservation and leads to unstable trajectories.

* <a name="rgbradii">`rgbradii`</a> (1.0) [nm]  
    Cut-off for the calculation of the Born radii. Currently must be equal to rlist

* <a name="gb-epsilon-solvent">`gb-epsilon-solvent`</a> (80)  
    Dielectric constant for the implicit solvent

* <a name="gb-saltconc">`gb-saltconc`</a> (0) [M]  
    Salt concentration for implicit solvent models, currently not used

* <a name="gb-obc-alpha">`gb-obc-alpha`</a> (1); gb-obc-beta (0.8); gb-obc-gamma (4.85);  
    Scale factors for the OBC model. Default values are OBC(II). Values for OBC(I) are 0.8, 0 and 2.91 respectively

* <a name="gb-dielectric-offset">`gb-dielectric-offset`</a> (0.009) [nm]  
    Distance for the di-electric offset when calculating the Born radii. This is the offset between the center of each atom the center of the polarization energy for the corresponding atom

* <a name="sa-algorithm">`sa-algorithm`</a>  

    + `Ace-approximation`    
    Use an Ace-type approximation (default)

    + `None`    
    No non-polar solvation calculation done. For GBSA only the polar part gets calculated

* <a name="sa-surface-tension">`sa-surface-tension`</a> [kJ mol-1 nm-2]  
    Default value for surface tension with SA algorithms. The default value is -1; Note that if this default value is not changed it will be overridden by [gmx grompp] using values that are specific for the choice of radii algorithm (0.0049 kcal/mol/Angstrom2 for Still, 0.0054 kcal/mol/Angstrom2 for HCT/OBC) Setting it to 0 will while using an sa-algorithm other than None means no non-polar calculations are done. 
  

* * *

### Adaptive Resolution Simulation

* <a name="adress">`adress`</a> (no)  
    Decide whether the AdResS feature is turned on.

* <a name="adress-type">`adress-type`</a> (Off)  

    + `Off`    
    Do an AdResS simulation with weight equal 1, which is equivalent to an explicit (normal) MD simulation. The difference to disabled AdResS is that the AdResS variables are still read-in and hence are defined.

    + `Constant`    
    Do an AdResS simulation with a constant weight, [`adress-const-wf`](#adress-const-wf) defines the value of the weight

    + `XSplit`    
    Do an AdResS simulation with simulation box split in x-direction, so basically the weight is only a function of the x coordinate and all distances are measured using the x coordinate only.

    + `Sphere`    
    Do an AdResS simulation with spherical explicit zone.

* <a name="adress-const-wf">`adress-const-wf`</a> (1)  
    Provides the weight for a constant weight simulation ([`adress-type`](#adress-type) `Constant`)

* <a name="adress-ex-width">`adress-ex-width`</a> (0)  
    Width of the explicit zone, measured from [`adress-reference-coords`](#adress-reference-coords).

* <a name="adress-hy-width">`adress-hy-width`</a> (0)  
    Width of the hybrid zone.

* <a name="adress-reference-coords">`adress-reference-coords`</a> (0,0,0)  
    Position of the center of the explicit zone. Periodic boundary conditions apply for measuring the distance from it.

* <a name="adress-cg-grp-names">`adress-cg-grp-names`</a>  
    The names of the coarse-grained energy groups. All other energy groups are considered explicit and their interactions will be automatically excluded with the coarse-grained groups.
* <a name="adress-site">`adress-site`</a> (COM)  The mapping point from which the weight is calculated.

    + `COM`    
    The weight is calculated from the center of mass of each charge group.

    + `COG`    
    The weight is calculated from the center of geometry of each charge group.

    + `Atom`    
    The weight is calculated from the position of 1st atom of each charge group.

    + `AtomPerAtom`    
    The weight is calculated from the position of each individual atom.

* <a name="adress-interface-correction">`adress-interface-correction`</a> (Off)  

    + `Off`    
    Do not a apply any interface correction.

    + `thermoforce`    
    Apply thermodynamic force interface correction. The table can be specified using the `-tabletf` option of [mdrun]. The table should contain the potential and force (acting on molecules) as function of the distance from [`adress-reference-coords`](#adress-reference-coords).

* <a name="adress-tf-grp-names">`adress-tf-grp-names`</a>  
    The names of the energy groups to which the `thermoforce` is applied if enabled in [`adress-interface-correction`](#adress-interface-correction). If no group is given the default table is applied.

* <a name="adress-ex-forcecap">`adress-ex-forcecap`</a> (0)  
    Cap the force in the hybrid region, useful for big molecules. 0 disables force capping.
  

* * *

### User defined thingies

* <a name="user1-grps">`user1-grps`</a>  
  <a name="user2-grps">`user2-grps`</a>  
  <a name="userint1">`userint1`</a> (0)  
  <a name="userint2">`userint2`</a> (0)  
  <a name="userint3">`userint3`</a> (0)  
  <a name="userint4">`userint4`</a> (0)  
  <a name="userreal1">`userreal1`</a> (0)  
  <a name="userreal2">`userreal2`</a> (0)  
  <a name="userreal3">`userreal3`</a> (0)  
  <a name="userreal4">`userreal4`</a> (0)  
    These you can use if you modify code. You can pass integers and reals to your subroutine. Check the inputrec definition in `src/include/types/inputrec.h`

