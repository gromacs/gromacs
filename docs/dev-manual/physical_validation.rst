Physical validation
===================

Physical validation tests check whether simulation results correspond
to physical (or mathematical) expectations.

Unlike the existing tests, we are not be able to keep these tests in
the "seconds, not minutes" time frame, rather aiming for "hours, not
days".  They should therefore be ran periodically, but probably not
for every build.

Also, given the long run time, it will in many cases be necessary to
separate running of the systems (e.g. to run it at a specific time, or
on a different resource), such that the make script does give the
option to

* prepare run files and an execution script,
* analyze already present simulations,
* or prepare, run and analyze in one go.


Test description
----------------

Currently, simulation results are tested against three physically /
mathematically expected results:

* *Integrator convergence*: A symplectic integrator can be shown to
  conserve a constant of motion (such as the energy in a
  micro-canonical simulation) up to a fluctuation that is quadratic in
  time step chosen. Comparing two or more constant-of-motion
  trajectories realized using different time steps (but otherwise
  unchanged simulation parameters) allows a check of the symplecticity
  of the integration. Note that lack of symplecticity does not
  necessarily imply an error in the integration algorithm, it can also
  hint at physical violations in other parts of the model, such as
  non-continuous potential functions, imprecise handling of
  constraints, etc.
* *Kinetic energy distribution*: The kinetic energy trajectory of a
  (equilibrated) system sampling a canonical or an isothermal-isobaric
  ensemble is expected to be Maxwell-Boltzmann distributed. The
  similarity between the physically expected and the observed
  distribution allows to validate the sampled kinetic energy ensemble.
* *Distribution of configurational quantities*: As the distribution of
  configurational quantities like the potential energy or the volume
  are in general not known analytically, testing the likelihood of a
  trajectory sampling a given ensemble is less straightforward than
  for the kinetic energy. However, generally, the ratio of the
  probability distribution between samples of the same ensemble at
  different state points (e.g. at different temperatures, different
  pressures) is known. Comparing two simulations at different state
  points therefore allows a validation of the sampled ensemble.

The physical validation included in |Gromacs| tests a range of the
most-used settings on several systems. The general philosophy is to
leave most settings to default values with the exception of the ones
explicitly tested in order to be sensitive to changes in the default
values. The test set will be enlarged as we discover interesting test
systems and corner cases. Under double precision, some additional
tests are ran, and some other tests are ran using a lower tolerance.


Integrator convergence
^^^^^^^^^^^^^^^^^^^^^^

All simulations performed under NVE on Argon (1000 atoms) and water
(900 molecules) systems. As these tests are very sensitive to
numerical imprecision, they are performed with long-range corrections
for both Lennard-Jones and electrostatic interactions, with a very low
pair-list tolerance (``verlet-buffer-tolerance = 1e-10``), and high
LINCS settings where applicable.

**Argon**:

* *Integrators*:
  - ``integrator = md``
  - ``integrator = md-vv``
* *Long-range corrections LJ*:
  - ``vdwtype = PME``
  - ``vdwtype = cut-off``, ``vdw-modifier = force-switch``, ``rvdw-switch = 0.8``

**Water**:

* *Integrators*:
  - ``integrator = md``
  - ``integrator = md-vv``
* *Long-range corrections LJ*:
  - ``vdwtype = PME``
  - ``vdwtype = cut-off``, ``vdw-modifier = force-switch``, ``rvdw-switch = 0.8``
* *Long-range corrections electrostatics*:
  - ``coulombtype = PME``, ``fourierspacing = 0.05``
* *Constraint algorithms*:
  - ``constraint-algorithm = lincs``, ``lincs-order = 6``, ``lincs-iter = 2``
  - ``constraint-algorithm = none``
  - SETTLE


Ensemble tests
^^^^^^^^^^^^^^

The generated ensembles are tested with Argon (1000 atoms) and water
(900 molecules, with SETTLE and PME) systems, in the following
combinations:

* ``integrator = md``, ``tcoupl = v-rescale``, ``tau-t = 0.1``,
  ``ref-t = 87.0`` (Argon) or ``ref-t = 298.15`` (Water)
* ``integrator = md``, ``tcoupl = v-rescale``, ``tau-t = 0.1``,
  ``ref-t = 87.0`` (Argon) or ``ref-t = 298.15`` (Water), ``pcoupl =
  parrinello-rahman``, ``ref-p = 1.0``, ``compressibility = 4.5e-5``
* ``integrator = md-vv``, ``tcoupl = v-rescale``, ``tau-t = 0.1``,
  ``ref-t = 87.0`` (Argon) or ``ref-t = 298.15`` (Water)
* ``integrator = md-vv``, ``tcoupl = nose-hoover``, ``tau-t = 1.0``,
  ``ref-t = 87.0`` (Argon) or ``ref-t = 298.15`` (Water), ``pcoupl =
  mttk``, ``ref-p = 1.0``, ``compressibility = 4.5e-5``

All thermostats are applied to the entire system (``tc-grps =
system``). The simulations run for 1ns at 2fs time step with Verlet
cut-off. All other settings left to default values.


Building and testing using the build system
-------------------------------------------

Since these tests can not be ran at the same frequency as the current
tests, they are kept strictly opt-in via
``-DGMX_PHYSICAL_VALIDATION=ON``, with
``-DGMX_PHYSICAL_VALIDATION=OFF`` being the default. Independently of
that, all previously existing build targets are unchanged, including
``make check``.

If physical validation is turned on, a number of additional make
targets can be used:

* ``make check`` is unchanged, it builds the main binaries and the unit
  tests, then runs the unit tests and, if available, the regression
  tests.
* ``make check-phys`` builds the main binaries, then runs the physical
  validation tests. **Warning**: This requires to simulate all systems
  and might take several hours on a average machine!
* ``make check-all`` combines ``make check`` and ``make check-phys``.

As the simulations needed to perform the physical validation tests may
take long, it might be advantageous to run them on an external
resource. To enable this, two additional make targets are present:

* ``make check-phys-prepare`` prepares all simulation files under
  ``tests/physicalvalidation`` of the build directory, as well as a
  rudimentary run script in the same directory.
* ``make check-phys-analyze`` runs the same tests as ``make
  check-phys``, but does not simulate the systems. Instead, this
  target assumes that the results can be found under
  ``tests/physicalvalidation`` of the build directory.

The intended usage of these additional targets is to prepare the
simulation files, then run them on a different resource or at a
different time, and later analyze them. If you want to use this, be
aware *(i)* that the run script generated is very simple and might
need (considerable) tuning to work with your setup, and *(ii)* that
the analysis script is sensitive to the folder structure, so make sure
to preserve it when copying the results to / from another resource.

Additionally to the mentioned make targets, a number of internal make
targets are defined. These are not intended to be used directly, but
are necessary to support the functionality described above, especially
the complex dependencies. These internal targets include
``run-ctest``, ``run-ctest-nophys``, ``run-ctest-phys`` and
``run-ctest-phys-analyze`` running the different tests,
``run-physval-sims`` running the simulations for physical validation,
and ``missing-tests-notice``, ``missing-tests-notice-all``,
``missing-phys-val-phys``, ``missing-phys-val-phys-analyze`` and
``missing-phys-val-all`` notifying users about missing tests.


Direct usage of the python script
---------------------------------

The ``make`` commands mentioned above are calling the python script
``tests/physicalvalidation/gmx_physicalvalidation.py``, which can be
used independently of the make system. Use the ``-h`` flag for the
general usage information, and the ``--tests`` for more details on the
available physical validations.

The script requires a ``json`` file defining the tests as an input.
Among other options, it allows to define the |Gromacs| binary and the
working directory to be used, and to decide whether to only prepare
the simulations, prepare and run the simulations, only analyze the
simulations, or do all three steps at once.


Adding new tests
----------------

The available tests are listed in the ``systems.json`` (tests
standardly used for single precision builds) and ``systems_d.json``
(tests standardly used for double precision builds) files in the same
directory, the |Gromacs| files are in the folder ``systems/``.

The ``json`` files lists the different test. Each test has a
``"name"`` attribute, which needs to be unique, a ``"dir"`` attribute,
which denotes the directory of the system (inside the ``systems/``
directory) to be tested, and a ``"test"`` attribute which lists the
validations to be performed on the system. Additionally, the optional
``"grompp_args"`` and ``"mdrun_args"`` attributes allow to pass
specific arguments to ``gmx grompp`` or ``gmx mdrun``, respectively. A
single test can contain several validations, and several independent
tests can be performed on the same input files.

To add a new test to a present system, add the test name and the
arguments to the ``json`` file(s). To use a new system, add a
subfolder in the ``systems/`` directory containing
``input/system.{gro,mdp,top}`` files defining your system.
