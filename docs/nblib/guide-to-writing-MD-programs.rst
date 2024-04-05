Guide to Writing MD Programs
============================

The goal of NB-LIB’s is to enable researchers to programmatically define molecular simulations.
Traditionally these have been performed using a collection of executables and a manual workflow followed by a “black-box” simulation engine.
NB-LIB allows users to script a variety of novel simulation and analysis workflows at a more granular level.

Many possible use cases are facilitated by the flexibility that NB-LIB allows.
These include customized update rules, defining custom forces, or orchestrating swarms of simulations.
NB-LIB also allows for writing conventional MD simulations and analysis.

This document goes over the steps to write MD programs using the API in NB-LIB that exposes features that are a part of the |Gromacs| package.

Global Definitions
------------------

NB-LIB programs are written in C++ so its headers for I/O or advanced tasks must be included.
In addition, one must include the headers for various capabilities and abstractions NB-LIB exposes as well.
This can be directly copied from here.
Finally, we use the namespace ``nblib`` for the data structures defined in the library.
The last line in the block allows one to skip this specifier each time a function or a data structure is used.

.. code:: cpp

   #include <cstdio>

   #include "nblib/box.h"
   #include "nblib/forcecalculator.h"
   #include "nblib/integrator.h"
   #include "nblib/molecules.h"
   #include "nblib/nbkerneloptions.h"
   #include "nblib/particletype.h"
   #include "nblib/simulationstate.h"
   #include "nblib/topology.h"

   using namespace nblib;

Define Particle Data
--------------------

.. code:: cpp

   // Parameters from a GROMOS compatible force-field 2016H66

   struct OWaterAtom
   {
       ParticleName         name = "Ow";
       Mass                 mass = 15.999;
       C6                   c6   = 0.0026173456;
       C12                  c12  = 2.634129e-06;
   };

   struct HwAtom
   {
       ParticleName         name = "Hw";
       Mass                 mass = 1.00784;
       C6                   c6   = 0.0;
       C12                  c12  = 0.0;  
   };

   struct CMethAtom
   {
       ParticleName         name = "Cm";
       Mass                 mass = 12.0107;
       C6                   c6   = 0.01317904;
       C12                  c12  = 34.363044e-06;
   };

   struct HcAtom
   {
       ParticleName         name = "Hc";
       Mass                 mass = 1.00784;
       C6                   c6   = 8.464e-05;
       C12                  c12  = 15.129e-09;  
   };

There can be as many structs of this kind as there are particle types in the system.
Organizing the data like this is not strictly necessary, but is shown for the purpose of clarity.
As shown here, there can be multiple particles that correspond to a single element as atomic mass can vary by molecular context.
For example, the carbon atom in a carboxyl group would have different parameters from one in the methyl group.
We can obtain the parameter set from any standard force-field, or generate new parameters to study new compounds or force fields.
This example comes from the `2016H66 Parameter Set <https://pubs.acs.org/doi/10.1021/acs.jctc.6b00187>`__.

Defining Coordinates, Velocities and Force Buffers
--------------------------------------------------

.. code:: cpp

   std::vector<gmx::RVec> coordinates = {
       { 0.794, 1.439, 0.610 }, { 1.397, 0.673, 1.916 }, { 0.659, 1.080, 0.573 },
       { 1.105, 0.090, 3.431 }, { 1.741, 1.291, 3.432 }, { 1.936, 1.441, 5.873 },
       { 0.960, 2.246, 1.659 }, { 0.382, 3.023, 2.793 }, { 0.053, 4.857, 4.242 },
       { 2.655, 5.057, 2.211 }, { 4.114, 0.737, 0.614 }, { 5.977, 5.104, 5.217 },
   };

   std::vector<gmx::RVec> velocities = {
       { 0.0055, -0.1400, 0.2127 }, { 0.0930, -0.0160, -0.0086 }, { 0.1678, 0.2476, -0.0660 },
       { 0.1591, -0.0934, -0.0835 }, { -0.0317, 0.0573, 0.1453 }, { 0.0597, 0.0013, -0.0462 },
       { 0.0484, -0.0357, 0.0168 }, { 0.0530, 0.0295, -0.2694 }, { -0.0550, -0.0896, 0.0494 },
       { -0.0799, -0.2534, -0.0079 }, { 0.0436, -0.1557, 0.1849 }, { -0.0214, 0.0446, 0.0758},
   };

   std::vector<gmx::RVec> forces = {
       { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
       { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
       { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
       { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
   };

We can initialize coordinates for our particles using ``std::vector`` of ``gmx::RVec`` which is a specific data type for holding 3D vector quantities. `Doxygen page on RVec here`__.

  __ doxygen-ref-rvec_

Writing the MD Program
----------------------

As with as any basic C++ program, there needs to be a ``main()`` function.


Define ParticleTypes
~~~~~~~~~~~~~~~~~~~~

.. code:: cpp

   int main()
   {
       // Bring the parameter structs to scope
       OwAtom      owAtom;
       HwAtom      hwAtom;
       CMethAtom   cmethAtom;
       HcAtom      hcAtom;
     
       // Create the particles
       ParticleType Ow(owAtom.name, owAtom.mass);
       ParticleType Hw(hwAtom.name, hwAtom.mass);
       ParticleType Cm(cmethAtom.name, cmethAtom.mass);
       ParticleType Hc(hcAtom.name, hcAtom.mass);

As before, the helper struct to define ``ParticleType`` data is not strictly needed, but is shown for clarity.
The line ``ParticleType CMethAtom(ParticleName("Cm"), Mass(12.0107));`` would be sufficient.

Define Non-Bonded Interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: cpp

   ParticleTypeInteractions interactions(CombinationRule::Geometric);

   // add non-bonded interactions for the particle types
   interactions.add(owAtom.name, owAtom.c6, owAtom.c12);
   interactions.add(hwAtom.name, hwAtom.c6, hwAtom.c12);
   interactions.add(cmethAtom.name, cmethAtom.c6, cmethAtom.c12);
   interactions.add(hcAtom.name, hcAtom.c6, hcAtom.c12);

For the Lennard-Jones interactions, we define a ``ParticleTypeInteractions`` object.
Each particle of the ``ParticleType`` interacts with each other based on the ``C6`` and ``C12`` parameters.
These parameters of the two different particles are averaged using ``Geometric`` or ``LorentzBerthelot`` ``CombinationRule``.
More details `here <http://manual.gromacs.org/documentation/2019/reference-manual/functions/nonbonded-interactions.html#the-lennard-jones-interaction>`__.
By default ``CombinationRule::Geometric`` is selected.

We add the interaction parameters of each of the particle types into the ``ParticleTypeInteractions`` object.
The result is a table that has interactions specified for all ``ParticleType`` pairs.
The following matrix describes the pair-wise C6 parameter created using ``CombinationRule::Geometric``.

== ====== === ======= =======
#  Ow     Hw  Cm      Hc
== ====== === ======= =======
Ow 0.0026 0.0 0.42    4.7e-4
Hw 0.0    0.0 0.0     0.0
Cm 0.42   0.0 0.013   1.05e-3
Hc 4.7e-4 0.0 1.05e-3 8.5e-5
== ====== === ======= =======

For a particular interaction pair, the user can also override the specified ``CombinationRule`` with custom parameters.
The following overload would replace the parameters computed from a ``CombinationRule``  between ``Ow`` and ``Cm`` particle types.

.. code:: cpp

   interactions.add("Ow", "Cm", 0.42, 42e-6);

To facilitate modular, reusable code, it is possible to combine multiple ``ParticleTypeInteractions`` objects.
Assuming ``otherInteractions`` is defined, this can be done with ``interactions.merge(otherInteractions)``

Define Molecules
~~~~~~~~~~~~~~~~

.. code:: cpp

   Molecule water("Water");
   Molecule methane("Methane");

   water.addParticle(ParticleName("O"), Ow);
   water.addParticle(ParticleName("H1"), Hw);
   water.addParticle(ParticleName("H2"), Hw);

   water.addExclusion("H1", "O");
   water.addExclusion("H2", "O");

   methane.addParticle(ParticleName("C"), Cm);
   methane.addParticle(ParticleName("H1"), Hc);
   methane.addParticle(ParticleName("H2"), Hc);
   methane.addParticle(ParticleName("H3"), Hc);
   methane.addParticle(ParticleName("H4"), Hc);

   methane.addExclusion("H1", "C");
   methane.addExclusion("H2", "C");
   methane.addExclusion("H3", "C");
   methane.addExclusion("H4", "C");

We begin declaring molecules with their constituent particles.
A string identifier must uniquely identify a specific particle within the molecule.
It is also possible to define partial charges on each particle for the computation of Coulomb interactions.
``water.addParticle(ParticleName("O"), Charge(-0.04), Ow);``

Adding exclusions ensures that non-bonded interactions are only computed when necessary.
For example, if two  particles share a bond, the potential energy of the bond makes the non-bonded term negligible.
Particle self-exclusions are enabled by default.
We use the unique identifiers specified during ``addParticle()`` for this and the listed interactions later.

Define Listed Interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~

Within a molecule, one can define interactions such as bonds, angles and dihedrals between the constituent particles.
NB-LIB provides concrete implementations of several commonly used 2, 3 and 4 center interactions.

.. code:: cpp

   HarmonicBondType ohHarmonicBond(1, 1);
   HarmonicBondType hcHarmonicBond(2, 1);

   DefaultAngle hohAngle(Degrees(120), 1);
   DefaultAngle hchAngle(Degrees(109.5), 1);

   //add harmonic bonds for water
   water.addInteraction("O", "H1", ohHarmonicBond);
   water.addInteraction("O", "H2", ohHarmonicBond);

   // add the angle for water
   water.addInteraction("H1", "O", "H2", hohAngle);

   // add harmonic bonds for methane
   methane.addInteraction("H1", "C", hcHarmonicBond);
   methane.addInteraction("H2", "C", hcHarmonicBond);
   methane.addInteraction("H3", "C", hcHarmonicBond);
   methane.addInteraction("H4", "C", hhcHarmonicBondc);

   // add the angles for methane
   methane.addInteraction("H1", "C", "H2", hchAngle);
   methane.addInteraction("H1", "C", "H3", hchAngle);
   methane.addInteraction("H1", "C", "H4", hchAngle);
   methane.addInteraction("H2", "C", "H3", hchAngle);
   methane.addInteraction("H2", "C", "H4", hchAngle);
   methane.addInteraction("H3", "C", "H4", hchAngle);

Define Options for the Simulation and Non-Bonded Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: cpp

   // Define a box for the simulation
   Box box(6.05449);

   // Define options for the non-bonded kernels
   NBKernelOptions options;

One can define the bounding box either with a single argument for a cube and 3 arguments to specify length, breadth and height separately.

``NBKernelOptions`` contains a set of flags and configuration options for both hardware context and the relevant calculations for the simulation.
The following table describes the possible options that can be set.

+----------------------+------+---------------------------------------+
| Flag or Config       | Type | Implications                          |
| Option               |      |                                       |
+======================+======+=======================================+
| ``useGpu``           | Bool | Use GPU for non-bonded computations   |
|                      | ean  |                                       |
+----------------------+------+---------------------------------------+
| ``numThreads``       | Inte | Number of CPU threads to use          |
|                      | ger  |                                       |
+----------------------+------+---------------------------------------+
| ``nbnxmSimd``        | Enum | Kernel SIMD type                      |
|                      |      | (``SimdAuto``/``SimdNo``/``Simd4XM``/ |
|                      |      | ``Simd2XMM``)                         |
+----------------------+------+---------------------------------------+
| ``useHalfLJOptimizat | Bool | Enable i-cluster half-LJ optimization |
| ion``                | ean  |                                       |
+----------------------+------+---------------------------------------+
| ``pairlistCutoff``   | Real | Specify pairlist and interaction      |
|                      |      | cut-off                               |
+----------------------+------+---------------------------------------+
| ``computeVirialAndEn | Bool | Enable energy computations            |
| ergy``               | ean  |                                       |
+----------------------+------+---------------------------------------+
| ``coulombType``      | Enum | Coulomb interaction function          |
|                      |      | (``Pme``/``Cutoff``/``ReactionField`` |
|                      |      | )                                     |
+----------------------+------+---------------------------------------+
| ``useTabulatedEwaldC | Bool | Use tabulated PME grid correction     |
| orr``                | ean  | instead of analytical                 |
+----------------------+------+---------------------------------------+
| ``numIterations``    | Inte | Specify number of iterations for each |
|                      | ger  | kernel                                |
+----------------------+------+---------------------------------------+
| ``cyclesPerPair``    | Bool | Enable printing cycles/pair instead   |
|                      | ean  | of pairs/cycle                        |
+----------------------+------+---------------------------------------+
| ``timestep``         | Real | Specify the time step                 |
+----------------------+------+---------------------------------------+

Define Topology and Simulation State
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We build the system topology using the ``TopologyBuilder`` class.
We add the ``Molecule`` objects that we defined previously along with the ``ParticleTypesInteractions`` using its public functions.
We get the actual ``Topology`` object complete with all exclusions, interaction maps and listed interaction data constructed based on the defined entities using the ``buildTopology()``\ function.

.. code:: cpp

   TopologyBuilder topologyBuilder;

   // add molecules
   topologyBuilder.addMolecule(water, 10);
   topologyBuilder.addMolecule(methane, 10);

   // add non-bonded interaction map
   topologyBuilder.addParticleTypesInteractions(interactions);

   Topology topology = topologyBuilder.buildTopology();

We now have all we need to fully describe our system using the ``SimulationState`` object.
This is built using the topology, the box, and the particle coordinates and velocities.
This object serves as a snapshot of the system that can be used for analysis or to start simulations from known states.

.. code:: cpp

   SimulationState simulationState(coordinates, velocities, forces, box, topology);



Writing the MD Loop
~~~~~~~~~~~~~~~~~~~

Now that we have fully described our system and the problem, we need two entities to write an MD loop.
The first is the ``ForceCalculator`` and the second is an Integrator.
NB-LIB comes with a ``LeapFrog`` integrator but it is also possible for users to write custom integrators.

.. code:: cpp

   // The force calculator contains all the data needed to compute forces
   ForceCalculator forceCalculator(simulationState, options);

   // Integration requires masses, positions, and forces
   LeapFrog integrator(simulationState);

   // Allocate a force buffer
   gmx::ArrayRef<gmx::RVec> userForces(topology.numParticles());

   // MD Loop
   int numSteps = 100;

   for (i = 0; i < numSteps; i++)
   {
     userForces = forceCalculator.compute();

     // The forces are not automatically updated in case the user wants to add their own
     std::copy(userForces.begin(), userForces.end(), begin(simulationState.forces()));

     // Integrate with a time step of 1 fs
     integrator.integrate(1.0);
   }

   return 0;
   } // main

.. _doxygen-ref-rvec: ../doxygen/html-lib/namespacegmx.xhtml#a139c1919a9680de4ad1450f42e37d33b