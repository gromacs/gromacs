General information:

Defaults for tests unless otherwise stated:
integrator               = md
dt                       = 0.001
sc-alpha                 = 0.5
sc-power                 = 1
sc-r-power               = 6
sc-sigma                 = 0.3
sc-coul                  = no
separate-dhdl-file       = no
dhdl-print-energy        = no

Test topology:

  The molecule is acetamide in a water box, but far less than liquid
  density (56 for a box of ~25 nm^3) to reduce expense. It would
  likely collapse to a droplet eventually, but we run for far less
  time.  There are 9 heavy molecules, 9 bonds, 10 1-4 pairs, 12
  angles, 10 proper dihedrals, and 2 improper dihedrals defined.
  Defaults use function 1 for bonds, pairs, angles, function 3 for
  dihedrals, and function 1 for impropers.
  
  
We test the separate-dhdl-file and dhdl-print-energy separately; we
check with these test that the data is put in the correct place in the
data structures, and then test with separate tests that it is exported
correctly into the .dhdl.xvg files.

Tests in this folder
  - coulandvdwintramol
    Tests turning off both vdw and coulomb while coupling intramolecular unteractions,
    using only the mdp options; no explicitly defined A and B states in the topology. 
    only a single init-lambda float is used, there is no array of lambda values, and therefore no init-lambda-state.
    Other relevant settings that could affect free energy:
    	  Constraints: h-bond with Lincs,
	  Coul: PME (Potential-shift-Verlet)
	  Vdw: Cutoff (Force-switch)
	  Dispersion: energy and pressure
	  Ensemble: constant T, berendsen
	  Frequencies: nstcalcenergy=1,nstenergy=1,nstdhdl=5,nsttcouple=1
  - coulandvdwsequential_coul
    Same as coulandvdwintramol, except it tests turning off both vdw and q, without
    coupling intramolecular interactions.  A lambda
    array is used, and the lambda value is in the region where vdw is turned on, and only coul
    is changing.
  - coulandvdwsequential_vdw
    Same as coulandvdwintramol, except it tests turning off both vdw and q, without
    coupling intramolecular interactions.  A lambda array is used, and the lambda value is
    in the region where coul is turned off, and only vdw is changing.
  - coulandvdwtogether
    Same as coulandvdwintermol EXCEPT without coupling intramolecular interactions.
  - expanded
    Tests expanded ensemble, with nstexpanded = 5.  Tests turning on vdw and coulomb, using a
    lambda array.  Only coul-lambda and vdw-lambda are specified, with the couple-intramol=no.
    coul and vdew are turned off sequentially.  Energies at all lambda states are calculated.
    Topology does not explicitly state A and B states.  Integrator is md-vv.
    Other relevant settings that could affect free energy:
    	  Constraints: h-bond with Lincs,
	  Coul: PME (Potential-shift-Verlet)
	  Vdw: shift (Potential-shift-Verlet)
	  Dispersion: energy and pressure
	  Ensemble: constant T, berendsen
	  Frequencies: nstcalcenergy=1,nstenergy=1,nstdhdl=10,nsttcouple=1  
  - relative
    Free energy changesare controlled by the .top file, not the .mdp file, as there is no couple_moltype,
    and A and B states are defined explicitly in the topol.top file.
    coul-lambdas, vdw-lambdas, and bonded-lambdas arrays vary from 0 to 1, with coulomb first being turned off.
    - Topolgy changes:
      - Atomtype opls_240 has a vdw parameter defined to accentuate the changes in lambda.
      - improper dihedrals have A and B states.
      - charges change on all atoms from state A to B.
      - masses stay the same
      - atom types from state A to B only change on atom 240.
    TODO: if couple-moltype is zero, why are the couple-lambda0 defined?
    Other relevant settings that could affect free energy:
    	  Same as coulandvdwintramol

  - relative-position-restraints
      Like relative, but includes posres.itp, which have A and B states on the position restraints, so
      should have lambda dependence. (also, does not have COM removal set, since it's unneeded,
      and the restraint-lambdas array is used.
  - restraints
      Very imilar to relative, except for the free energy parameters.
      Restraints are turned on. A and B states are the same except for
      the addition of dihedral_restraints, angle_restraints, and distance_restraints,
      each of which have A and B states.  The free energies are controlled in the topologies, but
      the _restraint terms are the only ones with A and B states.
  - simtemp
      Similar to relative, except only simulated tempering and expanded ensemble options are set.
      No topology B states are set.
  - testcoul
      uses md-vv, and uses the mdp to turn off both coul and vdw.
      Uses lambda array, turns off vdw and coul sequentially, and tests the coul region.
      otherwise fairly similar to coulandvdwtogether, and could possibly be merged later.
  - transformAtoB
      uses md-vv. Topology changes all atom types and charges with explicitly stated A and B states.
      Bonded remain the same, lambda array used, set using fep-lambdas. 
      Other relevant settings that could affect free energy: same as coulandvdwtogether.
  - vdwalone
      uses md-vv.  Turns off the interactions by defining B states with zero charge and vdw
      epsilon in the topology.  ALSO uses couple-lambda0 and couple-lambda1=0, which is weird, but
      there is no couple-moltype. Does not use a lambda array, but uses a single init-lambda.
      TODO: Results seem to be the same with and without couple-lambda0 and couple-lamda1, but there is an
      extra warning when they are not specified.

TODOS:
    - put in tests for individual bonded terms alone, of different
functional types.  These tests can be done without a water box.
    - put in tests for different types of long range cutoffs
    - put in tests that check for free energy changes for different combinations
      of nstenergy/nstcalcenergy/nstdhdl/nstexpanded/nsttemp/nstpres
    - tests with different thermostats/barostats
    - a bunch more ways to change from A to B
    - tests without constraints
    - are the energies printed out correctly at different dhdls.
    - what about different sc-powers/sc-rpower
    - add a separate test that makes sure the contents of nstdhdl get copied into the files correctly
    - tests of expanded ensemble options (mostly movement options, some testing for weight accumulation options)
    - tests of test of simulated tempering options
    - tests with different settings for lambdas (including changing mass)
