This folder contains forcefield parameters for general (G)MINFF, as well as tailored (T)MINFF parameters for specific mineral (see list below).

**Lennard-Jones parameters for the general MINFF**

Lennard-Jones parameters for all metal sites optimized for over 30 mineral synchronously, 
using the Lennard-Jones parameters for the oxygen atomtypes taken from OPC3 water model fixed.  

* Filename ffnonbonded_gminff.itp:
	- General MINFF for k=0, 250, 500, 1500 kJ/mol/rad2
	- CLAYFF
	- Several ion-pair potentials and water models

	The mineral set and the angle force constant are selected with two separate
	-D flags in the .mdp file:

	    define = -DMINFF_k500 -DOPC3 -DOPC3_HFE_LM

	where MINFF_k0, MINFF_k250, MINFF_k500 or MINFF_k1500 picks the O-M-O angle
	force constant in kJ/mol/rad2. No separate keyword for the parameter set is
	needed here; including ffnonbonded_gminff.itp (or ffnonbonded.itp) is what
	selects the general parameters. The older form -DGMINFF_k500 no longer works
	and must be replaced by -DMINFF_k500.

* Filename ffbonded_gminff.itp
	- bonded parameters for the general MINFF. Identical in content to
	  ffbonded.itp; kept under its own name so that each parameter set has a
	  matching forcefield/ffnonbonded/ffbonded trio.

* Filename forcefield_gminff.itp
	- include this to use the general parameters; it pulls in
	  ffnonbonded_gminff.itp and ffbonded_gminff.itp. forcefield.itp is
	  equivalent and includes ffnonbonded.itp and ffbonded.itp instead.


**Lennard-Jones parameters for the tailored MINFF**

Note that the Lennard-Jones parameters inn step 1 for each mineral were first optimized for only the metal sites only.
In a second step the oxygen Lennard-Jones parameters were optimized, keeping the metal parameters from step 1 fixed. 
These new oxygen parameters optimized in step 2 are commented by a ; and not necessarily better than the parameters optmized in step 1.

* Filename ffnonbonded_tminff.itp
	- tailored MINFF for all four angle force constants, sorted by mineral
	- CLAYFF
	- Several ion-pair potentials and water models

	This one file replaces the former ffnonbonded_tminff_k0|k250|k500|k1500.itp
	and ffnonbonded_tminff_all_k_sorted_by_mineral.itp, which have been removed.

	The mineral and the angle force constant are now selected independently,
	with two separate -D flags in the .mdp file. Both are required:

	    define = -DMontmorillonite -DMINFF_k500 -DOPC3 -DOPC3_HFE_LM

	where MINFF_k0, MINFF_k250, MINFF_k500 or MINFF_k1500 picks the O-M-O angle
	force constant in kJ/mol/rad2. The older combined form (-DMontmorillonite_k500)
	is no longer used. For the general sets, see above: -DGMINFF_k500 has likewise
	been replaced by -DMINFF_k500.

* Filename ffbonded_tminff.itp
	- bonded parameters for the tailored MINFF, guarded by mineral name alone
	  since they do not depend on the angle force constant: the O-M-O stiffness
	  is applied through KANGLE in each mineral .itp file, not here.
	- Only the terms whose atomtypes the chosen mineral actually declares are
	  emitted, so grompp never sees a bondtype referring to an undeclared type.
	  Minerals without hydroxyl groups contribute no bonded terms.

* Filename forcefield_tminff.itp
	- include this instead of forcefield.itp to use the tailored parameters;
	  it pulls in ffnonbonded_tminff.itp and ffbonded_tminff.itp.
	
Mineral list (# 13 and 22 missing, and that's ok)
1	Kaolinite
2	Pyrophyllite
3	Talc
4	Forsterite
5	Brucite
6	Corundum
7	Quartz
8	Gibbsite
9	Li2O
10	Coesite
11	Cristobalite
12	Maghemite
14	Akdalaite
15	Boehmite
16	Diaspore
17	Periclase
18	Goethite
19	Hematite
20	Lepidocrocite
21	Wustite
23	CaO
24	Portlandite
23	CaF2
26	Nontronite
27	Montmorillonite
28	Dickite
29	Hectorite-F
30	Hectorite-H
31	Nacrite
32	Imogolite
33	Anatase
34	Rutile
35	cis_Oct_Fe2_cis
36	cis_Oct_Fe2_trans
37	cis_Oct_Mg2cis_Fe3cis
38	cis_Oct_Mg2cis_Fe3trans
39	cis_Oct_Mg2trans_Fe3cis
40	cis_Oct_Mg2trans_Fe3trans
41	cis_Tet_Fe3
42	trans_Oct_Fe2_cis
43	trans_Oct_Mg2cis_Fe3cis
44	trans_Tet_Fe3
45	Muscovite
