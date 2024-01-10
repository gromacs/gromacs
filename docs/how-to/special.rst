.. _reference manual: gmx-manual-parent-dir_

.. _gmx-membrane:

Running membrane simulations in |Gromacs|
-----------------------------------------

Running Membrane Simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users frequently encounter problems when running simulations of lipid bilayers, especially
when a protein is involved. Users seeking to simulate membrane proteins may find this
`tutorial <https://tutorials.gromacs.org/membrane-protein.html>`__ useful.

One protocol for the simulation of membrane proteins consists of the following steps:

#. Choose a force field for which you have parameters for the protein and lipids.
#. Insert the protein into the membrane. (For instance, use g_membed on a pre-formed bilayer or do a
   coarse-grained self-assembly simulation and then convert back to the atomistic representation.)
#. Solvate the system and add ions to neutralize excess charges and adjust the final ion concentration.
#. Energy minimize.
#. Let the membrane adjust to the protein. Typically run MD for ~5-10ns with restraints (1000 kJ/(mol nm2) on all protein heavy atoms.
#. Equilibrate without restraints.
#. Run production MD.

Adding waters with genbox
^^^^^^^^^^^^^^^^^^^^^^^^^

When generating waters around a pre-formed lipid membrane with :ref:`solvate <gmx solvate>` you may find that
water molecules get introduced into interstices in the membrane. There are several approaches to removing these, including

* a short MD run to get the hydrophobic effect to exclude these waters. In general this
  is sufficient to reach a water-free hydrophobic phase, as the molecules are usually
  expelled quickly and without disrupting the general structure. If your setup relies
  on a completely water-free hydrophobic phase at the start, you can try to follow
  the advice below:
* Set the ``-radius`` option in :ref:`gmx solvate` to change the water exclusion radius,
* copy ``vdwradii.dat`` from your ``$GMXLIB`` location to the working directory, and edit it to
  increase the radii of your lipid atoms (between 0.35 and 0.5nm is suggested for carbon) to
  prevent :ref:`solvate <gmx solvate>` from seeing interstices large enough for water insertion,
* editing your structure by hand to delete them (remembering to adjust your atom count for :ref:`gro` files
  and to account for any changes in the :ref:`topology <top>`), or
* use a script someone wrote to remove them.

External material
^^^^^^^^^^^^^^^^^

* `Membrane simulations slides <https://extras.csc.fi/chem/courses/gmx2007/Erik_Talks/membrane_simulations.pdf>`_ ,
  `membrane simulations video <https://video.csc.fi/playlist/dedicated/0_7z3nas0q/0_0tr9yd2p>`_ - (Erik Lindahl).
* `tutorial for membrane protein simulations
  <http://www.mdtutorials.com/gmx/membrane_protein/index.html>`__ - designed to demonstrate what sorts of
  questions and problems occur when simulating proteins that are embedded within a lipid bilayer.
* `Combining the OPLS-AA forcefield with the Berger lipids <http://pomes.biochemistry.utoronto.ca/files/lipidCombinationRules.pdf>`_
  A detailed description of the motivation, method, and testing.

* Several Topologies for membrane proteins with different force fields gaff, charmm berger
  Shirley W. I. Siu, Robert Vacha, Pavel Jungwirth, Rainer A. Böckmann: Biomolecular simulations of membranes:
  `Physical properties from different force fields <https://doi.org/10.1063/1.2897760>`_.
* `Lipidbook <https://www.lipidbook.org/>`_ is a public repository for force-field parameters of lipids,
  detergents and other molecules that are used in
  the simulation of membranes and membrane proteins. It is described in: J. Domański, P. Stansfeld, M.S.P. Sansom,
  and O. Beckstein. J. Membrane Biol. 236 (2010), 255—258. `doi:10.1007/s00232-010-9296-8 <http://dx.doi.org/10.1007/s00232-010-9296-8>`_.


Parameterization of novel molecules
-----------------------------------

Most of your parametrization questions/problems can be resolved very simply, by remembering the following two rules:

* **You should not mix and match force fields**. :ref:`Force fields <gmx-force-field>` are (at best) designed to be self-consistent,
  and will not typically work well with other force fields. If you simulate part of your system with one
  force field and another part with a different force field which is not parametrized with the first force
  field in mind, your results will probably be questionable, and hopefully reviewers will be concerned.
  Pick a force field. Use that force field.
* If you need to develop new parameters, derive them in a manner consistent with how the rest of the force field
  was originally derived, which means that you will need to review the original literature. There isn't a single
  right way to derive force field parameters; what you need is to derive parameters that are consistent with the rest
  of the force field. How you go about doing this depends on which force field you want to use. For example, with
  AMBER force fields, deriving parameters for a non-standard amino acid would probably involve doing a number of
  different quantum calculations, while deriving GROMOS or OPLS parameters might involve more (a) fitting various fluid
  and liquid-state properties, and (b) adjusting parameters based on experience/chemical intuition/analogy. Some
  suggestions for automated approaches can be found :doc:`here <../user-guide/system-preparation>`.

It would be wise to have a reasonable amount of simulation experience with |Gromacs| before
attempting to parametrize new force fields, or new molecules for existing force fields.
These are expert topics, and not suitable for giving to (say) undergraduate students for
a research project, unless you like expensive quasi-random number generators. A very thorough knowledge
of :ref:`Chapter 5: Interaction function and force fields <ff>` of the |Gromacs| Reference Manual will be required. If you haven't been warned
strongly enough, please read below about parametrization for exotic species.

Another bit of advice: Don't be more haphazard in obtaining parameters than you would be buying
fine jewellery. Just because the guy on the street offers to sell you a *diamond* necklace for $10
doesn't mean that's where you should buy one. Similarly, it isn't necessarily the best strategy
to just download parameters for your molecule of interest from the website of someone you've
never heard of, especially if they don't explain how they got the parameters.

Be forewarned about using `PRODRG <http://davapc1.bioch.dundee.ac.uk/cgi-bin/prodrg>`_ topologies
without verifying their contents: the artifacts of doing so are now `published <http://pubs.acs.org/doi/abs/10.1021/ci100335w>`_,
along with some tips for properly deriving parameters for the GROMOS family of force fields.

Exotic Species
^^^^^^^^^^^^^^

So, you want to simulate a protein/nucleic acid system, but it binds various exotic metal
ions (ruthenium?), or there is an iron-sulfur cluster essential for its functionality, or similar.
But, (unfortunately?) there aren't parameters available for these in the force field you want
to use. What should you do? You shoot an e-mail to the |Gromacs| `user discussion forum`_, and get referred to the FAQs.

If you really insist on simulating these in molecular dynamics, you'll need to obtain parameters
for them, either from the literature, or by doing your own parametrization. But before doing so,
it's probably important to stop and think, as sometimes there is a reason there may not already
be parameters for such atoms/clusters. In particular, here are a couple of basic questions you
can ask yourself to see whether it's reasonable to develop/obtain standard parameters for these and use them in molecular dynamics:

* Are quantum effects (i.e. charge transfer) likely to be important? (i.e., if you have a
  divalent metal ion in an enzyme active site and are interested in studying enzyme
  functionality, this is probably a huge issue).
* Are standard force field parametrization techniques used for my force field of choice
  likely to fail for an atom/cluster of this type? (i.e. because Hartree-Fock 6-31G* can't
  adequately describe transition metals, for example)

If the answer to either of these questions is "Yes", you may want to consider doing your
simulations with something other than classical molecular dynamics.

Even if the answer to both of these is "No", you probably want to consult with someone who
is an expert on the compounds you're interested in, before attempting your own parametrization.
Further, you probably want to try parametrizing something more straightforward before you embark on one of these.


Potential of Mean Force
-----------------------

The potential of mean force (PMF) is defined as the potential that gives an average force over all the
configurations of a given system.  There are several ways to calculate the PMF in |Gromacs|, probably
the most common of which is to make use of the pull code. The steps for obtaining a PMF using umbrella
sampling, which allows for sampling of statistically-improbable states, are:

* Generate a series of configurations along a reaction coordinate (from a steered MD simulation,
  a normal MD simulation, or from some arbitrarily-created configurations)
* Use umbrella sampling to restrain these configurations within sampling windows.
* Use :ref:`gmx wham` to make use of the WHAM algorithm to reconstruct a PMF curve.

A more detailed tutorial is linked `here for umbrella
sampling <https://tutorials.gromacs.org/umbrella-sampling.html>`__.


Single-Point Energy
-------------------

Computing the energy of a single configuration is an operation that is sometimes useful. The best
way to do this with |Gromacs| is with the :ref:`mdrun <gmx mdrun>` ``-rerun`` mechanism, which
applies the model physics in the :ref:`tpr` to the configuration in the trajectory or coordinate file supplied to mdrun.

::

    mdrun -s input.tpr -rerun configuration.pdb

Note that the configuration supplied must match the topology you used when generating the :ref:`tpr`
file with :ref:`grompp <gmx grompp>`. The configuration you supplied to :ref:`grompp <gmx grompp>`
is irrelevant, except perhaps for atom names. You can also use this feature with energy groups
(see the `Reference manual`_), or with a trajectory of multiple configurations (and in this case,
by default :ref:`mdrun <gmx mdrun>` will do neighbour searching for each configuration, because
it can make no assumptions about the inputs being similar).

A zero-step energy minimization does a step before reporting the energy, and a zero-step MD run
has (avoidable) complications related to catering to possible restarts in the presence of
constraints, so neither of those procedures are recommended.


Carbon Nanotube
---------------

Robert Johnson's Tips
^^^^^^^^^^^^^^^^^^^^^

Taken from Robert Johnson's posts on the `gmx-users mailing list archive`_.

* Be absolutely sure that the "terminal" carbon atoms are sharing a bond in the topology file.
* Use ``periodic_molecules = yes`` in your :ref:`mdp` file for input in :ref:`gmx grompp`.
* Even if the topology is correct, crumpling may occur if you place the nanotube in a box of wrong
  dimension, so use `VMD`_ to visualize the nanotube and its periodic images and make sure that the
  space between images is correct. If the spacing is too small or too big, there will be a large amount
  of stress induced in the tube which will lead to crumpling or stretching.
* Don't apply pressure coupling along the axis of the nanotube. In fact, for debugging purposes,
  it might be better to turn off pressure coupling altogether until you figure out if anything
  is going wrong, and if so, what.
* When using :ref:`x2top <gmx x2top>` with a specific force field, things are assumed about the
  connectivity of the molecule. The terminal carbon atoms of your nanotube will only be bonded to,
  at most, 2 other carbons, if periodic, or one if non-periodic and capped with hydrogens.
* You can generate an "infinite" nanotube with the ``-pbc`` option to :ref:`x2top <gmx x2top>`.
  Here, :ref:`x2top <gmx x2top>` will recognize that the terminal C atoms actually share a
  chemical bond. Thus, when you use :ref:`grompp <gmx grompp>` you won't get an error
  about a single bonded C.

 
Andrea Minoia's tutorial
^^^^^^^^^^^^^^^^^^^^^^^^

Modeling Carbon Nanotubes with |Gromacs| (also archived
as http://chembytes.wikidot.com/grocnt) contains
everything to set up simple simulations of a CNT using OPLS-AA
parameters. Structures of simple CNTs can
be easily generated e.g. by `buildCstruct`_ (Python script that also adds
terminal hydrogens) or `TubeGen Online`_ (just copy and paste the
PDB output into a file and name it cnt.pdb).

To make it work with modern |Gromacs| you'll probably want to do the following:

* make a directory cnt_oplsaa.ff
* In this directory, create the following files, using the data from the tutorial page:

  * forcefield.itp from the file in section :ref:`itp`
  * atomnames2types.n2t from the file in section :ref:`n2t`
  * aminoacids.rtp from the file in section :ref:`rtp`

* generate a topology with the custom forcefield (the cnt_oplsaa.ff directory must be in the same directory as where the :ref:`gmx x2top`
  command is run or it must be found on the GMXLIB path), ``-noparam`` instructs :ref:`gmx x2top` to not use
  bond/angle/dihedral force constants from the command line (-kb, -ka, -kd) but rely on the force field files;
  however, this necessitates the next step (fixing the dihedral functions)

::

    gmx x2top -f cnt.gro -o cnt.top -ff cnt_oplsaa -name CNT -noparam

The function type for the dihedrals is set to '1' by :ref:`gmx x2top` but the force field file specifies type '3'.
Therefore, replace func type  '1' with '3' in the ``[ dihedrals ]`` section of the topology file. A quick way
is to use sed (but you might have to adapt this to your operating system; also manually look at the top file
and check that you only changed the dihedral func types):

::

    sed -i~ '/\[ dihedrals \]/,/\[ system \]/s/1 *$/3/' cnt.top

Once you have the topology you can set up your system. For instance, a simple in-vacuo simulation (using your
favourite parameters in em.\ :ref:`mdp` and md.\ :ref:`mdp`):

Put into a slightly bigger box:

::

    gmx editconf -f cnt.gro -o boxed.gro -bt dodecahedron -d 1

Energy minimise in vacuuo:

::

    gmx grompp -f em.mdp -c boxed.gro -p cnt.top -o em.tpr
    gmx mdrun -v -deffnm em

MD in vacuuo:

::

    gmx grompp -f md.mdp -c em.gro -p cnt.top -o md.tpr
    gmx mdrun -v -deffnm md

Look at trajectory:

::

    gmx trjconv -f md.xtc -s md.tpr -o md_centered.xtc -pbc mol -center
    gmx trjconv -s md.tpr -f md_centered.xtc -o md_fit.xtc -fit rot+trans
    vmd em.gro md_fit.xtc

.. _buildCstruct: http://chembytes.wikidot.com/buildcstruct
.. _TubeGen Online: http://turin.nss.udel.edu/research/tubegenonline.html


