.. _mimic:

MiMiC Hybrid Quantum Mechanical/Molecular Mechanical simulations
----------------------------------------------------------------

This section describes the coupling to a novel QM/MM interface.
The Multiscale Modeling in Computational Chemistry (MiMiC) interface
combines |Gromacs| with the `CPMD QM code <http://cpmd.org/>`__.
To find information about other QM/MM implementations in
|Gromacs| please refer to the section :ref:`qmmm`.
Within a QM/MM approach, typically a small part of the system
(e.g. active site of an enzyme where a chemical reaction can take place)
is treated at the QM level of theory (as we cannot neglect electronic
degrees of freedom while describing some processes e.g.  chemical
reactions), while the rest of the system (remainder of the
protein, solvent, etc.) is described by the classical forcefield (MM).

Overview
^^^^^^^^
MiMiC implements the  QM/MM coupling scheme developed by the group
of Prof. U. Roethlisberger described in
\ :ref:`180 <refRoethlisbergerQMMM>`. This additive
scheme uses electrostatic embedding of the classical system within
the quantum Hamiltonian. The total QM/MM energy is calculated as
a sum of subsystem contributions:

   .. math::

      E_{tot} = E_{QM}+E_{MM}+E_{QM/MM}

The QM contribution is computed by CPMD, while the MM part is
processed by |Gromacs| and the cross terms are treated by the
MiMiC interface. Cross terms, i.e. the terms involving simultaneously
atoms from the QM region and atoms from the MM region consist of
both bonded and non-bonded interactions. 

The bonded interactions are taken from the forcefield used to
describe the MM part. Whenever there is a chemical bond crossing
the QM/MM boundary additional care has to be taken to handle this
situation correctly. Otherwise the QM atom involved in the cut bond
is left with an unsaturated electronic orbital leading to
unphysical system behaviour. Therefore, the dangling bond has to be capped
with another QM atom. There are two different options available
in CPMD for bond capping:

#. Hydrogen capping - the simplest approach is to cap the bond with a
   hydrogen atom, constraining its relative position
   
#. Link atom pseudo-potential - this strategy uses an ad-hoc pseudo-potential
   developed to cap the bond. This pseudo-potential would represent the real
   atom and, thus, will not require the bond constraint.
   
As in standard forcefields, the non-bonded contributions to :math:`E_{QM/MM}`
can be separated into van der Waals and electrostatic contributions.
The first contribution is again taken from the MM forcefield. The second
part of non-bonded interactions is handled by MiMiC within the
electrostatic embedding approach. This adds additional terms to the
Hamiltonian of the system:

   .. math::

      E_{QM/MM}^{es} = -\sum_a^{N_{mm}}Q_a\int\rho(\mathbf{r})\frac{r_{c,a}^4 
      - |\mathbf{R_a} - \mathbf{r}|^4}{r_{c,a}^5 - |\mathbf{R_a} - \mathbf{r}|^5}d\mathbf{r} 
      + \sum_a^{N_{mm}}\sum_n^{N_{qm}}Q_aZ_n
      \frac{r_{c,a}^4 - |\mathbf{R_a} - \mathbf{R_n}|^4}
      {r_{c,a}^5 - |\mathbf{R_a} - \mathbf{R_n}|^5}

where :math:`N_{mm}` is a number of MM atoms :math:`N_{qm}`, is the number of QM atoms
and :math:`r_{c,a}` is the covalent radius of the MM atoms. The first
term above corresponds to the damped Coulomb interaction between the
eletronic density :math:`\rho(\mathbf{r})` of the QM region and the MM
atoms. The damping is needed due to the fact that CPMD uses a plane-wave
basis set to expand the electronic wavefunction. Unlike localized
basis sets, plane waves are delocalized and this may give a rise to
the so-called electron spill-out problem: positively charged MM atoms
may artificially overpolarize the electronic cloud due to the absence
of quantum mechanical effects (e.g. Pauli repusion) that would normally
prevent it (in a fully quantum system). This functional form of the
damped Coulomb potential from the equation above was introduced in
\ :ref:`180 <refRoethlisbergerQMMM>`.

Since computing the integrals in the first term above can be computational
extremely expensive, MiMiC also implements hierarchical electrostatic
embedding scheme in order to mitigate the enormous computational effort
needed to compute :math:`N_{mm}` integrals over the electronic grid.
Within this scheme the MM atoms are grouped into two shells according
to the distance from the QM region: the short-ranged and long-ranged one.
For the MM atoms in the short-ranged shell the QM/MM interactions are
calculated using the equation above. In contrast to that, the interactions
involving MM atoms from the long-ranged shell are computed using
the multipolar expansion of the QM electrostatic potential.
More details about it can be found in \ :ref:`180 <refRoethlisbergerQMMM>`.


Application coupling model
^^^^^^^^^^^^^^^^^^^^^^^^^^

Unlike the majority of QM/MM interfaces, MiMiC uses a loose coupling between
partner codes. This means that instead of compiling both codes into a
single binary MiMiC builds separate executables for CPMD and |Gromacs|.
The user will then prepare the input for both codes and run them simultaneously.
Each of the codes is running using a separate pool of MPI processes and 
communicate the necessary data (e.g. coordinates, energies and forces) 
through MPI client-server mechanism. Within MiMiC framework CPMD acts 
as a server and |Gromacs| becomes the client.

Software prerequisites
^^^^^^^^^^^^^^^^^^^^^^

#. |Gromacs| version 2019+. Newer major releases may support multiple versions of
   MiMiC.
#. CPMD version 4.1+.

Usage
^^^^^

After :ref:`installing with MiMiC`, to run a MiMiC QM/MM simulation
one needs to:

#. Get and compile CPMD with MiMiC support.
#. Do a normal classical equilibration with |Gromacs|.
#. Create an index group representing QM atoms within |Gromacs|.
   Keep in mind that this group should also include link atoms
   bound to atoms in the QM region, as they have to be treated
   at quantum level.
#. Prepare input for CPMD and |Gromacs| according to the recommendations
   below.
#. Run both CPMD and |Gromacs| as two independent instances within
   a single batch job.

Preparing the input file for |Gromacs|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to setup the :ref:`mdp` file for a MiMiC simulation one needs
to add two options:

#. :mdp-value:`integrator=mimic` to enable MiMiC workflow within GROMACS.
#. ``QMMM-grps=<name_of_qm_index_group>`` to indicate all the atoms
   that are going to be handled by CPMD.

Since CPMD is going to perform the MD integration, only :ref:`mdp`
options relating to force calculation and output are active.

After setting up the :ref:`mdp` file one can run :ref:`grompp <gmx
grompp>` as usual. :ref:`grompp <gmx grompp>` will set the charges of
all the QM atoms to zero to avoid double-counting of Coulomb
interactions. Moreover, it will update non-bonded exclusion lists to
exclude LJ interactions between QM atoms (since they will be described
by CPMD). Finally, it will remove bonds between QM atoms (if
present). We recommend to output the preprocessed topology file using
``gmx grompp -pp <preprocessed_topology_file>`` as it will help to
prepare the input for CPMD in an automated way.

Preparing the input file for CPMD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This section will only describe the MiMiC-related input in CPMD - for the
configuration of a DFT-related options - please refer to the `CPMD manual
<https://www.cpmd.org/>`__.
After preparing the input for GROMACS and having obtained the
preprocessed topology file, simply run the Python
preprocessor script provided within the MiMiC distribution to obtain
MiMiC-related part of the CPMD input file. The usage of the script is simple:

::

    prepare-qmmm.py <index_file> <gro_file> <preprocessed_topology_file> <qm_group_name>

Be advised that for MiMiC it is crucial that the forcefield contains the data about
the element number of each atom type! If it does not provide it, the preprocessor
will fail with the error:

::

    It looks like the forcefield that you are using has no information about the element number.
    The element number is needed to run QM/MM simulations.

Given all the relevant information the script will print the part of the CPMD
input that is related to MiMiC. Here is the sample output with the short
descriptions of keywords that can be found in this part of CPMD input:

::

    &MIMIC
    PATHS
    1
    <some_absolute_path>
    BOX
    35.77988547402689 35.77988547402689 35.77988547402689
    OVERLAPS
    3
    2 13 1 1
    2 14 1 2
    2 15 1 3
    &END
    
    &ATOMS
    O
    1
    17.23430225802002 17.76342557295923 18.576007806615877
    H
    2
    18.557110545368047 19.086233860307257 18.727185896598506
    17.57445296048094 16.705178943080806 17.06422690678956
    &END
    Suggested QM box size [12.661165036045407, 13.71941166592383, 13.00131573850633]

``&MIMIC`` section contains MiMiC settings:

    ``PATHS`` indicates number of MM client codes involved in the simulation
    and the absolute path to each of their respective folder. Keep in mind
    that this path has to point to the folder, where |Gromacs| is going to
    be run -- otherwise it will cause a deadlock in CPMD! The next line
    contains the number of MM codes (1 in this case) and next :math:`N`
    lines contain paths to their respective working directories
    
    ``BOX`` indicates the size of the whole simulation box in Bohr in
    an ``X Y Z`` format

    ``OVERLAPS`` - sets the number and IDs of atoms within |Gromacs| that are going to be 
    treated by CPMD. The format is the following:

    ::

        <code_id> <atom_id_in_code> <host_code_id> <atom_id_in_that_code>
    
    CPMD host code id is always ID 1. Therefore, in a QM/MM simulation
    |Gromacs| will have code ID 2.

    (OPTIONAL) ``LONG-RANGE COUPLING`` - enables the faster multipole coupling for
    atoms located at a certain distance from the QM box

    (OPTIONAL) ``CUTOFF DISTANCE`` - the next line contains the cutoff for
    explicit Coulomb coupling  (20 Bohr by default if ``LONG-RANGE COUPLING``
    is present)

    (OPTIONAL) ``MULTIPOLE ORDER`` - The next line will contain the order at which
    the multipolar expansion will be truncated (default 2, maximum 20).

The ``&ATOMS`` section of CPMD input file contains all the QM atoms
within the system and has a default CPMD formatting. Please refer
to the `CPMD manual <https://www.cpmd.org/>`__
to adjust it to your needs(one will need to set the correct pseudo-potential
for each atom species).

Finally, the preprocessor suggests the size of the QM box where the electronic
density is going to be contained. The suggested value is not final
- further adjustment by user may be required.

Running a MiMiC QM/MM simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to run the simulation, one will need to run both |Gromacs| and CPMD within one job.
This is easily done within the vast majority of queueing systems. For example in
case of SLURM queue system one can use two job steps within one job. Here is
the example job script running a 242-node slurm job, allocating 2 nodes to |Gromacs|
and 240 nodes to CPMD (both codes are launched in the same folder):

::

    #!/bin/bash -x
    #SBATCH --nodes=242
    #SBATCH --output=mpi-out.%j
    #SBATCH --error=mpi-err.%j
    #SBATCH --time=00:25:00
    #SBATCH --partition=batch
    
    # *** start of job script ***

    srun -N2 --ntasks-per-node=6 --cpus-per-task=4 -r0 gmx_mpi_d mdrun -deffnm mimic -ntomp 4 &
    srun -N240 --ntasks-per-node=6 --cpus-per-task=4 -r2 cpmd.x benchmark.inp <path_to_pp_folder> > benchmark-240-4.out &
    wait


Known Issues
^^^^^^^^^^^^

OpenMPI prior to version 3.x.x has a bug preventing the usage of MiMiC
completely - please use newer versions or other MPI distributions.

With IntelMPI communication between CPMD and |Gromacs| may result
in a deadlock in some situations. If it happens, setting an
IntelMPI-related environment variable may help:

::

    export FI_OFI_RXM_USE_SRX=1
