MiMiC Hybrid Quantum Mechanical/Molecular Mechanical simulations
----------------------------------------------------------------

This section desctribes the coupling to a novel MiMiC QM/MM interface
that enables the usage of CPMD QM code to drive the simulation.
Within QM/MM approach part of the system (e.g. active site of an enzyme)
is being treated by the QM level of theory (as we cannot neglect electronic
degrees of freedom while descibing some processes e.g.  chemical 
reactions). In the same time the rest of the system (remainder of the 
protein, solvent, etc.) is treated with GROMACS.

Overview
^^^^^^^^
MiMiC implements the coupling scheme developed by the group of Prof. U. Roethlisberger
described in \ :ref:`180 <refRoethlisbergerQMMM>`. This scheme 
uses an additive scheme with electrostatic embedding
of the classical system into the quantum Hamiltonian. The total QM/MM energy 
is calculated as a sum of subsytem contributions with added crossterms:

   .. math::

      E_{tot} = E_{QM}
      E_{QM}+E_{MM}+E_{QM/MM}

The QM contribution is computed by CPMD, MM part is processed by GROMACS 
and crossterms are treated by the MiMiC interface. Crossterms consist of
both bonded and non-bonded interactions. 

Whenever there is a chemical bond crossing the QM/MM boundary additional
care has to be taken to handle this situation correctly. Otherwise the QM
atom is left with  an unsaturated electronic orbital leading to unphysical
system behaviour. Therefore, the bond has to be capped with another QM
atom. There are two different options available in CPMD for bond capping:

#. Hydrogen capping - the simplest approach is to cap a bond with the
   hydrogen atom, constrining its relative position
   
#. Link atom pseudopotential - this strategy uses an ad-hoc pseudopotential
   developed to cap the bond. This pseudopotential would reprsent the real
   atom and, thus, will not require the bond constraint.
   
The non-bonded interactions can be separeted into two parts as well.
The first part includes electrostatic interaction. Since in MiMiC 
we use the electrostatic embedding this adds additional term to the 
Hamiltionian of the system:

   .. math::

      E_{QM/MM}^{es} = -\sum_a^{N_{mm}}Q_a\int\rho(\mathbf{r})\frac{r_{c,a}^4 
      - |\mathbf{R_a} - \mathbf{r}|^4}{r_{c,a}^5 - |\mathbf{R_a} - \mathbf{r}|^5}d\mathbf{r} 
      + \sum_a^{N_{mm}}\sum_n^{N_{qm}}Q_aZ_n
      \frac{r_{c,a}^4 - |\mathbf{R_a} - \mathbf{R_n}|^4}
      {r_{c,a}^5 - |\mathbf{R_a} - \mathbf{R_n}|^5}

Where :math:`N_{mm}` is a number of MM atoms :math:`N_{qm}` is the number of QM atoms
and :math:`r_{c,a}` is the covalent radius of the MM atoms. One can see that here the 
Coulomb interaction is damped. This is needed due to the fact that CPMD uses plane-wave 
basis set which is not localized. This may give a rise to the so-called electron 
spill-out problem. Positively charged MM atoms may overpolarize the electronic cloud
due to the absence of quantum mechanical effects (e.g. Pauli repusion) that would normally
prevent it. This form of damped Coulomb potential was introduced in
\ :ref:`180 <refRoethlisbergerQMMM>`.

MiMiC also implements hierarchical elestrostatic embedding scheme in order to mitigate
the enormous computational effort needed to compute :math:`N_mm` integrals over the electronic
grid. Within this scheme the interactions are separated into two shells: the short-ranged (where
the above-mentioned energy functional is used) and long-ranged. Within the long-ranged part
the multipolar expansion of the QM electrostatic potential is used to compute forces and energies.
More details about it can be found in \ :ref:`180 <refRoethlisbergerQMMM>`.

The second part of non-bonded interaction contains LJ energy. These are computed solely
by GROMACS.

Running model
^^^^^^^^^^^^^

Unlike the majority of QM/MM interfaces MiMiC uses a loose coupling between
partner codes. This means that instead of compiling both codes into a
single binary MiMiC builds separate executables for CPMD and GROMACS.
The user will then prepare the input for both codes and run them simultaneously.
Each of the codes is running using a separate pool of MPI processes and 
communicate the necessary data (e.g. coordinates, energies and forces) 
through MPI client-server mechanism. Within MiMiC framework CPMD acts 
as a server and GROMACS becomes the client.

Usage
^^^^^
To use MiMiC QM/MM coupling one would need to:

#. Get MiMiC Communication library can be downloaded `here
   <https://gitlab.com/MiMiC-projects/MiMiCCommLib>`__
   and add its install location to ``CMAKE_PREFIX_PATH``
#. Get CPMD with MiMiC support
#. Compile GROMACS with ``-DGMX_DOUBLE=ON -DGMX_MPI=ON -DGMX_MIMIC=ON``
#. Do a normal classical equilibration with GROMACS
#. Create an index group representing QM atoms within GROMACS
#. Prepare input for CPMD and GROMACS
#. Run both CPMD and GROMACS within one batch job

Preparing the input file for GROMACS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Since CPMD is going to do the MD itegration all :ref:`mdp` file options
controlling integrator parameters are going to be ignored. All other :ref:`mdp`
options can still be relevant and may affect your simulation (especially electrostatics
parameters)! In order to setup the :ref:`mdp` file for MiMiC simulation one needs
to add two options:

#. ``integrator = mimic`` - this will enable MiMiC workflow within GROMACS.
#. ``QMMM-grps = <name_of_qm_index_group>`` - this will indicate all the atoms
   that are going to be handled by CPMD. Keep in mind that link atoms should be
   in the group of QM atoms as they are going to be treated by CPMD

After setting up the file one can run the preprocessor as usual.
GROMACS preprocessor will zero charges of all QM atoms to avoid double-counting
of Coulomb interactions. Moreover, it will update non-bonded exclusion lists to exclude
LJ interactions between QM atoms (since they are taken care of by CPMD). Finally,
it will remove bonds between QM atoms (if present). We recommend to output also 
the preprocessed topology file using ``-pp <file_name>`` as it will help to prepare
input for CPMD in an automated way.

Preparing the input file for CPMD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This section will only touch the MiMiC related input in CPMD - for the
configuration of a DFT-related options - please refer to CPMD manual.
After preparing the input for GROMACS and having obtained the preprocessed topology
file the user can simply run the Python preprocessor script provided within 
MiMiC distribution to obtain MiMiC-related input in CPMD. The usage of the 
script is simple:

::

    prepare-qmmm.py <index_file> <gro_file> <preprocessed_topology_file> <qm_group_name>

Be advised that for MiMiC it is crucial that the forcefield contains the data about
the element number of each atom type! If it does not provide it - the preprocessor
will fail with the error.

Given all the relevant information the script will generate the part of the CPMD
input that is related to MiMiC. The sample output can be found here:

::

    &MIMIC
    PATHS
    1
    <some_absoulte_path>
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

    ``PATHS`` indicates number of MM client codes involved in the simulation and the absolute
    path to each of their respective folder. Keep in mind that this path has to point
    to the folder, where GROMACS is going to be run - otherwize it will cause a deadlock in CPMD!
    The next line contains the number of 
    MM codes (1 in this case) and next math:`N` lines contain paths to the respective folders
    
    ``BOX`` indicates the size of the whole simulation box in Bohr

    ``OVERLAPS`` - sets the number and IDs of atoms within GROMACS that are going to be 
    treated by CPMD. The format is the following:

    ::

        <code_id> <atom_id_in_code> <host_code_id> <atom_id_in_that_code>
    
    CPMD will always have ID 1 and GROMACS will have ID 2!.

    (OPTIONAL)``LONG-RANGE COUPLING`` - enables the faster multipole coupling for
    atom located at a certain distance from the QM box (20 Bohr by default)

    (OPTIONAL)``CUTOFF DISTANCE`` - the next line will contain the cutoff for
    explicit Coulomb coupling

    (OPTIONAL)``MULTIPOLE ORDER`` - The next line will contain the order at which
    the multipolar exansion will be truncated (default 2, maximum 20).

``&ATOMS`` section of CPMD input contains all the QM atoms within the system
and has a default CPMD formatting. Please refer to CPMD manual to adjust it to
your needs(one will need to set the correct pseudopotential for each atom species).

Finally, the preprocessor suggests the size of the QM box where the electronic
density is going to be contained. The choice is not final - further adjustment by
user may be required.

Running MiMiC QM/MM simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to run the simulation one will need to run both GROAMCS and CPMD within one job.
This is easily done within the vast majority of queueing systems. For example in
case of SLURM queue system one can use two job steps within one job. Here is
the example job script running a 242-node slurm job, allocating 2 nodes to GROMACS
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
    srun -N240 --ntasks-per-node=6 --cpus-per-task=4 -r2 cpmd.x benchmark.inp /homea/ias-5/bolnykh/PP > benchmark-240-4.out &
    wait


Known Issues
^^^^^^^^^^^^

OpenMPI prior to version 3.x.x has a bug preventing usage of MiMiC completely - please use
newer versions or other MPI distributions.

With IntelMPI communication between CPMD and GROMACS may result in a deadlock in
some situations. The way to avoid it is to use IntelMPI-related environment variable:
``export FI_OFI_RXM_USE_SRX=1``
