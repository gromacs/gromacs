/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "mdrun_main.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>

#include "gromacs/legacyheaders/checkpoint.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/main.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filenm.h"

int gmx_mdrun(int argc, char *argv[])
{
    const char   *desc[] = {
        "[THISMODULE] is the main computational chemistry engine",
        "within GROMACS. Obviously, it performs Molecular Dynamics simulations,",
        "but it can also perform Stochastic Dynamics, Energy Minimization,",
        "test particle insertion or (re)calculation of energies.",
        "Normal mode analysis is another option. In this case [TT]mdrun[tt]",
        "builds a Hessian matrix from single conformation.",
        "For usual Normal Modes-like calculations, make sure that",
        "the structure provided is properly energy-minimized.",
        "The generated matrix can be diagonalized by [gmx-nmeig].[PAR]",
        "The [TT]mdrun[tt] program reads the run input file ([TT]-s[tt])",
        "and distributes the topology over ranks if needed.",
        "[TT]mdrun[tt] produces at least four output files.",
        "A single log file ([TT]-g[tt]) is written, unless the option",
        "[TT]-seppot[tt] is used, in which case each rank writes a log file.",
        "The trajectory file ([TT]-o[tt]), contains coordinates, velocities and",
        "optionally forces.",
        "The structure file ([TT]-c[tt]) contains the coordinates and",
        "velocities of the last step.",
        "The energy file ([TT]-e[tt]) contains energies, the temperature,",
        "pressure, etc, a lot of these things are also printed in the log file.",
        "Optionally coordinates can be written to a compressed trajectory file",
        "([TT]-x[tt]).[PAR]",
        "The option [TT]-dhdl[tt] is only used when free energy calculation is",
        "turned on.[PAR]",
        "A simulation can be run in parallel using two different parallelization",
        "schemes: MPI parallelization and/or OpenMP thread parallelization.",
        "The MPI parallelization uses multiple processes when [TT]mdrun[tt] is",
        "compiled with a normal MPI library or threads when [TT]mdrun[tt] is",
        "compiled with the GROMACS built-in thread-MPI library. OpenMP threads",
        "are supported when [TT]mdrun[tt] is compiled with OpenMP. Full OpenMP support",
        "is only available with the Verlet cut-off scheme, with the (older)",
        "group scheme only PME-only ranks can use OpenMP parallelization.",
        "In all cases [TT]mdrun[tt] will by default try to use all the available",
        "hardware resources. With a normal MPI library only the options",
        "[TT]-ntomp[tt] (with the Verlet cut-off scheme) and [TT]-ntomp_pme[tt],",
        "for PME-only ranks, can be used to control the number of threads.",
        "With thread-MPI there are additional options [TT]-nt[tt], which sets",
        "the total number of threads, and [TT]-ntmpi[tt], which sets the number",
        "of thread-MPI threads.",
        "The number of OpenMP threads used by [TT]mdrun[tt] can also be set with",
        "the standard environment variable, [TT]OMP_NUM_THREADS[tt].",
        "The [TT]GMX_PME_NUM_THREADS[tt] environment variable can be used to specify",
        "the number of threads used by the PME-only ranks.[PAR]",
        "Note that combined MPI+OpenMP parallelization is in many cases",
        "slower than either on its own. However, at high parallelization, using the",
        "combination is often beneficial as it reduces the number of domains and/or",
        "the number of MPI ranks. (Less and larger domains can improve scaling,",
        "with separate PME ranks, using fewer MPI ranks reduces communication costs.)",
        "OpenMP-only parallelization is typically faster than MPI-only parallelization",
        "on a single CPU(-die). Since we currently don't have proper hardware",
        "topology detection, [TT]mdrun[tt] compiled with thread-MPI will only",
        "automatically use OpenMP-only parallelization when you use up to 4",
        "threads, up to 12 threads with Intel Nehalem/Westmere, or up to 16",
        "threads with Intel Sandy Bridge or newer CPUs. Otherwise MPI-only",
        "parallelization is used (except with GPUs, see below).",
        "[PAR]",
        "To quickly test the performance of the new Verlet cut-off scheme",
        "with old [TT].tpr[tt] files, either on CPUs or CPUs+GPUs, you can use",
        "the [TT]-testverlet[tt] option. This should not be used for production,",
        "since it can slightly modify potentials and it will remove charge groups",
        "making analysis difficult, as the [TT].tpr[tt] file will still contain",
        "charge groups. For production simulations it is highly recommended",
        "to specify [TT]cutoff-scheme = Verlet[tt] in the [TT].mdp[tt] file.",
        "[PAR]",
        "With GPUs (only supported with the Verlet cut-off scheme), the number",
        "of GPUs should match the number of particle-particle ranks, i.e.",
        "excluding PME-only ranks. With thread-MPI, unless set on the command line, the number",
        "of MPI threads will automatically be set to the number of GPUs detected.",
        "To use a subset of the available GPUs, or to manually provide a mapping of",
        "GPUs to PP ranks, you can use the [TT]-gpu_id[tt] option. The argument of [TT]-gpu_id[tt] is",
        "a string of digits (without delimiter) representing device id-s of the GPUs to be used.",
        "For example, \"[TT]02[tt]\" specifies using GPUs 0 and 2 in the first and second PP ranks per compute node",
        "respectively. To select different sets of GPU-s",
        "on different nodes of a compute cluster, use the [TT]GMX_GPU_ID[tt] environment",
        "variable instead. The format for [TT]GMX_GPU_ID[tt] is identical to ",
        "[TT]-gpu_id[tt], with the difference that an environment variable can have",
        "different values on different compute nodes. Multiple MPI ranks on each node",
        "can share GPUs. This is accomplished by specifying the id(s) of the GPU(s)",
        "multiple times, e.g. \"[TT]0011[tt]\" for four ranks sharing two GPUs in this node.",
        "This works within a single simulation, or a multi-simulation, with any form of MPI.",
        "[PAR]",
        "With the Verlet cut-off scheme and verlet-buffer-tolerance set,",
        "the pair-list update interval nstlist can be chosen freely with",
        "the option [TT]-nstlist[tt]. [TT]mdrun[tt] will then adjust",
        "the pair-list cut-off to maintain accuracy, and not adjust nstlist.",
        "Otherwise, by default, [TT]mdrun[tt] will try to increase the",
        "value of nstlist set in the [TT].mdp[tt] file to improve the",
        "performance. For CPU-only runs, nstlist might increase to 20, for",
        "GPU runs up to 40. For medium to high parallelization or with",
        "fast GPUs, a (user-supplied) larger nstlist value can give much",
        "better performance.",
        "[PAR]",
        "When using PME with separate PME ranks or with a GPU, the two major",
        "compute tasks, the non-bonded force calculation and the PME calculation",
        "run on different compute resources. If this load is not balanced,",
        "some of the resources will be idle part of time. With the Verlet",
        "cut-off scheme this load is automatically balanced when the PME load",
        "is too high (but not when it is too low). This is done by scaling",
        "the Coulomb cut-off and PME grid spacing by the same amount. In the first",
        "few hundred steps different settings are tried and the fastest is chosen",
        "for the rest of the simulation. This does not affect the accuracy of",
        "the results, but it does affect the decomposition of the Coulomb energy",
        "into particle and mesh contributions. The auto-tuning can be turned off",
        "with the option [TT]-notunepme[tt].",
        "[PAR]",
        "[TT]mdrun[tt] pins (sets affinity of) threads to specific cores,",
        "when all (logical) cores on a compute node are used by [TT]mdrun[tt],",
        "even when no multi-threading is used,",
        "as this usually results in significantly better performance.",
        "If the queuing systems or the OpenMP library pinned threads, we honor",
        "this and don't pin again, even though the layout may be sub-optimal.",
        "If you want to have [TT]mdrun[tt] override an already set thread affinity",
        "or pin threads when using less cores, use [TT]-pin on[tt].",
        "With SMT (simultaneous multithreading), e.g. Intel Hyper-Threading,",
        "there are multiple logical cores per physical core.",
        "The option [TT]-pinstride[tt] sets the stride in logical cores for",
        "pinning consecutive threads. Without SMT, 1 is usually the best choice.",
        "With Intel Hyper-Threading 2 is best when using half or less of the",
        "logical cores, 1 otherwise. The default value of 0 do exactly that:",
        "it minimizes the threads per logical core, to optimize performance.",
        "If you want to run multiple [TT]mdrun[tt] jobs on the same physical node,"
        "you should set [TT]-pinstride[tt] to 1 when using all logical cores.",
        "When running multiple [TT]mdrun[tt] (or other) simulations on the same physical",
        "node, some simulations need to start pinning from a non-zero core",
        "to avoid overloading cores; with [TT]-pinoffset[tt] you can specify",
        "the offset in logical cores for pinning.",
        "[PAR]",
        "When [TT]mdrun[tt] is started with more than 1 rank,",
        "parallelization with domain decomposition is used.",
        "[PAR]",
        "With domain decomposition, the spatial decomposition can be set",
        "with option [TT]-dd[tt]. By default [TT]mdrun[tt] selects a good decomposition.",
        "The user only needs to change this when the system is very inhomogeneous.",
        "Dynamic load balancing is set with the option [TT]-dlb[tt],",
        "which can give a significant performance improvement,",
        "especially for inhomogeneous systems. The only disadvantage of",
        "dynamic load balancing is that runs are no longer binary reproducible,",
        "but in most cases this is not important.",
        "By default the dynamic load balancing is automatically turned on",
        "when the measured performance loss due to load imbalance is 5% or more.",
        "At low parallelization these are the only important options",
        "for domain decomposition.",
        "At high parallelization the options in the next two sections",
        "could be important for increasing the performace.",
        "[PAR]",
        "When PME is used with domain decomposition, separate ranks can",
        "be assigned to do only the PME mesh calculation;",
        "this is computationally more efficient starting at about 12 ranks,",
        "or even fewer when OpenMP parallelization is used.",
        "The number of PME ranks is set with option [TT]-npme[tt],",
        "but this cannot be more than half of the ranks.",
        "By default [TT]mdrun[tt] makes a guess for the number of PME",
        "ranks when the number of ranks is larger than 16. With GPUs,",
        "using separate PME ranks is not selected automatically,",
        "since the optimal setup depends very much on the details",
        "of the hardware. In all cases, you might gain performance",
        "by optimizing [TT]-npme[tt]. Performance statistics on this issue",
        "are written at the end of the log file.",
        "For good load balancing at high parallelization, the PME grid x and y",
        "dimensions should be divisible by the number of PME ranks",
        "(the simulation will run correctly also when this is not the case).",
        "[PAR]",
        "This section lists all options that affect the domain decomposition.",
        "[PAR]",
        "Option [TT]-rdd[tt] can be used to set the required maximum distance",
        "for inter charge-group bonded interactions.",
        "Communication for two-body bonded interactions below the non-bonded",
        "cut-off distance always comes for free with the non-bonded communication.",
        "Atoms beyond the non-bonded cut-off are only communicated when they have",
        "missing bonded interactions; this means that the extra cost is minor",
        "and nearly indepedent of the value of [TT]-rdd[tt].",
        "With dynamic load balancing option [TT]-rdd[tt] also sets",
        "the lower limit for the domain decomposition cell sizes.",
        "By default [TT]-rdd[tt] is determined by [TT]mdrun[tt] based on",
        "the initial coordinates. The chosen value will be a balance",
        "between interaction range and communication cost.",
        "[PAR]",
        "When inter charge-group bonded interactions are beyond",
        "the bonded cut-off distance, [TT]mdrun[tt] terminates with an error message.",
        "For pair interactions and tabulated bonds",
        "that do not generate exclusions, this check can be turned off",
        "with the option [TT]-noddcheck[tt].",
        "[PAR]",
        "When constraints are present, option [TT]-rcon[tt] influences",
        "the cell size limit as well.",
        "Atoms connected by NC constraints, where NC is the LINCS order plus 1,",
        "should not be beyond the smallest cell size. A error message is",
        "generated when this happens and the user should change the decomposition",
        "or decrease the LINCS order and increase the number of LINCS iterations.",
        "By default [TT]mdrun[tt] estimates the minimum cell size required for P-LINCS",
        "in a conservative fashion. For high parallelization it can be useful",
        "to set the distance required for P-LINCS with the option [TT]-rcon[tt].",
        "[PAR]",
        "The [TT]-dds[tt] option sets the minimum allowed x, y and/or z scaling",
        "of the cells with dynamic load balancing. [TT]mdrun[tt] will ensure that",
        "the cells can scale down by at least this factor. This option is used",
        "for the automated spatial decomposition (when not using [TT]-dd[tt])",
        "as well as for determining the number of grid pulses, which in turn",
        "sets the minimum allowed cell size. Under certain circumstances",
        "the value of [TT]-dds[tt] might need to be adjusted to account for",
        "high or low spatial inhomogeneity of the system.",
        "[PAR]",
        "The option [TT]-gcom[tt] can be used to only do global communication",
        "every n steps.",
        "This can improve performance for highly parallel simulations",
        "where this global communication step becomes the bottleneck.",
        "For a global thermostat and/or barostat the temperature",
        "and/or pressure will also only be updated every [TT]-gcom[tt] steps.",
        "By default it is set to the minimum of nstcalcenergy and nstlist.[PAR]",
        "With [TT]-rerun[tt] an input trajectory can be given for which ",
        "forces and energies will be (re)calculated. Neighbor searching will be",
        "performed for every frame, unless [TT]nstlist[tt] is zero",
        "(see the [TT].mdp[tt] file).[PAR]",
        "ED (essential dynamics) sampling and/or additional flooding potentials",
        "are switched on by using the [TT]-ei[tt] flag followed by an [TT].edi[tt]",
        "file. The [TT].edi[tt] file can be produced with the [TT]make_edi[tt] tool",
        "or by using options in the essdyn menu of the WHAT IF program.",
        "[TT]mdrun[tt] produces a [TT].xvg[tt] output file that",
        "contains projections of positions, velocities and forces onto selected",
        "eigenvectors.[PAR]",
        "When user-defined potential functions have been selected in the",
        "[TT].mdp[tt] file the [TT]-table[tt] option is used to pass [TT]mdrun[tt]",
        "a formatted table with potential functions. The file is read from",
        "either the current directory or from the [TT]GMXLIB[tt] directory.",
        "A number of pre-formatted tables are presented in the [TT]GMXLIB[tt] dir,",
        "for 6-8, 6-9, 6-10, 6-11, 6-12 Lennard-Jones potentials with",
        "normal Coulomb.",
        "When pair interactions are present, a separate table for pair interaction",
        "functions is read using the [TT]-tablep[tt] option.[PAR]",
        "When tabulated bonded functions are present in the topology,",
        "interaction functions are read using the [TT]-tableb[tt] option.",
        "For each different tabulated interaction type the table file name is",
        "modified in a different way: before the file extension an underscore is",
        "appended, then a 'b' for bonds, an 'a' for angles or a 'd' for dihedrals",
        "and finally the table number of the interaction type.[PAR]",
        "The options [TT]-px[tt] and [TT]-pf[tt] are used for writing pull COM",
        "coordinates and forces when pulling is selected",
        "in the [TT].mdp[tt] file.[PAR]",
        "With [TT]-multi[tt] or [TT]-multidir[tt], multiple systems can be ",
        "simulated in parallel.",
        "As many input files/directories are required as the number of systems. ",
        "The [TT]-multidir[tt] option takes a list of directories (one for each ",
        "system) and runs in each of them, using the input/output file names, ",
        "such as specified by e.g. the [TT]-s[tt] option, relative to these ",
        "directories.",
        "With [TT]-multi[tt], the system number is appended to the run input ",
        "and each output filename, for instance [TT]topol.tpr[tt] becomes",
        "[TT]topol0.tpr[tt], [TT]topol1.tpr[tt] etc.",
        "The number of ranks per system is the total number of ranks",
        "divided by the number of systems.",
        "One use of this option is for NMR refinement: when distance",
        "or orientation restraints are present these can be ensemble averaged",
        "over all the systems.[PAR]",
        "With [TT]-replex[tt] replica exchange is attempted every given number",
        "of steps. The number of replicas is set with the [TT]-multi[tt] or ",
        "[TT]-multidir[tt] option, described above.",
        "All run input files should use a different coupling temperature,",
        "the order of the files is not important. The random seed is set with",
        "[TT]-reseed[tt]. The velocities are scaled and neighbor searching",
        "is performed after every exchange.[PAR]",
        "Finally some experimental algorithms can be tested when the",
        "appropriate options have been given. Currently under",
        "investigation are: polarizability.",
        "[PAR]",
        "The option [TT]-membed[tt] does what used to be g_membed, i.e. embed",
        "a protein into a membrane. The data file should contain the options",
        "that where passed to g_membed before. The [TT]-mn[tt] and [TT]-mp[tt]",
        "both apply to this as well.",
        "[PAR]",
        "The option [TT]-pforce[tt] is useful when you suspect a simulation",
        "crashes due to too large forces. With this option coordinates and",
        "forces of atoms with a force larger than a certain value will",
        "be printed to stderr.",
        "[PAR]",
        "Checkpoints containing the complete state of the system are written",
        "at regular intervals (option [TT]-cpt[tt]) to the file [TT]-cpo[tt],",
        "unless option [TT]-cpt[tt] is set to -1.",
        "The previous checkpoint is backed up to [TT]state_prev.cpt[tt] to",
        "make sure that a recent state of the system is always available,",
        "even when the simulation is terminated while writing a checkpoint.",
        "With [TT]-cpnum[tt] all checkpoint files are kept and appended",
        "with the step number.",
        "A simulation can be continued by reading the full state from file",
        "with option [TT]-cpi[tt]. This option is intelligent in the way that",
        "if no checkpoint file is found, Gromacs just assumes a normal run and",
        "starts from the first step of the [TT].tpr[tt] file. By default the output",
        "will be appending to the existing output files. The checkpoint file",
        "contains checksums of all output files, such that you will never",
        "loose data when some output files are modified, corrupt or removed.",
        "There are three scenarios with [TT]-cpi[tt]:[PAR]",
        "[TT]*[tt] no files with matching names are present: new output files are written[PAR]",
        "[TT]*[tt] all files are present with names and checksums matching those stored",
        "in the checkpoint file: files are appended[PAR]",
        "[TT]*[tt] otherwise no files are modified and a fatal error is generated[PAR]",
        "With [TT]-noappend[tt] new output files are opened and the simulation",
        "part number is added to all output file names.",
        "Note that in all cases the checkpoint file itself is not renamed",
        "and will be overwritten, unless its name does not match",
        "the [TT]-cpo[tt] option.",
        "[PAR]",
        "With checkpointing the output is appended to previously written",
        "output files, unless [TT]-noappend[tt] is used or none of the previous",
        "output files are present (except for the checkpoint file).",
        "The integrity of the files to be appended is verified using checksums",
        "which are stored in the checkpoint file. This ensures that output can",
        "not be mixed up or corrupted due to file appending. When only some",
        "of the previous output files are present, a fatal error is generated",
        "and no old output files are modified and no new output files are opened.",
        "The result with appending will be the same as from a single run.",
        "The contents will be binary identical, unless you use a different number",
        "of ranks or dynamic load balancing or the FFT library uses optimizations",
        "through timing.",
        "[PAR]",
        "With option [TT]-maxh[tt] a simulation is terminated and a checkpoint",
        "file is written at the first neighbor search step where the run time",
        "exceeds [TT]-maxh[tt]*0.99 hours.",
        "[PAR]",
        "When [TT]mdrun[tt] receives a TERM signal, it will set nsteps to the current",
        "step plus one. When [TT]mdrun[tt] receives an INT signal (e.g. when ctrl+C is",
        "pressed), it will stop after the next neighbor search step ",
        "(with nstlist=0 at the next step).",
        "In both cases all the usual output will be written to file.",
        "When running with MPI, a signal to one of the [TT]mdrun[tt] ranks",
        "is sufficient, this signal should not be sent to mpirun or",
        "the [TT]mdrun[tt] process that is the parent of the others.",
        "[PAR]",
        "Interactive molecular dynamics (IMD) can be activated by using at least one",
        "of the three IMD switches: The [TT]-imdterm[tt] switch allows to terminate the",
        "simulation from the molecular viewer (e.g. VMD). With [TT]-imdwait[tt],",
        "[TT]mdrun[tt] pauses whenever no IMD client is connected. Pulling from the",
        "IMD remote can be turned on by [TT]-imdpull[tt].",
        "The port [TT]mdrun[tt] listens to can be altered by [TT]-imdport[tt].The",
        "file pointed to by [TT]-if[tt] contains atom indices and forces if IMD",
        "pulling is used."
        "[PAR]",
        "When [TT]mdrun[tt] is started with MPI, it does not run niced by default."
    };
    t_commrec    *cr;
    t_filenm      fnm[] = {
        { efTPX, NULL,      NULL,       ffREAD },
        { efTRN, "-o",      NULL,       ffWRITE },
        { efCOMPRESSED, "-x", NULL,     ffOPTWR },
        { efCPT, "-cpi",    NULL,       ffOPTRD },
        { efCPT, "-cpo",    NULL,       ffOPTWR },
        { efSTO, "-c",      "confout",  ffWRITE },
        { efEDR, "-e",      "ener",     ffWRITE },
        { efLOG, "-g",      "md",       ffWRITE },
        { efXVG, "-dhdl",   "dhdl",     ffOPTWR },
        { efXVG, "-field",  "field",    ffOPTWR },
        { efXVG, "-table",  "table",    ffOPTRD },
        { efXVG, "-tabletf", "tabletf",    ffOPTRD },
        { efXVG, "-tablep", "tablep",   ffOPTRD },
        { efXVG, "-tableb", "table",    ffOPTRD },
        { efTRX, "-rerun",  "rerun",    ffOPTRD },
        { efXVG, "-tpi",    "tpi",      ffOPTWR },
        { efXVG, "-tpid",   "tpidist",  ffOPTWR },
        { efEDI, "-ei",     "sam",      ffOPTRD },
        { efXVG, "-eo",     "edsam",    ffOPTWR },
        { efXVG, "-devout", "deviatie", ffOPTWR },
        { efXVG, "-runav",  "runaver",  ffOPTWR },
        { efXVG, "-px",     "pullx",    ffOPTWR },
        { efXVG, "-pf",     "pullf",    ffOPTWR },
        { efXVG, "-ro",     "rotation", ffOPTWR },
        { efLOG, "-ra",     "rotangles", ffOPTWR },
        { efLOG, "-rs",     "rotslabs", ffOPTWR },
        { efLOG, "-rt",     "rottorque", ffOPTWR },
        { efMTX, "-mtx",    "nm",       ffOPTWR },
        { efNDX, "-dn",     "dipole",   ffOPTWR },
        { efRND, "-multidir", NULL,      ffOPTRDMULT},
        { efDAT, "-membed", "membed",   ffOPTRD },
        { efTOP, "-mp",     "membed",   ffOPTRD },
        { efNDX, "-mn",     "membed",   ffOPTRD },
        { efXVG, "-if",     "imdforces", ffOPTWR },
        { efXVG, "-swap",   "swapions", ffOPTWR }
    };
#define NFILE asize(fnm)

    /* Command line options ! */
    gmx_bool        bDDBondCheck  = TRUE;
    gmx_bool        bDDBondComm   = TRUE;
    gmx_bool        bTunePME      = TRUE;
    gmx_bool        bTestVerlet   = FALSE;
    gmx_bool        bVerbose      = FALSE;
    gmx_bool        bCompact      = TRUE;
    gmx_bool        bSepPot       = FALSE;
    gmx_bool        bRerunVSite   = FALSE;
    gmx_bool        bConfout      = TRUE;
    gmx_bool        bReproducible = FALSE;
    gmx_bool        bIMDwait      = FALSE;
    gmx_bool        bIMDterm      = FALSE;
    gmx_bool        bIMDpull      = FALSE;

    int             npme          = -1;
    int             nstlist       = 0;
    int             nmultisim     = 0;
    int             nstglobalcomm = -1;
    int             repl_ex_nst   = 0;
    int             repl_ex_seed  = -1;
    int             repl_ex_nex   = 0;
    int             nstepout      = 100;
    int             resetstep     = -1;
    gmx_int64_t     nsteps        = -2;   /* the value -2 means that the mdp option will be used */
    int             imdport       = 8888; /* can be almost anything, 8888 is easy to remember */

    rvec            realddxyz          = {0, 0, 0};
    const char     *ddno_opt[ddnoNR+1] =
    { NULL, "interleave", "pp_pme", "cartesian", NULL };
    const char     *dddlb_opt[] =
    { NULL, "auto", "no", "yes", NULL };
    const char     *thread_aff_opt[threadaffNR+1] =
    { NULL, "auto", "on", "off", NULL };
    const char     *nbpu_opt[] =
    { NULL, "auto", "cpu", "gpu", "gpu_cpu", NULL };
    real            rdd                   = 0.0, rconstr = 0.0, dlb_scale = 0.8, pforce = -1;
    char           *ddcsx                 = NULL, *ddcsy = NULL, *ddcsz = NULL;
    real            cpt_period            = 15.0, max_hours = -1;
    gmx_bool        bAppendFiles          = TRUE;
    gmx_bool        bKeepAndNumCPT        = FALSE;
    gmx_bool        bResetCountersHalfWay = FALSE;
    output_env_t    oenv                  = NULL;
    const char     *deviceOptions         = "";

    /* Non transparent initialization of a complex gmx_hw_opt_t struct.
     * But unfortunately we are not allowed to call a function here,
     * since declarations follow below.
     */
    gmx_hw_opt_t    hw_opt = {
        0, 0, 0, 0, threadaffSEL, 0, 0,
        { NULL, FALSE, 0, NULL }
    };

    t_pargs         pa[] = {

        { "-dd",      FALSE, etRVEC, {&realddxyz},
          "Domain decomposition grid, 0 is optimize" },
        { "-ddorder", FALSE, etENUM, {ddno_opt},
          "DD rank order" },
        { "-npme",    FALSE, etINT, {&npme},
          "Number of separate ranks to be used for PME, -1 is guess" },
        { "-nt",      FALSE, etINT, {&hw_opt.nthreads_tot},
          "Total number of threads to start (0 is guess)" },
        { "-ntmpi",   FALSE, etINT, {&hw_opt.nthreads_tmpi},
          "Number of thread-MPI threads to start (0 is guess)" },
        { "-ntomp",   FALSE, etINT, {&hw_opt.nthreads_omp},
          "Number of OpenMP threads per MPI rank to start (0 is guess)" },
        { "-ntomp_pme", FALSE, etINT, {&hw_opt.nthreads_omp_pme},
          "Number of OpenMP threads per MPI rank to start (0 is -ntomp)" },
        { "-pin",     FALSE, etENUM, {thread_aff_opt},
          "Set thread affinities" },
        { "-pinoffset", FALSE, etINT, {&hw_opt.core_pinning_offset},
          "The starting logical core number for pinning to cores; used to avoid pinning threads from different mdrun instances to the same core" },
        { "-pinstride", FALSE, etINT, {&hw_opt.core_pinning_stride},
          "Pinning distance in logical cores for threads, use 0 to minimize the number of threads per physical core" },
        { "-gpu_id",  FALSE, etSTR, {&hw_opt.gpu_opt.gpu_id},
          "List of GPU device id-s to use, specifies the per-node PP rank to GPU mapping" },
        { "-ddcheck", FALSE, etBOOL, {&bDDBondCheck},
          "Check for all bonded interactions with DD" },
        { "-ddbondcomm", FALSE, etBOOL, {&bDDBondComm},
          "HIDDENUse special bonded atom communication when [TT]-rdd[tt] > cut-off" },
        { "-rdd",     FALSE, etREAL, {&rdd},
          "The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates" },
        { "-rcon",    FALSE, etREAL, {&rconstr},
          "Maximum distance for P-LINCS (nm), 0 is estimate" },
        { "-dlb",     FALSE, etENUM, {dddlb_opt},
          "Dynamic load balancing (with DD)" },
        { "-dds",     FALSE, etREAL, {&dlb_scale},
          "Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased in order to "
          "provide a margin in which dynamic load balancing can act while preserving the minimum cell size." },
        { "-ddcsx",   FALSE, etSTR, {&ddcsx},
          "HIDDENA string containing a vector of the relative sizes in the x "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-ddcsy",   FALSE, etSTR, {&ddcsy},
          "HIDDENA string containing a vector of the relative sizes in the y "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-ddcsz",   FALSE, etSTR, {&ddcsz},
          "HIDDENA string containing a vector of the relative sizes in the z "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-gcom",    FALSE, etINT, {&nstglobalcomm},
          "Global communication frequency" },
        { "-nb",      FALSE, etENUM, {&nbpu_opt},
          "Calculate non-bonded interactions on" },
        { "-nstlist", FALSE, etINT, {&nstlist},
          "Set nstlist when using a Verlet buffer tolerance (0 is guess)" },
        { "-tunepme", FALSE, etBOOL, {&bTunePME},
          "Optimize PME load between PP/PME ranks or GPU/CPU" },
        { "-testverlet", FALSE, etBOOL, {&bTestVerlet},
          "Test the Verlet non-bonded scheme" },
        { "-v",       FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy" },
        { "-compact", FALSE, etBOOL, {&bCompact},
          "Write a compact log file" },
        { "-seppot",  FALSE, etBOOL, {&bSepPot},
          "Write separate V and dVdl terms for each interaction type and rank to the log file(s)" },
        { "-pforce",  FALSE, etREAL, {&pforce},
          "Print all forces larger than this (kJ/mol nm)" },
        { "-reprod",  FALSE, etBOOL, {&bReproducible},
          "Try to avoid optimizations that affect binary reproducibility" },
        { "-cpt",     FALSE, etREAL, {&cpt_period},
          "Checkpoint interval (minutes)" },
        { "-cpnum",   FALSE, etBOOL, {&bKeepAndNumCPT},
          "Keep and number checkpoint files" },
        { "-append",  FALSE, etBOOL, {&bAppendFiles},
          "Append to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names" },
        { "-nsteps",  FALSE, etINT64, {&nsteps},
          "Run this number of steps, overrides .mdp file option" },
        { "-maxh",   FALSE, etREAL, {&max_hours},
          "Terminate after 0.99 times this time (hours)" },
        { "-multi",   FALSE, etINT, {&nmultisim},
          "Do multiple simulations in parallel" },
        { "-replex",  FALSE, etINT, {&repl_ex_nst},
          "Attempt replica exchange periodically with this period (steps)" },
        { "-nex",  FALSE, etINT, {&repl_ex_nex},
          "Number of random exchanges to carry out each exchange interval (N^3 is one suggestion).  -nex zero or not specified gives neighbor replica exchange." },
        { "-reseed",  FALSE, etINT, {&repl_ex_seed},
          "Seed for replica exchange, -1 is generate a seed" },
        { "-imdport",    FALSE, etINT, {&imdport},
          "HIDDENIMD listening port" },
        { "-imdwait",  FALSE, etBOOL, {&bIMDwait},
          "HIDDENPause the simulation while no IMD client is connected" },
        { "-imdterm",  FALSE, etBOOL, {&bIMDterm},
          "HIDDENAllow termination of the simulation from IMD client" },
        { "-imdpull",  FALSE, etBOOL, {&bIMDpull},
          "HIDDENAllow pulling in the simulation from IMD client" },
        { "-rerunvsite", FALSE, etBOOL, {&bRerunVSite},
          "HIDDENRecalculate virtual site coordinates with [TT]-rerun[tt]" },
        { "-confout", FALSE, etBOOL, {&bConfout},
          "HIDDENWrite the last configuration with [TT]-c[tt] and force checkpointing at the last step" },
        { "-stepout", FALSE, etINT, {&nstepout},
          "HIDDENFrequency of writing the remaining wall clock time for the run" },
        { "-resetstep", FALSE, etINT, {&resetstep},
          "HIDDENReset cycle counters after these many time steps" },
        { "-resethway", FALSE, etBOOL, {&bResetCountersHalfWay},
          "HIDDENReset the cycle counters after half the number of steps or halfway [TT]-maxh[tt]" }
    };
    unsigned long   Flags, PCA_Flags;
    ivec            ddxyz;
    int             dd_node_order;
    gmx_bool        bAddPart;
    FILE           *fplog, *fpmulti;
    int             sim_part, sim_part_fn;
    const char     *part_suffix = ".part";
    char            suffix[STRLEN];
    int             rc;
    char          **multidir = NULL;


    cr = init_commrec();

    PCA_Flags = (PCA_CAN_SET_DEFFNM | (MASTER(cr) ? 0 : PCA_QUIET));

    /* Comment this in to do fexist calls only on master
     * works not with rerun or tables at the moment
     * also comment out the version of init_forcerec in md.c
     * with NULL instead of opt2fn
     */
    /*
       if (!MASTER(cr))
       {
       PCA_Flags |= PCA_NOT_READ_NODE;
       }
     */

    if (!parse_common_args(&argc, argv, PCA_Flags, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }


    /* we set these early because they might be used in init_multisystem()
       Note that there is the potential for npme>nnodes until the number of
       threads is set later on, if there's thread parallelization. That shouldn't
       lead to problems. */
    dd_node_order = nenum(ddno_opt);
    cr->npmenodes = npme;

    hw_opt.thread_affinity = nenum(thread_aff_opt);

    /* now check the -multi and -multidir option */
    if (opt2bSet("-multidir", NFILE, fnm))
    {
        if (nmultisim > 0)
        {
            gmx_fatal(FARGS, "mdrun -multi and -multidir options are mutually exclusive.");
        }
        nmultisim = opt2fns(&multidir, "-multidir", NFILE, fnm);
    }


    if (repl_ex_nst != 0 && nmultisim < 2)
    {
        gmx_fatal(FARGS, "Need at least two replicas for replica exchange (option -multi)");
    }

    if (repl_ex_nex < 0)
    {
        gmx_fatal(FARGS, "Replica exchange number of exchanges needs to be positive");
    }

    if (nmultisim >= 1)
    {
#ifndef GMX_THREAD_MPI
        gmx_bool bParFn = (multidir == NULL);
        init_multisystem(cr, nmultisim, multidir, NFILE, fnm, bParFn);
#else
        gmx_fatal(FARGS, "mdrun -multi or -multidir are not supported with the thread-MPI library. "
                  "Please compile GROMACS with a proper external MPI library.");
#endif
    }

    bAddPart = !bAppendFiles;

    /* Check if there is ANY checkpoint file available */
    sim_part    = 1;
    sim_part_fn = sim_part;
    if (opt2bSet("-cpi", NFILE, fnm))
    {
        if (bSepPot && bAppendFiles)
        {
            gmx_fatal(FARGS, "Output file appending is not supported with -seppot");
        }

        bAppendFiles =
            read_checkpoint_simulation_part(opt2fn_master("-cpi", NFILE,
                                                          fnm, cr),
                                            &sim_part_fn, NULL, cr,
                                            bAppendFiles, NFILE, fnm,
                                            part_suffix, &bAddPart);
        if (sim_part_fn == 0 && MULTIMASTER(cr))
        {
            fprintf(stdout, "No previous checkpoint file present, assuming this is a new run.\n");
        }
        else
        {
            sim_part = sim_part_fn + 1;
        }

        if (MULTISIM(cr) && MASTER(cr))
        {
            if (MULTIMASTER(cr))
            {
                /* Log file is not yet available, so if there's a
                 * problem we can only write to stderr. */
                fpmulti = stderr;
            }
            else
            {
                fpmulti = NULL;
            }
            check_multi_int(fpmulti, cr->ms, sim_part, "simulation part", TRUE);
        }
    }
    else
    {
        bAppendFiles = FALSE;
    }

    if (!bAppendFiles)
    {
        sim_part_fn = sim_part;
    }

    if (bAddPart)
    {
        /* Rename all output files (except checkpoint files) */
        /* create new part name first (zero-filled) */
        sprintf(suffix, "%s%04d", part_suffix, sim_part_fn);

        add_suffix_to_output_names(fnm, NFILE, suffix);
        if (MULTIMASTER(cr))
        {
            fprintf(stdout, "Checkpoint file is from part %d, new output files will be suffixed '%s'.\n", sim_part-1, suffix);
        }
    }

    Flags = opt2bSet("-rerun", NFILE, fnm) ? MD_RERUN : 0;
    Flags = Flags | (bSepPot       ? MD_SEPPOT       : 0);
    Flags = Flags | (bDDBondCheck  ? MD_DDBONDCHECK  : 0);
    Flags = Flags | (bDDBondComm   ? MD_DDBONDCOMM   : 0);
    Flags = Flags | (bTunePME      ? MD_TUNEPME      : 0);
    Flags = Flags | (bTestVerlet   ? MD_TESTVERLET   : 0);
    Flags = Flags | (bConfout      ? MD_CONFOUT      : 0);
    Flags = Flags | (bRerunVSite   ? MD_RERUN_VSITE  : 0);
    Flags = Flags | (bReproducible ? MD_REPRODUCIBLE : 0);
    Flags = Flags | (bAppendFiles  ? MD_APPENDFILES  : 0);
    Flags = Flags | (opt2parg_bSet("-append", asize(pa), pa) ? MD_APPENDFILESSET : 0);
    Flags = Flags | (bKeepAndNumCPT ? MD_KEEPANDNUMCPT : 0);
    Flags = Flags | (sim_part > 1    ? MD_STARTFROMCPT : 0);
    Flags = Flags | (bResetCountersHalfWay ? MD_RESETCOUNTERSHALFWAY : 0);
    Flags = Flags | (bIMDwait      ? MD_IMDWAIT      : 0);
    Flags = Flags | (bIMDterm      ? MD_IMDTERM      : 0);
    Flags = Flags | (bIMDpull      ? MD_IMDPULL      : 0);

    /* We postpone opening the log file if we are appending, so we can
       first truncate the old log file and append to the correct position
       there instead.  */
    if ((MASTER(cr) || bSepPot) && !bAppendFiles)
    {
        gmx_log_open(ftp2fn(efLOG, NFILE, fnm), cr,
                     !bSepPot, Flags & MD_APPENDFILES, &fplog);
        please_cite(fplog, "Abraham2015");
        please_cite(fplog, "Pall2015");
        please_cite(fplog, "Pronk2013");
        please_cite(fplog, "Hess2008b");
        please_cite(fplog, "Spoel2005a");
        please_cite(fplog, "Lindahl2001a");
        please_cite(fplog, "Berendsen95a");
    }
    else if (!MASTER(cr) && bSepPot)
    {
        gmx_log_open(ftp2fn(efLOG, NFILE, fnm), cr, !bSepPot, Flags, &fplog);
    }
    else
    {
        fplog = NULL;
    }

    ddxyz[XX] = (int)(realddxyz[XX] + 0.5);
    ddxyz[YY] = (int)(realddxyz[YY] + 0.5);
    ddxyz[ZZ] = (int)(realddxyz[ZZ] + 0.5);

    rc = mdrunner(&hw_opt, fplog, cr, NFILE, fnm, oenv, bVerbose, bCompact,
                  nstglobalcomm, ddxyz, dd_node_order, rdd, rconstr,
                  dddlb_opt[0], dlb_scale, ddcsx, ddcsy, ddcsz,
                  nbpu_opt[0], nstlist,
                  nsteps, nstepout, resetstep,
                  nmultisim, repl_ex_nst, repl_ex_nex, repl_ex_seed,
                  pforce, cpt_period, max_hours, deviceOptions, imdport, Flags);

    /* Log file has to be closed in mdrunner if we are appending to it
       (fplog not set here) */
    if (MASTER(cr) && !bAppendFiles)
    {
        gmx_log_close(fplog);
    }

    return rc;
}
