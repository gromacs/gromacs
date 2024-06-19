/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \defgroup module_mdrun Implementation of mdrun
 * \ingroup group_mdrun
 *
 * \brief This module contains code that implements mdrun.
 */
/*! \internal \file
 *
 * \brief This file implements mdrun
 *
 * \author Berk Hess <hess@kth.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "config.h"

#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/compat/pointers.h"
#include "gromacs/domdec/options.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdrun/legacymdrunoptions.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdrun/runner.h"
#include "gromacs/mdrun/simulationcontext.h"
#include "gromacs/mdrun/simulationinputhandle.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/logging.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/unique_cptr.h"

#include "mdrun_main.h"

namespace gmx
{

int gmx_mdrun(int argc, char* argv[])
{
    // Set up the communicator, where possible (see docs for
    // SimulationContext).
    MPI_Comm                 communicator = GMX_LIB_MPI ? MPI_COMM_WORLD : MPI_COMM_NULL;
    PhysicalNodeCommunicator physicalNodeCommunicator(communicator, gmx_physicalnode_id_hash());
    std::unique_ptr<gmx_hw_info_t> hwinfo = gmx_detect_hardware(physicalNodeCommunicator, communicator);
    return gmx_mdrun(communicator, *hwinfo, argc, argv);
}

int gmx_mdrun(MPI_Comm communicator, const gmx_hw_info_t& hwinfo, int argc, char* argv[])
{
    auto mdModules = std::make_unique<MDModules>();

    std::vector<const char*> desc = {
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
        "A single log file ([TT]-g[tt]) is written.",
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
        "Running mdrun efficiently in parallel is a complex topic,",
        "many aspects of which are covered in the online User Guide. You",
        "should look there for practical advice on using many of the options",
        "available in mdrun.[PAR]",
        "ED (essential dynamics) sampling and/or additional flooding potentials",
        "are switched on by using the [TT]-ei[tt] flag followed by an [REF].edi[ref]",
        "file. The [REF].edi[ref] file can be produced with the [TT]make_edi[tt] tool",
        "or by using options in the essdyn menu of the WHAT IF program.",
        "[TT]mdrun[tt] produces a [REF].xvg[ref] output file that",
        "contains projections of positions, velocities and forces onto selected",
        "eigenvectors.[PAR]",
        "When user-defined potential functions have been selected in the",
        "[REF].mdp[ref] file the [TT]-table[tt] option is used to pass [TT]mdrun[tt]",
        "a formatted table with potential functions. The file is read from",
        "either the current directory or from the [TT]GMXLIB[tt] directory.",
        "A number of pre-formatted tables are presented in the [TT]GMXLIB[tt] dir,",
        "for 6-8, 6-9, 6-10, 6-11, 6-12 Lennard-Jones potentials with",
        "normal Coulomb.",
        "When pair interactions are present, a separate table for pair interaction",
        "functions is read using the [TT]-tablep[tt] option.[PAR]",
        "When tabulated bonded functions are present in the topology,",
        "interaction functions are read using the [TT]-tableb[tt] option.",
        "For each different tabulated interaction type used, a table file name must",
        "be given. For the topology to work, a file name given here must match a",
        "character sequence before the file extension. That sequence is: an underscore,",
        "then a 'b' for bonds, an 'a' for angles or a 'd' for dihedrals,",
        "and finally the matching table number index used in the topology. Note that,",
        "these options are deprecated, and in future will be available via grompp.[PAR]",
        "The options [TT]-px[tt] and [TT]-pf[tt] are used for writing pull COM",
        "coordinates and forces when pulling is selected",
        "in the [REF].mdp[ref] file.",
        "[PAR]",
        "The option [TT]-membed[tt] does what used to be g_membed, i.e. embed",
        "a protein into a membrane. This module requires a number of settings",
        "that are provided in a data file that is the argument of this option.",
        "For more details in membrane embedding, see the documentation in the",
        "user guide. The options [TT]-mn[tt] and [TT]-mp[tt] are used to provide",
        "the index and topology files used for the embedding.",
        "[PAR]",
        "The option [TT]-pforce[tt] is useful when you suspect a simulation",
        "crashes due to too large forces. With this option coordinates and",
        "forces of atoms with a force larger than a certain value will",
        "be printed to stderr. It will also terminate the run when non-finite",
        "forces are present.",
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
        "if no checkpoint file is found, GROMACS just assumes a normal run and",
        "starts from the first step of the [REF].tpr[ref] file. By default the output",
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
        "exceeds [TT]-maxh[tt]\\*0.99 hours. This option is particularly useful in",
        "combination with setting [TT]nsteps[tt] to -1 either in the mdp or using the",
        "similarly named command line option (although the latter is deprecated).",
        "This results in an infinite run,",
        "terminated only when the time limit set by [TT]-maxh[tt] is reached (if any)",
        "or upon receiving a signal.",
        "[PAR]",
        "Interactive molecular dynamics (IMD) can be activated by using at least one",
        "of the three IMD switches: The [TT]-imdterm[tt] switch allows one to terminate",
        "the simulation from the molecular viewer (e.g. VMD). With [TT]-imdwait[tt],",
        "[TT]mdrun[tt] pauses whenever no IMD client is connected. Pulling from the",
        "IMD remote can be turned on by [TT]-imdpull[tt].",
        "The port [TT]mdrun[tt] listens to can be altered by [TT]-imdport[tt].The",
        "file pointed to by [TT]-if[tt] contains atom indices and forces if IMD",
        "pulling is used."
    };

    LegacyMdrunOptions options;

    if (options.updateFromCommandLine(argc, argv, desc) == 0)
    {
        return 0;
    }

    ArrayRef<const std::string> multiSimDirectoryNames =
            opt2fnsIfOptionSet("-multidir", gmx::ssize(options.filenames), options.filenames.data());

    // The SimulationContext is necessary with gmxapi so that
    // resources owned by the client code can have suitable
    // lifetime. The gmx wrapper binary uses the same infrastructure,
    // but the lifetime is now trivially that of the invocation of the
    // wrapper binary.
    SimulationContext simulationContext(communicator, multiSimDirectoryNames);

    StartingBehavior startingBehavior        = StartingBehavior::NewSimulation;
    LogFilePtr       logFileGuard            = nullptr;
    gmx_multisim_t*  ms                      = simulationContext.multiSimulation_.get();
    std::tie(startingBehavior, logFileGuard) = handleRestart(findIsSimulationMainRank(ms, communicator),
                                                             communicator,
                                                             ms,
                                                             options.mdrunOptions.appendingBehavior,
                                                             gmx::ssize(options.filenames),
                                                             options.filenames.data());

    /* The named components for the builder exposed here are descriptive of the
     * state of mdrun at implementation and are not intended to be prescriptive
     * of future design. (Note the ICommandLineOptions... framework used elsewhere.)
     * The modules should ultimately take part in composing the Director code
     * for an extensible Builder.
     *
     * In the near term, we assume that resources like domain decomposition and
     * neighbor lists must be reinitialized between simulation segments.
     * We would prefer to rebuild resources only as necessary, but we defer such
     * details to future optimizations.
     */
    auto builder = MdrunnerBuilder(std::move(mdModules),
                                   compat::not_null<SimulationContext*>(&simulationContext));
    builder.addHardwareDetectionResult(&hwinfo);
    builder.addSimulationMethod(options.mdrunOptions, options.pforce, startingBehavior);
    builder.addDomainDecomposition(options.domdecOptions);
    // \todo pass by value
    builder.addNonBonded(options.nbpu_opt_choices[0]);
    // \todo pass by value
    builder.addElectrostatics(options.pme_opt_choices[0], options.pme_fft_opt_choices[0]);
    builder.addBondedTaskAssignment(options.bonded_opt_choices[0]);
    builder.addUpdateTaskAssignment(options.update_opt_choices[0]);
    builder.addNeighborList(options.nstlist_cmdline);
    builder.addReplicaExchange(options.replExParams);
    // Need to establish run-time values from various inputs to provide a resource handle to Mdrunner
    builder.addHardwareOptions(options.hw_opt);
    // \todo File names are parameters that should be managed modularly through further factoring.
    builder.addFilenames(options.filenames);
    builder.addInput(makeSimulationInput(options));
    // Note: The gmx_output_env_t life time is not managed after the call to parse_common_args.
    // \todo Implement lifetime management for gmx_output_env_t.
    // \todo Output environment should be configured outside of Mdrunner and provided as a resource.
    builder.addOutputEnvironment(options.oenv);
    builder.addLogFile(logFileGuard.get());

    auto runner = builder.build();

    return runner.mdrunner();
}

} // namespace gmx
