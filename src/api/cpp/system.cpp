#include <stdio.h>

#include <array>

#include "gmxapi/gmxapi.h"
#include "gmxapi/md.h"
#include "gmxapi/runner.h"
#include "gmxapi/system.h"
#include "system-builder.h"

//#include "atoms.h"
#include "gromacs/compat/make_unique.h"
#include "md-impl.h"
#include "system-impl.h"
#include "gromacs/utility.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/smalloc.h"
#include "programs/mdrun/repl_ex.h"
#include "programs/mdrun/runner.h"
#include "api/cpp/gmxapi/md/runnerstate.h"

namespace gmxapi
{

System::Impl::Impl() :
    md_ {std::make_shared<EmptyMD>()},
runner_ {
    std::make_shared<EmptyMDRunnerState>()
}
{
}

System::Impl::~Impl() = default;

std::shared_ptr<MDEngine> System::Impl::md()
{
    assert(md_ != nullptr);
    return md_;
}

void System::Impl::md(std::shared_ptr<MDEngine> md)
{
    md_ = std::move(md);
}

std::shared_ptr<IMDRunner> &System::Impl::runner()
{
    assert(runner_ != nullptr);
    return runner_;
}

void System::Impl::runner(std::shared_ptr<IMDRunner> runner)
{
    runner_ = std::move(runner);
}

// Constructor and destructor needs to be defined after Impl is defined so that we can
// use unique_ptr
System::System() :
    impl_ {gmx::compat::make_unique<System::Impl>()}
{}

System::~System() = default;

std::shared_ptr<MDEngine> System::md()
{
    assert(impl_ != nullptr);
    return impl_->md();
}

void System::md(std::shared_ptr<MDEngine> md)
{
    impl_->md(std::move(md));
}

std::shared_ptr<IMDRunner> System::runner()
{
    assert(impl_ != nullptr);
    return impl_->runner();
}

void System::runner(std::shared_ptr<IMDRunner> runner)
{
    impl_->runner(std::move(runner));
}

System::System(System &&) noexcept = default;

System &System::operator=(System &&) noexcept = default;


System::Builder::Builder() :
    system_ {gmx::compat::make_unique<System>()}
{}

System::Builder::~Builder() = default;

//System::Builder& System::Builder::defaultContext()
//{
//    return *this;
//}
//
//System::Builder& System::Builder::runner(const MDInput& inputrec)
//{
//    auto mdModules = gmx::compat::make_unique<gmx::MDModules>();
//
//    mdModules->assignOptionsToModules(inputrec.params(), nullptr);
//
//    return *this;
//}

//System::Builder& System::Builder::structure(const MDInput& inputrec)
//{
//    auto atoms = gmx::compat::make_unique<gmxapi::Atoms>(*inputrec.state());
//
//    system_->impl_->setAtoms(*atoms);
//    return *this;
//}

//System::Builder& System::Builder::topology(std::shared_ptr<Topology> topology)
//{
//
//    return *this;
//}

System::Builder &System::Builder::mdEngine(std::shared_ptr<MDEngine> md)
{
    system_->md(std::move(md));
    return *this;
}

std::unique_ptr<System> System::Builder::build()
{
    // To reproduce the CLI control flow:
    // Context:
    //  * configure environment
    //  * initialize communicator
    //  * open log file
    // IRunner:
    //  * prepare mdrunner_arglist
    //  * configure mdModules
    //  * configure parallelism
    // MDEngine:
    //  * prepare arguments for integrator_t function

    // if (system_ == nullptr)
    // {
    //     throw();
    // };
    std::unique_ptr<System> product {
        nullptr
    };
    product.swap(system_);
    assert(product != nullptr);
    return product;
}

System::Builder &System::Builder::runner(std::shared_ptr<IMDRunner> runner)
{
    system_->runner(std::move(runner));
    return *this;
}

std::unique_ptr<gmxapi::System> fromTprFile(std::string filename)
{
    // The TPR file has enough information for us to
    //  1. choose an MD engine
    //  2. Get structure information
    //  3. Get topology information
    //  4. Get a lot of simulation and runtime parameters, but not all.
    // It does not have enough information on its own to determine much about the
    // necessary computation environment. That comes from environment
    // introspection and user runtime options.

    // for what it's worth, we can choose/configure a builder based
    // on the sort of system we are building.

    // For now have very limited execution environment abstraction
    // or flexibility.
//    builder->defaultContext();

    // Build MDEngine member
    auto mdBuilder = MDProxy().builder();
    assert(mdBuilder != nullptr);
    std::shared_ptr<MDEngine> md = mdBuilder->build();
    assert(md != nullptr);
    assert(!md->info().empty());

    auto runnerBuilder = UninitializedMDRunnerState::Builder();
    runnerBuilder.mdEngine(md);
    auto tpxState = gmx::TpxState::initializeFromFile(filename);
    assert(!tpxState->isDirty());
    assert(tpxState->isInitialized());
    runnerBuilder.tpxState(std::move(tpxState));
    auto runner = runnerBuilder.build();
    assert(runner != nullptr);

    // One way to get around the required CLI argument is to make sure topol.tpr is in the working directory...


    auto builder = gmx::compat::make_unique<System::Builder>();
    builder->mdEngine(md);
    builder->runner(std::move(runner));

    auto system = builder->build();

    // initialize libgromacs. parameters are passed to MPI initialization funcs.
    // gmx::init(int *argc, char ***argv)
    // noop if no MPI compiled in
    // TODO: handle MPI, avoid bare pointers in external API
//    gmx::init(nullptr, nullptr);

    //
    // Call broadcastArguments(int *argc, char ***argv)
    // returns immediately if gmx_node_num() <= 1
    // otherwise rewrites argv on non-master ranks
    // TODO: handle MPI

    // g_commandLineContext.reset(new CommandLineProgramContext(*argc, *argv));
    // creates a CommandLineProgramContext::CommandLineProgramContext which
    // is a IProgramContext, which seems to be concerned with introspection of
    // process name, path, and working directory.
    // Gets a pointer to new DefaultExecutableEnvironment(), which gets the
    // current working directory and is a IExecutableEnvironment.
    // Hold off on implementatin until the requirements dictate.

    // call setProgramContext(const IProgramContext *programContext)

    // g_libFileFinder.reset(new DataFileFinder());
    // g_libFileFinder->setSearchPathFromEnv("GMXLIB");
    // setLibraryFileFinder(g_libFileFinder.get());
    // TODO: Currently not necessary to run, but can be added for runtime
    // configurability.

    // Retain a handle to the IProgramContext

    // Create the module manager
    // gmx::CommandLineModuleManager manager("gmx", &context);
    // which stashes the name string and a pointer to the program context.
    // Module groups are maintained to organize help output. Various generations
    // of module implementations can be added with different methods and then
    // run(). e.g. manager->addModule(std::move(module));
    // The Module manager maintains a map of names to modules and their "help"
    // machinery.
    // mdrun is registered with
    // manager->addModuleCMainWithSettings(name, shortDescription, mainFunction,  &initSettingsNoNice);
    // where mainFunction is gmx_mdrun(int, char**) from mdrun.cpp
    // from which a wrapper is created with
    // CMainCommandLineModule(const char *name, const char *shortDescription,
    //                               CMainFunction mainFunction,
    //                               InitSettingsFunction settingsFunction)

    // Initialize options holder, handling various universal CLI flags.
    // Create a CommandLineParser with the optionsHolder and parse a subset of argv
    // if the gmx wrapper had CLI options.

    // Get a handle to the CMainCommandLineModule

    // Create a parser, which owns an OptionsAssigner that can act on the Options,
    // to which it has a handle.
    // The parser's .parse() method invokes the OptionsAssigner protocol for
    // impl_->assigner_, with provisions for various option-management protocols.
    // Common options are handled here.
    // optionsHolder.startupInfoFile() decides what the status output filehandle is.

    // Create a CommandLineModuleSettings and allow it to be initialized by the
    // module. For mdrun this just does gmx::CommandLineModuleSettings::setDefaultNiceLevel(0)

    // Then the optionsHolder gets a call optionsHolder.adjustFromSettings(settings);
    // and gmx::CommandLineCommonOptionsHolder updates its niceLevel_

    // gmx_set_max_backup_count(optionsHolder.shouldBackup() ? -1 : 0);
    // allows backup count to be set by optionsHolder, $GMX_MAXBACKUP, or default 99.
//    gmx_set_max_backup_count(-1);

    // optionally open debug file if (optionsHolder.debugLevel() > 0)

    // if (optionsHolder.enableFPExceptions()) { gmx_feenableexcept(); }
//    gmx_feenableexcept();

    // CMainCommandLineModule::run calls gmx_mdrun(argc, argv) with remaining args.
    // Possible mdrun filesystem interactions are defined by the t_filenm fnm[] member.
    // CLI options behavior is defined in the t_pargs pa[] member, which includes
    // many "hidden" options.
    // init_commrec() is called.
    // parse_common_args(...) creates a gmx::Options, which is given a gmx::FileNameOptionManager
    // A gmx::TimeUnitBehavior does timeUnitBehavior->addTimeUnitOption(&options, "tu");
    // and then the behavior is added to a gmx::OptionsBehaviorCollection initialized
    // with &options.
    // A gmx::OptionsAdapter receives the CLI arguments and processes fnm[] and pa[]
    // into &options.
    // All of the options map to gmx::AbstractOptionStorage pointers.
    // Finally, a gmx::CommandLineParser is created for &options.
    // '-s' option is recognized and filename gets passed to
    // gmx::CommandLineParser().impl_->assigner_.appendValue(arg)

    //Note that the
    // OptionsBehavior has its own protocol, and after parsing the following are called:
    // behaviors.optionsFinishing();
    // options.finish();
    // Note that OptionsAssigner is a state machine that relies on a protocol to
    // move it forward.

    // Some MPI awareness and threading stuff is set up and sanity checked.

    // Handle restarts and appending.

    // open log file, set a few more things, then call gmx::mdrunner()
    // TODO: follow up with Mark and Shirts group re reimplementation of mdrunner().
    // instance of t_inputrec is constructed.
    // t_state statInstances is created at runner.cpp:775
    // TPR is read by master with read_tpx_state(ftp2fn(efTPR, nfile, fnm), inputrec, state, mtop);

    // runner.cpp:900 mdModules.assignOptionsToModules(*inputrec->params, nullptr);
    // will allow the electric field module to receive parameters.

    // Initialize state with set_state_entries(state, inputrec);
    // Initialize more stuff... inputrec starts to be modified, e.g.
    // inputrec->pull_work = init_pull(...)

    // my_integrator(inputrec->eI) (...) causes the appropriate do_md(...) to be
    // called.
    // then finalize rot and pull and clean up stuff.

    return system;
}



} // end namespace gmxapi
