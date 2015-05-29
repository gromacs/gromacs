Wrapper binary implementation {#page_wrapperbinary}
=============================

This page mainly describes the implementation of the `gmx` wrapper binary.
Many of the details are not visible to most of the code, but this documentation
is included as part of the library API documentation to make it easier to
understand the overall implementation without reading extensive documentation.

%main() implementation
======================

The %main() method for the wrapper binary is implemented in
`src/programs/gmx.cpp`.  This is a very simple code that does these basic
tasks:
 1. Initializes \Gromacs using gmx::initForCommandLine()
    (see \ref page_usinglibrary).
 2. Creates a gmx::CommandLineModuleManager instance for the wrapper binary.
 3. Calls various methods to add modules to the manager and initialize it
    otherwise.  Many of the pre-5.0 binaries are added from
    `src/programs/legacymodules.cpp`.  New C++ tools are added from
    `src/gromacs/trajectoryanalysis/modules.cpp`.
 4. Passes control to the manager (see below).
 5. On successful return, deinitializes \Gromacs and returns the exit code from
    the manager.
The %main() method also catches all exceptions, and if one is caught, prints an
error message and terminates the program cleanly.

Command line modules
====================

All modules within the wrapper binary are implemented as classes that implement
the gmx::CommandLineModuleInterface interface.  There is generally some helper
class in between:
 * General C++ modules typically use gmx::Options for their command-line
   handling.  Instead of each module implementing parsing and help separately
   with identical code, they implement gmx::CommandLineOptionsModuleInterface
   instead.  The framework then provides a bridge class that contains the
   common code and wraps gmx::CommandLineOptionsModuleInterface into a
   gmx::CommandLineModuleInterface.
 * For C++ trajectory analysis modules, there is a general implementation for
   running the gmx::TrajectoryAnalysisModule subclasses in cmdlinerunner.cpp.
 * For old C-style %main() functions, see \ref section_wrapperbinary_cmain.

Command line manager {#section_wrapperbinary_manager}
====================

The core of the wrapper binary is the gmx::CommandLineModuleManager::run()
method.  This method:
 1. Parses the command line arguments before the module name as arguments to
    the wrapper binary.  Some arguments such as `-h` and `-version` cause rest
    of the command (the module name and all that follows) to be ignored.
 2. If a module is specified, also checks the command line arguments after the
    module name for the options understood by the wrapper binary, such as `-h`
    and `-version` (see below for details of how `-h` works).  Any such options
    are handled by the manager and removed from the command line for further
    processing.
 3. Print the startup header (contents of which can be controlled by the
    command line options).
 4. If a command line option requests termination after the startup header
    (such as `-version`), return.
 5. Passes control to the selected module.  If there is no module specified,
    the help module is invoked (see below).
 6. Print a quote at the end, and return the exit code from the module.

Command line help
-----------------

To handle the `gmx help ...` command, as well as for `gmx -h` and for
`gmx` _module_ `-h`, the command line manager internally creates a module that
handles the `help` command.  All command lines containing the `-h`, as well as
invocation of `gmx` without any arguments, are translated to corresponding
`gmx help` commands.  For example, `gmx` _module_ `-h` is handled exactly like
`gmx help` _module_.  Note that if `-h` is specified for a module, the command
line manager throws away all the other arguments before passing control to the
module.

After the above translations, the internal help module handles all the help
output.  All the help is organized into a hierarchy of gmx::HelpTopicInterface
instances.  The help module internally creates a root help topic that is
printed with `gmx help`.  If there are additional words after the `gmx help`
command, then those are taken to specify the topic to show in the hierarchy.

gmx::CommandLineModuleManager internally creates a help topic for each added
module.  These topics are shown when `gmx help` _module_ is invoked.
They forward the request to the actual module (to
gmx::CommandLineModuleInterface::writeHelp()).

In addition to the topics created internally, gmx::CommandLineModuleManager
provides methods to add additional help topics.  Currently, this is used to
expose some reference material for the selections (the same content that is
accessible using `help` in the selection prompt).

Help in other formats
---------------------

The build system provides a target, `make sphinx-programs`, that generates
reStructuredText help for the commands, which in turn is used to generate man
and HTML help.  Internally, this executes `gmx help -export rst`, which
triggers special handling in the internal help module.
See documentation for
\linktodevmanual{build-system,special targets in the build system} for details
of which targets to use for generating the documentation..

If this option is set, the help module loops through all the modules in the
binary, writing help for each into a separate file.  The help module writes
common headers and footers, and asks the actual module to write the
module-specific content (with gmx::CommandLineModuleInterface::writeHelp(),
using a different help context than for console output).

Additionally, a list of all the modules is generated (`gromacs.7` for man
pages, and alphabetical and by-topic lists for the HTML pages).

Handling C %main() functions {#section_wrapperbinary_cmain}
----------------------------

Many pre-5.0 modules are still implemented as a function with a C %main()
signature.  All these binaries call parse_common_args() as more or less the
first thing in their processing.  In order to implement the above approach, the
module manager internally creates a command line module for these (in
gmx::CommandLineModuleManager::addModuleCMain()).  The created module
collaborates with parse_common_args() to achieve the same functionality as for
the new C++ modules.

Running the module simply executes the provided %main() method.
Help writing is more complex, as it requires the help context to be passed from
the module to parse_common_args().  This is handled using a global instance of
the context (see gmx::GlobalCommandLineHelpContext).  This context is set in
the module, and if parse_common_args() detects it, it prints out the help and
returns `false` to indicate to the caller that it should immediately return.
