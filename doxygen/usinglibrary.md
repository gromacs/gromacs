Using \Gromacs as a library {#page_usinglibrary}
===========================

Getting started
===============

The \Gromacs library (`libgromacs`) provides a few different alternatives for
using it.  These are listed here from the highest level of abstraction to the
low-level functions.
 - If you are writing a trajectory analysis tool, please see
   \ref page_analysisframework.  \ref page_analysistemplate should contain
   all the ingredients to get started.
   If you have an existing tool written using the analysis template from 4.5 or
   4.6 (using the selection engine added in 4.5), you need to do some
   conversion work to get this work with the new template.  This is mostly
   straightforward, but requires some studying to understand the new framework.
 - If you are writing a command line tool for some other purpose, you can use
   the facilities provided by \ref module_commandline.  There are a few
   different alternatives, depending on how much control you want to give
   \Gromacs:
    - For C++ code, you can implement gmx::CommandLineModuleInterface, and
      use gmx::runCommandLineModule() to execute it.  This requires using some
      additional \Gromacs classes (in particular, for implementing
      gmx::CommandLineModuleInterface::writeHelp(), if you want to support the
      `-h` option).
    - For C code, you can use gmx_run_cmain() to wrap an existing C main
      method.  The only constraint on the provided main method is that it
      should use parse_common_args() for argument processing.
      This approach should allow you to convert existing C tools written
      against pre-5.0 \Gromacs (e.g., using the analysis template from 4.0 or
      earlier) to the new version.
    - If you want more control (for example, you do not want the default
      command line options added by \Gromacs), you can directly initialize
      \Gromacs using gmx::initForCommandLine() before calling other \Gromacs
      routines.  This allows you to write your own handling for command line
      options from scratch.  This is also discussed in \ref module_commandline.
 - For most control, you can use gmx::init() to do basic initialization, create
   your own implementation for gmx::ProgramContextInterface, and set that using
   gmx::setProgramContext().  This allows you to customize how the \Gromacs
   library shows the name of the program in messages, as well as how it locates
   its own data files.

If these do not fit your needs, you may need to modify the \Gromacs source code
yourself.  In particular, it is currently relatively difficult to extend the
functionality of `mdrun` without modifying the source code directly.
If you think that some particular API would be necessary for your work, and
think that it would be easy to expose, please drop a line on the
`gmx-developers` mailing list, or contribute the necessary changes on
http://gerrit.gromacs.org/.

Linking against `libgromacs`
============================

\Gromacs is a bit picky on how the headers need to be used: depending on
compilation options used for \Gromacs, some preprocessor defines may need to be
set, the required include path may also depend on compilation options, and some
extra libraries may need to be linked.  You will also likely need to use the
same compiler (or sufficiently similar one that uses the same standard library)
that was used to compile \Gromacs.

To manage this more easily, \Gromacs provides two mechanisms for getting the
correct flags for compilation and linking against the \Gromacs library:
 - `pkg-config`: \Gromacs installs `libgromacs.pc` file (suffixed with the
   library suffix) for use with `pkg-config` if that is present on the system.
   Sourcing `GMXRC` adjusts the `pkg-config` search path such that these files
   are found automatically.
   See `Makefile.pkg` installed with the analysis template for one example of
   how to use it (to use it with a differently suffixed \Gromacs, just replace
   `libgromacs` with `libgromacs`<em>_suffix</em> in the `pkg-config` calls).
 - CMake package configuration files that allow `find_package(GROMACS)` to work.
   See below for details about how to use this in CMake.
   Sourcing `GMXRC` sets an environment variable that allows CMake to find the
   configuration file automatically.
   See `CMakeLists.txt` installed with the analysis template for one example of
   how to use it.

These mechanisms are currently provided on a best-effort basis, but are not
routinely tested on a wide range of configurations.  Please report any issues
with details of how \Gromacs was built so that the mechanism can be improved.

CMake `find_package(GROMACS)` details
-------------------------------------

TODO: Describe FindGROMACS.cmake

Input options:

<dl>
<dt>`GROMACS_SUFFIX`</dt>
<dd>This CMake variable can be set before calling `find_package(GROMACS)` to
specify the \Gromacs suffix to search for.  If not set, an unsuffixed version
is searched for.</dd>
</dl>

Output variables:

<dl>
<dt>`GROMACS_INCLUDE_DIRS`</dt>
<dd>List of include directories necessary to compile against the \Gromacs
headers.</dd>
<dt>`GROMACS_LIBRARIES`</dt>
<dd>List of libraries to link with to link against \Gromacs.
Under the hood, this uses imported CMake targets to represent `libgromacs`.</dd>
<dt>`GROMACS_DEFINITIONS`</dt>
<dd>List of compile definitions (with `-D` in front) that are required to
compile the \Gromacs headers.</dd>
<dt>`GROMACS_IS_DOUBLE`</dt>
<dd>Whether the found \Gromacs was compiled in double precision.</dd>
</dl>

Declared macros/functions:

<dl>
<dt>`gromacs_check_double(GMX_DOUBLE)`</dt>
<dd>Checks that the found \Gromacs is in the expected precision.
The parameter `GMX_DOUBLE` should be the name of a cache variable that
specified whether double-precision was requested.</dd>
<dt>`gromacs_check_compiler(LANG)`<dt>
<dd>Checks that the found \Gromacs was compiled with the same compiler
that is used by the current CMake system.
Currently only `LANG=CXX` is supported.</dd>
</dl>

Notes on \Gromacs API
=====================

The headers for the public \Gromacs API are installed in `include/gromacs/`
under the installation directory.  The layout reflects the source code layout
under the `src/gromacs/` directory (see \ref page_codelayout).  The headers
directly under `include/gromacs/` do not contain any declarations, but instead
include a collection of headers from subdirectories.
You should prefer to include these convenience headers instead of individual
headers from the subdirectories, since they are much more stable.  The
individual headers in the subdirectories can be renamed or moved, but the goal
is to only rarely change the name of these top-level headers.

Pre-5.0 versions of \Gromacs installed (nearly) all headers directly under
`include/gromacs/`.  Most of these headers still exist, but now under
`include/gromacs/legacyheaders/`.  The long-term goal is to move these to
proper module hierarchy or get rid of them, but unfortunately this can take a
long time.  Thus, you should not expect much stability from the API in these
headers.  Some have already been moved, so if you do not find your favorite
header there, try searching for a declaration from the other subdirectories.

For headers under other subdirectories, some effort has been put to design the
API for stability.  However, with limited development resources, and the focus
of \Gromacs being in high performance simulations, all the APIs are subject to
change without notice.  With each new release (with possible exception of patch
releases), you should expect incompatible API changes.  This is in particular
true until the planned reorganization of the `legacyheaders/` subdirectory is
complete.

The header version.h (installed as `gromacs/version.h`) provides defines that
calling code can use to check the exact (released) version of \Gromacs that
installed the headers.

This Doxygen documentation only covers part of the API.  In particular, nearly
all of `include/gromacs/legacyheaders/` is undocumented, as well as code
recently moved from there.
