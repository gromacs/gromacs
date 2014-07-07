Example code for writing trajectory analysis tools {#page_analysistemplate}
==================================================

\Gromacs installation includes a template for writing trajectory analysis
tools using \ref page_analysisframework.
It can be found from `share/gromacs/template/` under the installation
directory, and from `share/template/` in the source distribution.

The full source code for the file is also included in this documentation:
\ref template.cpp "template.cpp"
The rest of this page walks through the code to explain the different parts.

\dontinclude template.cpp

Global definitions
==================

We start by including some generic C++ headers:
\skip  <string>
\until <vector>
and continue by including the header for the analysis library:
\skipline <gromacs/trajectoryanalysis.h>
This header includes other headers that together define all the basic data
types needed for writing trajectory analysis tools.
For convenience, we also import all names from the ::gmx namespace into the
global scope to avoid repeating the name everywhere:
\skipline using namespace


Tool module class declaration
=============================

We then define a class that implements our analysis tool:
\skip  AnalysisTemplate
\until };
The analysis tool class inherits from gmx::TrajectoryAnalysisModule, which
is an interface with a few convenience functions for easier interfacing
with other code.
Below, we walk through the different methods as implemented in the
template (note that the template does not implement some of the virtual
methods because they are less often needed), discussing some issues that can
arise in more complex cases.
See documentation of gmx::TrajectoryAnalysisModule for a full description of
the available virtual methods and convenience functions.
The first block of member variables are used to contain values provided to
the different options.  They will vary depending on the needs of the
analysis tool.
The AnalysisNeighborhood object provides neighborhood searching that is used
in the analysis.
The final block of variables are used to process output data.
See initAnalysis() for details on how they are used.

For the template, we do not need any custom frame-local data.  If you think
you need some for more complex analysis needs, see documentation of
gmx::TrajectoryAnalysisModuleData for more details.
If you do not care about parallelization, you do not need to consider this
part.  You can simply declare all variables in the module class itself,
initialize them in gmx::TrajectoryAnalysisModule::initAnalysis(), and do any
postprocessing in gmx::TrajectoryAnalysisModule::finishAnalysis()).


Construction
============

The constructor (and possible destructor) of the analysis module should be
simple: the constructor should just initialize default values, and the
destructor should free any memory managed by the module.  For the template,
we have no attributes in our class that need to be explicitly freed, so we
declare only a constructor:
\skip  AnalysisTemplate
\until }
In addition to initializing local variables that don't have default
constructors, we also provide a title and one-line description of our module
to the \p options_ object.  These values will only affect the help output.


Input options
=============

Initialization of the module is split into a few methods, two of which are
used in the template.  gmx::TrajectoryAnalysisModule::initOptions() is used
to set up options understood by the module, as well as for setting up
different options through gmx::TrajectoryAnalysisSettings (see the
documentation of that class for more details):
\skip  void
\until settings->
\until }
For the template, we first set a description text for the tool (used for
help text).  Then we declare an option to specify the output file name,
followed by options that are used to set selections, and finally an option
to set a cutoff value.  For the cutoff, the default value will be the one
that was set in the constructor, but it would also be possible to explicitly
set it here.  The values provided by the user for the options will be stored
in member variables.  Finally, we indicate that the tool always requires
topology information.  This is done for demonstration purposes only; the
code in the template works even without a topology.

For additional documentation on how to define different kinds of options, see
gmx::Options, basicoptions.h, and gmx::SelectionOption.  You only need to
define options that are specific to the analysis; common options, e.g., for
specifying input topology and trajectories are added by the framework.

To adjust settings or selection options (e.g., the number of accepted
selections) based on option values, you need to override
gmx::TrajectoryAnalysisModule::optionsFinished().  For simplicity,
this is not done in the template.


Analysis initialization
=======================

The actual analysis is initialized in
gmx::TrajectoryAnalysisModule::initAnalysis():
\skip  void
\until }
\until }
Information about the topology is passed as a parameter.  The settings
object can also be used to access information about user input.

One of the main tasks of this method is to set up appropriate
gmx::AnalysisData objects and modules for them (see
gmx::TrajectoryAnalysisModule for the general approach).
These objects will be used to process output from the tool.  Their main
purpose is to support parallelization, but even if you don't care about
parallelism, they still provide convenient building blocks, e.g., for
histogramming and file output.

For the template, we first set the cutoff for the neighborhood search.

Then, we create and register one gmx::AnalysisData object
that will contain, for each frame, one column for each input selection.
This will contain the main output from the tool: minimum distance between
the reference selection and that particular selection.
We then create and setup a module that will compute the average distance
for each selection (see writeOutput() for how it is used).
Finally, if an output file has been provided, we create and setup a module
that will plot the per-frame distances to a file.

If the analysis module needs some temporary storage during processing of a
frame (i.e., it uses a custom class derived from
gmx::TrajectoryAnalysisModuleData), this should be allocated in
gmx::TrajectoryAnalysisModule::startFrames() (see below) if parallelization
is to be supported.

If you need to do initialization based on data from the first frame (most
commonly, based on the box size), you need to override
gmx::TrajectoryAnalysisModule::initAfterFirstFrame(), but this is not used
in the template.


Analyzing the frames
====================

There is one more initialization method that needs to be overridden to
support automatic parallelization: gmx::TrajectoryAnalysisModule::startFrames().
If you do not need custom frame-local data (or parallelization at all), you
can skip this method and ignore the last parameter to
gmx::TrajectoryAnalysisModule::analyzeFrame() to make things simpler.
In the template, this method is not necessary.

The main part of the analysis is (in most analysis codes) done in the
gmx::TrajectoryAnalysisModule::analyzeFrame() method, which is called once
for each frame:
\skip  void
\until {
The \p frnr parameter gives a zero-based index of the current frame
(mostly for use with gmx::AnalysisData), \p pbc contains the PBC
information for the current frame for distance calculations with,
e.g., pbc_dx(), and \p pdata points to a data structure created in
startFrames().
Although usually not necessary (except for the time field), raw frame
data can be accessed through \p fr.
In most cases, the analysis should be written such that it gets all
position data through selections, and does not assume a constant size for
them.  This is all that is required to support the full flexibility of the
selection engine.

For the template, we first get data from our custom data structure for
shorthand access (if you use a custom data object, you need a \c static_cast
here):
\skip  AnalysisDataHandle
\until parallelSelection

We then do a simple calculation and use the AnalysisDataHandle class to set
the per-frame output for the tool:
\skip  nb
\until finishFrame()

After all the frames have been processed,
gmx::TrajectoryAnalysisModule::finishAnalysis() is called once.  This is the
place to do any custom postprocessing of the data.  For the template, we do
nothing, because all necessary processing is done in the data modules:
\skip  void
\until }

If the data structure created in gmx::TrajectoryAnalysisModule::startFrames()
is used to aggregate data across frames, you need to override
gmx::TrajectoryAnalysisModule::finishFrames() to combine the data from the
data structures (see documentation of the method for details).
This is not necessary for the template, because the ModuleData structure
only contains data used during the analysis of a single frame.


Output
======

Finally, most programs need to write out some values after the analysis is
complete.  In some cases, this can be achieved with proper chaining of data
modules, but often it is necessary to do some custom processing.
All such activities should be done in
gmx::TrajectoryAnalysisModule::writeOutput().  This makes it easier to reuse
analysis modules in, e.g., scripting languages, where output into files
might not be desired.  The template simply prints out the average distances
for each analysis group:
\skip  void
\until }
\until }
Here, we use the \c avem_ module, which we initialized in initAnalysis() to
aggregate the averages of the computed distances.


Definition of main()
====================

Now, the only thing remaining is to define the main() function.
To implement a command-line tool, it should create a module and run it using
gmx::TrajectoryAnalysisCommandLineRunner using the boilerplate code below:
\skip  int
\until }


\if libapi
Tools within \Gromacs
====================

Analysis tools implemented using the template can also be easily included into
the \Gromacs library.  To do this, follow these steps:

 1. Put your tool source code into `src/gromacs/trajectoryanalysis/modules/`.
 2. Remove `using namespace gmx;` and enclose all the code into
    `gmx::analysismodules` namespace, and the tool class into an unnamed
    namespace within this.
 3. Create a header file corresponding to your tool and add the following class
    into it withing `gmx::analysismodules` namespace (replace `Template` with
    the name of your tool):
~~~~{.cpp}
    class TemplateInfo
    {
        public:
            static const char name[];
            static const char shortDescription[];
            static TrajectoryAnalysisModulePointer create();
    };
~~~~
 4. Add definition for these items in the source file, outside the unnamed
    namespace (replace `Template`, `AnalysisTemplate` and the strings with
    correct values):
~~~~{.cpp}
    const char TemplateInfo::name[]             = "template";
    const char TemplateInfo::shortDescription[] =
        "Compute something";

    TrajectoryAnalysisModulePointer TemplateInfo::create()
    {
        return TrajectoryAnalysisModulePointer(new AnalysisTemplate);
    }
~~~~
 5. Change the constructor of your tool to refer to the strings in the new
    class:
~~~~{.cpp}
    AnalysisTemplate::AnalysisTemplate()
        : TrajectoryAnalysisModule(TemplateInfo::name, TemplateInfo::shortDescription),
~~~~
 6. Register your module in `src/gromacs/trajectoryanalysis/modules.cpp`.
 7. Done.  Your tool can now be invoked as `gmx template`, using the string you
    specified as the name.

See existing tools within the `src/gromacs/trajectoryanalysis/modules/` for
concrete examples and preferred layout of the files.  Please document yourself
as the author of the files, using Doxygen comments like in the existing files.

\endif
