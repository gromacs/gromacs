Framework for trajectory analysis {#page_analysisframework}
=================================

\Gromacs provides a framework for implementing flexible trajectory analysis
routines.  It consists of a few components that can also be used individually,
but in most cases it is desirable to use features from all of them to get most
out of the framework.  The main features are:

 - Support for flexible selections that can be used to provide the set of
   coordinates to analyze.  They can be dynamic, i.e., select different atoms
   for different trajectory frames, and also support evaluation of
   center-of-mass/center-of-geometry coordinates for a group of atoms.
   The latter means that a tool written to use the framework can automatically
   analyze also center-of-mass positions (or a mixture of center-of-mass and
   atomic positions) in addition to real atomic positions.
 - Support for per-frame parallelization.  The framework is designed to
   support running the analysis in parallel for multiple frames for cases where
   different frames can be analyzed (mostly) independently.  At this time, the
   actual parallelization is not implemented, but tools written to use the
   framework should be able to take advantage of it as soon as it materializes
   with no or minimal changes.
 - Access to a library of basic analysis routines.  Things such as computing
   averages and histogramming are provided as reusable modules.
 - Tool code can focus on the actual analysis.  Tools are implemented by
   subclassing an abstract class and providing an implementation for selected
   pure virtual methods.  The framework takes care of initialization tasks,
   loading the trajectory and looping over it, evaluating selections, and also
   provides basic features like making molecules whole before passing the frame
   to the analysis code.
   This approach also makes it possible to reuse the same tool code from a
   scripting language such as Python simply by implementing general support for
   such language bindings in the framework (no such integration is implemented
   at this time, though).

For a crash course on how to implement an analysis tool using the framework, see
\subpage page_analysistemplate.


High-level framework
====================

The \ref module_trajectoryanalysis module provides the high-level framework
that integrates all the pieces together.
It provides the abstract base class for analysis tool modules
(gmx::TrajectoryAnalysisModule), and the code that runs such a module as a
command-line tool (gmx::TrajectoryAnalysisCommandLineRunner).
See the [analysis template](\ref page_analysistemplate) and the
[trajectoryanalysis module documentation](\ref module_trajectoryanalysis) for
more details.


Selections
==========

The \ref module_selection module provides the support for selections.
Most of the work of managing the selections is taken care by the command-line
runner and the framework, and the analysis tool code only sees two main
classes:

 - gmx::SelectionOption and associated classes are used to declare the
   number and type of selections the tool accepts (see below for
   [details of the option support](#section_analysisframework_options)).
 - The tool receives a set of gmx::Selection objects as a value of the
   selection option.  These classes provide the evaluated value of the
   selections for during the analysis.  The framework evaluates them for each
   frame such that when the tool is called, it can access the selections for
   the current frame in the gmx::Selection objects it owns.

More details of the selection engine is available in the
[selection module documentation](\ref module_selection).
This is useful in particular for understanding how the selections work in
detail, or if you want to use the selection code outside the trajectory
analysis framework.


Output data handling
====================

The \ref module_analysisdata module provides two things:

 - Support for uniformly providing output data from the tool.
   Typically, the tool provides its raw data output through one or more
   gmx::AnalysisData objects.
   These data objects subclass gmx::AbstractData, and represent simple data
   that has a set of frames (typically the same number as in the input
   trajectory), and can contain one or more data sets, each of which has a set
   of columns.  The tool computes the values and places them into this data
   object for further processing.  This allows two things:
     - Reusable data modules can be applied across different tools to do common
       post-processing.  See below.
     - The gmx::AnalysisData object takes care of parallelization support.
       When multiple frames are being processed in parallel, the data may
       become available out of order.  In particular for writing the per-frame
       data into a file, but also for other types of post-processing, it is
       necessary to reorder the data sequentially.  This is implemented once by
       the framework, and the tool code does not need to worry, other than
       using the provided API.
 - Set of reusable data modules, each of which implements
   gmx::AnalysisDataModuleInterface and can post-process data in a
   gmx::AbstractData object.  These include functionality like averaging data,
   computing histograms, and plotting the data into a file.  Many of these
   modules also subclass gmx::AbstractData to provide the processed data like
   the averages or the histograms, allowing further data modules to be attached
   to them.

See the [analysisdata module documentation](\ref module_analysisdata) for more
details.


Input options {#section_analysisframework_options}
=============

To declare input data for the tool (typically, command-line options, including
input files and selections), \ref module_options module is used.
The analysis tool code receives a pre-initialized gmx::Options object in one of
its initialization methods, and fills it with its input options.
Basic options are declared in basicoptions.h, and also gmx::SelectionOption is
used in the same manner.  For each option, the tool declares a local variable
that will receive the value for that option.  After the options are parsed from
the command line (by the framework), the tool code can read the values from
these variables.  The option declarations, and other information filled into
the gmx::Options object, are also used to provide help to the user (also
handled by the framework).
See the documentation for gmx::TrajectoryAnalysisModule and the
[options module documentation](\ref module_options) for more details.
