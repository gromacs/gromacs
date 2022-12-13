Framework for energy analysis {#page_energyanalysisframework}
=============================

\Gromacs provides a framework for implementing flexible energy data analysis
routines.  It consists of a few components that can also be used individually,
but in most cases it is desirable to use features from all of them to get most
out of the framework.  The main features are:

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

High-level framework
====================

The \ref module_energyanalysis module provides the high-level framework
that integrates all the pieces together.
It provides the abstract base class for analysis tool modules
(gmx::EnergyAnalysisModule), and the code that runs such a module as a
command-line tool (gmx::EnergyAnalysisCommandLineRunner).
See [energyanalysis module documentation](\ref module_energyanalysis) for
more details.


Output data handling
====================

The \ref module_analysisdata module provides two things:

 - Support for uniformly providing output data from analysis tools.
   Tools compute their output values and place them into a
   _data object_ for further processing.  This allows two things:
     - Reusable data modules can be applied across different tools to do common
       post-processing.
     - The data object provides parallelization support.
 - Set of reusable data modules for post-processing the data.  These include
   functionality like averaging data, computing histograms, and plotting the
   data into a file.  Many of these modules also provide their output as a data
   object, allowing further data modules to be attached to them.

The general concept is explained in more detail on a separate page:
\subpage page_analysisdata.
The [analysisdata module documentation](\ref module_analysisdata) provides more
technical details.


Input options {#section_energyanalysisframework_options}
=============

To declare input data for the tool (typically, command-line options, including
input files and selections), \ref module_options module is used.
The analysis tool code receives an instance of gmx::IOptionsContainer for one of
its initialization methods, and uses it to provide its input options.
Basic options are declared in basicoptions.h. 
For each option, the tool declares a local variable
that will receive the value for that option.  After the options are parsed from
the command line (by the framework), the tool code can read the values from
these variables.  The option declarations filled into the
gmx::IOptionsContainer object are also used to provide help to the user (also
handled by the framework).
See the documentation for gmx::EnergyAnalysisModule and the
[options module documentation](\ref module_options) for more details.

