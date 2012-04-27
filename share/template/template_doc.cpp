/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \dir share/template
 * \brief Template code for writing analysis programs.
 */
/*! \example template.cpp
 * \brief Template code for writing analysis programs.
 *
 * See \ref share/template/template.cpp "documentation" of the template for
 * more information.
 */
/*! \internal \file
 * \brief Doxygen documentation source for template.cpp.
 */
/*! \file template.cpp
 * \brief Template code for writing analysis programs.
 *
 * The full source code for the file: \ref template.cpp "template.cpp"
 *
 * \dontinclude template.cpp
 *
 * \section template_global Global definitions
 *
 * We start by including some generic C++ headers:
 * \skip  <string>
 * \until <vector>
 * and continue by including the header for the analysis library:
 * \until <gromacs/trajectoryanalysis.h>
 * This header includes other headers that together define all the basic data
 * types needed for writing trajectory analysis tools.
 * For convenience, we also import all names from the ::gmx namespace into the
 * global scope to avoid repeating the name everywhere:
 * \skipline using namespace
 *
 *
 * \section template_class Tool module class declaration
 *
 * We then define a class that implements our analysis tool:
 * \until };
 * The analysis tool class inherits from gmx::TrajectoryAnalysisModule, which
 * is an interface with a few convenience functions for easier interfacing
 * with other code.
 * Below, we walk through the different methods as implemented in the
 * template (note that the template does not implement some of the virtual
 * methods because they are less often needed), discussing some issues that can
 * arise in more complex cases.
 * See documentation of gmx::TrajectoryAnalysisModule for a full description of
 * the available virtual methods and convenience functions.
 * The member variable \c options_ will be used to specify command-line options
 * that the tool accepts, and such a variable needs to be present to properly
 * implement gmx::TrajectoryAnalysisModule::initOptions().
 * The next block of member variables are used to contain values provided to
 * the different options.  They will vary depending on the needs of the
 * analysis tool.
 * The final block of variables are used to process output data.
 * See initAnalysis() for details on how they are used.
 *
 * We also declare a helper class, AnalysisTemplate::ModuleData, that derives
 * from gmx::TrajectoryAnalysisModuleData and will contain any data that needs
 * to be frame-local in parallel analysis:
 * \until };
 * See documentation of gmx::TrajectoryAnalysisModuleData for more details of
 * how this data can be used.  Note that not all types of analysis will require
 * a custom type for this purpose.
 * Also, if you don't care about parallelization, you can just include these
 * variables in the module class itself, initialize them in
 * gmx::TrajectoryAnalysisModule::initAnalysis(), and do any postprocessing in
 * gmx::TrajectoryAnalysisModule::finishAnalysis()).
 *
 *
 * \section template_ctor Construction
 *
 * The constructor (and possible destructor) of the analysis module should be
 * simple: the constructor should just initialize default values, and the
 * destructor should free any memory managed by the module.  For the template,
 * we have no attributes in our class that need to be explicitly freed, so we
 * declare only a constructor:
 * \skip  AnalysisTemplate
 * \until }
 * In addition to initializing local variables that don't have default
 * constructors, we also provide a title and one-line description of our module
 * to the \p options_ object.  These values will only affect the help output.
 *
 *
 * \section template_options Input options
 *
 * Initialization of the module is split into a few methods, two of which are
 * used in the template.  gmx::TrajectoryAnalysisModule::initOptions() is used
 * to set up options understood by the module, as well as for setting up
 * different options through gmx::TrajectoryAnalysisSettings (see the
 * documentation of that class for more details):
 * \skip  Options &
 * \until return options_;
 * \until }
 * For the template, we first set a description text for the tool (used for
 * help text).  Then we declare an option to specify the output file name,
 * followed by options that are used to set selections, and finally an option
 * to set a cutoff value.  For the cutoff, the default value will be the one
 * that was set in the constructor, but it would also be possible to explicitly
 * set it here.  The values provided by the user for the options will be stored
 * in member variables.  Finally, we indicate that the tool always requires
 * topology information.  This is done for demonstration purposes only; the
 * code in the template works even without a topology.
 *
 * For additional documentation on how to define different kinds of options, see
 * gmx::Options, basicoptions.h, and gmx::SelectionOption.  You only need to
 * define options that are specific to the analysis; common options, e.g., for
 * specifying input topology and trajectories are added by the framework.
 *
 * To adjust settings or selection options (e.g., the number of accepted
 * selections) based on option values, you need to override
 * gmx::TrajectoryAnalysisModule::initOptionsDone().  For simplicity,
 * this is not done in the template.
 *
 *
 * \section template_initialization Analysis initialization
 *
 * The actual analysis is initialized in
 * gmx::TrajectoryAnalysisModule::initAnalysis():
 * \skip  void
 * \until }
 * \until }
 * Information about the topology is passed as a parameter.  The settings
 * object can also be used to access information about user input.
 *
 * One of the main tasks of this method is to set up appropriate
 * gmx::AnalysisData objects and modules for them (see
 * gmx::TrajectoryAnalysisModule for the general approach).
 * These objects will be used to process output from the tool.  Their main
 * purpose is to support parallelization, but even if you don't care about
 * parallelism, they still provide convenient building blocks, e.g., for
 * histogramming and file output.
 *
 * For the template, we first create and register one gmx::AnalysisData object
 * that will contain, for each frame, one column for each input selection.
 * This will contain the main output from the tool: minimum distance between
 * the reference selection and that particular selection.
 * We then create and setup a module that will compute the average distance
 * for each selection (see writeOutput() for how it is used).
 * Finally, if an output file has been provided, we create and setup a module
 * that will plot the per-frame distances to a file.
 *
 * If the analysis module needs some temporary storage during processing of a
 * frame (i.e., it uses a custom class derived from
 * gmx::TrajectoryAnalysisModuleData), this should be allocated in
 * gmx::TrajectoryAnalysisModule::startFrames() (see below) if parallelization
 * is to be supported.
 *
 * If you need to do initialization based on data from the first frame (most
 * commonly, based on the box size), you need to override
 * gmx::TrajectoryAnalysisModule::initAfterFirstFrame(), but this is not used
 * in the template.
 *
 *
 * \section template_analysis Actual trajectory analysis
 *
 * There is one more initialization method that needs to be overridden to
 * support automatic parallelization.  If you do not need custom data (or
 * parallelization at all), you can skip this method and ignore the last
 * parameter to gmx::TrajectoryAnalysisModule::analyzeFrame() to make things
 * simpler.  In the template, we do include it to initialize our custom
 * ModuleData object:
 * \skip  TrajectoryAnalysisModuleDataPointer
 * \until }
 *
 * The main part of the analysis is (in most analysis codes) done in the
 * gmx::TrajectoryAnalysisModule::analyzeFrame() method, which is called once
 * for each frame:
 * \skip  void
 * \until {
 * The \p frnr parameter gives a zero-based index of the current frame
 * (mostly for use with gmx::AnalysisData), \p pbc contains the PBC
 * information for the current frame for distance calculations with,
 * e.g., pbc_dx(), and \p pdata points to a data structure created in
 * startFrames().
 * Although usually not necessary (except for the time field), raw frame
 * data can be accessed through \p fr.
 * In most cases, the analysis should be written such that it gets all
 * position data through selections, and does not assume a constant size for
 * them.  This is all that is required to support the full flexibility of the
 * selection engine.
 *
 * For the template, we first get data from our custom data structure for
 * shorthand access (note the cast on the second line to access our custom
 * data):
 * \skip  AnalysisDataHandle
 * \until parallelSelection
 *
 * We then do a simple calculation and use the AnalysisDataHandle class to set
 * the per-frame output for the tool:
 * \skip  nb
 * \until finishFrame()
 *
 * After all the frames have been processed,
 * gmx::TrajectoryAnalysisModule::finishAnalysis() is called once.  This is the
 * place to do any custom postprocessing of the data.  For the template, we do
 * nothing, because all necessary processing is done in the data modules:
 * \skip  void
 * \until }
 *
 * If the data structure created in gmx::TrajectoryAnalysisModule::startFrames()
 * is used to aggregate data across frames, you need to override
 * gmx::TrajectoryAnalysisModule::finishFrames() to combine the data from the
 * data structures (see documentation of the method for details).
 * This is not necessary for the template, because the ModuleData structure
 * only contains data used during the analysis of a single frame.
 *
 *
 * \section template_output Output
 *
 * Finally, most programs need to write out some values after the analysis is
 * complete.  In some cases, this can be achieved with proper chaining of data
 * modules, but often it is necessary to do some custom processing.
 * All such activities should be done in
 * gmx::TrajectoryAnalysisModule::writeOutput().  This makes it easier to reuse
 * analysis modules in, e.g., scripting languages, where output into files
 * might not be desired.  The template simply prints out the average distances
 * for each analysis group:
 * \skip  void
 * \until }
 * \until }
 * Here, we use the \c avem_ module, which we initialized in initAnalysis() to
 * aggregate the averages of the computed distances.
 *
 *
 * \section template_main Definition of main()
 *
 * Now, the only thing remaining is to define the main() function.
 * To implement a command-line tool, it should create a module and run it using
 * gmx::TrajectoryAnalysisCommandLineRunner using the boilerplate code below:
 * \skip  int
 * \until return 1;
 * \until }
 * \until }
 *
 *
 * \section template_references Where to go from here?
 *
 * For more information about the topics discussed here, see the following
 * pages:
 *  - \ref module_trajectoryanalysis
 *  - \ref module_selection
 *  - \ref module_analysisdata
 */
