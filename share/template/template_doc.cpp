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
 * The analysis tool class inherits gmx::TrajectoryAnalysisModule, which is
 * an interface with a few convenience functions for easier interfacing
 * with other code.
 * Below, we walk through the different methods as implemented in the
 * template (note that the template does not implement some of the virtual
 * methods because they are less often needed), discussing some issues that can
 * arise in more complex cases.
 * See documentation of gmx::TrajectoryAnalysisModule for a full description of
 * the available virtual methods and convenience functions.
 *
 * We also declare a helper class, AnalysisTemplate::ModuleData, that derives
 * from gmx::TrajectoryAnalysisModuleData and will contain any data that needs
 * to be frame-local in parallel analysis (if you don't care about
 * parallelization, you can just include these variables in the module class
 * itself, initialize them in gmx::TrajectoryAnalysisModule::initAnalysis(),
 * and do any postprocessing in
 * gmx::TrajectoryAnalysisModule::finishAnalysis()).
 * \until };
 * See documentation of gmx::TrajectoryAnalysisModuleData for more details of
 * how this data can be used.
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
 * to the \p _options object.  These values will only affect the help output.
 *
 *
 * \section template_initialization Initialization methods
 *
 * Initialization of the module is split into a few methods, two of which are
 * used in the template.  gmx::TrajectoryAnalysisModule::initOptions() is used
 * to set up options understood by the module, as well as for setting up
 * different options through gmx::TrajectoryAnalysisSettings (see the
 * documentation of that class for more details):
 * \skip  Options *
 * \until return &_options;
 * \until }
 * For additional documentation on how to define different kinds of options, see
 * gmx::Options, basicoptions.h, and gmx::SelectionOption.
 *
 * If you need to set up settings based on option values, you can also override
 * gmx::TrajectoryAnalysisModule::initOptionsDone().  For simplicity,
 * this is not done in the template.
 *
 * The actual analysis is initialized in
 * gmx::TrajectoryAnalysisModule::initAnalysis():
 * \skip  int
 * \until return 0;
 * \until }
 * Information about the topology is passed as a parameter.
 * If the analysis module needs some temporary storage during processing of a
 * frame, this should be allocated in gmx::TrajectoryAnalysisModule::startFrames()
 * (see below) if parallelization is to be supported.
 *
 * There is also a gmx::TrajectoryAnalysisModule::initAfterFirstFrame() method
 * that can be overridden if the module needs to do some initialization based
 * on data from the first frame (most commonly, based on the box size), but it
 * is not used in the template.
 *
 *
 * \section template_analysis Actual trajectory analysis
 *
 * There is one more initialization method that needs to be overridden to
 * support automatic parallelization.  If you do not need custom data (or
 * parallelization at all), you can skip this method and ignore the last
 * parameter to gmx::TrajectoryAnalysisModule::analyzeFrame() to make things
 * simpler.  In the template, we do include it:
 * \skip  int
 * \until return rc;
 * \until }
 *
 * The main part of the analysis is (in most analysis codes) done in the
 * gmx::TrajectoryAnalysisModule::analyzeFrame() method, which is called once
 * for each frame:
 * \skip  int
 * \until {
 * The \p frnr parameter gives a zero-based index of the current frame
 * (mostly for use with gmx::AnalysisData), \p pbc contains the PBC
 * information for the current frame for distance calculations with,
 * e.g., pbc_dx(), and \p pdata points to a data structure created in
 * gmx::TrajectoryAnalysisModule::startFrames().
 * Although usually not necessary (except for the time field), raw frame
 * data can be accessed through \p fr.
 * In most cases, the analysis should be written such that it gets all
 * position data through selections, and does not assume a constant size for
 * them.  This is all that is required to support the full flexibility of the
 * selection engine.
 *
 * First, we get data from our custom data structure (note the cast on the
 * second line to access our custom data) for shorthand access:
 * \skip  AnalysisDataHandle
 * \until NeighborhoodSearch
 *
 * For the template, we do a simple calculation:
 * \skip  nb
 * \until finishFrame()
 * Finally, we return zero to indicate that all went well.
 * \skipline return 0;
 *
 * After all the frames have been processed,
 * gmx::TrajectoryAnalysisModule::finishAnalysis() is called once.  This is the
 * place to do any custom postprocessing of the data.  For the template, we do
 * nothing, because all necessary processing is done in the data modules:
 * \skip  int
 * \until }
 *
 * There is one additional method that can be implemented if custom data from
 * frames is aggregated in the frame-local data structures created in
 * gmx::TrajectoryAnalysisModule::startFrames().  In such cases,
 * gmx::TrajectoryAnalysisModule::finishFrames() can be overridden to combine
 * the data from the data structures (see documentation of the method for
 * details).
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
 * \skipline int
 * \until return 0;
 * \until }
 *
 *
 * \section template_main Definition of main()
 *
 * Now, the only thing remaining is to define the main() function.
 * To implement a command-line tool, it should create a module and run it using
 * gmx::TrajectoryAnalysisCommandLineRunner:
 * \skip  int
 * \until }
 */
