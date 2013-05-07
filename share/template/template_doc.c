/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \dir share/template
 * \brief Template code for writing analysis programs.
 */
/*! \example template.c
 * \brief Template code for writing analysis programs.
 *
 * See \ref share/template/template.c "documentation" of the template for
 * more information.
 */
/*! \file
 * \brief Doxygen documentation source for template.c.
 */
/*! \file template.c
 * \brief Template code for writing analysis programs.
 *
 * The full source code for the file: \ref template.c "template.c"
 *
 * \dontinclude template.c
 *
 * \section global Global definitions
 *
 * We start by including some generic Gromacs headers
 * \until <xvgr.h>
 * and continue by including headers for the analysis library:
 * \until <trajana.h>
 * Here, we include just nbsearch.h that contains routines for performing
 * neighborhood searches for a set of positions.
 * Possibile additional header files include displacement.h and histogram.h;
 * see \ref include "the documentation of user header files".
 *
 * We then define a structure that holds data required during the analysis
 * in analyze_frame():
 * \until t_analysisdata;
 *
 *
 * \section analysis Frame analysis function
 *
 * We're now ready to define the function that (in most analysis codes) does
 * the main part of the work:
 * \until {
 * The topology can be accessed through \p top and the coordinates, box
 * size etc. through \p fr.
 * \p pbc can be used for periodic boundary calculations with, e.g., pbc_dx().
 * \p sel contains the selections defined by the user, containing both the
 * positions and the index groups that have been used to evaluate those.
 * The number of selections is available in \p nr.
 * The analysis should be written such that does not assume the selections to
 * have constant size.
 * The last parameter \p data points to our custom data (t_analysisdata), so we
 * first define a pointer that can be used to access the data:
 * \line t_analysisdata
 * Any changes made to the data through \p d can be accessed in later calls to
 * analyze_frame() as well as in the main function (gmx_template()).
 *
 * Here, we can do whatever calculations our program requires for a frame.
 * For the template, we first print the time for the current frame:
 * \skip  if (d->fp)
 * \until }
 * Then, we do a simple calculation:
 * \skip  nbsearch
 * \until }
 * \until }
 * \until }
 * \until }
 * After all the selections are processed, we print a newline to the output
 * file:
 * \skip  if (d->fp)
 * \until }
 * Finally, we return zero to indicate that all went well.
 * \skipline return 0;
 *
 *
 * \section gmx_template The main analysis tool
 *
 * We then define our main function in the same style as in Gromacs analysis
 * programs. We also provide a help text that can be shown with the \p -h
 * command-line option:
 * \skip  gmx_template
 * \until };
 *
 * Next, we define program-specific command-line arguments as in standard
 * Gromacs analysis programs:
 * \until };
 * Options for controlling the begin/end/skip of the trajectory are
 * automatically added, as are options for selections.
 * In the template, the second argument is for demonstartion purposes only;
 * it is not actually used.
 *
 * We also define the command-line options for output files:
 * \skip  t_filenm
 * \until NFILE
 * There should be no need to specify standard input files (trajectory,
 * topology, index), because they are automatically added and read in
 * by parse_trjana_args().
 * If you, however, define these file types with non-standard values, they
 * override the default values in parse_trjana_args().
 *
 * We then define some local variables:
 * \until gmx_ana_selection_t
 * The \p trj pointer holds internal data required for the library, and
 * can also be used to access topology information during initialization.
 * We also declare a t_analysisdata structure to hold data required in
 * analyze_frame(), as well as a few variables that will store information
 * about the selection that the user has provided.
 *
 * The actual code begins next. To follow the style of Gromacs tools, we
 * first print some information about the Gromacs build:
 * \skipline CopyRight
 *
 * Next, we initialize the \p trj pointer:
 * \skipline  gmx_ana_traj_create
 * We can use different flags to specify requirements for the selections and/or
 * other features of the library.
 * See \ref analysis_flags "Flags for gmx_ana_traj_create()".
 *
 * Next, we set some options:
 * \until nanagrps
 *
 * After initializing \p trj, we can call parse_trjana_args() to parse
 * command-line arguments and initialise the rest of the \p trj structure:
 * \skip  parse_trjana_args(trj,
 * \until oenv);
 *
 * After the call to parse_trjana_args(), all the command-line arguments are
 * available in the variables pointed by the \p pa array and can be used in
 * initialization.
 *
 * If you need to initialize the number of analysis groups
 * based on the command-line arguments, you can set \ref ANA_USER_SELINIT
 * and call gmx_ana_init_selections() separately.
 * Currently, the template does not demonstrate this to keep
 * it as simple as possible.
 *
 * After these calls, the program should perform any initialization that
 * is required before starting the analysis.
 * This can include allocating memory, precalculating things, preparing
 * output files and so on.
 * Any data that is required in the analysis of a single frame should be
 * put into \p d.
 * Data related to the selections, topology, and/or the first frame
 * can be accessed through \p trj; see the functions trajana.h.
 * For the template, we get the selection used as the reference positions
 * and the analysis groups using the following code:
 * \skip  gmx_ana_get_refsel
 * \until gmx_ana_get_anagrps
 * Notice that the reference selection is not included in the information
 * returned by the last two calls.
 *
 * For the template, we want to search the neighborhood of positions in
 * the first selection, so we initialize a neighborhood search.
 * This can be achieved with the following code:
 * \skip  nbsearch
 * \until }
 * Note that the calculation data is stored in \p d for it to be accessible in
 * analyze_frame().
 *
 * For the template, we also need some memory for calculating the average
 * distance for each index group, so we need to allocate it:
 * \skip  snew(d.ave
 * \until snew(d.n
 * We also open and prepare an output file for the data that is
 * written out during each frame:
 * \skip  d.fp
 * \until }
 *
 * Now we should be ready to do the actual loop over the frames.
 * This is accomplished by the function call
 * \skipline gmx_ana_do
 * This reads in each frame, updates the selections and calls analyze_frame()
 * for each frame.
 * All the initialized data (counters, file pointers, parameters etc.)
 * that is needed by the analysis (analyze_frame() in this case)
 * should be included in \p d.
 * Note that only the analysis groups are passed to analyze_frame();
 * the reference groups should be retrieved separately using
 * gmx_ana_get_refsel().
 * If this is not desired, one can also specify \ref ANA_USE_FULLGRPS to
 * include the reference groups to those passed to the analysis function.
 * Notice that this flag also changes the behavior of
 * gmx_ana_get_anagrps() and gmx_ana_get_grpnames().
 *
 * Finally, the analysis program probably needs to calculate and write out
 * some values such as averaged properties.
 * The template first closes the output file if one was opened
 * \skip  if (d.fp)
 * \until }
 * and then calculates and prints out the average distances for each analysis
 * group:
 * \skip  for (
 * \until }
 *
 * To follow the conventions of Gromacs tools, we finally print a (funny)
 * line
 * \skipline thanx
 * and return zero to indicate success.
 * \skipline return
 *
 *
 * \section main Definition of main()
 *
 * Now, the only thing remaining is to define the main() function.
 * It should simply call our gmx_template() function:
 * \skip int
 * \until }
 */
