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
/*! \file
 * \brief Main API of the trajectory analysis library.
 *
 * Contains the API for the core analysis library.
 *
 * \todo
 * Better handling of reference groups.
 * It would be nice to be able to provide a string that would be used in
 * prompting the groups, and also in automatic reporting of what the tool
 * is about to do.
 *
 * Most analysis tools should include trajana.h
 * (which automatically includes indexutil.h, selection.h, position.h)
 * and possibly one or more of the following headers:
 * displacement.h, histogram.h, nbsearch.h.
 * If the tool implements custom selection methods, it should also include
 * selmethod.h (which automatically includes selparam.h and selvalue.h).
 *
 * Other headers (centerofmass.h, poscalc.h) are used internally by the
 * library to calculate positions.
 * Analysis tools should preferably not use the routines in these headers
 * directly, but instead get all positions through selections. This makes
 * them more flexible with a minimal amount of work.
 */
#ifndef TRAJANA_H
#define TRAJANA_H
#include "visibility.h"
#include "typedefs.h"
#include "filenm.h"
#include "readinp.h"

#include "indexutil.h"
#include "selection.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Data structure for trajectory analysis tools. */
typedef struct gmx_ana_traj_t gmx_ana_traj_t;

/*! \name Flags for gmx_ana_traj_create()
 * \anchor analysis_flags
 * These flags can be used to alter the behavior of the analysis library to
 * suit the analysis tool.
 * They are given to the gmx_ana_traj_create() when creating the
 * \c gmx_ana_traj_t data structure, and affect the behavior of the other
 * functions in this header.
 */
/*@{*/
/*! \brief
 * Force loading of a topology file.
 *
 * If this flag is not specified, the topology file is loaded only if it is
 * provided on the command line explicitly.
 *
 * \see gmx_ana_get_topology()
 */
#define ANA_REQUIRE_TOP      (1<<0)
/*! \brief
 * Do no free the coordinates loaded from the topology.
 *
 * If this flag is specified, the coordinates loaded from the topology can
 * be accessed with gmx_ana_get_topconf().
 *
 * \see gmx_ana_get_topconf()
 */
#define ANA_USE_TOPX         (1<<1)
/*! \brief
 * Disallows the user from changing PBC handling.
 *
 * If this option is not specified, the analysis program (see gmx_analysisfunc())
 * may be passed a NULL PBC structure, and it should be able to handle such a
 * situation.
 */
#define ANA_NOUSER_PBC       (1<<4)
/*! \brief
 * Disallows the user from changing PBC removal.
 */
#define ANA_NOUSER_RMPBC     (1<<5)
/*! \brief
 * Disallows dynamic selections.
 *
 * If this flag is specified, an error is reported if the user specifies
 * any dynamic selections.
 */
#define ANA_NO_DYNSEL        (1<<8)
/*! \brief
 * Disallows breaking of residues in dynamic selections.
 *
 * Makes it impossible for the user to select atom-based dynamic selections.
 *
 * Only has effect if \ref ANA_NO_DYNSEL is not specified.
 */
#define ANA_REQUIRE_WHOLE    (1<<9)
/*! \brief
 * Disables automatic initialization of selection groups.
 *
 * If this flag is specified, parse_trjana_args() does not call
 * gmx_ana_init_selections(), allowing the program to do some initialization
 * before the call.
 * In particular, the program can use gmx_ana_set_nrefgprs() and
 * gmx_ana_set_nanagrps() before calling gmx_ana_init_selections() to
 * control the number of selections to expect.
 * This is useful if the program requires a different number of index groups
 * with different command-line arguments.
 * If the flag is specified, the program should call gmx_ana_init_selections()
 * exactly once after the parse_trjana_args() call.
 */
#define ANA_USER_SELINIT     (1<<10)
/*! \brief
 * Allow only atomic positions to be selected.
 *
 * Note that this flag only applies to the analysis groups, not the reference
 * groups. The reference groups should be checked in the analysis program
 * if some limitations are imposed on them.
 */
#define ANA_ONLY_ATOMPOS     (1<<11)
/*! \brief
 * Use masks in the positions to convey dynamic information.
 *
 * If this flag is specified, the positions calculated for the selections
 * are calculated always for the same group of atoms, even if the selection is
 * dynamic.
 * Dynamic selections only affect the \p refid field of the indexgroup map
 * given in the positions.
 */
#define ANA_USE_POSMASK      (1<<12)
/*! \brief
 * Pass the reference groups to gmx_analysisfunc().
 *
 * If this flag is specified, the reference groups are passed on to
 * gmx_analysisfunc().
 * Similarly, the arrays returned by gmx_ana_get_anagrps() and
 * gmx_ana_get_grpnames() contain the reference groups in the beginning.
 */
#define ANA_USE_FULLGRPS     (1<<13)
/*! \brief
 * Dump the parsed and compiled selection trees.
 *
 * This flag is used by internal debugging tools to make
 * gmx_ana_init_selections() dump the selection trees to stderr.
 */
#define ANA_DEBUG_SELECTION  (1<<16)
/*@}*/


/*! \name Functions for initialization */
/*@{*/

/** Allocates and initializes data structure for trajectory analysis. */
GMX_LIBGMX_EXPORT
int
gmx_ana_traj_create(gmx_ana_traj_t **data, unsigned long flags);
/** Frees the memory allocated for trajectory analysis data. */
void
gmx_ana_traj_free(gmx_ana_traj_t *d);
/** Sets additional flags after gmx_ana_traj_create() has been called. */
int
gmx_ana_add_flags(gmx_ana_traj_t *d, unsigned long flags);
/** Sets the number of reference groups required. */
GMX_LIBGMX_EXPORT
int
gmx_ana_set_nrefgrps(gmx_ana_traj_t *d, int nrefgrps);
/** Sets the number of analysis groups required. */
GMX_LIBGMX_EXPORT
int
gmx_ana_set_nanagrps(gmx_ana_traj_t *d, int nanagrps);
/** Sets whether PBC are used. */
int
gmx_ana_set_pbc(gmx_ana_traj_t *d, gmx_bool bPBC);
/** Sets whether molecules are made whole. */
int
gmx_ana_set_rmpbc(gmx_ana_traj_t *d, gmx_bool bRmPBC);
/** Sets flags that determine what to read from the trajectory. */
int
gmx_ana_set_frflags(gmx_ana_traj_t *d, int frflags);
/** Parses command-line arguments and performs some initialization. */
GMX_LIBGMX_EXPORT
int
parse_trjana_args(gmx_ana_traj_t *d, int *argc, char *argv[],
                  unsigned long pca_flags, int nfile, t_filenm fnm[],
                  int npargs, t_pargs *pa,
                  int ndesc, const char **desc,
                  int nbugs, const char **bugs,
                  output_env_t *oenv);
/** Initializes selection information. */
int
gmx_ana_init_selections(gmx_ana_traj_t *d);
/** Initializes calculation of covered fractions for selections. */
GMX_LIBGMX_EXPORT
int
gmx_ana_init_coverfrac(gmx_ana_traj_t *d, e_coverfrac_t type);

/** Returns whether PBC should be used. */
gmx_bool
gmx_ana_has_pbc(gmx_ana_traj_t *d);
/** Gets the topology information. */
GMX_LIBGMX_EXPORT
int
gmx_ana_get_topology(gmx_ana_traj_t *d, gmx_bool bReq, t_topology **top, gmx_bool *bTop);
/** Gets the configuration from the topology. */
int
gmx_ana_get_topconf(gmx_ana_traj_t *d, rvec **x, matrix box, int *ePBC);
/** Gets the first frame to be analyzed. */
int
gmx_ana_get_first_frame(gmx_ana_traj_t *d, t_trxframe **fr);

/** Gets the total number of selections provided by the user. */
int
gmx_ana_get_ngrps(gmx_ana_traj_t *d, int *ngrps);
/** Gets the number of analysis groups provided by the user. */
GMX_LIBGMX_EXPORT
int
gmx_ana_get_nanagrps(gmx_ana_traj_t *d, int *nanagrps);
/** Gets the selection object for a reference selection. */
GMX_LIBGMX_EXPORT
int
gmx_ana_get_refsel(gmx_ana_traj_t *d, int i, gmx_ana_selection_t **sel);
/** Gets the selection object for a reference selection. */
GMX_LIBGMX_EXPORT
int
gmx_ana_get_anagrps(gmx_ana_traj_t *d, gmx_ana_selection_t ***sel);
/** Gets an array of names for the selections. */
GMX_LIBGMX_EXPORT
int
gmx_ana_get_grpnames(gmx_ana_traj_t *d, char ***grpnames);
/** Gets the selection collection object that contains all the selections. */
int
gmx_ana_get_selcollection(gmx_ana_traj_t *d, gmx_ana_selcollection_t **sc);
/** Prints the selection strings into an XVGR file as comments. */
GMX_LIBGMX_EXPORT
int
xvgr_selections(FILE *out, gmx_ana_traj_t *d);

/*@}*/


/*! \name Functions for reading and analyzing the trajectory
 */
/*@{*/

/*! \brief
 * Function pointer type for frame analysis functions.
 *
 * \param[in] top  Topology structure.
 * \param[in] fr   Current frame.
 * \param[in] pbc  Initialized PBC structure for this frame.
 * \param[in] nr   Number of selections in the \p sel array.
 * \param[in] sel  Array of selections.
 * \param     data User data as provided to gmx_ana_do().
 * \returns   0 on success, a non-zero error code on error.
 *
 * This function is called by gmx_ana_do() for each frame that
 * needs to be analyzed.
 * Positions to be analyzed can be found in the \p sel array.
 * The selection array \p sel also provides information about the atoms that
 * have been used to evaluate the positions.
 * If a non-zero value is returned, gmx_ana_do() immediately exits and returns
 * the same value to the caller.
 */
typedef int (*gmx_analysisfunc)(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                                int nr, gmx_ana_selection_t *sel[], void *data);

/** Loops through all frames in the trajectory. */
GMX_LIBGMX_EXPORT
int
gmx_ana_do(gmx_ana_traj_t *d, int flags, gmx_analysisfunc analyze, void *data);
/** Gets the total number of frames analyzed. */
int
gmx_ana_get_nframes(gmx_ana_traj_t *d, int *nframes);

/*@}*/

#ifdef __cplusplus
}
#endif

#endif
