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
 * \brief API for handling selection (the \c gmx_ana_selection_t structure and related functions).
 *
 * There should be no need to use the data structures or call the
 * functions in this file directly unless using the selection routines outside
 * the main trajectory analysis API.
 */
#ifndef SELECTION_H
#define SELECTION_H

#include "position.h"
#include "indexutil.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Information for a collection of selections. */
typedef struct gmx_ana_selcollection_t gmx_ana_selcollection_t;

struct gmx_ana_poscalc_coll_t;

/** Defines the type of covered fraction. */
typedef enum
{
    CFRAC_NONE,         /**< No covered fraction (everything covered). */
    CFRAC_SOLIDANGLE    /**< Fraction of a solid (3D) angle covered. */
} e_coverfrac_t;

/*! \brief
 * Describes a single selection.
 */
typedef struct gmx_ana_selection_t
{
    /** Name of the selection. */
    char                       *name;
    /** The actual selection string. */
    char                       *selstr;
    /** Selected positions. */
    gmx_ana_pos_t               p;
    /** Masses associated with the positions. */
    real                       *m;
    /** Charges associated with the positions. */
    real                       *q;
    /** Pointer to the index group that holds the selected atoms. */
    struct gmx_ana_index_t     *g;
    /** TRUE if the value can change as a function of time. */
    gmx_bool                    bDynamic;
    /** Type of the covered fraction. */
    e_coverfrac_t               cfractype;
    /** TRUE if the covered fraction depends on the frame. */
    gmx_bool                    bCFracDyn;
    /** Covered fraction of the selection for the current frame. */
    real                        cfrac;
    /** The average covered fraction (over the trajectory). */
    real                        avecfrac;

    /*! \brief
     * Pointer to the root of the selection element tree (internal use only).
     *
     * \internal
     * This field is NULL if the selection has been loaded directly from an
     * index file.
     */
    struct t_selelem       *selelem;
    /** Original masses of all possible positions (internal use only). */
    real                   *orgm;
    /** Original charges of all possible positions (internal use only). */
    real                   *orgq;
} gmx_ana_selection_t;

/** Frees the memory allocated for a selection. */
void
gmx_ana_selection_free(gmx_ana_selection_t *sel);
/** Returns the name of a selection. */
char *
gmx_ana_selection_name(gmx_ana_selection_t *sel);
/** Prints out the selection information. */
void
gmx_ana_selection_print_info(gmx_ana_selection_t *sel);
/** Initializes the information for covered fraction. */
gmx_bool
gmx_ana_selection_init_coverfrac(gmx_ana_selection_t *sel, e_coverfrac_t type);

/** Creates a new empty selection collection. */
int
gmx_ana_selcollection_create(gmx_ana_selcollection_t      **sc,
                             struct gmx_ana_poscalc_coll_t *pcc);
/** Frees the memory allocated for a selection collection. */
void
gmx_ana_selcollection_free(gmx_ana_selcollection_t *sc);
/** Sets the default reference position handling for a selection collection. */
void
gmx_ana_selcollection_set_refpostype(gmx_ana_selcollection_t *sc, const char *type);
/** Sets the default output position handling for a selection collection. */
void
gmx_ana_selcollection_set_outpostype(gmx_ana_selcollection_t *sc,
                                     const char *type, gmx_bool bMaskOnly);
/** Request evaluation of velocities for selections. */
void
gmx_ana_selcollection_set_veloutput(gmx_ana_selcollection_t *sc,
                                    gmx_bool                 bVelOut);
/** Request evaluation of forces for selections. */
void
gmx_ana_selcollection_set_forceoutput(gmx_ana_selcollection_t *sc,
                                      gmx_bool                 bForceOut);
/** Sets the topology for a selection collection. */
int
gmx_ana_selcollection_set_topology(gmx_ana_selcollection_t *sc, t_topology *top,
                                   int natoms);
/** Returns the number of selections specified by a selection collection. */
int
gmx_ana_selcollection_get_count(gmx_ana_selcollection_t *sc);
/** Returns a selection by index. */
gmx_ana_selection_t *
gmx_ana_selcollection_get_selection(gmx_ana_selcollection_t *sc, int i);
/** Returns TRUE if the collection requires topology information for evaluation. */
gmx_bool
gmx_ana_selcollection_requires_top(gmx_ana_selcollection_t *sc);
/** Prints a human-readable version of the internal selection element tree. */
void
gmx_ana_selcollection_print_tree(FILE *fp, gmx_ana_selcollection_t *sc, gmx_bool bValues);
/** Prints the selection strings into an XVGR file as comments. */
void
xvgr_selcollection(FILE *fp, gmx_ana_selcollection_t *sc,
                   const output_env_t oenv);

/* In parsetree.c */
/** Parses selection(s) from standard input. */
int
gmx_ana_selcollection_parse_stdin(gmx_ana_selcollection_t *sc, int nr,
                                  gmx_ana_indexgrps_t *grps,
                                  gmx_bool bInteractive);
/** Parses selection(s) from a file. */
int
gmx_ana_selcollection_parse_file(gmx_ana_selcollection_t *sc, const char *fnm,
                                 gmx_ana_indexgrps_t *grps);
/** Parses selection(s) from a string. */
int
gmx_ana_selcollection_parse_str(gmx_ana_selcollection_t *sc, const char *str,
                                gmx_ana_indexgrps_t *grps);

/* In compiler.c */
/** Set debugging flag for selection compilation. */
void
gmx_ana_selcollection_set_compile_debug(gmx_ana_selcollection_t *sc, gmx_bool bDebug);
/** Prepares the selections for evaluation and performs some optimizations. */
int
gmx_ana_selcollection_compile(gmx_ana_selcollection_t *sc);

/* In evaluate.c */
/** Evaluates the selection. */
int
gmx_ana_selcollection_evaluate(gmx_ana_selcollection_t *sc,
                               t_trxframe *fr, t_pbc *pbc);
/** Evaluates the largest possible index groups from dynamic selections. */
int
gmx_ana_selcollection_evaluate_fin(gmx_ana_selcollection_t *sc, int nframes);

#ifdef __cplusplus
}
#endif

#endif
