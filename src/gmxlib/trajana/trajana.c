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
/*! \page libtrajana Library for trajectory analysis
 *
 * This is a trajectory analysis library for Gromacs.
 *
 * The main features of the library are:
 *  - \subpage selengine "Flexible handling of textual selections" that can
 *    also be dynamic, i.e., depend of the trajectory frame through
 *    positions of the particles.
 *    Selections evaluate directly to positions, which can then be used in
 *    the analysis program.
 *  - \subpage selmethods "Custom selection keywords" can also be easily
 *    implemented, also on a tool-by-tool basis.
 *  - Efficient \subpage nbsearch "neighborhood searching"
 *    (currently not very efficient...).
 *  - \subpage displacements "On-the-fly calculation of displacement vectors"
 *    during a single pass over the trajectory.
 *  - \subpage histograms "Calculation of histograms" with error estimates
 *    through block averaging.
 *
 * The functions also automate several things common to most analysis programs
 * such as making molecules whole if required, reading input files, and
 * setting up index groups for analysis.
 * The library unifies the structure of analysis programs (at least a bit)
 * and makes it easier to add common functionality to all analysis programs.
 *
 *
 * \section main_using Using the library
 *
 * There is a \ref share/template/template.c "template" for writing
 * analysis programs, the documentation for it and links from there should
 * help getting started.
 *
 *
 * \internal
 * \section libtrajana_impl Implementation details
 *
 * Some internal implementation details of the library are documented on
 * separate pages:
 *  - \subpage poscalcengine
 *  - \subpage selparser
 *  - \subpage selcompiler
 */
/*! \page selengine Text-based selections
 *
 * \section selection_basics Basics
 *
 * Selections are enabled automatically for an analysis program that uses
 * the library. The selection syntax is described in an online help that is
 * accessible from all tools that use the library.
 * By default, dynamic selections are allowed, and the user can freely
 * choose whether to analyze atoms or centers of mass/geometry of
 * residues/molecules.
 * These defaults, as well as some others, can be changed by specifying
 * flags for gmx_ana_traj_create().
 *
 * The analysis program can then access the selected positions for each frame
 * through a \p gmx_ana_selection_t array that is passed to the frame
 * analysis function (see gmx_analysisfunc()).
 * As long as the analysis program is written such that it does not assume
 * that the number of positions or the atoms in the groups groups remain
 * constant, any kind of selection expression can be used.
 *
 * Some analysis programs may require a special structure for the index groups
 * (e.g., \c g_angle requires the index group to be made of groups of three or
 * four atoms).
 * For such programs, it is up to the user to provide a proper selection
 * expression that always returns such positions.
 * Such analysis program can define \ref ANA_REQUIRE_WHOLE to make the
 * default behavior appropriate for the most common uses where the groups
 * should consist of atoms within a single residue/molecule.
 *
 * \section selection_methods Implementing new keywords
 *
 * New selection keywords can be easily implemented, either directly into the
 * library or as part of analysis programs (the latter may be useful for
 * testing or methods very specific to some analysis).
 * For both cases, you should first create a \c gmx_ana_selmethod_t structure
 * and fill it with the necessary information.
 * Details can be found on a separate page: \ref selmethods.
 */
/*! \internal \file
 * \brief Implementation of functions in trajana.h.
 */
/*! \internal \dir src/gmxlib/trajana
 * \brief Source code for common trajectory analysis functions.
 *
 * Selection handling is found in \ref src/gmxlib/selection.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include <filenm.h>
#include <futil.h>
#include <macros.h>
#include <pbc.h>
#include <rmpbc.h>
#include <smalloc.h>
#include <statutil.h>
#include <typedefs.h>
#include <tpxio.h>
#include <vec.h>

#include <poscalc.h>
#include <selection.h>
#include <selmethod.h>
#include <trajana.h>

/*! \internal \brief
 * Data structure for trajectory analysis tools.
 */
struct gmx_ana_traj_t
{
    /*! \brief
     * Flags that alter the behavior of the analysis library.
     *
     * This variable stores the flags passed to gmx_ana_traj_create() for use
     * in the other functions.
     */
    unsigned long             flags;
    /** Number of input reference groups. */
    int                       nrefgrps;
    /*! \brief
     * Number of input analysis groups.
     *
     * This is the number of groups in addition to the reference groups
     * that are required.
     * If -1, any number of groups (at least one) is acceptable.
     */
    int                       nanagrps;
    /*! \brief
     * Flags for init_first_frame() to specify what to read from the
     * trajectory.
     */
    int                           frflags;
    /** TRUE if molecules should be made whole for each frame. */
    gmx_bool                      bRmPBC;
    /*! \brief
     * TRUE if periodic boundary conditions should be used.
     *
     * If the flag is FALSE, the \p pbc pointer passed to gmx_analysisfunc()
     * is NULL.
     */
    gmx_bool                      bPBC;

    /** Name of the trajectory file (NULL if not provided). */
    char                     *trjfile;
    /** Name of the topology file (NULL if no topology loaded). */
    char                     *topfile;
    /** Non-NULL name of the topology file. */
    char                     *topfile_notnull;
    /** Name of the index file (NULL if no index file provided). */
    char                     *ndxfile;
    /** Name of the selection file (NULL if no file provided). */
    char                     *selfile;
    /** The selection string (NULL if not provided). */
    char                     *selection;

    /** The topology structure, or \p NULL if no topology loaded. */
    t_topology                   *top;
    /** TRUE if full tpx file was loaded, FALSE otherwise. */
    gmx_bool                      bTop;
    /** Coordinates from the topology (see \p bTopX). */
    rvec                         *xtop;
    /** The box loaded from the topology file. */
    matrix                        boxtop;
    /** The ePBC field loaded from the topology file. */
    int                           ePBC;

    /** The current frame, or \p NULL if no frame loaded yet. */
    t_trxframe               *fr;
    /** Used to store the status variable from read_first_frame(). */
    t_trxstatus              *status;
    /** The number of frames read. */
    int                       nframes;

    /** Number of elements in the \p sel array. */
    int                       ngrps;
    /*! \brief
     * Array of selection information (one element for each index group).
     *
     * After the call to gmx_ana_init_selections(), this array contains
     * information about the selections the user has provided.
     * The array contains \p ngrps elements;
     * if \p nanagrps was -1, the number may vary.
     * See \c gmx_ana_selection_t for details of the contents.
     *
     * The largest possible index groups for dynamic selections can be found
     * in \p sel[i]->g, i.e., the program can assume that any index group
     * passed to gmx_analysisfunc() is a subset of the provided group.
     * After gmx_ana_do(), the same groups can be used in post-processing.
     */
    gmx_ana_selection_t     **sel;
    /** Array of names of the selections for convenience. */
    char                    **grpnames;
    /** Position calculation data. */
    gmx_ana_poscalc_coll_t   *pcc;
    /** Selection data. */
    gmx_ana_selcollection_t  *sc;

    /** Data for statutil.c utilities. */
    output_env_t              oenv;
};

/** Loads the topology. */
static int load_topology(gmx_ana_traj_t *d, gmx_bool bReq);
/** Loads the first frame and does some checks. */
static int init_first_frame(gmx_ana_traj_t *d);

static int add_fnmarg(int nfile, t_filenm *fnm, t_filenm *fnm_add)
{
    memcpy(&(fnm[nfile]), fnm_add, sizeof(*fnm_add));
    return nfile + 1;
}

/* Copied from src/gmxlib/statutil.c */
static int
add_parg(int npargs, t_pargs *pa, t_pargs *pa_add)
{
    memcpy(&(pa[npargs]), pa_add, sizeof(*pa_add));
    return npargs + 1;
}

/*!
 * \param[out] data  Trajectory analysis data structure poitner to initialize.
 * \param[in]  flags Combination of flags (see \ref analysis_flags).
 * \returns    0 on success.
 */
int
gmx_ana_traj_create(gmx_ana_traj_t **data, unsigned long flags)
{
    gmx_ana_traj_t     *d;
    int                 rc;

    snew(d, 1);

    d->nrefgrps        = 0;
    d->nanagrps        = 1;
    d->frflags         = TRX_NEED_X;
    d->bRmPBC          = TRUE;
    d->bPBC            = TRUE;

    d->trjfile         = NULL;
    d->topfile         = NULL;
    d->ndxfile         = NULL;
    d->selfile         = NULL;
    d->selection       = NULL;

    d->top             = NULL;
    d->bTop            = FALSE;
    d->xtop            = NULL;
    d->ePBC            = -1;
    d->fr              = NULL;
    d->nframes         = 0;

    d->ngrps           = 0;
    d->sel             = NULL;
    d->grpnames        = NULL;

    d->flags           = flags;
    d->topfile_notnull = NULL;
    rc                 = gmx_ana_poscalc_coll_create(&d->pcc);
    if (rc != 0)
    {
        sfree(d);
        *data = NULL;
        return rc;
    }
    rc = gmx_ana_selcollection_create(&d->sc, d->pcc);
    if (rc != 0)
    {
        gmx_ana_poscalc_coll_free(d->pcc);
        sfree(d);
        *data = NULL;
        return rc;
    }
    d->status          = NULL;
    d->oenv            = NULL;

    *data              = d;
    return 0;
}

/*!
 * \param   d  Trajectory analysis data to free.
 */
void
gmx_ana_traj_free(gmx_ana_traj_t *d)
{
    int                 i;

    sfree(d->trjfile);
    sfree(d->topfile);
    sfree(d->topfile_notnull);
    sfree(d->ndxfile);
    sfree(d->selfile);
    if (d->top)
    {
        done_top(d->top);
        sfree(d->top);
    }
    if (d->fr)
    {
        /* Gromacs does not seem to have a function for freeing frame data */
        sfree(d->fr->x);
        sfree(d->fr->v);
        sfree(d->fr->f);
        sfree(d->fr);
    }
    sfree(d->xtop);
    sfree(d->sel);
    gmx_ana_selcollection_free(d->sc);
    gmx_ana_poscalc_coll_free(d->pcc);
    sfree(d->grpnames);
    output_env_done(d->oenv);
    sfree(d);
}

/*!
 * \param[in,out] d      Trajectory analysis data structure.
 * \param[in]     flags  Additional flags to set.
 * \returns       0 on success, a non-zero error code on error.
 *
 * Currently there are no checks whether the flags make sense or not.
 */
int
gmx_ana_add_flags(gmx_ana_traj_t *d, unsigned long flags)
{
    d->flags |= flags;
    return 0;
}

/*!
 * \param[in,out] d      Trajectory analysis data structure.
 * \param[in]     bPBC   TRUE if periodic boundary conditions should be used.
 * \returns       0 on success.
 *
 * If this function is called before parse_trjana_args(), it sets the default
 * for whether PBC are used in the analysis or not.
 * If \ref ANA_NOUSER_PBC is not set, a command-line option is provided for the
 * user to override the default value.
 * If called after parse_trjana_args(), it overrides the setting provided by
 * the user or an earlier call.
 *
 * If this function is not called, the default is to use PBC.
 *
 * If PBC are not used, the \p pbc pointer passed to gmx_analysisfunc()
 * is NULL.
 *
 * \see \ref ANA_NOUSER_PBC
 */
int
gmx_ana_set_pbc(gmx_ana_traj_t *d, gmx_bool bPBC)
{
    d->bPBC = bPBC;
    return 0;
}

/*!
 * \param[in] d      Trajectory analysis data structure.
 * \returns   TRUE if periodic boundary conditions are set to be used.
 */
gmx_bool
gmx_ana_has_pbc(gmx_ana_traj_t *d)
{
    return d->bPBC;
}

/*!
 * \param[in,out] d      Trajectory analysis data structure.
 * \param[in]     bRmPBC TRUE if molecules should be made whole.
 * \returns       0 on success.
 *
 * If this function is called before parse_trjana_args(), it sets the default
 * for whether molecules are made whole.
 * If \ref ANA_NOUSER_RMPBC is not set, a command-line option is provided for
 * the user to override the default value.
 * If called after parse_trjana_args(), it overrides the setting provided by
 * the user or an earlier call.
 *
 * If this function is not called, the default is to make molecules whole.
 *
 * The main use of this function is to call it with FALSE if your analysis
 * program does not require whole molecules as this can increase the
 * performance.
 * In such a case, you can also specify \ref ANA_NOUSER_RMPBC to not to
 * confuse the user with an option that would only slow the program
 * down.
 *
 * \see \ref ANA_NOUSER_RMPBC
 */
int
gmx_ana_set_rmpbc(gmx_ana_traj_t *d, gmx_bool bRmPBC)
{
    d->bRmPBC = bRmPBC;
    return 0;
}

/*!
 * \param[in,out] d       Trajectory analysis data structure.
 * \param[in]     frflags Flags for what to read from the trajectory file.
 * \returns       0 on success, an error code on error.
 *
 * The TRX_NEED_X flag is always set.
 * If the analysis tools needs some other information (velocities, forces),
 * it can call this function to load additional information from the
 * trajectory.
 */
int
gmx_ana_set_frflags(gmx_ana_traj_t *d, int frflags)
{
    if (d->sel)
    {
        gmx_call("cannot set trajectory flags after initializing selections");
        return -1;
    }
    if (d->fr)
    {
        gmx_call("cannot set trajectory flags after the first frame has been read");
        return -1;
    }
    frflags   |= TRX_NEED_X;
    d->frflags = frflags;
    return 0;
}

/*!
 * \param[in,out] d        Trajectory analysis data structure.
 * \param[in]     nrefgrps Number of reference groups required.
 * \returns       0 on success, a non-zero error code on error.
 *
 * \p nrefgrps should be a non-negative integer.
 * If this function is not called (or \p nrefgrps is 0), all selections are
 * treated as reference groups.
 */
int
gmx_ana_set_nrefgrps(gmx_ana_traj_t *d, int nrefgrps)
{
    if (nrefgrps < 0)
    {
        d->nrefgrps = 0;
        gmx_incons("number of reference groups is negative");
        return EINVAL;
    }
    d->nrefgrps = nrefgrps;
    return 0;
}

/*!
 * \param[in,out] d        Trajectory analysis data structure.
 * \param[in]     nanagrps Number of analysis groups required
 *   (-1 stands for any number of analysis groups).
 * \returns       0 on success, a non-zero error code on error.
 *
 * \p nanagrps should be a positive integer or -1.
 * In the latter case, any number of groups (but at least one) is acceptable.
 * gmx_ana_get_nanagrps() can be used to access the actual value after
 * gmx_ana_init_selections() has been called.
 * If this function is not called, a single analysis group is expected.
 */
int
gmx_ana_set_nanagrps(gmx_ana_traj_t *d, int nanagrps)
{
    if (nanagrps <= 0 && nanagrps != -1)
    {
        d->nanagrps = 1;
        gmx_incons("number of analysis groups is invalid");
        return EINVAL;
    }
    d->nanagrps = nanagrps;
    return 0;
}

/*!
 * \param[in]  d     Trajectory analysis data structure.
 * \param[out] ngrps Total number of selections specified by the user.
 * \returns    0 on success, a non-zero error code on error.
 *
 * The value stored in \p *ngrps is the sum of the number of reference groups
 * and the number of analysis groups.
 * If a specific number (not -1) of analysis groups has been set with
 * gmx_ana_set_nanagrps(), the value is always the sum of the values provided
 * to gmx_ana_set_nrefgrps() and gmx_ana_set_nanagrps().
 *
 * Should only be called after gmx_ana_init_selections().
 */
int
gmx_ana_get_ngrps(gmx_ana_traj_t *d, int *ngrps)
{
    if (d->nanagrps == -1)
    {
        *ngrps = 0;
        gmx_call("gmx_ana_init_selections() not called");
        return EINVAL;
    }
    *ngrps = d->nrefgrps + d->nanagrps;
    return 0;
}

/*!
 * \param[in]  d        Trajectory analysis data structure.
 * \param[out] nanagrps Number of analysis groups specified by the user.
 * \returns    0 on success, a non-zero error code on error.
 *
 * If a specific number (not -1) of analysis groups has been set with
 * gmx_ana_set_nanagrps(), the value is always the same value.
 * Hence, you only need to call this function if gmx_ana_set_nanagrps() has
 * been called with \p nanagrps set to -1.
 *
 * Should only be called after gmx_ana_init_selections().
 */
int
gmx_ana_get_nanagrps(gmx_ana_traj_t *d, int *nanagrps)
{
    if (d->nanagrps == -1)
    {
        *nanagrps = 0;
        gmx_call("gmx_ana_init_selections() not called");
        return EINVAL;
    }
    *nanagrps = d->nanagrps;
    return 0;
}

/*!
 * \param[in]  d   Trajectory analysis data structure.
 * \param[in]  i   Ordinal number of the reference selection to get.
 * \param[out] sel Selection object for the \p i'th reference group.
 * \returns    0 on success, a non-zero error code on error.
 *
 * The pointer returned in \p *sel should not be freed.
 * Should only be called after gmx_ana_init_selections().
 */
int
gmx_ana_get_refsel(gmx_ana_traj_t *d, int i, gmx_ana_selection_t **sel)
{
    if (i < 0 || i >= d->nrefgrps)
    {
        *sel = NULL;
        gmx_call("invalid reference group number");
        return EINVAL;
    }
    *sel = gmx_ana_selcollection_get_selection(d->sc, i);
    if (!*sel)
    {
        gmx_incons("gmx_ana_init_selections() not called");
        return EINVAL;
    }
    return 0;
}

/*!
 * \param[in]  d   Trajectory analysis data structure.
 * \param[out] sel Array of analysis selections.
 * \returns    0 on success, a non-zero error code on error.
 *
 * The pointer returned in \p *sel should not be freed.
 * Should only be called after gmx_ana_init_selections().
 */
int
gmx_ana_get_anagrps(gmx_ana_traj_t *d, gmx_ana_selection_t ***sel)
{
    if (!d->sel)
    {
        *sel = NULL;
        gmx_incons("gmx_ana_init_selections() not called");
        return EINVAL;
    }
    *sel = d->sel;
    return 0;
}

/*!
 * \param[in]  d        Trajectory analysis data structure.
 * \param[out] grpnames Array of selection names.
 * \returns    0 on success, a non-zero error code on error.
 *
 * The pointer returned in \p *grpnames should not be freed.
 * Should only be called after gmx_ana_init_selections().
 */
int
gmx_ana_get_grpnames(gmx_ana_traj_t *d, char ***grpnames)
{
    if (!d->grpnames)
    {
        *grpnames = NULL;
        gmx_call("gmx_ana_init_selections() not called");
        return EINVAL;
    }
    *grpnames = d->grpnames;
    return 0;
}

/*!
 * \param[in]  d   Trajectory analysis data structure.
 * \param[out] sc  Selection collection object.
 * \returns    0 on success.
 *
 * This function is mostly useful for debugging purposes.
 * The information commonly required in analysis programs is accessible
 * more conveniently through other means.
 *
 * The pointer returned in \p *sc should not be freed.
 * Can be called at any point.
 */
int
gmx_ana_get_selcollection(gmx_ana_traj_t *d, gmx_ana_selcollection_t **sc)
{
    *sc = d->sc;
    return 0;
}

/*!
 * \param[in,out] d     Trajectory analysis data structure.
 * \returns    0 on success, a non-zero error code on error.
 *
 * This function should be called first in the analysis program, much in
 * the same way as parse_common_args() in traditional Gromacs analysis
 * programs. It adds some command-line arguments of its own, and uses
 * parse_common_args() to parse the command-line arguments.
 * It also loads topology information if required or if a topology is
 * provided on the command line.
 * Selection handling is also initialized if it is enabled and
 * the user has selected it on the command line.
 *
 * The rest of the parameters are passed on to the Gromacs routine
 * parse_common_args(), and the use of this function should be identical
 * to parse_common_args(), with the exception that for \p pca_flags,
 * \p PCA_CAN_TIME and \p PCA_BE_NICE flags are automatically added.
 * \param      argc
 * \param      argv
 * \param      pca_flags
 * \param      nfile
 * \param      fnm
 * \param      npargs
 * \param      pa
 * \param      ndesc
 * \param      desc
 * \param      nbugs
 * \param      bugs
 * \param      oenv
 */
int
parse_trjana_args(gmx_ana_traj_t *d,
                  int *argc, char *argv[], unsigned long pca_flags,
                  int nfile, t_filenm fnm[], int npargs, t_pargs *pa,
                  int ndesc, const char **desc,
                  int nbugs, const char **bugs,
                  output_env_t *oenv)
{
    t_filenm               *all_fnm = NULL;
    int                     max_fnm, nfall;
    int                    *fnm_map;
    t_pargs                *all_pa = NULL;
    int                     max_pa, npall;
    size_t                  i;
    int                     k;
    int                     rc;
    const char             *tmp_fnm;

    t_filenm                def_fnm[] = {
        {efTRX, NULL,  NULL,        ffOPTRD},
        {efTPS, NULL,  NULL,        ffREAD},
        {efDAT, "-sf", "selection", ffOPTRD},
        {efNDX, NULL,  NULL,        ffOPTRD},
    };
    gmx_bool                bPBC     = TRUE;
    t_pargs                 pbc_pa[] = {
        {"-pbc",      FALSE, etBOOL, {&bPBC},
         "Use periodic boundary conditions for distance calculation"},
    };
    gmx_bool                bRmPBC     = TRUE;
    t_pargs                 rmpbc_pa[] = {
        {"-rmpbc",    FALSE, etBOOL, {&bRmPBC},
         "Make molecules whole for each frame"},
    };
    char                   *selection = NULL;
    const char            **rpost     = NULL;
    gmx_bool                bSelDump  = FALSE;
    t_pargs                 sel_pa[]  = {
        {"-select",   FALSE, etSTR,  {&selection},
         "Selection string (use 'help' for help). Note that the "
         "whole selection string will need to be quoted so that "
         "your shell will pass it in as a string. Example: "
         "[TT]g_select -select '\"Nearby water\" resname SOL "
         "and within 0.25 of group Protein'[tt]"},
        {"-seldebug", FALSE, etBOOL, {&bSelDump},
         "HIDDENPrint out the parsed and compiled selection trees"},
    };
    t_pargs                 dsel_pa[] = {
        {"-selrpos",  FALSE, etENUM, {NULL},
         "Selection reference position"},
    };
    const char            **spost      = NULL;
    t_pargs                 selpt_pa[] = {
        {"-seltype",  FALSE, etENUM, {NULL},
         "Default analysis positions"},
    };
#define MAX_PA asize(sel_pa)+asize(dsel_pa)+5

    if (d->nrefgrps < 0)
    {
        gmx_incons("number of reference groups is negative");
        return EINVAL;
    }

    if (d->flags & ANA_DEBUG_SELECTION)
    {
        bSelDump = TRUE;
    }

    rpost = gmx_ana_poscalc_create_type_enum(!(d->flags & ANA_REQUIRE_WHOLE));
    if (rpost == NULL)
    {
        return ENOMEM;
    }
    spost = gmx_ana_poscalc_create_type_enum(TRUE);
    if (spost == NULL)
    {
        sfree(rpost);
        return ENOMEM;
    }

    /* Construct the file name argument array */
    max_fnm = nfile + asize(def_fnm);
    snew(all_fnm, max_fnm);
    nfall = 0;
    if (!(d->flags & ANA_REQUIRE_TOP))
    {
        def_fnm[1].flag |= ffOPT;
    }
    snew(fnm_map, nfile);
    for (k = 0; k < nfile; ++k)
    {
        fnm_map[k] = -1;
    }

    for (i = 0; i < asize(def_fnm); ++i)
    {
        for (k = 0; k < nfile; ++k)
        {
            if (fnm_map[k] == -1 && def_fnm[i].opt == NULL
                && fnm[k].opt == NULL && fnm[k].ftp == def_fnm[i].ftp)
            {
                break;
            }
        }
        if (k < nfile)
        {
            fnm_map[k] = nfall;
            nfall      = add_fnmarg(nfall, all_fnm, &(fnm[k]));
        }
        else
        {
            nfall = add_fnmarg(nfall, all_fnm, &(def_fnm[i]));
        }
    }

    for (k = 0; k < nfile; ++k)
    {
        if (fnm_map[k] == -1)
        {
            fnm_map[k] = nfall;
            nfall      = add_fnmarg(nfall, all_fnm, &(fnm[k]));
        }
    }

    /* Construct the argument array */
    max_pa = npargs + MAX_PA;
    snew(all_pa, max_pa);
    npall = 0;

    if (!(d->flags & ANA_NOUSER_RMPBC))
    {
        for (i = 0; i < asize(rmpbc_pa); ++i)
        {
            npall = add_parg(npall, all_pa, &(rmpbc_pa[i]));
        }
    }
    if (!(d->flags & ANA_NOUSER_PBC))
    {
        for (i = 0; i < asize(pbc_pa); ++i)
        {
            npall = add_parg(npall, all_pa, &(pbc_pa[i]));
        }
    }

    for (i = 0; i < asize(sel_pa); ++i)
    {
        npall = add_parg(npall, all_pa, &(sel_pa[i]));
    }
    if (!(d->flags & ANA_NO_DYNSEL))
    {
        dsel_pa[0].u.c = rpost;
        for (i = 0; i < asize(dsel_pa); ++i)
        {
            npall = add_parg(npall, all_pa, &(dsel_pa[i]));
        }
    }

    if (!(d->flags & ANA_ONLY_ATOMPOS))
    {
        selpt_pa[0].u.c = spost;
        for (i = 0; i < asize(selpt_pa); ++i)
        {
            npall = add_parg(npall, all_pa, &(selpt_pa[i]));
        }
    }

    for (k = 0; k < npargs; ++k)
    {
        npall = add_parg(npall, all_pa, &(pa[k]));
    }

    pca_flags |= PCA_CAN_TIME | PCA_BE_NICE;
    parse_common_args(argc, argv, pca_flags, nfall, all_fnm, npall, all_pa,
                      ndesc, desc, nbugs, bugs, oenv);
    d->oenv = *oenv;

    /* Process our own options.
     * Make copies of file names for easier memory management. */
    tmp_fnm            = ftp2fn_null(efTRX, nfall, all_fnm);
    d->trjfile         = tmp_fnm ? strdup(tmp_fnm) : NULL;
    tmp_fnm            = ftp2fn_null(efTPS, nfall, all_fnm);
    d->topfile         = tmp_fnm ? strdup(tmp_fnm) : NULL;
    d->topfile_notnull = strdup(ftp2fn(efTPS, nfall, all_fnm));
    tmp_fnm            = ftp2fn_null(efNDX, nfall, all_fnm);
    d->ndxfile         = tmp_fnm ? strdup(tmp_fnm) : NULL;
    if (!(d->flags & ANA_NOUSER_RMPBC))
    {
        d->bRmPBC      = bRmPBC;
    }
    if (!(d->flags & ANA_NOUSER_PBC))
    {
        d->bPBC        = bPBC;
    }
    d->selection       = selection;
    tmp_fnm            = opt2fn_null("-sf", nfall, all_fnm);
    d->selfile         = tmp_fnm ? strdup(tmp_fnm) : NULL;

    /* Copy the results back */
    for (k = 0; k < nfile; ++k)
    {
        memcpy(&(fnm[k]), &(all_fnm[fnm_map[k]]), sizeof(fnm[k]));
        /* Delegate responsibility of freeing the file names to caller. */
        all_fnm[fnm_map[k]].nfiles = 0;
        all_fnm[fnm_map[k]].fns    = NULL;
    }
    for (i = 0, k = npall - npargs; i < (size_t)npargs; ++i, ++k)
    {
        memcpy(&(pa[i]), &(all_pa[k]), sizeof(pa[i]));
    }

    /* Free memory we have allocated. */
    done_filenms(nfall, all_fnm);
    sfree(all_fnm);
    sfree(fnm_map);
    sfree(all_pa);

    if (!(d->flags & ANA_NO_DYNSEL))
    {
        gmx_ana_selcollection_set_refpostype(d->sc, rpost[0]);
    }
    else
    {
        gmx_ana_selcollection_set_refpostype(d->sc, rpost[1]);
    }
    sfree(rpost);
    if (bSelDump)
    {
        d->flags |= ANA_DEBUG_SELECTION;
    }
    else
    {
        d->flags &= ~ANA_DEBUG_SELECTION;
    }

    if (!(d->flags & ANA_ONLY_ATOMPOS))
    {
        gmx_ana_selcollection_set_outpostype(d->sc, spost[0], d->flags & ANA_USE_POSMASK);
    }
    else
    {
        gmx_ana_selcollection_set_outpostype(d->sc, spost[1], d->flags & ANA_USE_POSMASK);
    }
    sfree(spost);

    /* Check if the user requested help on selections.
     * If so, call gmx_ana_init_selections() to print the help and exit. */
    if (selection && strncmp(selection, "help", 4) == 0)
    {
        gmx_ana_init_selections(d);
    }

    /* If no trajectory file is given, we need to set some flags to be able
     * to prepare a frame from the loaded topology information. Also, check
     * that a topology is provided. */
    if (!d->trjfile)
    {
        if (!d->topfile)
        {
            gmx_input("No trajectory or topology provided, nothing to do!");
            return -1;
        }
        d->flags |= ANA_REQUIRE_TOP;
        d->flags |= ANA_USE_TOPX;
    }

    /* Load the topology if so requested. */
    rc = load_topology(d, (d->flags & ANA_REQUIRE_TOP));
    if (rc != 0)
    {
        return rc;
    }

    /* Initialize the selections/index groups */
    if (!(d->flags & ANA_USER_SELINIT))
    {
        rc = gmx_ana_init_selections(d);
    }

    return rc;
}

/*!
 * \param[in,out] d     Trajectory analysis data structure.
 * \param[in]     bReq  If TRUE, topology loading is forced.
 * \returns       0 on success, a non-zero error code on error.
 *
 * Initializes the \c gmx_ana_traj_t::top, \c gmx_ana_traj_t::bTop,
 * \c gmx_ana_traj_t::boxtop and \c gmx_ana_traj_t::ePBC fields of the
 * analysis structure.
 * If \p bReq is TRUE, the topology is loaded even if it is not given on
 * the command line.
 *
 * The function can be called multiple times safely; extra calls are
 * ignored.
 */
static int load_topology(gmx_ana_traj_t *d, gmx_bool bReq)
{
    char                title[STRLEN];

    /* Return if topology already loaded */
    if (d->top)
    {
        return 0;
    }

    if (d->topfile || bReq)
    {
        snew(d->top, 1);
        d->bTop = read_tps_conf(d->topfile_notnull, title, d->top,
                                &d->ePBC, &d->xtop, NULL, d->boxtop, TRUE);
        if (!(d->flags & ANA_USE_TOPX))
        {
            sfree(d->xtop);
            d->xtop = NULL;
        }
    }
    return 0;
}

/*!
 * \param[in]  d     Trajectory analysis data structure.
 * \param[in]  bReq  If TRUE, topology loading is forced.
 * \param[out] top   Topology data pointer to initialize.
 * \param[out] bTop  TRUE if a full tpx file was loaded, FALSE otherwise
 *   (can be NULL, in which case it is not used).
 * \returns    0 on success, a non-zero error code on error.
 *
 * If \ref ANA_REQUIRE_TOP has not been specified and \p bReq is FALSE,
 * the pointer stored in \p *top may be NULL if no topology has been provided
 * on the command line.
 *
 * The pointer returned in \p *top should not be freed.
 */
int
gmx_ana_get_topology(gmx_ana_traj_t *d, gmx_bool bReq, t_topology **top, gmx_bool *bTop)
{
    int rc;

    rc = load_topology(d, bReq);
    if (rc != 0)
    {
        *top = NULL;
        return rc;
    }
    *top = d->top;
    if (bTop)
    {
        *bTop = d->bTop;
    }
    return 0;
}

/*!
 * \param[in]  d     Trajectory analysis data structure.
 * \param[out] x     Topology data pointer to initialize.
 *   (can be NULL, in which case it is not used).
 * \param[out] box   Box size from the topology file
 *   (can be NULL, in which case it is not used).
 * \param[out] ePBC  The ePBC field from the topology
 *   (can be NULL, in which case it is not used).
 * \returns    0 on success, a non-zero error code on error.
 *
 * If \ref ANA_USE_TOPX has not been specified, the \p x parameter should be
 * NULL.
 *
 * The pointer returned in \p *x should not be freed.
 */
int
gmx_ana_get_topconf(gmx_ana_traj_t *d, rvec **x, matrix box, int *ePBC)
{
    if (box)
    {
        copy_mat(d->boxtop, box);
    }
    if (ePBC)
    {
        *ePBC = d->ePBC;
    }
    if (x)
    {
        if (!(d->flags & ANA_USE_TOPX))
        {
            gmx_incons("topology coordinates requested by ANA_USE_TOPX not set");
            *x = NULL;
            return EINVAL;
        }
        *x = d->xtop;
    }
    return 0;
}

/*! \brief
 * Loads default index groups from a selection file.
 *
 * \param[in,out] d     Trajectory analysis data structure.
 * \param[out]    grps  Pointer to receive the default groups.
 * \returns       0 on success, a non-zero error code on error.
 */
static int
init_default_selections(gmx_ana_traj_t *d, gmx_ana_indexgrps_t **grps)
{
    gmx_ana_selcollection_t  *sc;
    char                     *fnm;
    int                       nr, nr_notempty, i;
    int                       rc;

    /* If an index file is provided, just load it and exit. */
    if (d->ndxfile)
    {
        gmx_ana_indexgrps_init(grps, d->top, d->ndxfile);
        return 0;
    }
    /* Initialize groups to NULL if we return prematurely. */
    *grps = NULL;
    /* Return immediately if no topology provided. */
    if (!d->top)
    {
        return 0;
    }

    /* Find the default selection file, return if none found. */
    fnm = low_gmxlibfn("defselection.dat", TRUE, FALSE);
    if (fnm == NULL)
    {
        return 0;
    }

    /* Create a temporary selection collection. */
    rc = gmx_ana_selcollection_create(&sc, d->pcc);
    if (rc != 0)
    {
        sfree(fnm);
        return rc;
    }
    rc = gmx_ana_selmethod_register_defaults(sc);
    if (rc != 0)
    {
        gmx_ana_selcollection_free(sc);
        sfree(fnm);
        gmx_fatal(FARGS, "default selection method registration failed");
        return rc;
    }
    /* FIXME: It would be better to not have the strings here hard-coded. */
    gmx_ana_selcollection_set_refpostype(sc, "atom");
    gmx_ana_selcollection_set_outpostype(sc, "atom", FALSE);

    /* Parse and compile the file with no external groups. */
    rc = gmx_ana_selcollection_parse_file(sc, fnm, NULL);
    sfree(fnm);
    if (rc != 0)
    {
        gmx_ana_selcollection_free(sc);
        fprintf(stderr, "\nWARNING: default selection(s) could not be parsed\n");
        return rc;
    }
    gmx_ana_selcollection_set_topology(sc, d->top, -1);
    rc = gmx_ana_selcollection_compile(sc);
    if (rc != 0)
    {
        gmx_ana_selcollection_free(sc);
        fprintf(stderr, "\nWARNING: default selection(s) could not be compiled\n");
        return rc;
    }

    /* Count the non-empty groups and check that there are no dynamic
     * selections. */
    nr          = gmx_ana_selcollection_get_count(sc);
    nr_notempty = 0;
    for (i = 0; i < nr; ++i)
    {
        gmx_ana_selection_t  *sel;

        sel = gmx_ana_selcollection_get_selection(sc, i);
        if (sel->bDynamic)
        {
            fprintf(stderr, "\nWARNING: dynamic default selection ignored\n");
        }
        else if (sel->g->isize > 0)
        {
            ++nr_notempty;
        }
    }

    /* Copy the groups to the output structure */
    gmx_ana_indexgrps_alloc(grps, nr_notempty);
    nr_notempty = 0;
    for (i = 0; i < nr; ++i)
    {
        gmx_ana_selection_t  *sel;

        sel = gmx_ana_selcollection_get_selection(sc, i);
        if (!sel->bDynamic && sel->g->isize > 0)
        {
            gmx_ana_index_t  *g;

            g = gmx_ana_indexgrps_get_grp(*grps, nr_notempty);
            gmx_ana_index_copy(g, sel->g, TRUE);
            g->name = strdup(sel->name);
            ++nr_notempty;
        }
    }

    gmx_ana_selcollection_free(sc);
    return 0;
}

/*!
 * \param[in,out] d     Trajectory analysis data structure.
 * \returns       0 on success, a non-zero error code on error.
 *
 * Initializes the selection data in \c gmx_ana_traj_t based on
 * the selection options and/or index files provided on the command line.
 *
 * This function is called automatically by parse_trjana_args() and should
 * not be called directly unless \ref ANA_USER_SELINIT is specified.
 *
 * \see ANA_USER_SELINIT
 */
int
gmx_ana_init_selections(gmx_ana_traj_t *d)
{
    int                      rc;
    int                      i;
    int                      nr;
    gmx_ana_indexgrps_t     *grps;
    int                      natoms;
    gmx_bool                 bStdIn;
    gmx_bool                 bInteractive;
    gmx_bool                 bOk;

    if (d->sel)
    {
        gmx_call("init_selections called more than once\n"
                 "perhaps you forgot ANA_USER_SELINIT");
        return -1;
    }

    gmx_ana_selcollection_set_veloutput(d->sc,
                                        d->frflags & (TRX_READ_V | TRX_NEED_V));
    gmx_ana_selcollection_set_forceoutput(d->sc,
                                          d->frflags & (TRX_READ_F | TRX_NEED_F));
    /* Check if we need some information from the topology */
    if (gmx_ana_selcollection_requires_top(d->sc))
    {
        rc = load_topology(d, TRUE);
        if (rc != 0)
        {
            return rc;
        }
    }
    /* Initialize the default selection methods */
    rc = gmx_ana_selmethod_register_defaults(d->sc);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "default selection method registration failed");
        return rc;
    }
    /* Initialize index groups.
     * We ignore the return value to continue without the default groups if
     * there is an error there. */
    init_default_selections(d, &grps);
    /* Parse the selections */
    bStdIn = (d->selfile && d->selfile[0] == '-' && d->selfile[1] == 0)
        || (d->selection && d->selection[0] == 0)
        || (!d->selfile && !d->selection);
    /* Behavior is not very pretty if we cannot check for interactive input,
     * but at least it should compile and work in most cases. */
#ifdef HAVE_UNISTD_H
    bInteractive = bStdIn && isatty(fileno(stdin));
#else
    bInteractive = bStdIn;
#endif
    if (bStdIn && bInteractive)
    {
        /* Parse from stdin */
        /* First we parse the reference groups if there are any */
        if (d->nrefgrps > 0)
        {
            fprintf(stderr, "\nSpecify ");
            if (d->nrefgrps == 1)
            {
                fprintf(stderr, "a reference selection");
            }
            else
            {
                fprintf(stderr, "%d reference selections", d->nrefgrps);
            }
            fprintf(stderr, ":\n");
            fprintf(stderr, "(one selection per line, 'help' for help)\n");
            rc = gmx_ana_selcollection_parse_stdin(d->sc, d->nrefgrps, grps, TRUE);
            nr = gmx_ana_selcollection_get_count(d->sc);
            if (rc != 0 || nr != d->nrefgrps)
            {
                gmx_ana_traj_free(d);
                gmx_input("unrecoverable error in selection parsing");
                return rc;
            }
        }
        /* Then, we parse the analysis groups */
        fprintf(stderr, "\nSpecify ");
        if (d->nanagrps == 1)
        {
            fprintf(stderr, "a selection");
        }
        else if (d->nanagrps == -1)
        {
            fprintf(stderr, "any number of selections");
        }
        else
        {
            fprintf(stderr, "%d selections", d->nanagrps);
        }
        fprintf(stderr, " for analysis:\n");
        fprintf(stderr, "(one selection per line, 'help' for help%s)\n",
                d->nanagrps == -1 ? ", Ctrl-D to end" : "");
        rc = gmx_ana_selcollection_parse_stdin(d->sc, d->nanagrps, grps, TRUE);
        fprintf(stderr, "\n");
    }
    else if (bStdIn)
    {
        rc = gmx_ana_selcollection_parse_stdin(d->sc, -1, grps, FALSE);
    }
    else if (d->selection)
    {
        rc = gmx_ana_selcollection_parse_str(d->sc, d->selection, grps);
    }
    else
    {
        rc = gmx_ana_selcollection_parse_file(d->sc, d->selfile, grps);
    }
    if (grps)
    {
        gmx_ana_indexgrps_free(grps);
    }
    if (rc != 0)
    {
        /* Free memory for memory leak checking */
        gmx_ana_traj_free(d);
        gmx_input("selection(s) could not be parsed");
        return rc;
    }

    /* Check the number of groups */
    nr = gmx_ana_selcollection_get_count(d->sc);
    if (nr == 0)
    {
        /* TODO: Don't print this if the user has requested help */
        fprintf(stderr, "Nothing selected, finishing up.\n");
        gmx_ana_traj_free(d);
        /* TODO: It would be better to return some code that tells the caller
         * that one should exit. */
        exit(0);
    }
    if (nr <= d->nrefgrps)
    {
        gmx_input("selection does not specify enough index groups");
        return -1;
    }
    if (d->nanagrps <= 0)
    {
        d->nanagrps = nr - d->nrefgrps;
    }
    else if (nr != d->nrefgrps + d->nanagrps)
    {
        gmx_input("selection does not specify the correct number of index groups");
        return -1;
    }

    if (d->flags & ANA_DEBUG_SELECTION)
    {
        gmx_ana_selcollection_print_tree(stderr, d->sc, FALSE);
    }
    if (gmx_ana_selcollection_requires_top(d->sc))
    {
        rc = load_topology(d, TRUE);
        if (rc != 0)
        {
            return rc;
        }
    }
    if (d->top)
    {
        natoms = -1;
    }
    else
    {
        rc = init_first_frame(d);
        if (rc != 0)
        {
            return rc;
        }
        natoms = d->fr->natoms;
    }
    gmx_ana_selcollection_set_topology(d->sc, d->top, natoms);
    rc = gmx_ana_selcollection_compile(d->sc);
    if (rc != 0)
    {
        /* Free memory for memory leak checking */
        gmx_ana_traj_free(d);
        gmx_input("selection could not be compiled");
        return rc;
    }
    /* Create the selection array */
    d->ngrps = gmx_ana_selcollection_get_count(d->sc);
    if (!(d->flags & ANA_USE_FULLGRPS))
    {
        d->ngrps -= d->nrefgrps;
    }
    snew(d->sel, d->ngrps);
    for (i = 0; i < d->ngrps; ++i)
    {
        if (d->flags & ANA_USE_FULLGRPS)
        {
            d->sel[i] = gmx_ana_selcollection_get_selection(d->sc, i);
        }
        else
        {
            d->sel[i] = gmx_ana_selcollection_get_selection(d->sc, i + d->nrefgrps);
        }
    }
    if (d->flags & ANA_DEBUG_SELECTION)
    {
        fprintf(stderr, "\n");
        gmx_ana_selcollection_print_tree(stderr, d->sc, FALSE);
        fprintf(stderr, "\n");
        gmx_ana_poscalc_coll_print_tree(stderr, d->pcc);
        fprintf(stderr, "\n");
    }

    /* Initialize the position evaluation */
    gmx_ana_poscalc_init_eval(d->pcc);
    if (d->flags & ANA_DEBUG_SELECTION)
    {
        gmx_ana_poscalc_coll_print_tree(stderr, d->pcc);
        fprintf(stderr, "\n");
    }

    /* Check that dynamic selections are not provided if not allowed */
    if (d->flags & ANA_NO_DYNSEL)
    {
        for (i = 0; i < d->nrefgrps + d->nanagrps; ++i)
        {
            gmx_ana_selection_t *sel;

            sel = gmx_ana_selcollection_get_selection(d->sc, i);
            if (sel->bDynamic)
            {
                gmx_fatal(FARGS, "%s does not support dynamic selections",
                          ShortProgram());
                return -1;
            }
        }
    }
    /* Check that non-atom positions are not provided if not allowed.
     * TODO: It would be better to have these checks in the parser. */
    if (d->flags & ANA_ONLY_ATOMPOS)
    {
        for (i = 0; i < d->nanagrps; ++i)
        {
            gmx_ana_selection_t *sel;

            sel = gmx_ana_selcollection_get_selection(d->sc, i + d->nrefgrps);
            if (sel->p.m.type != INDEX_ATOM)
            {
                gmx_fatal(FARGS, "%s does not support non-atom positions",
                          ShortProgram());
                return -1;
            }
        }
    }
    /* Create the names array */
    snew(d->grpnames, d->ngrps);
    for (i = 0; i < d->ngrps; ++i)
    {
        d->grpnames[i] = gmx_ana_selection_name(d->sel[i]);
    }

    return 0;
}

/*!
 * \param[in,out] d       Trajectory analysis data structure.
 * \param[in]     type    Type of covered fractions to calculate.
 * \returns       0 on success, a non-zero error code on error.
 *
 * By default, covered fractions are not calculated.
 * If this function is called, the covered fraction calculation is
 * initialize to calculate the fractions of type \p type for each selection.
 * The function must be called after gmx_ana_init_selections() has been
 * called.
 *
 * For more fine-grained control of the calculation, you can use
 * gmx_ana_selection_init_coverfrac(): if you initialize some selections
 * this function to something else than CFRAC_NONE before calling
 * gmx_ana_init_coverfrac(), these settings are not overwritten.
 *
 * You only need to call this function if your program needs to have
 * information about the covered fractions, e.g., for normalization.
 *
 * \see gmx_ana_selection_init_coverfrac()
 */
int
gmx_ana_init_coverfrac(gmx_ana_traj_t *d, e_coverfrac_t type)
{
    int                 g;

    if (type == CFRAC_NONE)
    {
        return 0;
    }

    for (g = 0; g < d->ngrps; ++g)
    {
        if (d->sel[g]->cfractype == CFRAC_NONE)
        {
            gmx_ana_selection_init_coverfrac(d->sel[g], type);
        }
    }
    return 0;
}

/*!
 * \param[in] out  Output file.
 * \param[in] d    Trajectory analysis data structure.
 * \returns   0 on success, a non-zero error code on error.
 */
int xvgr_selections(FILE *out, gmx_ana_traj_t *d)
{
    xvgr_selcollection(out, d->sc, d->oenv);
    return 0;
}

/*!
 * \param[in,out] d       Trajectory analysis data structure.
 * \returns       0 on success, a non-zero error code on error.
 */
static int init_first_frame(gmx_ana_traj_t *d)
{
    int                 i;

    /* Return if we have already initialized the trajectory */
    if (d->fr)
    {
        return 0;
    }

    d->frflags |= TRX_NEED_X;

    snew(d->fr, 1);

    if (d->trjfile)
    {
        if (!read_first_frame(d->oenv, &d->status, d->trjfile, d->fr, d->frflags))
        {
            gmx_input("could not read coordinates from trajectory");
            return EIO;
        }

        if (d->top && d->fr->natoms > d->top->atoms.nr)
        {
            gmx_fatal(FARGS, "Trajectory (%d atoms) does not match topology (%d atoms)",
                      d->fr->natoms, d->top->atoms.nr);
            return -1;
        }
        /* check index groups */
        for (i = 0; i < d->ngrps; ++i)
        {
            gmx_ana_index_check(d->sel[i]->g, d->fr->natoms);
        }
    }
    else
    {
        /* Prepare a frame from topology information */
        /* TODO: Initialize more of the fields */
        if (d->frflags & (TRX_NEED_V))
        {
            gmx_impl("Velocity reading from a topology not implemented");
            return -1;
        }
        if (d->frflags & (TRX_NEED_F))
        {
            gmx_input("Forces cannot be read from a topology");
            return -1;
        }
        d->fr->flags  = d->frflags;
        d->fr->natoms = d->top->atoms.nr;
        d->fr->bX     = TRUE;
        snew(d->fr->x, d->fr->natoms);
        memcpy(d->fr->x, d->xtop, sizeof(*d->fr->x)*d->fr->natoms);
        d->fr->bBox   = TRUE;
        copy_mat(d->boxtop, d->fr->box);
    }

    set_trxframe_ePBC(d->fr, d->ePBC);

    return 0;
}

/*!
 * \param[in,out] d       Trajectory analysis data structure.
 * \param[out]    fr      First frame in the trajectory.
 * \returns       0 on success, a non-zero error code on error.
 *
 * The pointer stored in \p *fr should not be freed by the caller.
 *
 * You only need to call this function if your program needs to do some
 * initialization for which it requires data from the first frame.
 *
 * \see gmx_ana_do()
 */
int gmx_ana_get_first_frame(gmx_ana_traj_t *d, t_trxframe **fr)
{
    int rc;

    rc = init_first_frame(d);
    if (rc != 0)
    {
        *fr = NULL;
        return rc;
    }
    *fr = d->fr;
    return 0;
}

/*!
 * \param[in,out] d   Trajectory analysis data structure.
 * \param[in] flags   Combination of flags
 *      (currently, there are no flags defined).
 * \param[in] analyze Pointer to frame analysis function.
 * \param     data    User data to be passed to \p analyze.
 * \returns   0 on success, a non-zero error code on error.
 *
 * This function performs the actual analysis of the trajectory.
 * It loops through all the frames in the trajectory, and calls
 * \p analyze for each frame to perform the actual analysis.
 * Before calling \p analyze, selections are updated (if there are any).
 * See gmx_analysisfunc() for description of \p analyze parameters.
 *
 * This function also calculates the number of frames during the run.
 */
int gmx_ana_do(gmx_ana_traj_t *d, int flags, gmx_analysisfunc analyze, void *data)
{
    t_pbc               pbc;
    t_pbc              *ppbc;
    int                 rc;
    gmx_rmpbc_t         gpbc = NULL;

    rc = init_first_frame(d);
    if (rc != 0)
    {
        return rc;
    }

    ppbc = d->bPBC ? &pbc : 0;
    if (!d->top)
    {
        d->bRmPBC = FALSE;
    }
    if (d->bRmPBC)
    {
        gpbc = gmx_rmpbc_init(&d->top->idef, d->ePBC, d->fr->natoms, d->fr->box);
    }
    d->nframes = 0;
    do
    {
        if (d->bRmPBC)
        {
            gmx_rmpbc_trxfr(gpbc, d->fr);
        }
        if (ppbc)
        {
            set_pbc(&pbc, d->ePBC, d->fr->box);
        }

        gmx_ana_poscalc_init_frame(d->pcc);
        rc = gmx_ana_selcollection_evaluate(d->sc, d->fr, ppbc);
        if (rc != 0)
        {
            close_trj(d->status);
            gmx_fatal(FARGS, "selection evaluation failed");
            return rc;
        }
        rc = analyze(d->top, d->fr, ppbc, d->ngrps, d->sel, data);
        if (rc != 0)
        {
            close_trj(d->status);
            return rc;
        }

        d->nframes++;
    }
    while (d->trjfile && read_next_frame(d->oenv, d->status, d->fr));
    if (d->bRmPBC)
    {
        gmx_rmpbc_done(gpbc);
    }
    if (d->trjfile)
    {
        close_trj(d->status);
        fprintf(stderr, "Analyzed %d frames, last time %.3f\n",
                d->nframes, d->fr->time);
    }
    else
    {
        fprintf(stderr, "Analyzed topology coordinates\n");
    }

    /* Restore the maximal groups for dynamic selections */
    rc = gmx_ana_selcollection_evaluate_fin(d->sc, d->nframes);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "selection evaluation failed");
    }

    return rc;
}

/*!
 * \param[in,out] d       Trajectory analysis data structure.
 * \param[out]    nframes Number of frames.
 * \returns   0 on success, a non-zero error code on error.
 *
 * If called before gmx_ana_do(), the behavior is currently undefined.
 */
extern int
gmx_ana_get_nframes(gmx_ana_traj_t *d, int *nframes)
{
    *nframes = d->nframes;
    return 0;
}
