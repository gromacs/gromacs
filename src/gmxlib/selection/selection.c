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
/*! \page selengine Text-based selections
 *
 * \section selection_basics Basics
 *
 * Selections are enabled for an analysis program automatically.
 * By default, dynamic selections are allowed, and the user can freely
 * choose whether to analyze atoms or centers of mass/geometry of
 * residues/molecules.
 * These defaults, as well as some others, can be changed by specifying
 * flags for gmx_ana_traj_create().
 *
 * Two command-line options are always added:
 *  - \p -select can be used to specify the selection on the command line.
 *  - \p -sf can be used to specify a file name from which the selection is
 *    read.
 *    If both options are specified, \p -select takes precedence.
 *  - If neither of the above is present, the user is prompted to type
 *    selections.
 *    This is also done if an empty string is passed to \p -select.
 *    An index file can be provided with \p -n, and the user can also select
 *    groups from it instead of writing a full selection.
 *    If no index file is provided, the default groups are generated.
 * 
 * \todo
 * Improve the user interface for interactive selection input.
 *
 * If \ref ANA_ONLY_ATOMPOS is not specified, an additional command-line
 * option is added:
 *  - \p -seltype can be used to specify the default type of positions to
 *    calculate for each selection.
 *
 * An additional command-line option is added if dynamic selections are
 * allowed (i.e., if \ref ANA_NO_DYNSEL is not specified).
 * This options control how position keywords work within selections:
 *  - By default, each atom can be selected separately based on its
 *    coordinates (unless \ref ANA_REQUIRE_WHOLE is specified), but this
 *    can be changed with \p -selrpos to select all the atoms in a
 *    residue/molecule based on the center of mass/geometry of the
 *    residue/molecule.
 *
 * The two options that specify the position types can take any of the
 * following values:
 *  - \p atom uses individual atom positions.
 *  - \p res_com and \p res_cog use the centers of mass/geometry of residues.
 *  - \p mol_com and \p mol_cog use the centers of mass/geometry of molecules.
 *  - The way the COM/COG are calculated can be altered with a one of the
 *    following prefixes (e.g., \p part_res_com):
 *     - \p whole_ prefix calcualtes the centers for the whole
 *       residue/molecule, even if only part of it is selected.
 *     - \p part_ prefix calculates the centers for the selected atoms,
 *       but uses always the same atoms for the same residue/molecule.
 *       The used atoms are determined from the the largest group allowed by
 *       the selection.
 *     - \p dyn_ prefix calculates the centers strictly only for the
 *       selected atoms.
 *     .
 *  - If no prefix is specified, \p -seltype defaults to \p part_ and
 *    \p -selrpos defaults to \p whole_.
 *    The latter is often desirable to select the same molecules in different
 *    tools, while the first is a compromise between speed (\p dyn_ positions
 *    can be slower to evaluate than \p part_) and intuitive behavior.
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
 * The analysis program can define \ref ANA_REQUIRE_WHOLE to make the
 * default behavior appropriate for the most common uses where the groups
 * should consist of atoms within a single residue/molecule.
 *
 * \note
 * Currently, it is not possible to write selection expressions that would
 * evaluate to index groups where the same atom occurs more than once.
 * Such index groups may be useful with tools like \c g_angle to calculate
 * the averages over several angles.
 * It is possible to implement such a feature, but currently the best solution is to
 * write/modify the tool such that it can deal with multiple index groups,
 * each of which then defines angles that don't have overlapping atoms.
 *
 * 
 * \section selection_syntax Selection syntax
 * 
 * A set of selections consists of one or more selections, separated by
 * semicolons. Each selection defines a set of positions for the analysis.
 * Each selection can also be preceded by a string that gives a name for
 * the selection for use in, e.g., graph legends.
 * If no name is provided, a name like "Selection 3" is generated automatically.
 * It is also possible to use variables to store selection expressions.
 * A variable is defined with the following syntax:
\verbatim
VARNAME = EXPR ;
\endverbatim
 * After this, you can use VARNAME anywhere where EXPR would be valid.
 *
 * No semicolons should be used when providing the selections interactively;
 * in the interactive prompt, each selection should appear on a line of its
 * own. Lines can be continued with \\ if necessary.
 * Notice that the above only applies to real interactive input,
 * not if you provide the selections, e.g., from a pipe.
 *
 * Positions can be defined in several different ways (ATOM_EXPR stands for any
 * expression that evaluates to a set of atoms):
 *  - <tt>(X, Y, Z)</tt>: a fixed position (X, Y, Z should be real numbers).
 *  - <tt>cog of ATOM_EXPR [pbc]</tt>: center of geometry of ATOM_EXPR.
 *  - <tt>com of EXPR [pbc]</tt>: center of mass of ATOM_EXPR.
 *  - <tt>res_com of ATOM_EXPR</tt>: a set of positions consisting of centers
 *    of mass of the residues in ATOM_EXPR.
 *    <tt>res_com</tt> can be replaced with any of the position types
 *    acceptable for the \p -seltype and \p -selrpos options.
 *  - <tt>ATOM_EXPR</tt>: the default positions provided with \p -seltype are
 *    evaluated for ATOM_EXPR.
 *
 * If you specify \p pbc with the \p cog or \p com expressions, 
 * the center is calculated iteratively to try to deal with cases where
 * ATOM_EXPR wraps around periodic boundary conditions.
 * 
 * To select atoms (ATOM_EXPR expressions above), the following keywords are
 * currently available:
 *  - <tt>atomnr [INT|INT to INT] ... </tt>: selects atoms by number
 *    (first atom = 1)
 *  - <tt>resnr [INT|INT to INT] ... </tt>: selects residues by number
 *    (first residue = 1)
 *  - <tt>name STR ... </tt>: selects atoms by name
 *  - <tt>type STR ... </tt>: selects atoms by type
 *  - <tt>resname STR ... </tt>: selects residues by name
 *  - <tt>mass</tt>: returns the mass of the atom, use, e.g., like
 *    <tt>mass < 15</tt>
 *  - <tt>charge</tt>: returns the charge of the atom
 *  - <tt>occupancy</tt>: returns the occupancy of the atom from the PDB file
 *  - <tt>betafactor</tt>: returns the B-factor of the atom from the PDB file
 *  - <tt>within RADIUS of POS_EXPR</tt>: Selects atoms that are within
 *    RADIUS of any positions defined by POS_EXPR. POS_EXPR can be defined
 *    as in selections above (\p -selrpos is used instead of \p -seltype to
 *    convert atom expressions to positions).
 *  - <tt>insolidangle span POS_EXPR center POS [cutoff ANGLE]</tt>:
 *    Selects atoms that are within ANGLE degrees (default=5) of any position
 *    in POS_EXPR as seen from POS (a position expression that evaluates to
 *    a single position), i.e., atoms in the solid angle spanned by the
 *    positions in POS_EXPR and centered at POS
 *    (see \subpage sm_insolidangle "detailed explanation").
 *
 * For string-based selection keywords, it is possible to use wildcards
 * (e.g., <tt>name "C*"</tt>) or regular expressions
 * (e.g., <tt>resname "R[AB]"</tt>).
 * Strings that contain non-alphanumeric characters should be enclosed in
 * double quotes as in the examples.
 * The match type is automatically guessed from the string: if it contains any
 * other characters than letters, numbers, '*', or '?', it is interpreted
 * as a regular expression.
 *
 * In addition, the following keywords yield numeric values that can
 * be compared with each other or constant values to select a subset of
 * the particles (\p resnr can also be used in this way):
 *  - <tt>distance from POS [cutoff REAL]</tt>: calculates the distance from
 *    POS
 *  - <tt>mindistance from POS_EXPR [cutoff REAL]</tt>: calculates the minimum
 *    distance from the positions in POS_EXPR
 *  - <tt>x</tt>, <tt>y</tt>, <tt>z</tt>: returns the respective coordinate
 *    of the particle
 *
 * For the distance selection keywords, it is possible to specify a cutoff
 * to possibly speed up the evaluation: all distances above the specified
 * cutoff are returned as equal to the cutoff.
 *
 * For all selection keywords that select by position, the position specified
 * by \p -selrpos is used. This can be overridden by prepending a reference
 * position specifier to the keyword. For example,
 * <tt>res_com dist from POS</tt>
 * evaluates the residue center of mass distances irrespective of the value of
 * \p -selrpos.
 *
 * Arithmetic expressions are not implemented.
 *
 * It is also possible to combine expressions through boolean operations:
 *  - <tt>not ATOM_EXPR</tt>
 *  - <tt>ATOM_EXPR and ATOM_EXPR</tt>
 *  - <tt>ATOM_EXPR or ATOM_EXPR</tt>
 *
 * are all valid expressions.
 * Parenthesis can also be used to alter the evaluation order.
 *
 * All selections are by default evaluated such that the atom indices are
 * returned in ascending order. This can be changed by appending
 * <tt>permute P1 P2 ... Pn</tt> to an expression.
 * The \c Pi should form a permutation of the numbers 1 to n.
 * This keyword permutes each n-position block in the selection such that the
 * i'th position in the block becomes Pi'th.
 * Note that it is the positions that are permuted, not individual atoms.
 * A fatal error occurs if the size of the selection is not a multiple of n.
 * It is only possible to permute the whole selection expression, not any
 * subexpressions, i.e., the \c permute keyword should appear last in a
 * selection.
 *
 *
 * \section selection_eval Selection evaluation and optimization
 *
 * Boolean evaluation proceeds from left to right and is short-circuiting,
 * i.e., as soon as it is known whether an atom will be selected, the
 * remaining expressions are not evaluated at all.
 * This can be used to optimize the selections: you should write the
 * most restrictive and/or the most inexpensive expressions first in
 * boolean expressions.
 * The relative ordering between dynamic and static expressions does not
 * matter: all static expressions are evaluated only once, before the first
 * frame, and the result becomes the leftmost expression.
 *
 * Another point for optimization is in common subexpressions: they are not
 * automatically recognized, but can be manually optimized by the use of
 * variables. This can have a big impact on the performance of complex
 * selections, in particular if you define several index groups like this:
\verbatim
rdist = distance from com of resnr 1 to 5;
resname RES and rdist < 2;
resname RES and rdist < 4;
resname RES and rdist < 6;
\endverbatim
 * Without the variable assignment, the distances would be evaluated three
 * times, although they are exactly the same within each selection.
 * Anything assigned into a variable becomes a common subexpression that
 * is evaluated only once during a frame.
 * Currently, in some cases the use of variables can actually lead to a small
 * performance loss because of the checks necessary to determine for which
 * atoms the expression has already been evaluated, but this should not be
 * a major problem.
 *
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
/*! \file
 * \brief
 * Implementation of functions in selection.h.
 */
/*! \dir src/gmxlib/selection
 * \brief
 * Source code for selection-related routines.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <smalloc.h>
#include <statutil.h>
#include <string2.h>
#include <xvgr.h>

#include <poscalc.h>
#include <selection.h>
#include <selmethod.h>

#include "selcollection.h"
#include "selelem.h"
#include "symrec.h"

/*!
 * \param[out] scp Pointer to a newly allocated empty selection collection.
 * \param[in]  pcc Position calculation data structure to use for selection
 *   position evaluation.
 * \returns    0 on success.
 */
int
gmx_ana_selcollection_create(gmx_ana_selcollection_t **scp,
                             gmx_ana_poscalc_coll_t *pcc)
{
    gmx_ana_selcollection_t *sc;
    
    snew(sc, 1);
    sc->rpost     = NULL;
    sc->spost     = NULL;
    sc->selstr    = NULL;
    sc->root      = NULL;
    sc->nr        = 0;
    sc->sel       = NULL;
    sc->top       = NULL;
    gmx_ana_index_clear(&sc->gall);
    sc->pcc       = pcc;
    _gmx_sel_symtab_create(&sc->symtab);
    *scp = sc;
    return 0;
}

/*!
 * \param[in,out] sc Selection collection to free.
 *
 * The pointer \p sc is invalid after the call.
 */
void
gmx_ana_selcollection_free(gmx_ana_selcollection_t *sc)
{
    int  i;

    sfree(sc->selstr);
    _gmx_selelem_free_chain(sc->root);
    if (sc->sel)
    {
        for (i = 0; i < sc->nr; ++i)
        {
            gmx_ana_selection_free(sc->sel[i]);
        }
    }
    sfree(sc->sel);
    gmx_ana_index_deinit(&sc->gall);
    _gmx_selcollection_clear_symtab(sc);
    sfree(sc);
}

/*!
 * \param[in,out] sc Selection collection.
 */
void
_gmx_selcollection_clear_symtab(gmx_ana_selcollection_t *sc)
{
    if (sc->symtab)
    {
        _gmx_sel_symtab_free(sc->symtab);
        sc->symtab = NULL;
    }
}

/*!
 * \param[in,out] sc        Selection collection to modify.
 * \param[in]     type      Default selection reference position type
 *   (one of the strings acceptable for gmx_ana_poscalc_type_from_enum()).
 *
 * Should be called before calling gmx_ana_selcollection_requires_top() or
 * gmx_ana_selcollection_parse_*().
 */
void
gmx_ana_selcollection_set_refpostype(gmx_ana_selcollection_t *sc,
                                     const char *type)
{
    sc->rpost     = type;
}

/*!
 * \param[in,out] sc        Selection collection to modify.
 * \param[in]     type      Default selection output position type
 *   (one of the strings acceptable for gmx_ana_poslcalc_type_from_enum()).
 * \param[in]     bMaskOnly If TRUE, the output positions are initialized
 *   using \ref POS_MASKONLY.
 *
 * If \p type is NULL, the default type is not modified.
 * Should be called before calling gmx_ana_selcollection_requires_top() or
 * gmx_ana_selcollection_parse_*().
 */
void
gmx_ana_selcollection_set_outpostype(gmx_ana_selcollection_t *sc,
                                     const char *type, bool bMaskOnly)
{
    if (type)
    {
        sc->spost     = type;
    }
    sc->bMaskOnly = bMaskOnly;
}

/*!
 * \param[in,out] sc        Selection collection to set the topology for.
 * \param[in]     top       Topology data.
 * \param[in]     natoms    Number of atoms. If <=0, the number of atoms in the
 *   topology is used.
 * \returns       0 on success, EINVAL if \p top is NULL and \p natoms <= 0.
 *
 * The topology is also set for the position calculation collection
 * associated with \p sc.
 *
 * \p natoms determines the largest atom index that can be selected by the
 * selection: even if the topology contains more atoms, they will not be
 * selected.
 */
int
gmx_ana_selcollection_set_topology(gmx_ana_selcollection_t *sc, t_topology *top,
                                   int natoms)
{
    gmx_ana_poscalc_coll_set_topology(sc->pcc, top);
    sc->top = top;

    /* Get the number of atoms from the topology if it is not given */
    if (natoms <= 0)
    {
        if (!sc->top)
        {
            gmx_incons("selections need either the topology or the number of atoms");
            return EINVAL;
        }
        natoms = sc->top->atoms.nr;
    }
    gmx_ana_index_init_simple(&sc->gall, natoms, NULL);
    return 0;
}

/*!
 * \param[in]  sc  Selection collection to query.
 * \returns    Number of selections in \p sc.
 *
 * If gmx_ana_selcollection_parse_*() has not been called, returns 0.
 *
 * \see gmx_ana_selcollection_get_selection()
 */
int
gmx_ana_selcollection_get_count(gmx_ana_selcollection_t *sc)
{
    return sc->nr;
}

/*!
 * \param[in]  sc  Selection collection to query.
 * \param[in]  i   Number of the selection.
 * \returns    Pointer to the \p i'th selection in \p sc,
 *   or NULL if there is no such selection.
 *
 * \p i should be between 0 and the value returned by
 * gmx_ana_selcollection_get_count().
 * The returned pointer should not be freed.
 * If gmx_ana_selcollection_compile() has not been called, the returned
 * selection is not completely initialized (but the returned pointer will be
 * valid even after compilation, and will point to the initialized selection).
 *
 * \see gmx_ana_selcollection_get_count()
 */
gmx_ana_selection_t *
gmx_ana_selcollection_get_selection(gmx_ana_selcollection_t *sc, int i)
{
    if (i < 0 || i >= sc->nr || !sc->sel)
        return NULL;
    return sc->sel[i];
}

/*!
 * \param[in]  sc  Selection collection to query.
 * \returns    TRUE if any selection in \p sc requires topology information,
 *   FALSE otherwise.
 *
 * Before gmx_ana_selcollection_parse_*(), the return value is based just on
 * the position types set.
 * After gmx_ana_selcollection_parse_*(), the return value also takes into account the
 * selection keywords used.
 */
bool
gmx_ana_selcollection_requires_top(gmx_ana_selcollection_t *sc)
{
    t_selelem   *sel;
    e_poscalc_t  type;
    int          flags;
    int          rc;

    if (sc->rpost)
    {
        flags = 0;
        rc = gmx_ana_poscalc_type_from_enum(sc->rpost, &type, &flags);
        if (rc == 0 && type != POS_ATOM)
        {
            return TRUE;
        }
    }
    if (sc->spost)
    {
        flags = 0;
        rc = gmx_ana_poscalc_type_from_enum(sc->spost, &type, &flags);
        if (rc == 0 && type != POS_ATOM)
        {
            return TRUE;
        }
    }

    sel = sc->root;
    while (sel)
    {
        if (_gmx_selelem_requires_top(sel))
        {
            return TRUE;
        }
        sel = sel->next;
    }
    return FALSE;
}

/*!
 * \param[in] fp      File handle to receive the output.
 * \param[in] sc      Selection collection to print.
 * \param[in] bValues If TRUE, the evaluated values of selection elements
 *   are printed as well.
 */
void
gmx_ana_selcollection_print_tree(FILE *fp, gmx_ana_selcollection_t *sc, bool bValues)
{
    t_selelem *sel;

    sel = sc->root;
    while (sel)
    {
        _gmx_selelem_print_tree(fp, sel, bValues, 0);
        sel = sel->next;
    }
}

/*!
 * \param[in] sel  Selection to free.
 *
 * After the call, the pointer \p sel is invalid.
 */
void
gmx_ana_selection_free(gmx_ana_selection_t *sel)
{
    sfree(sel->name);
    gmx_ana_pos_deinit(&sel->p);
    if (sel->m != sel->orgm)
    {
        sfree(sel->m);
    }
    if (sel->q != sel->orgq)
    {
        sfree(sel->q);
    }
    sfree(sel->orgm);
    sfree(sel->orgq);
    sfree(sel);
}

/*!
 * \param[in] sel  Selection whose name is needed.
 * \returns   Pointer to the name of the selection.
 *
 * The return value should not be freed by the caller.
 */
char *
gmx_ana_selection_name(gmx_ana_selection_t *sel)
{
    return sel->name;
}

/*!
 * \param[in] sel  Selection for which information should be printed.
 */
void
gmx_ana_selection_print_info(gmx_ana_selection_t *sel)
{
    fprintf(stderr, "\"%s\" (%d position%s, %d atom%s%s)", sel->name,
            sel->p.nr,     sel->p.nr     == 1 ? "" : "s",
            sel->g->isize, sel->g->isize == 1 ? "" : "s",
            sel->bDynamic ? ", dynamic" : "");
    fprintf(stderr, "\n");
}

/*!
 * \param[in] sel  Selection to initialize.
 * \param[in] type Type of covered fraction required.
 * \returns   TRUE if the covered fraction can be calculated for the selection,
 *   FALSE otherwise.    
 */
bool
gmx_ana_selection_init_coverfrac(gmx_ana_selection_t *sel, e_coverfrac_t type)
{
    sel->cfractype = type;
    if (type == CFRAC_NONE || !sel->selelem)
    {
        sel->bCFracDyn = FALSE;
    }
    else if (!_gmx_selelem_can_estimate_cover(sel->selelem))
    {
        sel->cfractype = CFRAC_NONE;
        sel->bCFracDyn = FALSE;
    }
    else
    {
        sel->bCFracDyn = TRUE;
    }
    sel->cfrac     = sel->bCFracDyn ? 0.0 : 1.0;
    sel->avecfrac  = sel->cfrac;
    return type == CFRAC_NONE || sel->cfractype != CFRAC_NONE;
}

/*!
 * \param[in] out  Output file.
 * \param[in] sc   Selection collection which should be written.
 */
void
xvgr_selcollection(FILE *out, gmx_ana_selcollection_t *sc)
{
    char  *buf;
    char  *p, *nl;

    if (bPrintXvgrCodes() && sc && sc->selstr)
    {
        fprintf(out, "# Selections:\n");
        buf = strdup(sc->selstr);
        p = buf;
        while (p && p[0] != 0)
        {
            nl = strchr(p, '\n');
            if (nl)
            {
                *nl = 0;
            }
            fprintf(out, "#   %s\n", p);
            p  = nl;
            if (nl)
            {
                ++p;
            }
        }
        fprintf(out, "#\n");
        sfree(buf);
    }
}
