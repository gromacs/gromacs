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
/*! \internal \file
 * \brief
 * Implementation of functions in selection.h.
 */
/*! \internal \dir src/gmxlib/selection
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
#include <gmx_fatal.h>

#include <poscalc.h>
#include <selection.h>
#include <selmethod.h>

#include "mempool.h"
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
                             gmx_ana_poscalc_coll_t   *pcc)
{
    gmx_ana_selcollection_t *sc;

    snew(sc, 1);
    sc->rpost         = NULL;
    sc->spost         = NULL;
    sc->bMaskOnly     = FALSE;
    sc->bVelocities   = FALSE;
    sc->bForces       = FALSE;
    sc->bDebugCompile = FALSE;
    sc->root          = NULL;
    sc->nr            = 0;
    sc->sel           = NULL;
    sc->nvars         = 0;
    sc->varstrs       = NULL;
    sc->top           = NULL;
    gmx_ana_index_clear(&sc->gall);
    sc->pcc       = pcc;
    sc->mempool   = NULL;
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

    _gmx_selelem_free_chain(sc->root);
    if (sc->sel)
    {
        for (i = 0; i < sc->nr; ++i)
        {
            gmx_ana_selection_free(sc->sel[i]);
        }
    }
    sfree(sc->sel);
    for (i = 0; i < sc->nvars; ++i)
    {
        sfree(sc->varstrs[i]);
    }
    sfree(sc->varstrs);
    gmx_ana_index_deinit(&sc->gall);
    if (sc->mempool)
    {
        _gmx_sel_mempool_destroy(sc->mempool);
    }
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
                                     const char              *type)
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
                                     const char *type, gmx_bool bMaskOnly)
{
    if (type)
    {
        sc->spost     = type;
    }
    sc->bMaskOnly = bMaskOnly;
}

/*!
 * \param[in,out] sc        Selection collection to modify.
 * \param[in]     bVelOut   If TRUE, selections will also evaluate
 *      velocities.
 */
void
gmx_ana_selcollection_set_veloutput(gmx_ana_selcollection_t *sc,
                                    gmx_bool                 bVelOut)
{
    sc->bVelocities = bVelOut;
}

/*!
 * \param[in,out] sc        Selection collection to modify.
 * \param[in]     bForceOut If TRUE, selections will also evaluate
 *      forces.
 */
void
gmx_ana_selcollection_set_forceoutput(gmx_ana_selcollection_t *sc,
                                      gmx_bool                 bForceOut)
{
    sc->bForces = bForceOut;
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
    {
        return NULL;
    }
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
gmx_bool
gmx_ana_selcollection_requires_top(gmx_ana_selcollection_t *sc)
{
    t_selelem   *sel;
    e_poscalc_t  type;
    int          flags;
    int          rc;

    if (sc->rpost)
    {
        flags = 0;
        rc    = gmx_ana_poscalc_type_from_enum(sc->rpost, &type, &flags);
        if (rc == 0 && type != POS_ATOM)
        {
            return TRUE;
        }
    }
    if (sc->spost)
    {
        flags = 0;
        rc    = gmx_ana_poscalc_type_from_enum(sc->spost, &type, &flags);
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
gmx_ana_selcollection_print_tree(FILE *fp, gmx_ana_selcollection_t *sc, gmx_bool bValues)
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
    sfree(sel->selstr);
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
gmx_bool
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
 * \param[in] oenv Output options structure.
 */
void xvgr_selcollection(FILE *out, gmx_ana_selcollection_t *sc,
                        const output_env_t oenv)
{
    int  i;

    if (output_env_get_xvg_format(oenv) != exvgNONE && sc)
    {
        fprintf(out, "# Selections:\n");
        for (i = 0; i < sc->nvars; ++i)
        {
            fprintf(out, "#   %s\n", sc->varstrs[i]);
        }
        for (i = 0; i < sc->nr; ++i)
        {
            fprintf(out, "#   %s\n", sc->sel[i]->selstr);
        }
        fprintf(out, "#\n");
    }
}
