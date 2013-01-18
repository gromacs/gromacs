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
/*! \internal
 * \page poscalcengine Position calculation engine
 *
 * The header file \ref poscalc.h defines an API for calculating positions
 * in an automated way. This is useful mostly in the selection engine, in
 * particular with dynamic selections, because the same COM/COG positions
 * may be needed in several contexts. The API makes it possible to
 * optimize the evaluation such that any heavy calculation is only done once,
 * and the results just copied if needed more than once.
 * The functions also provide a convenient interface for keeping the whole
 * \c gmx_ana_pos_t structure up-to-date.
 *
 * A new collection of position calculations is allocated with
 * gmx_ana_poscalc_coll_create().
 * Calculations within one collection should share the same topology, and
 * they are optimized. Calculations in different collections do not interact.
 * The topology for a collection can be set with
 * gmx_ana_poscalc_coll_set_topology().
 * This needs to be done before calling gmx_ana_poscalc_set_maxindex() for
 * any calculation in the collection, unless that calculation does not
 * require topology information.
 * All memory allocated for a collection and the calculations in it can be
 * freed with gmx_ana_poscalc_coll_free().
 *
 * A new calculation is created with gmx_ana_poscalc_create().
 * If flags need to be adjusted later, gmx_ana_poscalc_set_flags() can be
 * used.
 * After the flags are final, the largest possible index group for which the
 * positions are needed has to be set with gmx_ana_poscalc_set_maxindex().
 * gmx_ana_poscalc_coll_set_topology() should have been called before this
 * function is called.
 * After the above calls, gmx_ana_poscalc_init_pos() can be used to initialize
 * output to a \c gmx_ana_pos_t structure. Several different structures can be
 * initialized for the same calculation; the only requirement is that the
 * structure passed later to gmx_ana_poscalc_update() has been initialized
 * properly.
 * The memory allocated for a calculation can be freed with
 * gmx_ana_poscalc_free().
 *
 * The position evaluation is simple: gmx_ana_poscalc_init_frame() should be
 * called once for each frame, and gmx_ana_poscalc_update() can then be called
 * for each calculation that is needed for that frame.
 *
 * It is also possible to initialize the calculations based on a type provided
 * as a string.
 * The possible strings are returned by gmx_ana_poscalc_create_type_enum(),
 * and the string can be converted to the parameters for
 * gmx_ana_poscalc_create() using gmx_ana_poscalc_type_from_enum().
 * gmx_ana_poscalc_create_enum() is also provided for convenience.
 */
/*! \internal \file
 * \brief Implementation of functions in poscalc.h.
 *
 * \todo
 * There is probably some room for optimization in the calculation of
 * positions with bases.
 * In particular, the current implementation may do a lot of unnecessary
 * copying.
 * The interface would need to be changed to make it possible to use the
 * same output positions for several calculations.
 *
 * \todo
 * The current algorithm for setting up base calculations could be improved
 * in cases when there are calculations that cannot use a common base but
 * still overlap partially (e.g., with three calculations A, B, and C
 * such that A could use both B and C as a base, but B and C cannot use the
 * same base).
 * Setting up the bases in an optimal manner in every possible situation can
 * be quite difficult unless several bases are allowed for one calculation,
 * but better heuristics could probably be implemented.
 * For best results, the setup should probably be postponed (at least
 * partially) to gmx_ana_poscalc_init_eval().
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include <macros.h>
#include <smalloc.h>
#include <typedefs.h>
#include <pbc.h>
#include <vec.h>

#include <centerofmass.h>
#include <indexutil.h>
#include <poscalc.h>
#include <position.h>

/*! \internal \brief
 * Collection of \c gmx_ana_poscalc_t structures for the same topology.
 *
 * Calculations within the same structure are optimized to eliminate duplicate
 * calculations.
 */
struct gmx_ana_poscalc_coll_t
{
    /*! \brief
     * Topology data.
     *
     * Can be NULL if none of the calculations require topology data or if
     * gmx_ana_poscalc_coll_set_topology() has not been called.
     */
    t_topology                   *top;
    /** Pointer to the first data structure. */
    gmx_ana_poscalc_t            *first;
    /** Pointer to the last data structure. */
    gmx_ana_poscalc_t            *last;
    /** Whether the collection has been initialized for evaluation. */
    gmx_bool                      bInit;
};

/*! \internal \brief
 * Data structure for position calculation.
 */
struct gmx_ana_poscalc_t
{
    /*! \brief
     * Type of calculation.
     *
     * This field may differ from the type requested by the user, because
     * it is changed internally to the most effective calculation.
     * For example, if the user requests a COM calculation for residues
     * consisting of single atoms, it is simply set to POS_ATOM.
     * To provide a consistent interface to the user, the field \p itype
     * should be used when information should be given out.
     */
    e_poscalc_t               type;
    /*! \brief
     * Flags for calculation options.
     *
     * See \ref poscalc_flags "documentation of the flags".
     */
    int                       flags;

    /*! \brief
     * Type for the created indices.
     *
     * This field always agrees with the type that the user requested, but
     * may differ from \p type.
     */
    e_index_t                 itype;
    /*! \brief
     * Block data for the calculation.
     */
    t_blocka                  b;
    /*! \brief
     * Mapping from the blocks to the blocks of \p sbase.
     *
     * If \p sbase is NULL, this field is also.
     */
    int                      *baseid;
    /*! \brief
     * Maximum evaluation group.
     */
    gmx_ana_index_t           gmax;

    /** Position storage for calculations that are used as a base. */
    gmx_ana_pos_t            *p;

    /** TRUE if the positions have been evaluated for the current frame. */
    gmx_bool                      bEval;
    /*! \brief
     * Base position data for this calculation.
     *
     * If not NULL, the centers required by this calculation have
     * already been calculated in \p sbase.
     * The structure pointed by \p sbase is always a static calculation.
     */
    struct gmx_ana_poscalc_t *sbase;
    /** Next structure in the linked list of calculations. */
    struct gmx_ana_poscalc_t *next;
    /** Previous structure in the linked list of calculations. */
    struct gmx_ana_poscalc_t *prev;
    /** Number of references to this structure. */
    int                       refcount;
    /** Collection this calculation belongs to. */
    gmx_ana_poscalc_coll_t   *coll;
};

static const char *const poscalc_enum_strings[] = {
    "atom",
    "res_com",       "res_cog",
    "mol_com",       "mol_cog",
    "whole_res_com", "whole_res_cog",
    "whole_mol_com", "whole_mol_cog",
    "part_res_com",  "part_res_cog",
    "part_mol_com",  "part_mol_cog",
    "dyn_res_com",   "dyn_res_cog",
    "dyn_mol_com",   "dyn_mol_cog",
    NULL,
};
#define NENUM asize(poscalc_enum_strings)

/*! \brief
 * Returns the partition type for a given position type.
 *
 * \param [in] type  \c e_poscalc_t value to convert.
 * \returns    Corresponding \c e_indet_t.
 */
static e_index_t
index_type_for_poscalc(e_poscalc_t type)
{
    switch (type)
    {
        case POS_ATOM:    return INDEX_ATOM;
        case POS_RES:     return INDEX_RES;
        case POS_MOL:     return INDEX_MOL;
        case POS_ALL:     return INDEX_ALL;
        case POS_ALL_PBC: return INDEX_ALL;
    }
    return INDEX_UNKNOWN;
}

/*!
 * \param[in]     post  String (typically an enum command-line argument).
 *   Allowed values: 'atom', 'res_com', 'res_cog', 'mol_com', 'mol_cog',
 *   or one of the last four prepended by 'whole_', 'part_', or 'dyn_'.
 * \param[out]    type  \c e_poscalc_t corresponding to \p post.
 * \param[in,out] flags Flags corresponding to \p post.
 *   On input, the flags should contain the default flags.
 *   On exit, the flags \ref POS_MASS, \ref POS_COMPLMAX and
 *   \ref POS_COMPLWHOLE have been set according to \p post
 *   (the completion flags are left at the default values if no completion
 *   prefix is given).
 * \returns       0 if \p post is one of the valid strings, EINVAL otherwise.
 *
 * \attention
 * Checking is not complete, and other values than those listed above
 * may be accepted for \p post, but the results are undefined.
 */
int
gmx_ana_poscalc_type_from_enum(const char *post, e_poscalc_t *type, int *flags)
{
    const char *ptr;

    if (post[0] == 'a')
    {
        *type   = POS_ATOM;
        *flags &= ~(POS_MASS | POS_COMPLMAX | POS_COMPLWHOLE);
        return 0;
    }

    /* Process the prefix */
    ptr = post;
    if (post[0] == 'w')
    {
        *flags &= ~POS_COMPLMAX;
        *flags |= POS_COMPLWHOLE;
        ptr     = post + 6;
    }
    else if (post[0] == 'p')
    {
        *flags &= ~POS_COMPLWHOLE;
        *flags |= POS_COMPLMAX;
        ptr     = post + 5;
    }
    else if (post[0] == 'd')
    {
        *flags &= ~(POS_COMPLMAX | POS_COMPLWHOLE);
        ptr     = post + 4;
    }

    if (ptr[0] == 'r')
    {
        *type = POS_RES;
    }
    else if (ptr[0] == 'm')
    {
        *type = POS_MOL;
    }
    else
    {
        gmx_incons("unknown position calculation type");
        return EINVAL;
    }
    if (ptr[6] == 'm')
    {
        *flags |= POS_MASS;
    }
    else if (ptr[6] == 'g')
    {
        *flags &= ~POS_MASS;
    }
    else
    {
        gmx_incons("unknown position calculation type");
        return EINVAL;
    }
    return 0;
}

/*!
 * \param[in]  bAtom    If TRUE, the "atom" value is included.
 * \returns    NULL-terminated array of strings that contains the string
 *   values acceptable for gmx_ana_poscalc_type_from_enum().
 *
 * The first string in the returned list is always NULL to allow the list to
 * be used with Gromacs command-line parsing.
 */
const char **
gmx_ana_poscalc_create_type_enum(gmx_bool bAtom)
{
    const char **pcenum;
    size_t       i;

    if (bAtom)
    {
        snew(pcenum, NENUM+1);
        for (i = 0; i < NENUM; ++i)
        {
            pcenum[i+1] = poscalc_enum_strings[i];
        }
    }
    else
    {
        snew(pcenum, NENUM+1-1);
        for (i = 1; i < NENUM; ++i)
        {
            pcenum[i] = poscalc_enum_strings[i];
        }
    }
    pcenum[0] = NULL;
    return pcenum;
}

/*!
 * \param[out] pccp   Allocated position calculation collection.
 * \returns    0 for success.
 */
int
gmx_ana_poscalc_coll_create(gmx_ana_poscalc_coll_t **pccp)
{
    gmx_ana_poscalc_coll_t *pcc;

    snew(pcc, 1);
    pcc->top   = NULL;
    pcc->first = NULL;
    pcc->last  = NULL;
    pcc->bInit = FALSE;
    *pccp      = pcc;
    return 0;
}

/*!
 * \param[in,out] pcc   Position calculation collection data structure.
 * \param[in]     top   Topology data structure.
 *
 * This function should be called to set the topology before using
 * gmx_ana_poscalc_set_maxindex() for any calculation that requires
 * topology information.
 */
void
gmx_ana_poscalc_coll_set_topology(gmx_ana_poscalc_coll_t *pcc, t_topology *top)
{
    pcc->top = top;
}

/*!
 * \param[in] pcc   Position calculation collection to free.
 *
 * The pointer \p pcc is invalid after the call.
 * Any calculations in the collection are also freed, no matter how many
 * references to them are left.
 */
void
gmx_ana_poscalc_coll_free(gmx_ana_poscalc_coll_t *pcc)
{
    while (pcc->first)
    {
        gmx_ana_poscalc_free(pcc->first);
    }
    sfree(pcc);
}

/*!
 * \param[in] fp    File handle to receive the output.
 * \param[in] pcc   Position calculation collection to print.
 *
 * The output is very technical, making this function mainly useful for
 * debugging purposes.
 */
void
gmx_ana_poscalc_coll_print_tree(FILE *fp, gmx_ana_poscalc_coll_t *pcc)
{
    gmx_ana_poscalc_t *pc;
    int                i, j;

    fprintf(fp, "Position calculations:\n");
    i  = 1;
    pc = pcc->first;
    while (pc)
    {
        fprintf(fp, "%2d ", i);
        switch (pc->type)
        {
            case POS_ATOM:    fprintf(fp, "ATOM");    break;
            case POS_RES:     fprintf(fp, "RES");     break;
            case POS_MOL:     fprintf(fp, "MOL");     break;
            case POS_ALL:     fprintf(fp, "ALL");     break;
            case POS_ALL_PBC: fprintf(fp, "ALL_PBC"); break;
        }
        if (pc->itype != index_type_for_poscalc(pc->type))
        {
            fprintf(fp, " (");
            switch (pc->itype)
            {
                case INDEX_UNKNOWN: fprintf(fp, "???");  break;
                case INDEX_ATOM:    fprintf(fp, "ATOM"); break;
                case INDEX_RES:     fprintf(fp, "RES");  break;
                case INDEX_MOL:     fprintf(fp, "MOL");  break;
                case INDEX_ALL:     fprintf(fp, "ALL");  break;
            }
            fprintf(fp, ")");
        }
        fprintf(fp, " flg=");
        if (pc->flags & POS_MASS)
        {
            fprintf(fp, "M");
        }
        if (pc->flags & POS_DYNAMIC)
        {
            fprintf(fp, "D");
        }
        if (pc->flags & POS_MASKONLY)
        {
            fprintf(fp, "A");
        }
        if (pc->flags & POS_COMPLMAX)
        {
            fprintf(fp, "Cm");
        }
        if (pc->flags & POS_COMPLWHOLE)
        {
            fprintf(fp, "Cw");
        }
        if (!pc->flags)
        {
            fprintf(fp, "0");
        }
        fprintf(fp, " nr=%d nra=%d", pc->b.nr, pc->b.nra);
        fprintf(fp, " refc=%d", pc->refcount);
        fprintf(fp, "\n");
        if (pc->gmax.nalloc_index > 0)
        {
            fprintf(fp, "   Group: ");
            if (pc->gmax.isize > 20)
            {
                fprintf(fp, " %d atoms", pc->gmax.isize);
            }
            else
            {
                for (j = 0; j < pc->gmax.isize; ++j)
                {
                    fprintf(fp, " %d", pc->gmax.index[j] + 1);
                }
            }
            fprintf(fp, "\n");
        }
        if (pc->b.nalloc_a > 0)
        {
            fprintf(fp, "   Atoms: ");
            if (pc->b.nra > 20)
            {
                fprintf(fp, " %d atoms", pc->b.nra);
            }
            else
            {
                for (j = 0; j < pc->b.nra; ++j)
                {
                    fprintf(fp, " %d", pc->b.a[j] + 1);
                }
            }
            fprintf(fp, "\n");
        }
        if (pc->b.nalloc_index > 0)
        {
            fprintf(fp, "   Blocks:");
            if (pc->b.nr > 20)
            {
                fprintf(fp, " %d pcs", pc->b.nr);
            }
            else
            {
                for (j = 0; j <= pc->b.nr; ++j)
                {
                    fprintf(fp, " %d", pc->b.index[j]);
                }
            }
            fprintf(fp, "\n");
        }
        if (pc->sbase)
        {
            gmx_ana_poscalc_t *base;

            fprintf(fp, "   Base: ");
            j    = 1;
            base = pcc->first;
            while (base && base != pc->sbase)
            {
                ++j;
                base = base->next;
            }
            fprintf(fp, "%d", j);
            if (pc->baseid && pc->b.nr <= 20)
            {
                fprintf(fp, " id:");
                for (j = 0; j < pc->b.nr; ++j)
                {
                    fprintf(fp, " %d", pc->baseid[j]+1);
                }
            }
            fprintf(fp, "\n");
        }
        ++i;
        pc = pc->next;
    }
}

/*! \brief
 * Inserts a position calculation structure into its collection.
 *
 * \param pc     Data structure to insert.
 * \param before Data structure before which to insert
 *   (NULL = insert at end).
 *
 * Inserts \p pc to its collection before \p before.
 * If \p before is NULL, \p pc is appended to the list.
 */
static void
insert_poscalc(gmx_ana_poscalc_t *pc, gmx_ana_poscalc_t *before)
{
    if (before == NULL)
    {
        pc->next = NULL;
        pc->prev = pc->coll->last;
        if (pc->coll->last)
        {
            pc->coll->last->next = pc;
        }
        pc->coll->last = pc;
    }
    else
    {
        pc->prev     = before->prev;
        pc->next     = before;
        if (before->prev)
        {
            before->prev->next = pc;
        }
        before->prev = pc;
    }
    if (!pc->prev)
    {
        pc->coll->first = pc;
    }
}

/*! \brief
 * Removes a position calculation structure from its collection.
 *
 * \param pc    Data structure to remove.
 *
 * Removes \p pc from its collection.
 */
static void
remove_poscalc(gmx_ana_poscalc_t *pc)
{
    if (pc->prev)
    {
        pc->prev->next = pc->next;
    }
    else if (pc == pc->coll->first)
    {
        pc->coll->first = pc->next;
    }
    if (pc->next)
    {
        pc->next->prev = pc->prev;
    }
    else if (pc == pc->coll->last)
    {
        pc->coll->last = pc->prev;
    }
    pc->prev = pc->next = NULL;
}

/*! \brief
 * Initializes position calculation using the maximum possible input index.
 *
 * \param[in,out] pc  Position calculation data structure.
 * \param[in]     g   Maximum index group for the calculation.
 * \param[in]     bBase Whether \p pc will be used as a base or not.
 *
 * \p bBase affects on how the \p pc->gmax field is initialized.
 */
static void
set_poscalc_maxindex(gmx_ana_poscalc_t *pc, gmx_ana_index_t *g, gmx_bool bBase)
{
    gmx_ana_index_make_block(&pc->b, pc->coll->top, g, pc->itype, pc->flags & POS_COMPLWHOLE);
    /* Set the type to POS_ATOM if the calculation in fact is such. */
    if (pc->b.nr == pc->b.nra)
    {
        pc->type   = POS_ATOM;
        pc->flags &= ~(POS_MASS | POS_COMPLMAX | POS_COMPLWHOLE);
    }
    /* Set the POS_COMPLWHOLE flag if the calculation in fact always uses
     * complete residues and molecules. */
    if (!(pc->flags & POS_COMPLWHOLE)
        && (!(pc->flags & POS_DYNAMIC) || (pc->flags & POS_COMPLMAX))
        && (pc->type == POS_RES || pc->type == POS_MOL)
        && gmx_ana_index_has_complete_elems(g, pc->itype, pc->coll->top))
    {
        pc->flags &= ~POS_COMPLMAX;
        pc->flags |= POS_COMPLWHOLE;
    }
    /* Setup the gmax field */
    if ((pc->flags & POS_COMPLWHOLE) && !bBase && pc->b.nra > g->isize)
    {
        gmx_ana_index_copy(&pc->gmax, g, TRUE);
        sfree(pc->gmax.name);
        pc->gmax.name  = NULL;
    }
    else
    {
        gmx_ana_index_set(&pc->gmax, pc->b.nra, pc->b.a, NULL, 0);
    }
}

/*! \brief
 * Checks whether a position calculation should use a base at all.
 *
 * \param[in] pc   Position calculation data to check.
 * \returns   TRUE if \p pc can use a base and gets some benefit out of it,
 *   FALSE otherwise.
 */
static gmx_bool
can_use_base(gmx_ana_poscalc_t *pc)
{
    /* For atoms, it should be faster to do a simple copy, so don't use a
     * base. */
    if (pc->type == POS_ATOM)
    {
        return FALSE;
    }
    /* For dynamic selections that do not use completion, it is not possible
     * to use a base. */
    if ((pc->type == POS_RES || pc->type == POS_MOL)
        && (pc->flags & POS_DYNAMIC) && !(pc->flags & (POS_COMPLMAX | POS_COMPLWHOLE)))
    {
        return FALSE;
    }
    /* Dynamic calculations for a single position cannot use a base. */
    if ((pc->type == POS_ALL || pc->type == POS_ALL_PBC)
        && (pc->flags & POS_DYNAMIC))
    {
        return FALSE;
    }
    return TRUE;
}

/*! \brief
 * Checks whether two position calculations should use a common base.
 *
 * \param[in]     pc1 Calculation 1 to check for.
 * \param[in]     pc2 Calculation 2 to check for.
 * \param[in]     g1  Index group structure that contains the atoms from
 *   \p pc1.
 * \param[in,out] g   Working space, should have enough allocated memory to
 *   contain the intersection of the atoms in \p pc1 and \p pc2.
 * \returns   TRUE if the two calculations should be merged to use a common
 *   base, FALSE otherwise.
 */
static gmx_bool
should_merge(gmx_ana_poscalc_t *pc1, gmx_ana_poscalc_t *pc2,
             gmx_ana_index_t *g1, gmx_ana_index_t *g)
{
    gmx_ana_index_t  g2;

    /* Do not merge calculations with different mass weighting. */
    if ((pc1->flags & POS_MASS) != (pc2->flags & POS_MASS))
    {
        return FALSE;
    }
    /* Avoid messing up complete calculations. */
    if ((pc1->flags & POS_COMPLWHOLE) != (pc2->flags & POS_COMPLWHOLE))
    {
        return FALSE;
    }
    /* Find the overlap between the calculations. */
    gmx_ana_index_set(&g2, pc2->b.nra, pc2->b.a, NULL, 0);
    gmx_ana_index_intersection(g, g1, &g2);
    /* Do not merge if there is no overlap. */
    if (g->isize == 0)
    {
        return FALSE;
    }
    /* Full completion calculations always match if the type is correct. */
    if ((pc1->flags & POS_COMPLWHOLE) && (pc2->flags & POS_COMPLWHOLE)
        && pc1->type == pc2->type)
    {
        return TRUE;
    }
    /* The calculations also match if the intersection consists of full
     * blocks. */
    if (gmx_ana_index_has_full_ablocks(g, &pc1->b)
        && gmx_ana_index_has_full_ablocks(g, &pc2->b))
    {
        return TRUE;
    }
    return FALSE;
}

/*! \brief
 * Creates a static base for position calculation.
 *
 * \param     pc  Data structure to copy.
 * \returns   Pointer to a newly allocated base for \p pc.
 *
 * Creates and returns a deep copy of \p pc, but clears the
 * \ref POS_DYNAMIC and \ref POS_MASKONLY flags.
 * The newly created structure is set as the base (\c gmx_ana_poscalc_t::sbase)
 * of \p pc and inserted in the collection before \p pc.
 */
static gmx_ana_poscalc_t *
create_simple_base(gmx_ana_poscalc_t *pc)
{
    gmx_ana_poscalc_t *base;
    int                flags;
    int                rc;

    flags = pc->flags & ~(POS_DYNAMIC | POS_MASKONLY);
    rc    = gmx_ana_poscalc_create(&base, pc->coll, pc->type, flags);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "position calculation base creation failed");
    }
    set_poscalc_maxindex(base, &pc->gmax, TRUE);

    snew(base->p, 1);

    pc->sbase = base;
    remove_poscalc(base);
    insert_poscalc(base, pc);

    return base;
}

/*! \brief
 * Merges a calculation into another calculation such that the new calculation
 * can be used as a base.
 *
 * \param[in,out] base Base calculation to merge to.
 * \param[in,out] pc   Position calculation to merge to \p base.
 *
 * After the call, \p base can be used as a base for \p pc (or any calculation
 * that used it as a base).
 * It is assumed that any overlap between \p base and \p pc is in complete
 * blocks, i.e., that the merge is possible.
 */
static void
merge_to_base(gmx_ana_poscalc_t *base, gmx_ana_poscalc_t *pc)
{
    gmx_ana_index_t  gp, gb, g;
    int              isize, bnr;
    int              i, j, bi, bj, bo;

    base->flags |= pc->flags & (POS_VELOCITIES | POS_FORCES);
    gmx_ana_index_set(&gp, pc->b.nra, pc->b.a, NULL, 0);
    gmx_ana_index_set(&gb, base->b.nra, base->b.a, NULL, 0);
    isize = gmx_ana_index_difference_size(&gp, &gb);
    if (isize > 0)
    {
        gmx_ana_index_clear(&g);
        gmx_ana_index_reserve(&g, base->b.nra + isize);
        /* Find the new blocks */
        gmx_ana_index_difference(&g, &gp, &gb);
        /* Count the blocks in g */
        i = bi = bnr = 0;
        while (i < g.isize)
        {
            while (pc->b.a[pc->b.index[bi]] != g.index[i])
            {
                ++bi;
            }
            i += pc->b.index[bi+1] - pc->b.index[bi];
            ++bnr;
            ++bi;
        }
        /* Merge the atoms into a temporary structure */
        gmx_ana_index_merge(&g, &gb, &g);
        /* Merge the blocks */
        srenew(base->b.index, base->b.nr + bnr + 1);
        i                   = g.isize - 1;
        bi                  = base->b.nr - 1;
        bj                  = pc->b.nr - 1;
        bo                  = base->b.nr + bnr - 1;
        base->b.index[bo+1] = i + 1;
        while (bo >= 0)
        {
            if (bi < 0 || base->b.a[base->b.index[bi+1]-1] != g.index[i])
            {
                i -= pc->b.index[bj+1] - pc->b.index[bj];
                --bj;
            }
            else
            {
                if (bj >= 0 && pc->b.a[pc->b.index[bj+1]-1] == g.index[i])
                {
                    --bj;
                }
                i -= base->b.index[bi+1] - base->b.index[bi];
                --bi;
            }
            base->b.index[bo] = i + 1;
            --bo;
        }
        base->b.nr           += bnr;
        base->b.nalloc_index += bnr;
        sfree(base->b.a);
        base->b.nra      = g.isize;
        base->b.a        = g.index;
        base->b.nalloc_a = g.isize;
        /* Refresh the gmax field */
        gmx_ana_index_set(&base->gmax, base->b.nra, base->b.a, NULL, 0);
    }
}

/*! \brief
 * Merges two bases into one.
 *
 * \param[in,out] tbase Base calculation to merge to.
 * \param[in]     mbase Base calculation to merge to \p tbase.
 *
 * After the call, \p mbase has been freed and \p tbase is used as the base
 * for all calculations that previously had \p mbase as their base.
 * It is assumed that any overlap between \p tbase and \p mbase is in complete
 * blocks, i.e., that the merge is possible.
 */
static void
merge_bases(gmx_ana_poscalc_t *tbase, gmx_ana_poscalc_t *mbase)
{
    gmx_ana_poscalc_t *pc;

    merge_to_base(tbase, mbase);
    remove_poscalc(mbase);
    /* Set tbase as the base for all calculations that had mbase */
    pc = tbase->coll->first;
    while (pc)
    {
        if (pc->sbase == mbase)
        {
            pc->sbase = tbase;
            tbase->refcount++;
        }
        pc = pc->next;
    }
    /* Free mbase */
    mbase->refcount = 0;
    gmx_ana_poscalc_free(mbase);
}

/*! \brief
 * Setups the static base calculation for a position calculation.
 *
 * \param[in,out] pc   Position calculation to setup the base for.
 */
static void
setup_base(gmx_ana_poscalc_t *pc)
{
    gmx_ana_poscalc_t *base, *pbase, *next;
    gmx_ana_index_t    gp, g;

    /* Exit immediately if pc should not have a base. */
    if (!can_use_base(pc))
    {
        return;
    }

    gmx_ana_index_set(&gp, pc->b.nra, pc->b.a, NULL, 0);
    gmx_ana_index_clear(&g);
    gmx_ana_index_reserve(&g, pc->b.nra);
    pbase = pc;
    base  = pc->coll->first;
    while (base)
    {
        /* Save the next calculation so that we can safely delete base */
        next = base->next;
        /* Skip pc, calculations that already have a base (we should match the
         * base instead), as well as calculations that should not have a base.
         * If the above conditions are met, check whether we should do a
         * merge.
         */
        if (base != pc && !base->sbase && can_use_base(base)
            && should_merge(pbase, base, &gp, &g))
        {
            /* Check whether this is the first base found */
            if (pbase == pc)
            {
                /* Create a real base if one is not present */
                if (!base->p)
                {
                    pbase = create_simple_base(base);
                }
                else
                {
                    pbase = base;
                }
                /* Make it a base for pc as well */
                merge_to_base(pbase, pc);
                pc->sbase = pbase;
                pbase->refcount++;
            }
            else /* This was not the first base */
            {
                if (!base->p)
                {
                    /* If it is not a real base, just make the new base as
                     * the base for it as well. */
                    merge_to_base(pbase, base);
                    base->sbase = pbase;
                    pbase->refcount++;
                }
                else
                {
                    /* If base is a real base, merge it with the new base
                     * and delete it. */
                    merge_bases(pbase, base);
                }
            }
            gmx_ana_index_set(&gp, pbase->b.nra, pbase->b.a, NULL, 0);
            gmx_ana_index_reserve(&g, pc->b.nra);
        }
        /* Proceed to the next unchecked calculation */
        base = next;
    }

    gmx_ana_index_deinit(&g);

    /* If no base was found, create one if one is required */
    if (!pc->sbase && (pc->flags & POS_DYNAMIC)
        && (pc->flags & (POS_COMPLMAX | POS_COMPLWHOLE)))
    {
        create_simple_base(pc);
    }
}

/*!
 * \param[out] pcp   Position calculation data structure pointer to initialize.
 * \param[in,out] pcc   Position calculation collection.
 * \param[in]  type  Type of calculation.
 * \param[in]  flags Flags for setting calculation options
 *   (see \ref poscalc_flags "documentation of the flags").
 * \returns    0 on success.
 */
int
gmx_ana_poscalc_create(gmx_ana_poscalc_t **pcp, gmx_ana_poscalc_coll_t *pcc,
                       e_poscalc_t type, int flags)
{
    gmx_ana_poscalc_t *pc;

    snew(pc, 1);
    pc->type     = type;
    pc->itype    = index_type_for_poscalc(type);
    gmx_ana_poscalc_set_flags(pc, flags);
    pc->refcount = 1;
    pc->coll     = pcc;
    insert_poscalc(pc, NULL);
    *pcp = pc;
    return 0;
}

/*!
 * \param[out] pcp   Position calculation data structure pointer to initialize.
 * \param[in,out] pcc   Position calculation collection.
 * \param[in]  post  One of the strings acceptable for
 *   gmx_ana_poscalc_type_from_enum().
 * \param[in]  flags Flags for setting calculation options
 *   (see \ref poscalc_flags "documentation of the flags").
 * \returns    0 on success, a non-zero error value on error.
 *
 * This is a convenience wrapper for gmx_ana_poscalc_create().
 * \p flags sets the default calculation options if not overridden by \p post;
 * see gmx_ana_poscalc_type_from_enum().
 *
 * \see gmx_ana_poscalc_create(), gmx_ana_poscalc_type_from_enum()
 */
int
gmx_ana_poscalc_create_enum(gmx_ana_poscalc_t **pcp, gmx_ana_poscalc_coll_t *pcc,
                            const char *post, int flags)
{
    e_poscalc_t  type;
    int          cflags;
    int          rc;

    cflags = flags;
    rc     = gmx_ana_poscalc_type_from_enum(post, &type, &cflags);
    if (rc != 0)
    {
        *pcp = NULL;
        return rc;
    }
    return gmx_ana_poscalc_create(pcp, pcc, type, cflags);
}

/*!
 * \param[in,out] pc    Position calculation data structure.
 * \param[in]     flags New flags.
 *
 * \p flags are added to the old flags.
 * If calculation type is \ref POS_ATOM, \ref POS_MASS is automatically
 * cleared.
 * If both \ref POS_DYNAMIC and \ref POS_MASKONLY are provided,
 * \ref POS_DYNAMIC is cleared.
 * If calculation type is not \ref POS_RES or \ref POS_MOL,
 * \ref POS_COMPLMAX and \ref POS_COMPLWHOLE are automatically cleared.
 */
void
gmx_ana_poscalc_set_flags(gmx_ana_poscalc_t *pc, int flags)
{
    if (pc->type == POS_ATOM)
    {
        flags &= ~POS_MASS;
    }
    if (flags & POS_MASKONLY)
    {
        flags &= ~POS_DYNAMIC;
    }
    if (pc->type != POS_RES && pc->type != POS_MOL)
    {
        flags &= ~(POS_COMPLMAX | POS_COMPLWHOLE);
    }
    pc->flags |= flags;
}

/*!
 * \param[in,out] pc  Position calculation data structure.
 * \param[in]     g   Maximum index group for the calculation.
 *
 * Subsequent calls to gmx_ana_poscalc_update() should use only subsets of
 * \p g for evaluation.
 *
 * The topology should have been set for the collection of which \p pc is
 * a member.
 */
void
gmx_ana_poscalc_set_maxindex(gmx_ana_poscalc_t *pc, gmx_ana_index_t *g)
{
    set_poscalc_maxindex(pc, g, FALSE);
    setup_base(pc);
}

/*!
 * \param[in]  pc  Position calculation data structure.
 * \param[out] p   Output positions.
 *
 * Calls to gmx_ana_poscalc_update() using \p pc should use only positions
 * initialized with this function.
 * The \c p->g pointer is initialized to point to an internal group that
 * contains the maximum index group set with gmx_ana_poscalc_set_maxindex().
 */
void
gmx_ana_poscalc_init_pos(gmx_ana_poscalc_t *pc, gmx_ana_pos_t *p)
{
    gmx_ana_indexmap_init(&p->m, &pc->gmax, pc->coll->top, pc->itype);
    /* Only do the static optimization when there is no completion */
    if (!(pc->flags & POS_DYNAMIC) && pc->b.nra == pc->gmax.isize)
    {
        gmx_ana_indexmap_set_static(&p->m, &pc->b);
    }
    gmx_ana_pos_reserve(p, p->m.nr, 0);
    if (pc->flags & POS_VELOCITIES)
    {
        gmx_ana_pos_reserve_velocities(p);
    }
    if (pc->flags & POS_FORCES)
    {
        gmx_ana_pos_reserve_forces(p);
    }
    gmx_ana_pos_set_nr(p, p->m.nr);
    gmx_ana_pos_set_evalgrp(p, &pc->gmax);
}

/*!
 * \param  pc  Position calculation data to be freed.
 *
 * The \p pc pointer is invalid after the call.
 */
void
gmx_ana_poscalc_free(gmx_ana_poscalc_t *pc)
{
    if (!pc)
    {
        return;
    }

    pc->refcount--;
    if (pc->refcount > 0)
    {
        return;
    }

    remove_poscalc(pc);
    if (pc->b.nalloc_index > 0)
    {
        sfree(pc->b.index);
    }
    if (pc->b.nalloc_a > 0)
    {
        sfree(pc->b.a);
    }
    if (pc->flags & POS_COMPLWHOLE)
    {
        gmx_ana_index_deinit(&pc->gmax);
    }
    if (pc->p)
    {
        gmx_ana_pos_free(pc->p);
    }
    if (pc->sbase)
    {
        gmx_ana_poscalc_free(pc->sbase);
        sfree(pc->baseid);
    }
    sfree(pc);
}

/*!
 * \param[in] pc  Position calculation data to query.
 * \returns   TRUE if \p pc requires topology for initialization and/or
 *   evaluation, FALSE otherwise.
 */
gmx_bool
gmx_ana_poscalc_requires_top(gmx_ana_poscalc_t *pc)
{
    if ((pc->flags & POS_MASS) || pc->type == POS_RES || pc->type == POS_MOL)
    {
        return TRUE;
    }
    return FALSE;
}

/*!
 * \param[in,out] pcc Position calculation collection to initialize.
 *
 * This function does some final initialization of the data structures in the
 * collection to prepare them for evaluation.
 * After this function has been called, it is no longer possible to add new
 * calculations to the collection.
 *
 * This function is automatically called by gmx_ana_poscalc_init_frame()
 * if not called by the user earlier.
 * Multiple calls to the function are ignored.
 */
void
gmx_ana_poscalc_init_eval(gmx_ana_poscalc_coll_t *pcc)
{
    gmx_ana_poscalc_t      *pc;
    int                     bi, bj;

    if (pcc->bInit)
    {
        return;
    }
    pc = pcc->first;
    while (pc)
    {
        /* Initialize position storage for base calculations */
        if (pc->p)
        {
            gmx_ana_poscalc_init_pos(pc, pc->p);
        }
        /* Construct the mapping of the base positions */
        if (pc->sbase)
        {
            snew(pc->baseid, pc->b.nr);
            for (bi = bj = 0; bi < pc->b.nr; ++bi, ++bj)
            {
                while (pc->sbase->b.a[pc->sbase->b.index[bj]] != pc->b.a[pc->b.index[bi]])
                {
                    ++bj;
                }
                pc->baseid[bi] = bj;
            }
        }
        /* Free the block data for dynamic calculations */
        if (pc->flags & POS_DYNAMIC)
        {
            if (pc->b.nalloc_index > 0)
            {
                sfree(pc->b.index);
                pc->b.nalloc_index = 0;
            }
            if (pc->b.nalloc_a > 0)
            {
                sfree(pc->b.a);
                pc->b.nalloc_a = 0;
            }
        }
        pc = pc->next;
    }
    pcc->bInit = TRUE;
}

/*!
 * \param[in,out] pcc Position calculation collection to initialize.
 *
 * Clears the evaluation flag for all calculations.
 * Should be called for each frame before calling gmx_ana_poscalc_update().
 *
 * This function is automatically called by gmx_ana_do() for each
 * frame, and should not be called by the user unless gmx_ana_do() is
 * not being used.
 *
 * This function calls gmx_ana_poscalc_init_eval() automatically if it has
 * not been called earlier.
 */
void
gmx_ana_poscalc_init_frame(gmx_ana_poscalc_coll_t *pcc)
{
    gmx_ana_poscalc_t      *pc;

    if (!pcc->bInit)
    {
        gmx_ana_poscalc_init_eval(pcc);
    }
    /* Clear the evaluation flags */
    pc = pcc->first;
    while (pc)
    {
        pc->bEval = FALSE;
        pc        = pc->next;
    }
}

/*!
 * \param[in]     pc   Position calculation data.
 * \param[in,out] p    Output positions, initialized previously with
 *   gmx_ana_poscalc_init_pos() using \p pc.
 * \param[in]     g    Index group to use for the update.
 * \param[in]     fr   Current frame.
 * \param[in]     pbc  PBC data, or NULL if no PBC should be used.
 *
 * gmx_ana_poscalc_init_frame() should be called for each frame before calling
 * this function.
 */
void
gmx_ana_poscalc_update(gmx_ana_poscalc_t *pc, gmx_ana_pos_t *p,
                       gmx_ana_index_t *g, t_trxframe *fr, t_pbc *pbc)
{
    int  i, j, bi, bj;

    if (pc->bEval == TRUE && !(pc->flags & POS_MASKONLY))
    {
        return;
    }
    if (pc->sbase)
    {
        gmx_ana_poscalc_update(pc->sbase, NULL, NULL, fr, pbc);
    }
    if (!p)
    {
        p = pc->p;
    }
    if (!g)
    {
        g = &pc->gmax;
    }
    gmx_ana_pos_set_evalgrp(p, g);

    /* Update the index map */
    if (pc->flags & POS_DYNAMIC)
    {
        gmx_ana_indexmap_update(&p->m, g, FALSE);
        p->nr = p->m.nr;
    }
    else if (pc->flags & POS_MASKONLY)
    {
        gmx_ana_indexmap_update(&p->m, g, TRUE);
        if (pc->bEval)
        {
            return;
        }
    }
    if (!(pc->flags & POS_DYNAMIC))
    {
        pc->bEval = TRUE;
    }

    /* Evaluate the positions */
    if (pc->sbase)
    {
        /* TODO: It might be faster to evaluate the positions within this
         * loop instead of in the beginning. */
        if (pc->flags & POS_DYNAMIC)
        {
            for (bi = 0; bi < p->nr; ++bi)
            {
                bj = pc->baseid[p->m.refid[bi]];
                copy_rvec(pc->sbase->p->x[bj], p->x[bi]);
            }
            if (p->v)
            {
                for (bi = 0; bi < p->nr; ++bi)
                {
                    bj = pc->baseid[p->m.refid[bi]];
                    copy_rvec(pc->sbase->p->v[bj], p->v[bi]);
                }
            }
            if (p->f)
            {
                for (bi = 0; bi < p->nr; ++bi)
                {
                    bj = pc->baseid[p->m.refid[bi]];
                    copy_rvec(pc->sbase->p->f[bj], p->f[bi]);
                }
            }
        }
        else
        {
            for (bi = 0; bi < p->nr; ++bi)
            {
                bj = pc->baseid[bi];
                copy_rvec(pc->sbase->p->x[bj], p->x[bi]);
            }
            if (p->v)
            {
                for (bi = 0; bi < p->nr; ++bi)
                {
                    bj = pc->baseid[bi];
                    copy_rvec(pc->sbase->p->v[bj], p->v[bi]);
                }
            }
            if (p->f)
            {
                for (bi = 0; bi < p->nr; ++bi)
                {
                    bj = pc->baseid[bi];
                    copy_rvec(pc->sbase->p->f[bj], p->f[bi]);
                }
            }
        }
    }
    else /* pc->sbase is NULL */
    {
        if (pc->flags & POS_DYNAMIC)
        {
            pc->b.nr    = p->m.mapb.nr;
            pc->b.index = p->m.mapb.index;
            pc->b.nra   = g->isize;
            pc->b.a     = g->index;
        }
        if (p->v && !fr->bV)
        {
            for (i = 0; i < pc->b.nra; ++i)
            {
                clear_rvec(p->v[i]);
            }
        }
        if (p->f && !fr->bF)
        {
            for (i = 0; i < pc->b.nra; ++i)
            {
                clear_rvec(p->f[i]);
            }
        }
        /* Here, we assume that the topology has been properly initialized,
         * and do not check the return values of gmx_calc_comg*(). */
        switch (pc->type)
        {
            case POS_ATOM:
                for (i = 0; i < pc->b.nra; ++i)
                {
                    copy_rvec(fr->x[pc->b.a[i]], p->x[i]);
                }
                if (p->v && fr->bV)
                {
                    for (i = 0; i < pc->b.nra; ++i)
                    {
                        copy_rvec(fr->v[pc->b.a[i]], p->v[i]);
                    }
                }
                if (p->f && fr->bF)
                {
                    for (i = 0; i < pc->b.nra; ++i)
                    {
                        copy_rvec(fr->f[pc->b.a[i]], p->f[i]);
                    }
                }
                break;
            case POS_ALL:
                gmx_calc_comg(pc->coll->top, fr->x, pc->b.nra, pc->b.a,
                              pc->flags & POS_MASS, p->x[0]);
                if (p->v && fr->bV)
                {
                    gmx_calc_comg(pc->coll->top, fr->v, pc->b.nra, pc->b.a,
                                  pc->flags & POS_MASS, p->v[0]);
                }
                if (p->f && fr->bF)
                {
                    gmx_calc_comg_f(pc->coll->top, fr->f, pc->b.nra, pc->b.a,
                                    pc->flags & POS_MASS, p->f[0]);
                }
                break;
            case POS_ALL_PBC:
                gmx_calc_comg_pbc(pc->coll->top, fr->x, pbc, pc->b.nra, pc->b.a,
                                  pc->flags & POS_MASS, p->x[0]);
                if (p->v && fr->bV)
                {
                    gmx_calc_comg(pc->coll->top, fr->v, pc->b.nra, pc->b.a,
                                  pc->flags & POS_MASS, p->v[0]);
                }
                if (p->f && fr->bF)
                {
                    gmx_calc_comg_f(pc->coll->top, fr->f, pc->b.nra, pc->b.a,
                                    pc->flags & POS_MASS, p->f[0]);
                }
                break;
            default:
                gmx_calc_comg_blocka(pc->coll->top, fr->x, &pc->b,
                                     pc->flags & POS_MASS, p->x);
                if (p->v && fr->bV)
                {
                    gmx_calc_comg_blocka(pc->coll->top, fr->v, &pc->b,
                                         pc->flags & POS_MASS, p->v);
                }
                if (p->f && fr->bF)
                {
                    gmx_calc_comg_blocka(pc->coll->top, fr->f, &pc->b,
                                         pc->flags & POS_MASS, p->f);
                }
                break;
        }
    }
}
