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
 * \brief Implementation of functions in selelem.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <smalloc.h>
#include <gmx_fatal.h>

#include <indexutil.h>
#include <poscalc.h>
#include <position.h>
#include <selmethod.h>

#include "keywords.h"
#include "mempool.h"
#include "selelem.h"

/*!
 * \param[in] sel Selection for which the string is requested
 * \returns   Pointer to a string that corresponds to \p sel->type.
 *
 * The return value points to a string constant and should not be \p free'd.
 *
 * The function returns NULL if \p sel->type is not one of the valid values.
 */
const char *
_gmx_selelem_type_str(t_selelem *sel)
{
    switch (sel->type)
    {
        case SEL_CONST:      return "CONST";
        case SEL_EXPRESSION: return "EXPR";
        case SEL_BOOLEAN:    return "BOOL";
        case SEL_ARITHMETIC: return "ARITH";
        case SEL_ROOT:       return "ROOT";
        case SEL_SUBEXPR:    return "SUBEXPR";
        case SEL_SUBEXPRREF: return "REF";
        case SEL_MODIFIER:   return "MODIFIER";
    }
    return NULL;
}

/*!
 * \param[in] val Value structore for which the string is requested.
 * \returns   Pointer to a string that corresponds to \p val->type,
 *   NULL if the type value is invalid.
 *
 * The return value points to a string constant and should not be \p free'd.
 */
const char *
_gmx_sel_value_type_str(gmx_ana_selvalue_t *val)
{
    switch (val->type)
    {
        case NO_VALUE:       return "NONE";
        case INT_VALUE:      return "INT";
        case REAL_VALUE:     return "REAL";
        case STR_VALUE:      return "STR";
        case POS_VALUE:      return "VEC";
        case GROUP_VALUE:    return "GROUP";
    }
    return NULL;
}

/*! \copydoc _gmx_selelem_type_str() */
const char *
_gmx_selelem_gmx_boolean_type_str(t_selelem *sel)
{
    switch (sel->u.boolt)
    {
        case BOOL_NOT:  return "NOT"; break;
        case BOOL_AND:  return "AND"; break;
        case BOOL_OR:   return "OR";  break;
        case BOOL_XOR:  return "XOR"; break;
    }
    return NULL;
}

/*!
 * \param[in] type Type of selection element to allocate.
 * \returns   Pointer to the newly allocated and initialized element.
 *
 * \c t_selelem::type is set to \p type,
 * \c t_selelem::v::type is set to \ref GROUP_VALUE for gmx_boolean and comparison
 * expressions and \ref NO_VALUE for others,
 * \ref SEL_ALLOCVAL is set for non-root elements (\ref SEL_ALLOCDATA is also
 * set for \ref SEL_BOOLEAN elements),
 * and \c t_selelem::refcount is set to one.
 * All the pointers are set to NULL.
 */
t_selelem *
_gmx_selelem_create(e_selelem_t type)
{
    t_selelem *sel;

    snew(sel, 1);
    sel->name       = NULL;
    sel->type       = type;
    sel->flags      = (type != SEL_ROOT) ? SEL_ALLOCVAL : 0;
    if (type == SEL_BOOLEAN)
    {
        sel->v.type = GROUP_VALUE;
        sel->flags |= SEL_ALLOCDATA;
    }
    else
    {
        sel->v.type = NO_VALUE;
    }
    _gmx_selvalue_clear(&sel->v);
    sel->evaluate   = NULL;
    sel->mempool    = NULL;
    sel->child      = NULL;
    sel->next       = NULL;
    sel->refcount   = 1;

    return sel;
}

/*!
 * \param[in,out] sel   Selection element to set the type for.
 * \param[in]     vtype Value type for the selection element.
 * \returns       0 on success, EINVAL if the value type is invalid.
 *
 * If the new type is \ref GROUP_VALUE or \ref POS_VALUE, the
 * \ref SEL_ALLOCDATA flag is also set.
 *
 * This function should only be called at most once for each element,
 * preferably right after calling _gmx_selelem_create().
 */
int
_gmx_selelem_set_vtype(t_selelem *sel, e_selvalue_t vtype)
{
    if (sel->type == SEL_BOOLEAN && vtype != GROUP_VALUE)
    {
        gmx_bug("internal error");
        return EINVAL;
    }
    if (sel->v.type != NO_VALUE && vtype != sel->v.type)
    {
        gmx_call("_gmx_selelem_set_vtype() called more than once");
        return EINVAL;
    }
    sel->v.type = vtype;
    if (vtype == GROUP_VALUE || vtype == POS_VALUE)
    {
        sel->flags |= SEL_ALLOCDATA;
    }
    return 0;
}

/*!
 * \param[in,out] sel   Selection element to reserve.
 * \param[in]     count Number of values to reserve memory for.
 * \returns       0 on success or if no memory pool, non-zero on error.
 *
 * Reserves memory for the values of \p sel from the \p sel->mempool
 * memory pool. If no memory pool is set, nothing is done.
 */
int
_gmx_selelem_mempool_reserve(t_selelem *sel, int count)
{
    int rc = 0;

    if (!sel->mempool)
    {
        return 0;
    }
    switch (sel->v.type)
    {
        case INT_VALUE:
            rc = _gmx_sel_mempool_alloc(sel->mempool, (void **)&sel->v.u.i,
                                        sizeof(*sel->v.u.i)*count);
            break;

        case REAL_VALUE:
            rc = _gmx_sel_mempool_alloc(sel->mempool, (void **)&sel->v.u.r,
                                        sizeof(*sel->v.u.r)*count);
            break;

        case GROUP_VALUE:
            rc = _gmx_sel_mempool_alloc_group(sel->mempool, sel->v.u.g, count);
            break;

        default:
            gmx_incons("mem pooling not implemented for requested type");
            return -1;
    }
    return rc;
}

/*!
 * \param[in,out] sel   Selection element to release.
 *
 * Releases the memory allocated for the values of \p sel from the
 * \p sel->mempool memory pool. If no memory pool is set, nothing is done.
 */
void
_gmx_selelem_mempool_release(t_selelem *sel)
{
    if (!sel->mempool)
    {
        return;
    }
    switch (sel->v.type)
    {
        case INT_VALUE:
        case REAL_VALUE:
            _gmx_sel_mempool_free(sel->mempool, sel->v.u.ptr);
            _gmx_selvalue_setstore(&sel->v, NULL);
            break;

        case GROUP_VALUE:
            if (sel->v.u.g)
            {
                _gmx_sel_mempool_free_group(sel->mempool, sel->v.u.g);
            }
            break;

        default:
            gmx_incons("mem pooling not implemented for requested type");
            break;
    }
}

/*!
 * \param[in] sel Selection to free.
 */
void
_gmx_selelem_free_values(t_selelem *sel)
{
    int   i, n;

    _gmx_selelem_mempool_release(sel);
    if ((sel->flags & SEL_ALLOCDATA) && sel->v.u.ptr)
    {
        /* The number of position/group structures is constant, so the
         * backup of using sel->v.nr should work for them.
         * For strings, we report an error if we don't know the allocation
         * size here. */
        n = (sel->v.nalloc > 0) ? sel->v.nalloc : sel->v.nr;
        switch (sel->v.type)
        {
            case STR_VALUE:
                if (sel->v.nalloc == 0)
                {
                    gmx_bug("SEL_ALLOCDATA should only be set for allocated STR_VALUE values");
                    break;
                }
                for (i = 0; i < n; ++i)
                {
                    sfree(sel->v.u.s[i]);
                }
                break;
            case POS_VALUE:
                for (i = 0; i < n; ++i)
                {
                    gmx_ana_pos_deinit(&sel->v.u.p[i]);
                }
                break;
            case GROUP_VALUE:
                for (i = 0; i < n; ++i)
                {
                    gmx_ana_index_deinit(&sel->v.u.g[i]);
                }
                break;
            default: /* No special handling for other types */
                break;
        }
    }
    if (sel->flags & SEL_ALLOCVAL)
    {
        sfree(sel->v.u.ptr);
    }
    _gmx_selvalue_setstore(&sel->v, NULL);
    if (sel->type == SEL_SUBEXPRREF && sel->u.param)
    {
        sel->u.param->val.u.ptr = NULL;
    }
}

/*!
 * \param[in] method Method to free.
 * \param[in] mdata  Method data to free.
 */
void
_gmx_selelem_free_method(gmx_ana_selmethod_t *method, void *mdata)
{
    sel_freefunc free_func = NULL;

    /* Save the pointer to the free function. */
    if (method && method->free)
    {
        free_func = method->free;
    }

    /* Free the method itself.
     * Has to be done before freeing the method data, because parameter
     * values are typically stored in the method data, and here we may
     * access them. */
    if (method)
    {
        int  i, j;

        /* Free the memory allocated for the parameters that are not managed
         * by the selection method itself. */
        for (i = 0; i < method->nparams; ++i)
        {
            gmx_ana_selparam_t *param = &method->param[i];

            if (param->val.u.ptr)
            {
                if (param->val.type == GROUP_VALUE)
                {
                    for (j = 0; j < param->val.nr; ++j)
                    {
                        gmx_ana_index_deinit(&param->val.u.g[j]);
                    }
                }
                else if (param->val.type == POS_VALUE)
                {
                    for (j = 0; j < param->val.nr; ++j)
                    {
                        gmx_ana_pos_deinit(&param->val.u.p[j]);
                    }
                }

                if (param->val.nalloc > 0)
                {
                    sfree(param->val.u.ptr);
                }
            }
        }
        sfree(method->param);
        sfree(method);
    }
    /* Free method data. */
    if (mdata)
    {
        if (free_func)
        {
            free_func(mdata);
        }
        sfree(mdata);
    }
}

/*!
 * \param[in] sel Selection to free.
 */
void
_gmx_selelem_free_exprdata(t_selelem *sel)
{
    if (sel->type == SEL_EXPRESSION || sel->type == SEL_MODIFIER)
    {
        _gmx_selelem_free_method(sel->u.expr.method, sel->u.expr.mdata);
        sel->u.expr.mdata  = NULL;
        sel->u.expr.method = NULL;
        /* Free position data */
        if (sel->u.expr.pos)
        {
            gmx_ana_pos_free(sel->u.expr.pos);
            sel->u.expr.pos = NULL;
        }
        /* Free position calculation data */
        if (sel->u.expr.pc)
        {
            gmx_ana_poscalc_free(sel->u.expr.pc);
            sel->u.expr.pc = NULL;
        }
    }
    if (sel->type == SEL_ARITHMETIC)
    {
        sfree(sel->u.arith.opstr);
        sel->u.arith.opstr = NULL;
    }
    if (sel->type == SEL_SUBEXPR || sel->type == SEL_ROOT
        || (sel->type == SEL_CONST && sel->v.type == GROUP_VALUE))
    {
        gmx_ana_index_deinit(&sel->u.cgrp);
    }
}

/*!
 * \param[in] sel Selection to free.
 *
 * Decrements \ref t_selelem::refcount "sel->refcount" and frees the
 * memory allocated for \p sel and all its children if the reference count
 * reaches zero.
 */
void
_gmx_selelem_free(t_selelem *sel)
{
    /* Decrement the reference counter and do nothing if references remain */
    sel->refcount--;
    if (sel->refcount > 0)
    {
        return;
    }

    /* Free the children.
     * Must be done before freeing other data, because the children may hold
     * references to data in this element. */
    _gmx_selelem_free_chain(sel->child);

    /* Free value storage */
    _gmx_selelem_free_values(sel);

    /* Free other storage */
    _gmx_selelem_free_exprdata(sel);

    /* Free temporary compiler data if present */
    _gmx_selelem_free_compiler_data(sel);

    sfree(sel);
}

/*!
 * \param[in] first First selection to free.
 *
 * Frees \p first and all selections accessible through the
 * \ref t_selelem::next "first->next" pointer.
 */
void
_gmx_selelem_free_chain(t_selelem *first)
{
    t_selelem *child, *prev;

    child = first;
    while (child)
    {
        prev  = child;
        child = child->next;
        _gmx_selelem_free(prev);
    }
}

/*!
 * \param[in] fp      File handle to receive the output.
 * \param[in] sel     Root of the selection subtree to print.
 * \param[in] bValues If TRUE, the evaluated values of selection elements
 *   are printed as well.
 * \param[in] level   Indentation level, starting from zero.
 */
void
_gmx_selelem_print_tree(FILE *fp, t_selelem *sel, gmx_bool bValues, int level)
{
    t_selelem   *child;
    int          i;

    fprintf(fp, "%*c %s %s", level*2+1, '*',
            _gmx_selelem_type_str(sel), _gmx_sel_value_type_str(&sel->v));
    if (sel->name)
    {
        fprintf(fp, " \"%s\"", sel->name);
    }
    fprintf(fp, " flg=");
    if (sel->flags & SEL_FLAGSSET)
    {
        fprintf(fp, "s");
    }
    if (sel->flags & SEL_SINGLEVAL)
    {
        fprintf(fp, "S");
    }
    if (sel->flags & SEL_ATOMVAL)
    {
        fprintf(fp, "A");
    }
    if (sel->flags & SEL_VARNUMVAL)
    {
        fprintf(fp, "V");
    }
    if (sel->flags & SEL_DYNAMIC)
    {
        fprintf(fp, "D");
    }
    if (!(sel->flags & SEL_VALFLAGMASK))
    {
        fprintf(fp, "0");
    }
    if (sel->mempool)
    {
        fprintf(fp, "P");
    }
    if (sel->type == SEL_CONST)
    {
        if (sel->v.type == INT_VALUE)
        {
            fprintf(fp, " %d", sel->v.u.i[0]);
        }
        else if (sel->v.type == REAL_VALUE)
        {
            fprintf(fp, " %f", sel->v.u.r[0]);
        }
        else if (sel->v.type == GROUP_VALUE)
        {
            gmx_ana_index_t *g = sel->v.u.g;
            if (!g || g->isize == 0)
            {
                g = &sel->u.cgrp;
            }
            fprintf(fp, " (%d atoms)", g->isize);
        }
    }
    else if (sel->type == SEL_BOOLEAN)
    {
        fprintf(fp, " %s", _gmx_selelem_gmx_boolean_type_str(sel));
    }
    else if (sel->type == SEL_EXPRESSION
             && sel->u.expr.method->name == sm_compare.name)
    {
        _gmx_selelem_print_compare_info(fp, sel->u.expr.mdata);
    }
    if (sel->evaluate)
    {
        fprintf(fp, " eval=");
        _gmx_sel_print_evalfunc_name(fp, sel->evaluate);
    }
    if (sel->refcount > 1)
    {
        fprintf(fp, " refc=%d", sel->refcount);
    }
    if (!(sel->flags & SEL_ALLOCVAL))
    {
        fprintf(fp, " (ext. output)");
    }
    fprintf(fp, "\n");

    if ((sel->type == SEL_CONST && sel->v.type == GROUP_VALUE) || sel->type == SEL_ROOT)
    {
        gmx_ana_index_t *g = sel->v.u.g;
        if (!g || g->isize == 0 || sel->evaluate != NULL)
        {
            g = &sel->u.cgrp;
        }
        if (g->isize < 0)
        {
            fprintf(fp, "%*c group: (null)\n", level*2+1, ' ');
        }
        else if (g->isize > 0)
        {
            fprintf(fp, "%*c group:", level*2+1, ' ');
            if (g->isize <= 20)
            {
                for (i = 0; i < g->isize; ++i)
                {
                    fprintf(fp, " %d", g->index[i] + 1);
                }
            }
            else
            {
                fprintf(fp, " %d atoms", g->isize);
            }
            fprintf(fp, "\n");
        }
    }
    else if (sel->type == SEL_EXPRESSION)
    {
        if (sel->u.expr.pc)
        {
            fprintf(fp, "%*c COM", level*2+3, '*');
            fprintf(fp, "\n");
        }
    }

    if (sel->cdata)
    {
        _gmx_selelem_print_compiler_info(fp, sel, level);
    }

    if (bValues && sel->type != SEL_CONST && sel->type != SEL_ROOT && sel->v.u.ptr)
    {
        fprintf(fp, "%*c value: ", level*2+1, ' ');
        switch (sel->v.type)
        {
            case POS_VALUE:
                /* In normal use, the pointer should never be NULL, but it's
                 * useful to have the check for debugging to avoid accidental
                 * segfaults when printing the selection tree. */
                if (sel->v.u.p->x)
                {
                    fprintf(fp, "(%f, %f, %f)",
                            sel->v.u.p->x[0][XX], sel->v.u.p->x[0][YY],
                            sel->v.u.p->x[0][ZZ]);
                }
                else
                {
                    fprintf(fp, "(null)");
                }
                break;
            case GROUP_VALUE:
                fprintf(fp, "%d atoms", sel->v.u.g->isize);
                if (sel->v.u.g->isize < 20)
                {
                    if (sel->v.u.g->isize > 0)
                    {
                        fprintf(fp, ":");
                    }
                    for (i = 0; i < sel->v.u.g->isize; ++i)
                    {
                        fprintf(fp, " %d", sel->v.u.g->index[i] + 1);
                    }
                }
                break;
            default:
                fprintf(fp, "???");
                break;
        }
        fprintf(fp, "\n");
    }

    /* Print the subexpressions with one more level of indentation */
    child = sel->child;
    while (child)
    {
        if (!(sel->type == SEL_SUBEXPRREF && child->type == SEL_SUBEXPR))
        {
            _gmx_selelem_print_tree(fp, child, bValues, level+1);
        }
        child = child->next;
    }
}

/*!
 * \param[in] root Root of the subtree to query.
 * \returns TRUE if \p root or any any of its elements require topology
 *   information, FALSE otherwise.
 */
gmx_bool
_gmx_selelem_requires_top(t_selelem *root)
{
    t_selelem *child;

    if (root->type == SEL_EXPRESSION || root->type == SEL_MODIFIER)
    {
        if (root->u.expr.method && (root->u.expr.method->flags & SMETH_REQTOP))
        {
            return TRUE;
        }
        if (root->u.expr.pc && gmx_ana_poscalc_requires_top(root->u.expr.pc))
        {
            return TRUE;
        }
    }
    child = root->child;
    while (child)
    {
        if (_gmx_selelem_requires_top(child))
        {
            return TRUE;
        }
        child = child->next;
    }
    return FALSE;
}
