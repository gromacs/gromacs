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
/*! \internal \file
 * \brief
 * Implements functions in selelem.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#include <cstring>

#include "gromacs/legacyheaders/smalloc.h"

#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/poscalc.h"
#include "gromacs/selection/position.h"
#include "gromacs/selection/selmethod.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

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
_gmx_selelem_type_str(const gmx::SelectionTreeElement &sel)
{
    switch (sel.type)
    {
        case SEL_CONST:      return "CONST";
        case SEL_EXPRESSION: return "EXPR";
        case SEL_BOOLEAN:    return "BOOL";
        case SEL_ARITHMETIC: return "ARITH";
        case SEL_ROOT:       return "ROOT";
        case SEL_SUBEXPR:    return "SUBEXPR";
        case SEL_SUBEXPRREF: return "REF";
        case SEL_GROUPREF:   return "GROUPREF";
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
_gmx_sel_value_type_str(const gmx_ana_selvalue_t *val)
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
_gmx_selelem_boolean_type_str(const gmx::SelectionTreeElement &sel)
{
    switch (sel.u.boolt)
    {
        case BOOL_NOT:  return "NOT"; break;
        case BOOL_AND:  return "AND"; break;
        case BOOL_OR:   return "OR";  break;
        case BOOL_XOR:  return "XOR"; break;
    }
    return NULL;
}


namespace gmx
{

SelectionTreeElement::SelectionTreeElement(e_selelem_t type)
{
    this->type       = type;
    this->flags      = (type != SEL_ROOT) ? SEL_ALLOCVAL : 0;
    if (type == SEL_BOOLEAN)
    {
        this->v.type = GROUP_VALUE;
        this->flags |= SEL_ALLOCDATA;
    }
    else
    {
        this->v.type = NO_VALUE;
    }
    _gmx_selvalue_clear(&this->v);
    std::memset(&this->u, 0, sizeof(this->u));
    this->evaluate   = NULL;
    this->mempool    = NULL;
    this->cdata      = NULL;
}

SelectionTreeElement::~SelectionTreeElement()
{
    /* Free the children.
     * Must be done before freeing other data, because the children may hold
     * references to data in this element. */
    child.reset();

    freeValues();
    freeExpressionData();
    freeCompilerData();
}

void SelectionTreeElement::freeValues()
{
    mempoolRelease();
    if ((flags & SEL_ALLOCDATA) && v.u.ptr)
    {
        /* The number of position/group structures is constant, so the
         * backup of using sel->v.nr should work for them.
         * For strings, we report an error if we don't know the allocation
         * size here. */
        int n = (v.nalloc > 0) ? v.nalloc : v.nr;
        switch (v.type)
        {
            case STR_VALUE:
                GMX_RELEASE_ASSERT(v.nalloc != 0,
                        "SEL_ALLOCDATA should only be set for allocated "
                        "STR_VALUE values");
                for (int i = 0; i < n; ++i)
                {
                    sfree(v.u.s[i]);
                }
                break;
            case POS_VALUE:
                for (int i = 0; i < n; ++i)
                {
                    gmx_ana_pos_deinit(&v.u.p[i]);
                }
                break;
            case GROUP_VALUE:
                for (int i = 0; i < n; ++i)
                {
                    gmx_ana_index_deinit(&v.u.g[i]);
                }
                break;
            default: /* No special handling for other types */
                break;
        }
    }
    if (flags & SEL_ALLOCVAL)
    {
        sfree(v.u.ptr);
    }
    _gmx_selvalue_setstore(&v, NULL);
    if (type == SEL_SUBEXPRREF && u.param)
    {
        u.param->val.u.ptr = NULL;
    }
}

void
SelectionTreeElement::freeExpressionData()
{
    if (type == SEL_EXPRESSION || type == SEL_MODIFIER)
    {
        _gmx_selelem_free_method(u.expr.method, u.expr.mdata);
        u.expr.mdata = NULL;
        u.expr.method = NULL;
        /* Free position data */
        if (u.expr.pos)
        {
            gmx_ana_pos_free(u.expr.pos);
            u.expr.pos = NULL;
        }
        /* Free position calculation data */
        if (u.expr.pc)
        {
            gmx_ana_poscalc_free(u.expr.pc);
            u.expr.pc = NULL;
        }
    }
    if (type == SEL_ARITHMETIC)
    {
        sfree(u.arith.opstr);
        u.arith.opstr = NULL;
    }
    if (type == SEL_SUBEXPR || type == SEL_ROOT
        || (type == SEL_CONST && v.type == GROUP_VALUE))
    {
        gmx_ana_index_deinit(&u.cgrp);
    }
    if (type == SEL_GROUPREF)
    {
        sfree(u.gref.name);
    }
}

void SelectionTreeElement::mempoolReserve(int count)
{
    if (!mempool)
    {
        return;
    }
    switch (v.type)
    {
        case INT_VALUE:
            v.u.i = static_cast<int *>(
                    _gmx_sel_mempool_alloc(mempool, sizeof(*v.u.i)*count));
            break;

        case REAL_VALUE:
            v.u.r = static_cast<real *>(
                    _gmx_sel_mempool_alloc(mempool, sizeof(*v.u.r)*count));
            break;

        case GROUP_VALUE:
            _gmx_sel_mempool_alloc_group(mempool, v.u.g, count);
            break;

        default:
            GMX_THROW(gmx::InternalError("Memory pooling not implemented for requested type"));
    }
}

void SelectionTreeElement::mempoolRelease()
{
    if (!mempool)
    {
        return;
    }
    switch (v.type)
    {
        case INT_VALUE:
        case REAL_VALUE:
            _gmx_sel_mempool_free(mempool, v.u.ptr);
            _gmx_selvalue_setstore(&v, NULL);
            break;

        case GROUP_VALUE:
            if (v.u.g)
            {
                _gmx_sel_mempool_free_group(mempool, v.u.g);
            }
            break;

        default:
            GMX_THROW(gmx::InternalError("Memory pooling not implemented for requested type"));
    }
}

} // namespace gmx

/*!
 * \param[in,out] sel   Selection element to set the type for.
 * \param[in]     vtype Value type for the selection element.
 *
 * If the new type is \ref GROUP_VALUE or \ref POS_VALUE, the
 * \ref SEL_ALLOCDATA flag is also set.
 *
 * This function should only be called at most once for each element,
 * preferably right after calling _gmx_selelem_create().
 */
void
_gmx_selelem_set_vtype(const gmx::SelectionTreeElementPointer &sel,
                       e_selvalue_t vtype)
{
    GMX_RELEASE_ASSERT(sel->type != SEL_BOOLEAN || vtype == GROUP_VALUE,
                       "Boolean elements must have a group value");
    GMX_RELEASE_ASSERT(sel->v.type == NO_VALUE || vtype == sel->v.type,
                       "_gmx_selelem_set_vtype() called more than once");
    sel->v.type = vtype;
    if (vtype == GROUP_VALUE || vtype == POS_VALUE)
    {
        sel->flags |= SEL_ALLOCDATA;
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
        else
        {
            sfree(mdata);
        }
    }
}

/*!
 * \param[in] fp      File handle to receive the output.
 * \param[in] sel     Root of the selection subtree to print.
 * \param[in] bValues If true, the evaluated values of selection elements
 *   are printed as well.
 * \param[in] level   Indentation level, starting from zero.
 */
void
_gmx_selelem_print_tree(FILE *fp, const gmx::SelectionTreeElement &sel,
                        bool bValues, int level)
{
    int          i;

    fprintf(fp, "%*c %s %s", level*2+1, '*',
            _gmx_selelem_type_str(sel), _gmx_sel_value_type_str(&sel.v));
    if (!sel.name().empty())
    {
        fprintf(fp, " \"%s\"", sel.name().c_str());
    }
    fprintf(fp, " flg=");
    if (sel.flags & SEL_FLAGSSET)
    {
        fprintf(fp, "s");
    }
    if (sel.flags & SEL_SINGLEVAL)
    {
        fprintf(fp, "S");
    }
    if (sel.flags & SEL_ATOMVAL)
    {
        fprintf(fp, "A");
    }
    if (sel.flags & SEL_VARNUMVAL)
    {
        fprintf(fp, "V");
    }
    if (sel.flags & SEL_DYNAMIC)
    {
        fprintf(fp, "D");
    }
    if (!(sel.flags & SEL_VALFLAGMASK))
    {
        fprintf(fp, "0");
    }
    if (sel.mempool)
    {
        fprintf(fp, "P");
    }
    if (sel.type == SEL_CONST)
    {
        if (sel.v.type == INT_VALUE)
        {
            fprintf(fp, " %d", sel.v.u.i[0]);
        }
        else if (sel.v.type == REAL_VALUE)
        {
            fprintf(fp, " %f", sel.v.u.r[0]);
        }
        else if (sel.v.type == GROUP_VALUE)
        {
            const gmx_ana_index_t *g = sel.v.u.g;
            if (!g || g->isize == 0)
                g = &sel.u.cgrp;
            fprintf(fp, " (%d atoms)", g->isize);
        }
    }
    else if (sel.type == SEL_BOOLEAN)
    {
        fprintf(fp, " %s", _gmx_selelem_boolean_type_str(sel));
    }
    else if (sel.type == SEL_EXPRESSION
             && sel.u.expr.method->name == sm_compare.name)
    {
        _gmx_selelem_print_compare_info(fp, sel.u.expr.mdata);
    }
    if (sel.evaluate)
    {
        fprintf(fp, " eval=");
        _gmx_sel_print_evalfunc_name(fp, sel.evaluate);
    }
    if (!(sel.flags & SEL_ALLOCVAL))
    {
        fprintf(fp, " (ext. output)");
    }
    fprintf(fp, "\n");

    if ((sel.type == SEL_CONST && sel.v.type == GROUP_VALUE) || sel.type == SEL_ROOT)
    {
        const gmx_ana_index_t *g = sel.v.u.g;
        if (!g || g->isize == 0 || sel.evaluate != NULL)
        {
            g = &sel.u.cgrp;
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
    else if (sel.type == SEL_EXPRESSION)
    {
        if (sel.u.expr.pc)
        {
            fprintf(fp, "%*c COM", level*2+3, '*');
            fprintf(fp, "\n");
        }
    }

    if (sel.cdata)
    {
        _gmx_selelem_print_compiler_info(fp, sel, level);
    }

    if (bValues && sel.type != SEL_CONST && sel.type != SEL_ROOT && sel.v.u.ptr)
    {
        fprintf(fp, "%*c value: ", level*2+1, ' ');
        switch (sel.v.type)
        {
            case POS_VALUE:
                /* In normal use, the pointer should never be NULL, but it's
                 * useful to have the check for debugging to avoid accidental
                 * segfaults when printing the selection tree. */
                if (sel.v.u.p->x)
                {
                    fprintf(fp, "(%f, %f, %f)",
                            sel.v.u.p->x[0][XX], sel.v.u.p->x[0][YY],
                            sel.v.u.p->x[0][ZZ]);
                }
                else
                {
                    fprintf(fp, "(null)");
                }
                break;
            case GROUP_VALUE:
                fprintf(fp, "%d atoms", sel.v.u.g->isize);
                if (sel.v.u.g->isize < 20)
                {
                    if (sel.v.u.g->isize > 0)
                    {
                        fprintf(fp, ":");
                    }
                    for (i = 0; i < sel.v.u.g->isize; ++i)
                    {
                        fprintf(fp, " %d", sel.v.u.g->index[i] + 1);
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
    gmx::SelectionTreeElementPointer child = sel.child;
    while (child)
    {
        if (!(sel.type == SEL_SUBEXPRREF && child->type == SEL_SUBEXPR))
        {
            _gmx_selelem_print_tree(fp, *child, bValues, level+1);
        }
        child = child->next;
    }
}

/*!
 * \param[in] root Root of the subtree to query.
 * \returns true if \p root or any any of its elements require topology
 *   information, false otherwise.
 */
bool
_gmx_selelem_requires_top(const gmx::SelectionTreeElement &root)
{
    if (root.type == SEL_EXPRESSION || root.type == SEL_MODIFIER)
    {
        if (root.u.expr.method && (root.u.expr.method->flags & SMETH_REQTOP))
        {
            return true;
        }
        if (root.u.expr.pc && gmx_ana_poscalc_requires_top(root.u.expr.pc))
        {
            return true;
        }
    }
    gmx::SelectionTreeElementPointer child = root.child;
    while (child)
    {
        if (_gmx_selelem_requires_top(*child))
        {
            return true;
        }
        child = child->next;
    }
    return false;
}
