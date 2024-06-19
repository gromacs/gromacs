/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements functions in selelem.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "selelem.h"

#include <cstring>

#include <filesystem>

#include "gromacs/math/vectypes.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/selection/selectionenums.h"
#include "gromacs/selection/selparam.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "keywords.h"
#include "mempool.h"
#include "poscalc.h"
#include "selmethod.h"

struct gmx_ana_indexgrps_t;

/*!
 * \param[in] sel Selection for which the string is requested
 * \returns   Pointer to a string that corresponds to \p sel->type.
 *
 * The return value points to a string constant and should not be \p free'd.
 *
 * The function returns NULL if \p sel->type is not one of the valid values.
 */
const char* _gmx_selelem_type_str(const gmx::SelectionTreeElement& sel)
{
    const char* p = nullptr;
    switch (sel.type)
    {
        case SEL_CONST: p = "CONST"; break;
        case SEL_EXPRESSION: p = "EXPR"; break;
        case SEL_BOOLEAN: p = "BOOL"; break;
        case SEL_ARITHMETIC: p = "ARITH"; break;
        case SEL_ROOT: p = "ROOT"; break;
        case SEL_SUBEXPR: p = "SUBEXPR"; break;
        case SEL_SUBEXPRREF: p = "REF"; break;
        case SEL_GROUPREF: p = "GROUPREF"; break;
        case SEL_MODIFIER:
            p = "MODIFIER";
            break;
            // No default clause so we intentionally get compiler errors
            // if new selection choices are added later.
    }
    return p;
}

/*!
 * \param[in] val Value structore for which the string is requested.
 * \returns   Pointer to a string that corresponds to \p val->type,
 *   NULL if the type value is invalid.
 *
 * The return value points to a string constant and should not be \p free'd.
 */
const char* _gmx_sel_value_type_str(const gmx_ana_selvalue_t* val)
{
    const char* p = nullptr;
    switch (val->type)
    {
        case NO_VALUE: p = "NONE"; break;
        case INT_VALUE: p = "INT"; break;
        case REAL_VALUE: p = "REAL"; break;
        case STR_VALUE: p = "STR"; break;
        case POS_VALUE: p = "VEC"; break;
        case GROUP_VALUE:
            p = "GROUP";
            break;
            // No default clause so we intentionally get compiler errors
            // if new selection choices are added later.
    }
    return p;
}

/*! \copydoc _gmx_selelem_type_str() */
const char* _gmx_selelem_boolean_type_str(const gmx::SelectionTreeElement& sel)
{
    const char* p = nullptr;
    switch (sel.u.boolt)
    {
        case BOOL_NOT: p = "NOT"; break;
        case BOOL_AND: p = "AND"; break;
        case BOOL_OR: p = "OR"; break;
        case BOOL_XOR:
            p = "XOR";
            break;
            // No default clause so we intentionally get compiler errors
            // if new selection choices are added later.
    }
    return p;
}


namespace gmx
{

SelectionTreeElement::SelectionTreeElement(e_selelem_t elemType, const SelectionLocation& location) :
    location_(location)
{
    this->type  = elemType;
    this->flags = (elemType != SEL_ROOT) ? SEL_ALLOCVAL : 0;
    if (elemType == SEL_BOOLEAN)
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
    this->evaluate = nullptr;
    this->mempool  = nullptr;
    this->cdata    = nullptr;
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
            case GROUP_VALUE:
                for (int i = 0; i < n; ++i)
                {
                    gmx_ana_index_deinit(&v.u.g[i]);
                }
                break;
            default: /* No special handling for other types */ break;
        }
    }
    _gmx_selvalue_free(&v);
    if (type == SEL_SUBEXPRREF && u.param != nullptr)
    {
        // TODO: This is now called from two different locations.
        // It is likely that one of them is unnecessary, but that requires
        // extra analysis to clarify.
        _gmx_selelem_free_param(u.param);
    }
}

void SelectionTreeElement::freeExpressionData()
{
    if (type == SEL_EXPRESSION || type == SEL_MODIFIER)
    {
        _gmx_selelem_free_method(u.expr.method, u.expr.mdata);
        u.expr.mdata  = nullptr;
        u.expr.method = nullptr;
        /* Free position data */
        delete u.expr.pos;
        u.expr.pos = nullptr;
        /* Free position calculation data */
        if (u.expr.pc)
        {
            gmx_ana_poscalc_free(u.expr.pc);
            u.expr.pc = nullptr;
        }
    }
    if (type == SEL_SUBEXPR || type == SEL_ROOT || (type == SEL_CONST && v.type == GROUP_VALUE))
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
            v.u.i = static_cast<int*>(_gmx_sel_mempool_alloc(mempool, sizeof(*v.u.i) * count));
            break;

        case REAL_VALUE:
            v.u.r = static_cast<real*>(_gmx_sel_mempool_alloc(mempool, sizeof(*v.u.r) * count));
            break;

        case GROUP_VALUE: _gmx_sel_mempool_alloc_group(mempool, v.u.g, count); break;

        default: gmx_incons("Memory pooling not implemented for requested type");
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
            _gmx_selvalue_setstore(&v, nullptr);
            break;

        case GROUP_VALUE:
            if (v.u.g)
            {
                _gmx_sel_mempool_free_group(mempool, v.u.g);
            }
            break;

        default: gmx_incons("Memory pooling not implemented for requested type");
    }
}

void SelectionTreeElement::fillNameIfMissing(const char* selectionText)
{
    GMX_RELEASE_ASSERT(type == SEL_ROOT, "Should not be called for non-root elements");
    if (name().empty())
    {
        // Check whether the actual selection given was from an external group,
        // and if so, use the name of the external group.
        SelectionTreeElementPointer childElem = this->child;
        if (_gmx_selelem_is_default_kwpos(*childElem) && childElem->child
            && childElem->child->type == SEL_SUBEXPRREF && childElem->child->child)
        {
            if (childElem->child->child->type == SEL_CONST && childElem->child->child->v.type == GROUP_VALUE)
            {
                setName(childElem->child->child->name());
                return;
            }
            // If the group reference is still unresolved, leave the name empty
            // and fill it later.
            if (childElem->child->child->type == SEL_GROUPREF)
            {
                return;
            }
        }
        // If there still is no name, use the selection string.
        setName(selectionText);
    }
}

SelectionTopologyProperties SelectionTreeElement::requiredTopologyProperties() const
{
    SelectionTopologyProperties props;
    if (type == SEL_EXPRESSION || type == SEL_MODIFIER)
    {
        bool needsTop    = false;
        bool needsMasses = false;
        if (u.expr.method != nullptr)
        {
            needsTop    = ((u.expr.method->flags & SMETH_REQTOP) != 0);
            needsMasses = ((u.expr.method->flags & SMETH_REQMASS) != 0);
        }
        if (u.expr.pc != nullptr)
        {
            auto requiredTopologyInfo = gmx_ana_poscalc_required_topology_info(u.expr.pc);
            needsTop                  = needsTop
                       || (requiredTopologyInfo
                           != PositionCalculationCollection::RequiredTopologyInfo::None);
            needsMasses = needsMasses
                          || (requiredTopologyInfo
                              == PositionCalculationCollection::RequiredTopologyInfo::TopologyAndMasses);
        }
        if (needsTop)
        {
            props.merge(SelectionTopologyProperties::topology());
        }
        if (needsMasses)
        {
            props.merge(SelectionTopologyProperties::masses());
        }
    }
    SelectionTreeElementPointer childElem = this->child;
    while (childElem && !props.hasAll())
    {
        props.merge(childElem->requiredTopologyProperties());
        childElem = childElem->next;
    }
    return props;
}

void SelectionTreeElement::checkUnsortedAtoms(bool bUnsortedAllowed, ExceptionInitializer* errors) const
{
    const bool bUnsortedSupported =
            (type == SEL_CONST && v.type == GROUP_VALUE) || type == SEL_ROOT || type == SEL_SUBEXPR
            || type == SEL_SUBEXPRREF
            // TODO: Consolidate.
            || type == SEL_MODIFIER
            || (type == SEL_EXPRESSION && ((u.expr.method->flags & SMETH_ALLOW_UNSORTED) != 0));

    // TODO: For some complicated selections, this may result in the same
    // index group reference being flagged as an error multiple times for the
    // same selection.
    SelectionTreeElementPointer childElem = this->child;
    while (childElem)
    {
        childElem->checkUnsortedAtoms(bUnsortedAllowed && bUnsortedSupported, errors);
        childElem = childElem->next;
    }

    // The logic here is simplified by the fact that only constant groups can
    // currently be the root cause of SEL_UNSORTED being set, so only those
    // need to be considered in triggering the error.
    if (!bUnsortedAllowed && (flags & SEL_UNSORTED) && type == SEL_CONST && v.type == GROUP_VALUE)
    {
        std::string message = formatString(
                "Group '%s' cannot be used in selections except "
                "as a full value of the selection, "
                "because atom indices in it are not sorted and/or "
                "it contains duplicate atoms.",
                name().c_str());
        errors->addNested(InconsistentInputError(message));
    }
}

bool SelectionTreeElement::requiresIndexGroups() const
{
    if (type == SEL_GROUPREF)
    {
        return true;
    }
    SelectionTreeElementPointer childElem = this->child;
    while (childElem)
    {
        if (childElem->requiresIndexGroups())
        {
            return true;
        }
        childElem = childElem->next;
    }
    return false;
}

void SelectionTreeElement::resolveIndexGroupReference(gmx_ana_indexgrps_t* grps, int natoms)
{
    GMX_RELEASE_ASSERT(type == SEL_GROUPREF,
                       "Should only be called for index group reference elements");
    if (grps == nullptr)
    {
        std::string message = formatString(
                "Cannot match '%s', because index groups are not available.", name().c_str());
        GMX_THROW(InconsistentInputError(message));
    }

    gmx_ana_index_t foundGroup;
    std::string     foundName;
    if (u.gref.name != nullptr)
    {
        if (!gmx_ana_indexgrps_find(&foundGroup, &foundName, grps, u.gref.name))
        {
            std::string message = formatString(
                    "Cannot match '%s', because no such index group can be found.", name().c_str());
            GMX_THROW(InconsistentInputError(message));
        }
    }
    else
    {
        if (!gmx_ana_indexgrps_extract(&foundGroup, &foundName, grps, u.gref.id))
        {
            std::string message = formatString(
                    "Cannot match '%s', because no such index group can be found.", name().c_str());
            GMX_THROW(InconsistentInputError(message));
        }
    }

    if (!gmx_ana_index_check_sorted(&foundGroup))
    {
        flags |= SEL_UNSORTED;
    }

    sfree(u.gref.name);
    type = SEL_CONST;
    gmx_ana_index_set(&u.cgrp, foundGroup.isize, foundGroup.index, foundGroup.nalloc_index);
    setName(foundName);

    if (natoms > 0)
    {
        checkIndexGroup(natoms);
    }
}

void SelectionTreeElement::checkIndexGroup(int natoms)
{
    GMX_RELEASE_ASSERT(type == SEL_CONST && v.type == GROUP_VALUE,
                       "Should only be called for index group elements");
    if (!gmx_ana_index_check_range(&u.cgrp, natoms))
    {
        std::string message = formatString(
                "Group '%s' cannot be used in selections, because it "
                "contains negative atom indices and/or references atoms "
                "not present (largest allowed atom index is %d).",
                name().c_str(),
                natoms);
        GMX_THROW(InconsistentInputError(message));
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
void _gmx_selelem_set_vtype(const gmx::SelectionTreeElementPointer& sel, e_selvalue_t vtype)
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

void _gmx_selelem_free_param(gmx_ana_selparam_t* param)
{
    if (param->val.u.ptr != nullptr)
    {
        if (param->val.type == GROUP_VALUE)
        {
            for (int i = 0; i < param->val.nr; ++i)
            {
                gmx_ana_index_deinit(&param->val.u.g[i]);
            }
        }
        _gmx_selvalue_free(&param->val);
    }
}

void _gmx_selelem_free_method(gmx_ana_selmethod_t* method, void* mdata)
{
    sel_freefunc free_func = nullptr;

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
        /* Free the memory allocated for the parameters that are not managed
         * by the selection method itself. */
        for (int i = 0; i < method->nparams; ++i)
        {
            _gmx_selelem_free_param(&method->param[i]);
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
void _gmx_selelem_print_tree(FILE* fp, const gmx::SelectionTreeElement& sel, bool bValues, int level)
{
    int i;

    fprintf(fp, "%*c %s %s", level * 2 + 1, '*', _gmx_selelem_type_str(sel), _gmx_sel_value_type_str(&sel.v));
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
    if (sel.flags & SEL_ALLOCVAL)
    {
        fprintf(fp, "Av");
    }
    if (sel.flags & SEL_ALLOCDATA)
    {
        fprintf(fp, "Ad");
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
            const gmx_ana_index_t* g = sel.v.u.g;
            if (!g || g->isize == 0)
            {
                g = &sel.u.cgrp;
            }
            fprintf(fp, " (%d atoms)", g->isize);
        }
    }
    else if (sel.type == SEL_BOOLEAN)
    {
        fprintf(fp, " %s", _gmx_selelem_boolean_type_str(sel));
    }
    else if (sel.type == SEL_EXPRESSION && sel.u.expr.method->name == sm_compare.name)
    {
        _gmx_selelem_print_compare_info(fp, sel.u.expr.mdata);
    }
    if (sel.evaluate)
    {
        fprintf(fp, " eval=");
        _gmx_sel_print_evalfunc_name(fp, sel.evaluate);
    }
    if (sel.v.nalloc < 0)
    {
        fprintf(fp, " (ext)");
    }
    fprintf(fp, "\n");

    if ((sel.type == SEL_CONST && sel.v.type == GROUP_VALUE) || sel.type == SEL_ROOT)
    {
        const gmx_ana_index_t* g = sel.v.u.g;
        if (!g || g->isize == 0 || sel.evaluate != nullptr)
        {
            g = &sel.u.cgrp;
        }
        if (g->isize < 0)
        {
            fprintf(fp, "%*c group: (null)\n", level * 2 + 1, ' ');
        }
        else if (g->isize > 0)
        {
            fprintf(fp, "%*c group:", level * 2 + 1, ' ');
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
            fprintf(fp, "%*c COM", level * 2 + 3, '*');
            fprintf(fp, "\n");
        }
    }
    else if (sel.type == SEL_SUBEXPRREF && sel.u.param != nullptr)
    {
        fprintf(fp, "%*c param", level * 2 + 1, ' ');
        if (sel.u.param->name != nullptr)
        {
            fprintf(fp, " \"%s\"", sel.u.param->name);
        }
        if (sel.u.param->val.nalloc < 0)
        {
            fprintf(fp, " (ext)");
        }
        else
        {
            fprintf(fp, " nalloc: %d", sel.u.param->val.nalloc);
        }
        fprintf(fp, "\n");
    }

    if (sel.cdata)
    {
        _gmx_selelem_print_compiler_info(fp, sel, level);
    }

    if (bValues && sel.type != SEL_CONST && sel.type != SEL_ROOT && sel.v.u.ptr)
    {
        fprintf(fp, "%*c value: ", level * 2 + 1, ' ');
        switch (sel.v.type)
        {
            case POS_VALUE:
                /* In normal use, the pointer should never be NULL, but it's
                 * useful to have the check for debugging to avoid accidental
                 * segfaults when printing the selection tree. */
                if (sel.v.u.p->x)
                {
                    fprintf(fp, "(%f, %f, %f)", sel.v.u.p->x[0][XX], sel.v.u.p->x[0][YY], sel.v.u.p->x[0][ZZ]);
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
            default: fprintf(fp, "???"); break;
        }
        fprintf(fp, "\n");
    }

    /* Print the subexpressions with one more level of indentation */
    gmx::SelectionTreeElementPointer child = sel.child;
    while (child)
    {
        if (!(sel.type == SEL_SUBEXPRREF && child->type == SEL_SUBEXPR))
        {
            _gmx_selelem_print_tree(fp, *child, bValues, level + 1);
        }
        child = child->next;
    }
}
