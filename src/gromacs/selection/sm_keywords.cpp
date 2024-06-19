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
 * Implements internal selection methods for numeric and string keyword
 * evaluation.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include <cctype>
#include <cstring>

#include <memory>
#include <regex>
#include <string>
#include <vector>

#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/selection/selparam.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "keywords.h"
#include "parsetree.h"
#include "scanner.h"
#include "selelem.h"
#include "selmethod.h"

struct gmx_mtop_t;

/*! \brief
 * Allocates data for integer keyword evaluation.
 *
 * \param[in] npar  Not used.
 * \param     param Not used.
 * \returns   Pointer to the allocated data (\ref t_methoddata_kwint).
 *
 * Allocates memory for a \ref t_methoddata_kwint structure.
 */
static void* init_data_kwint(int npar, gmx_ana_selparam_t* param);
/*! \brief
 * Allocates data for real keyword evaluation.
 *
 * \param[in] npar  Not used.
 * \param     param Not used.
 * \returns   Pointer to the allocated data (\ref t_methoddata_kwreal).
 *
 * Allocates memory for a \ref t_methoddata_kwreal structure.
 */
static void* init_data_kwreal(int npar, gmx_ana_selparam_t* param);
/*! \brief
 * Allocates data for string keyword evaluation.
 *
 * \param[in] npar  Not used.
 * \param     param Not used.
 * \returns Pointer to the allocated data (t_methoddata_kwstr).
 *
 * Allocates memory for a t_methoddata_kwstr structure.
 */
static void* init_data_kwstr(int npar, gmx_ana_selparam_t* param);
/** /brief Initializes data for integer keyword evaluation.
 *
 * \param[in] top   Not used.
 * \param[in] npar  Not used (should be 2).
 * \param[in] param Method parameters (should point to \ref smparams_keyword_int).
 * \param[in] data  Should point to \ref t_methoddata_kwint.
 */
static void init_kwint(const gmx_mtop_t* top, int npar, gmx_ana_selparam_t* param, void* data);
/*! \brief
 * Initializes data for real keyword evaluation.
 *
 * \param[in] top   Not used.
 * \param[in] npar  Not used (should be 2).
 * \param[in] param Method parameters (should point to \ref smparams_keyword_real).
 * \param[in] data  Should point to \ref t_methoddata_kwreal.
 * \returns   0 (the initialization always succeeds).
 */
static void init_kwreal(const gmx_mtop_t* top, int npar, gmx_ana_selparam_t* param, void* data);
/*! \brief
 * Initializes data for string keyword evaluation.
 *
 * \param[in] top   Not used.
 * \param[in] npar  Not used (should be 2).
 * \param[in] param Method parameters (should point to \ref smparams_keyword_str).
 * \param[in] data  Should point to t_methoddata_kwstr.
 */
static void init_kwstr(const gmx_mtop_t* top, int npar, gmx_ana_selparam_t* param, void* data);
/** Frees the memory allocated for string keyword evaluation. */
static void free_data_kwstr(void* data);
/** Evaluates integer selection keywords. */
static void evaluate_keyword_int(const gmx::SelMethodEvalContext& context,
                                 gmx_ana_index_t*                 g,
                                 gmx_ana_selvalue_t*              out,
                                 void*                            data);
/** Evaluates real selection keywords. */
static void evaluate_keyword_real(const gmx::SelMethodEvalContext& context,
                                  gmx_ana_index_t*                 g,
                                  gmx_ana_selvalue_t*              out,
                                  void*                            data);
/** Evaluates string selection keywords. */
static void evaluate_keyword_str(const gmx::SelMethodEvalContext& context,
                                 gmx_ana_index_t*                 g,
                                 gmx_ana_selvalue_t*              out,
                                 void*                            data);

/*! \internal \brief
 * Data structure for integer keyword expression evaluation.
 */
typedef struct t_methoddata_kwint
{
    /** Array of values for the keyword. */
    int* v;
    /** Number of ranges in the \p r array. */
    int n;
    /*! \brief
     * Array of sorted integer ranges to match against.
     *
     * Each range is made of two integers, giving the endpoints (inclusive).
     * This field stores the pointer to the ranges allocated by the
     * parameter parser; see \ref SPAR_RANGES for more information.
     */
    int* r;
} t_methoddata_kwint;

/*! \internal \brief
 * Data structure for real keyword expression evaluation.
 */
typedef struct t_methoddata_kwreal
{
    /** Array of values for the keyword. */
    real* v;
    /** Number of ranges in the \p r array. */
    int n;
    /*! \brief
     * Array of sorted ranges to match against.
     *
     * Each range is made of two values, giving the endpoints (inclusive).
     * This field stores the pointer to the ranges allocated by the
     * parameter parser; see \ref SPAR_RANGES for more information.
     */
    real* r;
} t_methoddata_kwreal;

namespace
{

/*! \internal \brief
 * Single item in the list of strings/regular expressions to match.
 *
 * \ingroup module_selection
 */
class StringKeywordMatchItem
{
public:
    /*! \brief
     * Constructs a matcher from a string.
     *
     * \param[in] matchType String matching type.
     * \param[in] str       String to use for matching.
     */
    StringKeywordMatchItem(gmx::SelectionStringMatchType matchType, const char* str) : str_(str)
    {
        useRegExp_ = (matchType == gmx::eStringMatchType_RegularExpression);
        if (matchType == gmx::eStringMatchType_Auto)
        {
            for (size_t j = 0; j < std::strlen(str); ++j)
            {
                if (std::ispunct(str[j]) && str[j] != '?' && str[j] != '*')
                {
                    useRegExp_ = true;
                    break;
                }
            }
        }
        if (useRegExp_)
        {
            try
            {
                regex_.assign(str, std::regex::nosubs | std::regex::extended);
            }
            catch (const std::regex_error& /*ex*/)
            {
                // TODO: Better error messages.
                GMX_THROW(gmx::InvalidInputError(
                        gmx::formatString("Error in regular expression \"%s\"", str)));
            }
        }
    }

    /*! \brief
     * Checks whether this item matches a string.
     *
     * \param[in] matchType String matching type.
     * \param[in] value     String to match.
     * \returns   true if this item matches \p value.
     */
    bool match(gmx::SelectionStringMatchType matchType, const char* value) const
    {
        if (matchType == gmx::eStringMatchType_Exact)
        {
            return str_ == value;
        }
        else if (useRegExp_)
        {
            return std::regex_match(value, regex_);
        }
        else
        {
            return gmx_wcmatch(str_.c_str(), value) == 0;
        }
    }

private:
    //! The raw string passed for the matcher.
    std::string str_;
    //! Whether a regular expression match is used.
    bool useRegExp_;
    //! Regular expression compiled from \p str_, if applicable.
    std::regex regex_;
};

/*! \internal \brief
 * Data structure for string keyword expression evaluation.
 */
struct t_methoddata_kwstr
{
    /** Matching type for the strings. */
    gmx::SelectionStringMatchType matchType;
    /** Array of values for the keyword. */
    char** v;
    /** Array of strings/regular expressions to match against.*/
    std::vector<StringKeywordMatchItem> matches;
};

} // namespace

/** Parameters for integer keyword evaluation. */
static gmx_ana_selparam_t smparams_keyword_int[] = {
    { nullptr, { INT_VALUE, -1, { nullptr } }, nullptr, SPAR_ATOMVAL },
    { nullptr, { INT_VALUE, -1, { nullptr } }, nullptr, SPAR_RANGES | SPAR_VARNUM },
};

/** Parameters for real keyword evaluation. */
static gmx_ana_selparam_t smparams_keyword_real[] = {
    { nullptr, { REAL_VALUE, -1, { nullptr } }, nullptr, SPAR_ATOMVAL | SPAR_DYNAMIC },
    { nullptr, { REAL_VALUE, -1, { nullptr } }, nullptr, SPAR_RANGES | SPAR_VARNUM },
};

/** Parameters for string keyword evaluation. */
static gmx_ana_selparam_t smparams_keyword_str[] = {
    { nullptr, { STR_VALUE, -1, { nullptr } }, nullptr, SPAR_ATOMVAL },
    { nullptr, { STR_VALUE, -1, { nullptr } }, nullptr, SPAR_VARNUM },
};

/** Selection method data for integer keyword evaluation. */
gmx_ana_selmethod_t sm_keyword_int = {
    "kw_int",
    GROUP_VALUE,
    SMETH_SINGLEVAL,
    asize(smparams_keyword_int),
    smparams_keyword_int,
    &init_data_kwint,
    nullptr,
    &init_kwint,
    nullptr,
    nullptr,
    nullptr,
    &evaluate_keyword_int,
    nullptr,
    { nullptr, nullptr, 0, nullptr },
};

/** Selection method data for real keyword evaluation. */
gmx_ana_selmethod_t sm_keyword_real = {
    "kw_real",
    GROUP_VALUE,
    SMETH_SINGLEVAL,
    asize(smparams_keyword_real),
    smparams_keyword_real,
    &init_data_kwreal,
    nullptr,
    &init_kwreal,
    nullptr,
    nullptr,
    nullptr,
    &evaluate_keyword_real,
    nullptr,
    { nullptr, nullptr, 0, nullptr },
};

/** Selection method data for string keyword evaluation. */
gmx_ana_selmethod_t sm_keyword_str = {
    "kw_str",
    GROUP_VALUE,
    SMETH_SINGLEVAL,
    asize(smparams_keyword_str),
    smparams_keyword_str,
    &init_data_kwstr,
    nullptr,
    &init_kwstr,
    nullptr,
    &free_data_kwstr,
    nullptr,
    &evaluate_keyword_str,
    nullptr,
    { nullptr, nullptr, 0, nullptr },
};

/*! \brief
 * Initializes keyword evaluation for an arbitrary group.
 *
 * \param[in] top   Not used.
 * \param[in] npar  Not used.
 * \param[in] param Not used.
 * \param[in] data  Should point to \ref t_methoddata_kweval.
 * \returns   0 on success, a non-zero error code on return.
 *
 * Calls the initialization method of the wrapped keyword.
 */
static void init_kweval(const gmx_mtop_t* top, int npar, gmx_ana_selparam_t* param, void* data);
/*! \brief
 * Initializes output for keyword evaluation in an arbitrary group.
 *
 * \param[in]     top   Not used.
 * \param[in,out] out   Pointer to output data structure.
 * \param[in,out] data  Should point to \c t_methoddata_kweval.
 * \returns       0 for success.
 */
static void init_output_kweval(const gmx_mtop_t* top, gmx_ana_selvalue_t* out, void* data);
/** Frees the data allocated for keyword evaluation in an arbitrary group. */
static void free_data_kweval(void* data);
/** Initializes frame evaluation for keyword evaluation in an arbitrary group. */
static void init_frame_kweval(const gmx::SelMethodEvalContext& context, void* data);
/*! \brief
 * Evaluates keywords in an arbitrary group.
 *
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_kweval.
 *
 * Calls the evaluation function of the wrapped keyword with the given
 * parameters, with the exception of using \c t_methoddata_kweval::g for the
 * evaluation group.
 */
static void evaluate_kweval(const gmx::SelMethodEvalContext& context,
                            gmx_ana_index_t*                 g,
                            gmx_ana_selvalue_t*              out,
                            void*                            data);
/*! \brief
 * Evaluates keywords in an arbitrary set of positions.
 *
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_kweval.
 *
 * Calls the evaluation function of the wrapped keyword with the given
 * parameters, with the exception of using \c t_methoddata_kweval::p for the
 * evaluation positions.
 */
static void evaluate_kweval_pos(const gmx::SelMethodEvalContext& context,
                                gmx_ana_index_t*                 g,
                                gmx_ana_selvalue_t*              out,
                                void*                            data);

/*! \internal \brief
 * Data structure for keyword evaluation in arbitrary groups.
 */
struct t_methoddata_kweval
{
    //! Initialize new keyword evaluator for the given keyword.
    t_methoddata_kweval(gmx_ana_selmethod_t* method, void* data) : kwmethod(method), kwmdata(data)
    {
        gmx_ana_index_clear(&g);
    }
    ~t_methoddata_kweval() { _gmx_selelem_free_method(kwmethod, kwmdata); }

    //! Wrapped keyword method for evaluating the values.
    gmx_ana_selmethod_t* kwmethod;
    //! Method data for \p kwmethod.
    void* kwmdata;
    //! Group in which \p kwmethod should be evaluated.
    gmx_ana_index_t g;
    //! Positions for which \p kwmethod should be evaluated.
    gmx_ana_pos_t p;
};

/** Parameters for keyword evaluation in an arbitrary group. */
static gmx_ana_selparam_t smparams_kweval_group[] = {
    { nullptr, { GROUP_VALUE, 1, { nullptr } }, nullptr, SPAR_DYNAMIC },
};
/** Parameters for keyword evaluation from positions. */
static gmx_ana_selparam_t smparams_kweval_pos[] = {
    { nullptr, { POS_VALUE, 1, { nullptr } }, nullptr, SPAR_DYNAMIC },
};


/********************************************************************
 * INTEGER KEYWORD EVALUATION
 ********************************************************************/

static void* init_data_kwint(int /* npar */, gmx_ana_selparam_t* /* param */)
{
    t_methoddata_kwint* data;

    snew(data, 1);
    return data;
}

static void init_kwint(const gmx_mtop_t* /* top */, int /* npar */, gmx_ana_selparam_t* param, void* data)
{
    t_methoddata_kwint* d = static_cast<t_methoddata_kwint*>(data);

    d->v = param[0].val.u.i;
    d->n = param[1].val.nr;
    d->r = param[1].val.u.i;
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_kwint.
 *
 * Does a binary search to find which atoms match the ranges in the
 * \c t_methoddata_kwint structure for this selection.
 * Matching atoms are stored in \p out->u.g.
 */
static void evaluate_keyword_int(const gmx::SelMethodEvalContext& /*context*/,
                                 gmx_ana_index_t*    g,
                                 gmx_ana_selvalue_t* out,
                                 void*               data)
{
    t_methoddata_kwint* d = static_cast<t_methoddata_kwint*>(data);
    int                 n, i, j, jmin, jmax;
    int                 val;

    out->u.g->isize = 0;
    n               = d->n;
    for (i = 0; i < g->isize; ++i)
    {
        val = d->v[i];
        if (d->r[0] > val || d->r[2 * n - 1] < val)
        {
            continue;
        }
        jmin = 0;
        jmax = n;
        while (jmax - jmin > 1)
        {
            j = jmin + (jmax - jmin) / 2;
            if (val < d->r[2 * j])
            {
                jmax = j;
            }
            else
            {
                jmin = j;
                if (val <= d->r[2 * j + 1])
                {
                    break;
                }
                /* ++jmin;*/
            }
        }
        if (val <= d->r[2 * jmin + 1])
        {
            out->u.g->index[out->u.g->isize++] = g->index[i];
        }
    }
}


/********************************************************************
 * REAL KEYWORD EVALUATION
 ********************************************************************/

static void* init_data_kwreal(int /* npar */, gmx_ana_selparam_t* /* param */)
{
    t_methoddata_kwreal* data;

    snew(data, 1);
    return data;
}

static void init_kwreal(const gmx_mtop_t* /* top */, int /* npar */, gmx_ana_selparam_t* param, void* data)
{
    t_methoddata_kwreal* d = static_cast<t_methoddata_kwreal*>(data);

    d->v = param[0].val.u.r;
    d->n = param[1].val.nr;
    d->r = param[1].val.u.r;
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_kwreal.
 *
 * Does a binary search to find which atoms match the ranges in the
 * \c t_methoddata_kwreal structure for this selection.
 * Matching atoms are stored in \p out->u.g.
 */
static void evaluate_keyword_real(const gmx::SelMethodEvalContext& /*context*/,
                                  gmx_ana_index_t*    g,
                                  gmx_ana_selvalue_t* out,
                                  void*               data)
{
    t_methoddata_kwreal* d = static_cast<t_methoddata_kwreal*>(data);
    int                  n, i, j, jmin, jmax;
    real                 val;

    out->u.g->isize = 0;
    n               = d->n;
    for (i = 0; i < g->isize; ++i)
    {
        val = d->v[i];
        if (d->r[0] > val || d->r[2 * n - 1] < val)
        {
            continue;
        }
        jmin = 0;
        jmax = n;
        while (jmax - jmin > 1)
        {
            j = jmin + (jmax - jmin) / 2;
            if (val < d->r[2 * j])
            {
                jmax = j;
            }
            else
            {
                jmin = j;
                if (val <= d->r[2 * j + 1])
                {
                    break;
                }
                /* ++jmin;*/
            }
        }
        if (val <= d->r[2 * jmin + 1])
        {
            out->u.g->index[out->u.g->isize++] = g->index[i];
        }
    }
}


/********************************************************************
 * STRING KEYWORD EVALUATION
 ********************************************************************/

static void* init_data_kwstr(int /* npar */, gmx_ana_selparam_t* /* param */)
{
    t_methoddata_kwstr* data = new t_methoddata_kwstr();
    data->matchType          = gmx::eStringMatchType_Auto;
    return data;
}

/*!
 * \param[in,out] sel   Selection element to initialize.
 * \param[in]     matchType  Method to use to match string values.
 *
 * Sets the string matching method for string keyword matching.
 */
void _gmx_selelem_set_kwstr_match_type(const gmx::SelectionTreeElementPointer& sel,
                                       gmx::SelectionStringMatchType           matchType)
{
    t_methoddata_kwstr* d = static_cast<t_methoddata_kwstr*>(sel->u.expr.mdata);

    if (sel->type != SEL_EXPRESSION || !sel->u.expr.method
        || sel->u.expr.method->name != sm_keyword_str.name)
    {
        return;
    }
    d->matchType = matchType;
}

static void init_kwstr(const gmx_mtop_t* /* top */, int /* npar */, gmx_ana_selparam_t* param, void* data)
{
    t_methoddata_kwstr* d = static_cast<t_methoddata_kwstr*>(data);

    d->v = param[0].val.u.s;
    /* Return if this is not the first time */
    if (!d->matches.empty())
    {
        return;
    }
    int n = param[1].val.nr;
    d->matches.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        const char* s = param[1].val.u.s[i];
        d->matches.emplace_back(d->matchType, s);
    }
}

/*!
 * \param data Data to free (should point to a t_methoddata_kwstr).
 */
static void free_data_kwstr(void* data)
{
    t_methoddata_kwstr* d = static_cast<t_methoddata_kwstr*>(data);
    delete d;
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_kwstr.
 *
 * Does a linear search to find which atoms match the strings in the
 * \c t_methoddata_kwstr structure for this selection.
 * Wildcards are allowed in the strings.
 * Matching atoms are stored in \p out->u.g.
 */
static void evaluate_keyword_str(const gmx::SelMethodEvalContext& /*context*/,
                                 gmx_ana_index_t*    g,
                                 gmx_ana_selvalue_t* out,
                                 void*               data)
{
    t_methoddata_kwstr* d = static_cast<t_methoddata_kwstr*>(data);

    out->u.g->isize = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        for (size_t j = 0; j < d->matches.size(); ++j)
        {
            if (d->matches[j].match(d->matchType, d->v[i]))
            {
                out->u.g->index[out->u.g->isize++] = g->index[i];
                break;
            }
        }
    }
}


/********************************************************************
 * KEYWORD EVALUATION FOR ARBITRARY GROUPS
 ********************************************************************/

static void init_kweval(const gmx_mtop_t* top, int /* npar */, gmx_ana_selparam_t* /* param */, void* data)
{
    t_methoddata_kweval* d = static_cast<t_methoddata_kweval*>(data);

    d->kwmethod->init(top, 0, nullptr, d->kwmdata);
}

static void init_output_kweval(const gmx_mtop_t* /* top */, gmx_ana_selvalue_t* out, void* data)
{
    t_methoddata_kweval* d = static_cast<t_methoddata_kweval*>(data);

    out->nr = d->g.isize;
    _gmx_selvalue_reserve(out, out->nr);
}

/*!
 * \param data Data to free (should point to a \c t_methoddata_kweval).
 *
 * Frees the memory allocated for all the members of \c t_methoddata_kweval.
 */
static void free_data_kweval(void* data)
{
    t_methoddata_kweval* d = static_cast<t_methoddata_kweval*>(data);

    delete d;
}

/*!
 * \param[in]  context  Evaluation context.
 * \param      data Should point to a \ref t_methoddata_kweval.
 */
static void init_frame_kweval(const gmx::SelMethodEvalContext& context, void* data)
{
    t_methoddata_kweval* d = static_cast<t_methoddata_kweval*>(data);

    d->kwmethod->init_frame(context, d->kwmdata);
}

static void evaluate_kweval(const gmx::SelMethodEvalContext& context,
                            gmx_ana_index_t* /* g */,
                            gmx_ana_selvalue_t* out,
                            void*               data)
{
    t_methoddata_kweval* d = static_cast<t_methoddata_kweval*>(data);

    d->kwmethod->update(context, &d->g, out, d->kwmdata);
}

static void evaluate_kweval_pos(const gmx::SelMethodEvalContext& context,
                                gmx_ana_index_t* /* g */,
                                gmx_ana_selvalue_t* out,
                                void*               data)
{
    t_methoddata_kweval* d = static_cast<t_methoddata_kweval*>(data);

    d->kwmethod->pupdate(context, &d->p, out, d->kwmdata);
}

/*! \brief
 * Initializes a selection element for evaluating a keyword in a given group.
 *
 * \param[in]   method  Keyword selection method to evaluate.
 * \param[in]   params  Parameters to pass to initialization (the child group).
 * \param[in]   scanner Scanner data structure.
 * \returns     Pointer to the created selection element.
 *
 * Implements _gmx_sel_init_keyword_evaluator() for \ref GROUP_VALUE input
 * selections.
 */
static gmx::SelectionTreeElementPointer init_evaluator_group(gmx_ana_selmethod_t* method,
                                                             const gmx::SelectionParserParameterList& params,
                                                             void* scanner)
{
    if ((method->flags & (SMETH_SINGLEVAL | SMETH_VARNUMVAL)) || method->outinit || method->pupdate)
    {
        std::string message =
                gmx::formatString("Keyword '%s' cannot be evaluated in this context", method->name);
        GMX_THROW(gmx::InvalidInputError(message));
    }

    // TODO: For same ... as ..., some other location could be more intuitive.
    gmx::SelectionTreeElementPointer sel(new gmx::SelectionTreeElement(
            SEL_EXPRESSION, _gmx_sel_lexer_get_current_location(scanner)));
    _gmx_selelem_set_method(sel, method, scanner);

    t_methoddata_kweval* data = new t_methoddata_kweval(sel->u.expr.method, sel->u.expr.mdata);

    snew(sel->u.expr.method, 1);
    sel->u.expr.method->name        = data->kwmethod->name;
    sel->u.expr.method->type        = data->kwmethod->type;
    sel->u.expr.method->flags       = data->kwmethod->flags | SMETH_VARNUMVAL;
    sel->u.expr.method->init_data   = nullptr;
    sel->u.expr.method->set_poscoll = nullptr;
    sel->u.expr.method->init        = method->init ? &init_kweval : nullptr;
    sel->u.expr.method->outinit     = &init_output_kweval;
    sel->u.expr.method->free        = &free_data_kweval;
    sel->u.expr.method->init_frame  = method->init_frame ? &init_frame_kweval : nullptr;
    sel->u.expr.method->update      = &evaluate_kweval;
    sel->u.expr.method->pupdate     = nullptr;
    sel->u.expr.method->nparams     = asize(smparams_kweval_group);
    sel->u.expr.method->param       = smparams_kweval_group;
    _gmx_selelem_init_method_params(sel, scanner);
    sel->u.expr.mdata = data;

    sel->u.expr.method->param[0].val.u.g = &data->g;

    _gmx_sel_parse_params(params, sel->u.expr.method->nparams, sel->u.expr.method->param, sel, scanner);
    return sel;
}

/*! \brief
 * Initializes a selection element for evaluating a keyword in given positions.
 *
 * \param[in]   method  Keyword selection method to evaluate.
 * \param[in]   params  Parameters to pass to initialization (the child positions).
 * \param[in]   scanner Scanner data structure.
 * \returns     Pointer to the created selection element.
 *
 * Implements _gmx_sel_init_keyword_evaluator() for \ref POS_VALUE input
 * selections.
 */
static gmx::SelectionTreeElementPointer init_evaluator_pos(gmx_ana_selmethod_t* method,
                                                           const gmx::SelectionParserParameterList& params,
                                                           void* scanner)
{
    if ((method->flags & (SMETH_SINGLEVAL | SMETH_VARNUMVAL)) || method->outinit || method->pupdate == nullptr)
    {
        std::string message =
                gmx::formatString("Keyword '%s' cannot be evaluated in this context", method->name);
        GMX_THROW(gmx::InvalidInputError(message));
    }

    gmx::SelectionTreeElementPointer sel(new gmx::SelectionTreeElement(
            SEL_EXPRESSION, _gmx_sel_lexer_get_current_location(scanner)));
    _gmx_selelem_set_method(sel, method, scanner);

    t_methoddata_kweval* data = new t_methoddata_kweval(sel->u.expr.method, sel->u.expr.mdata);

    snew(sel->u.expr.method, 1);
    sel->u.expr.method->name        = data->kwmethod->name;
    sel->u.expr.method->type        = data->kwmethod->type;
    sel->u.expr.method->flags       = data->kwmethod->flags | SMETH_SINGLEVAL;
    sel->u.expr.method->init_data   = nullptr;
    sel->u.expr.method->set_poscoll = nullptr;
    sel->u.expr.method->init        = method->init ? &init_kweval : nullptr;
    sel->u.expr.method->outinit     = nullptr;
    sel->u.expr.method->free        = &free_data_kweval;
    sel->u.expr.method->init_frame  = method->init_frame ? &init_frame_kweval : nullptr;
    sel->u.expr.method->update      = &evaluate_kweval_pos;
    sel->u.expr.method->pupdate     = nullptr;
    sel->u.expr.method->nparams     = asize(smparams_kweval_pos);
    sel->u.expr.method->param       = smparams_kweval_pos;
    _gmx_selelem_init_method_params(sel, scanner);
    sel->u.expr.mdata = data;

    sel->u.expr.method->param[0].val.u.p = &data->p;

    _gmx_sel_parse_params(params, sel->u.expr.method->nparams, sel->u.expr.method->param, sel, scanner);
    return sel;
}

gmx::SelectionTreeElementPointer _gmx_sel_init_keyword_evaluator(gmx_ana_selmethod_t* method,
                                                                 const gmx::SelectionTreeElementPointer& child,
                                                                 void* scanner)
{
    gmx::SelectionParserParameterList params;
    params.push_back(gmx::SelectionParserParameter::createFromExpression(nullptr, child));
    if (child->v.type == GROUP_VALUE)
    {
        return init_evaluator_group(method, params, scanner);
    }
    else if (child->v.type == POS_VALUE)
    {
        return init_evaluator_pos(method, params, scanner);
    }
    else
    {
        std::string text(_gmx_sel_lexer_get_text(scanner, child->location()));
        std::string message =
                gmx::formatString("Expression '%s' cannot be used to evaluate keywords", text.c_str());
        GMX_THROW(gmx::InvalidInputError(message));
    }
}
