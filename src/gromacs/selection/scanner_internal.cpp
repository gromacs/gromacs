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
 * \brief Helper functions for the selection tokenizer.
 *
 * This file implements the functions in the headers scanner.h and
 * scanner_internal.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
/*! \cond
 * \internal \file scanner_flex.h
 * \brief Generated (from scanner.l) header file by Flex.
 *
 * This file contains definitions of functions that are needed in
 * scanner_internal.cpp.
 *
 * \ingroup module_selection
 * \endcond
 */
#include "gmxpre.h"

#include "scanner_internal.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <memory>
#include <string>

#include "gromacs/selection/scanner_flex.h"
#include "gromacs/selection/selparam.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "parser.h"
#include "parsetree.h"
#include "scanner.h"
#include "selectioncollection_impl.h"
#include "selelem.h"
#include "selmethod.h"
#include "symrec.h"

/* These are defined as macros in the generated scanner_flex.h.
 * We undefine them here to have them as variable names in the subroutines.
 * There are other ways of doing this, but this is probably the easiest. */
#undef yylval
#undef yytext
#undef yyleng

/*! \brief
 * Handles initialization of method parameter token.
 */
static int init_param_token(YYSTYPE* yylval, gmx_ana_selparam_t* param, bool bBoolNo)
{
    if (bBoolNo)
    {
        GMX_RELEASE_ASSERT(param->name != nullptr,
                           "bBoolNo should only be set for a parameters with a name");
        snew(yylval->str, strlen(param->name) + 3);
        yylval->str[0] = 'n';
        yylval->str[1] = 'o';
        strcpy(yylval->str + 2, param->name);
    }
    else
    {
        yylval->str = param->name ? gmx_strdup(param->name) : nullptr;
    }
    return PARAM;
}

/*! \brief
 * Processes a selection method token.
 */
static int init_method_token(YYSTYPE*                          yylval,
                             YYLTYPE*                          yylloc,
                             const gmx::SelectionParserSymbol* symbol,
                             bool                              bPosMod,
                             gmx_sel_lexer_t*                  state)
{
    gmx_ana_selmethod_t* method = symbol->methodValue();
    /* If the previous token was not KEYWORD_POS, return EMPTY_POSMOD
     * before the actual method to work around a limitation in Bison. */
    if (!bPosMod && method->type != POS_VALUE)
    {
        state->nextMethodSymbol = symbol;
        _gmx_sel_lexer_add_token(yylloc, nullptr, 0, state);
        return EMPTY_POSMOD;
    }
    _gmx_sel_lexer_add_token(yylloc, symbol->name().c_str(), -1, state);
    yylval->meth = method;
    if (!(method->flags & SMETH_MODIFIER) && method->nparams == 0)
    {
        /* Keyword */
        switch (method->type)
        {
            case INT_VALUE:
            case REAL_VALUE: state->bMatchOf = true; return KEYWORD_NUMERIC;
            case STR_VALUE: return KEYWORD_STR;
            case GROUP_VALUE: return KEYWORD_GROUP;
            default: GMX_THROW(gmx::InternalError("Unsupported keyword type"));
        }
    }
    else
    {
        /* Method with parameters or a modifier */
        if (method->flags & SMETH_MODIFIER)
        {
            /* Remove all methods from the stack */
            state->msp = -1;
            if (method->param[1].name == nullptr)
            {
                state->nextparam = &method->param[1];
            }
        }
        else
        {
            if (method->param[0].name == nullptr)
            {
                state->nextparam = &method->param[0];
            }
        }
        ++state->msp;
        if (state->msp >= state->mstack_alloc)
        {
            state->mstack_alloc += 10;
            srenew(state->mstack, state->mstack_alloc);
        }
        state->mstack[state->msp] = method;
        if (method->flags & SMETH_MODIFIER)
        {
            return MODIFIER;
        }
        switch (method->type)
        {
            case INT_VALUE: // Intended fall through
            case REAL_VALUE: return METHOD_NUMERIC;
            case POS_VALUE: return METHOD_POS;
            case GROUP_VALUE: return METHOD_GROUP;
            default: --state->msp; GMX_THROW(gmx::InternalError("Unsupported method type"));
        }
    }
}

int _gmx_sel_lexer_process_pending(YYSTYPE* yylval, YYLTYPE* yylloc, gmx_sel_lexer_t* state)
{
    if (state->nextparam)
    {
        gmx_ana_selparam_t* param   = state->nextparam;
        bool                bBoolNo = state->bBoolNo;

        if (state->neom > 0)
        {
            --state->neom;
            _gmx_sel_lexer_add_token(yylloc, nullptr, 0, state);
            return END_OF_METHOD;
        }
        state->nextparam = nullptr;
        state->bBoolNo   = false;
        _gmx_sel_lexer_add_token(yylloc, param->name, -1, state);
        return init_param_token(yylval, param, bBoolNo);
    }
    if (state->prev_pos_kw > 0)
    {
        --state->prev_pos_kw;
    }
    if (state->nextMethodSymbol)
    {
        const gmx::SelectionParserSymbol* symbol = state->nextMethodSymbol;
        state->nextMethodSymbol                  = nullptr;
        return init_method_token(yylval, yylloc, symbol, true, state);
    }
    return 0;
}

int _gmx_sel_lexer_process_identifier(YYSTYPE* yylval, YYLTYPE* yylloc, char* yytext, size_t yyleng, gmx_sel_lexer_t* state)
{
    /* Check if the identifier matches with a parameter name */
    if (state->msp >= 0)
    {
        gmx_ana_selparam_t* param   = nullptr;
        bool                bBoolNo = false;
        int                 sp      = state->msp;
        while (!param && sp >= 0)
        {
            int i;
            for (i = 0; i < state->mstack[sp]->nparams; ++i)
            {
                /* Skip NULL parameters and too long parameters */
                if (state->mstack[sp]->param[i].name == nullptr
                    || strlen(state->mstack[sp]->param[i].name) > yyleng)
                {
                    continue;
                }
                if (!strncmp(state->mstack[sp]->param[i].name, yytext, yyleng))
                {
                    param = &state->mstack[sp]->param[i];
                    break;
                }
                /* Check separately for a 'no' prefix on boolean parameters */
                if (state->mstack[sp]->param[i].val.type == NO_VALUE && yyleng > 2
                    && yytext[0] == 'n' && yytext[1] == 'o'
                    && !strncmp(state->mstack[sp]->param[i].name, yytext + 2, yyleng - 2))
                {
                    param   = &state->mstack[sp]->param[i];
                    bBoolNo = true;
                    break;
                }
            }
            if (!param)
            {
                --sp;
            }
        }
        if (param)
        {
            if (param->val.type == NO_VALUE && !bBoolNo)
            {
                state->bMatchBool = true;
            }
            if (sp < state->msp)
            {
                state->neom      = state->msp - sp - 1;
                state->nextparam = param;
                state->bBoolNo   = bBoolNo;
                return END_OF_METHOD;
            }
            _gmx_sel_lexer_add_token(yylloc, param->name, -1, state);
            return init_param_token(yylval, param, bBoolNo);
        }
    }

    /* Check if the identifier matches with a symbol */
    const gmx::SelectionParserSymbol* symbol = state->sc->symtab->findSymbol(std::string(yytext, yyleng));
    /* If there is no match, return the token as a string */
    if (!symbol)
    {
        yylval->str = gmx_strndup(yytext, yyleng);
        _gmx_sel_lexer_add_token(yylloc, yytext, yyleng, state);
        return IDENTIFIER;
    }
    gmx::SelectionParserSymbol::SymbolType symtype = symbol->type();
    /* For method symbols, we need some extra processing. */
    if (symtype == gmx::SelectionParserSymbol::MethodSymbol)
    {
        return init_method_token(yylval, yylloc, symbol, state->prev_pos_kw > 0, state);
    }
    _gmx_sel_lexer_add_token(yylloc, symbol->name().c_str(), -1, state);
    /* Reserved symbols should have been caught earlier */
    if (symtype == gmx::SelectionParserSymbol::ReservedSymbol)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString(
                "Mismatch between tokenizer and reserved symbol table (for '%s')", symbol->name().c_str())));
    }
    /* For variable symbols, return the type of the variable value */
    if (symtype == gmx::SelectionParserSymbol::VariableSymbol)
    {
        const gmx::SelectionTreeElementPointer& var = symbol->variableValue();
        /* Return simple tokens for constant variables */
        if (var->type == SEL_CONST)
        {
            switch (var->v.type)
            {
                case INT_VALUE: yylval->i = var->v.u.i[0]; return TOK_INT;
                case REAL_VALUE: yylval->r = var->v.u.r[0]; return TOK_REAL;
                case POS_VALUE: break;
                default: GMX_THROW(gmx::InternalError("Unsupported variable type"));
            }
        }
        yylval->sel = new gmx::SelectionTreeElementPointer(var);
        switch (var->v.type)
        {
            case INT_VALUE: // Intended fall through
            case REAL_VALUE: return VARIABLE_NUMERIC;
            case POS_VALUE: return VARIABLE_POS;
            case GROUP_VALUE: return VARIABLE_GROUP;
            default: delete yylval->sel; GMX_THROW(gmx::InternalError("Unsupported variable type"));
        }
        /* This position should not be reached. */
    }
    /* For position symbols, we need to return KEYWORD_POS, but we also need
     * some additional handling. */
    if (symtype == gmx::SelectionParserSymbol::PositionSymbol)
    {
        state->bMatchOf    = true;
        yylval->str        = gmx_strdup(symbol->name().c_str());
        state->prev_pos_kw = 2;
        return KEYWORD_POS;
    }
    /* Should not be reached */
    return INVALID;
}

void _gmx_sel_lexer_add_token(YYLTYPE* yylloc, const char* str, int len, gmx_sel_lexer_t* state)
{
    yylloc->startIndex = yylloc->endIndex = state->pselstr.size();
    /* Do nothing if the string is empty, or if it is a space and there is
     * no other text yet, or if there already is a space. */
    if (!str || len == 0 || strlen(str) == 0
        || (str[0] == ' ' && str[1] == 0 && (state->pselstr.empty() || state->pselstr.back() == ' ')))
    {
        return;
    }
    if (len < 0)
    {
        len = strlen(str);
    }
    /* Append the token to the stored string */
    state->pselstr.append(str, len);
    yylloc->endIndex = state->pselstr.size();
}

void _gmx_sel_init_lexer(yyscan_t*                       scannerp,
                         struct gmx_ana_selcollection_t* sc,
                         gmx::TextWriter*                statusWriter,
                         int                             maxnr,
                         bool                            bGroups,
                         struct gmx_ana_indexgrps_t*     grps)
{
    int rc = _gmx_sel_yylex_init(scannerp);
    if (rc != 0)
    {
        // TODO: Throw a more representative exception.
        GMX_THROW(gmx::InternalError("Lexer initialization failed"));
    }

    gmx_sel_lexer_t* state = new gmx_sel_lexer_t;

    state->sc      = sc;
    state->bGroups = bGroups;
    state->grps    = grps;
    state->nexpsel = (maxnr > 0 ? gmx::ssize(sc->sel) + maxnr : -1);

    state->statusWriter = statusWriter;

    state->currentLocation.startIndex = 0;
    state->currentLocation.endIndex   = 0;

    snew(state->mstack, 20);
    state->mstack_alloc     = 20;
    state->msp              = -1;
    state->neom             = 0;
    state->nextparam        = nullptr;
    state->nextMethodSymbol = nullptr;
    state->prev_pos_kw      = 0;
    state->bBoolNo          = false;
    state->bMatchOf         = false;
    state->bMatchBool       = false;
    state->bCmdStart        = true;
    state->bBuffer          = false;

    _gmx_sel_yyset_extra(state, *scannerp);
}

void _gmx_sel_free_lexer(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);

    sfree(state->mstack);
    if (state->bBuffer)
    {
        _gmx_sel_yy_delete_buffer(state->buffer, scanner);
    }
    delete state;
    _gmx_sel_yylex_destroy(scanner);
}

void _gmx_sel_lexer_set_exception(yyscan_t scanner, const std::exception_ptr& ex)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    state->exception       = ex;
}

void _gmx_sel_lexer_rethrow_exception_if_occurred(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    if (state->exception)
    {
        std::exception_ptr ex = state->exception;
        state->exception      = std::exception_ptr();
        std::rethrow_exception(ex);
    }
}

gmx::TextWriter* _gmx_sel_lexer_get_status_writer(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    return state->statusWriter;
}

struct gmx_ana_selcollection_t* _gmx_sel_lexer_selcollection(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    return state->sc;
}

bool _gmx_sel_lexer_has_groups_set(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    return state->bGroups;
}

struct gmx_ana_indexgrps_t* _gmx_sel_lexer_indexgrps(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    return state->grps;
}

int _gmx_sel_lexer_exp_selcount(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    return state->nexpsel;
}

const char* _gmx_sel_lexer_pselstr(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    return state->pselstr.c_str();
}

void _gmx_sel_lexer_set_current_location(yyscan_t scanner, const gmx::SelectionLocation& location)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    state->currentLocation = location;
}

const gmx::SelectionLocation& _gmx_sel_lexer_get_current_location(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    return state->currentLocation;
}

std::string _gmx_sel_lexer_get_current_text(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    return _gmx_sel_lexer_get_text(scanner, state->currentLocation);
}

std::string _gmx_sel_lexer_get_text(yyscan_t scanner, const gmx::SelectionLocation& location)
{
    gmx_sel_lexer_t* state      = _gmx_sel_yyget_extra(scanner);
    const int        startIndex = location.startIndex;
    const int        endIndex   = location.endIndex;
    if (startIndex >= endIndex)
    {
        return std::string();
    }
    return state->pselstr.substr(startIndex, endIndex - startIndex);
}

void _gmx_sel_lexer_clear_pselstr(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);
    state->pselstr.clear();
}

void _gmx_sel_lexer_clear_method_stack(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);

    state->msp = -1;
}

void _gmx_sel_finish_method(yyscan_t scanner)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);

    if (state->msp >= 0)
    {
        --state->msp;
    }
}

void _gmx_sel_set_lex_input_file(yyscan_t scanner, FILE* fp)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);

    state->bBuffer = true;
    state->buffer  = _gmx_sel_yy_create_buffer(fp, YY_BUF_SIZE, scanner);
    _gmx_sel_yy_switch_to_buffer(state->buffer, scanner);
}

void _gmx_sel_set_lex_input_str(yyscan_t scanner, const char* str)
{
    gmx_sel_lexer_t* state = _gmx_sel_yyget_extra(scanner);

    if (state->bBuffer)
    {
        _gmx_sel_yy_delete_buffer(state->buffer, scanner);
    }
    state->bBuffer = true;
    state->buffer  = _gmx_sel_yy_scan_string(str, scanner);
}
