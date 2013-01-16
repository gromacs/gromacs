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
 * \brief Implementation of functions in symrec.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <macros.h>
#include <smalloc.h>
#include <string2.h>
#include <typedefs.h>
#include <gmx_fatal.h>

#include <poscalc.h>

#include "selelem.h"
#include "symrec.h"

/*! \internal \brief
 * Symbol table for the selection parser.
 */
struct gmx_sel_symtab_t
{
    /** Pointer to the first symbol in the linked list of symbols. */
    gmx_sel_symrec_t *first;
};

/*! \internal \brief
 * Single symbol for the selection parser.
 */
struct gmx_sel_symrec_t
{
    /** Name of the symbol. */
    char                           *name;
    /** Type of the symbol. */
    e_symbol_t                      type;
    /** Value of the symbol. */
    union {
        /** Pointer to the method structure (\ref SYMBOL_METHOD). */
        struct gmx_ana_selmethod_t *meth;
        /** Pointer to the variable value (\ref SYMBOL_VARIABLE). */
        struct t_selelem           *var;
    }                               u;
    /** Pointer to the next symbol. */
    struct gmx_sel_symrec_t        *next;
};

/** List of reserved symbols to register in add_reserved_symbols(). */
static const char *const sym_reserved[] = {
    "group",
    "to",
    "not",
    "and",
    "or",
    "xor",
    "yes",
    "no",
    "on",
    "off",
    "help",
};

/*!
 * \param[in] sym Symbol to query.
 * \returns   The name of \p sym.
 *
 * The returned pointer should not be free'd.
 */
char *
_gmx_sel_sym_name(gmx_sel_symrec_t *sym)
{
    return sym->name;
}

/*!
 * \param[in] sym Symbol to query.
 * \returns   The type of \p sym.
 */
e_symbol_t
_gmx_sel_sym_type(gmx_sel_symrec_t *sym)
{
    return sym->type;
}

/*!
 * \param[in] sym Symbol to query.
 * \returns   The method associated with \p sym, or NULL if \p sym is not a
 *   \ref SYMBOL_METHOD symbol.
 */
struct gmx_ana_selmethod_t *
_gmx_sel_sym_value_method(gmx_sel_symrec_t *sym)
{
    if (sym->type != SYMBOL_METHOD)
    {
        gmx_call("symbol is not a method symbol");
        return NULL;
    }
    return sym->u.meth;
}

/*!
 * \param[in] sym Symbol to query.
 * \returns   The variable expression associated with \p sym, or NULL if
 *   \p sym is not a \ref SYMBOL_VARIABLE symbol.
 */
struct t_selelem *
_gmx_sel_sym_value_var(gmx_sel_symrec_t *sym)
{
    if (sym->type != SYMBOL_VARIABLE)
    {
        gmx_call("symbol is not a variable symbol");
        return NULL;
    }
    return sym->u.var;
}

/*! \brief
 * Adds the reserved symbols to a symbol table.
 *
 * \param[in,out] tab  Symbol table to which the symbols are added.
 *
 * Assumes that the symbol table is empty.
 */
static void
add_reserved_symbols(gmx_sel_symtab_t *tab)
{
    gmx_sel_symrec_t *sym;
    gmx_sel_symrec_t *last;
    size_t            i;

    last = NULL;
    for (i = 0; i < asize(sym_reserved); ++i)
    {
        snew(sym, 1);
        sym->name = strdup(sym_reserved[i]);
        sym->type = SYMBOL_RESERVED;
        sym->next = NULL;
        if (last)
        {
            last->next = sym;
        }
        else
        {
            tab->first = sym;
        }
        last = sym;
    }
}

/*! \brief
 * Adds the position symbols to the symbol list.
 *
 * \param[in,out] tab  Symbol table to which the symbols are added.
 */
static void
add_position_symbols(gmx_sel_symtab_t *tab)
{
    const char       **postypes;
    gmx_sel_symrec_t  *sym;
    gmx_sel_symrec_t  *last;
    int                i;

    postypes = gmx_ana_poscalc_create_type_enum(TRUE);
    last     = tab->first;
    while (last && last->next)
    {
        last = last->next;
    }
    for (i = 1; postypes[i] != NULL; ++i)
    {
        snew(sym, 1);
        sym->name = strdup(postypes[i]);
        sym->type = SYMBOL_POS;
        sym->next = NULL;
        if (last)
        {
            last->next = sym;
        }
        else
        {
            tab->first = sym;
        }
        last = sym;
    }
    sfree(postypes);
}

/*!
 * \param[out] tabp Symbol table pointer to initialize.
 *
 * Reserved and position symbols are added to the created table.
 */
int
_gmx_sel_symtab_create(gmx_sel_symtab_t **tabp)
{
    gmx_sel_symtab_t *tab;

    snew(tab, 1);
    add_reserved_symbols(tab);
    add_position_symbols(tab);
    *tabp = tab;
    return 0;
}

/*!
 * \param[in] tab Symbol table to free.
 *
 * The pointer \p tab is invalid after the call.
 */
void
_gmx_sel_symtab_free(gmx_sel_symtab_t *tab)
{
    gmx_sel_symrec_t *sym;

    while (tab->first)
    {
        sym        = tab->first;
        tab->first = sym->next;
        if (sym->type == SYMBOL_VARIABLE)
        {
            _gmx_selelem_free(sym->u.var);
        }
        sfree(sym->name);
        sfree(sym);
    }
    sfree(tab);
}

/*!
 * \param[in] tab    Symbol table to search.
 * \param[in] name   Symbol name to find.
 * \param[in] bExact If FALSE, symbols that begin with \p name are also
 *   considered.
 * \returns   Pointer to the symbol with name \p name, or NULL if not found.
 *
 * If no exact match is found and \p bExact is FALSE, returns a symbol that
 * begins with \p name if a unique matching symbol is found.
 */
gmx_sel_symrec_t *
_gmx_sel_find_symbol(gmx_sel_symtab_t *tab, const char *name, gmx_bool bExact)
{
    return _gmx_sel_find_symbol_len(tab, name, strlen(name), bExact);
}

/*!
 * \param[in] tab    Symbol table to search.
 * \param[in] name   Symbol name to find.
 * \param[in] len    Only consider the first \p len characters of \p name.
 * \param[in] bExact If FALSE, symbols that begin with \p name are also
 *   considered.
 * \returns   Pointer to the symbol with name \p name, or NULL if not found.
 *
 * If no exact match is found and \p bExact is FALSE, returns a symbol that
 * begins with \p name if a unique matching symbol is found.
 *
 * The parameter \p len is there to allow using this function from scanner.l
 * without modifying the text to be scanned or copying it.
 */
gmx_sel_symrec_t *
_gmx_sel_find_symbol_len(gmx_sel_symtab_t *tab, const char *name, size_t len,
                         gmx_bool bExact)
{
    gmx_sel_symrec_t     *sym;
    gmx_sel_symrec_t     *match;
    gmx_bool              bUnique;
    gmx_bool              bMatch;

    match   = NULL;
    bUnique = TRUE;
    bMatch  = FALSE;
    sym     = tab->first;
    while (sym)
    {
        if (!strncmp(sym->name, name, len))
        {
            if (strlen(sym->name) == len)
            {
                return sym;
            }
            if (bMatch)
            {
                bUnique = FALSE;
            }
            bMatch = TRUE;
            if (sym->type == SYMBOL_METHOD)
            {
                match = sym;
            }
        }
        sym = sym->next;
    }
    if (bExact)
    {
        return NULL;
    }

    if (!bUnique)
    {
        fprintf(stderr, "parse error: ambiguous symbol\n");
        return NULL;
    }
    return match;
}

/*!
 * \param[in] tab   Symbol table to search.
 * \param[in] type  Type of symbol to find.
 * \returns   The first symbol in \p tab with type \p type,
 *   or NULL if there are no such symbols.
 */
gmx_sel_symrec_t *
_gmx_sel_first_symbol(gmx_sel_symtab_t *tab, e_symbol_t type)
{
    gmx_sel_symrec_t *sym;

    sym = tab->first;
    while (sym)
    {
        if (sym->type == type)
        {
            return sym;
        }
        sym = sym->next;
    }
    return NULL;
}

/*!
 * \param[in] after Start the search after this symbol.
 * \param[in] type  Type of symbol to find.
 * \returns   The next symbol after \p after with type \p type,
 *   or NULL if there are no more symbols.
 */
gmx_sel_symrec_t *
_gmx_sel_next_symbol(gmx_sel_symrec_t *after, e_symbol_t type)
{
    gmx_sel_symrec_t *sym;

    sym = after->next;
    while (sym)
    {
        if (sym->type == type)
        {
            return sym;
        }
        sym = sym->next;
    }
    return NULL;
}

/*! \brief
 * Internal utility function used in adding symbols to a symbol table.
 *
 * \param[in,out] tab   Symbol table to add the symbol to.
 * \param[in]     name  Name of the symbol to add.
 * \param[out]    ctype On error, the type of the conflicting symbol is
 *   written to \p *ctype.
 * \returns       Pointer to the new symbol record, or NULL if \p name
 *   conflicts with an existing symbol.
 */
static gmx_sel_symrec_t *
add_symbol(gmx_sel_symtab_t *tab, const char *name, e_symbol_t *ctype)
{
    gmx_sel_symrec_t *sym, *psym;
    int               len;

    /* Check if there is a conflicting symbol */
    psym = NULL;
    sym  = tab->first;
    while (sym)
    {
        if (!gmx_strcasecmp(sym->name, name))
        {
            *ctype = sym->type;
            return NULL;
        }
        psym = sym;
        sym  = sym->next;
    }

    /* Create a new symbol record */
    if (psym == NULL)
    {
        snew(tab->first, 1);
        sym = tab->first;
    }
    else
    {
        snew(psym->next, 1);
        sym = psym->next;
    }
    sym->name = strdup(name);
    return sym;
}

/*!
 * \param[in,out] tab    Symbol table to add the symbol to.
 * \param[in]     name   Name of the new symbol.
 * \param[in]     sel    Value of the variable.
 * \returns       Pointer to the created symbol record, or NULL if there was a
 *   symbol with the same name.
 */
gmx_sel_symrec_t *
_gmx_sel_add_var_symbol(gmx_sel_symtab_t *tab, const char *name,
                        struct t_selelem *sel)
{
    gmx_sel_symrec_t *sym;
    e_symbol_t        ctype;

    sym = add_symbol(tab, name, &ctype);
    if (!sym)
    {
        fprintf(stderr, "parse error: ");
        switch (ctype)
        {
            case SYMBOL_RESERVED:
            case SYMBOL_POS:
                fprintf(stderr, "variable name (%s) conflicts with a reserved keyword\n",
                        name);
                break;
            case SYMBOL_VARIABLE:
                fprintf(stderr, "duplicate variable name (%s)\n", name);
                break;
            case SYMBOL_METHOD:
                fprintf(stderr, "variable name (%s) conflicts with a selection keyword\n",
                        name);
                break;
        }
        return NULL;
    }

    sym->type  = SYMBOL_VARIABLE;
    sym->u.var = sel;
    sel->refcount++;
    return sym;
}

/*!
 * \param[in,out] tab    Symbol table to add the symbol to.
 * \param[in]     name   Name of the new symbol.
 * \param[in]     method Method that this symbol represents.
 * \returns       Pointer to the created symbol record, or NULL if there was a
 *   symbol with the same name.
 */
gmx_sel_symrec_t *
_gmx_sel_add_method_symbol(gmx_sel_symtab_t *tab, const char *name,
                           struct gmx_ana_selmethod_t *method)
{
    gmx_sel_symrec_t *sym;
    e_symbol_t        ctype;

    sym = add_symbol(tab, name, &ctype);
    if (!sym)
    {
        fprintf(stderr, "parse error: ");
        switch (ctype)
        {
            case SYMBOL_RESERVED:
            case SYMBOL_POS:
                fprintf(stderr, "method name (%s) conflicts with a reserved keyword\n",
                        name);
                break;
            case SYMBOL_VARIABLE:
                fprintf(stderr, "method name (%s) conflicts with a variable name\n",
                        name);
                break;
            case SYMBOL_METHOD:
                fprintf(stderr, "duplicate method name (%s)\n", name);
                break;
        }
        return NULL;
    }

    sym->type   = SYMBOL_METHOD;
    sym->u.meth = method;
    return sym;
}
