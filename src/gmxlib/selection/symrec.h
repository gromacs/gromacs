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
 * \brief Handling of selection parser symbol table.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 */
#ifndef SELECTION_SYMREC_H
#define SELECTION_SYMREC_H

struct t_selelem;
struct gmx_ana_selmethod_t;

/** Defines the type of the symbol. */
typedef enum
{
    SYMBOL_RESERVED,    /**< The symbol is a reserved keyword. */
    SYMBOL_VARIABLE,    /**< The symbol is a variable. */
    SYMBOL_METHOD,      /**< The symbol is a selection method. */
    SYMBOL_POS          /**< The symbol is a position keyword. */
} e_symbol_t;

/** Symbol table for the selection parser. */
typedef struct gmx_sel_symtab_t gmx_sel_symtab_t;
/** Single symbol for the selection parser. */
typedef struct gmx_sel_symrec_t gmx_sel_symrec_t;

/** Returns the name of a symbol. */
char *
_gmx_sel_sym_name(gmx_sel_symrec_t *sym);
/** Returns the type of a symbol. */
e_symbol_t
_gmx_sel_sym_type(gmx_sel_symrec_t *sym);
/** Returns the method associated with a \ref SYMBOL_METHOD symbol. */
struct gmx_ana_selmethod_t *
_gmx_sel_sym_value_method(gmx_sel_symrec_t *sym);
/** Returns the method associated with a \ref SYMBOL_VARIABLE symbol. */
struct t_selelem *
_gmx_sel_sym_value_var(gmx_sel_symrec_t *sym);

/** Creates a new symbol table. */
int
_gmx_sel_symtab_create(gmx_sel_symtab_t **tabp);
/** Frees all memory allocated for a symbol table. */
void
_gmx_sel_symtab_free(gmx_sel_symtab_t *tab);
/** Finds a symbol by name. */
gmx_sel_symrec_t *
_gmx_sel_find_symbol(gmx_sel_symtab_t *tab, const char *name, gmx_bool bExact);
/** Finds a symbol by name. */
gmx_sel_symrec_t *
_gmx_sel_find_symbol_len(gmx_sel_symtab_t *tab, const char *name, size_t len,
                         gmx_bool bExact);
/** Returns the first symbol of a given type. */
gmx_sel_symrec_t *
_gmx_sel_first_symbol(gmx_sel_symtab_t *tab, e_symbol_t type);
/** Returns the next symbol of a given type. */
gmx_sel_symrec_t *
_gmx_sel_next_symbol(gmx_sel_symrec_t *after, e_symbol_t type);
/** Adds a new variable symbol. */
gmx_sel_symrec_t *
_gmx_sel_add_var_symbol(gmx_sel_symtab_t *tab, const char *name,
                        struct t_selelem *sel);
/** Adds a new method symbol. */
gmx_sel_symrec_t *
_gmx_sel_add_method_symbol(gmx_sel_symtab_t *tab, const char *name,
                           struct gmx_ana_selmethod_t *method);

#endif
