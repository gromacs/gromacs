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
 * \brief Help implementation for selections.
 */
#include <macros.h>
#include <string2.h>
#include <wman.h>

#include "selcollection.h"
#include "selhelp.h"
#include "symrec.h"

typedef struct {
    const char  *topic;
    int          nl;
    const char **text;
} t_selection_help_item;

static const char *help_common[] = {
    "SELECTION HELP[PAR]",

    "Please read the subtopic pages (available through \"help topic\") for",
    "more information.",
};

static const char *help_cmdline[] = {
    "SELECTION COMMAND-LINE ARGUMENTS[PAR]",

    "There are two alternative command-line arguments for specifying",
    "selections:[BR]",
    "1. [TT]-select[tt] can be used to specify the complete selection as a",
    "string on the command line.[BR]",
    "2. [TT]-sf[tt] can be used to specify a file name from which the",
    "selection is read.[BR]",
    "If both options are specified, [TT]-select[tt] takes precedence.",
    "If neither of the above is present, the user is prompted to type the",
    "selection on the standard input (a pipe can also be used to provide",
    "the selections in this case.",
    "This is also done if an empty string is passed to [TT]-select[tt].[PAR]",

    "Option [TT]-n[tt] can be used to provide an index file.",
    "If no index file is provided, the default groups are generated.",
    "In both cases, the user can also select an index group instead of",
    "writing a full selection.[PAR]",

    "Depending on the tool, two additional command-line arguments may be",
    "provided:[BR]",
    "1. [TT]-seltype[tt] can be used to specify the default type of",
    "positions to calculate for each selection.[BR]",
    "2. [TT]-selrpos[tt] can be used to specify the default type of",
    "positions used in selecting atoms by coordinates.[BR]",
    "See \"help positions\" for more information on these options.",
};

static const char *help_eval[] = {
    "SELECTION EVALUATION AND OPTIMIZATION[PAR]",

    "Boolean evaluation proceeds from left to right and is short-circuiting",
    "i.e., as soon as it is known whether an atom will be selected, the",
    "remaining expressions are not evaluated at all.",
    "This can be used to optimize the selections: you should write the",
    "most restrictive and/or the most inexpensive expressions first in",
    "boolean expressions.",
    "The relative ordering between dynamic and static expressions does not",
    "matter: all static expressions are evaluated only once, before the first",
    "frame, and the result becomes the leftmost expression.[PAR]",

    "Another point for optimization is in common subexpressions: they are not",
    "automatically recognized, but can be manually optimized by the use of",
    "variables. This can have a big impact on the performance of complex",
    "selections, in particular if you define several index groups like this:",
    "  [TT]rdist = distance from com of resnr 1 to 5;[tt][BR]",
    "  [TT]resname RES and rdist < 2;[tt][BR]",
    "  [TT]resname RES and rdist < 4;[tt][BR]",
    "  [TT]resname RES and rdist < 6;[tt][BR]",
    "Without the variable assignment, the distances would be evaluated three",
    "times, although they are exactly the same within each selection.",
    "Anything assigned into a variable becomes a common subexpression that",
    "is evaluated only once during a frame.",
    "Currently, in some cases the use of variables can actually lead to a small",
    "performance loss because of the checks necessary to determine for which",
    "atoms the expression has already been evaluated, but this should not be",
    "a major problem.",
};

static const char *help_keywords[] = {
    "SELECTION KEYWORDS[PAR]",

    "The following selection keywords are currently available.",
};

static const char *help_limits[] = {
    "SELECTION LIMITATIONS[PAR]",

    "Arithmetic expressions are not implemented.[PAR]",

    "Some analysis programs may require a special structure for the input",
    "selections (e.g., [TT]g_angle[tt] requires the index group to be made",
    "of groups of three or four atoms).",
    "For such programs, it is up to the user to provide a proper selection",
    "expression that always returns such positions.[PAR]",

    "Currently, it is not possible to write selection expressions that would",
    "evaluate to index groups where the same atom occurs more than once.",
};

static const char *help_positions[] = {
    "SPECIFYING POSITIONS[PAR]",

    "Possible ways of specifying positions in selections are:[PAR]",

    "1. A constant position can be defined as [TT](XX, YY, ZZ)[tt], where",
    "[TT]XX[tt], [TT]YY[tt] and [TT]ZZ[tt] are real numbers.[PAR]",

    "2. [TT]com of ATOM_EXPR [pbc][tt] or [TT]cog of ATOM_EXPR [pbc][tt]",
    "calculate the center of mass/geometry of [TT]ATOM_EXPR[tt]. If",
    "[TT]pbc[tt] is specified, the center is calculated iteratively to try",
    "to deal with cases where [TT]ATOM_EXPR[tt] wraps around periodic",
    "boundary conditions.[PAR]",

    "3. [TT]POSTYPE of ATOM_EXPR[tt] calculates the specified positions for",
    "the atoms in [TT]ATOM_EXPR[tt].",
    "[TT]POSTYPE[tt] can be [TT]atom[tt], [TT]res_com[tt], [TT]res_cog[tt],",
    "[TT]mol_com[tt] or [TT]mol_cog[tt], with an optional prefix [TT]whole_[tt]",
    "[TT]part_[tt] or [TT]dyn_[tt].",
    "[TT]whole_[tt] calcualtes the centers for the whole residue/molecule,",
    "even if only part of it is selected.",
    "[TT]part_[tt] prefix calculates the centers for the selected atoms, but",
    "uses always the same atoms for the same residue/molecule. The used atoms",
    "are determined from the the largest group allowed by the selection.",
    "[TT]dyn_[tt] calculates the centers strictly only for the selected atoms.",
    "If no prefix is specified, whole selections default to [TT]part_[tt] and",
    "other places default to [TT]whole_[tt].",
    "The latter is often desirable to select the same molecules in different",
    "tools, while the first is a compromise between speed ([TT]dyn_[tt]",
    "positions can be slower to evaluate than [TT]part_[tt]) and intuitive",
    "behavior.[PAR]",

    "4. [TT]ATOM_EXPR[tt], when given for whole selections, is handled as 3.",
    "above, using the position type from the command-line argument",
    "[TT]-seltype[tt].[PAR]",

    "Selection keywords that select atoms based on their positions, such as",
    "[TT]dist from[tt], use by default the positions defined by the",
    "[TT]-selrpos[tt] command-line option.",
    "This can be overridden by prepending a [TT]POSTYPE[tt] specifier to the",
    "keyword. For example, [TT]res_com dist from POS[tt] evaluates the",
    "residue center of mass distances. In the example, all atoms of a residue",
    "are either selected or not, based on the single distance calculated.",
};

static const char *help_syntax[] = {
    "SELECTION SYNTAX[PAR]",

    "A set of selections consists of one or more selections, separated by",
    "semicolons. Each selection defines a set of positions for the analysis.",
    "Each selection can also be preceded by a string that gives a name for",
    "the selection for use in, e.g., graph legends.",
    "If no name is provided, the string used for the selection is used",
    "automatically as the name.[PAR]",

    "It is possible to use variables to store selection expressions.",
    "A variable is defined with the following syntax:[BR]",
    "[TT]VARNAME = EXPR ;[tt][BR]",
    "where [TT]EXPR[tt] is any valid selection expression.",
    "After this, [TT]VARNAME[tt] can be used anywhere where [TT]EXPR[tt]",
    "would be valid.[PAR]",

    "For interactive input, the syntax is slightly altered: line breaks can",
    "also be used to separate selections. \\ followed by a line break can",
    "be used to continue a line if necessary.",
    "Notice that the above only applies to real interactive input,",
    "not if you provide the selections, e.g., from a pipe.[PAR]",

    "Selections are composed of three main types of expressions, those that",
    "define atoms ([TT]ATOM_EXPR[tt]s), those that define positions",
    "([TT]POS_EXPR[tt]s), and those that evaluate to numeric values",
    "([TT]NUM_EXPR[tt]s). Each selection should be a [TT]POS_EXPR[tt]",
    "or a [TT]ATOM_EXPR[tt] (the latter is automatically converted to",
    "positions). The basic rules are as follows:[BR]",
    "1. An expression like [TT]NUM_EXPR1 < NUM_EXPR2[tt] evaluates to an",
    "[TT]ATOM_EXPR[tt] that selects all the atoms for which the comparison",
    "is true.[BR]",
    "2. Atom expressions can be combined with boolean operations such as",
    "[TT]not ATOM_EXPR[tt], [TT]ATOM_EXPR and ATOM_EXPR[tt], or",
    "[TT]ATOM_EXPR or ATOM_EXPR[tt]. Parentheses can be used to alter the",
    "evaluation order.[BR]",
    "3. [TT]ATOM_EXPR[tt] expressions can be converted into [TT]POS_EXPR[tt]",
    "expressions in various ways, see \"help positions\" for more details.[PAR]",

    "Some keywords select atoms based on string values such as the atom name.",
    "For these keywords, it is possible to use wildcards ([TT]name \"C*\"[tt])",
    "or regular expressions (e.g., [TT]resname \"R[AB]\"[tt]).",
    "The match type is automatically guessed from the string: if it contains",
    "other characters than letters, numbers, '*', or '?', it is interpreted",
    "as a regular expression.",
    "Strings that contain non-alphanumeric characters should be enclosed in",
    "double quotes as in the examples. For other strings, the quotes are",
    "optional, but if the value conflicts with a reserved keyword, a syntax",
    "error will occur. If your strings contain uppercase letters, this should",
    "not happen."
};

static const t_selection_help_item helpitems[] = {
    {NULL,          asize(help_common),    help_common},
    {"cmdline",     asize(help_cmdline),   help_cmdline},
    {"evaluation",  asize(help_eval),      help_eval},
    {"keywords",    asize(help_keywords),  help_keywords},
    {"limitations", asize(help_limits),    help_limits},
    {"positions",   asize(help_positions), help_positions},
    {"syntax",      asize(help_syntax),    help_syntax},
};

/*!
 * \param[in]  sc    Selection collection for which help should be printed.
 * \param[in]  topic Topic to print help on, or NULL for general help.
 *
 * \p sc is used to get information on which keywords are available in the
 * present context.
 */
void
_gmx_sel_print_help(struct gmx_ana_selcollection_t *sc, const char *topic)
{
    const t_selection_help_item *item = 0;
    int  i;

    /* Find the item for the topic */
    if (!topic)
    {
        item = &helpitems[0];
    }
    else
    {
        for (i = 1; i < asize(helpitems); ++i)
        {
            if (strncmp(helpitems[i].topic, topic, strlen(topic)) == 0)
            {
                item = &helpitems[i];
            }
        }
    }
    /* If the topic is not found, tell the user and exit. */
    if (!item)
    {
        fprintf(stderr, "No help available for '%s'.\n", topic);
        return;
    }
    /* Print the help */
    print_tty_formatted(stderr, item->nl, item->text, 0, NULL, NULL, FALSE);
    /* Special handling of certain pages */
    if (!topic)
    {
        /* Print the subtopics on the main page */
        fprintf(stderr, "\nAvailable subtopics:\n");
        for (i = 1; i < asize(helpitems); ++i)
        {
            fprintf(stderr, "  %s", helpitems[i].topic);
        }
        fprintf(stderr, "\n");
    }
    else if (strcmp(item->topic, "keywords") == 0)
    {
        /* Print the list of keywords */
        gmx_sel_symrec_t *symbol;

        fprintf(stderr, "\nAvailable keywords:\n");
        symbol = _gmx_sel_first_symbol(sc->symtab, SYMBOL_METHOD);
        while (symbol)
        {
            fprintf(stderr, "  %s\n", _gmx_sel_sym_name(symbol));
            symbol = _gmx_sel_next_symbol(symbol, SYMBOL_METHOD);
        }
    }
}
