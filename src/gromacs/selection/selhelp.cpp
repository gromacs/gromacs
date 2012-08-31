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
 * Implements functions in selhelp.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#include <string>
#include <vector>
#include <utility>

#include <boost/scoped_ptr.hpp>

#include "gromacs/onlinehelp/helptopic.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/stringutil.h"

#include "selhelp.h"
#include "selmethod.h"
#include "symrec.h"

namespace
{

struct CommonHelpText
{
    static const char name[];
    static const char title[];
    static const char *const text[];
};

const char CommonHelpText::name[] = "selections";
const char CommonHelpText::title[] =
    "Selection syntax and usage";
const char *const CommonHelpText::text[] = {
    "Selections are used to select atoms/molecules/residues for analysis.",
    "In contrast to traditional index files, selections can be dynamic, i.e.,",
    "select different atoms for different trajectory frames.",
    "",
    "Each analysis tool requires a different number of selections and the",
    "selections are interpreted differently. The general idea is still the",
    "same: each selection evaluates to a set of positions, where a position",
    "can be an atom position or center-of-mass or center-of-geometry of",
    "a set of atoms. The tool then uses these positions for its analysis to",
    "allow very flexible processing. Some analysis tools may have limitations",
    "on the types of selections allowed.",
    "",
    "To get started with selections, run, e.g., ``[PROGRAM] select``",
    "without specifying selections on the command-line and use the interactive",
    "prompt to try out different selections.",
    "This tool provides output options that allow one to see what is actually",
    "selected by the given selections, and the interactive prompt reports",
    "syntax errors immediately, allowing one to try again.",
    "The subtopics listed below give more details on different aspects of",
    "selections.",
};

struct ArithmeticHelpText
{
    static const char name[];
    static const char title[];
    static const char *const text[];
};

const char ArithmeticHelpText::name[] = "arithmetic";
const char ArithmeticHelpText::title[] =
    "Arithmetic expressions in selections";
const char *const ArithmeticHelpText::text[] = {
    "Basic arithmetic evaluation is supported for numeric expressions.",
    "Supported operations are addition, subtraction, negation, multiplication,",
    "division, and exponentiation (using ^).",
    "Result of a division by zero or other illegal operations is undefined.",
};

struct CmdLineHelpText
{
    static const char name[];
    static const char title[];
    static const char *const text[];
};

const char CmdLineHelpText::name[] = "cmdline";
const char CmdLineHelpText::title[] =
    "Specifying selections from command line";
const char *const CmdLineHelpText::text[] = {
    "If no selections are provided on the command line, you are prompted to",
    "type the selections interactively (a pipe can also be used to provide",
    "the selections in this case for most tools). While this works well for",
    "testing, it is easier to provide the selections from the command line",
    "if they are complex or for scripting.",
    "",
    "Each tool has different command-line arguments for specifying selections",
    "(listed by ``[PROGRAM] help <tool>``).",
    "You can either pass a single string containing all selections (separated",
    "by semicolons), or multiple strings, each containing one selection.",
    "Note that you need to quote the selections to protect them from the",
    "shell.",
    "",
    "If you set a selection command-line argument, but do not provide any",
    "selections, you are prompted to type the selections for that argument",
    "interactively. This is useful if that selection argument is optional,",
    "in which case it is not normally prompted for.",
    "",
    "To provide selections from a file, use ``-sf file.dat`` in the place",
    "of the selection for a selection argument (e.g.,",
    "``-select -sf file.dat``). In general, the ``-sf`` argument reads",
    "selections from the provided file and assigns them to selection arguments",
    "that have been specified up to that point, but for which no selections",
    "have been provided.",
    "As a special case, ``-sf`` provided on its own, without preceding",
    "selection arguments, assigns the selections to all (yet unset) required",
    "selections (i.e., those that would be promted interactively if no",
    "selections are provided on the command line).",
    "",
    "To use groups from a traditional index file, use argument ``-n``",
    "to provide a file. See the \"syntax\" subtopic for how to use them.",
    "If this option is not provided, default groups are generated.",
    "The default groups are generated by reading selections from a file",
    "``defselection.dat``. If such a file is found in the current",
    "directory, it is used instead of the one provided by default.",
    "",
    "Depending on the tool, two additional command-line arguments may be",
    "available to control the behavior:",
    "",
    "1. ``-seltype`` can be used to specify the default type of",
    "positions to calculate for each selection.\n",
    "2. ``-selrpos`` can be used to specify the default type of",
    "positions used in selecting atoms by coordinates.",
    "",
    "See the \"positions\" subtopic for more information on these options.",
};

struct EvaluationHelpText
{
    static const char name[];
    static const char title[];
    static const char *const text[];
};

const char EvaluationHelpText::name[] = "evaluation";
const char EvaluationHelpText::title[] =
    "Selection evaluation and optimization";
const char *const EvaluationHelpText::text[] = {
    "Boolean evaluation proceeds from left to right and is short-circuiting",
    "i.e., as soon as it is known whether an atom will be selected, the",
    "remaining expressions are not evaluated at all.",
    "This can be used to optimize the selections: you should write the",
    "most restrictive and/or the most inexpensive expressions first in",
    "boolean expressions.",
    "The relative ordering between dynamic and static expressions does not",
    "matter: all static expressions are evaluated only once, before the first",
    "frame, and the result becomes the leftmost expression.",
    "",
    "Another point for optimization is in common subexpressions: they are not",
    "automatically recognized, but can be manually optimized by the use of",
    "variables. This can have a big impact on the performance of complex",
    "selections, in particular if you define several index groups like this::",
    "",
    "  rdist = distance from com of resnr 1 to 5;\n",
    "  resname RES and rdist < 2;\n",
    "  resname RES and rdist < 4;\n",
    "  resname RES and rdist < 6;",
    "",
    "Without the variable assignment, the distances would be evaluated three",
    "times, although they are exactly the same within each selection.",
    "Anything assigned into a variable becomes a common subexpression that",
    "is evaluated only once during a frame.",
    "Currently, in some cases the use of variables can actually lead to a small",
    "performance loss because of the checks necessary to determine for which",
    "atoms the expression has already been evaluated, but this should not be",
    "a major problem.",
};

struct ExamplesHelpText
{
    static const char name[];
    static const char title[];
    static const char *const text[];
};

const char ExamplesHelpText::name[] = "examples";
const char ExamplesHelpText::title[] =
    "Selection examples";
const char *const ExamplesHelpText::text[] = {
    // TODO: Once there are more tools available, use examples that invoke
    // tools and explain what the selections do in those tools.
    "Below, examples of increasingly complex selections are given.",
    "",
    "Selection of all water oxygens::",
    "",
    "  resname SOL and name OW",
    "",
    "Centers of mass of residues 1 to 5 and 10::",
    "",
    "  res_com of resnr 1 to 5 10",
    "",
    "All atoms farther than 1 nm of a fixed position::",
    "",
    "  not within 1 of (1.2, 3.1, 2.4)",
    "",
    "All atoms of a residue LIG within 0.5 nm of a protein (with a custom name)::",
    "",
    "  \"Close to protein\" resname LIG and within 0.5 of group \"Protein\"",
    "",
    "All protein residues that have at least one atom within 0.5 nm of a residue LIG::",
    "",
    "  group \"Protein\" and same residue as within 0.5 of resname LIG",
    "",
    "All RES residues whose COM is between 2 and 4 nm from the COM of all of them::",
    "",
    "  rdist = res_com distance from com of resname RES\n",
    "  resname RES and rdist >= 2 and rdist <= 4",
    "",
    "Selection like ``C1 C2 C2 C3 C3 C4 ... C8 C9`` (e.g., for ``g_bond``)::",
    "",
    "  name \"C[1-8]\" merge name \"C[2-9]\"",
};

struct KeywordsHelpText
{
    static const char name[];
    static const char title[];
    static const char *const text[];
};

const char KeywordsHelpText::name[] = "keywords";
const char KeywordsHelpText::title[] =
    "Selection keywords";
const char *const KeywordsHelpText::text[] = {
    "The following selection keywords are currently available.",
    "For keywords marked with a star, additional help is available through",
    "a subtopic KEYWORD, where KEYWORD is the name of the keyword.",
};

struct LimitationsHelpText
{
    static const char name[];
    static const char title[];
    static const char *const text[];
};

const char LimitationsHelpText::name[] = "limitations";
const char LimitationsHelpText::title[] =
    "Selection limitations";
const char *const LimitationsHelpText::text[] = {
    "Some analysis programs may require a special structure for the input",
    "selections (e.g., ``g_angle`` requires the index group to be made",
    "of groups of three or four atoms).",
    "For such programs, it is up to the user to provide a proper selection",
    "expression that always returns such positions.",
    "",
    "Due to technical reasons, having a negative value as the first value in",
    "expressions like ::",
    "",
    "  charge -1 to -0.7",
    "",
    "result in a syntax error. A workaround is to write ::",
    "",
    "  charge {-1 to -0.7}",
    "",
    "instead.",
};

struct PositionsHelpText
{
    static const char name[];
    static const char title[];
    static const char *const text[];
};

const char PositionsHelpText::name[] = "positions";
const char PositionsHelpText::title[] =
    "Specifying positions in selections";
const char *const PositionsHelpText::text[] = {
    "Possible ways of specifying positions in selections are:",
    "",
    "1. A constant position can be defined as ``[XX, YY, ZZ]``, where",
    "``XX``, ``YY`` and ``ZZ`` are real numbers.",
    "",
    "2. ``com of ATOM_EXPR [pbc]`` or ``cog of ATOM_EXPR [pbc]``",
    "calculate the center of mass/geometry of ``ATOM_EXPR``. If",
    "``pbc`` is specified, the center is calculated iteratively to try",
    "to deal with cases where ``ATOM_EXPR`` wraps around periodic",
    "boundary conditions.",
    "",
    "3. ``POSTYPE of ATOM_EXPR`` calculates the specified positions for",
    "the atoms in ``ATOM_EXPR``.",
    "``POSTYPE`` can be ``atom``, ``res_com``, ``res_cog``,",
    "``mol_com`` or ``mol_cog``, with an optional prefix ``whole_``",
    "``part_`` or ``dyn_``.",
    "``whole_`` calculates the centers for the whole residue/molecule,",
    "even if only part of it is selected.",
    "``part_`` prefix calculates the centers for the selected atoms, but",
    "uses always the same atoms for the same residue/molecule. The used atoms",
    "are determined from the the largest group allowed by the selection.",
    "``dyn_`` calculates the centers strictly only for the selected atoms.",
    "If no prefix is specified, whole selections default to ``part_`` and",
    "other places default to ``whole_``.",
    "The latter is often desirable to select the same molecules in different",
    "tools, while the first is a compromise between speed (``dyn_``",
    "positions can be slower to evaluate than ``part_``) and intuitive",
    "behavior.",
    "",
    "4. ``ATOM_EXPR``, when given for whole selections, is handled as 3.",
    "above, using the position type from the command-line argument",
    "``-seltype``.",
    "",
    "Selection keywords that select atoms based on their positions, such as",
    "``dist from``, use by default the positions defined by the",
    "``-selrpos`` command-line option.",
    "This can be overridden by prepending a ``POSTYPE`` specifier to the",
    "keyword. For example, ``res_com dist from POS`` evaluates the",
    "residue center of mass distances. In the example, all atoms of a residue",
    "are either selected or not, based on the single distance calculated.",
};

struct SyntaxHelpText
{
    static const char name[];
    static const char title[];
    static const char *const text[];
};

const char SyntaxHelpText::name[] = "syntax";
const char SyntaxHelpText::title[] =
    "Selection syntax";
const char *const SyntaxHelpText::text[] = {
    "A set of selections consists of one or more selections, separated by",
    "semicolons. Each selection defines a set of positions for the analysis.",
    "Each selection can also be preceded by a string that gives a name for",
    "the selection for use in, e.g., graph legends.",
    "If no name is provided, the string used for the selection is used",
    "automatically as the name.",
    "",
    "For interactive input, the syntax is slightly altered: line breaks can",
    "also be used to separate selections. \\ followed by a line break can",
    "be used to continue a line if necessary.",
    "Notice that the above only applies to real interactive input,",
    "not if you provide the selections, e.g., from a pipe.",
    "",
    "It is possible to use variables to store selection expressions.",
    "A variable is defined with the following syntax::",
    "",
    "  VARNAME = EXPR ;",
    "",
    "where ``EXPR`` is any valid selection expression.",
    "After this, ``VARNAME`` can be used anywhere where ``EXPR``",
    "would be valid.",
    "",
    "Selections are composed of three main types of expressions, those that",
    "define atoms (``ATOM_EXPR``), those that define positions",
    "(``POS_EXPR``), and those that evaluate to numeric values",
    "(``NUM_EXPR``). Each selection should be a ``POS_EXPR``",
    "or a ``ATOM_EXPR`` (the latter is automatically converted to",
    "positions). The basic rules are as follows:",
    "",
    "1. An expression like ``NUM_EXPR1 < NUM_EXPR2`` evaluates to an",
    "``ATOM_EXPR`` that selects all the atoms for which the comparison",
    "is true.\n",
    "2. Atom expressions can be combined with boolean operations such as",
    "``not ATOM_EXPR``, ``ATOM_EXPR and ATOM_EXPR``, or",
    "``ATOM_EXPR or ATOM_EXPR``. Parentheses can be used to alter the",
    "evaluation order.\n",
    "3. ``ATOM_EXPR`` expressions can be converted into ``POS_EXPR``",
    "expressions in various ways, see the \"positions\" subtopic for more",
    "details.",
    "",
    "Some keywords select atoms based on string values such as the atom name.",
    "For these keywords, it is possible to use wildcards (``name \"C*\"``)",
    "or regular expressions (e.g., ``resname \"R[AB]\"``).",
    "The match type is automatically guessed from the string: if it contains",
    "other characters than letters, numbers, '*', or '?', it is interpreted",
    "as a regular expression.",
    "To force the matching to use literal string matching, use",
    "''name = \"C*\"'' to match a literal C*.",
    "To force other type of matching, use '?' or '~' in place of '=' to force",
    "wildcard or regular expression matching, respectively.",
    "",
    "Strings that contain non-alphanumeric characters should be enclosed in",
    "double quotes as in the examples. For other strings, the quotes are",
    "optional, but if the value conflicts with a reserved keyword, a syntax",
    "error will occur. If your strings contain uppercase letters, this should",
    "not happen.",
    "",
    "Index groups provided with the ``-n`` command-line option or",
    "generated by default can be accessed with ``group NR`` or",
    "``group NAME``, where ``NR`` is a zero-based index of the group",
    "and ``NAME`` is part of the name of the desired group.",
    "The keyword ``group`` is optional if the whole selection is",
    "provided from an index group.",
    "To see a list of available groups in the interactive mode, press enter",
    "in the beginning of a line.",
};

} // namespace

namespace gmx
{

namespace
{

/*! \internal \brief
 * Help topic implementation for an individual selection method.
 *
 * \ingroup module_selection
 */
class KeywordDetailsHelpTopic : public AbstractSimpleHelpTopic
{
    public:
        //! Initialize help topic for the given selection method.
        explicit KeywordDetailsHelpTopic(const gmx_ana_selmethod_t &method)
            : method_(method)
        {
        }

        virtual const char *name() const
        {
            return method_.name;
        }
        virtual const char *title() const
        {
            return NULL;
        }

    protected:
        virtual std::string helpText() const
        {
            // FIXME: Reformat the selection method help texts.
            return concatenateStrings(method_.help.help, method_.help.nlhelp);
        }

    private:
        const gmx_ana_selmethod_t &method_;

        GMX_DISALLOW_COPY_AND_ASSIGN(KeywordDetailsHelpTopic);
};

/*! \internal \brief
 * Custom help topic for printing a list of selection keywords.
 *
 * \ingroup module_selection
 */
class KeywordsHelpTopic : public CompositeHelpTopic<KeywordsHelpText>
{
    public:
        KeywordsHelpTopic();

        virtual void writeHelp(const HelpWriterContext &context) const;

    private:
        /*! \brief
         * Container for known selection methods.
         *
         * The first item in the pair is the name of the selection method, and
         * the second points to the static data structure that describes the
         * method.
         * The name in the first item may differ from the name of the static
         * data structure if an alias is defined for that method.
         */
        typedef std::vector<std::pair<std::string,
                                      const gmx_ana_selmethod_t *> >
                MethodList;

        /*! \brief
         * Formats a brief list of keywords (selection methods) available.
         *
         * \param[in,out] text     Output is appended to this string.
         * \param[in]     type     Only methods of this type are printed.
         * \param[in]     bModifiers  If false, \ref SMETH_MODIFIER methods are
         *      excluded, otherwise only they are printed.
         */
        void formatKeywordList(std::string *result,
                               e_selvalue_t type, bool bModifiers) const;

        MethodList              methods_;
};

KeywordsHelpTopic::KeywordsHelpTopic()
{
    // TODO: This is not a very elegant way of getting the list of selection
    // methods, but this needs to be rewritten in any case if/when #652 is
    // implemented.
    boost::scoped_ptr<SelectionParserSymbolTable> symtab(
            new SelectionParserSymbolTable);
    gmx_ana_selmethod_register_defaults(symtab.get());

    SelectionParserSymbolIterator symbol
        = symtab->beginIterator(SelectionParserSymbol::MethodSymbol);
    while (symbol != symtab->endIterator())
    {
        const std::string &symname = symbol->name();
        const gmx_ana_selmethod_t *method = symbol->methodValue();
        methods_.push_back(std::make_pair(std::string(symname), method));
        if (method->help.nlhelp > 0 && method->help.help != NULL)
        {
            addSubTopic(HelpTopicPointer(new KeywordDetailsHelpTopic(*method)));
        }
        ++symbol;
    }
}

void KeywordsHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    context.writeTextBlock(helpText());
    context.writeTextBlock("\n");

    std::string text;
    // Print the list of keywords
    text = "| Keywords that select atoms by an integer property:\n"
           "| (use in expressions or like \"atomnr 1 to 5 7 9\")\n\n::\n\n";
    formatKeywordList(&text, INT_VALUE, false);
    text.append("\n");
    context.writeTextBlock(text);

    text = "| Keywords that select atoms by a numeric property:\n"
           "| (use in expressions or like \"occupancy 0.5 to 1\")\n\n::\n\n";
    formatKeywordList(&text, REAL_VALUE, false);
    text.append("\n");
    context.writeTextBlock(text);

    text = "| Keywords that select atoms by a string property:\n"
           "| (use like \"name PATTERN [PATTERN] ...\")\n\n::\n\n";
    formatKeywordList(&text, STR_VALUE, false);
    text.append("\n");
    context.writeTextBlock(text);

    text = "Additional keywords that directly select atoms::\n\n";
    formatKeywordList(&text, GROUP_VALUE, false);
    text.append("\n");
    context.writeTextBlock(text);

    text = "| Keywords that directly evaluate to positions:\n"
           "| (see also \"positions\" subtopic)\n\n::\n\n";
    formatKeywordList(&text, POS_VALUE, false);
    text.append("\n");
    context.writeTextBlock(text);

    text = "Additional keywords::\n\n";
    formatKeywordList(&text, POS_VALUE, true);
    formatKeywordList(&text, NO_VALUE, true);
    context.writeTextBlock(text);

    // TODO: For exported help, the extra help for keywords should also be
    // printed here.  It may be better to strip the stars from the above lists
    // in that case.
}

void KeywordsHelpTopic::formatKeywordList(
        std::string *result, e_selvalue_t type, bool bModifiers) const
{
    MethodList::const_iterator iter;
    for (iter = methods_.begin(); iter != methods_.end(); ++iter)
    {
        const gmx_ana_selmethod_t &method = *iter->second;
        bool bIsModifier = (method.flags & SMETH_MODIFIER) != 0;
        if (method.type == type && bModifiers == bIsModifier)
        {
            bool bHasHelp = (method.help.nlhelp > 0 && method.help.help != NULL);
            result->append("  ");
            if (method.help.syntax != NULL)
            {
                result->append(method.help.syntax);
            }
            else
            {
                std::string symname = iter->first;
                result->append(symname);
                if (symname != method.name)
                {
                    result->append(formatString(" (synonym for %s)", method.name));
                }
            }
            if (bHasHelp)
            {
                result->append(" (*)");
            }
            result->append("\n");
        }
    }
}

} // namespace

/*! \cond internal */
HelpTopicPointer createSelectionHelpTopic()
{
    CompositeHelpTopicPointer root(new CompositeHelpTopic<CommonHelpText>);
    root->registerSubTopic<SimpleHelpTopic<ArithmeticHelpText> >();
    root->registerSubTopic<SimpleHelpTopic<CmdLineHelpText> >();
    root->registerSubTopic<SimpleHelpTopic<EvaluationHelpText> >();
    root->registerSubTopic<SimpleHelpTopic<ExamplesHelpText> >();
    root->registerSubTopic<KeywordsHelpTopic>();
    root->registerSubTopic<SimpleHelpTopic<LimitationsHelpText> >();
    root->registerSubTopic<SimpleHelpTopic<PositionsHelpText> >();
    root->registerSubTopic<SimpleHelpTopic<SyntaxHelpText> >();
    return move(root);
}
//! \endcond

} // namespace gmx
