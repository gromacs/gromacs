/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * \brief
 * Implements functions in selhelp.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
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
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

const char        CommonHelpText::name[]  = "selections";
const char        CommonHelpText::title[] =
    "Selection syntax and usage";
const char *const CommonHelpText::text[] = {
    "Selections are used to select atoms/molecules/residues for analysis.",
    "In contrast to traditional index files, selections can be dynamic, i.e.,",
    "select different atoms for different trajectory frames.[PAR]",

    "Each analysis tool requires a different number of selections and the",
    "selections are interpreted differently. The general idea is still the",
    "same: each selection evaluates to a set of positions, where a position",
    "can be an atom position or center-of-mass or center-of-geometry of",
    "a set of atoms. The tool then uses these positions for its analysis to",
    "allow very flexible processing. Some analysis tools may have limitations",
    "on the types of selections allowed.[PAR]",

    "To get started with selections, run, e.g., [TT][PROGRAM] select[tt]",
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
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

const char        ArithmeticHelpText::name[]  = "arithmetic";
const char        ArithmeticHelpText::title[] =
    "Arithmetic expressions in selections";
const char *const ArithmeticHelpText::text[] = {
    "Basic arithmetic evaluation is supported for numeric expressions.",
    "Supported operations are addition, subtraction, negation, multiplication,",
    "division, and exponentiation (using ^).",
    "Result of a division by zero or other illegal operations is undefined.",
};

struct CmdLineHelpText
{
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

const char        CmdLineHelpText::name[]  = "cmdline";
const char        CmdLineHelpText::title[] =
    "Specifying selections from command line";
const char *const CmdLineHelpText::text[] = {
    "If no selections are provided on the command line, you are prompted to",
    "type the selections interactively (a pipe can also be used to provide",
    "the selections in this case for most tools). While this works well for",
    "testing, it is easier to provide the selections from the command line",
    "if they are complex or for scripting.[PAR]",

    "Each tool has different command-line arguments for specifying selections",
    "(listed by [TT][PROGRAM] help <tool>[tt]).",
    "You can either pass a single string containing all selections (separated",
    "by semicolons), or multiple strings, each containing one selection.",
    "Note that you need to quote the selections to protect them from the",
    "shell.[PAR]",

    "If you set a selection command-line argument, but do not provide any",
    "selections, you are prompted to type the selections for that argument",
    "interactively. This is useful if that selection argument is optional,",
    "in which case it is not normally prompted for.[PAR]",

    "To provide selections from a file, use [TT]-sf file.dat[tt] in the place",
    "of the selection for a selection argument (e.g.,",
    "[TT]-select -sf file.dat[tt]). In general, the [TT]-sf[tt] argument reads",
    "selections from the provided file and assigns them to selection arguments",
    "that have been specified up to that point, but for which no selections",
    "have been provided.",
    "As a special case, [TT]-sf[tt] provided on its own, without preceding",
    "selection arguments, assigns the selections to all (yet unset) required",
    "selections (i.e., those that would be promted interactively if no",
    "selections are provided on the command line).[PAR]",

    "To use groups from a traditional index file, use argument [TT]-n[tt]",
    "to provide a file. See the \"syntax\" subtopic for how to use them.",
    "If this option is not provided, default groups are generated.",
    "The default groups are generated by reading selections from a file",
    "[TT]defselection.dat[tt]. If such a file is found in the current",
    "directory, it is used instead of the one provided by default.[PAR]",

    "Depending on the tool, two additional command-line arguments may be",
    "available to control the behavior:[BR]",
    "1. [TT]-seltype[tt] can be used to specify the default type of",
    "positions to calculate for each selection.[BR]",
    "2. [TT]-selrpos[tt] can be used to specify the default type of",
    "positions used in selecting atoms by coordinates.[BR]",
    "See the \"positions\" subtopic for more information on these options.",
};

struct EvaluationHelpText
{
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

const char        EvaluationHelpText::name[]  = "evaluation";
const char        EvaluationHelpText::title[] =
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

struct ExamplesHelpText
{
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

const char        ExamplesHelpText::name[]  = "examples";
const char        ExamplesHelpText::title[] =
    "Selection examples";
const char *const ExamplesHelpText::text[] = {
    // TODO: Once there are more tools available, use examples that invoke
    // tools and explain what the selections do in those tools.
    "Below, examples of increasingly complex selections are given.[PAR]",

    "Selection of all water oxygens:[BR]",
    "  resname SOL and name OW",
    "[PAR]",

    "Centers of mass of residues 1 to 5 and 10:[BR]",
    "  res_com of resnr 1 to 5 10",
    "[PAR]",

    "All atoms farther than 1 nm of a fixed position:[BR]",
    "  not within 1 of [1.2, 3.1, 2.4]",
    "[PAR]",

    "All atoms of a residue LIG within 0.5 nm of a protein (with a custom name):[BR]",
    "  \"Close to protein\" resname LIG and within 0.5 of group \"Protein\"",
    "[PAR]",

    "All protein residues that have at least one atom within 0.5 nm of a residue LIG:[BR]",
    "  group \"Protein\" and same residue as within 0.5 of resname LIG",
    "[PAR]",

    "All RES residues whose COM is between 2 and 4 nm from the COM of all of them:[BR]",
    "  rdist = res_com distance from com of resname RES[BR]",
    "  resname RES and rdist >= 2 and rdist <= 4",
    "[PAR]",

    "Selection like C1 C2 C2 C3 C3 C4 ... C8 C9 (e.g., for g_bond):[BR]",
    "  name \"C[1-8]\" merge name \"C[2-9]\"",
};

struct KeywordsHelpText
{
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

const char        KeywordsHelpText::name[]  = "keywords";
const char        KeywordsHelpText::title[] =
    "Selection keywords";
const char *const KeywordsHelpText::text[] = {
    "The following selection keywords are currently available.",
    "For keywords marked with a star, additional help is available through",
    "a subtopic KEYWORD, where KEYWORD is the name of the keyword.",
};

struct LimitationsHelpText
{
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

const char        LimitationsHelpText::name[]  = "limitations";
const char        LimitationsHelpText::title[] =
    "Selection limitations";
const char *const LimitationsHelpText::text[] = {
    "Some analysis programs may require a special structure for the input",
    "selections (e.g., [TT]gmx angle[tt] requires the index group to be made",
    "of groups of three or four atoms).",
    "For such programs, it is up to the user to provide a proper selection",
    "expression that always returns such positions.",
    "[PAR]",

    "Due to technical reasons, having a negative value as the first value in",
    "expressions like[BR]",
    "[TT]charge -1 to -0.7[tt][BR]",
    "result in a syntax error. A workaround is to write[BR]",
    "[TT]charge {-1 to -0.7}[tt][BR]",
    "instead.[PAR]",

    "When [TT]name[tt] selection keyword is used together with PDB input",
    "files, the behavior may be unintuitive. When Gromacs reads in a PDB",
    "file, 4 character atom names that start with a digit are transformed",
    "such that, e.g., 1HG2 becomes HG21, and the latter is what is matched",
    "by the [TT]name[tt] keyword. Use [TT]pdbname[tt] to match the atom name",
    "as it appears in the input PDB file.",
};

struct PositionsHelpText
{
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

const char        PositionsHelpText::name[]  = "positions";
const char        PositionsHelpText::title[] =
    "Specifying positions in selections";
const char *const PositionsHelpText::text[] = {
    "Possible ways of specifying positions in selections are:[PAR]",

    "1. A constant position can be defined as [TT][XX, YY, ZZ][tt], where",
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
    "[TT]whole_[tt] calculates the centers for the whole residue/molecule,",
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

struct SyntaxHelpText
{
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

const char        SyntaxHelpText::name[]  = "syntax";
const char        SyntaxHelpText::title[] =
    "Selection syntax";
const char *const SyntaxHelpText::text[] = {
    "A set of selections consists of one or more selections, separated by",
    "semicolons. Each selection defines a set of positions for the analysis.",
    "Each selection can also be preceded by a string that gives a name for",
    "the selection for use in, e.g., graph legends.",
    "If no name is provided, the string used for the selection is used",
    "automatically as the name.[PAR]",

    "For interactive input, the syntax is slightly altered: line breaks can",
    "also be used to separate selections. \\ followed by a line break can",
    "be used to continue a line if necessary.",
    "Notice that the above only applies to real interactive input,",
    "not if you provide the selections, e.g., from a pipe.[PAR]",

    "It is possible to use variables to store selection expressions.",
    "A variable is defined with the following syntax:[BR]",
    "[TT]VARNAME = EXPR ;[tt][BR]",
    "where [TT]EXPR[tt] is any valid selection expression.",
    "After this, [TT]VARNAME[tt] can be used anywhere where [TT]EXPR[tt]",
    "would be valid.[PAR]",

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
    "expressions in various ways, see the \"positions\" subtopic for more",
    "details.[PAR]",

    "Some keywords select atoms based on string values such as the atom name.",
    "For these keywords, it is possible to use wildcards ([TT]name \"C*\"[tt])",
    "or regular expressions (e.g., [TT]resname \"R[AB]\"[tt]).",
    "The match type is automatically guessed from the string: if it contains",
    "other characters than letters, numbers, '*', or '?', it is interpreted",
    "as a regular expression.",
    "To force the matching to use literal string matching, use",
    "[TT]name = \"C*\"[tt] to match a literal C*.",
    "To force other type of matching, use '?' or '~' in place of '=' to force",
    "wildcard or regular expression matching, respectively.[PAR]",

    "Strings that contain non-alphanumeric characters should be enclosed in",
    "double quotes as in the examples. For other strings, the quotes are",
    "optional, but if the value conflicts with a reserved keyword, a syntax",
    "error will occur. If your strings contain uppercase letters, this should",
    "not happen.[PAR]",

    "Index groups provided with the [TT]-n[tt] command-line option or",
    "generated by default can be accessed with [TT]group NR[tt] or",
    "[TT]group NAME[tt], where [TT]NR[tt] is a zero-based index of the group",
    "and [TT]NAME[tt] is part of the name of the desired group.",
    "The keyword [TT]group[tt] is optional if the whole selection is",
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
        KeywordDetailsHelpTopic(const std::string         &name,
                                const gmx_ana_selmethod_t &method)
            : name_(name), method_(method)
        {
        }

        virtual const char *name() const
        {
            return name_.c_str();
        }
        virtual const char *title() const
        {
            return NULL;
        }

    protected:
        virtual std::string helpText() const
        {
            return concatenateStrings(method_.help.help, method_.help.nlhelp);
        }

    private:
        std::string                name_;
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
         * Prints a brief list of keywords (selection methods) available.
         *
         * \param[in] context  Context for printing the help.
         * \param[in] type     Only methods that return this type are printed.
         * \param[in] bModifiers  If false, \ref SMETH_MODIFIER methods are
         *      excluded, otherwise only them are printed.
         */
        void printKeywordList(const HelpWriterContext &context,
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
        const std::string         &symname = symbol->name();
        const gmx_ana_selmethod_t *method  = symbol->methodValue();
        methods_.push_back(std::make_pair(std::string(symname), method));
        if (method->help.nlhelp > 0 && method->help.help != NULL)
        {
            addSubTopic(HelpTopicPointer(
                                new KeywordDetailsHelpTopic(symname, *method)));
        }
        ++symbol;
    }
}

void KeywordsHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        GMX_THROW(NotImplementedError(
                          "Selection help is not implemented for this output format"));
    }
    // TODO: The markup here is not really appropriate, and printKeywordList()
    // still prints raw text, but these are waiting for discussion of the
    // markup format in #969.
    writeBasicHelpTopic(context, *this, helpText());
    context.writeTextBlock("[BR]");

    // Print the list of keywords
    context.writeTextBlock(
            "Keywords that select atoms by an integer property:[BR]"
            "(use in expressions or like \"atomnr 1 to 5 7 9\")[BR]");
    printKeywordList(context, INT_VALUE, false);
    context.writeTextBlock("[BR]");

    context.writeTextBlock(
            "Keywords that select atoms by a numeric property:[BR]"
            "(use in expressions or like \"occupancy 0.5 to 1\")[BR]");
    printKeywordList(context, REAL_VALUE, false);
    context.writeTextBlock("[BR]");

    context.writeTextBlock(
            "Keywords that select atoms by a string property:[BR]"
            "(use like \"name PATTERN [PATTERN] ...\")[BR]");
    printKeywordList(context, STR_VALUE, false);
    context.writeTextBlock("[BR]");

    context.writeTextBlock(
            "Additional keywords that directly select atoms:[BR]");
    printKeywordList(context, GROUP_VALUE, false);
    context.writeTextBlock("[BR]");

    context.writeTextBlock(
            "Keywords that directly evaluate to positions:[BR]"
            "(see also \"positions\" subtopic)[BR]");
    printKeywordList(context, POS_VALUE, false);
    context.writeTextBlock("[BR]");

    context.writeTextBlock("Additional keywords:[BR]");
    printKeywordList(context, POS_VALUE, true);
    printKeywordList(context, NO_VALUE, true);
}

void KeywordsHelpTopic::printKeywordList(const HelpWriterContext &context,
                                         e_selvalue_t             type,
                                         bool                     bModifiers) const
{
    File &file = context.outputFile();
    MethodList::const_iterator iter;
    for (iter = methods_.begin(); iter != methods_.end(); ++iter)
    {
        const gmx_ana_selmethod_t &method = *iter->second;
        bool bIsModifier                  = (method.flags & SMETH_MODIFIER) != 0;
        if (method.type == type && bModifiers == bIsModifier)
        {
            bool bHasHelp = (method.help.nlhelp > 0 && method.help.help != NULL);
            file.writeString(formatString(" %c ", bHasHelp ? '*' : ' '));
            if (method.help.syntax != NULL)
            {
                file.writeLine(method.help.syntax);
            }
            else
            {
                std::string symname = iter->first;
                if (symname != method.name)
                {
                    symname.append(formatString(" (synonym for %s)", method.name));
                }
                file.writeLine(symname);
            }
        }
    }
}

}   // namespace

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
