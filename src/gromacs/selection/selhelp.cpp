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
 * Implements functions in selhelp.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "selhelp.h"

#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/onlinehelp/helptopic.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "selmethod.h"
#include "symrec.h"

namespace
{

struct CommonHelpText
{
    static const char        name[];
    static const char        title[];
    static const char* const text[];
};

const char        CommonHelpText::name[]  = "selections";
const char        CommonHelpText::title[] = "Selection syntax and usage";
const char* const CommonHelpText::text[]  = {
    "Selections are used to select atoms/molecules/residues for analysis.",
    "In contrast to traditional index files, selections can be dynamic, i.e.,",
    "select different atoms for different trajectory frames. The GROMACS",
    "manual contains a short introductory section to selections in the",
    "Analysis chapter, including suggestions on how to get familiar with",
    "selections if you are new to the concept. The subtopics listed below",
    "provide more details on the technical and syntactic aspects of",
    "selections.[PAR]",

    "Each analysis tool requires a different number of selections and the",
    "selections are interpreted differently. The general idea is still the",
    "same: each selection evaluates to a set of positions, where a position",
    "can be an atom position or center-of-mass or center-of-geometry of",
    "a set of atoms. The tool then uses these positions for its analysis to",
    "allow very flexible processing. Some analysis tools may have limitations",
    "on the types of selections allowed."
};

struct ArithmeticHelpText
{
    static const char        name[];
    static const char        title[];
    static const char* const text[];
};

const char        ArithmeticHelpText::name[]  = "arithmetic";
const char        ArithmeticHelpText::title[] = "Arithmetic expressions in selections";
const char* const ArithmeticHelpText::text[]  = {
    "Basic arithmetic evaluation is supported for numeric expressions.",
    "Supported operations are addition, subtraction, negation, multiplication,",
    "division, and exponentiation (using ^).",
    "Result of a division by zero or other illegal operations is undefined.",
};

struct CmdLineHelpText
{
    static const char        name[];
    static const char        title[];
    static const char* const text[];
};

const char        CmdLineHelpText::name[]  = "cmdline";
const char        CmdLineHelpText::title[] = "Specifying selections from command line";
const char* const CmdLineHelpText::text[]  = {
    "If no selections are provided on the command line, you are prompted to",
    "type the selections interactively (a pipe can also be used to provide",
    "the selections in this case for most tools). While this works well for",
    "testing, it is easier to provide the selections from the command line",
    "if they are complex or for scripting.[PAR]",

    "Each tool has different command-line arguments for specifying selections",
    "(see the help for the individual tools).",
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
    "The default groups are generated with the same logic as for",
    "non-selection tools.",
    "",
    "Depending on the tool, two additional command-line arguments may be",
    "available to control the behavior:",
    "",
    "* [TT]-seltype[tt] can be used to specify the default type of",
    "  positions to calculate for each selection.",
    "* [TT]-selrpos[tt] can be used to specify the default type of",
    "  positions used in selecting atoms by coordinates.",
    "",
    "See the \"positions\" subtopic for more information on these options.",
    "",
    "Tools that take selections apply them to a structure/topology and/or",
    "a trajectory file. If the tool takes both (typically as [TT]-s[tt]",
    "for structure/topology and [TT]-f[tt] for trajectory), then the",
    "trajectory file is only used for coordinate information, and all other",
    "information, such as atom names and residue information, is read from",
    "the structure/topology file. If the tool only takes a structure file,",
    "or if only that input parameter is provided, then also the coordinates",
    "are taken from that file.",
    "For example, to select atoms from a [TT].pdb[tt]/[TT].gro[tt] file in",
    "a tool that provides both options, pass it as [TT]-s[tt] (only).",
    "There is no warning if the trajectory file specifies, e.g., different",
    "atom names than the structure file. Only the number of atoms is checked.",
    "Many selection-enabled tools also provide an [TT]-fgroup[tt] option",
    "to specify the atom indices that are present in the trajectory for cases",
    "where the trajectory only has a subset of atoms from the",
    "topology/structure file."
};

struct EvaluationHelpText
{
    static const char        name[];
    static const char        title[];
    static const char* const text[];
};

const char        EvaluationHelpText::name[]  = "evaluation";
const char        EvaluationHelpText::title[] = "Selection evaluation and optimization";
const char* const EvaluationHelpText::text[]  = {
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
    "selections, in particular if you define several index groups like this::",
    "",
    "  rdist = distance from com of resnr 1 to 5;",
    "  resname RES and rdist < 2;",
    "  resname RES and rdist < 4;",
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
    static const char        name[];
    static const char        title[];
    static const char* const text[];
};

const char        ExamplesHelpText::name[]  = "examples";
const char        ExamplesHelpText::title[] = "Selection examples";
const char* const ExamplesHelpText::text[]  = {
    "Below, examples of different types of selections are given.",
    "",
    "* Selection of all water oxygens::",
    "",
    "    resname SOL and name OW",
    "",
    "* Centers of mass of residues 1 to 5 and 10::",
    "",
    "    res_com of resnr 1 to 5 10",
    "",
    "* All atoms farther than 1 nm of a fixed position::",
    "",
    "    not within 1 of [1.2, 3.1, 2.4]",
    "",
    "* All atoms of a residue LIG within 0.5 nm of a protein (with a custom name)::",
    "",
    "    \"Close to protein\" resname LIG and within 0.5 of group \"Protein\"",
    "",
    "* All protein residues that have at least one atom within 0.5 nm of a residue LIG::",
    "",
    "    group \"Protein\" and same residue as within 0.5 of resname LIG",
    "",
    "* All RES residues whose COM is between 2 and 4 nm from the COM of all of them::",
    "",
    "    rdist = res_com distance from com of resname RES",
    "    resname RES and rdist >= 2 and rdist <= 4",
    "",
    // TODO: Make it possible to use links below.
    "* Selection like with duplicate atoms like C1 C2 C2 C3 C3 C4 ... C8 C9::",
    "",
    "    name \"C[1-8]\" merge name \"C[2-9]\"",
    "",
    "  This can be used with [TT]gmx distance[tt] to compute C1-C2, C2-C3 etc.",
    "  distances.",
    "",
    "* Selection with atoms in order C2 C1::",
    "",
    "    name C1 C2 permute 2 1",
    "",
    "  This can be used with [TT]gmx gangle[tt] to get C2->C1 vectors instead of",
    "  C1->C2.",
    "",
    "* Selection with COMs of two index groups::",
    "",
    "    com of group 1 plus com of group 2",
    "",
    "  This can be used with [TT]gmx distance[tt] to compute the distance between",
    "  these two COMs.",
    "",
    "* Fixed vector along x (can be used as a reference with [TT]gmx gangle[tt])::",
    "",
    "    [0, 0, 0] plus [1, 0, 0]",
    "",
    "* The following examples explain the difference between the various",
    "  position types.  This selection selects a position for each residue",
    "  where any of the three atoms C[123] has [TT]x < 2[tt].  The positions",
    "  are computed as the COM of all three atoms.",
    "  This is the default behavior if you just write [TT]res_com of[tt]. ::",
    "",
    "    part_res_com of name C1 C2 C3 and x < 2",
    "",
    "  This selection does the same, but the positions are computed as COM",
    "  positions of whole residues::",
    "",
    "    whole_res_com of name C1 C2 C3 and x < 2",
    "",
    "  Finally, this selection selects the same residues, but the positions",
    "  are computed as COM of exactly those atoms atoms that match the",
    "  [TT]x < 2[tt] criterion::",
    "",
    "    dyn_res_com of name C1 C2 C3 and x < 2",
    "",
    "* Without the [TT]of[tt] keyword, the default behavior is different from",
    "  above, but otherwise the rules are the same::",
    "",
    "    name C1 C2 C3 and res_com x < 2",
    "",
    "  works as if [TT]whole_res_com[tt] was specified, and selects the three",
    "  atoms from residues whose COM satisfiex [TT]x < 2[tt].",
    "  Using ::",
    "",
    "    name C1 C2 C3 and part_res_com x < 2",
    "",
    "  instead selects residues based on the COM computed from the C[123] atoms.",
};

struct KeywordsHelpText
{
    static const char        name[];
    static const char        title[];
    static const char* const text[];
};

const char        KeywordsHelpText::name[]  = "keywords";
const char        KeywordsHelpText::title[] = "Selection keywords";
const char* const KeywordsHelpText::text[]  = {
    "The following selection keywords are currently available.",
    "For keywords marked with a plus, additional help is available through",
    "a subtopic KEYWORD, where KEYWORD is the name of the keyword.",
};

struct LimitationsHelpText
{
    static const char        name[];
    static const char        title[];
    static const char* const text[];
};

const char        LimitationsHelpText::name[]  = "limitations";
const char        LimitationsHelpText::title[] = "Selection limitations";
const char* const LimitationsHelpText::text[]  = {
    "* Some analysis programs may require a special structure for the input",
    "  selections (e.g., some options of [TT]gmx gangle[tt] require the index",
    "  group to be made of groups of three or four atoms).",
    "  For such programs, it is up to the user to provide a proper selection",
    "  expression that always returns such positions.",
    "",
    "* All selection keywords select atoms in increasing order, i.e., you can",
    "  consider them as set operations that in the end return the atoms in",
    "  sorted numerical order.  For example, the following selections select",
    "  the same atoms in the same order::",
    "",
    "    resname RA RB RC",
    "    resname RB RC RA",
    "",
    "  ::",
    "",
    "    atomnr 10 11 12 13",
    "    atomnr 12 13 10 11",
    "    atomnr 10 to 13",
    "    atomnr 13 to 10",
    "",
    "  If you need atoms/positions in a different order, you can:",
    "",
    "  * use external index groups (for some static selections),",
    "  * use the [TT]permute[tt] keyword to change the final order, or",
    "  * use the [TT]merge[tt] or [TT]plus[tt] keywords to compose the",
    "    final selection from multiple distinct selections.",
    "",
    "* Due to technical reasons, having a negative value as the first value in",
    "  expressions like ::",
    "",
    "    charge -1 to -0.7",
    "",
    "  result in a syntax error. A workaround is to write ::",
    "",
    "    charge {-1 to -0.7}",
    "",
    "  instead.",
    "",
    "* When [TT]name[tt] selection keyword is used together with PDB input",
    "  files, the behavior may be unintuitive. When GROMACS reads in a PDB",
    "  file, 4 character atom names that start with a digit are transformed",
    "  such that, e.g., 1HG2 becomes HG21, and the latter is what is matched",
    "  by the [TT]name[tt] keyword. Use [TT]pdbname[tt] to match the atom name",
    "  as it appears in the input PDB file.",
};

struct PositionsHelpText
{
    static const char        name[];
    static const char        title[];
    static const char* const text[];
};

const char        PositionsHelpText::name[]  = "positions";
const char        PositionsHelpText::title[] = "Specifying positions in selections";
const char* const PositionsHelpText::text[]  = {
    "Possible ways of specifying positions in selections are:",
    "",
    "1. A constant position can be defined as [TT][XX, YY, ZZ][tt], where",
    "   [TT]XX[tt], [TT]YY[tt] and [TT]ZZ[tt] are real numbers.[PAR]",
    "",
    "2. [TT]com of ATOM_EXPR [pbc][tt] or [TT]cog of ATOM_EXPR [pbc][tt]",
    "   calculate the center of mass/geometry of [TT]ATOM_EXPR[tt]. If",
    "   [TT]pbc[tt] is specified, the center is calculated iteratively to try",
    "   to deal with cases where [TT]ATOM_EXPR[tt] wraps around periodic",
    "   boundary conditions.",
    "",
    "3. [TT]POSTYPE of ATOM_EXPR[tt] calculates the specified positions for",
    "   the atoms in [TT]ATOM_EXPR[tt].",
    "   [TT]POSTYPE[tt] can be [TT]atom[tt], [TT]res_com[tt], [TT]res_cog[tt],",
    "   [TT]mol_com[tt] or [TT]mol_cog[tt], with an optional prefix [TT]whole_[tt]",
    "   [TT]part_[tt] or [TT]dyn_[tt].",
    "   [TT]whole_[tt] calculates the centers for the whole residue/molecule,",
    "   even if only part of it is selected.",
    "   [TT]part_[tt] prefix calculates the centers for the selected atoms, but",
    "   uses always the same atoms for the same residue/molecule. The used atoms",
    "   are determined from the largest group allowed by the selection.",
    "   [TT]dyn_[tt] calculates the centers strictly only for the selected atoms.",
    "   If no prefix is specified, whole selections default to [TT]part_[tt] and",
    "   other places default to [TT]whole_[tt].",
    "   The latter is often desirable to select the same molecules in different",
    "   tools, while the first is a compromise between speed ([TT]dyn_[tt]",
    "   positions can be slower to evaluate than [TT]part_[tt]) and intuitive",
    "   behavior.",
    "",
    "4. [TT]ATOM_EXPR[tt], when given for whole selections, is handled as 3.",
    "   above, using the position type from the command-line argument",
    "   [TT]-seltype[tt].",
    "",
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
    static const char* const text[];
};

const char        SyntaxHelpText::name[]  = "syntax";
const char        SyntaxHelpText::title[] = "Selection syntax";
const char* const SyntaxHelpText::text[]  = {
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
    "A variable is defined with the following syntax::",
    "",
    "  VARNAME = EXPR ;",
    "",
    "where [TT]EXPR[tt] is any valid selection expression.",
    "After this, [TT]VARNAME[tt] can be used anywhere where [TT]EXPR[tt]",
    "would be valid.[PAR]",

    "Selections are composed of three main types of expressions, those that",
    "define atoms ([TT]ATOM_EXPR[tt]), those that define positions",
    "([TT]POS_EXPR[tt]), and those that evaluate to numeric values",
    "([TT]NUM_EXPR[tt]). Each selection should be a [TT]POS_EXPR[tt]",
    "or a [TT]ATOM_EXPR[tt] (the latter is automatically converted to",
    "positions). The basic rules are as follows:",
    "",
    "* An expression like [TT]NUM_EXPR1 < NUM_EXPR2[tt] evaluates to an",
    "  [TT]ATOM_EXPR[tt] that selects all the atoms for which the comparison",
    "  is true.",
    "* Atom expressions can be combined with boolean operations such as",
    "  [TT]not ATOM_EXPR[tt], [TT]ATOM_EXPR and ATOM_EXPR[tt], or",
    "  [TT]ATOM_EXPR or ATOM_EXPR[tt]. Parentheses can be used to alter the",
    "  evaluation order.",
    "* [TT]ATOM_EXPR[tt] expressions can be converted into [TT]POS_EXPR[tt]",
    "  expressions in various ways, see the \"positions\" subtopic for more",
    "  details.",
    "* [TT]POS_EXPR[tt] can be converted into [TT]NUM_EXPR[tt] using syntax",
    "  like \"[TT]x of POS_EXPR[tt]\". Currently, this is only supported for single",
    "  positions like in expression \"[TT]x of cog of ATOM_EXPR[tt]\".",
    "",

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
    KeywordDetailsHelpTopic(const std::string& name, const gmx_ana_selmethod_t& method) :
        name_(name), method_(method)
    {
    }

    const char* name() const override { return name_.c_str(); }
    const char* title() const override { return method_.help.helpTitle; }

protected:
    std::string helpText() const override
    {
        return joinStrings(method_.help.help, method_.help.help + method_.help.nlhelp, "\n");
    }

private:
    std::string                name_;
    const gmx_ana_selmethod_t& method_;

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

    void writeHelp(const HelpWriterContext& context) const override;

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
    typedef std::vector<std::pair<std::string, const gmx_ana_selmethod_t*>> MethodList;

    /*! \brief
     * Prints markup for starting a list of keywords.
     */
    static void writeKeywordListStart(const HelpWriterContext& context, const char* heading);
    /*! \brief
     * Prints markup for ending a list of keywords.
     */
    static void writeKeywordListEnd(const HelpWriterContext& context, const char* extraInfo);

    /*! \brief
     * Prints a brief list of keywords (selection methods) available.
     *
     * \param[in] context  Context for printing the help.
     * \param[in] type     Only methods that return this type are printed.
     * \param[in] bModifiers  If false, \ref SMETH_MODIFIER methods are
     *      excluded, otherwise only them are printed.
     */
    void printKeywordList(const HelpWriterContext& context, e_selvalue_t type, bool bModifiers) const;

    /*! \brief
     * Prints the detailed help for keywords for rst export.
     */
    void writeKeywordSubTopics(const HelpWriterContext& context) const;

    MethodList methods_;
};

KeywordsHelpTopic::KeywordsHelpTopic()
{
    // TODO: This is not a very elegant way of getting the list of selection
    // methods, but this needs to be rewritten in any case if/when #652 is
    // implemented.
    const std::unique_ptr<SelectionParserSymbolTable> symtab(new SelectionParserSymbolTable);
    gmx_ana_selmethod_register_defaults(symtab.get());

    SelectionParserSymbolIterator symbol = symtab->beginIterator(SelectionParserSymbol::MethodSymbol);
    while (symbol != symtab->endIterator())
    {
        const std::string&         symname = symbol->name();
        const gmx_ana_selmethod_t* method  = symbol->methodValue();
        methods_.push_back(std::make_pair(std::string(symname), method));
        if (method->help.nlhelp > 0 && method->help.help != nullptr)
        {
            addSubTopic(HelpTopicPointer(new KeywordDetailsHelpTopic(symname, *method)));
        }
        ++symbol;
    }
}

void KeywordsHelpTopic::writeHelp(const HelpWriterContext& context) const
{
    context.writeTextBlock(helpText());

    // Print the list of keywords
    writeKeywordListStart(context, "Keywords that select atoms by an integer property:");
    printKeywordList(context, INT_VALUE, false);
    writeKeywordListEnd(context, "(use in expressions or like \"atomnr 1 to 5 7 9\")");

    writeKeywordListStart(context, "Keywords that select atoms by a numeric property:");
    printKeywordList(context, REAL_VALUE, false);
    writeKeywordListEnd(context, "(use in expressions or like \"occupancy 0.5 to 1\")");

    writeKeywordListStart(context, "Keywords that select atoms by a string property:");
    printKeywordList(context, STR_VALUE, false);
    writeKeywordListEnd(context, "(use like \"name PATTERN [PATTERN] ...\")");

    writeKeywordListStart(context, "Additional keywords that directly select atoms:");
    printKeywordList(context, GROUP_VALUE, false);
    writeKeywordListEnd(context, nullptr);

    writeKeywordListStart(context, "Keywords that directly evaluate to positions:");
    printKeywordList(context, POS_VALUE, false);
    writeKeywordListEnd(context, "(see also \"positions\" subtopic)");

    writeKeywordListStart(context, "Additional keywords:");
    printKeywordList(context, POS_VALUE, true);
    printKeywordList(context, NO_VALUE, true);
    writeKeywordListEnd(context, nullptr);

    writeKeywordSubTopics(context);
}

void KeywordsHelpTopic::writeKeywordListStart(const HelpWriterContext& context, const char* heading)
{
    context.paragraphBreak();
    std::string fullHeading("* ");
    fullHeading.append(heading);
    context.writeTextBlock(fullHeading);
    if (context.outputFormat() == eHelpOutputFormat_Rst)
    {
        context.paragraphBreak();
        context.writeTextBlock("  ::");
        context.paragraphBreak();
    }
}

void KeywordsHelpTopic::writeKeywordListEnd(const HelpWriterContext& context, const char* extraInfo)
{
    if (context.outputFormat() == eHelpOutputFormat_Rst)
    {
        context.paragraphBreak();
    }
    if (!isNullOrEmpty(extraInfo))
    {
        std::string fullInfo("  ");
        fullInfo.append(extraInfo);
        context.writeTextBlock(fullInfo);
    }
}

void KeywordsHelpTopic::printKeywordList(const HelpWriterContext& context, e_selvalue_t type, bool bModifiers) const
{
    TextWriter&                file = context.outputFile();
    MethodList::const_iterator iter;
    for (iter = methods_.begin(); iter != methods_.end(); ++iter)
    {
        const gmx_ana_selmethod_t& method      = *iter->second;
        const bool                 bIsModifier = (method.flags & SMETH_MODIFIER) != 0;
        if (method.type == type && bModifiers == bIsModifier)
        {
            const bool bHasHelp       = (method.help.nlhelp > 0 && method.help.help != nullptr);
            const bool bPrintHelpMark = bHasHelp && context.outputFormat() == eHelpOutputFormat_Console;
            file.writeString(formatString("   %c ", bPrintHelpMark ? '+' : ' '));
            if (method.help.syntax != nullptr)
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

void KeywordsHelpTopic::writeKeywordSubTopics(const HelpWriterContext& context) const
{
    if (context.outputFormat() != eHelpOutputFormat_Rst)
    {
        return;
    }
    std::set<std::string>      usedSymbols;
    MethodList::const_iterator iter;
    for (iter = methods_.begin(); iter != methods_.end(); ++iter)
    {
        const gmx_ana_selmethod_t& method = *iter->second;
        const bool bHasHelp               = (method.help.nlhelp > 0 && method.help.help != nullptr);
        if (!bHasHelp || usedSymbols.count(iter->first) > 0)
        {
            continue;
        }

        std::string title;
        if (method.help.helpTitle != nullptr)
        {
            title = method.help.helpTitle;
            title.append(" - ");
        }
        title.append(iter->first);
        MethodList::const_iterator mergeIter = iter;
        for (++mergeIter; mergeIter != methods_.end(); ++mergeIter)
        {
            if (mergeIter->second->help.help == method.help.help)
            {
                title.append(", ");
                title.append(mergeIter->first);
                usedSymbols.insert(mergeIter->first);
            }
        }

        const IHelpTopic* subTopic = findSubTopic(iter->first.c_str());
        GMX_RELEASE_ASSERT(subTopic != nullptr, "Keyword subtopic no longer exists");
        HelpWriterContext subContext(context);
        subContext.enterSubSection(title);
        subTopic->writeHelp(subContext);
    }
}

} // namespace

//! \cond libapi */
HelpTopicPointer createSelectionHelpTopic()
{
    CompositeHelpTopicPointer root(new CompositeHelpTopic<CommonHelpText>);
    root->registerSubTopic<SimpleHelpTopic<CmdLineHelpText>>();
    root->registerSubTopic<SimpleHelpTopic<SyntaxHelpText>>();
    root->registerSubTopic<SimpleHelpTopic<PositionsHelpText>>();
    root->registerSubTopic<SimpleHelpTopic<ArithmeticHelpText>>();
    root->registerSubTopic<KeywordsHelpTopic>();
    root->registerSubTopic<SimpleHelpTopic<EvaluationHelpText>>();
    root->registerSubTopic<SimpleHelpTopic<LimitationsHelpText>>();
    root->registerSubTopic<SimpleHelpTopic<ExamplesHelpText>>();
    return HelpTopicPointer(root.release());
}
//! \endcond

} // namespace gmx
