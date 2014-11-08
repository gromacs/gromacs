/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements gmx::SelectionCollection.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "selectioncollection.h"

#include <cctype>
#include <cstdio>

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "gromacs/fileio/trx.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/onlinehelp/helpmanager.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selhelp.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "compiler.h"
#include "mempool.h"
#include "parser.h"
#include "poscalc.h"
#include "scanner.h"
#include "selectioncollection-impl.h"
#include "selelem.h"
#include "selmethod.h"
#include "symrec.h"

namespace gmx
{

/********************************************************************
 * SelectionCollection::Impl
 */

SelectionCollection::Impl::Impl()
    : maxAtomIndex_(0), debugLevel_(0), bExternalGroupsSet_(false), grps_(NULL)
{
    sc_.nvars     = 0;
    sc_.varstrs   = NULL;
    sc_.top       = NULL;
    gmx_ana_index_clear(&sc_.gall);
    sc_.mempool   = NULL;
    sc_.symtab.reset(new SelectionParserSymbolTable);
    gmx_ana_selmethod_register_defaults(sc_.symtab.get());
}


SelectionCollection::Impl::~Impl()
{
    clearSymbolTable();
    // The tree must be freed before the SelectionData objects, since the
    // tree may hold references to the position data in SelectionData.
    sc_.root.reset();
    sc_.sel.clear();
    for (int i = 0; i < sc_.nvars; ++i)
    {
        sfree(sc_.varstrs[i]);
    }
    sfree(sc_.varstrs);
    gmx_ana_index_deinit(&sc_.gall);
    if (sc_.mempool)
    {
        _gmx_sel_mempool_destroy(sc_.mempool);
    }
}


void
SelectionCollection::Impl::clearSymbolTable()
{
    sc_.symtab.reset();
}


namespace
{

/*! \brief
 * Reads a single selection line from stdin.
 *
 * \param[in]  infile        File to read from (typically File::standardInput()).
 * \param[in]  bInteractive  Whether to print interactive prompts.
 * \param[out] line          The read line in stored here.
 * \returns true if something was read, false if at end of input.
 *
 * Handles line continuation, reading also the continuing line(s) in one call.
 */
bool promptLine(File *infile, bool bInteractive, std::string *line)
{
    if (bInteractive)
    {
        fprintf(stderr, "> ");
    }
    if (!infile->readLineWithTrailingSpace(line))
    {
        return false;
    }
    while (endsWith(*line, "\\\n"))
    {
        line->resize(line->length() - 2);
        if (bInteractive)
        {
            fprintf(stderr, "... ");
        }
        std::string buffer;
        // Return value ignored, buffer remains empty and works correctly
        // if there is nothing to read.
        infile->readLineWithTrailingSpace(&buffer);
        line->append(buffer);
    }
    if (endsWith(*line, "\n"))
    {
        line->resize(line->length() - 1);
    }
    else if (bInteractive)
    {
        fprintf(stderr, "\n");
    }
    return true;
}

/*! \brief
 * Helper function for tokenizing the input and pushing them to the parser.
 *
 * \param     scanner       Tokenizer data structure.
 * \param     parserState   Parser data structure.
 * \param[in] bInteractive  Whether to operate in interactive mode.
 *
 * Repeatedly reads tokens using \p scanner and pushes them to the parser with
 * \p parserState until there is no more input, or until enough input is given
 * (only in interactive mode).
 */
int runParserLoop(yyscan_t scanner, _gmx_sel_yypstate *parserState,
                  bool bInteractive)
{
    int status    = YYPUSH_MORE;
    int prevToken = 0;
    do
    {
        YYSTYPE value;
        YYLTYPE location;
        int     token = _gmx_sel_yylex(&value, &location, scanner);
        if (bInteractive)
        {
            if (token == 0)
            {
                break;
            }
            // Empty commands cause the interactive parser to print out
            // status information. This avoids producing those unnecessarily,
            // e.g., from "resname RA;;".
            if (prevToken == CMD_SEP && token == CMD_SEP)
            {
                continue;
            }
            prevToken = token;
        }
        status = _gmx_sel_yypush_parse(parserState, token, &value, &location, scanner);
    }
    while (status == YYPUSH_MORE);
    _gmx_sel_lexer_rethrow_exception_if_occurred(scanner);
    return status;
}

/*! \brief
 * Print current status in response to empty line in interactive input.
 *
 * \param[in] sc             Selection collection data structure.
 * \param[in] grps           Available index groups.
 * \param[in] firstSelection Index of first selection from this interactive
 *     session.
 * \param[in] maxCount       Maximum number of selections.
 * \param[in] context        Context to print for what the selections are for.
 * \param[in] bFirst         Whether this is the header that is printed before
 *     any user input.
 *
 * Prints the available index groups and currently provided selections.
 */
void printCurrentStatus(gmx_ana_selcollection_t *sc, gmx_ana_indexgrps_t *grps,
                        size_t firstSelection, int maxCount,
                        const std::string &context, bool bFirst)
{
    if (grps != NULL)
    {
        std::fprintf(stderr, "Available static index groups:\n");
        gmx_ana_indexgrps_print(stderr, grps, 0);
    }
    std::fprintf(stderr, "Specify ");
    if (maxCount < 0)
    {
        std::fprintf(stderr, "any number of selections");
    }
    else if (maxCount == 1)
    {
        std::fprintf(stderr, "a selection");
    }
    else
    {
        std::fprintf(stderr, "%d selections", maxCount);
    }
    std::fprintf(stderr, "%s%s:\n",
                 context.empty() ? "" : " ", context.c_str());
    std::fprintf(stderr,
                 "(one per line, <enter> for status/groups, 'help' for help%s)\n",
                 maxCount < 0 ? ", Ctrl-D to end" : "");
    if (!bFirst && (sc->nvars > 0 || sc->sel.size() > firstSelection))
    {
        std::fprintf(stderr, "Currently provided selections:\n");
        for (int i = 0; i < sc->nvars; ++i)
        {
            std::fprintf(stderr, "     %s\n", sc->varstrs[i]);
        }
        for (size_t i = firstSelection; i < sc->sel.size(); ++i)
        {
            std::fprintf(stderr, " %2d. %s\n",
                         static_cast<int>(i - firstSelection + 1),
                         sc->sel[i]->selectionText());
        }
        if (maxCount > 0)
        {
            const int remaining
                = maxCount - static_cast<int>(sc->sel.size() - firstSelection);
            std::fprintf(stderr, "(%d more selection%s required)\n",
                         remaining, remaining > 1 ? "s" : "");
        }
    }
}

/*! \brief
 * Prints selection help in interactive selection input.
 *
 * \param[in] sc    Selection collection data structure.
 * \param[in] line  Line of user input requesting help (starting with `help`).
 *
 * Initializes the selection help if not yet initialized, and finds the help
 * topic based on words on the input line.
 */
void printHelp(gmx_ana_selcollection_t *sc, const std::string &line)
{
    if (sc->rootHelp.get() == NULL)
    {
        sc->rootHelp = createSelectionHelpTopic();
    }
    HelpWriterContext context(&File::standardError(),
                              eHelpOutputFormat_Console);
    HelpManager       manager(*sc->rootHelp, context);
    try
    {
        std::vector<std::string>                 topic = splitString(line);
        std::vector<std::string>::const_iterator value;
        // First item in the list is the 'help' token.
        for (value = topic.begin() + 1; value != topic.end(); ++value)
        {
            manager.enterTopic(*value);
        }
    }
    catch (const InvalidInputError &ex)
    {
        fprintf(stderr, "%s\n", ex.what());
        return;
    }
    manager.writeCurrentTopic();
}

/*! \brief
 * Helper function that runs the parser once the tokenizer has been
 * initialized.
 *
 * \param[in,out] scanner Scanner data structure.
 * \param[in]     bStdIn  Whether to use a line-based reading
 *      algorithm designed for interactive input.
 * \param[in]     maxnr   Maximum number of selections to parse
 *      (if -1, parse as many as provided by the user).
 * \param[in]     context Context to print for what the selections are for.
 * \returns       Vector of parsed selections.
 * \throws        std::bad_alloc if out of memory.
 * \throws        InvalidInputError if there is a parsing error.
 *
 * Used internally to implement parseFromStdin(), parseFromFile() and
 * parseFromString().
 */
SelectionList runParser(yyscan_t scanner, bool bStdIn, int maxnr,
                        const std::string &context)
{
    boost::shared_ptr<void>  scannerGuard(scanner, &_gmx_sel_free_lexer);
    gmx_ana_selcollection_t *sc   = _gmx_sel_lexer_selcollection(scanner);
    gmx_ana_indexgrps_t     *grps = _gmx_sel_lexer_indexgrps(scanner);

    size_t                   oldCount = sc->sel.size();
    {
        boost::shared_ptr<_gmx_sel_yypstate> parserState(
                _gmx_sel_yypstate_new(), &_gmx_sel_yypstate_delete);
        if (bStdIn)
        {
            File       &stdinFile(File::standardInput());
            const bool  bInteractive = _gmx_sel_is_lexer_interactive(scanner);
            if (bInteractive)
            {
                printCurrentStatus(sc, grps, oldCount, maxnr, context, true);
            }
            std::string line;
            int         status;
            while (promptLine(&stdinFile, bInteractive, &line))
            {
                if (bInteractive)
                {
                    line = stripString(line);
                    if (line.empty())
                    {
                        printCurrentStatus(sc, grps, oldCount, maxnr, context, false);
                        continue;
                    }
                    if (startsWith(line, "help")
                        && (line[4] == 0 || std::isspace(line[4])))
                    {
                        printHelp(sc, line);
                        continue;
                    }
                }
                line.append("\n");
                _gmx_sel_set_lex_input_str(scanner, line.c_str());
                status = runParserLoop(scanner, parserState.get(), true);
                if (status != YYPUSH_MORE)
                {
                    // TODO: Check if there is more input, and issue an
                    // error/warning if some input was ignored.
                    goto early_termination;
                }
            }
            {
                YYLTYPE location;
                status = _gmx_sel_yypush_parse(parserState.get(), 0, NULL,
                                               &location, scanner);
            }
            // TODO: Remove added selections from the collection if parsing failed?
            _gmx_sel_lexer_rethrow_exception_if_occurred(scanner);
early_termination:
            GMX_RELEASE_ASSERT(status == 0,
                               "Parser errors should have resulted in an exception");
        }
        else
        {
            int status = runParserLoop(scanner, parserState.get(), false);
            GMX_RELEASE_ASSERT(status == 0,
                               "Parser errors should have resulted in an exception");
        }
    }
    scannerGuard.reset();
    int nr = sc->sel.size() - oldCount;
    if (maxnr > 0 && nr != maxnr)
    {
        std::string message
            = formatString("Too few selections provided; got %d, expected %d",
                           nr, maxnr);
        GMX_THROW(InvalidInputError(message));
    }

    SelectionList                     result;
    SelectionDataList::const_iterator i;
    result.reserve(nr);
    for (i = sc->sel.begin() + oldCount; i != sc->sel.end(); ++i)
    {
        result.push_back(Selection(i->get()));
    }
    return result;
}

/*! \brief
 * Checks that index groups have valid atom indices.
 *
 * \param[in]    root    Root of selection tree to process.
 * \param[in]    natoms  Maximum number of atoms that the selections are set
 *     to evaluate.
 * \param        errors  Object for reporting any error messages.
 * \throws std::bad_alloc if out of memory.
 *
 * Recursively checks the selection tree for index groups.
 * Each found group is checked that it only contains atom indices that match
 * the topology/maximum number of atoms set for the selection collection.
 * Any issues are reported to \p errors.
 */
void checkExternalGroups(const SelectionTreeElementPointer &root,
                         int                                natoms,
                         ExceptionInitializer              *errors)
{
    if (root->type == SEL_CONST && root->v.type == GROUP_VALUE)
    {
        try
        {
            root->checkIndexGroup(natoms);
        }
        catch (const UserInputError &)
        {
            errors->addCurrentExceptionAsNested();
        }
    }

    SelectionTreeElementPointer child = root->child;
    while (child)
    {
        checkExternalGroups(child, natoms, errors);
        child = child->next;
    }
}

}   // namespace


void SelectionCollection::Impl::resolveExternalGroups(
        const SelectionTreeElementPointer &root,
        ExceptionInitializer              *errors)
{

    if (root->type == SEL_GROUPREF)
    {
        try
        {
            root->resolveIndexGroupReference(grps_, sc_.gall.isize);
        }
        catch (const UserInputError &)
        {
            errors->addCurrentExceptionAsNested();
        }
    }

    SelectionTreeElementPointer child = root->child;
    while (child)
    {
        resolveExternalGroups(child, errors);
        root->flags |= (child->flags & SEL_UNSORTED);
        child        = child->next;
    }
}


/********************************************************************
 * SelectionCollection
 */

SelectionCollection::SelectionCollection()
    : impl_(new Impl)
{
}


SelectionCollection::~SelectionCollection()
{
}


void
SelectionCollection::initOptions(Options *options)
{
    const char * const debug_levels[]
        = { "no", "basic", "compile", "eval", "full" };

    bool bAllowNonAtomOutput = false;
    SelectionDataList::const_iterator iter;
    for (iter = impl_->sc_.sel.begin(); iter != impl_->sc_.sel.end(); ++iter)
    {
        const internal::SelectionData &sel = **iter;
        if (!sel.hasFlag(efSelection_OnlyAtoms))
        {
            bAllowNonAtomOutput = true;
        }
    }

    const char *const *postypes = PositionCalculationCollection::typeEnumValues;
    options->addOption(StringOption("selrpos")
                           .enumValueFromNullTerminatedArray(postypes)
                           .store(&impl_->rpost_).defaultValue(postypes[0])
                           .description("Selection reference positions"));
    if (bAllowNonAtomOutput)
    {
        options->addOption(StringOption("seltype")
                               .enumValueFromNullTerminatedArray(postypes)
                               .store(&impl_->spost_).defaultValue(postypes[0])
                               .description("Default selection output positions"));
    }
    else
    {
        impl_->spost_ = postypes[0];
    }
    GMX_RELEASE_ASSERT(impl_->debugLevel_ >= 0 && impl_->debugLevel_ <= 4,
                       "Debug level out of range");
    options->addOption(StringOption("seldebug").hidden(impl_->debugLevel_ == 0)
                           .enumValue(debug_levels)
                           .defaultValue(debug_levels[impl_->debugLevel_])
                           .storeEnumIndex(&impl_->debugLevel_)
                           .description("Print out selection trees for debugging"));
}


void
SelectionCollection::setReferencePosType(const char *type)
{
    GMX_RELEASE_ASSERT(type != NULL, "Cannot assign NULL position type");
    // Check that the type is valid, throw if it is not.
    e_poscalc_t  dummytype;
    int          dummyflags;
    PositionCalculationCollection::typeFromEnum(type, &dummytype, &dummyflags);
    impl_->rpost_ = type;
}


void
SelectionCollection::setOutputPosType(const char *type)
{
    GMX_RELEASE_ASSERT(type != NULL, "Cannot assign NULL position type");
    // Check that the type is valid, throw if it is not.
    e_poscalc_t  dummytype;
    int          dummyflags;
    PositionCalculationCollection::typeFromEnum(type, &dummytype, &dummyflags);
    impl_->spost_ = type;
}


void
SelectionCollection::setDebugLevel(int debugLevel)
{
    impl_->debugLevel_ = debugLevel;
}


void
SelectionCollection::setTopology(t_topology *top, int natoms)
{
    GMX_RELEASE_ASSERT(natoms > 0 || top != NULL,
                       "The number of atoms must be given if there is no topology");
    // Get the number of atoms from the topology if it is not given.
    if (natoms <= 0)
    {
        natoms = top->atoms.nr;
    }
    if (impl_->bExternalGroupsSet_)
    {
        ExceptionInitializer        errors("Invalid index group references encountered");
        SelectionTreeElementPointer root = impl_->sc_.root;
        while (root)
        {
            checkExternalGroups(root, natoms, &errors);
            root = root->next;
        }
        if (errors.hasNestedExceptions())
        {
            GMX_THROW(InconsistentInputError(errors));
        }
    }
    gmx_ana_selcollection_t *sc = &impl_->sc_;
    // Do this first, as it allocates memory, while the others don't throw.
    gmx_ana_index_init_simple(&sc->gall, natoms);
    sc->pcc.setTopology(top);
    sc->top = top;
}


void
SelectionCollection::setIndexGroups(gmx_ana_indexgrps_t *grps)
{
    GMX_RELEASE_ASSERT(grps == NULL || !impl_->bExternalGroupsSet_,
                       "Can only set external groups once or clear them afterwards");
    impl_->grps_               = grps;
    impl_->bExternalGroupsSet_ = true;

    ExceptionInitializer        errors("Invalid index group reference(s)");
    SelectionTreeElementPointer root = impl_->sc_.root;
    while (root)
    {
        impl_->resolveExternalGroups(root, &errors);
        root->checkUnsortedAtoms(true, &errors);
        root = root->next;
    }
    if (errors.hasNestedExceptions())
    {
        GMX_THROW(InconsistentInputError(errors));
    }
    for (size_t i = 0; i < impl_->sc_.sel.size(); ++i)
    {
        impl_->sc_.sel[i]->refreshName();
    }
}


bool
SelectionCollection::requiresTopology() const
{
    e_poscalc_t  type;
    int          flags;

    if (!impl_->rpost_.empty())
    {
        flags = 0;
        // Should not throw, because has been checked earlier.
        PositionCalculationCollection::typeFromEnum(impl_->rpost_.c_str(),
                                                    &type, &flags);
        if (type != POS_ATOM)
        {
            return true;
        }
    }
    if (!impl_->spost_.empty())
    {
        flags = 0;
        // Should not throw, because has been checked earlier.
        PositionCalculationCollection::typeFromEnum(impl_->spost_.c_str(),
                                                    &type, &flags);
        if (type != POS_ATOM)
        {
            return true;
        }
    }

    SelectionTreeElementPointer sel = impl_->sc_.root;
    while (sel)
    {
        if (_gmx_selelem_requires_top(*sel))
        {
            return true;
        }
        sel = sel->next;
    }
    return false;
}


SelectionList
SelectionCollection::parseFromStdin(int nr, bool bInteractive,
                                    const std::string &context)
{
    yyscan_t scanner;

    _gmx_sel_init_lexer(&scanner, &impl_->sc_, bInteractive, nr,
                        impl_->bExternalGroupsSet_,
                        impl_->grps_);
    return runParser(scanner, true, nr, context);
}


SelectionList
SelectionCollection::parseFromFile(const std::string &filename)
{

    try
    {
        yyscan_t scanner;
        File     file(filename, "r");
        // TODO: Exception-safe way of using the lexer.
        _gmx_sel_init_lexer(&scanner, &impl_->sc_, false, -1,
                            impl_->bExternalGroupsSet_,
                            impl_->grps_);
        _gmx_sel_set_lex_input_file(scanner, file.handle());
        return runParser(scanner, false, -1, std::string());
    }
    catch (GromacsException &ex)
    {
        ex.prependContext(formatString(
                                  "Error in parsing selections from file '%s'",
                                  filename.c_str()));
        throw;
    }
}


SelectionList
SelectionCollection::parseFromString(const std::string &str)
{
    yyscan_t scanner;

    _gmx_sel_init_lexer(&scanner, &impl_->sc_, false, -1,
                        impl_->bExternalGroupsSet_,
                        impl_->grps_);
    _gmx_sel_set_lex_input_str(scanner, str.c_str());
    return runParser(scanner, false, -1, std::string());
}


void
SelectionCollection::compile()
{
    if (impl_->sc_.top == NULL && requiresTopology())
    {
        GMX_THROW(InconsistentInputError("Selection requires topology information, but none provided"));
    }
    if (!impl_->bExternalGroupsSet_)
    {
        setIndexGroups(NULL);
    }
    if (impl_->debugLevel_ >= 1)
    {
        printTree(stderr, false);
    }

    SelectionCompiler compiler;
    compiler.compile(this);

    if (impl_->debugLevel_ >= 1)
    {
        std::fprintf(stderr, "\n");
        printTree(stderr, false);
        std::fprintf(stderr, "\n");
        impl_->sc_.pcc.printTree(stderr);
        std::fprintf(stderr, "\n");
    }
    impl_->sc_.pcc.initEvaluation();
    if (impl_->debugLevel_ >= 1)
    {
        impl_->sc_.pcc.printTree(stderr);
        std::fprintf(stderr, "\n");
    }

    // TODO: It would be nicer to associate the name of the selection option
    // (if available) to the error message.
    SelectionDataList::const_iterator iter;
    for (iter = impl_->sc_.sel.begin(); iter != impl_->sc_.sel.end(); ++iter)
    {
        const internal::SelectionData &sel = **iter;
        if (sel.hasFlag(efSelection_OnlyAtoms))
        {
            if (!sel.hasOnlyAtoms())
            {
                std::string message = formatString(
                            "Selection '%s' does not evaluate to individual atoms. "
                            "This is not allowed in this context.",
                            sel.selectionText());
                GMX_THROW(InvalidInputError(message));
            }
        }
        if (sel.hasFlag(efSelection_DisallowEmpty))
        {
            if (sel.posCount() == 0)
            {
                std::string message = formatString(
                            "Selection '%s' never matches any atoms.",
                            sel.selectionText());
                GMX_THROW(InvalidInputError(message));
            }
        }
    }
}


void
SelectionCollection::evaluate(t_trxframe *fr, t_pbc *pbc)
{
    if (fr->natoms <= impl_->maxAtomIndex_)
    {
        std::string message = formatString(
                    "Trajectory has less atoms (%d) than what is required for "
                    "evaluating the provided selections (atoms up to index %d "
                    "are required).", fr->natoms, impl_->maxAtomIndex_ + 1);
        GMX_THROW(InconsistentInputError(message));
    }
    impl_->sc_.pcc.initFrame();

    SelectionEvaluator evaluator;
    evaluator.evaluate(this, fr, pbc);

    if (impl_->debugLevel_ >= 3)
    {
        std::fprintf(stderr, "\n");
        printTree(stderr, true);
    }
}


void
SelectionCollection::evaluateFinal(int nframes)
{
    SelectionEvaluator evaluator;
    evaluator.evaluateFinal(this, nframes);
}


void
SelectionCollection::printTree(FILE *fp, bool bValues) const
{
    SelectionTreeElementPointer sel = impl_->sc_.root;
    while (sel)
    {
        _gmx_selelem_print_tree(fp, *sel, bValues, 0);
        sel = sel->next;
    }
}


void
SelectionCollection::printXvgrInfo(FILE *out, output_env_t oenv) const
{
    if (output_env_get_xvg_format(oenv) != exvgNONE)
    {
        const gmx_ana_selcollection_t &sc = impl_->sc_;
        std::fprintf(out, "# Selections:\n");
        for (int i = 0; i < sc.nvars; ++i)
        {
            std::fprintf(out, "#   %s\n", sc.varstrs[i]);
        }
        for (size_t i = 0; i < sc.sel.size(); ++i)
        {
            std::fprintf(out, "#   %s\n", sc.sel[i]->selectionText());
        }
        std::fprintf(out, "#\n");
    }
}

} // namespace gmx
