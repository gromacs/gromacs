/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
 * Implements gmx::SelectionCollection.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/selection/selectioncollection.h"

#include <cctype>
#include <cstdio>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "gromacs/onlinehelp/helpmanager.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/onlinehelp/ihelptopic.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/parsetree.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionenums.h"
#include "gromacs/selection/selhelp.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textstream.h"
#include "gromacs/utility/textwriter.h"

#include "compiler.h"
#include "mempool.h"
#include "parser.h"
#include "poscalc.h"
#include "scanner.h"
#include "selectioncollection_impl.h"
#include "selelem.h"
#include "selmethod.h"
#include "symrec.h"

struct gmx_ana_indexgrps_t;
struct t_pbc;

namespace gmx
{

/********************************************************************
 * SelectionCollection::Impl
 */

SelectionCollection::Impl::Impl() :
    debugLevel_(DebugLevel::None), bExternalGroupsSet_(false), grps_(nullptr)
{
    sc_.nvars   = 0;
    sc_.varstrs = nullptr;
    sc_.top     = nullptr;
    gmx_ana_index_clear(&sc_.gall);
    sc_.mempool = nullptr;
    sc_.symtab  = std::make_unique<SelectionParserSymbolTable>();
    gmx_ana_index_clear(&requiredAtoms_);
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
    gmx_ana_index_deinit(&requiredAtoms_);
}


void SelectionCollection::Impl::clearSymbolTable()
{
    sc_.symtab.reset();
}


namespace
{

/*! \brief
 * Reads a single selection line from stdin.
 *
 * \param[in]  inputStream   Stream to read from (typically the StandardInputStream).
 * \param[in]  statusWriter  Stream to print prompts to (if NULL, no output is done).
 * \param[out] line          The read line in stored here.
 * \returns true if something was read, false if at end of input.
 *
 * Handles line continuation, reading also the continuing line(s) in one call.
 */
bool promptLine(TextInputStream* inputStream, TextWriter* statusWriter, std::string* line)
{
    if (statusWriter != nullptr)
    {
        statusWriter->writeString("> ");
    }
    if (!inputStream->readLine(line))
    {
        return false;
    }
    while (endsWith(*line, "\\\n"))
    {
        line->resize(line->length() - 2);
        if (statusWriter != nullptr)
        {
            statusWriter->writeString("... ");
        }
        std::string buffer;
        // Return value ignored, buffer remains empty and works correctly
        // if there is nothing to read.
        inputStream->readLine(&buffer);
        line->append(buffer);
    }
    if (endsWith(*line, "\n"))
    {
        line->resize(line->length() - 1);
    }
    else if (statusWriter != nullptr)
    {
        statusWriter->writeLine();
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
int runParserLoop(yyscan_t scanner, _gmx_sel_yypstate* parserState, bool bInteractive)
{
    int status = YYPUSH_MORE;
    do
    {
        YYSTYPE value;
        YYLTYPE location;
        int     token = _gmx_sel_yylex(&value, &location, scanner);
        if (bInteractive && token == 0)
        {
            break;
        }
        status = _gmx_sel_yypush_parse(parserState, token, &value, &location, scanner);
    } while (status == YYPUSH_MORE);
    _gmx_sel_lexer_rethrow_exception_if_occurred(scanner);
    return status;
}

/*! \brief
 * Print current status in response to empty line in interactive input.
 *
 * \param[in] writer         Writer to use for the output.
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
void printCurrentStatus(TextWriter*              writer,
                        gmx_ana_selcollection_t* sc,
                        gmx_ana_indexgrps_t*     grps,
                        size_t                   firstSelection,
                        int                      maxCount,
                        const std::string&       context,
                        bool                     bFirst)
{
    if (grps != nullptr)
    {
        writer->writeLine("Available static index groups:");
        gmx_ana_indexgrps_print(writer, grps, 0);
    }
    writer->writeString("Specify ");
    if (maxCount < 0)
    {
        writer->writeString("any number of selections");
    }
    else if (maxCount == 1)
    {
        writer->writeString("a selection");
    }
    else
    {
        writer->writeString(formatString("%d selections", maxCount));
    }
    writer->writeString(formatString("%s%s:\n", context.empty() ? "" : " ", context.c_str()));
    writer->writeString(
            formatString("(one per line, <enter> for status/groups, 'help' for help%s)\n",
                         maxCount < 0 ? ", Ctrl-D to end" : ""));
    if (!bFirst && (sc->nvars > 0 || sc->sel.size() > firstSelection))
    {
        writer->writeLine("Currently provided selections:");
        for (int i = 0; i < sc->nvars; ++i)
        {
            writer->writeString(formatString("     %s\n", sc->varstrs[i]));
        }
        for (size_t i = firstSelection; i < sc->sel.size(); ++i)
        {
            writer->writeString(formatString(
                    " %2d. %s\n", static_cast<int>(i - firstSelection + 1), sc->sel[i]->selectionText()));
        }
        if (maxCount > 0)
        {
            const int remaining = maxCount - static_cast<int>(sc->sel.size() - firstSelection);
            writer->writeString(formatString(
                    "(%d more selection%s required)\n", remaining, remaining > 1 ? "s" : ""));
        }
    }
}

/*! \brief
 * Prints selection help in interactive selection input.
 *
 * \param[in] writer Writer to use for the output.
 * \param[in] sc    Selection collection data structure.
 * \param[in] line  Line of user input requesting help (starting with `help`).
 *
 * Initializes the selection help if not yet initialized, and finds the help
 * topic based on words on the input line.
 */
void printHelp(TextWriter* writer, gmx_ana_selcollection_t* sc, const std::string& line)
{
    if (sc->rootHelp.get() == nullptr)
    {
        sc->rootHelp = createSelectionHelpTopic();
    }
    HelpWriterContext context(writer, eHelpOutputFormat_Console);
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
    catch (const InvalidInputError& ex)
    {
        writer->writeLine(ex.what());
        return;
    }
    manager.writeCurrentTopic();
}

/*! \brief
 * Helper function that runs the parser once the tokenizer has been
 * initialized.
 *
 * \param[in,out] scanner       Scanner data structure.
 * \param[in]     inputStream   Stream to use for input (currently only with
 *      `bInteractive==true`).
 * \param[in]     bInteractive  Whether to use a line-based reading
 *      algorithm designed for interactive input.
 * \param[in]     maxnr   Maximum number of selections to parse
 *      (if -1, parse as many as provided by the user).
 * \param[in]     context Context to print for what the selections are for.
 * \returns       Vector of parsed selections.
 * \throws        std::bad_alloc if out of memory.
 * \throws        InvalidInputError if there is a parsing error.
 *
 * Used internally to implement parseInteractive(), parseFromFile() and
 * parseFromString().
 */
SelectionList runParser(yyscan_t           scanner,
                        TextInputStream*   inputStream,
                        bool               bInteractive,
                        int                maxnr,
                        const std::string& context)
{
    std::shared_ptr<void>    scannerGuard(scanner, &_gmx_sel_free_lexer);
    gmx_ana_selcollection_t* sc   = _gmx_sel_lexer_selcollection(scanner);
    gmx_ana_indexgrps_t*     grps = _gmx_sel_lexer_indexgrps(scanner);

    size_t oldCount = sc->sel.size();
    {
        std::shared_ptr<_gmx_sel_yypstate> parserState(_gmx_sel_yypstate_new(), &_gmx_sel_yypstate_delete);
        if (bInteractive)
        {
            TextWriter* statusWriter = _gmx_sel_lexer_get_status_writer(scanner);
            if (statusWriter != nullptr)
            {
                printCurrentStatus(statusWriter, sc, grps, oldCount, maxnr, context, true);
            }
            std::string line;
            int         status;
            while (promptLine(inputStream, statusWriter, &line))
            {
                if (statusWriter != nullptr)
                {
                    line = stripString(line);
                    if (line.empty())
                    {
                        printCurrentStatus(statusWriter, sc, grps, oldCount, maxnr, context, false);
                        continue;
                    }
                    if (startsWith(line, "help") && (line[4] == 0 || (std::isspace(line[4]) != 0)))
                    {
                        printHelp(statusWriter, sc, line);
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
                status = _gmx_sel_yypush_parse(parserState.get(), 0, nullptr, &location, scanner);
            }
            // TODO: Remove added selections from the collection if parsing failed?
            _gmx_sel_lexer_rethrow_exception_if_occurred(scanner);
        early_termination:
            GMX_RELEASE_ASSERT(status == 0, "Parser errors should have resulted in an exception");
        }
        else
        {
            int status = runParserLoop(scanner, parserState.get(), false);
            GMX_RELEASE_ASSERT(status == 0, "Parser errors should have resulted in an exception");
        }
    }
    scannerGuard.reset();
    int nr = sc->sel.size() - oldCount;
    if (maxnr > 0 && nr != maxnr)
    {
        std::string message = formatString("Too few selections provided; got %d, expected %d", nr, maxnr);
        GMX_THROW(InvalidInputError(message));
    }

    SelectionList                     result;
    SelectionDataList::const_iterator i;
    result.reserve(nr);
    for (i = sc->sel.begin() + oldCount; i != sc->sel.end(); ++i)
    {
        result.emplace_back(i->get());
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
void checkExternalGroups(const SelectionTreeElementPointer& root, int natoms, ExceptionInitializer* errors)
{
    if (root->type == SEL_CONST && root->v.type == GROUP_VALUE)
    {
        try
        {
            root->checkIndexGroup(natoms);
        }
        catch (const UserInputError&)
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

//! Checks whether the given topology properties are available.
void checkTopologyProperties(const gmx_mtop_t* top, const SelectionTopologyProperties& props)
{
    if (top == nullptr)
    {
        if (props.hasAny())
        {
            GMX_THROW(InconsistentInputError(
                    "Selection requires topology information, but none provided"));
        }
        return;
    }
    if (props.needsMasses_ && !gmx_mtop_has_masses(top))
    {
        GMX_THROW(InconsistentInputError(
                "Selection requires mass information, but it is not available in the topology"));
    }
}

} // namespace


void SelectionCollection::Impl::resolveExternalGroups(const SelectionTreeElementPointer& root,
                                                      ExceptionInitializer*              errors)
{

    if (root->type == SEL_GROUPREF)
    {
        try
        {
            root->resolveIndexGroupReference(grps_, sc_.gall.isize);
        }
        catch (const UserInputError&)
        {
            errors->addCurrentExceptionAsNested();
        }
    }

    SelectionTreeElementPointer child = root->child;
    while (child)
    {
        resolveExternalGroups(child, errors);
        root->flags |= (child->flags & SEL_UNSORTED);
        child = child->next;
    }
}


bool SelectionCollection::Impl::areForcesRequested() const
{
    for (const auto& sel : sc_.sel)
    {
        if (sel->hasFlag(gmx::efSelection_EvaluateForces))
        {
            return true;
        }
    }
    return false;
}


SelectionTopologyProperties
SelectionCollection::Impl::requiredTopologyPropertiesForPositionType(const std::string& post, bool forces)
{
    SelectionTopologyProperties props;
    if (!post.empty())
    {
        switch (PositionCalculationCollection::requiredTopologyInfoForType(post.c_str(), forces))
        {
            case PositionCalculationCollection::RequiredTopologyInfo::None: break;
            case PositionCalculationCollection::RequiredTopologyInfo::Topology:
                props.merge(SelectionTopologyProperties::topology());
                break;
            case PositionCalculationCollection::RequiredTopologyInfo::TopologyAndMasses:
                props.merge(SelectionTopologyProperties::masses());
                break;
        }
    }
    return props;
}


/********************************************************************
 * SelectionCollection
 */

SelectionCollection::SelectionCollection() : impl_(new Impl) {}


SelectionCollection::~SelectionCollection() = default;

SelectionCollection::SelectionCollection(const SelectionCollection& rhs) :
    impl_(std::make_unique<Impl>())
{
    setReferencePosType(rhs.impl_->rpost_.empty() ? PositionCalculationCollection::typeEnumValues[0]
                                                  : rhs.impl_->rpost_.c_str());
    setOutputPosType(rhs.impl_->spost_.empty() ? PositionCalculationCollection::typeEnumValues[0]
                                               : rhs.impl_->spost_.c_str());
    setDebugLevel(static_cast<int>(rhs.impl_->debugLevel_));

    for (size_t i = 0; i < rhs.impl_->sc_.sel.size(); i++)
    {
        const auto& selectionOption = rhs.impl_->sc_.sel[i];
        parseFromString(selectionOption->selectionText());
        impl_->sc_.sel[i]->setFlags(selectionOption->flags());
    }

    // Topology has been initialized in rhs if top is non-null or natoms is set.
    // Note this needs to be set after selections are parsed to register topology requirements properly.
    if (rhs.impl_->sc_.top != nullptr || rhs.impl_->sc_.gall.isize > 0)
    {
        setTopology(rhs.impl_->sc_.top, rhs.impl_->sc_.gall.isize);
        gmx_ana_index_copy(&impl_->sc_.gall, &rhs.impl_->sc_.gall, /*balloc=*/false);
    }

    if (rhs.impl_->grps_ != nullptr)
    {
        setIndexGroups(rhs.impl_->grps_);
    }

    // Only compile the selection if rhs is compiled.
    if (rhs.impl_->sc_.mempool != nullptr)
    {
        compile();
    }
}

SelectionCollection& SelectionCollection::operator=(SelectionCollection rhs)
{
    rhs.swap(*this);
    return *this;
}

void SelectionCollection::swap(SelectionCollection& rhs) noexcept
{
    using std::swap;
    swap(impl_, rhs.impl_);
}

void SelectionCollection::initOptions(IOptionsContainer* options, SelectionTypeOption selectionTypeOption)
{
    static const EnumerationArray<Impl::DebugLevel, const char*> s_debugLevelNames = {
        { "no", "basic", "compile", "eval", "full" }
    };

    const char* const* postypes = PositionCalculationCollection::typeEnumValues;
    options->addOption(StringOption("selrpos")
                               .enumValueFromNullTerminatedArray(postypes)
                               .store(&impl_->rpost_)
                               .defaultValue(postypes[0])
                               .description("Selection reference positions"));
    if (selectionTypeOption == IncludeSelectionTypeOption)
    {
        options->addOption(StringOption("seltype")
                                   .enumValueFromNullTerminatedArray(postypes)
                                   .store(&impl_->spost_)
                                   .defaultValue(postypes[0])
                                   .description("Default selection output positions"));
    }
    else
    {
        impl_->spost_ = postypes[0];
    }
    GMX_RELEASE_ASSERT(impl_->debugLevel_ != Impl::DebugLevel::Count, "Debug level out of range");
    options->addOption(EnumOption<Impl::DebugLevel>("seldebug")
                               .hidden(impl_->debugLevel_ == Impl::DebugLevel::None)
                               .enumValue(s_debugLevelNames)
                               .store(&impl_->debugLevel_)
                               .description("Print out selection trees for debugging"));
}


void SelectionCollection::setReferencePosType(const char* type)
{
    GMX_RELEASE_ASSERT(type != nullptr, "Cannot assign NULL position type");
    // Check that the type is valid, throw if it is not.
    e_poscalc_t dummytype;
    int         dummyflags;
    PositionCalculationCollection::typeFromEnum(type, &dummytype, &dummyflags);
    impl_->rpost_ = type;
}


void SelectionCollection::setOutputPosType(const char* type)
{
    GMX_RELEASE_ASSERT(type != nullptr, "Cannot assign NULL position type");
    // Check that the type is valid, throw if it is not.
    e_poscalc_t dummytype;
    int         dummyflags;
    PositionCalculationCollection::typeFromEnum(type, &dummytype, &dummyflags);
    impl_->spost_ = type;
}


void SelectionCollection::setDebugLevel(int debugLevel)
{
    impl_->debugLevel_ = Impl::DebugLevel(debugLevel);
}


void SelectionCollection::setTopology(const gmx_mtop_t* top, int natoms)
{
    GMX_RELEASE_ASSERT(natoms > 0 || top != nullptr,
                       "The number of atoms must be given if there is no topology");
    checkTopologyProperties(top, requiredTopologyProperties());
    // Get the number of atoms from the topology if it is not given.
    if (natoms <= 0)
    {
        natoms = top->natoms;
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
    gmx_ana_selcollection_t* sc = &impl_->sc_;
    // Do this first, as it allocates memory, while the others don't throw.
    gmx_ana_index_init_simple(&sc->gall, natoms);
    sc->top = top;
    sc->pcc.setTopology(top);
}


void SelectionCollection::setIndexGroups(gmx_ana_indexgrps_t* grps)
{
    GMX_RELEASE_ASSERT(grps == nullptr || !impl_->bExternalGroupsSet_,
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

SelectionTopologyProperties SelectionCollection::requiredTopologyProperties() const
{
    SelectionTopologyProperties props;

    // These should not throw, because has been checked earlier.
    props.merge(impl_->requiredTopologyPropertiesForPositionType(impl_->rpost_, false));
    const bool forcesRequested = impl_->areForcesRequested();
    props.merge(impl_->requiredTopologyPropertiesForPositionType(impl_->spost_, forcesRequested));

    SelectionTreeElementPointer sel = impl_->sc_.root;
    while (sel && !props.hasAll())
    {
        props.merge(sel->requiredTopologyProperties());
        sel = sel->next;
    }
    return props;
}


bool SelectionCollection::requiresIndexGroups() const
{
    SelectionTreeElementPointer sel = impl_->sc_.root;
    while (sel)
    {
        if (sel->requiresIndexGroups())
        {
            return true;
        }
        sel = sel->next;
    }
    return false;
}


SelectionList SelectionCollection::parseFromStdin(int count, bool bInteractive, const std::string& context)
{
    StandardInputStream inputStream;
    return parseInteractive(
            count, &inputStream, bInteractive ? &TextOutputFile::standardError() : nullptr, context);
}

namespace
{

//! Helper function to initialize status writer for interactive selection parsing.
std::unique_ptr<TextWriter> initStatusWriter(TextOutputStream* statusStream)
{
    std::unique_ptr<TextWriter> statusWriter;
    if (statusStream != nullptr)
    {
        statusWriter = std::make_unique<TextWriter>(statusStream);
        statusWriter->wrapperSettings().setLineLength(78);
    }
    return statusWriter;
}

} // namespace

SelectionList SelectionCollection::parseInteractive(int                count,
                                                    TextInputStream*   inputStream,
                                                    TextOutputStream*  statusStream,
                                                    const std::string& context)
{
    yyscan_t scanner;

    const std::unique_ptr<TextWriter> statusWriter(initStatusWriter(statusStream));
    _gmx_sel_init_lexer(
            &scanner, &impl_->sc_, statusWriter.get(), count, impl_->bExternalGroupsSet_, impl_->grps_);
    return runParser(scanner, inputStream, true, count, context);
}


SelectionList SelectionCollection::parseFromFile(const std::filesystem::path& filename)
{

    try
    {
        yyscan_t      scanner;
        TextInputFile file(filename);
        // TODO: Exception-safe way of using the lexer.
        _gmx_sel_init_lexer(&scanner, &impl_->sc_, nullptr, -1, impl_->bExternalGroupsSet_, impl_->grps_);
        _gmx_sel_set_lex_input_file(scanner, file.handle());
        return runParser(scanner, nullptr, false, -1, std::string());
    }
    catch (GromacsException& ex)
    {
        ex.prependContext(formatString("Error in parsing selections from file '%s'",
                                       filename.string().c_str()));
        throw;
    }
}


SelectionList SelectionCollection::parseFromString(const std::string& str)
{
    yyscan_t scanner;

    _gmx_sel_init_lexer(&scanner, &impl_->sc_, nullptr, -1, impl_->bExternalGroupsSet_, impl_->grps_);
    _gmx_sel_set_lex_input_str(scanner, str.c_str());
    return runParser(scanner, nullptr, false, -1, std::string());
}


void SelectionCollection::compile()
{
    checkTopologyProperties(impl_->sc_.top, requiredTopologyProperties());
    if (!impl_->bExternalGroupsSet_)
    {
        setIndexGroups(nullptr);
    }
    if (impl_->debugLevel_ != Impl::DebugLevel::None)
    {
        printTree(stderr, false);
    }

    compileSelection(this);

    if (impl_->debugLevel_ != Impl::DebugLevel::None)
    {
        std::fprintf(stderr, "\n");
        printTree(stderr, false);
        std::fprintf(stderr, "\n");
        impl_->sc_.pcc.printTree(stderr);
        std::fprintf(stderr, "\n");
    }
    impl_->sc_.pcc.initEvaluation();
    if (impl_->debugLevel_ != Impl::DebugLevel::None)
    {
        impl_->sc_.pcc.printTree(stderr);
        std::fprintf(stderr, "\n");
    }

    // TODO: It would be nicer to associate the name of the selection option
    // (if available) to the error message.
    SelectionDataList::const_iterator iter;
    for (iter = impl_->sc_.sel.begin(); iter != impl_->sc_.sel.end(); ++iter)
    {
        const internal::SelectionData& sel = **iter;
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
            if (sel.hasFlag(efSelection_OnlySorted))
            {
                if (!sel.hasSortedAtomIndices())
                {
                    const std::string message = formatString(
                            "Selection '%s' does not evaluate to atoms in an "
                            "ascending (sorted) order. "
                            "This is not allowed in this context.",
                            sel.selectionText());
                    GMX_THROW(InvalidInputError(message));
                }
            }
        }
        if (sel.hasFlag(efSelection_DisallowEmpty))
        {
            if (sel.posCount() == 0)
            {
                std::string message =
                        formatString("Selection '%s' never matches any atoms.", sel.selectionText());
                GMX_THROW(InvalidInputError(message));
            }
        }
    }
    impl_->rpost_.clear();
    impl_->spost_.clear();
}


void SelectionCollection::evaluate(t_trxframe* fr, t_pbc* pbc)
{
    checkTopologyProperties(impl_->sc_.top, requiredTopologyProperties());
    if (fr->bIndex)
    {
        gmx_ana_index_t g;
        gmx_ana_index_set(&g, fr->natoms, fr->index, 0);
        GMX_RELEASE_ASSERT(gmx_ana_index_check_sorted(&g),
                           "Only trajectories with atoms in ascending order "
                           "are currently supported");
        if (!gmx_ana_index_contains(&g, &impl_->requiredAtoms_))
        {
            const std::string message = formatString(
                    "Trajectory does not contain all atoms required for "
                    "evaluating the provided selections.");
            GMX_THROW(InconsistentInputError(message));
        }
    }
    else
    {
        const int maxAtomIndex = gmx_ana_index_get_max_index(&impl_->requiredAtoms_);
        if (fr->natoms <= maxAtomIndex)
        {
            const std::string message = formatString(
                    "Trajectory has less atoms (%d) than what is required for "
                    "evaluating the provided selections (atoms up to index %d "
                    "are required).",
                    fr->natoms,
                    maxAtomIndex + 1);
            GMX_THROW(InconsistentInputError(message));
        }
    }
    impl_->sc_.pcc.initFrame(fr);

    SelectionEvaluator evaluator;
    evaluator.evaluate(this, fr, pbc);

    if (impl_->debugLevel_ == Impl::DebugLevel::Evaluated || impl_->debugLevel_ == Impl::DebugLevel::Full)
    {
        std::fprintf(stderr, "\n");
        printTree(stderr, true);
    }
}


void SelectionCollection::evaluateFinal(int nframes)
{
    SelectionEvaluator evaluator;
    evaluator.evaluateFinal(this, nframes);
}

std::optional<Selection> SelectionCollection::selection(std::string_view selName) const
{
    const auto& selections = impl_->sc_.sel;
    if (const auto foundIter = std::find_if(
                selections.cbegin(),
                selections.cend(),
                [selName](const auto& selection) { return selection->name() == selName; });
        foundIter != selections.end())
    {
        return Selection(foundIter->get());
    }
    return std::nullopt;
}


void SelectionCollection::printTree(FILE* fp, bool bValues) const
{
    SelectionTreeElementPointer sel = impl_->sc_.root;
    while (sel)
    {
        _gmx_selelem_print_tree(fp, *sel, bValues, 0);
        sel = sel->next;
    }
}


void SelectionCollection::printXvgrInfo(FILE* out) const
{
    const gmx_ana_selcollection_t& sc = impl_->sc_;
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

void swap(SelectionCollection& lhs, SelectionCollection& rhs) noexcept
{
    lhs.swap(rhs);
}

} // namespace gmx
