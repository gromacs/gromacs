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
/*! \file
 * \brief
 * Declares gmx::SelectionCollection.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONCOLLECTION_H
#define GMX_SELECTION_SELECTIONCOLLECTION_H

#include <cstdio>

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/selection/selection.h" // For gmx::SelectionList

struct gmx_ana_indexgrps_t;
struct gmx_mtop_t;
struct t_pbc;
struct t_trxframe;

namespace gmx
{

class IOptionsContainer;
class SelectionEvaluator;
class TextInputStream;
class TextOutputStream;
struct SelectionTopologyProperties;

/*! \brief
 * Collection of selections.
 *
 * This class is the main interface to the core of the selection engine.
 * It is used to initialize and manage a collection of selections that share
 * the same topology.  Selections within one collection can share variables and
 * can be optimized together.  Selections from two different collections do not
 * interact.
 *
 * The constructor creates an empty selection collection object.  To initialize
 * the object, either call initOptions(), or both setReferencePosType() and
 * setOutputPosType().  See these methods for more details on the
 * initialization options.
 *
 * SelectionCollections can be copied. Copies retain the same pointers to external indices (if set)
 * and the topology (if set), and are compiled if the copied collection is compiled. Selection
 * objects created from a given SelectionCollection are tied only to the original collection, so
 * a copy of a SelectionCollection will not update pre-existing Selections on evaluate() calls.
 *
 * After setting the default values, one or more selections can be parsed with
 * one or more calls to parseInteractive(), parseFromStdin(), parseFromFile(), and/or
 * parseFromString().  After all selections are parsed, the topology must be
 * set with setTopology() unless requiresTopology() returns false (the topology
 * can also be set earlier).
 * setIndexGroups() must also be called if external index group references are
 * used in the selections; it can be called at any point before compile().
 * Once all selections are parsed, they must be compiled all at once using
 * compile().
 *
 * After compilation, dynamic selections have the maximum number of atoms they
 * can evaluate to, but positions have undefined values (see \ref Selection and
 * SelectionPosition).  evaluate() can be used to update the selections for a
 * new frame.  evaluateFinal() can be called after all the frames have been
 * processed to restore the selection values back to the ones they were after
 * compile().
 *
 * At any point, requiresTopology() can be called to see whether the
 * information provided so far requires loading the topology.
 * Similarly, requiresIndexGroups() tells whether external index groups are
 * requires.
 * printTree() can be used to print the internal representation of the
 * selections (mostly useful for debugging).
 *
 * Note that for trajectory analysis using TrajectoryAnalysisModule, the
 * SelectionCollection object is managed by Gromacs, and \ref Selection objects
 * are obtained from SelectionOption.
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class SelectionCollection
{
public:
    //! Flag for initOptions() to select how to behave with -seltype option.
    enum SelectionTypeOption
    {
        /*! \brief
         * Add the option for the user to select default value for
         * setOutputPosType().
         */
        IncludeSelectionTypeOption,
        /*! \brief
         * Do not add the option, selections will always select atoms by
         * default.
         */
        AlwaysAtomSelections
    };

    /*! \brief
     * Creates an empty selection collection.
     *
     * \throws  std::bad_alloc if out of memory.
     */
    SelectionCollection();
    ~SelectionCollection();

    SelectionCollection(const SelectionCollection& rhs);
    SelectionCollection& operator=(SelectionCollection rhs);
    void                 swap(SelectionCollection& rhs) noexcept;

    /*! \brief
     * Initializes options for setting global properties on the collection.
     *
     * \param[in,out] options Options object to initialize.
     * \param[in]     selectionTypeOption
     *     Whether to add option to influence setOutputPosType().
     * \throws        std::bad_alloc if out of memory.
     *
     * Adds options to \p options that can be used to set the default
     * position types (see setReferencePosType() and setOutputPosType())
     * and debugging flags.
     */
    void initOptions(IOptionsContainer* options, SelectionTypeOption selectionTypeOption);

    /*! \brief
     * Sets the default reference position handling for a selection
     * collection.
     *
     * \param[in]     type      Default selection reference position type
     *     (one of the strings acceptable for
     *     PositionCalculationCollection::typeFromEnum()).
     * \throws  InternalError if \p type is invalid.
     *
     * Should be called before calling the parser functions, unless
     * initOptions() has been called.  In the latter case, can still be
     * used to override the default value (before initOptions() is called)
     * and/or the value provided through the Options object.
     *
     * Strong exception safety.
     */
    void setReferencePosType(const char* type);
    /*! \brief
     * Sets the default reference position handling for a selection
     * collection.
     *
     * \param[in]     type      Default selection output position type
     *     (one of the strings acceptable for
     *     PositionCalculationCollection::typeFromEnum()).
     * \throws  InternalError if \p type is invalid.
     *
     * Should be called before calling the parser functions, unless
     * initOptions() has been called.  In the latter case, can still be
     * used to override the default value (before initOptions() is called)
     * and/or the value provided through the Options object.
     *
     * Strong exception safety.
     */
    void setOutputPosType(const char* type);
    /*! \brief
     * Sets the debugging level for the selection collection.
     *
     * \param[in]   debugLevel  Debug level to set (0 = no debug
     *      information).
     *
     * initOptions() creates debugging options that can also be used to set
     * the debug level.  These are normally hidden, but if this method is
     * called before initOptions() with a non-zero \p debugLevel, they are
     * made visible.
     *
     * Mostly useful for debugging tools.
     *
     * Does not throw.
     */
    void setDebugLevel(int debugLevel);

    /*! \brief
     * Returns what topology information is required for evaluation.
     *
     * \returns What topology information is required for compiling and/or
     *     evaluating the selections in the collection.
     *
     * Before the parser functions have been called, the return value is
     * based just on the position types set.
     * After parser functions have been called, the return value also takes
     * into account the selection keywords used.
     * After the compiler has been called, the return value is final and
     * also considers possible force evaluation requested for the
     * selections.
     *
     * Does not throw.
     */
    SelectionTopologyProperties requiredTopologyProperties() const;
    /*! \brief
     * Returns true if the collection requires external index groups.
     *
     * \returns true if any selection has an unresolved index group reference.
     *
     * The return value is `false` after setIndexGroups() has been called.
     *
     * Does not throw.
     */
    bool requiresIndexGroups() const;
    /*! \brief
     * Sets the topology for the collection.
     *
     * \param[in]     top       Topology data.
     * \param[in]     natoms    Number of atoms. If <=0, the number of
     *      atoms in the topology is used.
     *
     * Either the topology must be provided, or \p natoms must be > 0.
     *
     * \p natoms determines the largest atom index that can be selected by
     * the selection: even if the topology contains more atoms, they will
     * not be selected.
     *
     * Does not throw currently, but this is subject to change when more
     * underlying code is converted to C++.
     */
    void setTopology(const gmx_mtop_t* top, int natoms);
    /*! \brief
     * Sets the external index groups to use for the selections.
     *
     * \param[in]  grps  Index groups to use for the selections.
     * \throws  std::bad_alloc if out of memory.
     * \throws  InconsistentInputError if a group reference cannot be resolved.
     *
     * Only the first call to this method can have a non-NULL \p grps.
     * At this point, any selections that have already been provided are
     * searched for references to external groups, and the references are
     * replaced by the contents of the groups.  If any referenced group
     * cannot be found in \p grps (or if \p grps is NULL and there are any
     * references), InconsistentInputError is thrown.
     *
     * The selection collection keeps a reference to \p grps until this
     * method is called with a NULL \p grps.
     * If this method is not called before compile(), it is automatically
     * called as setIndexGroups(NULL).
     */
    void setIndexGroups(gmx_ana_indexgrps_t* grps);
    /*! \brief
     * Parses selection(s) from standard input.
     *
     * \param[in]  count    Number of selections to parse
     *      (if -1, parse as many as provided by the user).
     * \param[in]  bInteractive Whether the parser should behave
     *      interactively.
     * \param[in]  context  Context to print for interactive input.
     * \returns    Vector of parsed selections.
     * \throws     std::bad_alloc if out of memory.
     * \throws     InvalidInputError if there is a parsing error
     *      (an interactive parser only throws this if too few selections
     *      are provided and the user forced the end of input).
     *
     * The returned objects remain valid for the lifetime of
     * the selection collection.
     * Some information about the selections only becomes available once
     * compile() has been called; see \ref Selection.
     *
     * The string provided to \p context should be such that it can replace
     * the three dots in "Specify selections ...:".  It can be empty.
     */
    SelectionList parseFromStdin(int count, bool bInteractive, const std::string& context);
    /*! \brief
     * Parses selection(s) interactively using provided streams.
     *
     * \param[in]  count    Number of selections to parse
     *      (if -1, parse as many as provided by the user).
     * \param[in]  inputStream  Stream to use for input.
     * \param[in]  outputStream Stream to use for output
     *      (if NULL, the parser runs non-interactively and does not
     *      produce any status messages).
     * \param[in]  context  Context to print for interactive input.
     * \returns    Vector of parsed selections.
     * \throws     std::bad_alloc if out of memory.
     * \throws     InvalidInputError if there is a parsing error
     *      (an interactive parser only throws this if too few selections
     *      are provided and the user forced the end of input).
     *
     * Works the same as parseFromStdin(), except that the caller can
     * provide streams to use instead of `stdin` and `stderr`.
     *
     * Mainly usable for unit testing interactive input.
     */
    SelectionList parseInteractive(int                count,
                                   TextInputStream*   inputStream,
                                   TextOutputStream*  outputStream,
                                   const std::string& context);
    /*! \brief
     * Parses selection(s) from a file.
     *
     * \param[in]  filename Name of the file to parse selections from.
     * \returns    Vector of parsed selections.
     * \throws     std::bad_alloc if out of memory.
     * \throws     InvalidInputError if there is a parsing error.
     *
     * The returned objects remain valid for the lifetime of
     * the selection collection.
     * Some information about the selections only becomes available once
     * compile() has been called; see \ref Selection.
     */
    SelectionList parseFromFile(const std::filesystem::path& filename);
    /*! \brief
     * Parses selection(s) from a string.
     *
     * \param[in]  str      String to parse selections from.
     * \returns    Vector of parsed selections.
     * \throws     std::bad_alloc if out of memory.
     * \throws     InvalidInputError if there is a parsing error.
     *
     * The returned objects remain valid for the lifetime of
     * the selection collection.
     * Some information about the selections only becomes available once
     * compile() has been called; see \ref Selection.
     */
    SelectionList parseFromString(const std::string& str);
    /*! \brief
     * Prepares the selections for evaluation and performs optimizations.
     *
     * \throws  InconsistentInputError if topology is required but not set.
     * \throws  InvalidInputError if setIndexGroups() has not been called
     *      and there are index group references.
     * \throws  unspecified if compilation fails (TODO: list/reduce these).
     *
     * Before compilation, selections should have been added to the
     * collection using the parseFrom*() functions.
     * The compiled selection collection can be passed to evaluate() to
     * evaluate the selection for a frame.
     * Before the compiled selection is evaluated, the selections indicate
     * the maximal set of atoms/positions to which they can be evaluated;
     * see \ref Selection.
     *
     * If an error occurs, the collection is cleared.
     *
     * The covered fraction information is initialized to ::CFRAC_NONE for
     * all selections.
     */
    void compile();
    /*! \brief
     * Evaluates selections in the collection.
     *
     * \param[in] fr  Frame for which the evaluation should be carried out.
     * \param[in] pbc PBC data, or NULL if no PBC should be used.
     * \throws    unspeficied  Multiple possible exceptions to indicate
     *      evaluation failure (TODO: enumerate).
     */
    void evaluate(t_trxframe* fr, t_pbc* pbc);
    /*! \brief
     * Evaluates the largest possible index groups from dynamic selections.
     *
     * \param[in] nframes Total number of frames.
     *
     * This method restores the selections to the state they were after
     * compile().
     *
     * \p nframes should equal the number of times evaluate() has been
     * called.
     *
     * Does not throw.
     */
    void evaluateFinal(int nframes);
    /*! \brief
     * Retrieves a Selection handle for the selection with the given name
     *
     * @param selName name of the selection to return
     * @return The selection with the given name, or nullopt if no such selection exists.
     */
    [[nodiscard]] std::optional<Selection> selection(std::string_view selName) const;
    /*! \brief
     * Prints a human-readable version of the internal selection element
     * tree.
     *
     * \param[in] fp      File handle to receive the output.
     * \param[in] bValues If true, the evaluated values of selection
     *      elements are printed as well.
     *
     * The output is very techical, and intended for debugging purposes.
     *
     * Does not throw.
     */
    void printTree(FILE* fp, bool bValues) const;
    /*! \brief
     * Prints the selection strings into an XVGR file as comments.
     *
     * \param[in] fp   Output file.
     *
     * Does not throw.
     */
    void printXvgrInfo(FILE* fp) const;

private:
    class Impl;

    std::unique_ptr<Impl> impl_;

    // Needed for the compiler to freely modify the collection.
    friend void compileSelection(SelectionCollection* coll);
    // Needed for the evaluator to freely modify the collection.
    friend class SelectionEvaluator;
};

void swap(SelectionCollection& lhs, SelectionCollection& rhs) noexcept;

} // namespace gmx

#endif
