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
/*! \file
 * \brief
 * Declares gmx::SelectionCollection.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONCOLLECTION_H
#define GMX_SELECTION_SELECTIONCOLLECTION_H

#include <string>
#include <vector>

#include "../legacyheaders/typedefs.h"

#include "../utility/common.h"
#include "selection.h" // For gmx::SelectionList

struct gmx_ana_indexgrps_t;

namespace gmx
{

class Options;
class SelectionCompiler;
class SelectionEvaluator;
class SelectionOptionStorage;

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
 * After setting the default values, one or more selections can be parsed with
 * one or more calls to parseFromStdin(), parseFromFile(), and/or
 * parseFromString().  parseRequestedFromStdin() and parseRequestedFromString()
 * are provided for integration with SelectionOption.  After all selections are
 * parsed, the topology must be set with setTopology() unless
 * requiresTopology() returns false (the topology can also be set earlier).
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
        /*! \brief
         * Creates an empty selection collection.
         *
         * \throws  std::bad_alloc if out of memory.
         */
        SelectionCollection();
        ~SelectionCollection();

        /*! \brief
         * Initializes options for setting global properties on the collection.
         *
         * \returns Initialized options object.
         * \throws  std::bad_alloc if out of memory.
         *
         * The returned options can be used to set the default position types
         * (see setReferencePosType() and setOutputPosType()) and debugging
         * options.
         */
        Options &initOptions();

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
        void setReferencePosType(const char *type);
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
        void setOutputPosType(const char *type);
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
         * Returns true if the collection requires topology information for
         * evaluation.
         *
         * \returns true if any selection in the collection requires topology
         *      information.
         *
         * Before the parser functions have been called, the return value is
         * based just on the position types set.
         * After parser functions have been called, the return value also takes
         * into account the selection keywords used.
         *
         * Does not throw.
         */
        bool requiresTopology() const;
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
        void setTopology(t_topology *top, int natoms);
        /*! \brief
         * Sets the external index groups to use for the selections.
         *
         * \param[in]  grps  Index groups to use for the selections.
         * \throws  std::bad_alloc if out of memory.
         * \throws  InvalidInputError if a group reference cannot be resolved.
         *
         * Only the first call to this method can have a non-NULL \p grps.
         * At this point, any selections that have already been provided are
         * searched for references to external groups, and the references are
         * replaced by the contents of the groups.  If any referenced group
         * cannot be found in \p grps (or if \p grps is NULL and there are any
         * references), InvalidInputError is thrown.
         *
         * The selection collection keeps a reference to \p grps until this
         * method is called with a NULL \p grps.
         * If this method is not called before compile(), it is automatically
         * called as setIndexGroups(NULL).
         */
        void setIndexGroups(gmx_ana_indexgrps_t *grps);
        /*! \brief
         * Parses selection(s) from standard input for options not yet
         * provided.
         *
         * \param[in]  bInteractive Whether the parser should behave
         *      interactively.
         * \throws     unspecified  Can throw any exception thrown by
         *      parseFromStdin().
         * \throws     std::bad_alloc if out of memory.
         *
         * This method cooperates with SelectionOption to allow interactive
         * input of missing selections after all options have been processed.
         * It should be called after the Options::finish() method has been
         * called on all options that add selections to this collection.
         * For each required selection option that has not been given, as well
         * as for optional selection options that have been specified without
         * values, it will prompt the user to input the necessary selections.
         */
        void parseRequestedFromStdin(bool bInteractive);
        /*! \brief
         * Parses selection(s) from a string for options not yet provided.
         *
         * \param[in]  str     String to parse.
         * \throws     unspecified  Can throw any exception thrown by
         *      parseFromString().
         * \throws     std::bad_alloc if out of memory.
         * \throws     InvalidInputError if
         *      - the number of selections in \p str doesn't match the number
         *        requested.
         *      - any selection uses a feature that is not allowed for the
         *        corresponding option.
         * \throws     APIError if there is a request for any number of
         *      selections that is not the last (in which case it is not
         *      possible to determine which selections belong to which
         *      request).
         *
         * This method behaves as parseRequestedFromStdin(), but reads the
         * selections from a string instead of standard input.
         * This method is mainly used for testing.
         *
         * \see parseRequestedFromStdin()
         */
        void parseRequestedFromString(const std::string &str);
        /*! \brief
         * Parses selection(s) from standard input.
         *
         * \param[in]  count    Number of selections to parse
         *      (if -1, parse as many as provided by the user).
         * \param[in]  bInteractive Whether the parser should behave
         *      interactively.
         * \param[out] output   Vector to which parsed selections are appended.
         * \throws     std::bad_alloc if out of memory.
         * \throws     InvalidInputError if there is a parsing error
         *      (an interactive parser only throws this if too few selections
         *      are provided and the user forced the end of input).
         *
         * Parsed selections are appended to \p output without clearing it
         * first.  If parsing fails, \p output is not modified.
         *
         * The objects returned in \p output remain valid for the lifetime of
         * the selection collection.
         * Some information about the selections only becomes available once
         * compile() has been called; see \ref Selection.
         */
        void parseFromStdin(int count, bool bInteractive,
                            SelectionList *output);
        /*! \brief
         * Parses selection(s) from a file.
         *
         * \param[in]  filename Name of the file to parse selections from.
         * \param[out] output   Vector to which parsed selections are appended.
         * \throws     std::bad_alloc if out of memory.
         * \throws     InvalidInputError if there is a parsing error.
         *
         * Parsed selections are appended to \p output without clearing it
         * first.  If parsing fails, \p output is not modified.
         *
         * The objects returned in \p output remain valid for the lifetime of
         * the selection collection.
         * Some information about the selections only becomes available once
         * compile() has been called; see \ref Selection.
         */
        void parseFromFile(const std::string &filename,
                           SelectionList *output);
        /*! \brief
         * Parses selection(s) from a string.
         *
         * \param[in]  str      String to parse selections from.
         * \param[out] output   Vector to which parsed selections are appended.
         * \throws     std::bad_alloc if out of memory.
         * \throws     InvalidInputError if there is a parsing error.
         *
         * Parsed selections are appended to \p output without clearing it
         * first.  If parsing fails, \p output is not modified.
         *
         * The objects returned in \p output remain valid for the lifetime of
         * the selection collection.
         * Some information about the selections only becomes available once
         * compile() has been called; see \ref Selection.
         */
        void parseFromString(const std::string &str,
                             SelectionList *output);
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
        void evaluate(t_trxframe *fr, t_pbc *pbc);
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
        void printTree(FILE *fp, bool bValues) const;
        /*! \brief
         * Prints the selection strings into an XVGR file as comments.
         *
         * \param[in] fp   Output file.
         * \param[in] oenv Output options structure.
         *
         * Does not throw.
         */
        void printXvgrInfo(FILE *fp, output_env_t oenv) const;

    private:
        class Impl;

        PrivateImplPointer<Impl> _impl;

        /*! \brief
         * Needed for the compiler to freely modify the collection.
         */
        friend class SelectionCompiler;
        /*! \brief
         * Needed for the evaluator to freely modify the collection.
         */
        friend class SelectionEvaluator;
        /*! \brief
         * Needed for handling delayed selection parsing requests.
         */
        friend class SelectionOptionStorage;
};

} // namespace gmx

#endif
