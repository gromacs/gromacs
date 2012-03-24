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
/*! \libinternal \file
 * \brief
 * Declares gmx::SelectionCollection.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
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
struct gmx_ana_poscalc_coll_t;

namespace gmx
{

class Options;
class SelectionCompiler;
class SelectionEvaluator;
class SelectionOptionStorage;

/*! \libinternal \brief
 * Collection of selections.
 *
 * Some default values must then be set with
 * gmx_ana_selcollection_set_refpostype() and
 * gmx_ana_selcollection_set_outpostype().
 *
 * After setting the default values, one or more selections can be parsed
 * with one or more calls to parseFromStdin(), parseFromFile(), and/or
 * parseFromString().  After all selections are parsed, the topology must be
 * set with setTopology() unless requiresTopology() returns false (the topology
 * can also be set earlier).  Once all selections are parsed, they must be
 * compiled all at once using compile().
 * After compilation, dynamic selections have the maximum number of atoms they
 * can evaluate to, but positions have undefined values.  evaluate() can be
 * used to update the selections for a new frame.
 * evaluateFinal() can be called after all the frames have been processed to
 * restore the selection values back to the ones they were after compile().
 *
 * At any point, requiresTopology() can be called to see whether the
 * information provided so far requires loading the topology.
 * printTree() can be used to print the internal representation of the
 * selections (mostly useful for debugging).
 *
 * \inlibraryapi
 * \ingroup module_selection
 */
class SelectionCollection
{
    public:
        /*! \brief
         * Creates an empty selection collection.
         *
         * \param[in] pcc  Position calculation collection to use for selection
         *      evaluation.
         *
         * If \p pcc is NULL, an internal collection is created and managed by
         * the object.
         */
        explicit SelectionCollection(gmx_ana_poscalc_coll_t *pcc);
        ~SelectionCollection();

        /*! \brief
         * Initializes options for setting global properties on the collection.
         */
        Options &initOptions();

        /*! \brief
         * Sets the default reference position handling for a selection
         * collection.
         *
         * \param[in]     type      Default selection reference position type
         *   (one of the strings acceptable for gmx_ana_poscalc_type_from_enum()).
         *
         * Should be called before calling the parser functions, unless
         * initOptions() has been called.  In the latter case, can still be
         * used to override the default value and/or the value provided through
         * the Options object.
         */
        void setReferencePosType(const char *type);
        /*! \brief
         * Sets the default reference position handling for a selection
         * collection.
         *
         * \param[in]     type      Default selection output position type
         *   (one of the strings acceptable for gmx_ana_poscalc_type_from_enum()).
         *
         * Should be called before calling the parser functions, unless
         * initOptions() has been called.  In the latter case, can still be
         * used to override the default value and/or the value provided through
         * the Options object.
         */
        void setOutputPosType(const char *type);
        /*! \brief
         * Sets the debugging level for the selection collection.
         */
        void setDebugLevel(int debuglevel);

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
         */
        bool requiresTopology() const;
        /*! \brief
         * Sets the topology for the collection.
         *
         * \param[in]     top       Topology data.
         * \param[in]     natoms    Number of atoms. If <=0, the number of
         *      atoms in the topology is used.
         * \retval  0 on success.
         * \retval  ::eeInvalidValue if \p top is NULL and \p natoms <= 0.
         *
         * The topology is also set for the position calculation collection
         * associated with the collection.
         *
         * \p natoms determines the largest atom index that can be selected by
         * the selection: even if the topology contains more atoms, they will
         * not be selected.
         */
        void setTopology(t_topology *top, int natoms);
        /*! \brief
         * Sets the external index groups to use for the selections.
         *
         * Can be called only once with non-NULL \p grps.
         */
        void setIndexGroups(gmx_ana_indexgrps_t *grps);
        /*! \brief
         * Parses selection(s) from standard input for options not yet
         * provided.
         *
         * \param[in]  bInteractive Whether the parser should behave
         *      interactively.
         *
         * This method cooperates with SelectionOption to allow interactive
         * input of missing selections after all options have been processed.
         * It should be called after the Options::finish() method has been
         * called on all options that add selections to this collection.
         */
        void parseRequestedFromStdin(bool bInteractive);
        /*! \brief
         * Parses selection(s) from a string for options not yet provided.
         *
         * \param[in]  str     String to parse.
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
         * \retval     0 on success.
         * \retval     ::eeInvalidInput on syntax error (an interactive parser
         *      only returns this if an incorrect number of selections is
         *      provided).
         *
         * Parsed selections are appended to \p output without clearing it
         * first.  If parsing fails, \p output is not modified.
         *
         * The objects returned in \p output remain valid for the lifetime of
         * the selection collection, and should not be freed by the user.
         * Some information about the selections only becomes available once
         * compile() has been called.
         */
        void parseFromStdin(int count, bool bInteractive,
                            SelectionList *output);
        /*! \brief
         * Parses selection(s) from a file.
         *
         * \param[in]  filename Name of the file to parse selections from.
         * \param[out] output   Vector to which parsed selections are appended.
         * \retval     0 on success.
         * \retval     ::eeInvalidInput on syntax error.
         *
         * Parsed selections are appended to \p output without clearing it
         * first.  If parsing fails, \p output is not modified.
         *
         * The objects returned in \p output remain valid for the lifetime of
         * the selection collection, and should not be freed by the user.
         * Some information about the selections only becomes available once
         * compile() has been called.
         */
        void parseFromFile(const std::string &filename,
                           SelectionList *output);
        /*! \brief
         * Parses selection(s) from a string.
         *
         * \param[in]  str      String to parse selections from.
         * \param[out] output   Vector to which parsed selections are appended.
         * \retval     0 on success.
         * \retval     ::eeInvalidInput on syntax error.
         *
         * Parsed selections are appended to \p output without clearing it
         * first.  If parsing fails, \p output is not modified.
         *
         * The objects returned in \p output remain valid for the lifetime of
         * the selection collection, and should not be freed by the user.
         * Some information about the selections only becomes available once
         * compile() has been called.
         */
        void parseFromString(const std::string &str,
                             SelectionList *output);
        /*! \brief
         * Prepares the selections for evaluation and performs optimizations.
         *
         * \retval  0 on successful compilation, a non-zero error code on error.
         *
         * Before compilation, selections should have been added to the
         * collection using the parseFrom*() functions.
         * The compiled selection collection can be passed to evaluate() to
         * evaluate the selection for a frame.
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
         * \returns   0 on successful evaluation, a non-zero error code on error.
         */
        void evaluate(t_trxframe *fr, t_pbc *pbc);
        /*! \brief
         * Evaluates the largest possible index groups from dynamic selections.
         *
         * \param[in] nframes Total number of frames.
         *
         * This function does not throw.
         */
        void evaluateFinal(int nframes);

        /*! \brief
         * Prints a human-readable version of the internal selection element
         * tree.
         *
         * \param[in] fp      File handle to receive the output.
         * \param[in] bValues If true, the evaluated values of selection
         *      elements are printed as well.
         */
        void printTree(FILE *fp, bool bValues) const;
        /*! \brief
         * Prints the selection strings into an XVGR file as comments.
         *
         * \param[in] fp   Output file.
         * \param[in] oenv Output options structure.
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
