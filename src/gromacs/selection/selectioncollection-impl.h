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
 * Declares private implementation class for gmx::SelectionCollection.
 *
 * This header also defines ::gmx_ana_selcollection_t, which is used in the old
 * C code for handling selection collections.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONCOLLECTION_IMPL_H
#define GMX_SELECTION_SELECTIONCOLLECTION_IMPL_H

#include <string>
#include <vector>

#include "../legacyheaders/typedefs.h"

#include "../options/options.h"
#include "../utility/uniqueptr.h"
#include "indexutil.h"
#include "poscalc.h"
#include "selection.h" // For gmx::SelectionList
#include "selectioncollection.h"

namespace gmx
{
//! Smart pointer for managing an internal selection data object.
typedef gmx_unique_ptr<internal::SelectionData>::type SelectionDataPointer;
//! Container for storing a list of selections internally.
typedef std::vector<SelectionDataPointer> SelectionDataList;
}

/*! \internal \brief
 * Information for a collection of selections.
 *
 * \ingroup module_selection
 */
struct gmx_ana_selcollection_t
{
    /** Root of the selection element tree. */
    struct t_selelem           *root;
    /*! \brief
     * Array of compiled selections.
     *
     * Has the responsibility of managing the memory for the contained objects,
     * but note that gmx::Selection instances also hold pointers to the
     * objects.
     */
    gmx::SelectionDataList         sel;
    /** Number of variables defined. */
    int                            nvars;
    /** Selection strings for variables. */
    char                         **varstrs;

    /** Topology for the collection. */
    t_topology                    *top;
    /** Index group that contains all the atoms. */
    struct gmx_ana_index_t         gall;
    /** Position calculation collection used for selection position evaluation. */
    gmx::PositionCalculationCollection  pcc;
    /** Memory pool used for selection evaluation. */
    struct gmx_sel_mempool_t      *mempool;
    /** Parser symbol table. */
    struct gmx_sel_symtab_t     *symtab;
};

namespace gmx
{

class MessageStringCollector;
class SelectionOptionStorage;

/*! \internal \brief
 * Private implemention class for SelectionCollection.
 *
 * \ingroup module_selection
 */
class SelectionCollection::Impl
{
    public:
        /*! \brief
         * Request for postponed parsing of selections.
         *
         * Used to communicate what needs to be parsed with
         * parseRequestedFromStdin() or parseRequstedFromString().
         */
        struct SelectionRequest
        {
            //! Initializes a request for the given option.
            SelectionRequest(const std::string &name, const std::string &descr,
                             SelectionOptionStorage *storage)
                : name(name), descr(descr), storage(storage)
            {
            }

            /*! \brief
             * Returns the number of selections requested in this request.
             *
             * -1 indicates no upper limit.
             */
            int count() const;

            //! Name of the option to which this request relates to.
            std::string                 name;
            //! Description of the option to which this request relates to.
            std::string                 descr;
            //! Storage object to which the selections will be added.
            SelectionOptionStorage     *storage;
        };

        //! Collection for a list of selection requests.
        typedef std::vector<SelectionRequest> RequestList;

        /*! \brief
         * Helper class that clears a request list on scope exit.
         *
         * Methods in this class do not throw.
         */
        class RequestsClearer
        {
            public:
                //! Constructs an object that clears given list on scope exit.
                explicit RequestsClearer(RequestList *requests)
                    : requests_(requests)
                {
                }
                //! Clears the request list given to the constructor.
                ~RequestsClearer()
                {
                    requests_->clear();
                }

            private:
                RequestList    *requests_;
        };

        /*! \brief
         * Creates a new selection collection.
         *
         * \throws  std::bad_alloc if out of memory.
         */
        Impl();
        ~Impl();

        /*! \brief
         * Clears the symbol table of the selection collection.
         *
         * Does not throw.
         */
        void clearSymbolTable();
        /*! \brief
         * Helper function that runs the parser once the tokenizer has been
         * initialized.
         *
         * \param[in,out] scanner Scanner data structure.
         * \param[in]     maxnr   Maximum number of selections to parse
         *      (if -1, parse as many as provided by the user).
         * \param[out]    output  Vector to which parsed selections are
         *      appended.
         * \throws        std::bad_alloc if out of memory.
         * \throws        InvalidInputError if there is a parsing error.
         *
         * Parsed selections are appended to \p output without clearing it
         * first.  If parsing fails, \p output is not modified.
         *
         * Used internally to implement parseFromStdin(), parseFromFile() and
         * parseFromString().
         */
        void runParser(void *scanner, int maxnr,
                       SelectionList *output);
        /*! \brief
         * Adds a selection request for delayed user input.
         *
         * \param[in] name    Name for the requested selections.
         * \param[in] descr   Description of the requested selections.
         * \param     storage Storage object to receive the selections.
         * \throws    std::bad_alloc if out of memory.
         *
         * Strong exception safety.
         *
         * \see parseRequestedFromStdin()
         */
        void requestSelections(const std::string &name,
                               const std::string &descr,
                               SelectionOptionStorage *storage);
        /*! \brief
         * Assign selections from a list to pending requests.
         *
         * \param[in] selections  List of selections to assign.
         * \throws    std::bad_alloc if out of memory.
         * \throws    InvalidInputError if the assignment cannot be done
         *      (see parseRequestedFromFile() for documented conditions).
         *
         * Loops through \p selections and the pending requests lists in order,
         * and for each requests, assigns the first yet unassigned selections
         * from the list.
         *
         * Used to implement parseRequestedFromFile() and
         * parseRequestedFromStdin().
         */
        void placeSelectionsInRequests(const SelectionList &selections);
        /*! \brief
         * Replace group references by group contents.
         *
         * \param[in]    root    Root of selection tree to process.
         * \param        errors  Object for reporting any error messages.
         *
         * Recursively searches the selection tree for unresolved external
         * references.  If one is found, finds the referenced group in
         * \a _grps and replaces the reference with a constant element that
         * contains the atoms from the referenced group.  Any failures to
         * resolve references are reported to \p errors.
         *
         * Does not throw currently, but this is subject to change when more
         * underlying code is converted to C++.
         */
        void resolveExternalGroups(struct t_selelem *root,
                                   MessageStringCollector *errors);

        //! Internal data, used for interfacing with old C code.
        gmx_ana_selcollection_t _sc;
        //! Options object for setting global properties on the collection.
        Options                 _options;
        //! Default reference position type for selections.
        std::string             _rpost;
        //! Default output position type for selections.
        std::string             _spost;
        /*! \brief
         * Debugging level for the collection.
         *
         * Possible values:
         *  - 0: no debugging
         *  - 1: print selection trees after parsing and compilation
         *  - 2: like 1, also print intermediate compilation trees
         *  - 3: like 1, also print the tree after evaluation
         *  - 4: combine 2 and 3
         */
        int                     _debugLevel;
        //! Whether setIndexGroups() has been called.
        bool                    _bExternalGroupsSet;
        //! External index groups (can be NULL).
        gmx_ana_indexgrps_t    *_grps;
        //! List of selections requested for later parsing.
        RequestList             _requests;
};

/*! \internal \brief
 * Implements selection evaluation.
 *
 * This class is used to implement SelectionCollection::evaluate() and
 * SelectionCollection::evaluateFinal().
 *
 * \ingroup module_selection
 */
class SelectionEvaluator
{
    public:
        SelectionEvaluator();

        /*! \brief
         * Evaluates selections in a collection.
         */
        void evaluate(SelectionCollection *sc, t_trxframe *fr, t_pbc *pbc);
        /*! \brief
         * Evaluates the final state for dynamic selections.
         */
        void evaluateFinal(SelectionCollection *sc, int nframes);
};

} // namespace gmx

#endif
