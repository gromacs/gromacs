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

#include <typedefs.h>

#include "../options/options.h"
#include "indexutil.h"
#include "selectioncollection.h"

namespace gmx
{
class Selection;
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
    /** Array of compiled selections. */
    std::vector<gmx::Selection *>  sel;
    /** Number of variables defined. */
    int                            nvars;
    /** Selection strings for variables. */
    char                         **varstrs;

    /** Topology for the collection. */
    t_topology                    *top;
    /** Index group that contains all the atoms. */
    struct gmx_ana_index_t         gall;
    /** Position calculation collection used for selection position evaluation. */
    struct gmx_ana_poscalc_coll_t *pcc;
    /** Memory pool used for selection evaluation. */
    struct gmx_sel_mempool_t      *mempool;
    /** Parser symbol table. */
    struct gmx_sel_symtab_t     *symtab;
};

namespace gmx
{

class SelectionOptionStorage;

/*! \internal \brief
 * Private implemention class for SelectionCollection.
 *
 * \ingroup module_selection
 */
class SelectionCollection::Impl
{
    public:
        struct SelectionRequest
        {
            SelectionRequest(const std::string &name, const std::string &descr,
                             int count, SelectionOptionStorage *storage)
                : name(name), descr(descr), count(count), storage(storage)
            { }

            std::string                 name;
            std::string                 descr;
            int                         count;
            SelectionOptionStorage     *storage;
        };

        //! Shorthand for a list of selections stored internally.
        typedef std::vector<Selection *> SelectionList;
        //! Shorthand for a list of selection requests.
        typedef std::vector<SelectionRequest> RequestList;

        //! Possibel flags for the selection collection.
        enum Flag
        {
            efOwnPositionCollection = 1<<0
        };

        //! Creates a new selection collection.
        explicit Impl(gmx_ana_poscalc_coll_t *pcc);
        ~Impl();

        //! Returns true if the given flag has been set.
        bool hasFlag(Flag flag) const { return _flags & flag; }
        //! Sets or clears the given flag.
        void setFlag(Flag flat, bool bSet);
        //! Clears the symbol table of the selection collection.
        void clearSymbolTable();
        //! Registers the default selection methods for the collection.
        int registerDefaultMethods();
        /*! \brief
         * Helper function that runs the parser once the tokenizer has been
         * initialized.
         *
         * \param[in,out] scanner Scanner data structure.
         * \param[in]     maxnr   Maximum number of selections to parse
         *      (if -1, parse as many as provided by the user).
         * \param[out]    output  Vector to which parsed selections are
         *      appended.
         * \retval        0 on success.
         * \retval        ::eeInvalidInput on error.
         *
         * Does not clear \p output.
         */
        int runParser(void *scanner, int maxnr,
                      std::vector<Selection *> *output);
        void requestSelections(const std::string &name,
                               const std::string &descr,
                               int count, SelectionOptionStorage *storage);

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
        //! Flags for various properties of the collection.
        unsigned long           _flags;
        //! External index groups (can be NULL).
        gmx_ana_indexgrps_t    *_grps;
        //! List of selections requested for later parsing.
        RequestList             _requests;
};

} // namespace gmx

/*! \addtogroup module_selection
 * \{
 */

/* In compiler.cpp */
/*! \internal \brief
 * Prepares the selections for evaluation and performs some optimizations.
 */
int
gmx_ana_selcollection_compile(gmx::SelectionCollection *coll);

/* In evaluate.cpp */
/*! \internal \brief
 * Evaluates the selection.
 */
int
gmx_ana_selcollection_evaluate(gmx_ana_selcollection_t *sc,
                               t_trxframe *fr, t_pbc *pbc);
/*! \internal \brief
 * Evaluates the largest possible index groups from dynamic selections.
 */
int
gmx_ana_selcollection_evaluate_fin(gmx_ana_selcollection_t *sc, int nframes);

/*!\}*/

#endif
