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
 * Implements gmx::SelectionOptionManager.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#include "selectionoptionmanager.h"

#include <cstdio>

#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/format.h"

#include "selectionoptionstorage.h"

namespace gmx
{

/********************************************************************
 * SelectionOptionManager::Impl
 */

/*! \internal \brief
 * Private implemention class for SelectionOptionManager.
 *
 * \ingroup module_selection
 */
class SelectionOptionManager::Impl
{
    public:
        /*! \brief
         * Request for postponed parsing of selections.
         */
        struct SelectionRequest
        {
            //! Initializes a request for the given option.
            explicit SelectionRequest(SelectionOptionStorage *storage)
                : storage_(storage)
            {
            }

            //! Returns name of the requested selection optin.
            const std::string &name() const
            {
                return storage_->name();
            }
            //! Returns description for the requested selection option.
            const std::string &description() const
            {
                return storage_->description();
            }
            /*! \brief
             * Returns the number of selections requested in this request.
             *
             * -1 indicates no upper limit.
             */
            int count() const
            {
                return storage_->maxValueCount();
            }

            //! Storage object to which the selections will be added.
            SelectionOptionStorage     *storage_;
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
        explicit Impl(SelectionCollection *collection);

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
         */
        void placeSelectionsInRequests(const SelectionList &selections);

        //! Selection collection to which selections are stored.
        SelectionCollection    &collection_;
        //! List of selections requested for later parsing.
        RequestList             requests_;
};

SelectionOptionManager::Impl::Impl(SelectionCollection *collection)
    : collection_(*collection)
{
}

void SelectionOptionManager::Impl::placeSelectionsInRequests(
        const SelectionList &selections)
{
    RequestsClearer clearRequestsOnExit(&requests_);

    SelectionList::const_iterator first = selections.begin();
    SelectionList::const_iterator last = first;
    RequestList::const_iterator i;
    // TODO: Improve error messages.
    for (i = requests_.begin(); i != requests_.end(); ++i)
    {
        const SelectionRequest &request = *i;
        if (request.count() > 0)
        {
            if (selections.end() - first < request.count())
            {
                GMX_THROW(InvalidInputError("Too few selections provided"));
            }
            last = first + request.count();
        }
        else
        {
            if (i != requests_.end() - 1)
            {
                GMX_THROW(InvalidInputError(
                            formatString("Request for selection '%s' must "
                                         "not be followed by others",
                                         request.name().c_str())));
            }
            last = selections.end();
        }
        SelectionList curr(first, last);
        request.storage_->addSelections(curr, true);
        first = last;
    }
    if (last != selections.end())
    {
        GMX_THROW(InvalidInputError("Too many selections provided"));
    }
}


/********************************************************************
 * SelectionOptionManager
 */

SelectionOptionManager::SelectionOptionManager(SelectionCollection *collection)
    : impl_(new Impl(collection))
{
}

SelectionOptionManager::~SelectionOptionManager()
{
}

SelectionCollection &
SelectionOptionManager::selectionCollection()
{
    return impl_->collection_;
}

void
SelectionOptionManager::requestDelayedParsing(SelectionOptionStorage *storage)
{
    impl_->requests_.push_back(Impl::SelectionRequest(storage));
}

void
SelectionOptionManager::parseRequestedFromStdin(bool bInteractive)
{
    Impl::RequestsClearer clearRequestsOnExit(&impl_->requests_);

    Impl::RequestList::const_iterator i;
    for (i = impl_->requests_.begin(); i != impl_->requests_.end(); ++i)
    {
        const Impl::SelectionRequest &request = *i;
        if (bInteractive)
        {
            std::fprintf(stderr, "\nSpecify ");
            if (request.count() < 0)
            {
                std::fprintf(stderr, "any number of selections");
            }
            else if (request.count() == 1)
            {
                std::fprintf(stderr, "a selection");
            }
            else
            {
                std::fprintf(stderr, "%d selections", request.count());
            }
            std::fprintf(stderr, " for option '%s' (%s):\n",
                         request.name().c_str(), request.description().c_str());
            std::fprintf(stderr, "(one selection per line, 'help' for help%s)\n",
                         request.count() < 0 ? ", Ctrl-D to end" : "");
        }
        SelectionList selections;
        impl_->collection_.parseFromStdin(request.count(), bInteractive, &selections);
        request.storage_->addSelections(selections, true);
    }
}

void
SelectionOptionManager::parseRequestedFromFile(const std::string &filename)
{
    SelectionList selections;
    impl_->collection_.parseFromFile(filename, &selections);
    impl_->placeSelectionsInRequests(selections);
}

void
SelectionOptionManager::parseRequestedFromString(const std::string &str)
{
    SelectionList selections;
    impl_->collection_.parseFromString(str, &selections);
    impl_->placeSelectionsInRequests(selections);
}

} // namespace gmx
