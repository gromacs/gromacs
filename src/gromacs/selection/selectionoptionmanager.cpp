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

#include "gromacs/options/optionsvisitor.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectionfileoption.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

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
        //! Collection for list of option storage objects.
        typedef std::vector<SelectionOptionStorage *> OptionList;

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
        /*! \brief
         * Adds a request for each required option that is not yet set.
         *
         * \throws    std::bad_alloc if out of memory.
         */
        void requestUnsetRequiredOptions();

        //! Selection collection to which selections are stored.
        SelectionCollection    &collection_;
        //! List of selection options (storage objects) this manager manages.
        OptionList              options_;
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
    if (requests_.empty())
    {
        requestUnsetRequiredOptions();
    }

    RequestsClearer               clearRequestsOnExit(&requests_);

    SelectionList::const_iterator first = selections.begin();
    SelectionList::const_iterator last  = first;
    RequestList::const_iterator   i;
    for (i = requests_.begin(); i != requests_.end(); ++i)
    {
        const SelectionRequest &request = *i;
        if (request.count() > 0)
        {
            int remaining = selections.end() - first;
            if (remaining < request.count())
            {
                int assigned = first - selections.begin();
                GMX_THROW(InvalidInputError(formatString(
                                                    "Too few selections provided for '%s': "
                                                    "Expected %d selections, but only %d left "
                                                    "after assigning the first %d to other selections.",
                                                    request.name().c_str(), request.count(),
                                                    remaining, assigned)));
            }
            last = first + request.count();
        }
        else
        {
            RequestList::const_iterator nextRequest = i;
            ++nextRequest;
            if (nextRequest != requests_.end())
            {
                const char *name         = request.name().c_str();
                const char *conflictName = nextRequest->name().c_str();
                GMX_THROW(InvalidInputError(formatString(
                                                    "Ambiguous selections for '%s' and '%s': "
                                                    "Any number of selections is acceptable for "
                                                    "'%s', but you have requested subsequent "
                                                    "selections to be assigned to '%s'. "
                                                    "Resolution for such cases is not implemented, "
                                                    "and may be impossible.",
                                                    name, conflictName, name, conflictName)));
            }
            last = selections.end();
        }
        SelectionList curr(first, last);
        request.storage_->addSelections(curr, true);
        first = last;
    }
    if (last != selections.end())
    {
        int count     = selections.end() - selections.begin();
        int remaining = selections.end() - last;
        int assigned  = last - selections.begin();
        GMX_THROW(InvalidInputError(formatString(
                                            "Too many selections provided: "
                                            "Expected %d selections, but %d provided. "
                                            "Last %d selections could not be assigned to any option.",
                                            assigned, count, remaining)));
    }
}

void SelectionOptionManager::Impl::requestUnsetRequiredOptions()
{
    OptionList::const_iterator i;
    for (i = options_.begin(); i != options_.end(); ++i)
    {
        SelectionOptionStorage &storage = **i;
        if (storage.isRequired() && !storage.isSet())
        {
            requests_.push_back(SelectionRequest(&storage));
        }
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

void
SelectionOptionManager::registerOption(SelectionOptionStorage *storage)
{
    impl_->requests_.reserve(impl_->options_.size() + 1);
    impl_->options_.push_back(storage);
}

void
SelectionOptionManager::convertOptionValue(SelectionOptionStorage *storage,
                                           const std::string      &value)
{
    SelectionList selections = impl_->collection_.parseFromString(value);
    storage->addSelections(selections, false);
}

void
SelectionOptionManager::requestOptionDelayedParsing(
        SelectionOptionStorage *storage)
{
    impl_->requests_.push_back(Impl::SelectionRequest(storage));
}

void
SelectionOptionManager::parseRequestedFromStdin(bool bInteractive)
{
    Impl::RequestsClearer             clearRequestsOnExit(&impl_->requests_);

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
        SelectionList selections
            = impl_->collection_.parseFromStdin(request.count(), bInteractive);
        request.storage_->addSelections(selections, true);
    }
}

void
SelectionOptionManager::parseRequestedFromFile(const std::string &filename)
{
    SelectionList selections = impl_->collection_.parseFromFile(filename);
    try
    {
        impl_->placeSelectionsInRequests(selections);
    }
    catch (GromacsException &ex)
    {
        ex.prependContext(formatString(
                                  "Error in adding selections from file '%s'",
                                  filename.c_str()));
        throw;
    }
}

void
SelectionOptionManager::parseRequestedFromString(const std::string &str)
{
    SelectionList selections = impl_->collection_.parseFromString(str);
    impl_->placeSelectionsInRequests(selections);
}

/********************************************************************
 * Global functions
 */

namespace
{

/*! \internal \brief
 * Visitor that sets the manager for each selection option.
 *
 * \ingroup module_selection
 */
class SelectionOptionManagerSetter : public OptionsModifyingVisitor
{
    public:
        //! Construct a visitor that sets given manager.
        explicit SelectionOptionManagerSetter(SelectionOptionManager *manager)
            : manager_(manager)
        {
        }

        void visitSubSection(Options *section)
        {
            OptionsModifyingIterator iterator(section);
            iterator.acceptSubSections(this);
            iterator.acceptOptions(this);
        }

        void visitOption(OptionInfo *option)
        {
            SelectionOptionInfo *selOption
                = option->toType<SelectionOptionInfo>();
            if (selOption != NULL)
            {
                selOption->setManager(manager_);
            }
            SelectionFileOptionInfo *selFileOption
                = option->toType<SelectionFileOptionInfo>();
            if (selFileOption != NULL)
            {
                selFileOption->setManager(manager_);
            }
        }

    private:
        SelectionOptionManager *manager_;
};

}   // namespace

void setManagerForSelectionOptions(Options                *options,
                                   SelectionOptionManager *manager)
{
    SelectionOptionManagerSetter(manager).visitSubSection(options);
}

} // namespace gmx
