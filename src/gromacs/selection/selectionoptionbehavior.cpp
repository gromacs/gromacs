/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Implements gmx::SelectionOptionBehavior.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "selectionoptionbehavior.h"

#include <cstdio>

#include <string>

#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionfileoption.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/utility/filestream.h"

namespace gmx
{

/********************************************************************
 * ITopologyProvider
 */

ITopologyProvider::~ITopologyProvider()
{
}

/********************************************************************
 * SelectionOptionBehavior
 */

class SelectionOptionBehavior::Impl
{
    public:
        Impl(SelectionCollection *selections,
             ITopologyProvider   *topologyProvider)
            : selections_(*selections), topologyProvider_(*topologyProvider),
              manager_(selections), grps_(NULL)
        {
        }
        ~Impl()
        {
            if (grps_ != NULL)
            {
                gmx_ana_indexgrps_free(grps_);
            }
        }

        void promptSelections()
        {
            const bool isInteractive = StandardInputStream::instance().isInteractive();
            initIndexGroups();
            manager_.parseRequestedFromStdin(isInteractive);
            doneIndexGroups();
        }
        void initIndexGroups()
        {
            if (!selections_.requiresIndexGroups()
                && !manager_.hasRequestedSelections())
            {
                if (!ndxfile_.empty())
                {
                    std::fprintf(stderr, "NOTE: You provided an index file\n"
                                 "  %s\n(with -n), but it was not used by any selection.\n",
                                 ndxfile_.c_str());
                }
                selections_.setIndexGroups(NULL);
                return;
            }
            if (ndxfile_.empty())
            {
                t_topology *top = topologyProvider_.getTopology(false);
                gmx_ana_indexgrps_init(&grps_, top, NULL);
            }
            else
            {
                gmx_ana_indexgrps_init(&grps_, NULL, ndxfile_.c_str());
            }
            selections_.setIndexGroups(grps_);
        }
        void doneIndexGroups()
        {
            if (grps_ != NULL)
            {
                selections_.setIndexGroups(NULL);
                gmx_ana_indexgrps_free(grps_);
                grps_ = NULL;
            }
        }

        void compileSelections()
        {
            t_topology *top    = topologyProvider_.getTopology(selections_.requiresTopology());
            int         natoms = -1;
            if (top == NULL)
            {
                natoms = topologyProvider_.getAtomCount();
            }
            selections_.setTopology(top, natoms);
            selections_.compile();
        }

        SelectionCollection    &selections_;
        ITopologyProvider      &topologyProvider_;
        SelectionOptionManager  manager_;
        //! Name of the index file (empty if no index file provided).
        std::string             ndxfile_;
        gmx_ana_indexgrps_t    *grps_;
};

SelectionOptionBehavior::SelectionOptionBehavior(
        SelectionCollection *selections,
        ITopologyProvider   *topologyProvider)
    : impl_(new Impl(selections, topologyProvider))
{
}

SelectionOptionBehavior::~SelectionOptionBehavior()
{
}

void
SelectionOptionBehavior::initOptions(IOptionsContainer *options)
{
    options->addOption(FileNameOption("n")
                           .filetype(eftIndex).inputFile()
                           .store(&impl_->ndxfile_)
                           .defaultBasename("index")
                           .description("Extra index groups"));
    options->addOption(SelectionFileOption("sf"));
    impl_->manager_.initOptions(options);
}

void
SelectionOptionBehavior::initBehavior(Options *options)
{
    options->addManager(&impl_->manager_);
}

void
SelectionOptionBehavior::optionsFinished()
{
    impl_->promptSelections();
    impl_->compileSelections();
}

} // namespace gmx
