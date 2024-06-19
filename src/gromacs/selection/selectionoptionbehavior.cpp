/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Implements gmx::SelectionOptionBehavior.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/selection/selectionoptionbehavior.h"

#include <cstdio>

#include <memory>
#include <string>
#include <vector>

#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionenums.h"
#include "gromacs/selection/selectionfileoption.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/gmxassert.h"

struct gmx_ana_indexgrps_t;

namespace gmx
{

/********************************************************************
 * ITopologyProvider
 */

ITopologyProvider::~ITopologyProvider() {}

/********************************************************************
 * SelectionOptionBehavior
 */

class SelectionOptionBehavior::Impl
{
public:
    Impl(SelectionCollection* selections, ITopologyProvider* topologyProvider) :
        selections_(*selections), topologyProvider_(*topologyProvider), manager_(selections), grps_(nullptr)
    {
    }
    ~Impl()
    {
        if (grps_ != nullptr)
        {
            gmx_ana_indexgrps_free(grps_);
        }
    }

    void promptSelections()
    {
        const bool isInteractive = StandardInputStream::isInteractive();
        initIndexGroups();
        manager_.parseRequestedFromStdin(isInteractive);
        doneIndexGroups();
    }
    void initIndexGroups()
    {
        if (!selections_.requiresIndexGroups() && !manager_.hasRequestedSelections())
        {
            if (!ndxfile_.empty())
            {
                std::fprintf(stderr,
                             "NOTE: You provided an index file\n"
                             "  %s\n(with -n), but it was not used by any selection.\n",
                             ndxfile_.c_str());
            }
            selections_.setIndexGroups(nullptr);
            return;
        }
        if (ndxfile_.empty())
        {
            gmx_mtop_t* top = topologyProvider_.getTopology(false);
            gmx_ana_indexgrps_init(&grps_, top, nullptr);
        }
        else
        {
            gmx_ana_indexgrps_init(&grps_, nullptr, ndxfile_.c_str());
        }
        selections_.setIndexGroups(grps_);
    }
    void doneIndexGroups()
    {
        if (grps_ != nullptr)
        {
            selections_.setIndexGroups(nullptr);
            gmx_ana_indexgrps_free(grps_);
            grps_ = nullptr;
        }
    }

    void compileSelections()
    {
        const bool  topRequired = selections_.requiredTopologyProperties().needsTopology_;
        gmx_mtop_t* top         = topologyProvider_.getTopology(topRequired);
        int         natoms      = -1;
        if (top == nullptr)
        {
            natoms = topologyProvider_.getAtomCount();
        }
        getMassesIfRequired(top);
        selections_.setTopology(top, natoms);
        selections_.compile();
        // Situation may have changed after compilation.
        getMassesIfRequired(top);
    }

    void getMassesIfRequired(gmx_mtop_t* top) const
    {
        const bool massRequired = selections_.requiredTopologyProperties().needsMasses_;
        if (!massRequired)
        {
            return;
        }
        // TODO: There can be some corner cases that still hit this assert
        // when the user has not provided the topology.
        GMX_RELEASE_ASSERT(top != nullptr, "Masses are required, but no topology is loaded");
        for (gmx_moltype_t& moltype : top->moltype)
        {
            if (!moltype.atoms.haveMass)
            {
                atomsSetMassesBasedOnNames(&moltype.atoms, TRUE);
                if (!moltype.atoms.haveMass)
                {
                    GMX_THROW(InconsistentInputError(
                            "Selections require mass information for evaluation, but it is not "
                            "available in the input and could not be determined for all atoms "
                            "based on atom names."));
                }
            }
        }
    }

    SelectionCollection&   selections_;
    ITopologyProvider&     topologyProvider_;
    SelectionOptionManager manager_;
    //! Name of the index file (empty if no index file provided).
    std::string          ndxfile_;
    gmx_ana_indexgrps_t* grps_;
};

SelectionOptionBehavior::SelectionOptionBehavior(SelectionCollection* selections,
                                                 ITopologyProvider*   topologyProvider) :
    impl_(new Impl(selections, topologyProvider))
{
}

SelectionOptionBehavior::~SelectionOptionBehavior() {}

void SelectionOptionBehavior::initOptions(IOptionsContainer* options)
{
    options->addOption(FileNameOption("n")
                               .filetype(OptionFileType::AtomIndex)
                               .inputFile()
                               .store(&impl_->ndxfile_)
                               .defaultBasename("index")
                               .description("Extra index groups"));
    options->addOption(SelectionFileOption("sf"));
    impl_->manager_.initOptions(options);
}

void SelectionOptionBehavior::initBehavior(Options* options)
{
    options->addManager(&impl_->manager_);
}

void SelectionOptionBehavior::optionsFinished()
{
    impl_->promptSelections();
    impl_->compileSelections();
}

} // namespace gmx
