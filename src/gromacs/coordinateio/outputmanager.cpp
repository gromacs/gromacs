/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*!\file
 * \internal
 * \brief
 * Implements methods from outputmanager.h
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */


#include "gmxpre.h"

#include "outputmanager.h"

#include <algorithm>

#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/*! \brief
 * Create a deep copy of a t_trxframe \p input into \p copy
 *
 * When running the analysis tools and changing values with the
 * outputadapters, a deep copy of the \p input coordinate frame has to be
 * created first to ensure that the data is not changed if it is needed for other
 * tools following with analysis later. Therefore, the data is passed
 * to \p copy by performing a deep copy first.
 *
 * The method allocates new storage for coordinates of the x, v, and f arrays
 * in the new coordinate frame. This means that those arrays need to be free'd
 * after the frame has been processed and been written to disk.
 *
 * \param[in]     input Reference input coordinate frame.
 * \param[in,out] copy  Pointer to new output frame that will receive the deep copy.
 */
static void deepCopy_t_trxframe(const t_trxframe &input, t_trxframe *copy)
{
    copy->not_ok    = input.not_ok;
    copy->bStep     = input.bStep;
    copy->bTime     = input.bTime;
    copy->bLambda   = input.bLambda;
    copy->bFepState = input.bFepState;
    copy->bAtoms    = input.bAtoms;
    copy->bPrec     = input.bPrec;
    copy->bX        = input.bX;
    copy->bV        = input.bV;
    copy->bF        = input.bF;
    copy->bBox      = input.bBox;
    copy->bDouble   = input.bDouble;
    copy->natoms    = input.natoms;
    copy->step      = input.step;
    copy->time      = input.time;
    copy->lambda    = input.lambda;
    copy->fep_state = input.fep_state;
    if (input.bAtoms)
    {
        copy->atoms = copy_t_atoms(input.atoms);
    }
    copy->atoms     = input.atoms;
    copy->prec      = input.prec;
    if (copy->bX)
    {
        snew(copy->x, copy->natoms);
    }
    if (copy->bV)
    {
        snew(copy->v, copy->natoms);
    }
    if (copy->bF)
    {
        snew(copy->f, copy->natoms);
    }
    for (int i = 0; i < copy->natoms; i++)
    {
        if (copy->bX)
        {
            copy_rvec(input.x[i], copy->x[i]);
        }
        if (copy->bV)
        {
            copy_rvec(input.v[i], copy->v[i]);
        }
        if (copy->bF)
        {
            copy_rvec(input.f[i], copy->f[i]);
        }
    }
    copy_mat(input.box, copy->box);
    copy->bPBC   = input.bPBC;
    copy->ePBC   = input.ePBC;
}

//! Clean up local storage of atom information.
static void clearLocalAtoms(t_atoms *atoms)
{
    done_atom(atoms);
    sfree(atoms);
}

/*! \brief
 * Method to open TNG file.
 *
 * Only need extra method to open this kind of file as it may need access to
 * a Selection \p sel, if it is valid. Otherwise atom indices will be taken
 * from the topology \p mtop.
 *
 * \param[in] name Name of the output file.
 * \param[in] sel Reference to selection.
 * \param[in] mtop Pointer to topology, tested before that it is valid.
 */
static t_trxstatus *openTNG(const std::string &name, const Selection &sel, const gmx_mtop_t *mtop)
{
    const char *filemode = "w";
    if (sel.isValid())
    {
        GMX_ASSERT(sel.hasOnlyAtoms(), "Can only work with selections consisting out of atoms");
        return trjtools_gmx_prepare_tng_writing_old(name.c_str(),
                                                    filemode[0],
                                                    nullptr, //infile_, //how to get the input file here?
                                                    nullptr,
                                                    sel.atomCount(),
                                                    const_cast<gmx_mtop_t *>(mtop),
                                                    sel.atomIndices().data(),
                                                    sel.name());
    }
    else
    {
        return trjtools_gmx_prepare_tng_writing(name.c_str(),
                                                filemode[0],
                                                nullptr, //infile_, //how to get the input file here?
                                                nullptr,
                                                mtop->natoms,
                                                const_cast<gmx_mtop_t *>(mtop),
                                                get_atom_index(mtop).data(),
                                                "System");
    }
}

void
OutputManager::canAddNewModule(unsigned long newModule)
{
    if ((registeredModules_ & newModule) == newModule)
    {
        GMX_THROW(InternalError("Trying to add module that has already been added"));
    }
    else if ( (newModule == efChangeAtomInformationModule) &&
              ((registeredModules_ & efChangeCoordinateSelectionModule)
               == efChangeCoordinateSelectionModule))
    {
        std::string message = formatString("Trying to add module that changes atom "
                                           "information after changing output selection, this is not supported.");
        GMX_THROW(InternalError(message));
    }
    if ((outputFlag_ & newModule) == 0u)
    {
        std::string message = formatString("Adding module will violate requirements for "
                                           "output to trajectory file.");
        GMX_THROW(InconsistentInputError(message));
    }
    registeredModules_ |= newModule;
}

void
OutputManager::addOutputAdapter(CoordinateOutputPointer module)
{
    // first check if the module dependencies are satisfied
    // by the OutputManager object that has been constructed.
    canAddNewModule(module->getModuleFlag());
    outputAdapters_.emplace_back(std::move(module));
}

void
OutputManager::closeFile()
{
    if (outputFile_ != nullptr)
    {
        close_trx(outputFile_);
    }
    outputFile_ = nullptr;
}

void
OutputManager::initOutput()
{
    if (outputFile_ == nullptr)
    {
        const char *filemode = "w";
        switch (filetype_)
        {
            case (efTNG):
                outputFile_ = openTNG(outputFileName_,
                                      sel_,
                                      getMtop());
                break;
            case (efPDB):
            case (efGRO):
            case (efTRR):
            case (efXTC):
            case (efG96):
                outputFile_ = open_trx(outputFileName_.c_str(), filemode);
                break;
            default:
                GMX_THROW(InvalidInputError("Invalid file type"));
        }
    }
}

void
OutputManager::prepareFrame(const int framenumber, const t_trxframe &input)
{
    if (outputFile_ == nullptr)
    {
        initOutput();
    }
    if (!outputAdapters_.empty())
    {
        t_trxframe local;
        clear_trxframe(&local, true);
        deepCopy_t_trxframe(input, &local);
        for (const auto &outputAdapter : outputAdapters_)
        {
            outputAdapter.module_->processFrame(framenumber, &local);
        }
        write_trxframe(outputFile_, &local, nullptr);
        clearLocalAtoms(local.atoms);
        clear_trxframe(&local, true);
    }
    else
    {
        write_trxframe(outputFile_, const_cast<t_trxframe *>(&input), nullptr);
    }
}

} // namespace gmx
