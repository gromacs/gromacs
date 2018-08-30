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
/*! \file
 * \brief
 * Declares gmx::ProcessFrameConversion
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinatedata
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_ANALYSE_H
#define GMX_TRAJECTORYANALYSIS_MODULES_ANALYSE_H

#include <algorithm>
#include <utility>

#include "gromacs/coordinatedata/frameconverters/frameconverter.h"

namespace gmx
{

/*!\brief
 * ProcessFrameConversion class for handling the running of several analysis steps
 *
 * This analysis module allows to register modules for coordinate frame manipulation
 * that are then run ones this modules convertFrame method is invoked.
 *
 * \inlibraryapi
 * \ingroup module_coordinatedata
 *
 */
class ProcessFrameConversion : public IFrameConverter
{
    public:
        /*! \brief
         * Default constructor for ProcessFrameConversion.
         */
        ProcessFrameConversion() {}

        virtual ~ProcessFrameConversion() {}

        /*! \brief
         * Change coordinate frame information for output.
         *
         * This method is used to perform the actual coordinate frame manipulation.
         * In this case, it acts as a wrapper that runs the method for all
         * modules that have been registered for the analysis chain.
         *
         * \param[in] input Coordinate frame to be modified later.
         */
        virtual void convertFrame(const t_trxframe &input);

        /*! \brief
         * Add framemodule to analysis chain.
         *
         * Other modules derived from IFrameConverter can be registered to be analysed
         * by this method instead of having to be analysed separately.
         *
         * \param[in] module Framemodification module to add to the chain.
         * \throws unspecified Any exception thrown by any of the \p module
         * objects in the chain during analysis.
         */
        void addFrameConverter(FrameConverterPointer module);

        /*! \libinternal \brief
         * Storage for info about different modules chained together.
         *
         * This is storing the pointers to the individual methods in the analysis chain.
         * For each method, one pointer to the method is stored to be used in the custom
         * convertFrame method.
         *
         */
        struct FrameModule
        {
            //! Initializes module, stolen from datamodulemanager.
            explicit FrameModule(FrameConverterPointer module)
                : module(std::move(module))
            {
            }
            //! Pointer to module.
            FrameConverterPointer module;
        };
        //! Shorthand for list of chained modules
        typedef std::vector<FrameModule> FrameModuleList;

        //! List of chained modules.
        FrameModuleList moduleChain_;

};

//! Smart pointer to manage the analyse object.
typedef std::unique_ptr<ProcessFrameConversion>
    ProcessFrameConversionPointer;

} // namespace gmx

#endif
