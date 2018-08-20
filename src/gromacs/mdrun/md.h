/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2018, by the GROMACS development team, led by
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
 *
 * \brief Declares the integrator for normal molecular dynamics simulations
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_MD_H
#define GMX_MDRUN_MD_H

#include <memory>

#include "gromacs/mdrun/integrator.h"

namespace gmx
{

/*! \internal
 * \brief MD simulations
 *
 */
class MDIntegrator : public IIntegrator
{

    private:
        class Impl;

        std::unique_ptr<Impl> impl_;

    public:
        ~MDIntegrator() override;

        void run() override;

        class Builder;

        /*!
         * \brief Private constructor for use by Builder.
         *
         * \param implementation implementation object with which to instantiate the integrator.
         *
         * Create a new integrator object by transfering ownership of a MDIntegrator::Impl.
         */
        explicit MDIntegrator(std::unique_ptr<Impl> implementation);
};

class MDIntegrator::Builder : public IntegratorBuilder::Base
{
    public:
        Builder();
        ~Builder() override;

        Base &addContext(const md::Context &context) override;

        std::unique_ptr<IIntegrator> build() override;

        IntegratorBuilder::DataSentry setAggregateAdapter(std::unique_ptr<IntegratorAggregateAdapter> container)
        override;

    private:
        std::unique_ptr<IntegratorAggregateAdapter> container_;
        std::unique_ptr<md::Context>                context_;
};

}      // namespace gmx

#endif // GMX_MDRUN_MD_H
