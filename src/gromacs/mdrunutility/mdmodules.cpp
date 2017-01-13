/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "mdmodules.h"

#include "gromacs/applied-forces/electricfield.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

class MDModules::Impl
{
public:
    Impl()
        : ir_(nullptr)
    {
    }
    ~Impl()
    {
        if (ir_ != nullptr)
        {
            done_inputrec(ir_);
            sfree(ir_);
        }
    }

    void ensureInputrecInitialized()
    {
        if (ir_ == nullptr)
        {
            field_ = createElectricFieldModule();
            snew(ir_, 1);
            snew(ir_->fepvals, 1);
            snew(ir_->expandedvals, 1);
            snew(ir_->simtempvals, 1);
            ir_->efield = field_.get();
        }
    }

    std::unique_ptr<IInputRecExtension> field_;
    t_inputrec *                        ir_;
};

MDModules::MDModules()
    : impl_(new Impl)
{
}

MDModules::~MDModules()
{
}

t_inputrec *MDModules::inputrec()
{
    impl_->ensureInputrecInitialized();
    return impl_->ir_;
}

} // namespace gmx
