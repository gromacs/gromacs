/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
/*! \internal
 * \brief Encapsulates membed methods
 *
 * \author Joe Jordan <ejjordan@kth.se>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "membedholder.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/mdlib/membed.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/real.h"

namespace gmx
{

MembedHolder::MembedHolder(int nfile, const t_filenm fnm[]) :
    doMembed_(opt2bSet("-membed", nfile, fnm))
{
}

MembedHolder::~MembedHolder()
{
    if (doMembed_)
    {
        free_membed(membed_);
    }
}

void MembedHolder::initializeMembed(FILE*          fplog,
                                    int            nfile,
                                    const t_filenm fnm[],
                                    gmx_mtop_t*    mtop,
                                    t_inputrec*    inputrec,
                                    t_state*       state,
                                    t_commrec*     cr,
                                    real*          cpt)
{
    if (doMembed_)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "Initializing membed");
        }
        /* Note that membed cannot work in parallel because mtop is
         * changed here. Fix this if we ever want to make it run with
         * multiple ranks. */
        membed_ = init_membed(fplog, nfile, fnm, mtop, inputrec, state, cr, cpt);
    }
}

gmx_membed_t* MembedHolder::membed()
{
    return membed_;
}

MembedHolder::MembedHolder(MembedHolder&& holder) noexcept
{
    doMembed_        = holder.doMembed_;
    membed_          = holder.membed_;
    holder.membed_   = nullptr;
    holder.doMembed_ = false;
}

MembedHolder& MembedHolder::operator=(MembedHolder&& holder) noexcept
{
    if (&holder != this)
    {
        doMembed_        = holder.doMembed_;
        membed_          = holder.membed_;
        holder.membed_   = nullptr;
        holder.doMembed_ = false;
    }
    return *this;
}

} // namespace gmx
