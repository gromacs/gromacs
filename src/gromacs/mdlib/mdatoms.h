/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_MDATOMS_H
#define GMX_MDLIB_MDATOMS_H

#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/unique_cptr.h"

struct gmx_mtop_t;
struct t_inputrec;

namespace gmx
{

/*! \libinternal
 * \brief Contains a C-style t_mdatoms while managing some of its
 * memory with C++ vectors with allocators.
 *
 * The group-scheme kernels need to use a plain C-style t_mdatoms, so
 * this type combines that with the memory management needed for
 * efficient PME on GPU transfers.
 *
 * \todo Refactor this class and rename MDAtoms once the group scheme
 * is removed. */
class MDAtoms
{
    //! C-style mdatoms struct.
    unique_cptr<t_mdatoms> mdatoms_;
    //! Memory for chargeA that can be set up for efficient GPU transfer.
    HostVector<real>       chargeA_;
    public:
        // TODO make this private
        //! Constructor.
        MDAtoms();
        //! Getter.
        t_mdatoms *mdatoms()
        {
            return mdatoms_.get();
        }
        /*! \brief Resizes memory.
         *
         * \throws std::bad_alloc  If out of memory.
         */
        void resize(int newSize);
        /*! \brief Reserves memory.
         *
         * \throws std::bad_alloc  If out of memory.
         */
        void reserve(int newCapacity);
        //! Builder function.
        friend std::unique_ptr<MDAtoms>
        makeMDAtoms(FILE *fp, const gmx_mtop_t &mtop, const t_inputrec &ir,
                    bool useGpuForPme);
};

//! Builder function for MdAtomsWrapper.
std::unique_ptr<MDAtoms>
makeMDAtoms(FILE *fp, const gmx_mtop_t &mtop, const t_inputrec &ir,
            bool useGpuForPme);

} // namespace

void atoms2md(const gmx_mtop_t *mtop, const t_inputrec *ir,
              int nindex, const int *index,
              int homenr,
              gmx::MDAtoms *mdAtoms);
/* This routine copies the atoms->atom struct into md.
 * If index!=NULL only the indexed atoms are copied.
 * For the masses the A-state (lambda=0) mass is used.
 * Sets md->lambda = 0.
 * In free-energy runs, update_mdatoms() should be called after atoms2md()
 * to set the masses corresponding to the value of lambda at each step.
 */

void update_mdatoms(t_mdatoms *md, real lambda);
/* When necessary, sets all the mass parameters to values corresponding
 * to the free-energy parameter lambda.
 * Sets md->lambda = lambda.
 */

#endif
