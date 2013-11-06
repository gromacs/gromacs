/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Declares code for calling GROMACS routines from an external program
 *
 * \author David van der Spoel <spoel@xray.bmc.uu.se>
 * \inpublicapi
 * \ingroup module_qmmm
 */
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/mm_slave.h"
#include "gromacs/mmslave/mm_slave.h"

//! Abstract type for the mm_slave code
typedef struct {
    //! Embedded C++ class
    MMSlave *mms;
} gmx_mm_slave;

/* Routines for C interface to the MMSlave class */
gmx_mm_slave_t mm_slave_init(void)
{
    gmx_mm_slave *gms;
    
    snew(gms,1);
    gms->mms = new MMSlave();
    
    return gms;
}
    
void mm_slave_done(gmx_mm_slave_t gms)
{
    delete gms->mms;
    sfree(gms);
}

gmx_bool mm_slave_read_tpr(const char *tpr, 
                           gmx_mm_slave_t gms)
{
    return (gmx_bool) gms->mms->readTpr(tpr);
}

gmx_bool mm_slave_set_q(gmx_mm_slave_t gms,
                        atom_id id,
                        double q)
{
    return (gmx_bool) gms->mms->setAtomQ(id, q);
}

gmx_bool mm_slave_calc_energy(gmx_mm_slave_t gms, 
                              const rvec *x,
                              rvec *f,
                              double *energy)
{
    return (gmx_bool) gms->mms->calcEnergy(x, f, energy);
}
    
namespace gmx
{

bool MMSlave::readTpr(const char *tpr)
{
    return true;
}

bool MMSlave::setAtomQ(atom_id id, double q)
{
    return true;
}

bool MMSlave::calcEnergy(const rvec *x,
                         rvec *f,
                         double *energy)
{
    return true;
}

}
    
#endif
