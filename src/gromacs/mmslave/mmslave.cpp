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
 * Contains code for calling GROMACS routines from an external program
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 * \ingroup module_qmmm
 */
#include <stdlib.h>
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/mtop_util.h"
#include "gromacs/legacyheaders/vec.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/mmslave.h"
#include "gromacs/mmslave/mmslave.h"

/* Note: the C-interface is all the way down in this file */

namespace gmx
{

MMSlave::MMSlave()
{
    x_      = NULL;
    v_      = NULL;
    f_      = NULL;
    natoms_ = 0;
}

bool MMSlave::readTpr(const char *tpr)
{
    t_tpxheader tpx;
    int         version, generation;
    
    read_tpxheader(tpr, &tpx, FALSE, &version, &generation);
    natoms_ = tpx.natoms;
    x_ = (rvec *)calloc(natoms_, sizeof(rvec));
    v_ = (rvec *)calloc(natoms_, sizeof(rvec));
    f_ = (rvec *)calloc(natoms_, sizeof(rvec));
    (void) read_tpx(tpr, &inputrec_, box_, &natoms_, x_, v_, f_, &mtop_);

    return true;
}

static bool copyIt(int natoms_dst, rvec *x_dst, int natoms_src, rvec *x_src)
{
    if ((natoms_dst < natoms_src) || (NULL == x_dst))
    {
        return false;
    }
    for(int i=0; (i < natoms_src); i++)
    {
        copy_rvec(x_src[i], x_dst[i]);
    }
    return true;
}

bool MMSlave::copyX(int natoms, rvec *x)
{
    return copyIt(natoms, x, natoms_, x_);
}

bool MMSlave::copyV(int natoms, rvec *v)
{
    return copyIt(natoms, v, natoms_, v_);
}

bool MMSlave::copyF(int natoms, rvec *f)
{
    return copyIt(natoms, f, natoms_, f_);
}

void MMSlave::cleanUp()
{
    if (NULL != x_)
    {
        free(x_);
    }
    if (NULL != v_)
    {
        free(v_);
    }
    if (NULL != f_)
    {
        free(f_);
    }
}

bool MMSlave::setAtomQ(atom_id id, double q)
{
    t_atom *atom;
    
    gmx_mtop_atomlookup_t alook = gmx_mtop_atomlookup_init(&mtop_);
    
    gmx_mtop_atomnr_to_atom(alook, id, &atom);
    atom->q = q;
    atom->qB = q;

    gmx_mtop_atomlookup_destroy(alook);

    return true;
}

bool MMSlave::calcEnergy(const rvec *x,
                         rvec *f,
                         double *energy)
{
    fprintf(stderr, "Warning: computing the force and energy is not implemented yet.\n");
    return true;
}

}

//! Abstract type for the mmslave code
typedef struct gmx_mmslave {
    //! Embedded C++ class
    gmx::MMSlave *mms;
} gmx_mmslave;

/* Routines for C interface to the MMSlave class */
gmx_mmslave_t mmslave_init(void)
{
    gmx_mmslave *gms;
    
    gms = (gmx_mmslave *) calloc(1,sizeof(gmx_mmslave));
    gms->mms = new gmx::MMSlave();
    
    return gms;
}
    
void mmslave_done(gmx_mmslave_t gms)
{
    delete gms->mms;
    free(gms);
}

gmx_bool mmslave_read_tpr(const char *tpr, 
                          gmx_mmslave_t gms)
{
    return (gmx_bool) gms->mms->readTpr(tpr);
}

int mmslave_natoms(gmx_mmslave_t gms)
{
    return gms->mms->nAtoms();
}

gmx_bool mmslave_copyX(gmx_mmslave_t gms, int natoms, rvec *x)
{
    return (gmx_bool) gms->mms->copyX(natoms, x);
}
    
gmx_bool mmslave_copyV(gmx_mmslave_t gms, int natoms, rvec *v)
{
    return (gmx_bool) gms->mms->copyX(natoms, v);
}
    
gmx_bool mmslave_copyF(gmx_mmslave_t gms, int natoms, rvec *f)
{
    return (gmx_bool) gms->mms->copyX(natoms, f);
}
    
void mmslave_clean(gmx_mmslave_t gms)
{
    gms->mms->cleanUp();
    delete gms->mms;
    free(gms);
}

gmx_bool mmslave_set_q(gmx_mmslave_t gms,
                       atom_id id,
                       double q)
{
    return (gmx_bool) gms->mms->setAtomQ(id, q);
}

gmx_bool mmslave_calc_energy(gmx_mmslave_t gms, 
                             const rvec *x,
                             rvec *f,
                             double *energy)
{
    return (gmx_bool) gms->mms->calcEnergy(x, f, energy);
}
    
