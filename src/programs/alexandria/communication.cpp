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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
 
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"

#include "communication.h"
#include "gmx_simple_comm.h"
    
    
const char *cs_name(CommunicationStatus cs)
{
    switch (cs)
    {
        case CS_OK:
            static const char *ok = "Communication OK";
            return ok;
        case CS_ERROR:
            static const char *err = "Communication Error";
            return err;
        case CS_SEND_DATA:
            static const char *sd = "Communication sent data";
            return sd;
        case CS_RECV_DATA:
            static const char *rd = "Communication OK";
            return rd;
        default:
            gmx_fatal(FARGS, "Unknown communication status %d", (int) cs);
    }
    return nullptr;
};

#define GMX_SEND_DATA 19823
#define GMX_SEND_DONE -666
CommunicationStatus gmx_send_data(t_commrec *cr, int dest)
{
    gmx_send_int(cr, dest, GMX_SEND_DATA);

    return CS_OK;
}

CommunicationStatus gmx_send_done(t_commrec *cr, int dest)
{
    gmx_send_int(cr, dest, GMX_SEND_DONE);

    return CS_OK;
}

static CommunicationStatus gmx_recv_data_(t_commrec *cr, int src, int line)
{
    int kk = gmx_recv_int(cr, src);

    if ((kk != GMX_SEND_DATA) && (kk != GMX_SEND_DONE))
    {
        gmx_fatal(FARGS, "Received %d in gmx_recv_data (line %d). Was expecting either %d or %d\n.", kk, line,
                  (int)GMX_SEND_DATA, (int)GMX_SEND_DONE);
    }
    return CS_OK;
}

CommunicationStatus gmx_recv_data(t_commrec *cr, int src) 
{
    return gmx_recv_data_(cr, src, __LINE__);
}
#undef GMX_SEND_DATA
#undef GMX_SEND_DONE

