/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */


#include "communication.h"

#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"

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
CommunicationStatus gmx_send_data(const t_commrec *cr, int dest)
{
    gmx_send_int(cr, dest, GMX_SEND_DATA);

    return CS_OK;
}

CommunicationStatus gmx_send_done(const t_commrec *cr, int dest)
{
    gmx_send_int(cr, dest, GMX_SEND_DONE);

    return CS_OK;
}

static CommunicationStatus gmx_recv_data_(const t_commrec *cr, int src, int line)
{
    int kk = gmx_recv_int(cr, src);

    if ((kk != GMX_SEND_DATA) && (kk != GMX_SEND_DONE))
    {
        gmx_fatal(FARGS, "Received %d in gmx_recv_data (line %d). Was expecting either %d or %d\n.", kk, line,
                  (int)GMX_SEND_DATA, (int)GMX_SEND_DONE);
    }
    return CS_OK;
}

CommunicationStatus gmx_recv_data(const t_commrec *cr, int src)
{
    return gmx_recv_data_(cr, src, __LINE__);
}
#undef GMX_SEND_DATA
#undef GMX_SEND_DONE
