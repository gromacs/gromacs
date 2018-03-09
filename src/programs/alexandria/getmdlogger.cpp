/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
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
 
 
#include "gmxpre.h"

#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/loggerbuilder.h"

gmx::MDLogger getMdLogger(const t_commrec *cr,
                          FILE            *fplog)
{
    gmx::LoggerBuilder builder;
    builder.addTargetFile(gmx::MDLogger::LogLevel::Info, fplog);
    if (cr == nullptr || SIMMASTER(cr))
    {
        builder.addTargetStream(gmx::MDLogger::LogLevel::Warning,
                                &gmx::TextOutputFile::standardError());
    }
    gmx::LoggerOwner logOwner(builder.build());
    gmx::MDLogger    mdlog(logOwner.logger());

    return mdlog;
}
