/*
 * This source file is part of the Aleandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"
#include <stdlib.h>
#include "gromacs/commandline/pargs.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/topology/atomprop.h"
#include "poldata_xml.h"

int alex_poldata_test(int argc, char*argv[])
{
    static const char               *desc[] = {
        "poldata_test reads a poldata (force field) file and writes a new one.",
    };
    output_env_t                     oenv;
    t_filenm                         fnm[] = {
        { efDAT, "-f", "pdin", ffREAD },
        { efDAT, "-o", "pdout", ffWRITE }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, 0, NULL,
                           asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    gmx_atomprop_t aps  = gmx_atomprop_init();
    gmx_poldata_t  pd   = gmx_poldata_read(opt2fn("-f", NFILE, fnm), aps);
    gmx_poldata_write(opt2fn("-o", NFILE, fnm), pd, 0);

    return 0;
}
