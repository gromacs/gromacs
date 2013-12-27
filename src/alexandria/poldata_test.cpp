/*
 * $Id: pdtest.c,v 1.5 2009/04/05 11:46:58 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#include <stdlib.h>
#include "gromacs/legacyheaders/statutil.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/legacyheaders/macros.h"
#include "poldata_xml.hpp"
#include "atomprop.h"

int alex_poldata_test(int argc,char*argv[])
{
    static const char               *desc[] = {
        "poldata_test reads a poldata (force field) file and writes a new one.",
    };
    output_env_t   oenv;    
    t_filenm        fnm[] = {
        { efDAT, "-f", "pdin", ffREAD },
        { efDAT, "-o", "pdout", ffWRITE }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, 0, NULL,
                           asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }
    
    gmx_atomprop_t aps = gmx_atomprop_init();
    gmx_poldata_t pd   = gmx_poldata_read(opt2fn("-f", NFILE, fnm), aps);
    gmx_poldata_write(opt2fn("-o", NFILE, fnm), pd, 0);
    
    return 0;
}
