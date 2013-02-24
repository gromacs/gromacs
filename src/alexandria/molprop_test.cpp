/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: mptest.c,v 1.5 2009/04/05 11:46:57 spoel Exp $
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
 *
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
#include <stdio.h>
#include <stdlib.h>
#include <gmx_fatal.h>
#include "atomprop.h"
#include "poldata_xml.h"

#include <vector>
#include "molprop.hpp"
#include "molprop_xml.hpp"

int main(int argc,char*argv[])
{
    std::vector<alexandria::MolProp> mpt;
    
    if (argc < 3) 
    {
        fprintf(stderr,"Usage: %s infile outfile\n",argv[0]);
        exit(1);
    }
    MolPropRead(argv[1],mpt);
    printf("Read %lu molecules from %s\n",mpt.size(),argv[1]);
    MolPropWrite(argv[2],mpt,FALSE);
    
    return 0;
}

