/* -*- mode: c; tab-width: 4; indent-tabs-mode: n; c-basic-offset: 4 -*- 
 *
 * $Id$
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#include <stdio.h>

#include <types/nrnb.h>

#include "../nb_kerneltype.h"

#include "nb_kernel.h"

/* Include standard kernel headers in local directory */
#include "nb_kernel010.h"
#include "nb_kernel020.h"
#include "nb_kernel030.h"
#include "nb_kernel100.h"
#include "nb_kernel101.h"
#include "nb_kernel102.h"
#include "nb_kernel103.h"
#include "nb_kernel104.h"
#include "nb_kernel110.h"
#include "nb_kernel111.h"
#include "nb_kernel112.h"
#include "nb_kernel113.h"
#include "nb_kernel114.h"
#include "nb_kernel120.h"
#include "nb_kernel121.h"
#include "nb_kernel122.h"
#include "nb_kernel123.h"
#include "nb_kernel124.h"
#include "nb_kernel130.h"
#include "nb_kernel131.h"
#include "nb_kernel132.h"
#include "nb_kernel133.h"
#include "nb_kernel134.h"
#include "nb_kernel200.h"
#include "nb_kernel201.h"
#include "nb_kernel202.h"
#include "nb_kernel203.h"
#include "nb_kernel204.h"
#include "nb_kernel210.h"
#include "nb_kernel211.h"
#include "nb_kernel212.h"
#include "nb_kernel213.h"
#include "nb_kernel214.h"
#include "nb_kernel220.h"
#include "nb_kernel221.h"
#include "nb_kernel222.h"
#include "nb_kernel223.h"
#include "nb_kernel224.h"
#include "nb_kernel230.h"
#include "nb_kernel231.h"
#include "nb_kernel232.h"
#include "nb_kernel233.h"
#include "nb_kernel234.h"
#include "nb_kernel300.h"
#include "nb_kernel301.h"
#include "nb_kernel302.h"
#include "nb_kernel303.h"
#include "nb_kernel304.h"
#include "nb_kernel310.h"
#include "nb_kernel311.h"
#include "nb_kernel312.h"
#include "nb_kernel313.h"
#include "nb_kernel314.h"
#include "nb_kernel320.h"
#include "nb_kernel321.h"
#include "nb_kernel322.h"
#include "nb_kernel323.h"
#include "nb_kernel324.h"
#include "nb_kernel330.h"
#include "nb_kernel331.h"
#include "nb_kernel332.h"
#include "nb_kernel333.h"
#include "nb_kernel334.h"
#include "nb_kernel400.h"
#include "nb_kernel410.h"
#include "nb_kernel420.h"
#include "nb_kernel430.h"


static nb_kernel_t *
kernellist[eNR_NBKERNEL_NR] = 
{
    F77_OR_C_FUNC_(nb_kernel010,NB_KERNEL010),
    F77_OR_C_FUNC_(nb_kernel020,NB_KERNEL020),
    F77_OR_C_FUNC_(nb_kernel030,NB_KERNEL030),
    F77_OR_C_FUNC_(nb_kernel100,NB_KERNEL100),
    F77_OR_C_FUNC_(nb_kernel101,NB_KERNEL101),
    F77_OR_C_FUNC_(nb_kernel102,NB_KERNEL102),
    F77_OR_C_FUNC_(nb_kernel103,NB_KERNEL103),
    F77_OR_C_FUNC_(nb_kernel104,NB_KERNEL104),
    F77_OR_C_FUNC_(nb_kernel110,NB_KERNEL110),
    F77_OR_C_FUNC_(nb_kernel111,NB_KERNEL111),
    F77_OR_C_FUNC_(nb_kernel112,NB_KERNEL112),
    F77_OR_C_FUNC_(nb_kernel113,NB_KERNEL113),
    F77_OR_C_FUNC_(nb_kernel114,NB_KERNEL114),
    F77_OR_C_FUNC_(nb_kernel120,NB_KERNEL120),
    F77_OR_C_FUNC_(nb_kernel121,NB_KERNEL121),
    F77_OR_C_FUNC_(nb_kernel122,NB_KERNEL122),
    F77_OR_C_FUNC_(nb_kernel123,NB_KERNEL123),
    F77_OR_C_FUNC_(nb_kernel124,NB_KERNEL124),
    F77_OR_C_FUNC_(nb_kernel130,NB_KERNEL130),
    F77_OR_C_FUNC_(nb_kernel131,NB_KERNEL131),
    F77_OR_C_FUNC_(nb_kernel132,NB_KERNEL132),
    F77_OR_C_FUNC_(nb_kernel133,NB_KERNEL133),
    F77_OR_C_FUNC_(nb_kernel134,NB_KERNEL134),
    F77_OR_C_FUNC_(nb_kernel200,NB_KERNEL200),
    F77_OR_C_FUNC_(nb_kernel201,NB_KERNEL201),
    F77_OR_C_FUNC_(nb_kernel202,NB_KERNEL202),
    F77_OR_C_FUNC_(nb_kernel203,NB_KERNEL203),
    F77_OR_C_FUNC_(nb_kernel204,NB_KERNEL204),
    F77_OR_C_FUNC_(nb_kernel210,NB_KERNEL210),
    F77_OR_C_FUNC_(nb_kernel211,NB_KERNEL211),
    F77_OR_C_FUNC_(nb_kernel212,NB_KERNEL212),
    F77_OR_C_FUNC_(nb_kernel213,NB_KERNEL213),
    F77_OR_C_FUNC_(nb_kernel214,NB_KERNEL214),
    F77_OR_C_FUNC_(nb_kernel220,NB_KERNEL220),
    F77_OR_C_FUNC_(nb_kernel221,NB_KERNEL221),
    F77_OR_C_FUNC_(nb_kernel222,NB_KERNEL222),
    F77_OR_C_FUNC_(nb_kernel223,NB_KERNEL223),
    F77_OR_C_FUNC_(nb_kernel224,NB_KERNEL224),
    F77_OR_C_FUNC_(nb_kernel230,NB_KERNEL230),
    F77_OR_C_FUNC_(nb_kernel231,NB_KERNEL231),
    F77_OR_C_FUNC_(nb_kernel232,NB_KERNEL232),
    F77_OR_C_FUNC_(nb_kernel233,NB_KERNEL233),
    F77_OR_C_FUNC_(nb_kernel234,NB_KERNEL234),
    F77_OR_C_FUNC_(nb_kernel300,NB_KERNEL300),
    F77_OR_C_FUNC_(nb_kernel301,NB_KERNEL301),
    F77_OR_C_FUNC_(nb_kernel302,NB_KERNEL302),
    F77_OR_C_FUNC_(nb_kernel303,NB_KERNEL303),
    F77_OR_C_FUNC_(nb_kernel304,NB_KERNEL304),
    F77_OR_C_FUNC_(nb_kernel310,NB_KERNEL310),
    F77_OR_C_FUNC_(nb_kernel311,NB_KERNEL311),
    F77_OR_C_FUNC_(nb_kernel312,NB_KERNEL312),
    F77_OR_C_FUNC_(nb_kernel313,NB_KERNEL313),
    F77_OR_C_FUNC_(nb_kernel314,NB_KERNEL314),
    F77_OR_C_FUNC_(nb_kernel320,NB_KERNEL320),
    F77_OR_C_FUNC_(nb_kernel321,NB_KERNEL321),
    F77_OR_C_FUNC_(nb_kernel322,NB_KERNEL322),
    F77_OR_C_FUNC_(nb_kernel323,NB_KERNEL323),
    F77_OR_C_FUNC_(nb_kernel324,NB_KERNEL324),
    F77_OR_C_FUNC_(nb_kernel330,NB_KERNEL330),
    F77_OR_C_FUNC_(nb_kernel331,NB_KERNEL331),
    F77_OR_C_FUNC_(nb_kernel332,NB_KERNEL332),
    F77_OR_C_FUNC_(nb_kernel333,NB_KERNEL333),
    F77_OR_C_FUNC_(nb_kernel334,NB_KERNEL334),
    F77_OR_C_FUNC_(nb_kernel400,NB_KERNEL400),
    F77_OR_C_FUNC_(nb_kernel410,NB_KERNEL410),
    F77_OR_C_FUNC_(nb_kernel430,NB_KERNEL430)
};


void
nb_kernel_setup(FILE *log,nb_kernel_t **list)
{
	int i;
    nb_kernel_t *p;

	for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist[i];
        if(p!=NULL)
            list[i] = p; 
	}
}    
