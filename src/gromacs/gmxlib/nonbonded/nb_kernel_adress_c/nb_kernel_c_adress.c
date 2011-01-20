
/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000_adress, University of Groningen_adress, The Netherlands.
 * Copyright (c) 2001-2009_adress, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis_adress,
 * written by Erik Lindahl_adress, David van der Spoel_adress, Berk Hess_adress, and others - for
 * a full list of developers and information_adress, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2 of the License_adress, or (at your option) any
 * later version.
 * As a special exception_adress, you may use this file as part of a free software
 * library without restriction.  Specifically_adress, if other files instantiate
 * templates or use macros or inline functions from this file_adress, or you compile
 * this file and link it with other files to produce an executable_adress, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL_adress, not an application linking with it.
 *
 * To help fund GROMACS development_adress, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include "types/nrnb.h"
#include "nb_kernel_c_adress.h"
#include "../nb_kerneltype.h"


/* Include standard kernel headers in local directory */
#include "nb_kernel010_adress.h"
#include "nb_kernel020_adress.h"
#include "nb_kernel030_adress.h"
#include "nb_kernel100_adress.h"
#include "nb_kernel101_adress.h"
#include "nb_kernel102_adress.h"
#include "nb_kernel103_adress.h"
#include "nb_kernel104_adress.h"
#include "nb_kernel110_adress.h"
#include "nb_kernel111_adress.h"
#include "nb_kernel112_adress.h"
#include "nb_kernel113_adress.h"
#include "nb_kernel114_adress.h"
#include "nb_kernel120_adress.h"
#include "nb_kernel121_adress.h"
#include "nb_kernel122_adress.h"
#include "nb_kernel123_adress.h"
#include "nb_kernel124_adress.h"
#include "nb_kernel130_adress.h"
#include "nb_kernel131_adress.h"
#include "nb_kernel132_adress.h"
#include "nb_kernel133_adress.h"
#include "nb_kernel134_adress.h"
#include "nb_kernel200_adress.h"
#include "nb_kernel201_adress.h"
#include "nb_kernel202_adress.h"
#include "nb_kernel203_adress.h"
#include "nb_kernel204_adress.h"
#include "nb_kernel210_adress.h"
#include "nb_kernel211_adress.h"
#include "nb_kernel212_adress.h"
#include "nb_kernel213_adress.h"
#include "nb_kernel214_adress.h"
#include "nb_kernel220_adress.h"
#include "nb_kernel221_adress.h"
#include "nb_kernel222_adress.h"
#include "nb_kernel223_adress.h"
#include "nb_kernel224_adress.h"
#include "nb_kernel230_adress.h"
#include "nb_kernel231_adress.h"
#include "nb_kernel232_adress.h"
#include "nb_kernel233_adress.h"
#include "nb_kernel234_adress.h"
#include "nb_kernel300_adress.h"
#include "nb_kernel301_adress.h"
#include "nb_kernel302_adress.h"
#include "nb_kernel303_adress.h"
#include "nb_kernel304_adress.h"
#include "nb_kernel310_adress.h"
#include "nb_kernel311_adress.h"
#include "nb_kernel312_adress.h"
#include "nb_kernel313_adress.h"
#include "nb_kernel314_adress.h"
#include "nb_kernel320_adress.h"
#include "nb_kernel321_adress.h"
#include "nb_kernel322_adress.h"
#include "nb_kernel323_adress.h"
#include "nb_kernel324_adress.h"
#include "nb_kernel330_adress.h"
#include "nb_kernel331_adress.h"
#include "nb_kernel332_adress.h"
#include "nb_kernel333_adress.h"
#include "nb_kernel334_adress.h"
#include "nb_kernel400_adress.h"
#include "nb_kernel410_adress.h"
#include "nb_kernel420_adress.h"
#include "nb_kernel430_adress.h"



static nb_adress_kernel_t *
kernellist_adress[eNR_NBKERNEL_NR] =
{
    nb_kernel010_adress_cg,
    nb_kernel020_adress_cg,
    nb_kernel030_adress_cg,
    nb_kernel100_adress_cg,
    nb_kernel101_adress_cg,
    nb_kernel102_adress_cg,
    nb_kernel103_adress_cg,
    nb_kernel104_adress_cg,
    nb_kernel110_adress_cg,
    nb_kernel111_adress_cg,
    nb_kernel112_adress_cg,
    nb_kernel113_adress_cg,
    nb_kernel114_adress_cg,
    nb_kernel120_adress_cg,
    nb_kernel121_adress_cg,
    nb_kernel122_adress_cg,
    nb_kernel123_adress_cg,
    nb_kernel124_adress_cg,
    nb_kernel130_adress_cg,
    nb_kernel131_adress_cg,
    nb_kernel132_adress_cg,
    nb_kernel133_adress_cg,
    nb_kernel134_adress_cg,
    nb_kernel200_adress_cg,
    nb_kernel201_adress_cg,
    nb_kernel202_adress_cg,
    nb_kernel203_adress_cg,
    nb_kernel204_adress_cg,
    nb_kernel210_adress_cg,
    nb_kernel211_adress_cg,
    nb_kernel212_adress_cg,
    nb_kernel213_adress_cg,
    nb_kernel214_adress_cg,
    nb_kernel220_adress_cg,
    nb_kernel221_adress_cg,
    nb_kernel222_adress_cg,
    nb_kernel223_adress_cg,
    nb_kernel224_adress_cg,
    nb_kernel230_adress_cg,
    nb_kernel231_adress_cg,
    nb_kernel232_adress_cg,
    nb_kernel233_adress_cg,
    nb_kernel234_adress_cg,
    nb_kernel300_adress_cg,
    nb_kernel301_adress_cg,
    nb_kernel302_adress_cg,
    nb_kernel303_adress_cg,
    nb_kernel304_adress_cg,
    nb_kernel310_adress_cg,
    nb_kernel311_adress_cg,
    nb_kernel312_adress_cg,
    nb_kernel313_adress_cg,
    nb_kernel314_adress_cg,
    nb_kernel320_adress_cg,
    nb_kernel321_adress_cg,
    nb_kernel322_adress_cg,
    nb_kernel323_adress_cg,
    nb_kernel324_adress_cg,
    nb_kernel330_adress_cg,
    nb_kernel331_adress_cg,
    nb_kernel332_adress_cg,
    nb_kernel333_adress_cg,
    nb_kernel334_adress_cg,
    nb_kernel400_adress_cg,
    nb_kernel410_adress_cg,
    nb_kernel430_adress_cg,
    nb_kernel010_adress_ex,
    nb_kernel020_adress_ex,
    nb_kernel030_adress_ex,
    nb_kernel100_adress_ex,
    nb_kernel101_adress_ex,
    nb_kernel102_adress_ex,
    nb_kernel103_adress_ex,
    nb_kernel104_adress_ex,
    nb_kernel110_adress_ex,
    nb_kernel111_adress_ex,
    nb_kernel112_adress_ex,
    nb_kernel113_adress_ex,
    nb_kernel114_adress_ex,
    nb_kernel120_adress_ex,
    nb_kernel121_adress_ex,
    nb_kernel122_adress_ex,
    nb_kernel123_adress_ex,
    nb_kernel124_adress_ex,
    nb_kernel130_adress_ex,
    nb_kernel131_adress_ex,
    nb_kernel132_adress_ex,
    nb_kernel133_adress_ex,
    nb_kernel134_adress_ex,
    nb_kernel200_adress_ex,
    nb_kernel201_adress_ex,
    nb_kernel202_adress_ex,
    nb_kernel203_adress_ex,
    nb_kernel204_adress_ex,
    nb_kernel210_adress_ex,
    nb_kernel211_adress_ex,
    nb_kernel212_adress_ex,
    nb_kernel213_adress_ex,
    nb_kernel214_adress_ex,
    nb_kernel220_adress_ex,
    nb_kernel221_adress_ex,
    nb_kernel222_adress_ex,
    nb_kernel223_adress_ex,
    nb_kernel224_adress_ex,
    nb_kernel230_adress_ex,
    nb_kernel231_adress_ex,
    nb_kernel232_adress_ex,
    nb_kernel233_adress_ex,
    nb_kernel234_adress_ex,
    nb_kernel300_adress_ex,
    nb_kernel301_adress_ex,
    nb_kernel302_adress_ex,
    nb_kernel303_adress_ex,
    nb_kernel304_adress_ex,
    nb_kernel310_adress_ex,
    nb_kernel311_adress_ex,
    nb_kernel312_adress_ex,
    nb_kernel313_adress_ex,
    nb_kernel314_adress_ex,
    nb_kernel320_adress_ex,
    nb_kernel321_adress_ex,
    nb_kernel322_adress_ex,
    nb_kernel323_adress_ex,
    nb_kernel324_adress_ex,
    nb_kernel330_adress_ex,
    nb_kernel331_adress_ex,
    nb_kernel332_adress_ex,
    nb_kernel333_adress_ex,
    nb_kernel334_adress_ex,
    nb_kernel400_adress_ex,
    nb_kernel410_adress_ex,
    nb_kernel430_adress_ex,
};

void
nb_kernel_setup_adress(FILE *log, nb_adress_kernel_t **list_adress)
{
  int i;
  nb_adress_kernel_t *p;

    if(NULL != log)
  fprintf(log,"AdResS simulation: Configuring adress C nonbonded kernels...\n");

  for(i=0;i<eNR_NBKERNEL_NR;i++)
  {
    p = kernellist_adress[i];
    if(p!=NULL)
      list_adress[i] = p;
  }
}
