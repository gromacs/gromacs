/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include "types/nrnb.h"
#include "nb_kernel_f77_double.h"
#include "../nb_kerneltype.h"


/* Include standard kernel headers in local directory */
#include "nbkernel010_f77_double.h"
#include "nbkernel020_f77_double.h"
#include "nbkernel030_f77_double.h"
#include "nbkernel100_f77_double.h"
#include "nbkernel101_f77_double.h"
#include "nbkernel102_f77_double.h"
#include "nbkernel103_f77_double.h"
#include "nbkernel104_f77_double.h"
#include "nbkernel110_f77_double.h"
#include "nbkernel111_f77_double.h"
#include "nbkernel112_f77_double.h"
#include "nbkernel113_f77_double.h"
#include "nbkernel114_f77_double.h"
#include "nbkernel120_f77_double.h"
#include "nbkernel121_f77_double.h"
#include "nbkernel122_f77_double.h"
#include "nbkernel123_f77_double.h"
#include "nbkernel124_f77_double.h"
#include "nbkernel130_f77_double.h"
#include "nbkernel131_f77_double.h"
#include "nbkernel132_f77_double.h"
#include "nbkernel133_f77_double.h"
#include "nbkernel134_f77_double.h"
#include "nbkernel200_f77_double.h"
#include "nbkernel201_f77_double.h"
#include "nbkernel202_f77_double.h"
#include "nbkernel203_f77_double.h"
#include "nbkernel204_f77_double.h"
#include "nbkernel210_f77_double.h"
#include "nbkernel211_f77_double.h"
#include "nbkernel212_f77_double.h"
#include "nbkernel213_f77_double.h"
#include "nbkernel214_f77_double.h"
#include "nbkernel220_f77_double.h"
#include "nbkernel221_f77_double.h"
#include "nbkernel222_f77_double.h"
#include "nbkernel223_f77_double.h"
#include "nbkernel224_f77_double.h"
#include "nbkernel230_f77_double.h"
#include "nbkernel231_f77_double.h"
#include "nbkernel232_f77_double.h"
#include "nbkernel233_f77_double.h"
#include "nbkernel234_f77_double.h"
#include "nbkernel300_f77_double.h"
#include "nbkernel301_f77_double.h"
#include "nbkernel302_f77_double.h"
#include "nbkernel303_f77_double.h"
#include "nbkernel304_f77_double.h"
#include "nbkernel310_f77_double.h"
#include "nbkernel311_f77_double.h"
#include "nbkernel312_f77_double.h"
#include "nbkernel313_f77_double.h"
#include "nbkernel314_f77_double.h"
#include "nbkernel320_f77_double.h"
#include "nbkernel321_f77_double.h"
#include "nbkernel322_f77_double.h"
#include "nbkernel323_f77_double.h"
#include "nbkernel324_f77_double.h"
#include "nbkernel330_f77_double.h"
#include "nbkernel331_f77_double.h"
#include "nbkernel332_f77_double.h"
#include "nbkernel333_f77_double.h"
#include "nbkernel334_f77_double.h"
#include "nbkernel400_f77_double.h"
#include "nbkernel410_f77_double.h"
#include "nbkernel420_f77_double.h"
#include "nbkernel430_f77_double.h"


static nb_kernel_t *
kernellist[eNR_NBKERNEL_NR] = 
{
    nbkernel010_f77_double,
    nbkernel020_f77_double,
    nbkernel030_f77_double,
    nbkernel100_f77_double,
    nbkernel101_f77_double,
    nbkernel102_f77_double,
    nbkernel103_f77_double,
    nbkernel104_f77_double,
    nbkernel110_f77_double,
    nbkernel111_f77_double,
    nbkernel112_f77_double,
    nbkernel113_f77_double,
    nbkernel114_f77_double,
    nbkernel120_f77_double,
    nbkernel121_f77_double,
    nbkernel122_f77_double,
    nbkernel123_f77_double,
    nbkernel124_f77_double,
    nbkernel130_f77_double,
    nbkernel131_f77_double,
    nbkernel132_f77_double,
    nbkernel133_f77_double,
    nbkernel134_f77_double,
    nbkernel200_f77_double,
    nbkernel201_f77_double,
    nbkernel202_f77_double,
    nbkernel203_f77_double,
    nbkernel204_f77_double,
    nbkernel210_f77_double,
    nbkernel211_f77_double,
    nbkernel212_f77_double,
    nbkernel213_f77_double,
    nbkernel214_f77_double,
    nbkernel220_f77_double,
    nbkernel221_f77_double,
    nbkernel222_f77_double,
    nbkernel223_f77_double,
    nbkernel224_f77_double,
    nbkernel230_f77_double,
    nbkernel231_f77_double,
    nbkernel232_f77_double,
    nbkernel233_f77_double,
    nbkernel234_f77_double,
    nbkernel300_f77_double,
    nbkernel301_f77_double,
    nbkernel302_f77_double,
    nbkernel303_f77_double,
    nbkernel304_f77_double,
    nbkernel310_f77_double,
    nbkernel311_f77_double,
    nbkernel312_f77_double,
    nbkernel313_f77_double,
    nbkernel314_f77_double,
    nbkernel320_f77_double,
    nbkernel321_f77_double,
    nbkernel322_f77_double,
    nbkernel323_f77_double,
    nbkernel324_f77_double,
    nbkernel330_f77_double,
    nbkernel331_f77_double,
    nbkernel332_f77_double,
    nbkernel333_f77_double,
    nbkernel334_f77_double,
    nbkernel400_f77_double,
    nbkernel410_f77_double,
    nbkernel430_f77_double,
	nbkernel010nf_f77_double,
    nbkernel020nf_f77_double,
    nbkernel030nf_f77_double,
    nbkernel100nf_f77_double,
    nbkernel101nf_f77_double,
    nbkernel102nf_f77_double,
    nbkernel103nf_f77_double,
    nbkernel104nf_f77_double,
    nbkernel110nf_f77_double,
    nbkernel111nf_f77_double,
    nbkernel112nf_f77_double,
    nbkernel113nf_f77_double,
    nbkernel114nf_f77_double,
    nbkernel120nf_f77_double,
    nbkernel121nf_f77_double,
    nbkernel122nf_f77_double,
    nbkernel123nf_f77_double,
    nbkernel124nf_f77_double,
    nbkernel130nf_f77_double,
    nbkernel131nf_f77_double,
    nbkernel132nf_f77_double,
    nbkernel133nf_f77_double,
    nbkernel134nf_f77_double,
    nbkernel200nf_f77_double,
    nbkernel201nf_f77_double,
    nbkernel202nf_f77_double,
    nbkernel203nf_f77_double,
    nbkernel204nf_f77_double,
    nbkernel210nf_f77_double,
    nbkernel211nf_f77_double,
    nbkernel212nf_f77_double,
    nbkernel213nf_f77_double,
    nbkernel214nf_f77_double,
    nbkernel220nf_f77_double,
    nbkernel221nf_f77_double,
    nbkernel222nf_f77_double,
    nbkernel223nf_f77_double,
    nbkernel224nf_f77_double,
    nbkernel230nf_f77_double,
    nbkernel231nf_f77_double,
    nbkernel232nf_f77_double,
    nbkernel233nf_f77_double,
    nbkernel234nf_f77_double,
    nbkernel300nf_f77_double,
    nbkernel301nf_f77_double,
    nbkernel302nf_f77_double,
    nbkernel303nf_f77_double,
    nbkernel304nf_f77_double,
    nbkernel310nf_f77_double,
    nbkernel311nf_f77_double,
    nbkernel312nf_f77_double,
    nbkernel313nf_f77_double,
    nbkernel314nf_f77_double,
    nbkernel320nf_f77_double,
    nbkernel321nf_f77_double,
    nbkernel322nf_f77_double,
    nbkernel323nf_f77_double,
    nbkernel324nf_f77_double,
    nbkernel330nf_f77_double,
    nbkernel331nf_f77_double,
    nbkernel332nf_f77_double,
    nbkernel333nf_f77_double,
    nbkernel334nf_f77_double,
    nbkernel400nf_f77_double,
    nbkernel410nf_f77_double,
    nbkernel430nf_f77_double,
};


void
nb_kernel_setup_f77_double(FILE *log,nb_kernel_t **list)
{
  int i;
  nb_kernel_t *p;
  
  if(log)
    fprintf(log,"Configuring double precision Fortran kernels...\n");
  
  for(i=0;i<eNR_NBKERNEL_NR;i++)
  {
      p = kernellist[i];
      if(p!=NULL)
	list[i] = p; 
  }
}    
