/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_testtab_c = "$Id$";
#include <stdio.h>
#include "typedefs.h"
#include "force.h"
#include "shift_util.h"

int main(int argc,char *argv[])
{
  t_forcerec *fr;
  rvec box;
  
  fr=mk_forcerec();
  fr->r1 = 0.6;
  fr->rc = 0.9;
  fr->eeltype = eelTWIN;
  box[XX]=box[YY]=box[ZZ]=1.0;
  
  set_shift_consts(stdout,fr->r1,fr->rc,box,fr);

  make_tables(fr);
  
  return 0;
}
