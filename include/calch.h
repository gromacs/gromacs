/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROningen MAchine for Chemical Simulation
 */

#ifndef _calch_h
#define _calch_h

static char *SRCID_calch_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) calch.h 1.8 2/2/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
	
extern void calc_h_pos(int nht,int nh[],int na[],rvec x[]);
/*
 *    w.f. van gunsteren, groningen, july 1981 
 *
 *    translated to c d. van der spoel groningen jun 1993
 *    added option 5 jan 95
 *
 *    subroutine genh (nht,nh,na,d,alfa,x)                          
 *                                                                  
 *    genh generates cartesian coordinates for hydrogen atoms   
 *    using the coordinates of neighbour atoms.                     
 *                                                                  
 *    nht = 1 : one hydrogen atom (n) is generated, lying in the plane 
 *              of atoms (i,j,k) on the line bisecting angle (j-i-k)
 *              at a distance d from atom i, such that the angles   
 *              (n-i-j) and (n-i-k) are > 90 degrees                
 *        = 2 : one hydrogen atom (n) is generated at a distance d  
 *              from atom i, such that angle (n-i-j)=alfa and dihedral 
 *              (n-i-j-k)=trans                                     
 *        = 3 : two hydrogens (n1,n2) are generated at a distance d 
 *              from atom i, such that angle (n1-i-j)=(n2-i-j)=alfa 
 *              and dihedral (n1-i-j-k)=trans and (n2-i-j-k)=cis    
 *        = 4 : three (n1,n2,n3) or two (n1,n2) hydrogens are generated
 *              at a distance d from atom i, such that angle (n1-i-j)= 
 *              (n2-i-j)=(n3-i-j)=alfa, dihedral (n1-i-j-k)=trans,  
 *              (n2-i-j-k)=trans+120 and (n3-i-j-k)=trans+240 degrees  
 *        = 5 : one hydrogen is generated connected to n1, such that it is
 *              tetrahedral configuration with n2,n3 and n4
 *    nh(1.. ) = sequence numbers of the hydrogen atoms that are to be 
 *               generated (see x)                                  
 *               if nht=4 and nh(3)=0, only two hydrogens are generated
 *    na(1..4) = sequence numbers of the atoms i, j and k and l
 *    x(1.. ) = atom cartesian coordinates                          
 *    default bond lengths and angles are defined internally
 */

#endif	/* _calch_h */
