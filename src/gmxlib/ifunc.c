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
static char *SRCID_ifunc_c = "$Id$";
#include "typedefs.h"
#include "bondf.h"
#include "disre.h"
#include "orires.h"

#define def_bonded(str,lstr,nra,nrpa,nrpb,ind,func)\
   {str,lstr,(nra),(nrpa),(nrpb),IF_BOND,                        (ind),(func)}
   
#define  def_angle(str,lstr,nra,nrpa,nrpb,ind,func)\
   {str,lstr,(nra),(nrpa),(nrpb),IF_BOND | IF_ATYPE,(ind),(func)}
   
#define   def_bond(str,lstr,nra,nrpa,nrpb,ind,func)\
   {str,lstr,(nra),(nrpa),(nrpb),IF_BOND | IF_CONNECT | IF_BTYPE,(ind),(func)}

#define def_bondnb(str,lstr,nra,nrpa,nrpb,ind,func)\
   {str,lstr,(nra),(nrpa),(nrpb),IF_BOND | IF_CONNECT,(ind),(func)}

#define  def_dummy(str,lstr,nra,nrpa)\
   {str,lstr,(nra),(nrpa),     0,IF_DUMMY | IF_CONNECT,     -1, unimplemented}
   
#define    def_shk(str,lstr,nra,nrpa,nrpb)\
   {str,lstr,(nra),(nrpa),(nrpb),IF_CONSTRAINT,             -1, unimplemented}

#define def_shkcon(str,lstr,nra,nrpa,nrpb)\
   {str,lstr,(nra),(nrpa),(nrpb),IF_CONSTRAINT | IF_CONNECT,-1, unimplemented}
   
#define     def_nb(str,lstr,nra, nrp)\
   {str,lstr,(nra), (nrp),     0,IF_NULL,                    -1,unimplemented}
   
#define   def_nofc(str,lstr)\
   {str,lstr,    0,     0,     0,IF_NULL,                    -1,unimplemented}

/* this MUST correspond to the enum in include/types/idef.h */
t_interaction_function interaction_function[F_NRE]=
{
  def_bond   ("BONDS",    "Bond",            2, 2, 2,  eNR_BONDS,  bonds    ),
  def_bond   ("G96BONDS", "G96Bond",         2, 2, 2,  eNR_BONDS,  g96bonds ),
  def_bond   ("MORSE",    "Morse",           2, 3, 0,  eNR_MORSE, morsebonds),
  def_bond   ("CUBICBONDS","Cubic Bonds",    2, 3, 0,  eNR_CUBICBONDS, cubicbonds),
  def_bondnb ("CONNBONDS","Connect Bonds",   2, 0, 0,  0,      unimplemented),
  def_bonded ("HARMONIC", "Harmonic Pot.",   2, 2, 2,  eNR_BONDS,  bonds    ),
  def_angle  ("ANGLES",   "Angle",           3, 2, 2,  eNR_ANGLES, angles   ),
  def_angle  ("G96ANGLES","G96Angle",        3, 2, 2,  eNR_ANGLES, g96angles),
  def_bonded ("PDIHS",    "Proper Dih.",     4, 3, 3,  eNR_PROPER, pdihs    ),
  def_bonded ("RBDIHS",   "Ryckaert-Bell.",  4, 6, 0,  eNR_RB, rbdihs       ),
  def_bonded ("IDIHS",    "Improper Dih.",   4, 2, 2,  eNR_IMPROPER,idihs   ),
  def_bonded ("LJ14",     "LJ-14",           2, 2, 2,  eNR_INL1100, do_14   ),
  def_nofc   ("COUL14",   "Coulomb-14"       ),
  def_nb     ("LJ",       "LJ (SR)",         2, 2      ),
  def_nb     ("BHAM",     "BuckingHam",      2, 3      ),
  def_nofc   ("LJLR",     "LJ (LR)"          ),
  def_nofc   ("DISPCORR", "Disper. corr."    ),
  def_nofc   ("SR",       "Coulomb (SR)"     ),
  def_nofc   ("LR",       "Coulomb (LR)"     ),
  def_bonded ("WATERPOL", "Water Pol.",      1, 6, 0,  eNR_WPOL,   water_pol),
  def_bonded ("POSRES",   "Position Rest.",  1, 3, 0,  eNR_POSRES, posres   ),
  def_bonded ("DISRES",   "Dis. Rest.",      2, 6, 0,  eNR_DISRES, ta_disres),
  def_nofc   ("DRVIOL",   "D. R. Viol. (nm)" ),    
  def_bonded ("ORIRES",   "Orient. Rest.",   2, 6, 0,  eNR_ORIRES, orires   ),
  def_nofc   ("ORDEV",    "Ori. R. RMSD"     ),  
  def_bonded ("ANGRES",   "Angle Rest.",     4, 3, 3,  eNR_ANGRES, angres   ),
  def_bonded ("ANGRESZ",  "Angle Rest. Z",   2, 3, 3,  eNR_ANGRESZ,angresz  ),
  def_shkcon ("CONSTR",   "Constraint",      2, 1, 1   ),
  def_shk    ("CONSTRNC", "Constr. No Conn.",2, 1, 1   ),
  def_shk    ("SETTLE",   "Settle",          1, 2, 0   ),
  def_dummy  ("DUMMY2",   "Dummy2",          3, 1      ),
  def_dummy  ("DUMMY3",   "Dummy3",          4, 2      ),
  def_dummy  ("DUMMY3FD", "Dummy3fd",        4, 2      ),
  def_dummy  ("DUMMY3FAD","Dummy3fad",       4, 2      ),
  def_dummy  ("DUMMY3OUT","Dummy3out",       4, 3      ),
  def_dummy  ("DUMMY4FD", "Dummy4fd",        5, 3      ),
  def_nofc   ("EQM",      "Quantum En."      ),
  def_nofc   ("EPOT",     "Potential"        ),
  def_nofc   ("EKIN",     "Kinetic En."      ),
  def_nofc   ("ETOT",     "Total Energy"     ),
  def_nofc   ("TEMP",     "Temperature"      ),
  def_nofc   ("PRES",     "Pressure (bar)"   ),
  def_nofc   ("DV/DL",    "dVpot/dlambda"    ),
  def_nofc   ("DK/DL",    "dEkin/dlambda"    )
};

bool have_interaction(t_idef *idef,int ftype)
{
  int i;
  
  for(i=0; (i<idef->ntypes); i++)
    if (idef->functype[i] == ftype)
      return TRUE;
  return FALSE;
}
