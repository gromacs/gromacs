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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
static char *SRCID_ifunc_c = "$Id$";

#include "typedefs.h"
#include "bondf.h"
#include "disre.h"

#define def_bond(str,lstr,nra,nrpa,nrpb,ind,func) \
   {str,lstr,(nra),(nrpa),(nrpb),IF_BOND,(ind),(func)}
   
#define def_connect(str,lstr,nra,nrpa,nrpb,ind,func) \
   {str,lstr,(nra),(nrpa),(nrpb),IF_BOND | IF_CONNECT,(ind),(func)}
   
#define def_dumm(str,lstr,nra,nrpa) \
   {str,lstr,(nra),(nrpa),0,IF_DUMMY | IF_CONNECT, -1, unimplemented}
   
#define def_shk(str,lstr,nra,nrpa,nrpb,ind,func)  \
   {str,lstr,(nra),(nrpa),(nrpb),IF_SHAKE,(ind),(func)}
   
#define def_shkcon(str,lstr,nra,nrpa,nrpb,ind,func)  \
   {str,lstr,(nra),(nrpa),(nrpb),IF_SHAKE | IF_CONNECT,(ind),(func)}
   
#define def_nb(str,lstr,nra,nrp)                  \
   {str,lstr,(nra),(nrp),0,IF_NULL,-1,unimplemented}
   
#define def_nofc(str,lstr)                        \
   {str,lstr,0,0,0,IF_NULL,-1,unimplemented}

/* this MUST correspond to the enum in include/types/idef.h */
t_interaction_function interaction_function[F_NRE]=
{
  def_bond   ("ANGLES",   "Angle",          3, 2, 2,  eNR_ANGLES, angles),
  def_nb     ("BHAM",     "BuckingHam",     2, 3),
  def_connect("BONDS",    "Bonds",          2, 2, 2,  eNR_BONDS,  bonds),
  def_connect("MORSE",    "Morse",          2, 3, 0,  eNR_MORSE,  morsebonds),
  def_bond   ("WATERPOL", "Water Pol.",     1, 6, 0,  eNR_WPOL,   water_pol),
  def_bond   ("DISRES",   "Dis. Res",       2, 6, 0,  eNR_DISRES, ta_disres),
  def_bond   ("IDIHS",    "Impropers",      4, 2, 2,  eNR_IMPROPER, idihs),
  def_nb     ("LJ",       "LJ",             2, 2),
  def_bond   ("LJ14",     "Coulomb+LJ-14",  2, 2, 2,  eNR_LJC, do_14),
  def_nofc   ("LR",       "Coulomb (LR)"),
  def_nofc   ("LJLR",     "Dispersion(LR)"),
  def_bond   ("PDIHS",    "Proper Dih.",    4, 3, 3,  eNR_PROPER, pdihs),
  def_bond   ("POSRES",   "Position Rest.", 1, 3, 0,  eNR_POSRES, posres),
  def_bond   ("RBDIHS",   "Ryckaert-Bell.", 4, 6, 0,  eNR_RB, rbdihs),
  def_shkcon ("SHAKE",    "Shake",          2, 1, 1,  -1, unimplemented),
  def_shk    ("SETTLE",   "Settle",         1, 2, 0,  -1, unimplemented),
  def_dumm   ("DUMMY1",   "Dummy1",         3, 1),
  def_dumm   ("DUMMY2",   "Dummy2",         4, 2),
  def_dumm   ("DUMMY2FD", "Dummy2'",        4, 2),
  def_dumm   ("DUMMY2FAD","Dummy2''",       4, 2),
  def_dumm   ("DUMMY3",   "Dummy3",         4, 3),
  def_nofc   ("SR",       "Coulomb (SR)"),
  def_nofc   ("EPOT",     "Potential"), 
  def_nofc   ("EKIN",     "Kinetic En."),
  def_nofc   ("ETOT",     "Total Energy"),
  def_nofc   ("TEMP",     "Temperature"),
  def_nofc   ("PRES",     "Pressure"),
  def_nofc   ("DV/DL",    "d Vpot/d mu"),
  def_nofc   ("DK/DL",    "d Ekin/d mu")
};
#undef def_bond
#undef def_connect
#undef def_nb
#undef def_dumm
#undef def_shk
#undef def_shkcon
#undef def_nofc

bool have_interaction(t_idef *idef,int ftype)
{
  int i;
  
  for(i=0; (i<idef->ntypes); i++)
    if (idef->functype[i] == ftype)
      return TRUE;
  return FALSE;
}
