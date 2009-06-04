/*
 *
 * Gromacs 4.0                         Copyright (c) 1991-2003
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

/* This file is NOT threadsafe, but it is only used to create
 * the nonbonded functions during the build process, so it will never be
 * executed by multiple threads.
 */

#include <mknb_common.h>
#include <mknb_metacode.h>

#include <string.h>
#include <stdlib.h>


/* fscal terms to be mult. by -tabscale/r */
static char 
fs_minus_tabscale_rinv[1024]; 


/* fscal terms to be mult. by -1/r */
static char 
fs_minus_rinv[1024];         


/* fscal terms to be mult. by 1/(r*r) */
static char 
fs_rinvsq[1024];             




/* UTILITY ROUTINES FOR TABLE INTERACTIONS: */

/* Calculate the table index from rsq and rinv */
int
mknb_calc_table_index(char *r)
{
  int nflops = 0;
  
  mknb_comment("Calculate table index");
  mknb_assign("rt","r*tabscale");
  nflops++;
  
  /* Truncate rt to an integer. */
  mknb_assign("n0","rt");

  mknb_assign("eps","rt-n0");
  mknb_assign("eps2","eps*eps");
  nflops += 2;

  mknb_assign("nnn","%d*n0%s",mknb_func.table_element_size, 
	      (mknb_fortran) ? "+1" : "");

  return nflops;

}

/* Calculate the table index from rsq and rinv */
int
mknb_calc_gbtable_index(char *r)
{
  int nflops = 0;
  
  mknb_comment("Calculate table index");
  mknb_assign("rt","r*gbscale");
  nflops ++;
  
  /* Truncate rt to an integer. */
  mknb_assign("n0","rt");
  
  mknb_assign("eps","rt-n0");
  mknb_assign("eps2","eps*eps");
  nflops += 2;

  mknb_assign("nnn","4*n0%s", (mknb_fortran) ? "+1" : "");
  return nflops;

}


/* Perform a table lookup and calculate VV and FF */
int
mknb_read_table(char *tabname)
{
  int nflops = 0;

  /* See the Gromacs manual for details on cubic spline table interpolation */
  mknb_assign("Y",mknb_array(tabname,"nnn"));
  mknb_assign("F",mknb_array(tabname,"nnn+1"));
  mknb_assign("Geps","eps*%s",mknb_array(tabname,"nnn+2"));
  mknb_assign("Heps2","eps2*%s",mknb_array(tabname,"nnn+3"));
  
  mknb_assign("Fp","F+Geps+Heps2");
  mknb_assign("VV","Y+eps*Fp");
  nflops += 6;

  if(mknb_func.do_force) {
    mknb_assign("FF","Fp+Geps+2.0*Heps2");
    nflops += 3;
  }
  return nflops;
}

/* COULOMB INTERACTIONS */

int
mknb_coul_normal(char *rinv)
{
  mknb_comment("Coulomb interaction");
  /* The Coulomb potential is simply charge/r.
   * vcoul is just a temporary variable defined in mknb_declarations.c */
  mknb_assign("vcoul","qq*%s",rinv);
  /* The Coulomb force is -charge/(r*r) = -vcoul/r.
   * We are going to multiply by the vector r later, and take care of
   * the sign when incrementing/decrementing forces, so the part we
   * want to add to fs is  vcoul/(r*r).
   */
  if(mknb_func.do_force)
    sprintf(fs_rinvsq,"vcoul");

  /* Update total Coulomb energy */
  mknb_assign("vctot","vctot+vcoul");
  /* Done. 2 flops */
  return 2;
}

int
mknb_coul_rf(char *rsq,char *rinv)
{
  mknb_comment("Coulomb reaction-field interaction");

  mknb_assign("krsq","krf*%s",rsq);
  mknb_assign("vcoul","qq*(%s+krsq-crf)",rinv);

  if(mknb_func.do_force)
    sprintf(fs_rinvsq,"qq*(%s-2.0*krsq)",rinv);

  mknb_assign("vctot","vctot+vcoul");
  /* Done. 8 flops with force, 5 for energy only */
  return mknb_func.do_force ? 8 : 5;
}

int
mknb_coul_tab(char *rsq, char *rinv)
{
  int nflops = 0;
  
  mknb_comment("Tabulated coulomb interaction");

  /* We already calculated the table index. Coulomb is always first in table,
   * so we just use it and do the lookup.
   */
  nflops += mknb_read_table("VFtab");
  mknb_assign("vcoul","qq*VV");
  nflops++;
  if(mknb_func.do_force) {
    mknb_assign("fijC","qq*FF");
    nflops++;
    /* fs_minusrinv is empty */
    sprintf(fs_minus_tabscale_rinv,"fijC");
  }
  mknb_assign("vctot","vctot + vcoul");
  nflops++;
  
  return nflops;
}
  

int
mknb_coul_gb(char *rsq, char *rinv)
{
  int nflops = 0;
  
  mknb_comment("Tabulated Generalized-Born interaction");
  if(mknb_func.do_force)
    mknb_assign("dvdaj",mknb_array("dvda","jnr"));
  /* Coulomb is always the first interaction, so we
   * have not calculated the table index yet.
   */
  mknb_assign("r","%s*%s",rsq,rinv);
  nflops++;
  nflops += mknb_calc_gbtable_index("r");
 
  nflops += mknb_read_table("GBtab");
  mknb_assign("vgb","qq*VV");
  nflops++;
  if(mknb_func.do_force) {
    mknb_assign("fijC","qq*FF*gbscale");
    mknb_assign("dvdatmp","-0.5*(vgb+fijC*r)");
    mknb_assign("dvdasum","dvdasum + dvdatmp");
    nflops+=6;
    /* fs_minusrinv is empty */
    sprintf(fs_minus_rinv,"fijC-fscal");
    /* Update j atom dvda */
    mknb_assign(mknb_array("dvda","jnr"),"dvdaj+dvdatmp*isaj*isaj");
  }
  /* This will only give thw Coulomb part back to the total potential */
  mknb_assign("vctot","vctot + vcoul"); 
  nflops++;
  
  return nflops;
}
  
/* MKNB_VDW INTERACTIONS */
int
mknb_vdw_lj(char *rsq, char *rinv)
{
  mknb_comment("Lennard-Jones interaction");

  mknb_assign("rinvsix","rinvsq*rinvsq*rinvsq");
  mknb_assign("Vvdw6","c6*rinvsix");
  mknb_assign("Vvdw12","c12*rinvsix*rinvsix");

  if(mknb_func.do_force) {
    if(strlen(fs_rinvsq)==0)
      sprintf(fs_rinvsq,"12.0*Vvdw12-6.0*Vvdw6");
    else /* append */
      strcat(fs_rinvsq,"+12.0*Vvdw12-6.0*Vvdw6");
  }

  mknb_assign("Vvdwtot","Vvdwtot+Vvdw12-Vvdw6");
  /* Done. 11 flops with force, 7 for energy only */
  return mknb_func.do_force ? 11 : 7;
}


int
mknb_vdw_bham(char *rsq, char *rinv)
{
  mknb_comment("Buckingham interaction");

  mknb_assign("rinvsix","rinvsq*rinvsq*rinvsq");
  mknb_assign("Vvdw6","c6*rinvsix");
  mknb_assign("br","cexp2*%s*%s",rsq,rinv); /* br=cexp2*r */
  mknb_assign("Vvdwexp","cexp1*exp(-br)");

  if(mknb_func.do_force) {
    if(strlen(fs_rinvsq)==0)
      sprintf(fs_rinvsq,"br*Vvdwexp-6.0*Vvdw6");
    else /* append */
      strcat(fs_rinvsq,"+br*Vvdwexp-6.0*Vvdw6");
  }

  mknb_assign("Vvdwtot","Vvdwtot+Vvdwexp-Vvdw6");
  /* exp() is expensive, 25 flops is a low estimate.
   * This gives about 37 flops , 34 for energy only
   */      
  return mknb_func.do_force ? 37 : 34;
}


int
mknb_vdw_tab(char *rsq, char *rinv)
{
    int nflops = 0;
    
    mknb_comment("Tabulated VdW interaction - dispersion");
    
    /* If we just did tabulated coulomb we already
     * have the lookup seed. Just add the offset (4)
     * to use the LJ instead of the coulomb table.
     */
    if(mknb_func.table_element_size==12) 
    {
        mknb_assign("nnn","nnn+4");
    }
    /* Without coulomb, dispersion is first element - nothing to do */
   
    nflops += mknb_read_table("VFtab");
    mknb_assign("Vvdw6","c6*VV");
    
    nflops++;
    if(mknb_func.do_force) 
    {
        mknb_assign("fijD","c6*FF");
        nflops++;
        if(strlen(fs_minus_tabscale_rinv)==0)
            sprintf(fs_minus_tabscale_rinv,"fijD");
        else 
        {
            strcat(fs_minus_tabscale_rinv,"+fijD");
            nflops++;
        }
    }

    mknb_comment("Tabulated VdW interaction - repulsion");
    mknb_assign("nnn", "nnn+4");
    nflops += mknb_read_table("VFtab");
    mknb_assign("Vvdw12","c12*VV");
    
    nflops++;
    if(mknb_func.do_force) 
    {
        mknb_assign("fijR","c12*FF");
        nflops++;
        /* string is never empty here */
        strcat(fs_minus_tabscale_rinv,"+fijR");
        nflops++;
    }
    mknb_assign("Vvdwtot","Vvdwtot+ Vvdw6 + Vvdw12");
    nflops += 2;
    
    return nflops;
}


int
mknb_calculate_interaction(char *rsq, char *rinv)
{
  int nflops = 0;
  char tmp[512];
  
  /* Time to do the actual work. All the extra code in the function
   * generation program only comes down to one thing: to make it
   * possible to execute these switch statements at compile time
   * instead of when running the code (as we do in the unoptimized
   * routine for nonbonded forces).
   */

  /* Available variables:
   * rsq: a textstring with the variable holding r*r
   * rinv: textstring with the variable for 1/r
   *
   * The only exception is LJ-only functions, where we only need 1/(r*r)
   * and calculated it directly as 1/(rsq). This value can
   * be accessed directly as "rinvsq".
   */

  /* Add your potential to vctot or Vvdwtot.
   * The force is a bit special:
   *
   * First - skip all force calculations when mknb_func.do_force is false!
   *
   * To calculate the vectorial force we need to divide by the scalar r
   * and then multiply by the vector r. To avoid a stupid division, we
   * calculate (f/r) directly and add it to the variable "fs" with the
   * same sign as the potential (i.e. fs=grad(V)*r).
   *
   * To avoid duplicate calculations, we use two intermediary strings:
   * "fs_rinv" and "fs_rinvsq".
   *
   * The contents of these strings will be
   * multiplied by by rinv and rinvsq respectively, and added to get fs.
   *
   * So, if your coulomb interaction writes "A" to fs_rinv, "B" to
   * fs_rinvsq, and the VdW routine adds "-C" to the end of fs_rinv,
   * the fs variable will be
   *
   * fs = (A-C)*rinv + (B)*rinvsq
   */

  /* A lot of routines need 1/(r*r). To avoid duplicating operations
   * we calculate it before calling the interactions, when necessary.
   *
   * For LJ-only calculation it has already been computed.
   * All coulomb calculations except tabulated ones need rinvsq for force.
   * All vdw calculations except tabulated ones need it (even for energy).
   */
  if(((mknb_func.coul==MKNB_COUL_NORMAL || mknb_func.coul==MKNB_COUL_RF) && mknb_func.do_force) ||
     (mknb_func.vdw==MKNB_VDW_LJ && mknb_func.coul) || mknb_func.vdw==MKNB_VDW_BHAM) {
    mknb_assign("rinvsq","%s*%s",rinv,rinv);
    nflops++;
  }
  
  strcpy(fs_minus_tabscale_rinv,"");
  strcpy(fs_minus_rinv,"");
  strcpy(fs_rinvsq,"");


  if(mknb_func.coul==MKNB_COUL_TAB)
  {
      mknb_comment("Calculate table index");
      mknb_assign("r","%s*%s",rsq,rinv);
      nflops++;
      nflops += mknb_calc_table_index("r");
  }
  
  switch(mknb_func.coul) {
  case MKNB_COUL_NO:
    break;
  case MKNB_COUL_NORMAL:
    nflops += mknb_coul_normal(rinv);
    break;
  case MKNB_COUL_RF:
    nflops += mknb_coul_rf(rsq,rinv);
    break;
  case MKNB_COUL_TAB:
    nflops += mknb_coul_tab(rsq,rinv);
    break;
  case MKNB_COUL_GB:
    nflops += mknb_coul_gb(rsq,rinv);
    break;
  default:
    fprintf(stderr,"Error: Coulomb type %d undefined (mknb_interactions.c)\n",mknb_func.coul);
    exit(0);    
  }

  if(mknb_func.vdw==MKNB_VDW_TAB && mknb_func.coul!=MKNB_COUL_TAB)
  {
      mknb_comment("Calculate table index");
      mknb_assign("r","%s*%s",rsq,rinv);
      nflops++;
      nflops += mknb_calc_table_index("r");
  }
  
  switch(mknb_func.vdw) {
  case MKNB_VDW_NO:
    break;
  case MKNB_VDW_LJ:
    nflops += mknb_vdw_lj(rsq,rinv);
    break;
  case MKNB_VDW_BHAM:
    nflops += mknb_vdw_bham(rsq,rinv);
    break;
  case MKNB_VDW_TAB:
    nflops += mknb_vdw_tab(rsq,rinv);
    break;
  default:
    fprintf(stderr,"Error: VdW type %d undefined (check mknb_interactions.c)",mknb_func.vdw);
  }

  /* Assembly the scalar force (actually f/dr, since we will multiply with vectorial dr) */
  if(strlen(fs_minus_tabscale_rinv)>0) {
    if(strlen(fs_minus_rinv)>0) {
      sprintf(tmp,"(%s)*tabscale+%s",fs_minus_tabscale_rinv,fs_minus_rinv);
      nflops+=2;
    } else {
      sprintf(tmp,"(%s)*tabscale",fs_minus_tabscale_rinv);
      nflops++;
    }
    strcpy(fs_minus_rinv,tmp);
  }
  
  if(strlen(fs_rinvsq)>0 && strlen(fs_minus_rinv)>0) {
	  mknb_assign("fscal","(%s)*rinvsq-(%s)*%s",fs_rinvsq,fs_minus_rinv,rinv);
	  nflops += 3;
  } else if(strlen(fs_rinvsq)>0) {
	  mknb_assign("fscal","(%s)*rinvsq",fs_rinvsq);
	  nflops++;
  } else if(strlen(fs_minus_rinv)>0) {
	  mknb_assign("fscal","-(%s)*%s",fs_minus_rinv,rinv);
	  nflops += 2;
  }
	
  return nflops;
}

