/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_pullio_c = "$Id$";

#include <string.h>
#include <stdlib.h>
#include "sysstuff.h"
#include "princ.h"
#include "futil.h"
#include "statutil.h"
#include "vec.h"
#include "smalloc.h"
#include "typedefs.h"
#include "names.h"
#include "fatal.h"
#include "macros.h"
#include "rdgroup.h"
#include "symtab.h"
#include "index.h"
#include "confio.h"
#include "pull.h"
#include "string.h"

void dump_conf(t_pull *pull,rvec x[],matrix box,t_topology *top, 
	       int nout, real time) 
{
  FILE *confout;
  char buf[128];
  rvec tmp,tmp1,tmp2;

  sprintf(buf,"out_%d.gro",nout);
  nout++;
  confout = ffopen(buf,"w");

  /* calculate the current positions of the center of mass of the grps 
     printed is pull - reference, so position with respect to reference
     group 
   */
  if (pull->pull.n == 2) {
    rvec_sub(pull->pull.x_unc[0],pull->ref.x_unc[0],tmp1);
    rvec_sub(pull->pull.x_unc[1],pull->ref.x_unc[0],tmp2);
    sprintf(buf,"grp1:%8.3f%8.3f%8.3f grp2:%8.3f%8.3f%8.3f t:%8.3f",
	    tmp1[0],tmp1[1],tmp1[2],tmp2[0],tmp2[1],tmp2[2],time);
  } else {
    rvec_sub(pull->pull.x_unc[0],pull->ref.x_unc[0],tmp1);
    sprintf(buf,"grp1:%8.3f%8.3f%8.3f t:%8.3f",
	    tmp1[XX],tmp1[YY],tmp1[ZZ],time);
  }
  write_hconf(confout,buf,&top->atoms,x,NULL,box);  
}

void print_start(t_pull *pull, int step) 
{
  int i;
  
  for (i=0;i<pull->pull.n;i++)
    fprintf(pull->out,"%d:%d x:%8.3f%8.3f%8.3f\n",
	    step,i,pull->pull.x_unc[i][XX],
	    pull->pull.x_unc[i][YY],pull->pull.x_unc[i][ZZ]);
}

void print_afm(t_pull *pull, int step)
{
  int i;

  for (i=0;i<pull->pull.n;i++)
    if (pull->bVerbose) 
      fprintf(pull->out,
	      "%d:%d: f:%8.3f s:%8.3f\n",step,i,pull->pull.f[i][ZZ],
	      pull->pull.spring[i][ZZ]);
    else
      fprintf(pull->out,"%8.3f\n",pull->pull.f[i][ZZ]);
}

void print_constraint(t_pull *pull, rvec *f, int step, matrix box) 
{
  int i,ii,m; 
  rvec tmp,tmp2,tmp3;

  for (i=0;i<pull->pull.n;i++) {
    if (pull->bCyl) 
      rvec_sub(pull->pull.x_con[i],pull->dyna.x_con[i],tmp);
    else 
      rvec_sub(pull->pull.x_con[i],pull->ref.x_con[0],tmp);
    for (m=0;m<DIM;m++) {
      tmp[m] *= pull->dims[m];
      if (tmp[m] < -0.5*box[m][m]) tmp[m] += box[m][m];
      if (tmp[m] > 0.5*box[m][m])  tmp[m] -= box[m][m];
    }
    if (pull->bVerbose) 
      fprintf(pull->out,"%d:%d ds:%8.3f f: %8.3f\n", step,i,norm(tmp),
	      pull->pull.f[i][ZZ]);
    else
      fprintf(pull->out,"%8.3f ",pull->pull.f[i][ZZ]);
  }
  
  fprintf(pull->out,"\n");

  /* DEBUG */ /* this code doesn't correct for pbc, needs improvement */
  if (pull->bVerbose) {
    for (i=0;i<pull->pull.n;i++) {
      if (pull->bCyl) 
	fprintf(pull->out,"eConstraint: step %d. Refgroup = dynamic (%f,%f\n"
		"Group %d (%s): ref. dist = %8.3f, unconstr. dist = %8.3f"
		" con. dist = %8.3f f_i = %8.3f\n", step, pull->r,pull->rc,
		i,pull->pull.grps[i],
		pull->dyna.x_ref[i][ZZ]-pull->pull.x_ref[i][ZZ],
		pull->dyna.x_unc[i][ZZ]-pull->pull.x_unc[i][ZZ],
		pull->dyna.x_con[i][ZZ]-pull->pull.x_con[i][ZZ],
		pull->pull.f[i][ZZ]);
      else {
	rvec_sub(pull->ref.x_ref[0],pull->pull.x_ref[i],tmp);
	rvec_sub(pull->ref.x_unc[0],pull->pull.x_unc[i],tmp2);
	rvec_sub(pull->ref.x_con[0],pull->pull.x_con[i],tmp3);
	fprintf(stderr,"grp %d:ref (%8.3f,%8.3f,%8.3f) unc(%8.3f%8.3f%8.3f\n"
		"con (%8.3f%8.3f%8.3f)\n",i, tmp[0],tmp[1],tmp[2],
		tmp2[0],tmp2[1],tmp2[2],tmp3[0],tmp3[1],tmp3[2]);
      }
    }
  } /* END DEBUG */

}

void print_umbrella(t_pull *pull, int step) 
{
  int i;
  for (i=0;i<pull->pull.n;i++)
    fprintf(pull->out,"=%d= %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ", step,
	    pull->pull.x_con[i][XX],pull->pull.x_con[i][YY],
	    pull->pull.x_con[i][ZZ],
	    pull->pull.f[i][XX],pull->pull.f[i][YY],
	    pull->pull.f[i][ZZ]);
  fprintf(pull->out,"\n");
}

void read_pullparams(t_pull *pull, char *infile, char *outfile) 
{
  t_inpfile *inp;
  int ninp,i;
  char *tmp;           /* for the input parsing macros */
  char dummy[STRLEN];  /* idem */
  char grp1buf[STRLEN], grp2buf[STRLEN], grp3buf[STRLEN], 
       bf1[STRLEN], bf2[STRLEN], dir[STRLEN], 
       refdir1[STRLEN],refdir2[STRLEN];

  int bReverse; int tmpref; int tmprun; 

  enum {erunSTART, erunAFM, erunConstraint, erunUmbrella, erunTest, erunNR};
  static char *runtypes[erunNR+1] = { 
    "start", "afm", "constraint", "umbrella", "test", NULL
  };
  enum {erefCom, erefComT0, erefDyn, erefDynT0, erefNR};
  static char *reftypes[erefNR+1] = {
    "com", "com_t0", "dynamic", "dynamic_t0", NULL
  };
  enum {ereverseTO_REF, ereverseFROM_REF, ereverseNR};
  static char *reversetypes[ereverseNR+1] = {
    "from_reference", "to_reference", NULL
  };
  enum {everboseYES, everboseNO, everboseNR};
  static char *verbosetypes[erefNR+1] = {
    "no", "yes", NULL
  };
  int nerror = 0;

  /* read input parameter file */
  fprintf(stderr,"Reading parameter file %s\n",infile);
  inp=read_inpfile(infile,&ninp);

  /* general options */
  CTYPE("GENERAL");
  EETYPE("verbose",         pull->bVerbose, verbosetypes, &nerror, TRUE);
  CTYPE("Runtype: start, afm, constraint, umbrella, test");
  EETYPE("runtype",         tmprun, runtypes, &nerror, TRUE);
  CTYPE("Groups to be pulled");
  STYPE("group_1",          grp1buf, "");
  STYPE("group_2",          grp2buf, "");
  CTYPE("The group for the reaction force.");
  STYPE("reference_group",  grp3buf, "");
  CTYPE("Ref. type: com, com_t0, dynamic, dynamic_t0");
  EETYPE("reftype",         tmpref, reftypes, &nerror, TRUE);
  CTYPE("Use running average for reflag steps for com calculation");
  ITYPE("reflag",           pull->reflag, 1);
  CTYPE("Select components for constraint vector. default: z-only");
  STYPE("direction",        dir, "0.0 0.0 1.0");
  CTYPE("Direction for start/afm: to_reference, from_reference");
  EETYPE("reverse",          pull->bReverse, reversetypes, &nerror, TRUE);

  /* options for dynamic reference groups */
  CTYPE("DYNAMIC REFERENCE GROUP OPTIONS");
  CTYPE("Cylinder radius for dynamic reaction force groups (nm)");
  RTYPE("r",                pull->r, 0.0);
  CTYPE("Switch from r to rc in case of dynamic reaction force");
  RTYPE("rc",   pull->rc,   0.0);
  CTYPE("Update frequency for dynamic reference groups (steps)");
  ITYPE("update",           pull->update, 1);

  /* options for AFM type pulling simulations */
  CCTYPE("AFM OPTIONS");
  CTYPE("pull rate in nm/timestep");
  RTYPE("pullrate",         pull->rate,    0.0);
  CTYPE("forceconstant in kJ/(mol*nm^2)");
  RTYPE("forceconstant",    pull->k, 0.0);

  /* umbrella sampling options */
  CCTYPE("UMBRELLA SAMPLING OPTIONS");
  CTYPE("Width of umbrella sampling potential in kJ/(mol*nm^2)");
  RTYPE("width",            pull->um_width, 0.0);

  /* options for making starting structures */
  CCTYPE("STARTING STRUCTURE OPTIONS");
  CTYPE("Start coord. for making starting struct, rel. to ref. grp.: x y z");
  STYPE("r0_group1",        bf1, "");
  STYPE("r0_group2",        bf2, ""); 
  /*  CTYPE("Constrain rotations around principle axes? Needs work");
      ITYPE("rotation_x",       pull->bRot[0], 0);
      ITYPE("rotation_y",       pull->bRot[1], 0);
      ITYPE("rotation_z",       pull->bRot[2], 0);
      CTYPE("Rate of rotation (degrees/step)");
      RTYPE("rotation_rate",    pull->rot_rate, 0.0);
  */
  RTYPE("tolerance",        pull->tolerance, 0.05);
  CTYPE("Rate of translation in all directions (nm/step)");
  RTYPE("translation_rate", pull->xlt_rate, 0.0); 
  CTYPE("Write out structure every ndegr degrees, transstep nm");
  ITYPE("ndegr",            pull->rot_incr, 0);
  RTYPE("transstep",        pull->xlt_incr, 0.001);

  write_inpfile(outfile,ninp,inp);
  for (i=0; (i<ninp); i++) {
    sfree(inp[i].name);
    sfree(inp[i].value);
  }
  sfree(inp);

  pull->runtype = (t_runtype)tmprun;
  pull->reftype = (t_reftype)tmpref;

  /* sort out the groups */
  fprintf(stderr,"Groups: %s %s %s\n",grp1buf,grp2buf,grp3buf);
  if (!strcmp(grp1buf,"") || !strcmp(grp3buf,"")) 
    fatal_error(0,"Need to specify at least group_1 and reference_group");
  pull->pull.n = !strcmp(grp2buf,"") ? 1 : 2;

  /* initialize the names of the groups */
  snew(pull->pull.grps,pull->pull.n);
  snew(pull->ref.grps,1);
  pull->pull.grps[0] = (char *)strdup(grp1buf);
  if (pull->pull.n == 2)
    pull->pull.grps[1] = (char *)strdup(grp2buf);
  pull->ref.grps[0]  = (char *)strdup(grp3buf);

  if (pull->runtype == eStart) {
    snew(pull->pull.xtarget,pull->pull.n);
    string2rvec(bf1,pull->pull.xtarget[0]);
    if (pull->pull.n == 2)
      string2rvec(bf2,pull->pull.xtarget[1]);
  }

  string2rvec(dir,pull->dims);
  fprintf(stderr,"Using distance components %2.1f %2.1f %2.1f\n",
	  pull->dims[0],pull->dims[1],pull->dims[2]);

  if (pull->r > 0.001) 
    pull->bCyl = TRUE;
  else
    pull->bCyl = FALSE;
}








