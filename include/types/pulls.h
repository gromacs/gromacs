/* It would be more logical to have a t_pullgrp and an array of those
   for each type of groups than a t_pullgrps for each type with an array
   of attributes for each group in that type. Oh well. PT 
*/

typedef enum {eStart, eAfm, eConstraint, eUmbrella, eTest} t_runtype;
typedef enum {eCom, eComT0, eDyn, eDynT0} t_reftype;

typedef struct {
  int        n;         /* number of groups */
  atom_id    **idx;     /* atoms indices in the groups */
  real       **weights; /* optional weight (used for switch function) */
  int        *ngx;      /* nr of atoms in the groups */
  char       **grps;    /* names of the groups */
  real       *tmass;    /* total mass of the groups */
  rvec       **x0;      /* coordinates at t=0 */
  rvec       **xp;      /* coordinates at previous step */
  rvec       *x_con;    /* center of mass, obeying constraints */
  rvec       *f;        /* forces due to the pulling/constraining */
  rvec       *spring;   /* coordinates of the springs (eAfm) */
  rvec       *x_unc;    /* center of mass before constraining */
  rvec       *x_ref;    /* reference positions */
  rvec       *dir;      /* direction of constraint */
  real       *d_ref;    /* reference distance  */
  rvec       *xstart;   /* starting point for structure generation */
  rvec       *xprev;    /* position of coms in last written structure */
  rvec       **comhist; /* com over the last nhist steps (for running aver) */
} t_pullgrps; 

typedef struct {
  t_pullgrps dyna;      /* dynamic groups for use with local constraints */
  t_pullgrps pull;      /* groups to pull/restrain/etc/ */
  t_pullgrps ref;       /* reference group, reaction force grps */
  t_runtype  runtype;   /* 0: startstructure 1: afmrun 2: constraint 
			   3: umbrella sampling 4: particle insertion */
  t_reftype  reftype;   /* 0: com 1: com at t=0 2: dynamic */
  rvec       dir;       /* used to select components for constraint */
  rvec       coor;      /* reaction coordinate */
  real       r;         /* radius of cylinder for dynamic COM */
  real       rc;        /* radius of cylinder including switch length */
  int        bRot[3];   /* rotation around x, y, z? */
  real       rot_rate;  /* rate of rotation, for startstructure run */
  real       xlt_rate;  /* rate of translation, for startstructure run */
  int        rot_incr;  /* write out structure every rot_incr degrees */
  real       xlt_incr;  /* write out structure every xlt_incr nm */
  real       tolerance; /* tolerance for reaching desired coordinates (nm) */
  bool       bPull;     /* true if we're doing any pulling */
  bool       bCyl;      /* true if we're using dynamic ref. groups */
  bool       bReverse;  /* reverse reference direction */
  FILE       *out;      /* output file for pull data */
  real       k;         /* force constant for atoms */
  real       rate;      /* pull rate, in nm/timestep */
  real       um_width;  /* width umbrella potential */  
  int        update;    /* update frequency for dynamic grps */
  int        reflag;    /* running average over reflag steps for com */
  bool       bVerbose;  /* be loud and noise */
} t_pull;



