/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
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
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <time.h>
#include "typedefs.h"
#include "string2.h"
#include "smalloc.h"
#include "names.h"
#include "confio.h"
#include "mvdata.h"
#include "txtdump.h"
#include "vec.h"
#include "time.h"
#include "nrnb.h"
#include "mshift.h"
#include "mdrun.h"
#include "update.h"
#include "physics.h"
#include "rmpbc.h"
#include "nrjac.h"
#include "edsam.h"
#include "mpelogging.h"

/* Defines used for flooding */
#define FIT 1
#define BLOWUP 2
#define DOFIT
#define NODEBUG BLOWUP
#ifdef DEBUG
#define DUMP_FORCES 
#define DEBUG_PRINT(X) fprintf(stderr,"%s\n",(X))
#else
#define DEBUG_PRINT(X) {}
#endif

#ifdef DUMP_FORCES
static FILE* logfile = NULL;
#endif
/* end of defines used for flooding */


/* We use the same defines as in mvdata.c here */
#define  block_bc(cr,   d) gmx_bcast(     sizeof(d),     &(d),(cr))
#define nblock_bc(cr,nr,d) gmx_bcast((nr)*sizeof((d)[0]), (d),(cr))
#define   snew_bc(cr,d,nr) { if (!MASTER(cr)) snew((d),(nr)); }


/* enum to identify the type of ED: none, normal ED, flooding */
enum {eEDnone, eEDedsam, eEDflood, eEDnr};

/* enum to identify operations on reference, average, origin, target structures */
enum {eedREF, eedAV, eedORI, eedTAR, eedNR};


typedef struct 
{
    rvec *x;
    rvec *transvec;
    rvec *forces_cartesian;
} t_edlocals;


typedef struct
{
    int    neig;     /* nr of eigenvectors             */
    int   *ieig;     /* index nrs of eigenvectors      */
    real  *stpsz;    /* stepsizes (per eigenvector)    */
    rvec  **vec;     /* eigenvector components         */
    real  *xproj;    /* instantaneous x projections    */
    real  *vproj;    /* instantaneous v projections    */
    real  *fproj;    /* instantaneous f projections    */
    real  *refproj;  /* starting or target projecions  */
    real  radius;    /* instantaneous radius           */
} t_eigvec;


typedef struct
{
    t_eigvec      mon;            /* only monitored, no constraints       */
    t_eigvec      linfix;         /* fixed linear constraints             */
    t_eigvec      linacc;         /* acceptance linear constraints        */
    t_eigvec      radfix;         /* fixed radial constraints (exp)       */
    t_eigvec      radacc;         /* acceptance radial constraints (exp)  */  
    t_eigvec      radcon;         /* acceptance rad. contraction constr.  */
} t_edvecs;


typedef struct
{ 
    real deltaF0;
    bool bHarmonic;
    real tau;
    real deltaF;
    real Efl;
    real kT; 
    real Vfl;
    real dt;
    real constEfl;
    real alpha2; 
    int flood_id;
    t_edlocals loc;
    t_eigvec vecs;   /* use flooding for these */
} t_edflood;


/* This type is for the average, reference, target, and origin structure    */
typedef struct gmx_edx
{
    int           nr;             /* number of atoms this structure contains  */
    int           nr_loc;         /* number of atoms on local node            */
    int           *anrs;          /* atom index numbers                       */
    int           *anrs_loc;      /* local atom index numbers                 */
    int           *c_ind;         /* at which position of the whole anrs 
                                   * array is a local atom?, i.e. 
                                   * c_ind[0...nr_loc-1] gives the atom index 
                                   * with respect to the collective 
                                   * anrs[0...nr-1] array                     */
    rvec          *x;             /* positions                                */
    real          *m;             /* masses                                   */
    real          mtot;           /* total mass (only used in sref)           */
    real          *sqrtm;         /* sqrt of the masses used for mass-
                                   * weighting of analysis (only used in sav) */
    int           c_pbcatom;      /* PBC: Put each atom of this structure at  
                                   * the periodic image which is closest to 
                                   * the coordinates of the pbcatom, 
                                   * coords are at c_ind[c_pbcatom]          */
} t_gmx_edx;


typedef struct edpar
{
    int            nini;           /* total Nr of atoms                    */
    bool           fitmas;         /* true if trans fit with cm            */
    bool           pcamas;         /* true if mass-weighted PCA            */
    int            presteps;       /* number of steps to run without any   
                                    *    perturbations ... just monitoring */
    int            outfrq;         /* freq (in steps) of writing to edo    */
    int            maxedsteps;     /* max nr of steps per cycle            */

    /* all gmx_edx datasets are copied to all nodes in the parallel case    */
    struct gmx_edx sref;           /* reference positions, to these fitting
                                    * will be done                         */
    bool           bRefEqAv;       /* If true, reference & average indices
                                    * are the same. Used for optimization  */
    struct gmx_edx sav;            /* average positions                    */
    struct gmx_edx star;           /* target positions                     */
    struct gmx_edx sori;           /* origin positions                     */

    t_edvecs       vecs;           /* eigenvectors                         */
    real           slope;          /* minimal slope in acceptance radexp   */

    /* the following will be obsolete once flooding is parallelized */
    int           ned;            /* Nr of atoms in essdyn buffer         */
    int           nmass;          /* Nr of masses                         */ /* = edi->sref.nr   */
    int           *masnrs;        /* index nrs for atoms with masses      */ /* = edi->sref.anrs */
    real          *mass;          /* atomic masses                        */ /* to be put in edi->sref.mass */
    real          tmass;          /* total mass                           */ /* to be put in edi->sref.tmass */
    int           nfit;           /* Number of atoms to use for rot fit   */ /* = edi->sref.nr */
    int           nfit_loc;       /* ... and how many of these are on           = edi->sref.nr_loc  
                                     the local node                       */
    int           *fitnrs;        /* index nrs of atoms to use for              =edi->sref.anrs
                                   * rotational fit, this is a copy of
                                   * the reference structure edi->sref    */
    int           *fitnrs_loc;    /* fitnrs_loc[0...nfit_loc-1] give the 
                                   * atom index with respect to the whole 
                                   * fitnrs[0...nfit-1] array             */
    /* end obsolete */

    bool           bNeedDoEdsam;   /* if any of the options mon, linfix, ...
                                    * is used (i.e. apart from flooding)   */
    t_edflood      flood;          /* parameters especially for flooding   */
    struct t_ed_buffer *buf;       /* handle to local buffers              */
    struct edpar   *next_edi;      /* Pointer to another ed dataset        */
} t_edpar;


typedef struct gmx_edsam
{
    int           eEDtype;        /* Type of ED: see enums above          */
    char          *edinam;        /* name of ED sampling input file       */
    char          *edonam;        /*                     output           */
    FILE          *edo;           /* output file pointer                  */
    t_edpar       *edpar;
    int           ePBC;           /* PBC type from inputrec               */
    t_pbc         *pbc;           /* PBC information                      */
} t_gmx_edsam;


struct t_do_edsam
{
    matrix old_rotmat;
    real oldrad;
    rvec old_transvec,older_transvec,transvec_compact;
    rvec *xc_dum;
    rvec *xcoll;         /* Coordinates from all nodes, this is the collective set of coords we work on.
                          * It coordinates of atoms with average structure indices */
    rvec *xc_ref;        /* same but with reference structure indices */
    ivec *shifts_xcoll;  /* Shifts for xcoll  */
    ivec *shifts_xc_ref; /* Shifts for xc_ref */
};


/* definition of ED buffer structure */
struct t_ed_buffer
{
    struct t_fitit *                fitit;
    struct t_do_edfit *             do_edfit;
    struct t_remove_pbc_effect *    remove_pbc_effect;
    struct t_do_edsam *             do_edsam;
    struct t_do_radcon *            do_radcon;
};



/* Does not subtract average positions, projection on single eigenvector is returned
 * used by: do_linfix, do_linacc, do_radfix, do_radacc, do_radcon
 * Average position is subtracted in ed_apply_constraints prior to calling projectx
 */
static real projectx(t_edpar *edi, rvec *xcoll, rvec *vec, t_commrec *cr)
{
    int  i;
    real proj=0.0;


    for (i=0; i<edi->sav.nr; i++)
        proj += edi->sav.sqrtm[i]*iprod(vec[i], xcoll[i]);

    return proj;
}


/* Same as projectx, but mass-weighting is applied differently -> forces */
static real projectf(t_edpar *edi, rvec *xcoll, rvec *vec, t_commrec *cr)
{
    int  i;
    real proj=0.0;


    for (i=0; i<edi->sav.nr; i++)
        proj += iprod(vec[i], xcoll[i])/edi->sav.sqrtm[i];

    return proj;
}


/* specialized: projection is stored in vec->refproj
 * ---> used for radacc, radfix, radcon  and center of flooding potential
 * subtracts average positions, projects vector x, uses atoms sav.anrs[i] of x */
static void rad_project(t_edpar *edi, rvec *x, t_eigvec *vec, t_commrec *cr)
{
    int i,j,k;
    real rad=0.0;

    /* subtract average positions */
    for (i = 0; i < edi->sav.nr; i++)
        rvec_dec(x[i], edi->sav.x[i]);

    for (i = 0; i < vec->neig; i++)
    {
        vec->refproj[i] = projectx(edi,x,vec->vec[i],cr);
        rad += pow((vec->refproj[i]-vec->xproj[i]),2);
    }
    vec->radius=sqrt(rad);

    /* add average positions */
    for (i = 0; i < edi->sav.nr; i++) 
        rvec_inc(x[i], edi->sav.x[i]);
}


/* New version of project_to_eigenvectors */
static void project_to_eigvectors_p(rvec       *x,    /* The coordinates to project to an eigenvector */ 
                                    t_eigvec   *vec,  /* The eigenvectors */
                                    t_edpar    *edi, 
                                    char       *mode, 
                                    t_commrec  *cr)
{
    int  i,j,k;
    real proj;


    if (!vec->neig) return;

    /* subtract average positions */
    if (strcmp(mode,"x") == 0)
    {
        for (i=0; i<edi->sav.nr; i++) 
            rvec_dec(x[i], edi->sav.x[i]);
    }

    for (i=0; i<vec->neig; i++)
    { 
        if      (strcmp(mode,"x") == 0) vec->xproj[i] = projectx(edi, x, vec->vec[i], cr);
        else if (strcmp(mode,"v") == 0) vec->vproj[i] = projectx(edi, x, vec->vec[i], cr);
        else if (strcmp(mode,"f") == 0) vec->fproj[i] = projectf(edi, x, vec->vec[i], cr);
        /* this has no influence on flooding forces, since this routine is called from the 
         * edsam branch completely disconnected from the do_flood branch */
    }

    /* add average positions */
    if (strcmp(mode,"x") == 0)
    {
        for (i=0; i<edi->sav.nr; i++) 
            rvec_inc(x[i], edi->sav.x[i]);
    }
}


/* project vector x, uses atoms sav.anrs[i] of x, 
 * if mode  is "x" it subtract average positions prior to projection
 * and add them afterwards to retain the unchanged vector x
 * mode = "x","v","f" -> store in xproj, vproj or fproj
 * XXX mass-weighting is applied
 */
static void project_to_eigvectors(rvec *x, t_eigvec *vec, t_edpar *edi, char *mode, t_commrec *cr)
{
    int i,j,k;
    real proj;
    
    
    if (!vec->neig) return;
    /* subtract average positions */
    if (strcmp(mode,"x") == 0)
    {
        for (i=0;(i<edi->sav.nr);i++) 
            rvec_dec(x[edi->sav.anrs[i]],edi->sav.x[i]);
    }

    for (i=0;(i<vec->neig);i++)
    {
        if      (strcmp(mode,"x") == 0) vec->xproj[i]=projectx(edi,x,vec->vec[i],cr);
        else if (strcmp(mode,"v") == 0) vec->vproj[i]=projectx(edi,x,vec->vec[i],cr);
        else if (strcmp(mode,"f") == 0) vec->fproj[i]=projectf(edi,x,vec->vec[i],cr);
        /* this has no influence on flooding forces, since this routine is called from the edsam branch
       completely disconnected from the do_flood branch */
    }

    /* add average positions */
    if (strcmp(mode,"x") == 0)
    {
        for (i=0;(i<edi->sav.nr);i++) 
            rvec_inc(x[edi->sav.anrs[i]],edi->sav.x[i]);
    }
}


/* New version of project */
static void project_p(rvec      *x,     /* coordinates to project */ 
                      t_edpar   *edi,   /* edi data set */
                      char      *mode, 
                      t_commrec *cr)
{
    int i;

    /* make the projections */
    /* it is not more work to subtract the average position in every subroutine again, 
     because these routines are rarely used simultanely */
    project_to_eigvectors_p(x, &edi->vecs.mon   , edi, mode, cr);
    project_to_eigvectors_p(x, &edi->vecs.linfix, edi, mode, cr);
    project_to_eigvectors_p(x, &edi->vecs.linacc, edi, mode, cr);
    project_to_eigvectors_p(x, &edi->vecs.radfix, edi, mode, cr);
    project_to_eigvectors_p(x, &edi->vecs.radacc, edi, mode, cr);
    project_to_eigvectors_p(x, &edi->vecs.radcon, edi, mode, cr);
}


/* wrapper: project vector x onto all edi->vecs (mon, linfix,...) */
static void project(rvec *x, t_edpar *edi, char *mode, t_commrec *cr)
{
    int i;
    /* make the projections */
    /* it is not more work to subtract the average position in every subroutine again, 
     because these routines are rarely used simultanely */
    project_to_eigvectors(x,&edi->vecs.mon   ,edi,mode, cr);
    project_to_eigvectors(x,&edi->vecs.linfix,edi,mode, cr);
    project_to_eigvectors(x,&edi->vecs.linacc,edi,mode, cr);
    project_to_eigvectors(x,&edi->vecs.radfix,edi,mode, cr);
    project_to_eigvectors(x,&edi->vecs.radacc,edi,mode, cr);
    project_to_eigvectors(x,&edi->vecs.radcon,edi,mode, cr);
}


static real calc_radius(t_eigvec *vec)
{
    int i;
    real rad=0.0;


    for (i=0; i<vec->neig; i++) 
        rad += pow((vec->refproj[i]-vec->xproj[i]),2);

    return rad=sqrt(rad);
}


/* Debug helper */
static void dump_edi_positions(FILE *out, struct gmx_edx *s, char name[])
{
    int i;

    
    fprintf(out, "#%s coordinates:\n%d\n", name, s->nr);
    if (s->nr == 0)
        return;

    fprintf(out, "#index, x, y, z");
    if (s->sqrtm)
        fprintf(out, ", sqrt(m)");
    for (i=0; i<s->nr; i++) 
    {
        fprintf(out, "\n%6d  %11.6f %11.6f %11.6f",s->anrs[i], s->x[i][XX], s->x[i][YY], s->x[i][ZZ]);
        if (s->sqrtm)
            fprintf(out,"%9.3f",s->sqrtm[i]);
    }
    fprintf(out, "\n");
}


/* Debug helper */
static void dump_edi_eigenvecs(FILE *out, t_eigvec *ev, char name[], int length)
{
    int i,j;

    
    fprintf(out, "#%s eigenvectors:\n%d\n", name, ev->neig);
    /* Dump the data for every eigenvector: */
    for (i=0; i<ev->neig; i++)
    {
        fprintf(out, "EV %4d\ncomponents %d\nstepsize %f\nxproj %f\nvproj %f\nfproj %f\nrefproj %f\nradius %f\nComponents:\n", 
                ev->ieig[i], length, ev->stpsz[i], ev->xproj[i], ev->vproj[i], ev->fproj[i], ev->refproj[i], ev->radius);
        for (j=0; j<length; j++)
            fprintf(out, "%11.6f %11.6f %11.6f\n", ev->vec[i][j][XX], ev->vec[i][j][YY], ev->vec[i][j][ZZ]);
    }
}


/* Debug helper */
static void dump_edi(t_edpar *edpars, t_commrec *cr)
{
    int   i;
    FILE  *out;
    char  fname[255];


    sprintf(fname, "EDdump%.2d", cr->nodeid);
    out = fopen(fname, "w");

    fprintf(out,"#NINI\n %d\n#SELMAS\n %d\n#ANALYSIS_MAS\n %d\n",
            edpars->nini,edpars->fitmas,edpars->pcamas);
    fprintf(out,"#OUTFRQ\n %d\n#MAXLEN\n %d\n#SLOPECRIT\n %f\n",
            edpars->outfrq,edpars->maxedsteps,edpars->slope);
    fprintf(out,"#PRESTEPS\n %d\n#DELTA_F0\n %f\n#TAU\n %f\n#EFL_NULL\n %f\n#ALPHA2\n %f\n",
            edpars->presteps,edpars->flood.deltaF0,edpars->flood.tau,edpars->flood.constEfl,edpars->flood.alpha2);

    /* Dump reference, average, target, origin positions */
    dump_edi_positions(out, &edpars->sref, "REFERENCE");
    dump_edi_positions(out, &edpars->sav , "AVERAGE"  );
    dump_edi_positions(out, &edpars->star, "TARGET"   );
    dump_edi_positions(out, &edpars->sori, "ORIGIN"   );

    /* Dump eigenvectors */
    dump_edi_eigenvecs(out, &edpars->vecs.mon   , "MONITORED", edpars->sav.nr);
    dump_edi_eigenvecs(out, &edpars->vecs.linfix, "LINFIX"   , edpars->sav.nr);
    dump_edi_eigenvecs(out, &edpars->vecs.linacc, "LINACC"   , edpars->sav.nr);
    dump_edi_eigenvecs(out, &edpars->vecs.radfix, "RADFIX"   , edpars->sav.nr);
    dump_edi_eigenvecs(out, &edpars->vecs.radacc, "RADACC"   , edpars->sav.nr);
    dump_edi_eigenvecs(out, &edpars->vecs.radcon, "RADCON"   , edpars->sav.nr);

    /* Dump flooding eigenvectors */
    dump_edi_eigenvecs(out, &edpars->flood.vecs, "FLOODING", edpars->sav.nr);

    /* Dump ed local buffer */
    fprintf(out, "buf->fitit            =%p\n", edpars->buf->fitit            );
    fprintf(out, "buf->do_edfit         =%p\n", edpars->buf->do_edfit         );
    fprintf(out, "buf->remove_pbc_effect=%p\n", edpars->buf->remove_pbc_effect);
    fprintf(out, "buf->do_edsam         =%p\n", edpars->buf->do_edsam         );
    fprintf(out, "buf->do_radcon        =%p\n", edpars->buf->do_radcon        );

    fclose(out);
}


/* Debug helper */
static void dump_rotmat(FILE* out,matrix rotmat)
{
    fprintf(out,"ROTMAT: %f %f %f\n",rotmat[XX][XX],rotmat[XX][YY],rotmat[XX][ZZ]);
    fprintf(out,"ROTMAT: %f %f %f\n",rotmat[YY][XX],rotmat[YY][YY],rotmat[YY][ZZ]);
    fprintf(out,"ROTMAT: %f %f %f\n",rotmat[ZZ][XX],rotmat[ZZ][YY],rotmat[ZZ][ZZ]);
}


/* Debug helper */
static void dump_rvec(FILE *out, int dim, rvec *x) 
{
    int i;

    
    for (i=0; i<dim;i++)
        fprintf(out,"%4d   %f %f %f\n",i,x[i][XX],x[i][YY],x[i][ZZ]);
}


/* Debug helper */
static void dump_mat(FILE* out, int dim, double** mat)
{
    int i,j;

    
    fprintf(out,"MATRIX:\n");
    for (i=0;i<dim;i++)
    {
        for (j=0;j<dim;j++)
            fprintf(out,"%f ",mat[i][j]);
        fprintf(out,"\n");
    }
}


struct t_do_edfit { 
    double **omega;
    double **om;
};

static void do_edfit(int natoms,rvec *xp,rvec *x,matrix R,t_edpar *edi)
{
    /* this is a copy of do_fit with some modifications */
    int    c,r,n,j,i,irot;
    double d[6],xnr,xpc;
    matrix vh,vk,u;
    int    index;
    real   max_d;

    struct t_do_edfit *loc;
    bool bFirst;

    if(edi->buf->do_edfit != NULL)
    {
        bFirst = FALSE;
        loc    = edi->buf->do_edfit;
    }
    else
    {
        bFirst = TRUE;
        snew(edi->buf->do_edfit,1);
    }
    loc = edi->buf->do_edfit;

    if (bFirst) 
    {
        snew(loc->omega,2*DIM);
        snew(loc->om,2*DIM);
        for(i=0; i<2*DIM; i++)
        {
            snew(loc->omega[i],2*DIM);
            snew(loc->om[i],2*DIM);
        }
    }

    for(i=0;(i<6);i++)
    {
        d[i]=0;
        for(j=0;(j<6);j++)
        {
            loc->omega[i][j]=0;
            loc->om[i][j]=0;
        }
    }

    /* calculate the matrix U */
    clear_mat(u);
    for(n=0;(n<natoms);n++)
    {
        for(c=0; (c<DIM); c++)
        {
            xpc=xp[n][c];
            for(r=0; (r<DIM); r++)
            {
                xnr=x[n][r];
                u[c][r]+=xnr*xpc;
            }
        }
    }

    /* construct loc->omega */
    /* loc->omega is symmetric -> loc->omega==loc->omega' */
    for(r=0;(r<6);r++)
        for(c=0;(c<=r);c++)
            if ((r>=3) && (c<3))
            {
                loc->omega[r][c]=u[r-3][c];
                loc->omega[c][r]=u[r-3][c];
            }
            else
            {
                loc->omega[r][c]=0;
                loc->omega[c][r]=0;
            }

    /* determine h and k */
    DEBUG_PRINT("call jacobi");

#ifdef DEBUG
    {
        int i;
        dump_mat(stderr,2*DIM,loc->omega);
        for (i=0; i<6; i++)
            fprintf(stderr,"d[%d] = %f\n",i,d[i]);
    }
#endif
    jacobi(loc->omega,6,d,loc->om,&irot);

    if (irot==0)
        fprintf(stderr,"IROT=0\n");

    index=0; /* For the compiler only */

    for(j=0;(j<3);j++)
    {
        max_d=-1000;
        for(i=0;(i<6);i++)
            if (d[i]>max_d)
            {
                max_d=d[i];
                index=i;
            }
        d[index]=-10000;
        for(i=0;(i<3);i++)
        {
            vh[j][i]=M_SQRT2*loc->om[i][index];
            vk[j][i]=M_SQRT2*loc->om[i+DIM][index];
        }
    }

    /* determine R */
    for(c=0;(c<3);c++)
        for(r=0;(r<3);r++)
            R[c][r]=vk[0][r]*vh[0][c]+
            vk[1][r]*vh[1][c]+
            vk[2][r]*vh[2][c];
    if (det(R) < 0)
        for(c=0;(c<3);c++)
            for(r=0;(r<3);r++)
                R[c][r]=vk[0][r]*vh[0][c]+
                vk[1][r]*vh[1][c]-
                vk[2][r]*vh[2][c];
}


static void rotate_x(int nr,rvec *x,matrix rmat)
{
    int i,j,k;
    rvec x_old;

    DEBUG_PRINT(" apply the rotation matrix \n");
    for(i=0;(i<nr);i++)
    {
        for(j=0;(j<3);j++)
            x_old[j]=x[i][j];
        for(j=0;(j<3);j++)
        {
            x[i][j]=0;
            for(k=0;(k<3);k++)
                x[i][j]+=rmat[k][j]*x_old[k];
        }
    }
}


struct t_fitit {
    rvec *xdum1;
    int  nr1;
    rvec *xdum2;
    int  nfit;
};

/* fits x[edi->fitnrs[i]] bzw x[edi->masnrs[i]] which should be equivalent
    X has to contain ALL Atoms */
static void fitit(int nr, rvec *x,t_edpar *edi,rvec *transvec,matrix rmat, t_commrec *cr)
{
    int i,j,k;
    bool bFirst;
    struct t_fitit *buf;

    
    if(edi->buf->fitit != NULL)
    {
        bFirst = FALSE;
        buf    = edi->buf->fitit;
    }
    else
    {
        bFirst = TRUE;
        snew(edi->buf->fitit,1);
    }
    buf = edi->buf->fitit;

    if (bFirst)
    {
        snew(buf->xdum1,nr);
        buf->nr1=nr;
    }
    if (nr>buf->nr1)
    {
        srenew(buf->xdum1,nr);
        buf->nr1=nr;
    }

    for(i=0; (i<nr); i++)
        copy_rvec(x[i],buf->xdum1[i]); 

    DEBUG_PRINT(" first do translational fit \n");
    /* put_in_origin(nr,x,edi->nmass,edi->masnrs,edi->mass,edi->tmass); */

    /* determine transvec from difference after translational fit */
    for(i=0; (i<nr); i++) 
        rvec_sub(x[i],buf->xdum1[i],transvec[i]); 


    DEBUG_PRINT(" now rotational fit \n");
    if (bFirst) {
        buf->nfit = edi->nfit;
        snew(buf->xdum2,edi->nfit);
    }
    if (buf->nfit < edi->nfit) { /* happens in flooding with more than one matrix */
        srenew(buf->xdum2,edi->nfit);
        buf->nfit=edi->nfit;
    }
    for(i=0; (i<edi->nfit); i++)
        copy_rvec(x[edi->fitnrs[i]],buf->xdum2[i]);
    DEBUG_PRINT("do_edfit..");
    do_edfit(edi->nfit,edi->sref.x,buf->xdum2,rmat,edi);
    rotate_x(nr,x,rmat);
}


static void rmrotfit(int nat, rvec *xcoll, matrix rotmat)
{
    int    i,j,k;
    matrix r_inv;
    rvec   xdum;


    /* invert the rotation matrix and apply */
    for (i=0; i<nat; i++)
    {
        for (j=0; j<3; j++)
            xdum[j]=xcoll[i][j];
        for (j=0; j<3; j++)
        {
            xcoll[i][j]=0;
            for (k=0; k<3; k++)
                xcoll[i][j] += rotmat[j][k]*xdum[k];
        }
    }
}


static void rmtransfit(int nat, rvec *xcoll, rvec transvec) 
{ 
    int i;

    
    /* subtract the translation vector */
    for(i=0; i<nat; i++)
        rvec_dec(xcoll[i], transvec);
}


static void rmfit(int nat, rvec *xcoll, rvec transvec, matrix rotmat) 
{
    rmrotfit(nat, xcoll, rotmat);
    rmtransfit(nat, xcoll, transvec);
}


/**********************************************************************************
 ******************** FLOODING ****************************************************
 **********************************************************************************
 
The flooding ability was added later to edsam. Many of the edsam functionality could be reused for that purpose. 
The flooding covariance matrix, i.e. the selected eigenvectors and their corresponding eigenvalues are 
read as 7th Component Group. The eigenvalues are coded into the stepsize parameter (as used by -linfix or -linacc). 

do_md clls right in the beginning the function init_edsam, which reads the edi file, saves all the necessary information in
the edi structure and calls init_flood, to initialise some extra fields in the edi->flood structure.

since the flooding acts on forces do_flood is called from the function force() (force.c), while the other 
edsam functionality is hooked into md via the update() (update.c) function acting as constraint on positions. 
flooding works correctly even if do_edsam() is not called.

do_flood makes a copy of the positions,
fits them, projects them computes flooding_energy, and flooding forces. The forces are computed in the 
space of the eigenvectors and are then blown up to the full cartesian space and rotated back to remove the
fit. Then do_flood adds these forces to the forcefield-forces
(given as parameter) and updates the adaptive flooding parameters Efl and deltaF.

To center the flooding potential at a different location one can use the -ori option in make_edi. The ori
structure is projected to the system of eigenvectors and then this position in the subspace is used as
center of the flooding potential.   If the option is not used, the center will be zero in the subspace,
i.e. the average structure as given in the make_edi file.

To use the flooding potential as restraint, make_edi has the option -restrain, which leads to inverted
signs of alpha2 and Efl, such that the sign in the exponential of Vfl is not inverted but the sign of
Vfl is inverted. Vfl = Efl * exp (- .../Efl/alpha2*x^2...) With tau>0 the negative Efl will grow slowly
so that the restraint is switched off slowly. When Efl==0 and inverted flooding is ON is reached no
 further adaption is applied, Efl will stay constant at zero. 

to use restraints with harmonic potentials switch -restrain and -harmonic. Then the eigenvalues are 
used as spring constants for the harmonic potential. 
Note that eq3 in the flooding paper (J. Comp. Chem. 2006, 27, 1693-1702) defines the parameter lambda \
as the inverse of the spring constant,
whereas the implementation uses lambda as the spring constant.

to use more than one flooding matrix just concatenate severale .edi files (cat flood1.edi flood2.edi > flood_all.edi )
the routine read_edi_file reads all of theses flooding files.
The structure t_edi is now organized as a list of t_edis  and the function do_flood cycles through the list 
calling the do_single_flood() routine for every single entry. Since every state variables have been kept in one 
edi there is no interdependence whatsoever. The forces are added together. 

  To write energies into the .edr file, call the function 
        get_flood_enx_names(char**, int *nnames) to get the Header (Vfl1 Vfl2... Vfln)
and call
        get_flood_energies(real Vfl[],int nnames); 

  TODO:
- one could program the whole thing such that Efl, Vfl and deltaF is written to the .edr file. -- i dont know how to do that, yet.

  At the moment one can have multiple flooding matrices, but only the first input is used for edsam. 
  Especially with the angular motion remover there might be a market for multiple edsam inputs as well. 

  Secondly: Maybe one should give a range of atoms for which to remove motion, so that motion is removed with 
  two edsam files from two peptide chains
*/

static void write_edo_flood(t_edpar *edi, FILE *fp_edo, int step) 
{
    int i;

    fprintf(fp_edo,"%d.th FL: %d %g %g %g\n",edi->flood.flood_id,step, edi->flood.Efl, edi->flood.Vfl, edi->flood.deltaF);
    fprintf(fp_edo,"FL_FORCES: ");
    for (i=0;i<edi->flood.vecs.neig;i++)
        fprintf(fp_edo," %f",edi->flood.vecs.fproj[i]);
    fprintf(fp_edo,"\n");
    fflush(fp_edo);
}


/* project fitted structure onto supbspace -> store in edi->vec.flood.xproj */
static void flood_project(rvec *x, t_edpar *edi, t_commrec *cr)
{
    /* projects the positions onto the subspace */
    int i;

    /* do projection */
    project_to_eigvectors(x,&edi->flood.vecs,edi,"x",cr);
}


/* from flood.xproj compute the Vfl(x) at this point -> store in edi->flood.Vfl*/
static real flood_energy(t_edpar *edi)
{
    /* compute flooding energy Vfl
     Vfl = Efl * exp( - \frac {kT} {2Efl alpha^2} * sum_i { \lambda_i c_i^2 } )
     \lambda_i is the reciproce eigenvalue 1/\sigma_i
         it is already computed by make_edi and stored in stpsz[i]
     bHarmonic:
       Vfl = - Efl * 1/2(sum _i {\frac 1{\lambda_i} c_i^2})
     */
    real summe;
    int i;

    summe=0.0;
    /*compute sum which will be the exponent of the exponential */
    if (edi->flood.bHarmonic)
        for (i=0;i<edi->flood.vecs.neig; i++)
            summe+=edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i])*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]);
    else
        for (i=0;i<edi->flood.vecs.neig; i++)
            summe+=edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i])*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]);
#ifdef DUMP_FORCES
    fprintf(logfile, "REFPROJ: ");
    for (i=0;i<edi->flood.vecs.neig; i++)
        fprintf(logfile, "%f ",edi->flood.vecs.refproj[i]);
    fprintf(logfile, "\n");

    fprintf(logfile, "XPROJ: ");
    for (i=0;i<edi->flood.vecs.neig; i++)
        fprintf(logfile, "%f ",edi->flood.vecs.xproj[i]);
    fprintf(logfile, "\n");

    fprintf(logfile, "SUMME: %f kT : %f alpha2 : %f Efl %f\n ", summe, edi->flood.kT, edi->flood.alpha2, edi->flood.Efl);
#endif

    /*compute the gauss function*/
    if (edi->flood.bHarmonic)
        edi->flood.Vfl=  - 0.5*edi->flood.Efl*summe;  /* minus sign because Efl is negativ, if restrain is on. */
    else
        edi->flood.Vfl= edi->flood.Efl!=0 ? edi->flood.Efl*exp(-edi->flood.kT/2/edi->flood.Efl/edi->flood.alpha2*summe) :0;
        return edi->flood.Vfl;
}


/* from the position and the Vfl compute forces in subspace -> stored in edi->vec.flood.fproj */
static void flood_forces(t_edpar *edi)
{
    /* compute the forces in the subspace of the flooding eigenvectors
   by the formula F_i= V_{fl}(c) * ( \frac {kT} {E_{fl}} \lambda_i c_i */
    int i;
    real energy=edi->flood.Vfl;
#ifdef DUMP_FORCES
    fprintf(logfile, "Vfl= %f, Efl= %f, xproj= %f, refproj= %f\n",energy, edi->flood.Efl,edi->flood.vecs.xproj[0],edi->flood.vecs.refproj[0]);
#endif
    if (edi->flood.bHarmonic)
        for (i=0; i<edi->flood.vecs.neig; i++)
        {
            edi->flood.vecs.fproj[i]= edi->flood.Efl* edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]);
#ifdef DUMP_FORCES
            fprintf(logfile, "%f ",edi->flood.vecs.fproj[i]);
#endif
        }
    else
        for (i=0; i<edi->flood.vecs.neig; i++)
        {
            /* if Efl is zero the forces are zero if not use the formula */
            edi->flood.vecs.fproj[i]= edi->flood.Efl!=0 ? edi->flood.kT/edi->flood.Efl/edi->flood.alpha2*energy*edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]) : 0;
#ifdef DUMP_FORCES
            fprintf(logfile, "force %f ",edi->flood.vecs.fproj[i]);
#endif
        }
#ifdef DUMP_FORCES
    fprintf(logfile,"\n");
#endif
}


/* raise forces from subspace into cartesian space */
static void flood_blowup(t_edpar *edi, rvec *forces_cart)
{
    /* this function lifts the forces from the subspace to the cartesian space
     all the values not contained in the subspace are assumed to be zero and then 
     a coordinate transformation from eigenvector to cartesian vectors is performed 
     The nonexistent values don't have to be set to zero explicitly, they would occur 
     as zero valued summands, hence we just stop to compute this part of the sum.

     for every atom we add all the contributions to this atom from all the different eigenvectors.

     NOTE: one could add directly to the forcefield forces, would mean we wouldn't have to clear the 
     field forces_cart prior the computation, but momentarily we want to compute the forces seperately 
     to have them accessible for diagnostics
     */
    int i,j,eig;
    rvec dum;
    real *forces_sub;
    forces_sub=edi->flood.vecs.fproj;
    /* clear forces first */
    for (j=0; j<edi->ned; j++) 
        clear_rvec(forces_cart[j]);
#if (DEBUG==BLOWUP)
    fprintf(stderr,"cleared cartesian force vector:");
    dump_rvec(stderr, edi->ned, forces_cart);
#endif
    /* now compute atomwise */
    for (j=0; j<edi->sav.nr; j++)
    {
        /* should this be sav.nr ??? */
        /* compute    forces_cart[edi->sav.anrs[j]] */
        for (eig=0; eig<edi->flood.vecs.neig; eig++)
        {
            /* force vector is force * eigenvector compute only atom j */
            svmul(forces_sub[eig],edi->flood.vecs.vec[eig][j],dum);
            /* add this vector to the cartesian forces */
            rvec_inc(forces_cart[edi->sav.anrs[j]],dum);
        }
#ifdef DUMP_FORCES
        fprintf(logfile,"%d %f %f %f\n",
                edi->sav.anrs[j], forces_cart[edi->sav.anrs[j]][XX],forces_cart[edi->sav.anrs[j]][YY],forces_cart[edi->sav.anrs[j]][ZZ]); 
    } 
    fprintf(logfile,"\n--------------------------------\n");
#if 0
    {
#endif
#else
    }
#endif
}


/* update the values of Efl, deltaF depending on tau and Vfl */
static void update_adaption(t_edpar *edi)
{
    /* this function updates the parameter Efl and deltaF according to the rules given in 
   'predicting unimolecular chemical reactions: chemical flooding' M Mueller et al, J. chem Phys.
     */

    if ((edi->flood.tau < 0 ? -edi->flood.tau : edi->flood.tau )>0.00000001)
    {
        edi->flood.Efl=edi->flood.Efl+edi->flood.dt/edi->flood.tau*(edi->flood.deltaF0-edi->flood.deltaF);
        /* check if restrain (inverted flooding) --> don't let EFL become positiv*/
        if (edi->flood.alpha2<0 && edi->flood.Efl>-0.00000001)
            edi->flood.Efl=0;

        edi->flood.deltaF=(1-edi->flood.dt/edi->flood.tau)*edi->flood.deltaF+edi->flood.dt/edi->flood.tau*edi->flood.Vfl;
    }
}


static void do_single_flood(FILE *log,FILE *edo,rvec x_orig[],rvec force[], t_edpar *edi, int step, int nr, t_commrec *cr) 
{  
    int i,j,ned=edi->ned,iupdate=500;
    matrix rotmat;
    real mas,rad;
    t_edpar *actual_edi;


    /*make copy of coordinates */
    for (i=0;i<edi->ned;i++)
    {
        copy_rvec(x_orig[i],edi->flood.loc.x[i]);
        if (!finite(edi->flood.loc.x[i][XX]))
        {
            fprintf(stderr,"Found invalid coordinate: \n");
            dump_rvec(stderr,edi->ned,edi->flood.loc.x);
            gmx_fatal(FARGS,"Something is wrong with the coordinates: a LINCS error? to much flooding strength? wrong flooding vectors?");
        }
    }

    /* fit the structure */
    DEBUG_PRINT("fit the structure...");
#ifdef DEBUG
    fprintf(stderr,"the structure to be fitted:\n");
    dump_rvec(stderr,edi->ned,edi->flood.loc.x);
#endif
#ifdef DOFIT
    fitit(ned,edi->flood.loc.x,edi,edi->flood.loc.transvec,rotmat,cr);
#endif
#ifdef DEBUG
    dump_rotmat(stderr,rotmat);
#endif
    /* put projected values into edi->flood.vecs.xproj */
    DEBUG_PRINT("put projected values into edi->flood.vecs.xproj");
    flood_project(edi->flood.loc.x,edi,cr);

    DEBUG_PRINT("compute_energy");
    flood_energy(edi);
    update_adaption(edi);


    DEBUG_PRINT("compute flooding forces");
    flood_forces(edi);

    /* translate them into cartesian coordinates */
    flood_blowup(edi, edi->flood.loc.forces_cartesian);

    /* rotate forces back so that they correspond to the given structure and not to the fitted one */
#ifdef DOFIT
    rmrotfit(ned, edi->flood.loc.forces_cartesian, rotmat);
#endif
    /* and finally: add forces to the master force variable */
    for (i=0; i<edi->ned; i++)
        rvec_inc(force[i],edi->flood.loc.forces_cartesian[i]);


    if (do_per_step(step,edi->outfrq))
        write_edo_flood(edi,edo,step);
}


extern void do_flood(FILE *log, t_commrec *cr, rvec x_orig[],rvec force[], gmx_edsam_t ed, int step)
{
    /* this is the hook, 
     do_flood is called by mdrun via the force() (force.c) routine.
     the parameters are
     logfile, commrec (to no if we are on the master node), x_orig (positions), 
     force ( forcefield forces , we add the flooding force to them), edi - all the parameters)
     */
    int i,j,ned,iupdate=500;
    int nr;
    matrix rotmat;
    real mas,rad;
    t_edpar *edi;
    t_edpar *actual_edi;

    if (!ed) 
        return;
    if (ed->eEDtype != eEDflood)
        return;

    edi = ed->edpar;
    if (!edi->flood.vecs.neig) return;
    ned=edi->ned;
    nr=0;
    
    if(!MASTER(cr))
        return;
    
    actual_edi=edi;
    while (actual_edi)
    {
        DEBUG_PRINT("call flooding for one matrix");
        do_single_flood(log,ed->edo,x_orig,force,actual_edi,step,nr++,cr);
        actual_edi=actual_edi->next_edi;
    }
}


/* called by init_edi, configure some flooding related variables and structures, 
 * print headers to output files */
static void init_flood(t_edpar *edi, FILE *fp_edo, real dt, t_commrec *cr)
{
    int i;
    matrix rotmat;


    edi->flood.Efl = edi->flood.constEfl;
    edi->flood.Vfl = 0;
    edi->flood.dt  = dt;

    if (edi->flood.vecs.neig)
    {
        gmx_fatal(FARGS, "Sorry - Flooding is not yet implemented in this version.");

#ifdef DUMP_FORCES
        logfile=ffopen("dump_force2.log","w");
        fprintf(logfile,"Vfl Efl\nForce in Subspace\nForces in cartesian space\n");
#endif
        fprintf(fp_edo,"FL_HEADER: Flooding of matrix %d is switched on! The flooding output will have the following format:\n",
                edi->flood.flood_id);
        if (edi->flood.flood_id<1)
            fprintf(fp_edo,"FL_HEADER: Step Efl Vfl deltaF \n");
        /* get memory only once, that should speed up computations */

        snew(edi->flood.loc.x,edi->ned);
        snew(edi->flood.loc.transvec,edi->ned);
        snew(edi->flood.loc.forces_cartesian,edi->ned);

        /* set center of flooding potential */
        if (edi->sori.nr > 0) {
            /* ... to the structure given with -ori */
            fitit(edi->ned, edi->sori.x, edi,edi->flood.loc.transvec, rotmat, cr);
            rad_project(edi, edi->sori.x, &edi->flood.vecs, cr);
        } else
        {
            /* ... to the center of the covariance matrix, i.e. the average structure, i.e. zero in the projected system */
            for (i=0;i<edi->flood.vecs.neig;i++)
                edi->flood.vecs.refproj[i]=0.0;
        }
    }
}


/*********** Energy book keeping ******/
static void get_flood_enx_names(t_edpar *edi, char** names, int *nnames)  /* get header of energies */
{
    t_edpar *actual;
    int count;
    char buf[STRLEN];
    actual=edi;
    count = 1;
    while (actual)
    {
        srenew(names,count);
        sprintf(buf,"Vfl_%d",count);
        names[count-1]=strdup(buf);
        actual=actual->next_edi;
        count++;
    }
    *nnames=count-1;
}


static void get_flood_energies(t_edpar *edi, real Vfl[],int nnames)
{
    /*fl has to be big enough to capture nnames-many entries*/
    t_edpar *actual;
    int count;
    char buf[STRLEN];
    
    actual=edi;
    count = 1;
    while (actual)
    {
        Vfl[count-1]=actual->flood.Vfl;
        actual=actual->next_edi;
        count++;
    }
    if (nnames!=count-1) 
        gmx_fatal(FARGS,"Number of energies is not consistent with t_edi structure");
}
/************* END of FLOODING IMPLEMENTATION ****************************/



static void get_COM(int   nat,   /* number of atoms in the coordinate buffer */
                    rvec  *x,    /* coordinate buffer */
                    real  *m,    /* buffer for the masses */
                    real  tmass, /* total mass */
                    rvec  com)   /* the center of mass */
{
    int  i;
    rvec xm, dum_com = {.0, .0, .0};


    /* calculate COM */
    for (i=0; i<nat; i++)
    {
        svmul(m[i], x[i], xm);
        rvec_inc(dum_com, xm);
    }
    svmul(1.0/tmass, dum_com, dum_com);

    com[XX] = dum_com[XX];
    com[YY] = dum_com[YY];
    com[ZZ] = dum_com[ZZ]; 
}


/* Put current coordinates into origin */
static void subtract_COM(int   nat,  /* number of atoms in the coordinate buffer */
                         rvec  *x,   /* coordinate buffer */
                         rvec  com)  /* the center of mass */
{
    int  i;  

    
    /* subtract COM */
    for (i=0; i<nat; i++) 
        rvec_dec(x[i], com);
}


gmx_edsam_t ed_open(int nfile,t_filenm fnm[],t_commrec *cr)
{   
    gmx_edsam_t ed;
    
    
    /* Allocate space for the ED data structure */
    snew(ed, 1);
    
    /* Set the flavour of ED that we want to do */
    ed->eEDtype = eEDedsam;

    if (MASTER(cr)) 
    {
        /* Open .edi input file: */
        ed->edinam=ftp2fn(efEDI,nfile,fnm);
        /* The master opens the .edo output file */
        fprintf(stderr,"ED sampling will be performed!\n");        
        ed->edonam = ftp2fn(efEDO,nfile,fnm);
        ed->edo    = ffopen(ed->edonam,"w");
    }
    return ed;
}


/* Broadcasts the structure data */
static void bc_ed_positions(t_commrec *cr, struct gmx_edx *s, int stype)
{
    snew_bc(cr, s->anrs, s->nr   );    /* Index numbers     */
    snew_bc(cr, s->x   , s->nr   );    /* Positions         */
    nblock_bc(cr, s->nr, s->anrs );
    nblock_bc(cr, s->nr, s->x    );

    /* For the average & reference structures we need an array for the collective indices,
     * and we need to broadcast the masses as well */
    if (stype == eedAV || stype == eedREF)
    {
        /* We need these additional variables in the parallel case: */
        snew(s->c_ind    , s->nr   );   /* Collective indices */
        snew(s->anrs_loc , s->nr   );   /* Local atom indices */
    }

    /* broadcast masses for the reference structure (for mass-weighted fitting) */
    if (stype == eedREF)
    {
        snew_bc(cr, s->m, s->nr);
        nblock_bc(cr, s->nr, s->m);
    }

    /* For the average structure we might need the masses for mass-weighting */
    if (stype == eedAV)
    {
        snew_bc(cr, s->sqrtm, s->nr);
        nblock_bc(cr, s->nr, s->sqrtm);
        snew_bc(cr, s->m, s->nr);
        nblock_bc(cr, s->nr, s->m);
    }
}


/* Broadcasts the eigenvector data */
static void bc_ed_vecs(t_commrec *cr, t_eigvec *ev, int length)
{
    int i;

    snew_bc(cr, ev->ieig   , ev->neig);  /* index numbers of eigenvector  */
    snew_bc(cr, ev->stpsz  , ev->neig);  /* stepsizes per eigenvector     */
    snew_bc(cr, ev->xproj  , ev->neig);  /* instantaneous x projection    */
    snew_bc(cr, ev->vproj  , ev->neig);  /* instantaneous v projection    */
    snew_bc(cr, ev->fproj  , ev->neig);  /* instantaneous f projection    */
    snew_bc(cr, ev->refproj, ev->neig);  /* starting or target projection */

    nblock_bc(cr, ev->neig, ev->ieig   );
    nblock_bc(cr, ev->neig, ev->stpsz  );
    nblock_bc(cr, ev->neig, ev->xproj  );
    nblock_bc(cr, ev->neig, ev->vproj  );
    nblock_bc(cr, ev->neig, ev->fproj  );
    nblock_bc(cr, ev->neig, ev->refproj);

    snew_bc(cr, ev->vec, ev->neig);      /* Eigenvector components        */
    for (i=0; i<ev->neig; i++)
    {
        snew_bc(cr, ev->vec[i], length);
        nblock_bc(cr, length, ev->vec[i]);
    }
}


/* Broadcasts the ED / flooding data to other nodes, 
 * also allocates memory */
static void broadcast_ed_data(t_commrec *cr, gmx_edsam_t ed, int numedis)
{
    int     i,nr;

    
    if (numedis > 1) 
        gmx_fatal(FARGS, "parallel ED not yet implemented with multiple datasets - sorry");

    /* First let everybody know how many ED data sets to expect */
    gmx_bcast(sizeof(numedis), &numedis, cr);

    /* Now transfer every single ED data set */
    for (i=0; i<numedis; i++)
    {
        /* Broadcast the edpar data structure */
        snew_bc(cr, ed->edpar,1);
        block_bc(cr, *ed->edpar);

        /* Broadcast positions */
        bc_ed_positions(cr, &(ed->edpar->sref), eedREF); /* reference positions (don't broadcast masses)    */
        bc_ed_positions(cr, &(ed->edpar->sav ), eedAV ); /* average positions (do broadcast masses as well) */
        bc_ed_positions(cr, &(ed->edpar->star), eedTAR); /* target positions                                */
        bc_ed_positions(cr, &(ed->edpar->sori), eedORI); /* origin positions                                */

        /* Broadcast eigenvectors */
        bc_ed_vecs(cr, &ed->edpar->vecs.mon   , ed->edpar->sav.nr);
        bc_ed_vecs(cr, &ed->edpar->vecs.linfix, ed->edpar->sav.nr);
        bc_ed_vecs(cr, &ed->edpar->vecs.linacc, ed->edpar->sav.nr);
        bc_ed_vecs(cr, &ed->edpar->vecs.radfix, ed->edpar->sav.nr);
        bc_ed_vecs(cr, &ed->edpar->vecs.radacc, ed->edpar->sav.nr);
        bc_ed_vecs(cr, &ed->edpar->vecs.radcon, ed->edpar->sav.nr);
        /* Broadcast flooding eigenvectors */
        bc_ed_vecs(cr, &ed->edpar->flood.vecs,  ed->edpar->sav.nr);

        /* Initialize local buffer structure, content is NULL anyway */
        snew_bc(cr, ed->edpar->buf, 1);
        /* not needed:  block_bc(cr, *ed->edpar->buf); */
    }
}


/* init-routine called for every *.edi-cycle, initialises t_edpar structure */
static void init_edi(t_topology *top,t_inputrec *ir,
        t_commrec *cr,gmx_edsam_t ed,t_edpar *edi) 
{
    int  i;
    rvec transvec;
    real totalmass = 0.0;
    rvec com;


    /* NOTE Init_edi is executed on the master process only 
     * The initialized data sets are then transmitted to the
     * other nodes in broadcast_ed_data */

    edi->bNeedDoEdsam = edi->vecs.mon.neig 
                     || edi->vecs.linfix.neig 
                     || edi->vecs.linacc.neig 
                     || edi->vecs.radfix.neig 
                     || edi->vecs.radacc.neig 
                     || edi->vecs.radcon.neig;

    /* evaluate masses (reference structure) */
    snew(edi->sref.m, edi->sref.nr);
    for (i = 0; i < edi->sref.nr; i++)
    {
        if (edi->fitmas)
            edi->sref.m[i] = top->atoms.atom[edi->sref.anrs[i]].m;
        else
            edi->sref.m[i] = 1.0;
        totalmass += edi->sref.m[i];
    }
    edi->sref.mtot = totalmass;

    /* Masses m and sqrt(m) for the average structure. Note that m 
     * is needed if forces have to be evaluated in do_edsam */
    snew(edi->sav.sqrtm, edi->sav.nr );
    snew(edi->sav.m    , edi->sav.nr );
    for (i = 0; i < edi->sav.nr; i++)
    {
        edi->sav.m[i] = top->atoms.atom[edi->sav.anrs[i]].m;
        if (edi->pcamas)
            edi->sav.sqrtm[i] = sqrt(top->atoms.atom[edi->sav.anrs[i]].m);
        else
            edi->sav.sqrtm[i] = 1.0;
    }

    /* put reference structure in origin */
    get_COM(edi->sref.nr, edi->sref.x, edi->sref.m, edi->sref.mtot, com);
    subtract_COM(edi->sref.nr, edi->sref.x, com);

    /* Init flooding parameters unless it's a TEE-REX simulation */
    init_flood(edi,ed->edo,ir->delta_t,cr);

    /* Init ED buffer */
    snew(edi->buf, 1);
}


/* Return the atom which is closest to the center of the structure */
static int ed_set_pbcatom(struct gmx_edx *s)
{
    int i, j, i_ind=0, j_ind=0;
    real d, dmax=0.0, dmin;
    int nat;
    rvec center;


    /* Find the two atoms with the largest mutual distance */
    nat = s->nr;
    for (i=0; i<nat; i++)
    {
        for (j=i+1; j<nat; j++)
        {
            d = distance2(s->x[i], s->x[j]);
            /* If we find a larger maximum, save atom pair indices */
            if (d > dmax)
            {
                i_ind = i;
                j_ind = j;
                dmax = d;
            }
        }
    }
    /* Largest distance is between atoms i_ind and j_ind,
     * now find the atom closest to the midpont between i and j */
    rvec_add(s->x[i_ind], s->x[j_ind], center);
    svmul(0.5, center, center);

    d    = distance2(s->x[0], center); /* to begin with */
    dmin = d;
    for (i=1; i<nat; i++)
    {
        d = distance2(s->x[i], center);
        /* Save index of atom closest to center */
        if (d < dmin)
        {
            i_ind = i;
            dmin = d;
        }
    }

    return i_ind;
}


static void check(char *line, char *label)
{
    if (!strstr(line,label)) 
        gmx_fatal(FARGS,"Could not find input parameter %s at expected position in edsam input-file (.edi)\nline read instead is %s",label,line);
}


static int read_checked_edint(FILE *file,char *label)
{
    char line[STRLEN+1];
    int idum;

    
    fgets2 (line,STRLEN,file);
    check(line,label);
    fgets2 (line,STRLEN,file);
    sscanf (line,"%d",&idum);
    return idum;
} 


static int read_edint(FILE *file,bool *bEOF)
{
    char line[STRLEN+1];
    int idum;
    char *eof;

    
    eof=fgets2 (line,STRLEN,file);
    if (eof==NULL)
    {
        *bEOF = TRUE;
        return -1;
    }
    eof=fgets2 (line,STRLEN,file);  
    if (eof==NULL)
    {
        *bEOF = TRUE;
        return -1;
    }
    sscanf (line,"%d",&idum);
    *bEOF = FALSE;
    return idum;
}


static real read_checked_edreal(FILE *file,char *label)
{
    char line[STRLEN+1];
    double rdum;

    
    fgets2 (line,STRLEN,file);
    check(line,label);
    fgets2 (line,STRLEN,file);
    sscanf (line,"%lf",&rdum);
    return (real) rdum; /* always read as double and convert to single */
}


static real read_edreal(FILE *file)
{
    char line[STRLEN+1];
    double rdum;

    
    fgets2 (line,STRLEN,file);
    fgets2 (line,STRLEN,file);
    sscanf (line,"%lf",&rdum);
    return (real) rdum; /* always read as double and convert to single */
}


static int read_edint2(FILE *file)
{
    char line[STRLEN+1];
    int idum;

    
    fgets2 (line,STRLEN,file);
    sscanf (line,"%d",&idum);
    return idum;
}


static void read_edx(FILE *file,int number,int *anrs,rvec *x)
{
    int i,j;
    char line[STRLEN+1];
    double d[3];

    
    for(i=0; (i < number); i++)
    {
        fgets2 (line,STRLEN,file);
        sscanf (line,"%d%lf%lf%lf",&anrs[i],&d[0],&d[1],&d[2]);
        anrs[i]--; /* we are reading FORTRAN indices */
        for(j=0; (j < 3); j++)
            x[i][j]=d[j]; /* always read as double and convert to single */
    }
}


static void scan_edvec(FILE *in,int nr,rvec *vec)
{
    char line[STRLEN+1];
    int i;
    double x,y,z;

    
    for(i=0; (i < nr); i++)
    {
        fgets2 (line,STRLEN,in);
        sscanf (line,"%le%le%le",&x,&y,&z);
        vec[i][XX]=x;
        vec[i][YY]=y;
        vec[i][ZZ]=z;
    }
}


static void read_edvec(FILE *in,int nr,t_eigvec *tvec)
{
    int i,idum;
    double rdum;
    char line[STRLEN+1];

    
    tvec->neig=read_checked_edint(in,"NUMBER OF EIGENVECTORS");
    if (tvec->neig >0)
    {
        snew(tvec->ieig,tvec->neig);
        snew(tvec->stpsz,tvec->neig);
        snew(tvec->vec,tvec->neig);
        snew(tvec->xproj,tvec->neig);
        snew(tvec->vproj,tvec->neig);
        snew(tvec->fproj,tvec->neig);
        snew(tvec->refproj,tvec->neig);
        for(i=0; (i < tvec->neig); i++)
        {
            fgets2 (line,STRLEN,in);
            sscanf (line,"%d%lf",&idum,&rdum);
            tvec->ieig[i]=idum;
            tvec->stpsz[i]=rdum;
        }
        for(i=0; (i < tvec->neig); i++)
        {
            snew(tvec->vec[i],nr);
            scan_edvec(in,nr,tvec->vec[i]);
        }
    }
}


/* calls read_edvec for the vector groups, only for flooding there is an extra call */
static void read_edvecs(FILE *in,int nr,t_edvecs *vecs)
{
    read_edvec(in,nr,&vecs->mon   );
    read_edvec(in,nr,&vecs->linfix);
    read_edvec(in,nr,&vecs->linacc);
    read_edvec(in,nr,&vecs->radfix);
    read_edvec(in,nr,&vecs->radacc);
    read_edvec(in,nr,&vecs->radcon);
}


/* Check if the same atom indices are used for reference and average positions */
static bool check_if_same(struct gmx_edx sref, struct gmx_edx sav)
{
    int i;

    
    /* If the number of atoms differs between the two structures,
     * they cannot be identical */
    if (sref.nr != sav.nr)
        return FALSE;

    /* Now that we know that both stuctures have the same number of atoms,
     * check if also the indices are identical */
    for (i=0; i < sav.nr; i++)
    {
        if (sref.anrs[i] != sav.anrs[i])
            return FALSE;
    }
    fprintf(stderr, "ED: Note: Reference and average structure are composed of the same atom indices.\n");

    return TRUE;
}


static int read_edi(FILE* in, gmx_edsam_t ed,t_edpar *edi,int nr_mdatoms, int edi_nr)
{
    int i,j,idum,readmagic;
    static const int magic=669;
    int ignore;
    rvec *xdum;
    bool bEOF;

    
    /* the edi file is not free format, so expect problems if the input is corrupt. */

    /* check the magic number */
    readmagic=read_edint(in,&bEOF);
    if (bEOF)
        return 0;
    if (readmagic != magic)
    {
        if (readmagic==666 || readmagic==667 || readmagic==668)
            gmx_fatal(FARGS,"wrong magic number: Use newest version of make_edi to produce edi file");
        else
            gmx_fatal(FARGS,"Wrong magic number %d in %s",readmagic,ed->edinam);
    }

    /* check the number of atoms */
    edi->nini=read_edint(in,&bEOF);
    if (edi->nini != nr_mdatoms)
        gmx_fatal(FARGS,"Nr of atoms in %s (%d) does not match nr of md atoms (%d)",
                ed->edinam,edi->nini,nr_mdatoms); 

    /* Done checking. For the rest we blindly trust the input */
    edi->fitmas          = read_checked_edint(in,"FITMAS");
    edi->pcamas          = read_checked_edint(in,"ANALYSIS_MAS");
    edi->outfrq          = read_checked_edint(in,"OUTFRQ");
    edi->maxedsteps      = read_checked_edint(in,"MAXLEN");
    edi->slope           = read_checked_edreal(in,"SLOPECRIT");

    edi->presteps        = read_checked_edint(in,"PRESTEPS");
    edi->flood.deltaF0   = read_checked_edreal(in,"DELTA_F0");
    edi->flood.deltaF    = read_checked_edreal(in,"INIT_DELTA_F");
    edi->flood.tau       = read_checked_edreal(in,"TAU");
    edi->flood.constEfl  = read_checked_edreal(in,"EFL_NULL");
    edi->flood.alpha2    = read_checked_edreal(in,"ALPHA2");
    edi->flood.kT        = read_checked_edreal(in,"KT");
    edi->flood.bHarmonic = read_checked_edint(in,"HARMONIC");
    edi->flood.flood_id  = edi_nr;
    edi->sref.nr         = read_checked_edint(in,"NREF");

    /* allocate space for reference positions and read them */
    snew(edi->sref.anrs,edi->sref.nr);
    snew(edi->sref.x   ,edi->sref.nr);
    edi->sref.sqrtm    =NULL;
    read_edx(in,edi->sref.nr,edi->sref.anrs,edi->sref.x);

    /* average positions. they define which atoms will be used for ED sampling */
    edi->sav.nr=read_checked_edint(in,"NAV");
    snew(edi->sav.anrs,edi->sav.nr);
    snew(edi->sav.x   ,edi->sav.nr);
    read_edx(in,edi->sav.nr,edi->sav.anrs,edi->sav.x);

    /* Check if the same atom indices are used for reference and average positions */
    edi->bRefEqAv = check_if_same(edi->sref, edi->sav);

    edi->ned=edi->sref.anrs[edi->sref.nr-1]+1;
    if (edi->sav.anrs[edi->sav.nr-1] > edi->ned)
        edi->ned=edi->sav.anrs[edi->sav.nr-1]+1;

    /* eigenvectors */
    read_edvecs(in,edi->sav.nr,&edi->vecs);
    read_edvec(in,edi->sav.nr,&edi->flood.vecs);

    /* target positions */
    edi->star.nr=read_edint(in,&bEOF);
    if (edi->star.nr > 0)
    {
        snew(edi->star.anrs,edi->star.nr);
        snew(xdum          ,edi->star.nr);
        edi->star.sqrtm    =NULL;
        read_edx(in,edi->star.nr,edi->star.anrs,xdum);
        snew(edi->star.x,edi->ned);

        for(j=0; (j < edi->star.nr); j++) 
            if (edi->star.anrs[j] < 0 || edi->star.anrs[j] > edi->ned)
                gmx_fatal(FARGS,"ED sampling target index out of bounds: %d\n",edi->star.anrs[j]);
        for(i=0; (i < edi->ned); i++) 
        {
            for(j=0; (j < edi->star.nr); j++) 
            {
                if (edi->star.anrs[j] == i)
                    copy_rvec(xdum[j],edi->star.x[i]);
            }
        }
        sfree(xdum);
    }

    /* positions defining origin of expansion circle */
    edi->sori.nr=read_edint(in,&bEOF);
    if (edi->sori.nr > 0)
    {
        snew(edi->sori.anrs,edi->sori.nr);
        snew(xdum          ,edi->sori.nr);
        edi->sori.sqrtm    =NULL;
        read_edx(in,edi->sori.nr,edi->sori.anrs,xdum);
        snew(edi->sori.x,edi->ned);

        for(j=0; (j < edi->sori.nr); j++) 
            if (edi->sori.anrs[j] < 0 || edi->sori.anrs[j] > edi->ned)
                gmx_fatal(FARGS,"ED sampling origin index out of bounds: %d\n",edi->sori.anrs[j]);
        for(i=0; (i < edi->ned); i++)
        {
            for(j=0; (j < edi->sori.nr); j++)
            {
                if (edi->sori.anrs[j] == i)
                    copy_rvec(xdum[j],edi->sori.x[i]);
            }
        }
        sfree(xdum);
    }
    /* all done */
    return 1;
}


static int read_edi_file(gmx_edsam_t ed, t_edpar *edi, int nr_mdatoms) 
{
    FILE    *in;
    t_edpar *actual_edi;
    t_edpar *edi_read;
    int     edi_nr;
    int     nat_ed;  /* Nr of atoms in essdyn buffer */

    
    /* this routine is called by the master only */

    /* Open the .edi parameter input file */
    in=ffopen(ed->edinam,"r");  
    fprintf(stderr, "ED: Opening edi file %s\n", ed->edinam);

    /* Now read a sequence of ED input parameter sets from the edi file */
    actual_edi=edi;
    edi_nr=0;
    read_edi(in,ed,actual_edi,nr_mdatoms,edi_nr++);
    nat_ed = actual_edi->ned;
    if (edi->nini!=nr_mdatoms)
        gmx_fatal(FARGS,"edi file %s was made for %d atoms, but the simulation contains %d atoms.", 
                ed->edinam, edi->nini, nr_mdatoms);
    snew(edi_read,1);
    while( read_edi(in, ed, edi_read, nr_mdatoms, edi_nr++))
    {
        actual_edi->next_edi=edi_read;
        actual_edi=edi_read;
        snew(edi_read,1);
    }
    sfree(edi_read);
    actual_edi->next_edi=NULL;

    /* Close the .edi file again */
    ffclose(in);

    return nat_ed;
}


/* Fit the current coordinates to the reference coordinates 
 * Do not actually do the fit, just return rotation and translation.
 * Note that the COM of the reference structure was already put into 
 * the origin by init_edi. New version of fitit */
static void fit_to_reference(rvec      *xcoll,    /* The coordinates to be fitted */
                             rvec      transvec,  /* The translation vector */ 
                             matrix    rotmat,    /* The rotation matrix */
                             t_edpar   *edi,      /* Just needed for do_edfit */
                             t_commrec *cr)
{
    static rvec *xcopy=NULL;  /* Working copy of the coordinates */
    static rvec com;          /* center of mass */
    int         i;


    GMX_MPE_LOG(ev_fit_to_reference_start);

    /* We do not touch the original coordinates but work on a copy */
    if (!xcopy)
        snew(xcopy, edi->sref.nr);

    for (i=0; i<edi->sref.nr; i++)
        copy_rvec(xcoll[i], xcopy[i]);

    /* Calculate the center of mass */
    get_COM(edi->sref.nr, xcopy, edi->sref.m, edi->sref.mtot, com);

    /* Subtract the center of mass from the copy */
    subtract_COM(edi->sref.nr, xcopy, com);

    /* Determine the rotation matrix */
    do_edfit(edi->sref.nr, edi->sref.x, xcopy, rotmat, edi);

    transvec[XX] = -com[XX];
    transvec[YY] = -com[YY];
    transvec[ZZ] = -com[ZZ];

    GMX_MPE_LOG(ev_fit_to_reference_finish);
}


static void translate_and_rotate(rvec *x,         /* The coordinates to be translated and rotated */
                                 int nat,         /* How many coordinates are there */
                                 rvec transvec,   /* The translation vector */ 
                                 matrix rotmat)   /* The rotation matrix */
{
    int i;

    
    /* Translation */
    for (i=0; i<nat; i++)
        rvec_inc(x[i], transvec);

    /* Rotation */
    rotate_x(nat, x, rotmat);
}


/* Gets the rms deviation of the x coordinates to the structure s */
/* fit_to_structure has to be called before calling this routine! */
static real rmsd_from_structure(rvec           *x,  /* The x coordinates under consideration */ 
                                struct gmx_edx *s)  /* The structure from which the rmsd shall be computed */
{
    real  rmsd=0.0;
    int   i;


    for (i=0; i < s->nr; i++)
        rmsd += distance2(s->x[i], x[i]);

    rmsd /= (real) s->nr;
    rmsd = sqrt(rmsd);

    return rmsd;
}


static real get_rmsd(t_edpar *edi, rvec *xcoll) 
{
    /* fit has to be done previously */
    real rmsd = 0.0;
    int  i;


    for (i=0; i<edi->sref.nr; i++)
        rmsd += distance2(edi->sref.x[i], xcoll[i]);

    rmsd /= (real) edi->sref.nr;
    rmsd = sqrt(rmsd);

    return rmsd;
}


/* select the indices of the ED atoms which are local 
 *
 * Only the indices that correspond to the structure s are
 * taken into account and saved in s->c_ind[] 
 */
static void dd_make_local_indices(gmx_domdec_t *dd, struct gmx_edx *s, t_mdatoms *md)
{
    int         i,ii;
    gmx_ga2la_t *ga2la=NULL;


    ga2la = dd->ga2la;

    /* we have not yet found a local atom */
    s->nr_loc = 0; 
    /* go through all the atom indices of the structure */
    for(i=0; i<s->nr; i++)
    {
        if (ga2la[s->anrs[i]].cell == 0)
        {
            ii = ga2la[s->anrs[i]].a;
            if (ii < md->start+md->homenr)
            {
                /* The atom with this index is a home atom, therefore 
                 * save its local index in local atom numbers array */
                s->anrs_loc[s->nr_loc] = ii;
                /* keep track of where this local atom is in the collective c_ind array: */
                s->c_ind[s->nr_loc] = i;
                /* add one to the local atom count */
                s->nr_loc++; 
            }
        }
    }
}


void dd_make_local_ed_indices(gmx_domdec_t *dd, struct gmx_edsam *ed,t_mdatoms *md)
{
    if (ed->eEDtype != eEDnone)
    {
        /* Local atoms of the reference structure (for fitting) */
        dd_make_local_indices(dd, &ed->edpar->sref, md);
        /* Local atoms of the average structure (on these ED will be performed) */
        dd_make_local_indices(dd, &ed->edpar->sav , md);
    }
}


static void remove_pbc_effect(int ePBC,rvec transvec_compact, matrix box,gmx_edsam_t ed) 
{
    bool bFirst;
    rvec null = {0.0,0.0,0.0};


    if(ed->pbc)
        bFirst = FALSE;
    else
    {
        bFirst = TRUE;
        snew(ed->pbc,1);
    }

    if (bFirst)
        set_pbc(ed->pbc,ePBC,box);

    pbc_dx(ed->pbc,null,transvec_compact,transvec_compact);
}


static void ed_get_shifts(int npbcdim,matrix box,
        rvec *xc, const rvec reference, ivec *shifts, int nat, t_commrec *cr)
{
    int  i,m,d;
    rvec dx;

    
    /* Get the shifts such that each atom is within closest
     * distance to the reference atom after shifting */
    for (i=0; i<nat; i++)
        clear_ivec(shifts[i]);

    for (i=0; i<nat; i++)
    {
        rvec_sub(xc[i],reference,dx); /* distance atom - reference */

        for(m=npbcdim-1; m>=0; m--)
        {
            while (dx[m] < -0.5*box[m][m])
            {
                for(d=0; d<DIM; d++)
                    dx[d] += box[m][d];
                shifts[i][m]++;
            }
            while (dx[m] >= 0.5*box[m][m])
            {
                for(d=0; d<DIM; d++)
                    dx[d] -= box[m][d];
                shifts[i][m]--;
            }
        }
    }
}


static void ed_shift_coords(matrix box, rvec x[], ivec *is, int nr)
{
    int      i,tx,ty,tz;


    GMX_MPE_LOG(ev_shift_start);

    /* Loop over the ED atoms */
    if(TRICLINIC(box)) 
    {
        for (i=0; i < nr; i++)
        {
            tx=is[i][XX];
            ty=is[i][YY];
            tz=is[i][ZZ];

            x[i][XX]=x[i][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
            x[i][YY]=x[i][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
            x[i][ZZ]=x[i][ZZ]+tz*box[ZZ][ZZ];
        }
    } else
    {
        for (i=0; i < nr; i++)
        {
            tx=is[i][XX];
            ty=is[i][YY];
            tz=is[i][ZZ];

            x[i][XX]=x[i][XX]+tx*box[XX][XX];
            x[i][YY]=x[i][YY]+ty*box[YY][YY];
            x[i][ZZ]=x[i][ZZ]+tz*box[ZZ][ZZ];
        }
    }    
    GMX_MPE_LOG(ev_shift_finish);
}


static void ed_unshift_coords(matrix box, rvec x[], ivec *is, int nr)
{
    int i,tx,ty,tz;


    GMX_MPE_LOG(ev_unshift_start);

    if(TRICLINIC(box))
    {
        for(i=0; i < nr; i++)
        {
            tx=is[i][XX];
            ty=is[i][YY];
            tz=is[i][ZZ];

            x[i][XX]=x[i][XX]-tx*box[XX][XX]-ty*box[YY][XX]-tz*box[ZZ][XX];
            x[i][YY]=x[i][YY]-ty*box[YY][YY]-tz*box[ZZ][YY];
            x[i][ZZ]=x[i][ZZ]-tz*box[ZZ][ZZ];
        }
    } else
    {
        for(i=0; i < nr; i++)
        {
            tx=is[i][XX];
            ty=is[i][YY];
            tz=is[i][ZZ];

            x[i][XX]=x[i][XX]-tx*box[XX][XX];
            x[i][YY]=x[i][YY]-ty*box[YY][YY];
            x[i][ZZ]=x[i][ZZ]-tz*box[ZZ][ZZ];
        }
    }

    GMX_MPE_LOG(ev_unshift_finish);
}


/* Assemble the coordinates such that every node has all of them. 
 * Get the indices from structure s */
static void get_coordinates(t_commrec      *cr, 
                            rvec           *xc,         /* Collective array of coordinates (write here) */
                            ivec           *shifts_xc,  /* Collective array of shifts */
                            rvec           *x_loc,      /* Local coordinates on this node (read coords from here) */ 
                            struct gmx_edx *s, 
                            matrix         box,
                            char           title[])
{
    int i;
    rvec reference_x;


    GMX_MPE_LOG(ev_get_coords_start);

    /* Zero out the collective coordinate array */
    clear_rvecs(s->nr, xc);

    /* Put the local coordinates that this node has into the right place of 
     * the collective array. Note that in the serial case, s->c_ind[i] = i */
    for (i=0; i<s->nr_loc; i++)
        copy_rvec(x_loc[s->anrs_loc[i]], xc[s->c_ind[i]]);

    if (PAR(cr))
        /* Add the arrays from all nodes together */
        gmx_sum(s->nr*3, xc[0], cr);

    /* We now need to move the assembled coordinates within closest distance 
     * to the reference atom */
    /* Get the current pbc reference coordinate */
    copy_rvec(xc[s->c_pbcatom], reference_x);

    /* Choose periodic images closest to pbcatom */
    ed_get_shifts(3, box, xc, reference_x, shifts_xc, s->nr, cr);
    ed_shift_coords(box, xc, shifts_xc, s->nr);

    GMX_MPE_LOG(ev_get_coords_finish);
}


static void rotate_vec(int nr, rvec *x, matrix rotmat)
{
    int i,j,k;
    rvec xdum;


    /* apply the rotation matrix */
    for (i=0; i<nr; i++)
    {
        for (j=0; j<3; j++)
            xdum[j]=x[i][j];
        for (j=0; j<3; j++)
        {
            x[i][j]=0;
            for (k=0; k<3; k++)
                x[i][j] += rotmat[k][j]*xdum[k];
        }
    }
}



static void do_linfix(rvec *xcoll, t_edpar *edi, int step, t_commrec *cr)
{
    int  i, j;
    real proj, add;
    rvec vec_dum;


    /* loop over linfix vectors */
    for (i=0; i<edi->vecs.linfix.neig; i++)
    {
        /* calculate the projection */
        proj = projectx(edi, xcoll, edi->vecs.linfix.vec[i], cr);

        /* calculate the correction */
        add = edi->vecs.linfix.refproj[i] + step*edi->vecs.linfix.stpsz[i] - proj;

        /* apply the correction */
        add /= edi->sav.sqrtm[i];
        for (j=0; j<edi->sav.nr; j++)
        {
            svmul(add, edi->vecs.linfix.vec[i][j], vec_dum);
            rvec_inc(xcoll[j], vec_dum);
        }
    }
}


static void do_linacc(rvec *xcoll, t_edpar *edi, t_commrec *cr)
{
    int  i, j;
    real proj, add;
    rvec vec_dum;


    /* loop over linacc vectors */
    for (i=0; i<edi->vecs.linacc.neig; i++)
    {
        /* calculate the projection */
        proj=projectx(edi, xcoll, edi->vecs.linacc.vec[i], cr);

        /* calculate the correction */
        add = 0.0;
        if (edi->vecs.linacc.stpsz[i] > 0.0)
        {
            if ((proj-edi->vecs.linacc.refproj[i]) < 0.0)
                add = edi->vecs.linacc.refproj[i] - proj;
        }
        if (edi->vecs.linacc.stpsz[i] < 0.0)
        {
            if ((proj-edi->vecs.linacc.refproj[i]) > 0.0)
                add = edi->vecs.linacc.refproj[i] - proj;
        }

        /* apply the correction */
        add /= edi->sav.sqrtm[i];
        for (j=0; j<edi->sav.nr; j++)
        {
            svmul(add, edi->vecs.linacc.vec[i][j], vec_dum);
            rvec_inc(xcoll[j], vec_dum);
        }

        /* new positions will act as reference */
        edi->vecs.linacc.refproj[i] = proj + add;
    }
}


static void do_radfix(rvec *xcoll, t_edpar *edi, int step, t_commrec *cr)
{
    int  i,j;
    real *proj, rad=0.0, ratio;
    rvec vec_dum;


    if (edi->vecs.radfix.neig == 0) 
        return;

    snew(proj, edi->vecs.radfix.neig);

    /* loop over radfix vectors */
    for (i=0; i<edi->vecs.radfix.neig; i++)
    {
        /* calculate the projections, radius */
        proj[i] = projectx(edi, xcoll, edi->vecs.radfix.vec[i], cr);
        rad += pow(proj[i] - edi->vecs.radfix.refproj[i], 2);
    }

    rad   = sqrt(rad);
    ratio = (edi->vecs.radfix.stpsz[0]+edi->vecs.radfix.radius)/rad - 1.0;
    edi->vecs.radfix.radius += edi->vecs.radfix.stpsz[0];

    /* loop over radfix vectors */
    for (i=0; i<edi->vecs.radfix.neig; i++)
    {
        proj[i] -= edi->vecs.radfix.refproj[i];

        /* apply the correction */
        proj[i] /= edi->sav.sqrtm[i];
        proj[i] *= ratio;
        for (j=0; j<edi->sav.nr; j++) {
            svmul(proj[i], edi->vecs.radfix.vec[i][j], vec_dum);
            rvec_inc(xcoll[j], vec_dum);
        }
    }

    sfree(proj);
}


static void do_radacc(rvec *xcoll, t_edpar *edi, t_commrec *cr)
{
    int  i,j;
    real *proj, rad=0.0, ratio=0.0;
    rvec vec_dum;


    if (edi->vecs.radacc.neig == 0) 
        return;

    snew(proj,edi->vecs.radacc.neig);

    /* loop over radacc vectors */
    for (i=0; i<edi->vecs.radacc.neig; i++)
    {
        /* calculate the projections, radius */
        proj[i] = projectx(edi, xcoll, edi->vecs.radacc.vec[i], cr);
        rad += pow(proj[i] - edi->vecs.radacc.refproj[i], 2);
    }
    rad = sqrt(rad);

    /* only correct when radius decreased */
    if (rad < edi->vecs.radacc.radius)
    {
        ratio = edi->vecs.radacc.radius/rad - 1.0;
        rad   = edi->vecs.radacc.radius;
    }
    else 
        edi->vecs.radacc.radius = rad;

    /* loop over radacc vectors */
    for (i=0; i<edi->vecs.radacc.neig; i++)
    {
        proj[i] -= edi->vecs.radacc.refproj[i];

        /* apply the correction */
        proj[i] /= edi->sav.sqrtm[i];
        proj[i] *= ratio;
        for (j=0; j<edi->sav.nr; j++) {
            svmul(proj[i], edi->vecs.radacc.vec[i][j], vec_dum);
            rvec_inc(xcoll[j], vec_dum);
        }
    }  
    sfree(proj);
}


struct t_do_radcon {
    real *proj;
};

static void do_radcon(rvec *xcoll, t_edpar *edi, t_commrec *cr)
{
    int  i,j,k;
    real rad=0.0, ratio=0.0;
    struct t_do_radcon *loc;
    bool bFirst;
    rvec vec_dum;


    if(edi->buf->do_radcon != NULL)
    {
        bFirst = FALSE;
        loc    = edi->buf->do_radcon;
    } 
    else
    {
        bFirst = TRUE;
        snew(edi->buf->do_radcon, 1);
    }
    loc = edi->buf->do_radcon;

    if (edi->vecs.radcon.neig == 0) 
        return;

    if (bFirst)
        snew(loc->proj, edi->vecs.radcon.neig);

    /* loop over radcon vectors */
    for (i=0; i<edi->vecs.radcon.neig; i++)
    {
        /* calculate the projections, radius */
        loc->proj[i] = projectx(edi, xcoll, edi->vecs.radcon.vec[i], cr);
        rad += pow(loc->proj[i] - edi->vecs.radcon.refproj[i], 2);
    }
    rad = sqrt(rad);
    /* only correct when radius increased */
    if (rad > edi->vecs.radcon.radius)
    {
        ratio = edi->vecs.radcon.radius/rad - 1.0;
    } 
    else
        edi->vecs.radcon.radius = rad;

    /* loop over radcon vectors */
    for (i=0; i<edi->vecs.radcon.neig; i++)
    {
        /* apply the correction */
        loc->proj[i] -= edi->vecs.radcon.refproj[i];
        loc->proj[i] /= edi->sav.sqrtm[i];
        loc->proj[i] *= ratio;
        for (j=0; j<edi->sav.nr; j++)
        {
            svmul(loc->proj[i], edi->vecs.radcon.vec[i][j], vec_dum);
            rvec_inc(xcoll[i], vec_dum);
        }
    }

    if (rad != edi->vecs.radcon.radius)
    {
        rad = 0.0;
        for (i=0; i<edi->vecs.radcon.neig; i++)
        {
            /* calculate the projections, radius */
            loc->proj[i] = projectx(edi, xcoll, edi->vecs.radcon.vec[i], cr);
            rad += pow(loc->proj[i] - edi->vecs.radcon.refproj[i], 2);
        }
        rad = sqrt(rad);
    }
}


static void ed_apply_constraints(rvec *xcoll, t_edpar *edi, int step, t_commrec *cr)
{
    int i;


    GMX_MPE_LOG(ev_ed_apply_cons_start);

    /* subtract the average positions */
    for (i=0; i<edi->sav.nr; i++)
        rvec_dec(xcoll[i], edi->sav.x[i]);

    /* apply the constraints */
    if (step >= 0) 
        do_linfix(xcoll, edi, step, cr);
    do_linacc(xcoll, edi, cr);
    if (step >= 0) 
        do_radfix(xcoll, edi, step, cr);
    do_radacc(xcoll, edi, cr);
    do_radcon(xcoll, edi, cr);

    /* add back the average positions */
    for (i=0; i<edi->sav.nr; i++) 
        rvec_inc(xcoll[i], edi->sav.x[i]);

    GMX_MPE_LOG(ev_ed_apply_cons_finish);
}


static void write_edo(t_edpar *edi, gmx_edsam_t ed, int step,real rmsd)
{  
    int i;


    if (step == -1)
        fprintf(ed->edo, "Initial projections:\n");
    else
    {
        fprintf(ed->edo,"Step %d",step);
        fprintf(ed->edo,"  RMSD %f nm\n",rmsd);
        if (ed->eEDtype == eEDflood)
          fprintf(ed->edo, "  Efl=%f  deltaF=%f  Vfl=%f\n",edi->flood.Efl,edi->flood.deltaF,edi->flood.Vfl);
    }

    if (edi->vecs.mon.neig)
    {
        fprintf(ed->edo,"Monitor eigenvectors");
        for (i=0; i<edi->vecs.mon.neig; i++)
            fprintf(ed->edo," %d: %12.6e ",edi->vecs.mon.ieig[i],edi->vecs.mon.xproj[i]);
        fprintf(ed->edo,"\n");
    }
    if (edi->vecs.linfix.neig)
    {
        fprintf(ed->edo,"Linfix  eigenvectors");
        for (i=0; i<edi->vecs.linfix.neig; i++)
            fprintf(ed->edo," %d: %12.6e ",edi->vecs.linfix.ieig[i],edi->vecs.linfix.xproj[i]);
        fprintf(ed->edo,"\n");
    }
    if (edi->vecs.linacc.neig)
    {
        fprintf(ed->edo,"Linacc  eigenvectors");
        for (i=0; i<edi->vecs.linacc.neig; i++)
            fprintf(ed->edo," %d: %12.6e ",edi->vecs.linacc.ieig[i],edi->vecs.linacc.xproj[i]);
        fprintf(ed->edo,"\n");
    }
    if (edi->vecs.radfix.neig)
    {
        fprintf(ed->edo,"Radfix  eigenvectors");
        for (i=0; i<edi->vecs.radfix.neig; i++)
            fprintf(ed->edo," %d: %12.6e ",edi->vecs.radfix.ieig[i],edi->vecs.radfix.xproj[i]);
        fprintf(ed->edo,"\n");
        fprintf(ed->edo,"fixed increment radius = %f\n", calc_radius(&edi->vecs.radfix));
    }
    if (edi->vecs.radacc.neig)
    {
        fprintf(ed->edo,"Radacc  eigenvectors");
        for (i=0; i<edi->vecs.radacc.neig; i++)
            fprintf(ed->edo," %d: %12.6e ",edi->vecs.radacc.ieig[i],edi->vecs.radacc.xproj[i]);
        fprintf(ed->edo,"\n");
        fprintf(ed->edo,"acceptance radius      = %f\n", calc_radius(&edi->vecs.radacc));
    }
    if (edi->vecs.radcon.neig)
    {
        fprintf(ed->edo,"Radcon  eigenvectors");
        for (i=0; i<edi->vecs.radcon.neig; i++)
            fprintf(ed->edo," %d: %12.6e ",edi->vecs.radcon.ieig[i],edi->vecs.radcon.xproj[i]);
        fprintf(ed->edo,"\n");
        fprintf(ed->edo,"contracting radius     = %f\n", calc_radius(&edi->vecs.radcon));
    }
}


static int ed_constraints(gmx_edsam_t ed)
{ 
    /* returns if any constraints are switched on */
    t_edpar *edi;

    if (ed->eEDtype == eEDedsam) 
    {
        edi=ed->edpar;
        return (edi->vecs.linfix.neig || edi->vecs.linacc.neig || 
                edi->vecs.radfix.neig || edi->vecs.radacc.neig ||  
                edi->vecs.radcon.neig);
    } 
    return 0;
}


void init_edsam(t_topology  *top,    /* global topology                    */
                t_inputrec  *ir,     /* input record                       */
                t_commrec   *cr,     /* communication record               */
                gmx_edsam_t ed,      /* contains all ED data               */
                rvec        x[],     /* coordinates of the whole MD system */
                matrix      box)     /* the box                            */
{
    t_edpar *edi;          /* points to a single edi data set */
    int     numedis=0;     /* keep track of the number of ED data sets in edi file */
    int     i,j;
    int     nat_ed=-1;     /* number of atoms in ED buffer */
    rvec    *x_pbc;        /* coordinates of the whole MD system with pbc removed  */
    rvec    *transvec;
    matrix  rotmat;
    rvec    *xfit;         /* the coordinates which will be fitted to the reference structure  */
    rvec    *xstart;       /* the coordinates which are subject to ED sampling */
    rvec    fit_transvec;  /* translation ... */
    matrix  fit_rotmat;    /* ... and rotation from fit to reference structure */
    char    fn_graph[255]; /* Filename for outputting graph information */

static bool bDebugWait = 0;

while (bDebugWait)
    ;


    if (PARTDECOMP(cr) && PAR(cr))
        gmx_fatal(FARGS, "Please switch on domain decomposition to use essential dynamics in parallel.");

    GMX_MPE_LOG(ev_edsam_start);

    if (MASTER(cr))
        fprintf(stderr, "ED: Initialzing essential dynamics constraints.\n");

    ed->ePBC = ir->ePBC;

    /* The input file is read by the master and the edi structures are
     * initialized here. Input is stored in ed->edpar. Then the edi
     * structures are transferred to the other nodes */
    if (MASTER(cr))
    {
        snew(ed->edpar,1);
        nat_ed = read_edi_file(ed,ed->edpar,top->atoms.nr);

        /* Initialization for every ED/flooding dataset */
        /* Flooding can use multiple edi datasets */
        edi=ed->edpar;
        while(edi != NULL)
        {
            init_edi(top,ir,cr,ed,edi);
            edi=edi->next_edi;
            numedis++;
        }
    }

    /* The master does the work here. The other nodes get the coordinates 
     * not before dd_partition_system which is called after init_edsam */
    if (MASTER(cr))
    {
        /* Reset pointer to first ED data set which contains the actual ED data */
        edi=ed->edpar;

        /* Remove pbc, make molecule whole */
        snew(x_pbc,top->atoms.nr);
        m_rveccopy(top->atoms.nr,x,x_pbc);
        rm_pbc(&(top->idef),ir->ePBC,top->atoms.nr,box,x_pbc,x_pbc);
        /* Extract the coordinates of the atoms to which will be fitted */
        snew(xfit, edi->sref.nr);
        for (i=0; i < edi->sref.nr; i++) {
            copy_rvec(x_pbc[edi->sref.anrs[i]], xfit[i]);
        }
        /* Extract the coordinates of the atoms subject to ED sampling */
        snew(xstart, edi->sav.nr);
        for (i=0; i < edi->sav.nr; i++)
            copy_rvec(x_pbc[edi->sav.anrs[i]], xstart[i]);

        /* Make the fit to the REFERENCE structure, get translation and rotation */
        fit_to_reference(xfit, fit_transvec, fit_rotmat, edi, cr);
        
        /* Output how well we fit to the reference at the start */
        translate_and_rotate(xfit, edi->sref.nr, fit_transvec, fit_rotmat);        
        fprintf(stderr, "ED: Initial RMSD from reference after fit = %f nm\n",
                rmsd_from_structure(xfit, &edi->sref));
        sfree(xfit);

        /* Now apply the translation and rotation to the atoms on which ED sampling will be performed */
        translate_and_rotate(xstart, edi->sav.nr, fit_transvec, fit_rotmat);
        
        /* calculate initial projections */
        project_p(xstart, edi, "x", cr);

        /* process target structure, if required */
        if (edi->star.nr > 0)
        {
            /* get translation & rotation for fit of target structure to reference structure */
            fit_to_reference(edi->star.x, fit_transvec, fit_rotmat, edi, cr);
            /* do the fit */
            translate_and_rotate(edi->star.x, edi->sav.nr, fit_transvec, fit_rotmat);
            rad_project(edi, edi->star.x, &edi->vecs.radcon, cr);
        } else
            rad_project(edi, xstart, &edi->vecs.radcon, cr);

        /* process structure that will serve as origin of expansion circle */
        if (edi->sori.nr > 0)
        {
            /* fit this structure to reference structure */
            fit_to_reference(edi->sori.x, fit_transvec, fit_rotmat, edi, cr);
            /* do the fit */
            translate_and_rotate(edi->sori.x, edi->sav.nr, fit_transvec, fit_rotmat);
            rad_project(edi, edi->sori.x, &edi->vecs.radacc, cr);
            rad_project(edi, edi->sori.x, &edi->vecs.radfix, cr);
        }
        else
        {
            rad_project(edi, xstart, &edi->vecs.radacc, cr);
            rad_project(edi, xstart, &edi->vecs.radfix, cr);
        }

        /* set starting projections for linsam */
        rad_project(edi, xstart, &edi->vecs.linacc, cr);
        rad_project(edi, xstart, &edi->vecs.linfix, cr);

        sfree(xstart);

        /* Define a reference atom of the ED structure */
        edi->sav.c_pbcatom  = ed_set_pbcatom(&(edi->sav ));
        edi->sref.c_pbcatom = ed_set_pbcatom(&(edi->sref));
        
        /* Output to file */
        if (ed->edo)
        {
            fprintf(ed->edo, "pbc atom of average structure: %d. pbc atom of reference structure: %d\n", 
                    edi->sav.c_pbcatom, edi->sref.c_pbcatom);
            /* Set the step to -1 so that write_edo knows it was called from init_edsam */
            write_edo(edi, ed, -1, 0);
        }
    } /* end of MASTER only section */

    /* Broadcast the essential dynamics / flooding data to all nodes */
    if (PAR(cr))
        broadcast_ed_data(cr, ed, numedis);

    edi=ed->edpar;

    if (!PAR(cr))
    {
        /* Point the local atom numbers pointers to the global one, so
         * that we can use the same notation in serial and parallel case: */
        edi=ed->edpar;
        edi->sref.anrs_loc = edi->sref.anrs;
        edi->sav.anrs_loc  = edi->sav.anrs;
        edi->star.anrs_loc = edi->star.anrs;
        edi->sori.anrs_loc = edi->sori.anrs;
        /* For the same reason as above, make a dummy c_ind array: */
        snew(edi->sav.c_ind, edi->sav.nr);
        /* Initialize the array */
        for (i=0; i<edi->sav.nr; i++)
            edi->sav.c_ind[i] = i;
        /* In the general case we will need a different-sized array for the reference indices: */
        if (!edi->bRefEqAv)
        {
            snew(edi->sref.c_ind, edi->sref.nr);
            for (i=0; i<edi->sref.nr; i++)
                edi->sref.c_ind[i] = i;
        }
        /* Point to the very same array in case of other structures: */
        edi->star.c_ind = edi->sav.c_ind;
        edi->sori.c_ind = edi->sav.c_ind;
        /* In the serial case, the local number of atoms is the global one: */
        edi->sref.nr_loc = edi->sref.nr;
        edi->sav.nr_loc  = edi->sav.nr;
        edi->star.nr_loc = edi->star.nr;
        edi->sori.nr_loc = edi->sori.nr;
    }

#ifdef DUMPEDI
    /* Dump it all into one file per process */
    dump_edi(ed->edpar, cr);
#endif

    /* Allocate space for ED buffer */
    snew(edi->buf, 1);
    snew(edi->buf->do_edsam, 1);

    /* Allocate space for collective ED buffer variables */

    /* Collective coordinates of atoms with the average indices */
    snew(edi->buf->do_edsam->xcoll            , edi->sav.nr);
    snew(edi->buf->do_edsam->shifts_xcoll     , edi->sav.nr); /* buffer for xcoll shifts */ 
    /* Collective coordinates of atoms with the reference indices */
    if (!edi->bRefEqAv)
    {
        snew(edi->buf->do_edsam->xc_ref       , edi->sref.nr);
        snew(edi->buf->do_edsam->shifts_xc_ref, edi->sref.nr); /* To store the shifts in */
    }
 
    /* Flush the edo file so that the user can check some things 
     * when the simulation has started */
    if (ed->edo)
        fflush(ed->edo);

    GMX_MPE_LOG(ev_edsam_finish);
}


void do_edsam(t_topology  *top,   /* Local topology  */
              t_inputrec  *ir,
              int         step,
              t_mdatoms   *md,
              t_commrec   *cr,
              rvec        xs[],   /* The local current coordinates on this processor */
              matrix      box,
              gmx_edsam_t ed)
{
    int     i,j,iupdate=500;
    matrix  rotmat;         /* rotation matrix */
    rvec    transvec;       /* translation vector */
    rvec    dx, null={.0, .0, .0}; /* distance */
    static bool  bFirst=TRUE;
    real    dt,dt_1,dt_2;
    struct t_do_edsam *buf;
    t_edpar *edi;
    real    rmsdev;        /* RMSD from reference structure prior to applying the constraints */


#ifdef DEBUGPRINT
    static FILE *fdebug;
    char fname[255];

    if (bFirst)
    {
        sprintf(fname, "debug%2.2d", cr->nodeid);
        fdebug = fopen(fname, "w");
    }
#endif

    /* Check if ED sampling has to be performed */
    if ( ed->eEDtype==eEDnone )
        return;
    
    edi=ed->edpar;  
    if (!edi->bNeedDoEdsam)
        return;

    /* Initialise some variables */
    dt   = ir->delta_t;
    dt_1 = 1.0/dt;
    dt_2 = 1.0/(dt*dt);
    buf=edi->buf->do_edsam;

    if (bFirst)
        /* initialise radacc radius for slope criterion */
        buf->oldrad=calc_radius(&edi->vecs.radacc);

    /* Copy the coordinates into buf->xc* arrays and after ED 
     * feed back corrections to the official coordinates */

    /* Broadcast the ED coordinates such that every node has all of them
     * Every node contributes its local coordinates xs and stores it in
     * the collective buf->xcoll array. */  
    get_coordinates(cr, buf->xcoll, buf->shifts_xcoll, xs, &edi->sav,  box, "XC_AVERAGE");  

    /* Only assembly reference coordinates if their indices differ from the average ones */
    if (!edi->bRefEqAv)
        get_coordinates(cr, buf->xc_ref, buf->shifts_xc_ref, xs, &edi->sref, box, "XC_REFERENCE");

    /* Now all nodes have all of the ED coordinates in edi->sav->xcoll,
     * as well as the indices in edi->sav.anrs */

    /* Fit the reference indices to the reference structure */
    if (edi->bRefEqAv)
        fit_to_reference(buf->xcoll , transvec, rotmat, edi, cr);
    else
        fit_to_reference(buf->xc_ref, transvec, rotmat, edi, cr);

    /* Now apply the translation and rotation to the ED structure */
    translate_and_rotate(buf->xcoll, edi->sav.nr, transvec, rotmat);

    /* Find out how well we fit to the reference (just for output steps) */
    if (do_per_step(step,edi->outfrq) && MASTER(cr))
    {
        if (edi->bRefEqAv)
        {
            /* Indices of reference and average structures are identical, 
             * thus we can calculate the rmsd to SREF using xcoll */
            rmsdev = rmsd_from_structure(buf->xcoll,&edi->sref);
        }
        else
        {
            /* We have to translate & rotate the reference atoms first */
            translate_and_rotate(buf->xc_ref, edi->sref.nr, transvec, rotmat);
            rmsdev = rmsd_from_structure(buf->xc_ref,&edi->sref);
        }
    }
    
    copy_rvec(transvec, buf->transvec_compact);

    /* this is to remove virtual jumps in translational velocity due to periodic boundaries */
    remove_pbc_effect(ed->ePBC,buf->transvec_compact,box,ed);

    if (bFirst)
    {
        copy_mat(rotmat, buf->old_rotmat);
        copy_rvec(buf->transvec_compact, buf->old_transvec  );
        copy_rvec(buf->old_transvec    , buf->older_transvec);
        dt_1 = 0; 
        dt_2 = 0;
    }

    /* update radsam references, when required */
    if (do_per_step(step,edi->maxedsteps) && step > edi->presteps)
    {
        project_p(buf->xcoll, edi, "x", cr);
        rad_project(edi, buf->xcoll, &edi->vecs.radacc, cr);
        rad_project(edi, buf->xcoll, &edi->vecs.radfix, cr);
        buf->oldrad=-1.e5;
    }

    /* update radacc references, when required */
    if (do_per_step(step,iupdate) && step > edi->presteps)
    {
        edi->vecs.radacc.radius = calc_radius(&edi->vecs.radacc);
        if (edi->vecs.radacc.radius - buf->oldrad < edi->slope)
        {
            project_p(buf->xcoll, edi, "x", cr);
            rad_project(edi, buf->xcoll, &edi->vecs.radacc, cr);
            buf->oldrad = 0.0;
        } else
            buf->oldrad = edi->vecs.radacc.radius;
    }

    /* apply the constraints */
    if (step > edi->presteps && ed_constraints(ed))
        ed_apply_constraints(buf->xcoll, edi, step, cr);

    /* write to edo, when required */
    if (do_per_step(step,edi->outfrq))
    {
        project_p(buf->xcoll, edi, "x", cr);
        if (MASTER(cr))
          write_edo(edi, ed, step, rmsdev);
    }

    /* remove fitting */
    rmfit(edi->sav.nr, buf->xcoll, transvec, rotmat);

    /* Unshift ED coordinates */
    ed_unshift_coords(box, buf->xcoll, buf->shifts_xcoll, edi->sav.nr);

    /* Copy the ED corrected coordinates into the coordinate array */
    /* Each node copies his local part. In the serial case, nat_loc is the 
     * total number of ED atoms */

    /* Copy back the coordinates unless monitoring only */
    if (ed_constraints(ed))
        for (i=0; i<edi->sav.nr_loc; i++)
            copy_rvec(buf->xcoll[edi->sav.c_ind[i]], xs[edi->sav.anrs_loc[i]]);

    bFirst = FALSE;

    GMX_MPE_LOG(ev_edsam_finish);
}
