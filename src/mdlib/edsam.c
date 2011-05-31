/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
#include "nrjac.h"
#include "mtop_util.h"
#include "edsam.h"
#include "mpelogging.h"
#include "gmxfio.h"
#include "groupcoord.h"


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
    int    neig;     /* nr of eigenvectors             */
    int   *ieig;     /* index nrs of eigenvectors      */
    real  *stpsz;    /* stepsizes (per eigenvector)    */
    rvec  **vec;     /* eigenvector components         */
    real  *xproj;    /* instantaneous x projections    */
    real  *fproj;    /* instantaneous f projections    */
    real  radius;    /* instantaneous radius           */
    real  *refproj;  /* starting or target projecions  */
    /* When using flooding as harmonic restraint: The current reference projection
     * is at each step calculated from the initial refproj0 and the slope. */
    real  *refproj0,*refprojslope;
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
    gmx_bool bHarmonic;           /* Use flooding for harmonic restraint on
                                     the eigenvector                          */
    gmx_bool bConstForce;         /* Do not calculate a flooding potential,
                                     instead flood with a constant force      */
    real tau;
    real deltaF;
    real Efl;
    real kT;
    real Vfl;
    real dt;
    real constEfl;
    real alpha2;
    int flood_id;
    rvec *forces_cartesian;
    t_eigvec vecs;         /* use flooding for these */
} t_edflood;


/* This type is for the average, reference, target, and origin structure    */
typedef struct gmx_edx
{
    int           nr;             /* number of atoms this structure contains  */
    int           nr_loc;         /* number of atoms on local node            */
    int           *anrs;          /* atom index numbers                       */
    int           *anrs_loc;      /* local atom index numbers                 */
    int           nalloc_loc;     /* allocation size of anrs_loc              */
    int           *c_ind;         /* at which position of the whole anrs
                                   * array is a local atom?, i.e.
                                   * c_ind[0...nr_loc-1] gives the atom index
                                   * with respect to the collective
                                   * anrs[0...nr-1] array                     */
    rvec          *x;             /* positions for this structure             */
    rvec          *x_old;         /* used to keep track of the shift vectors
                                     such that the ED molecule can always be
                                     made whole in the parallel case          */
    real          *m;             /* masses                                   */
    real          mtot;           /* total mass (only used in sref)           */
    real          *sqrtm;         /* sqrt of the masses used for mass-
                                   * weighting of analysis (only used in sav) */
} t_gmx_edx;


typedef struct edpar
{
    int            nini;           /* total Nr of atoms                    */
    gmx_bool       fitmas;         /* true if trans fit with cm            */
    gmx_bool       pcamas;         /* true if mass-weighted PCA            */
    int            presteps;       /* number of steps to run without any
                                    *    perturbations ... just monitoring */
    int            outfrq;         /* freq (in steps) of writing to edo    */
    int            maxedsteps;     /* max nr of steps per cycle            */

    /* all gmx_edx datasets are copied to all nodes in the parallel case   */
    struct gmx_edx sref;           /* reference positions, to these fitting
                                    * will be done                         */
    gmx_bool       bRefEqAv;       /* If true, reference & average indices
                                    * are the same. Used for optimization  */
    struct gmx_edx sav;            /* average positions                    */
    struct gmx_edx star;           /* target positions                     */
    struct gmx_edx sori;           /* origin positions                     */

    t_edvecs       vecs;           /* eigenvectors                         */
    real           slope;          /* minimal slope in acceptance radexp   */

    gmx_bool       bNeedDoEdsam;   /* if any of the options mon, linfix, ...
                                    * is used (i.e. apart from flooding)   */
    t_edflood      flood;          /* parameters especially for flooding   */
    struct t_ed_buffer *buf;       /* handle to local buffers              */
    struct edpar   *next_edi;      /* Pointer to another ed dataset        */
} t_edpar;


typedef struct gmx_edsam
{
    int           eEDtype;        /* Type of ED: see enums above          */
    const char    *edinam;        /* name of ED sampling input file       */
    const char    *edonam;        /*                     output           */
    FILE          *edo;           /* output file pointer                  */
    t_edpar       *edpar;
    gmx_bool      bFirst;
    gmx_bool      bStartFromCpt;
} t_gmx_edsam;


struct t_do_edsam
{
    matrix old_rotmat;
    real oldrad;
    rvec old_transvec,older_transvec,transvec_compact;
    rvec *xcoll;         /* Positions from all nodes, this is the
                            collective set we work on.
                            These are the positions of atoms with
                            average structure indices */
    rvec *xc_ref;        /* same but with reference structure indices */
    ivec *shifts_xcoll;        /* Shifts for xcoll  */
    ivec *extra_shifts_xcoll;  /* xcoll shift changes since last NS step */
    ivec *shifts_xc_ref;       /* Shifts for xc_ref */
    ivec *extra_shifts_xc_ref; /* xc_ref shift changes since last NS step */
    gmx_bool bUpdateShifts;    /* TRUE in NS steps to indicate that the
                                  ED shifts for this ED dataset need to
                                  be updated */
};


/* definition of ED buffer structure */
struct t_ed_buffer
{
    struct t_fit_to_ref *           fit_to_ref;
    struct t_do_edfit *             do_edfit;
    struct t_do_edsam *             do_edsam;
    struct t_do_radcon *            do_radcon;
};


/* Function declarations */
static void fit_to_reference(rvec *xcoll,rvec transvec,matrix rotmat,t_edpar *edi);

static void translate_and_rotate(rvec *x,int nat,rvec transvec,matrix rotmat);
/* End function declarations */


/* Does not subtract average positions, projection on single eigenvector is returned
 * used by: do_linfix, do_linacc, do_radfix, do_radacc, do_radcon
 * Average position is subtracted in ed_apply_constraints prior to calling projectx
 */
static real projectx(t_edpar *edi, rvec *xcoll, rvec *vec)
{
    int  i;
    real proj=0.0;


    for (i=0; i<edi->sav.nr; i++)
        proj += edi->sav.sqrtm[i]*iprod(vec[i], xcoll[i]);

    return proj;
}


/* Specialized: projection is stored in vec->refproj
 * -> used for radacc, radfix, radcon  and center of flooding potential
 * subtracts average positions, projects vector x */
static void rad_project(t_edpar *edi, rvec *x, t_eigvec *vec, t_commrec *cr)
{
    int i;
    real rad=0.0;

    /* Subtract average positions */
    for (i = 0; i < edi->sav.nr; i++)
        rvec_dec(x[i], edi->sav.x[i]);

    for (i = 0; i < vec->neig; i++)
    {
        vec->refproj[i] = projectx(edi,x,vec->vec[i]);
        rad += pow((vec->refproj[i]-vec->xproj[i]),2);
    }
    vec->radius=sqrt(rad);

    /* Add average positions */
    for (i = 0; i < edi->sav.nr; i++)
        rvec_inc(x[i], edi->sav.x[i]);
}


/* Project vector x, subtract average positions prior to projection and add
 * them afterwards to retain the unchanged vector. Store in xproj. Mass-weighting
 * is applied. */
static void project_to_eigvectors(rvec       *x,    /* The positions to project to an eigenvector */
                                  t_eigvec   *vec,  /* The eigenvectors */
                                  t_edpar    *edi)
{
    int  i;


    if (!vec->neig) return;

    /* Subtract average positions */
    for (i=0; i<edi->sav.nr; i++)
        rvec_dec(x[i], edi->sav.x[i]);

    for (i=0; i<vec->neig; i++)
        vec->xproj[i] = projectx(edi, x, vec->vec[i]);

    /* Add average positions */
    for (i=0; i<edi->sav.nr; i++)
        rvec_inc(x[i], edi->sav.x[i]);
}


/* Project vector x onto all edi->vecs (mon, linfix,...) */
static void project(rvec      *x,     /* positions to project */
                    t_edpar   *edi)   /* edi data set */
{
    /* It is not more work to subtract the average position in every
     * subroutine again, because these routines are rarely used simultanely */
    project_to_eigvectors(x, &edi->vecs.mon   , edi);
    project_to_eigvectors(x, &edi->vecs.linfix, edi);
    project_to_eigvectors(x, &edi->vecs.linacc, edi);
    project_to_eigvectors(x, &edi->vecs.radfix, edi);
    project_to_eigvectors(x, &edi->vecs.radacc, edi);
    project_to_eigvectors(x, &edi->vecs.radcon, edi);
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
#ifdef DEBUGHELPERS
static void dump_xcoll(t_edpar *edi, struct t_do_edsam *buf, t_commrec *cr,
                       int step)
{
    int i;
    FILE *fp;
    char fn[STRLEN];
    rvec *xcoll;
    ivec *shifts, *eshifts;


    if (!MASTER(cr))
        return;

    xcoll   = buf->xcoll;
    shifts  = buf->shifts_xcoll;
    eshifts = buf->extra_shifts_xcoll;

    sprintf(fn, "xcolldump_step%d.txt", step);
    fp = fopen(fn, "w");

    for (i=0; i<edi->sav.nr; i++)
        fprintf(fp, "%d %9.5f %9.5f %9.5f   %d %d %d   %d %d %d\n",
                edi->sav.anrs[i]+1,
                xcoll[i][XX]  , xcoll[i][YY]  , xcoll[i][ZZ],
                shifts[i][XX] , shifts[i][YY] , shifts[i][ZZ],
                eshifts[i][XX], eshifts[i][YY], eshifts[i][ZZ]);

    fclose(fp);
}


/* Debug helper */
static void dump_edi_positions(FILE *out, struct gmx_edx *s, const char name[])
{
    int i;


    fprintf(out, "#%s positions:\n%d\n", name, s->nr);
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
static void dump_edi_eigenvecs(FILE *out, t_eigvec *ev,
                               const char name[], int length)
{
    int i,j;


    fprintf(out, "#%s eigenvectors:\n%d\n", name, ev->neig);
    /* Dump the data for every eigenvector: */
    for (i=0; i<ev->neig; i++)
    {
        fprintf(out, "EV %4d\ncomponents %d\nstepsize %f\nxproj %f\nfproj %f\nrefproj %f\nradius %f\nComponents:\n",
                ev->ieig[i], length, ev->stpsz[i], ev->xproj[i], ev->fproj[i], ev->refproj[i], ev->radius);
        for (j=0; j<length; j++)
            fprintf(out, "%11.6f %11.6f %11.6f\n", ev->vec[i][j][XX], ev->vec[i][j][YY], ev->vec[i][j][ZZ]);
    }
}


/* Debug helper */
static void dump_edi(t_edpar *edpars, t_commrec *cr, int nr_edi)
{
    FILE  *out;
    char  fn[STRLEN];


    sprintf(fn, "EDdump_node%d_edi%d", cr->nodeid, nr_edi);
    out = ffopen(fn, "w");

    fprintf(out,"#NINI\n %d\n#FITMAS\n %d\n#ANALYSIS_MAS\n %d\n",
            edpars->nini,edpars->fitmas,edpars->pcamas);
    fprintf(out,"#OUTFRQ\n %d\n#MAXLEN\n %d\n#SLOPECRIT\n %f\n",
            edpars->outfrq,edpars->maxedsteps,edpars->slope);
    fprintf(out,"#PRESTEPS\n %d\n#DELTA_F0\n %f\n#TAU\n %f\n#EFL_NULL\n %f\n#ALPHA2\n %f\n",
            edpars->presteps,edpars->flood.deltaF0,edpars->flood.tau,
            edpars->flood.constEfl,edpars->flood.alpha2);

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
    dump_edi_eigenvecs(out, &edpars->flood.vecs, "FLOODING"  , edpars->sav.nr);

    /* Dump ed local buffer */
    fprintf(out, "buf->do_edfit         =%p\n", (void*)edpars->buf->do_edfit  );
    fprintf(out, "buf->do_edsam         =%p\n", (void*)edpars->buf->do_edsam  );
    fprintf(out, "buf->do_radcon        =%p\n", (void*)edpars->buf->do_radcon );

    ffclose(out);
}


/* Debug helper */
static void dump_rotmat(FILE* out,matrix rotmat)
{
    fprintf(out,"ROTMAT: %12.8f %12.8f %12.8f\n",rotmat[XX][XX],rotmat[XX][YY],rotmat[XX][ZZ]);
    fprintf(out,"ROTMAT: %12.8f %12.8f %12.8f\n",rotmat[YY][XX],rotmat[YY][YY],rotmat[YY][ZZ]);
    fprintf(out,"ROTMAT: %12.8f %12.8f %12.8f\n",rotmat[ZZ][XX],rotmat[ZZ][YY],rotmat[ZZ][ZZ]);
}


/* Debug helper */
static void dump_rvec(FILE *out, int dim, rvec *x)
{
    int i;


    for (i=0; i<dim; i++)
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
#endif


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
    gmx_bool bFirst;

    if(edi->buf->do_edfit != NULL)
        bFirst = FALSE;
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


static void rmfit(int nat, rvec *xcoll, rvec transvec, matrix rotmat)
{
    rvec vec;
    matrix tmat;


    /* Remove rotation.
     * The inverse rotation is described by the transposed rotation matrix */
    transpose(rotmat,tmat);
    rotate_x(xcoll, nat, tmat);

    /* Remove translation */
    vec[XX]=-transvec[XX];
    vec[YY]=-transvec[YY];
    vec[ZZ]=-transvec[ZZ];
    translate_x(xcoll, nat, vec);
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

To use restraints with harmonic potentials switch -restrain and -harmonic. Then the eigenvalues are
used as spring constants for the harmonic potential.
Note that eq3 in the flooding paper (J. Comp. Chem. 2006, 27, 1693-1702) defines the parameter lambda \
as the inverse of the spring constant, whereas the implementation uses lambda as the spring constant.

To use more than one flooding matrix just concatenate several .edi files (cat flood1.edi flood2.edi > flood_all.edi)
the routine read_edi_file reads all of theses flooding files.
The structure t_edi is now organized as a list of t_edis and the function do_flood cycles through the list
calling the do_single_flood() routine for every single entry. Since every state variables have been kept in one
edi there is no interdependence whatsoever. The forces are added together.

  To write energies into the .edr file, call the function
        get_flood_enx_names(char**, int *nnames) to get the Header (Vfl1 Vfl2... Vfln)
and call
        get_flood_energies(real Vfl[],int nnames);

  TODO:
- one could program the whole thing such that Efl, Vfl and deltaF is written to the .edr file. -- i dont know how to do that, yet.

  Maybe one should give a range of atoms for which to remove motion, so that motion is removed with
  two edsam files from two peptide chains
*/

static void write_edo_flood(t_edpar *edi, FILE *fp, gmx_large_int_t step)
{
    int i;
    char buf[22];
    gmx_bool bOutputRef=FALSE;


    fprintf(fp,"%d.th FL: %s %12.5e %12.5e %12.5e\n",
            edi->flood.flood_id, gmx_step_str(step,buf),
            edi->flood.Efl, edi->flood.Vfl, edi->flood.deltaF);


    /* Check whether any of the references changes with time (this can happen
     * in case flooding is used as harmonic restraint). If so, output all the
     * current reference projections. */
    if (edi->flood.bHarmonic)
    {
        for (i = 0; i < edi->flood.vecs.neig; i++)
        {
            if (edi->flood.vecs.refprojslope[i] != 0.0)
                bOutputRef=TRUE;
        }
        if (bOutputRef)
        {
            fprintf(fp, "Ref. projs.: ");
            for (i = 0; i < edi->flood.vecs.neig; i++)
            {
                fprintf(fp, "%12.5e ", edi->flood.vecs.refproj[i]);
            }
            fprintf(fp, "\n");
        }
    }
    fprintf(fp,"FL_FORCES: ");

    for (i=0; i<edi->flood.vecs.neig; i++)
        fprintf(fp," %12.5e",edi->flood.vecs.fproj[i]);

    fprintf(fp,"\n");
}


/* From flood.xproj compute the Vfl(x) at this point */
static real flood_energy(t_edpar *edi, gmx_large_int_t step)
{
    /* compute flooding energy Vfl
     Vfl = Efl * exp( - \frac {kT} {2Efl alpha^2} * sum_i { \lambda_i c_i^2 } )
     \lambda_i is the reciprocal eigenvalue 1/\sigma_i
         it is already computed by make_edi and stored in stpsz[i]
     bHarmonic:
       Vfl = - Efl * 1/2(sum _i {\frac 1{\lambda_i} c_i^2})
     */
    real sum;
    real Vfl;
    int i;


    /* Each time this routine is called (i.e. each time step), we add a small
     * value to the reference projection. This way a harmonic restraint towards
     * a moving reference is realized. If no value for the additive constant
     * is provided in the edi file, the reference will not change. */
    if (edi->flood.bHarmonic)
    {
        for (i=0; i<edi->flood.vecs.neig; i++)
        {
            edi->flood.vecs.refproj[i] = edi->flood.vecs.refproj0[i] + step * edi->flood.vecs.refprojslope[i];
        }
    }

    sum=0.0;
    /* Compute sum which will be the exponent of the exponential */
    for (i=0; i<edi->flood.vecs.neig; i++)
    {
        /* stpsz stores the reciprocal eigenvalue 1/sigma_i */
        sum += edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i])*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]);
    }

    /* Compute the Gauss function*/
    if (edi->flood.bHarmonic)
    {
        Vfl = -0.5*edi->flood.Efl*sum;  /* minus sign because Efl is negative, if restrain is on. */
    }
    else
    {
        Vfl = edi->flood.Efl!=0 ? edi->flood.Efl*exp(-edi->flood.kT/2/edi->flood.Efl/edi->flood.alpha2*sum) :0;
    }

    return Vfl;
}


/* From the position and from Vfl compute forces in subspace -> store in edi->vec.flood.fproj */
static void flood_forces(t_edpar *edi)
{
    /* compute the forces in the subspace of the flooding eigenvectors
     * by the formula F_i= V_{fl}(c) * ( \frac {kT} {E_{fl}} \lambda_i c_i */

    int i;
    real energy=edi->flood.Vfl;


    if (edi->flood.bHarmonic)
        for (i=0; i<edi->flood.vecs.neig; i++)
        {
            edi->flood.vecs.fproj[i] = edi->flood.Efl* edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]);
        }
    else
        for (i=0; i<edi->flood.vecs.neig; i++)
        {
            /* if Efl is zero the forces are zero if not use the formula */
            edi->flood.vecs.fproj[i] = edi->flood.Efl!=0 ? edi->flood.kT/edi->flood.Efl/edi->flood.alpha2*energy*edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]) : 0;
        }
}


/* Raise forces from subspace into cartesian space */
static void flood_blowup(t_edpar *edi, rvec *forces_cart)
{
    /* this function lifts the forces from the subspace to the cartesian space
     all the values not contained in the subspace are assumed to be zero and then
     a coordinate transformation from eigenvector to cartesian vectors is performed
     The nonexistent values don't have to be set to zero explicitly, they would occur
     as zero valued summands, hence we just stop to compute this part of the sum.

     for every atom we add all the contributions to this atom from all the different eigenvectors.

     NOTE: one could add directly to the forcefield forces, would mean we wouldn't have to clear the
     field forces_cart prior the computation, but we compute the forces separately
     to have them accessible for diagnostics
     */
    int  j,eig;
    rvec dum;
    real *forces_sub;


    forces_sub = edi->flood.vecs.fproj;


    /* Calculate the cartesian forces for the local atoms */

    /* Clear forces first */
    for (j=0; j<edi->sav.nr_loc; j++)
        clear_rvec(forces_cart[j]);

    /* Now compute atomwise */
    for (j=0; j<edi->sav.nr_loc; j++)
    {
        /* Compute forces_cart[edi->sav.anrs[j]] */
        for (eig=0; eig<edi->flood.vecs.neig; eig++)
        {
            /* Force vector is force * eigenvector (compute only atom j) */
            svmul(forces_sub[eig],edi->flood.vecs.vec[eig][edi->sav.c_ind[j]],dum);
            /* Add this vector to the cartesian forces */
            rvec_inc(forces_cart[j],dum);
        }
    }
}


/* Update the values of Efl, deltaF depending on tau and Vfl */
static void update_adaption(t_edpar *edi)
{
    /* this function updates the parameter Efl and deltaF according to the rules given in
     * 'predicting unimolecular chemical reactions: chemical flooding' M Mueller et al,
     * J. chem Phys. */

    if ((edi->flood.tau < 0 ? -edi->flood.tau : edi->flood.tau ) > 0.00000001)
    {
        edi->flood.Efl = edi->flood.Efl+edi->flood.dt/edi->flood.tau*(edi->flood.deltaF0-edi->flood.deltaF);
        /* check if restrain (inverted flooding) -> don't let EFL become positive */
        if (edi->flood.alpha2<0 && edi->flood.Efl>-0.00000001)
            edi->flood.Efl = 0;

        edi->flood.deltaF = (1-edi->flood.dt/edi->flood.tau)*edi->flood.deltaF+edi->flood.dt/edi->flood.tau*edi->flood.Vfl;
    }
}


static void do_single_flood(
        FILE *edo,
        rvec x[],
        rvec force[],
        t_edpar *edi,
        gmx_large_int_t step,
        matrix box,
        t_commrec *cr)
{
    int i;
    matrix  rotmat;         /* rotation matrix */
    matrix  tmat;           /* inverse rotation */
    rvec    transvec;       /* translation vector */
    struct t_do_edsam *buf;


    buf=edi->buf->do_edsam;

    /* Broadcast the positions of the AVERAGE structure such that they are known on
     * every processor. Each node contributes its local positions x and stores them in
     * the collective ED array buf->xcoll */
    communicate_group_positions(cr, buf->xcoll, buf->shifts_xcoll, buf->extra_shifts_xcoll, buf->bUpdateShifts, x,
                    edi->sav.nr, edi->sav.nr_loc, edi->sav.anrs_loc, edi->sav.c_ind, edi->sav.x_old, box);

    /* Only assembly REFERENCE positions if their indices differ from the average ones */
    if (!edi->bRefEqAv)
        communicate_group_positions(cr, buf->xc_ref, buf->shifts_xc_ref, buf->extra_shifts_xc_ref, buf->bUpdateShifts, x,
                edi->sref.nr, edi->sref.nr_loc, edi->sref.anrs_loc, edi->sref.c_ind, edi->sref.x_old, box);

    /* If bUpdateShifts was TRUE, the shifts have just been updated in get_positions.
     * We do not need to update the shifts until the next NS step */
    buf->bUpdateShifts = FALSE;

    /* Now all nodes have all of the ED/flooding positions in edi->sav->xcoll,
     * as well as the indices in edi->sav.anrs */

    /* Fit the reference indices to the reference structure */
    if (edi->bRefEqAv)
        fit_to_reference(buf->xcoll , transvec, rotmat, edi);
    else
        fit_to_reference(buf->xc_ref, transvec, rotmat, edi);

    /* Now apply the translation and rotation to the ED structure */
    translate_and_rotate(buf->xcoll, edi->sav.nr, transvec, rotmat);

    /* Project fitted structure onto supbspace -> store in edi->flood.vecs.xproj */
    project_to_eigvectors(buf->xcoll,&edi->flood.vecs,edi);

    if (FALSE == edi->flood.bConstForce)
    {
        /* Compute Vfl(x) from flood.xproj */
        edi->flood.Vfl = flood_energy(edi, step);

        update_adaption(edi);

        /* Compute the flooding forces */
        flood_forces(edi);
    }

    /* Translate them into cartesian positions */
    flood_blowup(edi, edi->flood.forces_cartesian);

    /* Rotate forces back so that they correspond to the given structure and not to the fitted one */
    /* Each node rotates back its local forces */
    transpose(rotmat,tmat);
    rotate_x(edi->flood.forces_cartesian, edi->sav.nr_loc, tmat);

    /* Finally add forces to the main force variable */
    for (i=0; i<edi->sav.nr_loc; i++)
        rvec_inc(force[edi->sav.anrs_loc[i]],edi->flood.forces_cartesian[i]);

    /* Output is written by the master process */
    if (do_per_step(step,edi->outfrq) && MASTER(cr))
        write_edo_flood(edi,edo,step);
}


/* Main flooding routine, called from do_force */
extern void do_flood(
        FILE            *log,    /* md.log file */
        t_commrec       *cr,     /* Communication record */
        rvec            x[],     /* Positions on the local processor */
        rvec            force[], /* forcefield forces, to these the flooding forces are added */
        gmx_edsam_t     ed,      /* ed data structure contains all ED and flooding datasets */
        matrix          box,     /* the box */
        gmx_large_int_t step)    /* The relative time step since ir->init_step is already subtracted */
{
    t_edpar *edi;


    if (ed->eEDtype != eEDflood)
        return;

    edi = ed->edpar;
    while (edi)
    {
        /* Call flooding for one matrix */
        if (edi->flood.vecs.neig)
            do_single_flood(ed->edo,x,force,edi,step,box,cr);
        edi = edi->next_edi;
    }
}


/* Called by init_edi, configure some flooding related variables and structures,
 * print headers to output files */
static void init_flood(t_edpar *edi, gmx_edsam_t ed, real dt, t_commrec *cr)
{
    int i;


    edi->flood.Efl = edi->flood.constEfl;
    edi->flood.Vfl = 0;
    edi->flood.dt  = dt;

    if (edi->flood.vecs.neig)
    {
        /* If in any of the datasets we find a flooding vector, flooding is turned on */
        ed->eEDtype = eEDflood;

        fprintf(stderr,"ED: Flooding of matrix %d is switched on.\n", edi->flood.flood_id);

        if (edi->flood.bConstForce)
        {
            /* We have used stpsz as a vehicle to carry the fproj values for constant
             * force flooding. Now we copy that to flood.vecs.fproj. Note that
             * in const force flooding, fproj is never changed. */
            for (i=0; i<edi->flood.vecs.neig; i++)
            {
                edi->flood.vecs.fproj[i] = edi->flood.vecs.stpsz[i];

                fprintf(stderr, "ED: applying on eigenvector %d a constant force of %g\n",
                        edi->flood.vecs.ieig[i], edi->flood.vecs.fproj[i]);
            }
        }
        fprintf(ed->edo,"FL_HEADER: Flooding of matrix %d is switched on! The flooding output will have the following format:\n",
                edi->flood.flood_id);
        fprintf(ed->edo,"FL_HEADER: Step     Efl          Vfl       deltaF\n");
    }
}


#ifdef DEBUGHELPERS
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
#endif


gmx_edsam_t ed_open(int nfile,const t_filenm fnm[],unsigned long Flags,t_commrec *cr)
{
    gmx_edsam_t ed;


    /* Allocate space for the ED data structure */
    snew(ed, 1);

    /* We want to perform ED (this switch might later be upgraded to eEDflood) */
    ed->eEDtype = eEDedsam;

    if (MASTER(cr))
    {
        /* Open .edi input file: */
        ed->edinam=ftp2fn(efEDI,nfile,fnm);
        /* The master opens the .edo output file */
        fprintf(stderr,"ED sampling will be performed!\n");
        ed->edonam = ftp2fn(efEDO,nfile,fnm);
        ed->edo    = gmx_fio_fopen(ed->edonam,(Flags & MD_APPENDFILES)? "a+" : "w+");
        ed->bStartFromCpt = Flags & MD_STARTFROMCPT;
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
        /* Local atom indices get assigned in dd_make_local_group_indices.
         * There, also memory is allocated */
        s->nalloc_loc = 0;              /* allocation size of s->anrs_loc */
        snew_bc(cr, s->x_old, s->nr);   /* To be able to always make the ED molecule whole, ...        */
        nblock_bc(cr, s->nr, s->x_old); /* ... keep track of shift changes with the help of old coords */
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
static void bc_ed_vecs(t_commrec *cr, t_eigvec *ev, int length, gmx_bool bHarmonic)
{
    int i;

    snew_bc(cr, ev->ieig   , ev->neig);  /* index numbers of eigenvector  */
    snew_bc(cr, ev->stpsz  , ev->neig);  /* stepsizes per eigenvector     */
    snew_bc(cr, ev->xproj  , ev->neig);  /* instantaneous x projection    */
    snew_bc(cr, ev->fproj  , ev->neig);  /* instantaneous f projection    */
    snew_bc(cr, ev->refproj, ev->neig);  /* starting or target projection */

    nblock_bc(cr, ev->neig, ev->ieig   );
    nblock_bc(cr, ev->neig, ev->stpsz  );
    nblock_bc(cr, ev->neig, ev->xproj  );
    nblock_bc(cr, ev->neig, ev->fproj  );
    nblock_bc(cr, ev->neig, ev->refproj);

    snew_bc(cr, ev->vec, ev->neig);      /* Eigenvector components        */
    for (i=0; i<ev->neig; i++)
    {
        snew_bc(cr, ev->vec[i], length);
        nblock_bc(cr, length, ev->vec[i]);
    }

    /* For harmonic restraints the reference projections can change with time */
    if (bHarmonic)
    {
        snew_bc(cr, ev->refproj0    , ev->neig);
        snew_bc(cr, ev->refprojslope, ev->neig);
        nblock_bc(cr, ev->neig, ev->refproj0    );
        nblock_bc(cr, ev->neig, ev->refprojslope);
    }
}


/* Broadcasts the ED / flooding data to other nodes
 * and allocates memory where needed */
static void broadcast_ed_data(t_commrec *cr, gmx_edsam_t ed, int numedis)
{
    int     nr;
    t_edpar *edi;


    /* Master lets the other nodes know if its ED only or also flooding */
    gmx_bcast(sizeof(ed->eEDtype), &(ed->eEDtype), cr);

    snew_bc(cr, ed->edpar,1);
    /* Now transfer the ED data set(s) */
    edi = ed->edpar;
    for (nr=0; nr<numedis; nr++)
    {
        /* Broadcast a single ED data set */
        block_bc(cr, *edi);

        /* Broadcast positions */
        bc_ed_positions(cr, &(edi->sref), eedREF); /* reference positions (don't broadcast masses)    */
        bc_ed_positions(cr, &(edi->sav ), eedAV ); /* average positions (do broadcast masses as well) */
        bc_ed_positions(cr, &(edi->star), eedTAR); /* target positions                                */
        bc_ed_positions(cr, &(edi->sori), eedORI); /* origin positions                                */

        /* Broadcast eigenvectors */
        bc_ed_vecs(cr, &edi->vecs.mon   , edi->sav.nr, FALSE);
        bc_ed_vecs(cr, &edi->vecs.linfix, edi->sav.nr, FALSE);
        bc_ed_vecs(cr, &edi->vecs.linacc, edi->sav.nr, FALSE);
        bc_ed_vecs(cr, &edi->vecs.radfix, edi->sav.nr, FALSE);
        bc_ed_vecs(cr, &edi->vecs.radacc, edi->sav.nr, FALSE);
        bc_ed_vecs(cr, &edi->vecs.radcon, edi->sav.nr, FALSE);
        /* Broadcast flooding eigenvectors and, if needed, values for the moving reference */
        bc_ed_vecs(cr, &edi->flood.vecs,  edi->sav.nr, edi->flood.bHarmonic);

        /* Set the pointer to the next ED dataset */
        if (edi->next_edi)
        {
          snew_bc(cr, edi->next_edi, 1);
          edi = edi->next_edi;
        }
    }
}


/* init-routine called for every *.edi-cycle, initialises t_edpar structure */
static void init_edi(gmx_mtop_t *mtop,t_inputrec *ir,
                     t_commrec *cr,gmx_edsam_t ed,t_edpar *edi)
{
    int  i;
    real totalmass = 0.0;
    rvec com;
    t_atom *atom;

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
        {
            gmx_mtop_atomnr_to_atom(mtop,edi->sref.anrs[i],&atom);
            edi->sref.m[i] = atom->m;
        }
        else
        {
            edi->sref.m[i] = 1.0;
        }

        /* Check that every m > 0. Bad things will happen otherwise. */
        if (edi->sref.m[i] <= 0.0)
        {
            gmx_fatal(FARGS, "Reference structure atom %d (sam.edi index %d) has a mass of %g.\n"
                             "For a mass-weighted fit, all reference structure atoms need to have a mass >0.\n"
                             "Either make the covariance analysis non-mass-weighted, or exclude massless\n"
                             "atoms from the reference structure by creating a proper index group.\n",
                      i, edi->sref.anrs[i]+1, edi->sref.m[i]);
        }

        totalmass += edi->sref.m[i];
    }
    edi->sref.mtot = totalmass;

    /* Masses m and sqrt(m) for the average structure. Note that m
     * is needed if forces have to be evaluated in do_edsam */
    snew(edi->sav.sqrtm, edi->sav.nr );
    snew(edi->sav.m    , edi->sav.nr );
    for (i = 0; i < edi->sav.nr; i++)
    {
        gmx_mtop_atomnr_to_atom(mtop,edi->sav.anrs[i],&atom);
        edi->sav.m[i] = atom->m;
        if (edi->pcamas)
        {
            edi->sav.sqrtm[i] = sqrt(atom->m);
        }
        else
        {
            edi->sav.sqrtm[i] = 1.0;
        }

        /* Check that every m > 0. Bad things will happen otherwise. */
        if (edi->sav.sqrtm[i] <= 0.0)
        {
            gmx_fatal(FARGS, "Average structure atom %d (sam.edi index %d) has a mass of %g.\n"
                             "For ED with mass-weighting, all average structure atoms need to have a mass >0.\n"
                             "Either make the covariance analysis non-mass-weighted, or exclude massless\n"
                             "atoms from the average structure by creating a proper index group.\n",
                      i, edi->sav.anrs[i]+1, atom->m);
        }
    }

    /* put reference structure in origin */
    get_center(edi->sref.x, edi->sref.m, edi->sref.nr, com);
    com[XX] = -com[XX];
    com[YY] = -com[YY];
    com[ZZ] = -com[ZZ];
    translate_x(edi->sref.x, edi->sref.nr, com);

    /* Init ED buffer */
    snew(edi->buf, 1);
}


static void check(const char *line, const char *label)
{
    if (!strstr(line,label))
        gmx_fatal(FARGS,"Could not find input parameter %s at expected position in edsam input-file (.edi)\nline read instead is %s",label,line);
}


static int read_checked_edint(FILE *file,const char *label)
{
    char line[STRLEN+1];
    int idum;


    fgets2 (line,STRLEN,file);
    check(line,label);
    fgets2 (line,STRLEN,file);
    sscanf (line,"%d",&idum);
    return idum;
}


static int read_edint(FILE *file,gmx_bool *bEOF)
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


static real read_checked_edreal(FILE *file,const char *label)
{
    char line[STRLEN+1];
    double rdum;


    fgets2 (line,STRLEN,file);
    check(line,label);
    fgets2 (line,STRLEN,file);
    sscanf (line,"%lf",&rdum);
    return (real) rdum; /* always read as double and convert to single */
}


static void read_edx(FILE *file,int number,int *anrs,rvec *x)
{
    int i,j;
    char line[STRLEN+1];
    double d[3];


    for(i=0; i<number; i++)
    {
        fgets2 (line,STRLEN,file);
        sscanf (line,"%d%lf%lf%lf",&anrs[i],&d[0],&d[1],&d[2]);
        anrs[i]--; /* we are reading FORTRAN indices */
        for(j=0; j<3; j++)
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


static void read_edvec(FILE *in,int nr,t_eigvec *tvec,gmx_bool bReadRefproj, gmx_bool *bHaveReference)
{
    int i,idum,nscan;
    double rdum,refproj_dum=0.0,refprojslope_dum=0.0;
    char line[STRLEN+1];


    tvec->neig=read_checked_edint(in,"NUMBER OF EIGENVECTORS");
    if (tvec->neig >0)
    {
        snew(tvec->ieig   ,tvec->neig);
        snew(tvec->stpsz  ,tvec->neig);
        snew(tvec->vec    ,tvec->neig);
        snew(tvec->xproj  ,tvec->neig);
        snew(tvec->fproj  ,tvec->neig);
        snew(tvec->refproj,tvec->neig);
        if (bReadRefproj)
        {
            snew(tvec->refproj0    ,tvec->neig);
            snew(tvec->refprojslope,tvec->neig);
        }

        for(i=0; (i < tvec->neig); i++)
        {
            fgets2 (line,STRLEN,in);
            if (bReadRefproj) /* ONLY when using flooding as harmonic restraint */
            {
                nscan = sscanf(line,"%d%lf%lf%lf",&idum,&rdum,&refproj_dum,&refprojslope_dum);
                /* Zero out values which were not scanned */
                switch(nscan)
                {
                    case 4:
                        /* Every 4 values read, including reference position */
                        *bHaveReference = TRUE;
                        break;
                    case 3:
                        /* A reference position is provided */
                        *bHaveReference = TRUE;
                        /* No value for slope, set to 0 */
                        refprojslope_dum = 0.0;
                        break;
                    case 2:
                        /* No values for reference projection and slope, set to 0 */
                        refproj_dum      = 0.0;
                        refprojslope_dum = 0.0;
                        break;
                    default:
                        gmx_fatal(FARGS,"Expected 2 - 4 (not %d) values for flooding vec: <nr> <spring const> <refproj> <refproj-slope>\n", nscan);
                        break;
                }
                tvec->refproj[i]=refproj_dum;
                tvec->refproj0[i]=refproj_dum;
                tvec->refprojslope[i]=refprojslope_dum;
            }
            else /* Normal flooding */
            {
                nscan = sscanf(line,"%d%lf",&idum,&rdum);
                if (nscan != 2)
                    gmx_fatal(FARGS,"Expected 2 values for flooding vec: <nr> <stpsz>\n");
            }
            tvec->ieig[i]=idum;
            tvec->stpsz[i]=rdum;
        } /* end of loop over eigenvectors */

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
	gmx_bool bHaveReference = FALSE;


    read_edvec(in, nr, &vecs->mon   , FALSE, &bHaveReference);
    read_edvec(in, nr, &vecs->linfix, FALSE, &bHaveReference);
    read_edvec(in, nr, &vecs->linacc, FALSE, &bHaveReference);
    read_edvec(in, nr, &vecs->radfix, FALSE, &bHaveReference);
    read_edvec(in, nr, &vecs->radacc, FALSE, &bHaveReference);
    read_edvec(in, nr, &vecs->radcon, FALSE, &bHaveReference);
}


/* Check if the same atom indices are used for reference and average positions */
static gmx_bool check_if_same(struct gmx_edx sref, struct gmx_edx sav)
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


static int read_edi(FILE* in, gmx_edsam_t ed,t_edpar *edi,int nr_mdatoms, int edi_nr, t_commrec *cr)
{
    int readmagic;
    const int magic=670;
    gmx_bool bEOF;

    /* Was a specific reference point for the flooding/umbrella potential provided in the edi file? */
    gmx_bool bHaveReference = FALSE;


    /* the edi file is not free format, so expect problems if the input is corrupt. */

    /* check the magic number */
    readmagic=read_edint(in,&bEOF);
    /* Check whether we have reached the end of the input file */
    if (bEOF)
        return 0;

    if (readmagic != magic)
    {
        if (readmagic==666 || readmagic==667 || readmagic==668)
            gmx_fatal(FARGS,"Wrong magic number: Use newest version of make_edi to produce edi file");
        else if (readmagic == 669)
            ;
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
    if (readmagic > 669)
        edi->flood.bConstForce = read_checked_edint(in,"CONST_FORCE_FLOODING");
    else
        edi->flood.bConstForce = FALSE;
    edi->flood.flood_id  = edi_nr;
    edi->sref.nr         = read_checked_edint(in,"NREF");

    /* allocate space for reference positions and read them */
    snew(edi->sref.anrs,edi->sref.nr);
    snew(edi->sref.x   ,edi->sref.nr);
    if (PAR(cr))
        snew(edi->sref.x_old,edi->sref.nr);
    edi->sref.sqrtm    =NULL;
    read_edx(in,edi->sref.nr,edi->sref.anrs,edi->sref.x);

    /* average positions. they define which atoms will be used for ED sampling */
    edi->sav.nr=read_checked_edint(in,"NAV");
    snew(edi->sav.anrs,edi->sav.nr);
    snew(edi->sav.x   ,edi->sav.nr);
    if (PAR(cr))
        snew(edi->sav.x_old,edi->sav.nr);
    read_edx(in,edi->sav.nr,edi->sav.anrs,edi->sav.x);

    /* Check if the same atom indices are used for reference and average positions */
    edi->bRefEqAv = check_if_same(edi->sref, edi->sav);

    /* eigenvectors */
    read_edvecs(in,edi->sav.nr,&edi->vecs);
    read_edvec(in,edi->sav.nr,&edi->flood.vecs,edi->flood.bHarmonic, &bHaveReference);

    /* target positions */
    edi->star.nr=read_edint(in,&bEOF);
    if (edi->star.nr > 0)
    {
        snew(edi->star.anrs,edi->star.nr);
        snew(edi->star.x   ,edi->star.nr);
        edi->star.sqrtm    =NULL;
        read_edx(in,edi->star.nr,edi->star.anrs,edi->star.x);
    }

    /* positions defining origin of expansion circle */
    edi->sori.nr=read_edint(in,&bEOF);
    if (edi->sori.nr > 0)
    {
    	if (bHaveReference)
    	{
    		/* Both an -ori structure and a at least one manual reference point have been
    		 * specified. That's ambiguous and probably not intentional. */
    		gmx_fatal(FARGS, "ED: An origin structure has been provided and a at least one (moving) reference\n"
    		                 "    point was manually specified in the edi file. That is ambiguous. Aborting.\n");
    	}
        snew(edi->sori.anrs,edi->sori.nr);
        snew(edi->sori.x   ,edi->sori.nr);
        edi->sori.sqrtm    =NULL;
        read_edx(in,edi->sori.nr,edi->sori.anrs,edi->sori.x);
    }

    /* all done */
    return 1;
}



/* Read in the edi input file. Note that it may contain several ED data sets which were
 * achieved by concatenating multiple edi files. The standard case would be a single ED
 * data set, though. */
static void read_edi_file(gmx_edsam_t ed, t_edpar *edi, int nr_mdatoms, t_commrec *cr)
{
    FILE    *in;
    t_edpar *curr_edi,*last_edi;
    t_edpar *edi_read;
    int     edi_nr = 0;


    /* This routine is executed on the master only */

    /* Open the .edi parameter input file */
    in = gmx_fio_fopen(ed->edinam,"r");
    fprintf(stderr, "ED: Reading edi file %s\n", ed->edinam);

    /* Now read a sequence of ED input parameter sets from the edi file */
    curr_edi=edi;
    last_edi=edi;
    while( read_edi(in, ed, curr_edi, nr_mdatoms, edi_nr, cr) )
    {
        edi_nr++;
        /* Make shure that the number of atoms in each dataset is the same as in the tpr file */
        if (edi->nini != nr_mdatoms)
            gmx_fatal(FARGS,"edi file %s (dataset #%d) was made for %d atoms, but the simulation contains %d atoms.",
                    ed->edinam, edi_nr, edi->nini, nr_mdatoms);
        /* Since we arrived within this while loop we know that there is still another data set to be read in */
        /* We need to allocate space for the data: */
        snew(edi_read,1);
        /* Point the 'next_edi' entry to the next edi: */
        curr_edi->next_edi=edi_read;
        /* Keep the curr_edi pointer for the case that the next dataset is empty: */
        last_edi = curr_edi;
        /* Let's prepare to read in the next edi data set: */
        curr_edi = edi_read;
    }
    if (edi_nr == 0)
        gmx_fatal(FARGS, "No complete ED data set found in edi file %s.", ed->edinam);

    /* Terminate the edi dataset list with a NULL pointer: */
    last_edi->next_edi = NULL;

    fprintf(stderr, "ED: Found %d ED dataset%s.\n", edi_nr, edi_nr>1? "s" : "");

    /* Close the .edi file again */
    gmx_fio_fclose(in);
}


struct t_fit_to_ref {
    rvec *xcopy;       /* Working copy of the positions in fit_to_reference */
};

/* Fit the current positions to the reference positions
 * Do not actually do the fit, just return rotation and translation.
 * Note that the COM of the reference structure was already put into
 * the origin by init_edi. */
static void fit_to_reference(rvec      *xcoll,    /* The positions to be fitted */
                             rvec      transvec,  /* The translation vector */
                             matrix    rotmat,    /* The rotation matrix */
                             t_edpar   *edi)      /* Just needed for do_edfit */
{
    rvec com;          /* center of mass */
    int  i;
    struct t_fit_to_ref *loc;


    GMX_MPE_LOG(ev_fit_to_reference_start);

    /* Allocate memory the first time this routine is called for each edi dataset */
    if (NULL == edi->buf->fit_to_ref)
    {
        snew(edi->buf->fit_to_ref, 1);
        snew(edi->buf->fit_to_ref->xcopy, edi->sref.nr);
    }
    loc = edi->buf->fit_to_ref;

    /* We do not touch the original positions but work on a copy. */
    for (i=0; i<edi->sref.nr; i++)
        copy_rvec(xcoll[i], loc->xcopy[i]);

    /* Calculate the center of mass */
    get_center(loc->xcopy, edi->sref.m, edi->sref.nr, com);

    transvec[XX] = -com[XX];
    transvec[YY] = -com[YY];
    transvec[ZZ] = -com[ZZ];

    /* Subtract the center of mass from the copy */
    translate_x(loc->xcopy, edi->sref.nr, transvec);

    /* Determine the rotation matrix */
    do_edfit(edi->sref.nr, edi->sref.x, loc->xcopy, rotmat, edi);

    GMX_MPE_LOG(ev_fit_to_reference_finish);
}


static void translate_and_rotate(rvec *x,         /* The positions to be translated and rotated */
                                 int nat,         /* How many positions are there? */
                                 rvec transvec,   /* The translation vector */
                                 matrix rotmat)   /* The rotation matrix */
{
    /* Translation */
    translate_x(x, nat, transvec);

    /* Rotation */
    rotate_x(x, nat, rotmat);
}


/* Gets the rms deviation of the positions to the structure s */
/* fit_to_structure has to be called before calling this routine! */
static real rmsd_from_structure(rvec           *x,  /* The positions under consideration */
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


void dd_make_local_ed_indices(gmx_domdec_t *dd, struct gmx_edsam *ed)
{
    t_edpar *edi;


    if (ed->eEDtype != eEDnone)
    {
        /* Loop over ED datasets (usually there is just one dataset, though) */
        edi=ed->edpar;
        while (edi)
        {
            /* Local atoms of the reference structure (for fitting), need only be assembled
             * if their indices differ from the average ones */
            if (!edi->bRefEqAv)
                dd_make_local_group_indices(dd->ga2la, edi->sref.nr, edi->sref.anrs,
                        &edi->sref.nr_loc, &edi->sref.anrs_loc, &edi->sref.nalloc_loc, edi->sref.c_ind);

            /* Local atoms of the average structure (on these ED will be performed) */
            dd_make_local_group_indices(dd->ga2la, edi->sav.nr, edi->sav.anrs,
                    &edi->sav.nr_loc, &edi->sav.anrs_loc, &edi->sav.nalloc_loc, edi->sav.c_ind);

            /* Indicate that the ED shift vectors for this structure need to be updated
             * at the next call to communicate_group_positions, since obviously we are in a NS step */
            edi->buf->do_edsam->bUpdateShifts = TRUE;

            /* Set the pointer to the next ED dataset (if any) */
            edi=edi->next_edi;
        }
    }
}


static inline void ed_unshift_single_coord(matrix box, const rvec x, const ivec is, rvec xu)
{
    int tx,ty,tz;


    GMX_MPE_LOG(ev_unshift_start);

    tx=is[XX];
    ty=is[YY];
    tz=is[ZZ];

    if(TRICLINIC(box))
    {
        xu[XX] = x[XX]-tx*box[XX][XX]-ty*box[YY][XX]-tz*box[ZZ][XX];
        xu[YY] = x[YY]-ty*box[YY][YY]-tz*box[ZZ][YY];
        xu[ZZ] = x[ZZ]-tz*box[ZZ][ZZ];
    } else
    {
        xu[XX] = x[XX]-tx*box[XX][XX];
        xu[YY] = x[YY]-ty*box[YY][YY];
        xu[ZZ] = x[ZZ]-tz*box[ZZ][ZZ];
    }

    GMX_MPE_LOG(ev_unshift_finish);
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
        proj = projectx(edi, xcoll, edi->vecs.linfix.vec[i]);

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
        proj=projectx(edi, xcoll, edi->vecs.linacc.vec[i]);

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
        proj[i] = projectx(edi, xcoll, edi->vecs.radfix.vec[i]);
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
        proj[i] = projectx(edi, xcoll, edi->vecs.radacc.vec[i]);
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
        for (j=0; j<edi->sav.nr; j++)
        {
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
    int  i,j;
    real rad=0.0, ratio=0.0;
    struct t_do_radcon *loc;
    gmx_bool bFirst;
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
        loc->proj[i] = projectx(edi, xcoll, edi->vecs.radcon.vec[i]);
        rad += pow(loc->proj[i] - edi->vecs.radcon.refproj[i], 2);
    }
    rad = sqrt(rad);
    /* only correct when radius increased */
    if (rad > edi->vecs.radcon.radius)
    {
        ratio = edi->vecs.radcon.radius/rad - 1.0;

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
                rvec_inc(xcoll[j], vec_dum);
            }
        }
    }
    else
        edi->vecs.radcon.radius = rad;

    if (rad != edi->vecs.radcon.radius)
    {
        rad = 0.0;
        for (i=0; i<edi->vecs.radcon.neig; i++)
        {
            /* calculate the projections, radius */
            loc->proj[i] = projectx(edi, xcoll, edi->vecs.radcon.vec[i]);
            rad += pow(loc->proj[i] - edi->vecs.radcon.refproj[i], 2);
        }
        rad = sqrt(rad);
    }
}


static void ed_apply_constraints(rvec *xcoll, t_edpar *edi, gmx_large_int_t step, t_commrec *cr)
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


/* Write out the projections onto the eigenvectors */
static void write_edo(int nr_edi, t_edpar *edi, gmx_edsam_t ed, gmx_large_int_t step,real rmsd)
{
    int i;
    char buf[22];


    if (edi->bNeedDoEdsam)
    {
        if (step == -1)
            fprintf(ed->edo, "Initial projections:\n");
        else
        {
            fprintf(ed->edo,"Step %s, ED #%d  ", gmx_step_str(step, buf), nr_edi);
            fprintf(ed->edo,"  RMSD %f nm\n",rmsd);
        }

        if (edi->vecs.mon.neig)
        {
            fprintf(ed->edo,"  Monitor eigenvectors");
            for (i=0; i<edi->vecs.mon.neig; i++)
                fprintf(ed->edo," %d: %12.5e ",edi->vecs.mon.ieig[i],edi->vecs.mon.xproj[i]);
            fprintf(ed->edo,"\n");
        }
        if (edi->vecs.linfix.neig)
        {
            fprintf(ed->edo,"  Linfix  eigenvectors");
            for (i=0; i<edi->vecs.linfix.neig; i++)
                fprintf(ed->edo," %d: %12.5e ",edi->vecs.linfix.ieig[i],edi->vecs.linfix.xproj[i]);
            fprintf(ed->edo,"\n");
        }
        if (edi->vecs.linacc.neig)
        {
            fprintf(ed->edo,"  Linacc  eigenvectors");
            for (i=0; i<edi->vecs.linacc.neig; i++)
                fprintf(ed->edo," %d: %12.5e ",edi->vecs.linacc.ieig[i],edi->vecs.linacc.xproj[i]);
            fprintf(ed->edo,"\n");
        }
        if (edi->vecs.radfix.neig)
        {
            fprintf(ed->edo,"  Radfix  eigenvectors");
            for (i=0; i<edi->vecs.radfix.neig; i++)
                fprintf(ed->edo," %d: %12.5e ",edi->vecs.radfix.ieig[i],edi->vecs.radfix.xproj[i]);
            fprintf(ed->edo,"\n");
            fprintf(ed->edo,"  fixed increment radius = %f\n", calc_radius(&edi->vecs.radfix));
        }
        if (edi->vecs.radacc.neig)
        {
            fprintf(ed->edo,"  Radacc  eigenvectors");
            for (i=0; i<edi->vecs.radacc.neig; i++)
                fprintf(ed->edo," %d: %12.5e ",edi->vecs.radacc.ieig[i],edi->vecs.radacc.xproj[i]);
            fprintf(ed->edo,"\n");
            fprintf(ed->edo,"  acceptance radius      = %f\n", calc_radius(&edi->vecs.radacc));
        }
        if (edi->vecs.radcon.neig)
        {
            fprintf(ed->edo,"  Radcon  eigenvectors");
            for (i=0; i<edi->vecs.radcon.neig; i++)
                fprintf(ed->edo," %d: %12.5e ",edi->vecs.radcon.ieig[i],edi->vecs.radcon.xproj[i]);
            fprintf(ed->edo,"\n");
            fprintf(ed->edo,"  contracting radius     = %f\n", calc_radius(&edi->vecs.radcon));
        }
    }
}

/* Returns if any constraints are switched on */
static int ed_constraints(gmx_bool edtype, t_edpar *edi)
{
    if (edtype == eEDedsam || edtype == eEDflood)
    {
        return (edi->vecs.linfix.neig || edi->vecs.linacc.neig ||
                edi->vecs.radfix.neig || edi->vecs.radacc.neig ||
                edi->vecs.radcon.neig);
    }
    return 0;
}


/* Copies reference projection 'refproj' to fixed 'refproj0' variable for flooding/
 * umbrella sampling simulations. */
static void copyEvecReference(t_eigvec* floodvecs)
{
	int i;


	for (i=0; i<floodvecs->neig; i++)
	{
		floodvecs->refproj0[i] = floodvecs->refproj[i];
	}
}


void init_edsam(gmx_mtop_t  *mtop,   /* global topology                    */
                t_inputrec  *ir,     /* input record                       */
                t_commrec   *cr,     /* communication record               */
                gmx_edsam_t ed,      /* contains all ED data               */
                rvec        x[],     /* positions of the whole MD system   */
                matrix      box)     /* the box                            */
{
    t_edpar *edi = NULL;    /* points to a single edi data set */
    int     numedis=0;      /* keep track of the number of ED data sets in edi file */
    int     i,nr_edi;
    rvec    *x_pbc  = NULL; /* positions of the whole MD system with pbc removed  */
    rvec    *xfit   = NULL; /* the positions which will be fitted to the reference structure  */
    rvec    *xstart = NULL; /* the positions which are subject to ED sampling */
    rvec    fit_transvec;   /* translation ... */
    matrix  fit_rotmat;     /* ... and rotation from fit to reference structure */


    if (!DOMAINDECOMP(cr) && PAR(cr) && MASTER(cr))
        gmx_fatal(FARGS, "Please switch on domain decomposition to use essential dynamics in parallel.");

    GMX_MPE_LOG(ev_edsam_start);

    if (MASTER(cr))
        fprintf(stderr, "ED: Initializing essential dynamics constraints.\n");

    /* Needed for initializing radacc radius in do_edsam */
    ed->bFirst = 1;

    /* The input file is read by the master and the edi structures are
     * initialized here. Input is stored in ed->edpar. Then the edi
     * structures are transferred to the other nodes */
    if (MASTER(cr))
    {
        snew(ed->edpar,1);
        /* Read the whole edi file at once: */
        read_edi_file(ed,ed->edpar,mtop->natoms,cr);

        /* Initialization for every ED/flooding dataset. Flooding uses one edi dataset per
         * flooding vector, Essential dynamics can be applied to more than one structure
         * as well, but will be done in the order given in the edi file, so
         * expect different results for different order of edi file concatenation! */
        edi=ed->edpar;
        while(edi != NULL)
        {
            init_edi(mtop,ir,cr,ed,edi);

            /* Init flooding parameters if needed */
            init_flood(edi,ed,ir->delta_t,cr);

            edi=edi->next_edi;
            numedis++;
        }
    }

    /* The master does the work here. The other nodes get the positions
     * not before dd_partition_system which is called after init_edsam */
    if (MASTER(cr))
    {
        /* Remove pbc, make molecule whole.
         * When ir->bContinuation=TRUE this has already been done, but ok.
         */
        snew(x_pbc,mtop->natoms);
        m_rveccopy(mtop->natoms,x,x_pbc);
        do_pbc_first_mtop(NULL,ir->ePBC,box,mtop,x_pbc);

        /* Reset pointer to first ED data set which contains the actual ED data */
        edi=ed->edpar;

        /* Loop over all ED/flooding data sets (usually only one, though) */
        for (nr_edi = 1; nr_edi <= numedis; nr_edi++)
        {
            /* We use srenew to allocate memory since the size of the buffers
             * is likely to change with every ED dataset */
            srenew(xfit  , edi->sref.nr );
            srenew(xstart, edi->sav.nr  );

            /* Extract the positions of the atoms to which will be fitted */
            for (i=0; i < edi->sref.nr; i++)
            {
                copy_rvec(x_pbc[edi->sref.anrs[i]], xfit[i]);

                /* Save the sref positions such that in the next time step the molecule can
                 * be made whole again (in the parallel case) */
                if (PAR(cr))
                    copy_rvec(xfit[i], edi->sref.x_old[i]);
            }

            /* Extract the positions of the atoms subject to ED sampling */
            for (i=0; i < edi->sav.nr; i++)
            {
                copy_rvec(x_pbc[edi->sav.anrs[i]], xstart[i]);

                /* Save the sav positions such that in the next time step the molecule can
                 * be made whole again (in the parallel case) */
                if (PAR(cr))
                    copy_rvec(xstart[i], edi->sav.x_old[i]);
            }

            /* Make the fit to the REFERENCE structure, get translation and rotation */
            fit_to_reference(xfit, fit_transvec, fit_rotmat, edi);

            /* Output how well we fit to the reference at the start */
            translate_and_rotate(xfit, edi->sref.nr, fit_transvec, fit_rotmat);
            fprintf(stderr, "ED: Initial RMSD from reference after fit = %f nm (dataset #%d)\n",
                    rmsd_from_structure(xfit, &edi->sref), nr_edi);

            /* Now apply the translation and rotation to the atoms on which ED sampling will be performed */
            translate_and_rotate(xstart, edi->sav.nr, fit_transvec, fit_rotmat);

            /* calculate initial projections */
            project(xstart, edi);

            /* process target structure, if required */
            if (edi->star.nr > 0)
            {
                fprintf(stderr, "ED: Fitting target structure to reference structure\n");
                /* get translation & rotation for fit of target structure to reference structure */
                fit_to_reference(edi->star.x, fit_transvec, fit_rotmat, edi);
                /* do the fit */
                translate_and_rotate(edi->star.x, edi->sav.nr, fit_transvec, fit_rotmat);
                rad_project(edi, edi->star.x, &edi->vecs.radcon, cr);
            } else
                rad_project(edi, xstart, &edi->vecs.radcon, cr);

            /* process structure that will serve as origin of expansion circle */
            if ( (eEDflood == ed->eEDtype) && (FALSE == edi->flood.bConstForce) )
                fprintf(stderr, "ED: Setting center of flooding potential (0 = average structure)\n");
            if (edi->sori.nr > 0)
            {
                fprintf(stderr, "ED: Fitting origin structure to reference structure\n");
                /* fit this structure to reference structure */
                fit_to_reference(edi->sori.x, fit_transvec, fit_rotmat, edi);
                /* do the fit */
                translate_and_rotate(edi->sori.x, edi->sav.nr, fit_transvec, fit_rotmat);
                rad_project(edi, edi->sori.x, &edi->vecs.radacc, cr);
                rad_project(edi, edi->sori.x, &edi->vecs.radfix, cr);
                if ( (eEDflood == ed->eEDtype) && (FALSE == edi->flood.bConstForce) )
                {
                    fprintf(stderr, "ED: The ORIGIN structure will define the flooding potential center.\n");
                    /* Set center of flooding potential to the ORIGIN structure */
                    rad_project(edi, edi->sori.x, &edi->flood.vecs, cr);
                    /* We already know that no (moving) reference position was provided,
                     * therefore we can overwrite refproj[0]*/
                    copyEvecReference(&edi->flood.vecs);
                }
            }
            else /* No origin structure given */
            {
                rad_project(edi, xstart, &edi->vecs.radacc, cr);
                rad_project(edi, xstart, &edi->vecs.radfix, cr);
                if ( (eEDflood == ed->eEDtype) && (FALSE == edi->flood.bConstForce) )
                {
                    if (edi->flood.bHarmonic)
                    {
                        fprintf(stderr, "ED: A (possibly changing) ref. projection will define the flooding potential center.\n");
                        for (i=0; i<edi->flood.vecs.neig; i++)
                            edi->flood.vecs.refproj[i] = edi->flood.vecs.refproj0[i];
                    }
                    else
                    {
                        fprintf(stderr, "ED: The AVERAGE structure will define the flooding potential center.\n");
                        /* Set center of flooding potential to the center of the covariance matrix,
                         * i.e. the average structure, i.e. zero in the projected system */
                        for (i=0; i<edi->flood.vecs.neig; i++)
                            edi->flood.vecs.refproj[i] = 0.0;
                    }
                }
            }
            /* For convenience, output the center of the flooding potential for the eigenvectors */
            if ( (eEDflood == ed->eEDtype) && (FALSE == edi->flood.bConstForce) )
            {
                for (i=0; i<edi->flood.vecs.neig; i++)
                {
                    fprintf(stdout, "ED: EV %d flooding potential center: %11.4e", i, edi->flood.vecs.refproj[i]);
                    if (edi->flood.bHarmonic)
                        fprintf(stdout, " (adding %11.4e/timestep)", edi->flood.vecs.refprojslope[i]);
                    fprintf(stdout, "\n");
                }
            }

            /* set starting projections for linsam */
            rad_project(edi, xstart, &edi->vecs.linacc, cr);
            rad_project(edi, xstart, &edi->vecs.linfix, cr);

            /* Output to file, set the step to -1 so that write_edo knows it was called from init_edsam */
            if (ed->edo && !(ed->bStartFromCpt))
                write_edo(nr_edi, edi, ed, -1, 0);

            /* Prepare for the next edi data set: */
            edi=edi->next_edi;
        }
        /* Cleaning up on the master node: */
        sfree(x_pbc);
        sfree(xfit);
        sfree(xstart);

    } /* end of MASTER only section */

    if (PAR(cr))
    {
        /* First let everybody know how many ED data sets to expect */
        gmx_bcast(sizeof(numedis), &numedis, cr);
        /* Broadcast the essential dynamics / flooding data to all nodes */
        broadcast_ed_data(cr, ed, numedis);
    }
    else
    {
        /* In the single-CPU case, point the local atom numbers pointers to the global
         * one, so that we can use the same notation in serial and parallel case: */

        /* Loop over all ED data sets (usually only one, though) */
        edi=ed->edpar;
        for (nr_edi = 1; nr_edi <= numedis; nr_edi++)
        {
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

            /* An on we go to the next edi dataset */
            edi=edi->next_edi;
        }
    }

    /* Allocate space for ED buffer variables */
    /* Again, loop over ED data sets */
    edi=ed->edpar;
    for (nr_edi = 1; nr_edi <= numedis; nr_edi++)
    {
        /* Allocate space for ED buffer */
        snew(edi->buf, 1);
        snew(edi->buf->do_edsam, 1);

        /* Space for collective ED buffer variables */

        /* Collective positions of atoms with the average indices */
        snew(edi->buf->do_edsam->xcoll                  , edi->sav.nr);
        snew(edi->buf->do_edsam->shifts_xcoll           , edi->sav.nr); /* buffer for xcoll shifts */
        snew(edi->buf->do_edsam->extra_shifts_xcoll     , edi->sav.nr);
        /* Collective positions of atoms with the reference indices */
        if (!edi->bRefEqAv)
        {
            snew(edi->buf->do_edsam->xc_ref             , edi->sref.nr);
            snew(edi->buf->do_edsam->shifts_xc_ref      , edi->sref.nr); /* To store the shifts in */
            snew(edi->buf->do_edsam->extra_shifts_xc_ref, edi->sref.nr);
        }

        /* Get memory for flooding forces */
        snew(edi->flood.forces_cartesian                , edi->sav.nr);

#ifdef DUMPEDI
        /* Dump it all into one file per process */
        dump_edi(edi, cr, nr_edi);
#endif

        /* An on we go to the next edi dataset */
        edi=edi->next_edi;
    }

    /* Flush the edo file so that the user can check some things
     * when the simulation has started */
    if (ed->edo)
        fflush(ed->edo);

    GMX_MPE_LOG(ev_edsam_finish);
}


void do_edsam(t_inputrec  *ir,
              gmx_large_int_t step,
              t_mdatoms   *md,
              t_commrec   *cr,
              rvec        xs[],   /* The local current positions on this processor */
              rvec        v[],    /* The velocities */
              matrix      box,
              gmx_edsam_t ed)
{
    int     i,edinr,iupdate=500;
    matrix  rotmat;         /* rotation matrix */
    rvec    transvec;       /* translation vector */
    rvec    dv,dx,x_unsh;   /* tmp vectors for velocity, distance, unshifted x coordinate */
    real    dt_1;           /* 1/dt */
    struct t_do_edsam *buf;
    t_edpar *edi;
    real    rmsdev=-1;      /* RMSD from reference structure prior to applying the constraints */
    gmx_bool bSuppress=FALSE; /* Write .edo file on master? */


    /* Check if ED sampling has to be performed */
    if ( ed->eEDtype==eEDnone )
        return;

    /* Suppress output on first call of do_edsam if
     * two-step sd2 integrator is used */
    if ( (ir->eI==eiSD2) && (v != NULL) )
        bSuppress = TRUE;

    dt_1 = 1.0/ir->delta_t;

    /* Loop over all ED datasets (usually one) */
    edi  = ed->edpar;
    edinr = 0;
    while (edi != NULL)
    {
        edinr++;
        if (edi->bNeedDoEdsam)
        {

            buf=edi->buf->do_edsam;

            if (ed->bFirst)
                /* initialise radacc radius for slope criterion */
                buf->oldrad=calc_radius(&edi->vecs.radacc);

            /* Copy the positions into buf->xc* arrays and after ED
             * feed back corrections to the official positions */

            /* Broadcast the ED positions such that every node has all of them
             * Every node contributes its local positions xs and stores it in
             * the collective buf->xcoll array. Note that for edinr > 1
             * xs could already have been modified by an earlier ED */

            communicate_group_positions(cr, buf->xcoll, buf->shifts_xcoll, buf->extra_shifts_xcoll, buf->bUpdateShifts, xs,
                    edi->sav.nr, edi->sav.nr_loc, edi->sav.anrs_loc, edi->sav.c_ind, edi->sav.x_old,  box);

#ifdef DEBUG_ED
            dump_xcoll(edi, buf, cr, step);
#endif
            /* Only assembly reference positions if their indices differ from the average ones */
            if (!edi->bRefEqAv)
                communicate_group_positions(cr, buf->xc_ref, buf->shifts_xc_ref, buf->extra_shifts_xc_ref, buf->bUpdateShifts, xs,
                        edi->sref.nr, edi->sref.nr_loc, edi->sref.anrs_loc, edi->sref.c_ind, edi->sref.x_old, box);

            /* If bUpdateShifts was TRUE then the shifts have just been updated in get_positions.
             * We do not need to uptdate the shifts until the next NS step */
            buf->bUpdateShifts = FALSE;

            /* Now all nodes have all of the ED positions in edi->sav->xcoll,
             * as well as the indices in edi->sav.anrs */

            /* Fit the reference indices to the reference structure */
            if (edi->bRefEqAv)
                fit_to_reference(buf->xcoll , transvec, rotmat, edi);
            else
                fit_to_reference(buf->xc_ref, transvec, rotmat, edi);

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

            /* update radsam references, when required */
            if (do_per_step(step,edi->maxedsteps) && step >= edi->presteps)
            {
                project(buf->xcoll, edi);
                rad_project(edi, buf->xcoll, &edi->vecs.radacc, cr);
                rad_project(edi, buf->xcoll, &edi->vecs.radfix, cr);
                buf->oldrad=-1.e5;
            }

            /* update radacc references, when required */
            if (do_per_step(step,iupdate) && step >= edi->presteps)
            {
                edi->vecs.radacc.radius = calc_radius(&edi->vecs.radacc);
                if (edi->vecs.radacc.radius - buf->oldrad < edi->slope)
                {
                    project(buf->xcoll, edi);
                    rad_project(edi, buf->xcoll, &edi->vecs.radacc, cr);
                    buf->oldrad = 0.0;
                } else
                    buf->oldrad = edi->vecs.radacc.radius;
            }

            /* apply the constraints */
            if (step >= edi->presteps && ed_constraints(ed->eEDtype, edi))
            {
                /* ED constraints should be applied already in the first MD step
                 * (which is step 0), therefore we pass step+1 to the routine */
                ed_apply_constraints(buf->xcoll, edi, step+1 - ir->init_step, cr);
            }

            /* write to edo, when required */
            if (do_per_step(step,edi->outfrq))
            {
                project(buf->xcoll, edi);
                if (MASTER(cr) && !bSuppress)
                    write_edo(edinr, edi, ed, step, rmsdev);
            }

            /* Copy back the positions unless monitoring only */
            if (ed_constraints(ed->eEDtype, edi))
            {
                /* remove fitting */
                rmfit(edi->sav.nr, buf->xcoll, transvec, rotmat);

                /* Copy the ED corrected positions into the coordinate array */
                /* Each node copies its local part. In the serial case, nat_loc is the
                 * total number of ED atoms */
                for (i=0; i<edi->sav.nr_loc; i++)
                {
                    /* Unshift local ED coordinate and store in x_unsh */
                    ed_unshift_single_coord(box, buf->xcoll[edi->sav.c_ind[i]],
                                            buf->shifts_xcoll[edi->sav.c_ind[i]], x_unsh);

                    /* dx is the ED correction to the positions: */
                    rvec_sub(x_unsh, xs[edi->sav.anrs_loc[i]], dx);

                    if (v != NULL)
                    {
                        /* dv is the ED correction to the velocity: */
                        svmul(dt_1, dx, dv);
                        /* apply the velocity correction: */
                        rvec_inc(v[edi->sav.anrs_loc[i]], dv);
                    }
                    /* Finally apply the position correction due to ED: */
                    copy_rvec(x_unsh, xs[edi->sav.anrs_loc[i]]);
                }
            }
        } /* END of if (edi->bNeedDoEdsam) */

        /* Prepare for the next ED dataset */
        edi = edi->next_edi;

    } /* END of loop over ED datasets */

    ed->bFirst = FALSE;

    GMX_MPE_LOG(ev_edsam_finish);
}
