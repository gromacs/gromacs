/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "domdec.h"
#include "gmx_wallcycle.h"
#include "trnio.h"
#include "smalloc.h"
#include "network.h"
#include "pbc.h"
#include "futil.h"
#include "mdrun.h"
#include "txtdump.h"
#include "names.h"
#include "mtop_util.h"
#include "names.h"
#include "nrjac.h"
#include "vec.h"
#include "gmx_ga2la.h"
#include "xvgr.h"
#include "gmxfio.h"
#include "mpelogging.h"
#include "groupcoord.h"



/* Helper structure for sorting positions along rotation vector             */
typedef struct {
    real xcproj;            /* Projection of xc on the rotation vector        */
    int ind;                /* Index of xc                                    */
} sort_along_vec_t;


/* Enforced rotation / flexible: determine the angle of each slab             */
typedef struct gmx_slabdata
{
    int  nat;               /* Number of atoms belonging to this slab         */
    rvec *x;                /* The positions belonging to this slab. In 
                               general, this should be all positions of the 
                               whole rotation group, but we leave those away 
                               that have a small enough weight                */
    rvec *ref;              /* Same for reference                             */
    real *weight;           /* The weight for each atom                       */
} t_gmx_slabdata;


/* Enforced rotation data for all groups                                      */
typedef struct gmx_enfrot
{
    FILE  *out_rot;         /* Output file for rotation data                  */
    FILE  *out_torque;      /* Output file for torque data                    */
    FILE  *out_angles;      /* Output file for slab angles for flexible type  */
    FILE  *out_slabs;       /* Output file for slab centers                   */
    int   bufsize;          /* Allocation size of buf                         */
    rvec  *buf;             /* Coordinate buffer variable                     */
    sort_along_vec_t *data; /* Buffer variable needed for position sorting    */
    real  *mpi_inbuf;       /* MPI buffer                                     */
    real  *mpi_outbuf;      /* MPI buffer                                     */
    int   mpi_bufsize;      /* Allocation size of in & outbuf                 */
    real  Vrot;             /* (Local) part of the enf. rotation potential    */
} t_gmx_enfrot;


/* Global enforced rotation data for a single rotation group                  */
typedef struct gmx_enfrotgrp
{
    real    degangle;       /* Rotation angle in degree                       */
    matrix  rotmat;         /* Rotation matrix                                */
    atom_id *ind_loc;       /* Local rotation indices                         */
    int     nat_loc;        /* Number of local group atoms                    */
    int     nalloc_loc;     /* Allocation size for ind_loc and weight_loc     */
    rvec    center;         /* Center of the rotation group positions if
                             * chose origin is COG or COM, otherwise (0,0,0)  */

    /* Collective coordinates for the whole rotation group */
    rvec  *xc_ref;          /* Reference (unrotated) positions                */
    real  *xc_ref_length;   /* Length of each x_rotref vector after x_rotref 
                               has been put into origin                       */
    int   *xc_ref_ind;      /* Position of each local atom in the collective
                               array                                          */
    rvec  xc_ref_center;    /* Center of the reference positions. May be 
                               mass-weighted                                  */
    rvec  *xc;              /* Current (collective) positions                 */
    ivec  *xc_shifts;       /* Current (collective) shifts                    */
    ivec  *xc_eshifts;      /* Extra shifts since last DD step                */
    rvec  *xc_old;          /* Old (collective) positions                     */
    rvec  *xc_norm;         /* Normalized form of the current positions       */
    rvec  *xc_ref_sorted;   /* Reference positions (sorted in the same order 
                               as xc when sorted)                             */
    int   *xc_sortind;      /* Indices of sorted positions                    */
    real  *mc;              /* Collective masses                              */
    /* Fixed rotation only */
    rvec  *xr_loc;          /* Local reference coords, correctly rotated      */
    rvec  *x_loc_pbc;       /* Local current coords, correct PBC image        */
    real  *m_loc;           /* Masses of the current local atoms              */
    real  fix_torque_v;     /* Torque in the direction of rotation vector     */
    real  fix_angles_v;
    real  fix_weight_v;
    /* Flexible rotation only */
    int   nslabs_alloc;     /* For this many slabs memory is allocated        */
    int   slab_first;       /* Lowermost slab for that the calculation needs 
                               to be performed at a given time step           */
    int   slab_last;        /* Uppermost slab ...                             */
    int   slab_first_ref;   /* First slab for which reference COG is stored   */
    int   slab_last_ref;    /* Last ...                                       */
    int   slab_buffer;      /* Slab buffer region around reference slabs      */
    int   *firstatom;       /* First relevant atom for a slab                 */
    int   *lastatom;        /* Last relevant atom for a slab                  */
    rvec  *slab_center;     /* Gaussian-weighted slab center (COG)            */
    rvec  *slab_center_ref; /* Gaussian-weighted slab COG for the 
                               reference positions                            */
    real  *slab_weights;    /* Sum of gaussian weights in a slab              */
    real  *slab_torque_v;   /* Torque T = r x f for each slab.                */
                            /* torque_v = m.v = angular momentum in the 
                               direction of v                                 */
    real  max_beta;         /* min_gaussian from inputrec->rotgrp is the
                               minimum value the gaussian must have so that 
                               the force is actually evaluated max_beta is 
                               just another way to put it                     */
    real  *gn_atom;         /* Precalculated gaussians for a single atom      */
    int   *gn_slabind;      /* Tells to which slab each precalculated gaussian 
                               belongs                                        */
    real  V;                /* Rotation potential for this rotation group     */
    rvec  *f_rot_loc;       /* Array to store the forces on the local atoms 
                               resulting from enforced rotation potential     */
    t_gmx_slabdata *slab_data; /* Holds atom positions and gaussian weights 
                               of atoms belonging to a slab                   */
} t_gmx_enfrotgrp;


static double** allocate_square_matrix(int dim)
{
    int i;
    double** mat = NULL; 
    
    
    snew(mat, dim);
    for(i=0; i<dim; i++)
        snew(mat[i], dim);

    return mat;
}


static void free_square_matrix(double** mat, int dim)
{
    int i;
    
    
    for (i=0; i<dim; i++)
        sfree(mat[i]);
    sfree(mat);
}


/* Output rotation energy and torque for each rotation group */
static void reduce_output(t_commrec *cr, t_rot *rot, real t)
{
    int      g,i,islab,nslabs=0;
    int      count;      /* MPI element counter                               */
    t_rotgrp *rotg;
    gmx_enfrot_t er;     /* Pointer to the enforced rotation buffer variables */
    gmx_enfrotgrp_t erg; /* Pointer to enforced rotation group data           */

    
    er=rot->enfrot;
    
    /* Fill the MPI buffer with stuff to reduce: */
    if (PAR(cr))
    {
        count=0;
        for (g=0; g < rot->ngrp; g++)
        {
            rotg = &rot->grp[g];
            erg=rotg->enfrotgrp;
            nslabs = erg->slab_last - erg->slab_first + 1;
            er->mpi_inbuf[count++] = erg->V;
            switch (rotg->eType)
            {
            case erotgFIXED:
            case erotgFIXED_PLANE:
                er->mpi_inbuf[count++] = erg->fix_torque_v;
                er->mpi_inbuf[count++] = erg->fix_angles_v;
                er->mpi_inbuf[count++] = erg->fix_weight_v;
                break;
            case erotgFLEX1:
            case erotgFLEX2:
                /* (Re-)allocate memory for MPI buffer: */
                if (er->mpi_bufsize < count+nslabs)
                {
                    er->mpi_bufsize = count+nslabs;
                    srenew(er->mpi_inbuf , er->mpi_bufsize);
                    srenew(er->mpi_outbuf, er->mpi_bufsize);
                }
                for (i=0; i<nslabs; i++)
                    er->mpi_inbuf[count++] = erg->slab_torque_v[i];
                break;
            default:
                break;
            }
        }
#ifdef GMX_MPI
        MPI_Reduce(er->mpi_inbuf, er->mpi_outbuf, count, GMX_MPI_REAL, MPI_SUM, MASTERRANK(cr), cr->mpi_comm_mygroup);
#endif
        /* Copy back the reduced data from the buffer on the master */
        if (MASTER(cr))
        {
            count=0;
            for (g=0; g < rot->ngrp; g++)
            {
                rotg = &rot->grp[g];
                erg=rotg->enfrotgrp;
                nslabs = erg->slab_last - erg->slab_first + 1;
                erg->V = er->mpi_outbuf[count++];
                switch (rotg->eType)
                {
                case erotgFIXED:
                case erotgFIXED_PLANE:
                    erg->fix_torque_v = er->mpi_outbuf[count++];
                    erg->fix_angles_v = er->mpi_outbuf[count++];
                    erg->fix_weight_v = er->mpi_outbuf[count++];
                    break;
                case erotgFLEX1:
                case erotgFLEX2:
                    for (i=0; i<nslabs; i++)
                        erg->slab_torque_v[i] = er->mpi_outbuf[count++];
                    break;
                default:
                    break;
                }
            }
        }
    }
    
    /* Output */
    if (MASTER(cr))
    {
        /* Av. angle and total torque for each rotation group */
        for (g=0; g < rot->ngrp; g++)
        {
            rotg=&rot->grp[g];
            erg=rotg->enfrotgrp;
            
            /* Output to main rotation log file: */
            if (rotg->eType == erotgFIXED || rotg->eType == erotgFIXED_PLANE)
            {
                fprintf(er->out_rot, "%12.4f%12.3e", 
                        (erg->fix_angles_v/erg->fix_weight_v)*180.0*M_1_PI,
                        erg->fix_torque_v);
            }
            fprintf(er->out_rot, "%12.3e", erg->V);
                        
            /* Output to torque log file: */
            if (rotg->eType == erotgFLEX1 || rotg->eType == erotgFLEX2)
            {
                fprintf(er->out_torque, "%12.3e%6d", t, g);
                for (i=erg->slab_first; i<=erg->slab_last; i++)
                {
                    islab = i - erg->slab_first;  /* slab index */
                    /* Only output if enough weight is in slab */
                    if (erg->slab_weights[islab] > rotg->min_gaussian)
                        fprintf(er->out_torque, "%6d%12.3e", i, erg->slab_torque_v[islab]);
                }
                fprintf(er->out_torque , "\n");
            }
        }
        fprintf(er->out_rot, "\n");
    }
}


/* Add the forces from enforced rotation potential to the local forces.
 * Should be called after the SR forces have been evaluated */
extern real add_rot_forces(t_rot *rot, rvec f[], t_commrec *cr, int step, real t)
{
    int g,l,ii;
    t_rotgrp *rotg;
    gmx_enfrot_t er;     /* Pointer to the enforced rotation buffer variables */
    gmx_enfrotgrp_t erg; /* Pointer to enforced rotation group data           */

    
    er=rot->enfrot;
    
    GMX_MPE_LOG(ev_add_rot_forces_start);
    
    /* Reduce energy,torque, angles etc. to get the sum values (per rotation group) 
     * on the master and output these values to file. */
    if (do_per_step(step, rot->nsttout))
        reduce_output(cr, rot, t);

    /* Total rotation potential is the sum over all rotation groups */
    er->Vrot = 0.0; 
        
    /* Loop over enforced rotation groups (usually 1, though)
     * Apply the forces from rotation potentials */
    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        erg=rotg->enfrotgrp;
        er->Vrot += erg->V;
        for (l=0; l<erg->nat_loc; l++)
        {
            /* Get the right index of the local force */
            ii = erg->ind_loc[l];
            /* Add */
            rvec_inc(f[ii],erg->f_rot_loc[l]);
        }
    }
    
    GMX_MPE_LOG(ev_add_rot_forces_finish);

    return (MASTER(cr)? er->Vrot : 0.0);
}


/* Calculate the maximum beta that leads to a gaussian larger min_gaussian,
 * also does some checks
 */
static double calc_beta_max(real min_gaussian, real slab_dist)
{
    const double norm = 0.5698457353514458216;  /* = 1/1.7548609 */
    double sigma;
    double arg;
    
    
    if (slab_dist <= 0)
        gmx_fatal(FARGS, "Slab distance of flexible rotation groups must be >=0 !");
    
    /* Define the sigma value */
    sigma = 0.7*slab_dist;

    /* Calculate the argument for the logarithm and check that the log() result is negative or 0 */
    arg = min_gaussian/norm;
    if (arg > 1.0)
        gmx_fatal(FARGS, "min_gaussian of flexible rotation groups must be <%g", norm);
    
    return sqrt(-2.0*sigma*sigma*log(min_gaussian/norm));
}


static inline real calc_beta(rvec curr_x, t_rotgrp *rotg, int n)
{
    return iprod(curr_x, rotg->vec) - rotg->slab_dist * n;
}


static inline real gaussian_weight(rvec curr_x, t_rotgrp *rotg, int n)
{
    /* norm is chosen such that the sum of the gaussians
     * over the slabs is approximately 1.0 everywhere */
    const real norm = 0.5698457353514458216;  /* = 1/1.7548609 */
    real       sigma;

    
    /* Define the sigma value */
    sigma = 0.7*rotg->slab_dist;
    /* Calculate the Gaussian value of slab n for position curr_x */
    return norm * exp( -0.5 * sqr( calc_beta(curr_x, rotg, n)/sigma ) );
}


static void get_slab_centers(
        t_rotgrp *rotg,       /* The rotation group information               */
        rvec      *xc,        /* The rotation group positions; will 
                                 typically be enfrotgrp->xc, but at first call 
                                 it is enfrotgrp->xc_ref                      */
        t_commrec *cr,        /* Communication record                         */
        int       g,          /* The number of the rotation group             */
        real      time,       /* Used for output only                         */
        FILE      *out_slabs, /* For outputting COG per slab information      */
        bool      bOutStep,   /* Is this an output step?                      */
        bool      bReference) /* If this routine is called from
                                 init_rot_group we need to store
                                 the reference slab COGs                      */
{
    rvec curr_x;              /* The position of an atom                 */
    rvec curr_x_weighted;     /* The gaussian-weighted position          */
    real gaussian;            /* The gaussian weight                     */
    int i,j,islab;
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;

    /* Loop over slabs */
    for (j = erg->slab_first; j <= erg->slab_last; j++)
    {
        islab = j - erg->slab_first;
        /* Initialize data for this slab: */
        clear_rvec(erg->slab_center[islab]);
        erg->slab_weights[islab] = 0.0;

        /* loop over all atoms in the rotation group */
        for(i=0; i<rotg->nat;i++)
        {
            copy_rvec(xc[i], curr_x);
            gaussian = gaussian_weight(curr_x, rotg, j);
            svmul(gaussian, curr_x, curr_x_weighted);
            rvec_add(erg->slab_center[islab], curr_x_weighted, erg->slab_center[islab]);
            erg->slab_weights[islab] += gaussian;
        } /* END of loop over rotation group atoms */

        /* Do the calculations ONLY if there is weight in the slab! */
        if (erg->slab_weights[islab] > 0.0)
        {
            svmul(1.0/erg->slab_weights[islab], erg->slab_center[islab], erg->slab_center[islab]);
        }
        else
        {
            /* We need to check this here, since we divide through slab_weights
             * in the flexible low-level routines! */
            gmx_fatal(FARGS, "Slab %d has weight 0, cannot determine its center!", j);
        }
        
        /* At first time step: save the COGs of the reference structure */
        if (bReference)
            copy_rvec(erg->slab_center[islab], erg->slab_center_ref[islab]);
    } /* END of loop over slabs */
    
    /* Output on the master */
    if (MASTER(cr) && bOutStep)
    {
        fprintf(out_slabs, "%12.3e%6d", time, g);
        for (j = erg->slab_first; j <= erg->slab_last; j++)
        {
            islab = j - erg->slab_first;
            fprintf(out_slabs, "%6d%12.3e%12.3e%12.3e",
                    j,erg->slab_center[islab][XX],erg->slab_center[islab][YY],erg->slab_center[islab][ZZ]);
        }
        fprintf(out_slabs, "\n");
    }
}


static void calc_rotmat(
        rvec vec,
        real degangle,  /* Angle alpha of rotation at time t in degrees       */
        matrix rotmat)  /* Rotation matrix                                    */
{
    real radangle;            /* Rotation angle in radians */
    real cosa;                /* cosine alpha              */
    real sina;                /* sine alpha                */
    real OMcosa;              /* 1 - cos(alpha)            */
    real dumxy, dumxz, dumyz; /* save computations         */
    rvec rot_vec;             /* Rotate around rot_vec ... */


    radangle = degangle * M_PI/180.0;
    copy_rvec(vec , rot_vec );

    /* Precompute some variables: */
    cosa   = cos(radangle);
    sina   = sin(radangle);
    OMcosa = 1.0 - cosa;
    dumxy  = rot_vec[XX]*rot_vec[YY]*OMcosa;
    dumxz  = rot_vec[XX]*rot_vec[ZZ]*OMcosa;
    dumyz  = rot_vec[YY]*rot_vec[ZZ]*OMcosa;

    /* Construct the rotation matrix for this rotation group: */
    /* 1st column: */
    rotmat[XX][XX] = cosa  + rot_vec[XX]*rot_vec[XX]*OMcosa;
    rotmat[YY][XX] = dumxy + rot_vec[ZZ]*sina;
    rotmat[ZZ][XX] = dumxz - rot_vec[YY]*sina;
    /* 2nd column: */
    rotmat[XX][YY] = dumxy - rot_vec[ZZ]*sina;
    rotmat[YY][YY] = cosa  + rot_vec[YY]*rot_vec[YY]*OMcosa;
    rotmat[ZZ][YY] = dumyz + rot_vec[XX]*sina;
    /* 3rd column: */
    rotmat[XX][ZZ] = dumxz + rot_vec[YY]*sina;
    rotmat[YY][ZZ] = dumyz - rot_vec[XX]*sina;
    rotmat[ZZ][ZZ] = cosa  + rot_vec[ZZ]*rot_vec[ZZ]*OMcosa;

#ifdef PRINTMATRIX
    int iii,jjj;

    for (iii=0; iii<3; iii++) {
        for (jjj=0; jjj<3; jjj++)
            fprintf(stderr, " %10.8f ",  rotmat[iii][jjj]);
        fprintf(stderr, "\n");
    }
#endif
}


/* Calculates torque on the rotation axis tau = coord x force */
static inline real torque(
        rvec rotvec,  /* rotation vector; MUST be normalized!                 */
        rvec force,   /* force                                                */
        rvec x,       /* position of atom on which the force acts             */
        rvec offset)  /* piercing point of rotation axis 
                         (or center of the slab for the flexible types)       */
{
    rvec vectmp, tau;

    
    /* Subtract offset */
    rvec_sub(x,offset,vectmp);
    
    /* coord x force */
    cprod(vectmp, force, tau);
    
    /* Return the part of the torque which is parallel to the rotation vector */
    return iprod(tau, rotvec);
}


/* Right-aligned output of value with standard width */
static void print_aligned(FILE *fp, char *str)
{
    fprintf(fp, "%12s", str);
}


/* Right-aligned output of value with standard short width */
static void print_aligned_short(FILE *fp, char *str)
{
    fprintf(fp, "%6s", str);
}


/* Right-aligned output of value with standard width */
static void print_aligned_group(FILE *fp, char *str, int g)
{
    char sbuf[STRLEN];
    
    
    sprintf(sbuf, "%s%d", str, g);
    fprintf(fp, "%12s", sbuf);
}


static FILE *open_output_file(const char *fn, int steps)
{
    FILE *fp;
    
    
    fp = ffopen(fn, "w");

    fprintf(fp, "# Output is written every %d time steps.\n\n", steps);
    
    return fp;
}


/* Open output file for slab COG data. Call on master only */
static FILE *open_slab_out(t_rot *rot, const char *fn)
{
    FILE      *fp=NULL;
    int       g;
    t_rotgrp  *rotg;


    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (rotg->eType == erotgFLEX1 || rotg->eType == erotgFLEX2)
        {
            if (NULL == fp)
                fp = open_output_file(fn, rot->nsttout);
            fprintf(fp, "# Rotation group %d (%s), slab distance %f nm\n", g, erotg_names[rotg->eType], rotg->slab_dist);
        }
    }
    
    if (fp != NULL)
    {
        fprintf(fp, "# The following columns will have the syntax: (COG = center of geometry, gaussian weighted)\n");
        fprintf(fp, "#     ");
        print_aligned_short(fp, "t");
        print_aligned_short(fp, "grp");
        print_aligned_short(fp, "slab");
        print_aligned(fp, "COG-X");
        print_aligned(fp, "COG-Y");
        print_aligned(fp, "COG-Z");
        print_aligned_short(fp, "slab");
        print_aligned(fp, "COG-X");
        print_aligned(fp, "COG-Y");
        print_aligned(fp, "COG-Z");
        print_aligned_short(fp, "slab");
        fprintf(fp, " ...\n");
        fflush(fp);
    }
    
    return fp;
}


/* Open output file and print some general information about the rotation groups.
 * Call on master only */
static FILE *open_rot_out(const char *fn, t_rot *rot, const output_env_t oenv, 
                          unsigned long Flags)
{
    FILE      *fp;
    int       g,nsets;
    t_rotgrp  *rotg;
    char      **setname,buf[50];
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */

    
    if (Flags & MD_APPENDFILES)
    {
        fp = gmx_fio_fopen(fn,"a");
    }
    else
    {
        fp = xvgropen(fn, "Rotation angles and energy", "Time (ps)", "angles and energies", oenv);
        fprintf(fp, "# The scalar tau is the torque in the direction of the rotation vector v.\n");
        fprintf(fp, "# To obtain the vectorial torque, multiply tau with the group's rot_vec.\n");
        fprintf(fp, "# Torques are given in [kJ/mol], angles in degrees, time in ps.\n");
        
        for (g=0; g<rot->ngrp; g++)
        {
            rotg = &rot->grp[g];
            erg=rotg->enfrotgrp;

            fprintf(fp, "# Rotation group %d (%s):\n", g, erotg_names[rotg->eType]);
            fprintf(fp, "# rot_vec%d            %12.5e %12.5e %12.5e\n", g, rotg->vec[XX], rotg->vec[YY], rotg->vec[ZZ]);
            fprintf(fp, "# rot_rate%d           %12.5e degree/ps\n",     g, rotg->rate);
            fprintf(fp, "# rot_k%d              %12.5e kJ/(mol*nm^2)\n", g, rotg->k);
            fprintf(fp, "# rot_origin%d          %s\n",                  g, erotg_originnames[rotg->eOrigin]);
            
            switch (rotg->eType)
            {
            case erotgFIXED:
            case erotgFIXED_PLANE:
                fprintf(fp, "# rot_pivot%d          %12.5e %12.5e %12.5e\n", g, rotg->pivot[XX], rotg->pivot[YY], rotg->pivot[ZZ]);
                break;
            case erotgFLEX1:
            case erotgFLEX2:
                fprintf(fp, "# rot_slab_distance%d   %f nm\n", g, rotg->slab_dist);
                fprintf(fp, "# rot_min_gaussian%d   %12.5e\n", g, rotg->min_gaussian);
                fprintf(fp, "# COG of ref. grp. %d  %12.5e %12.5e %12.5e  (only needed for fit)\n", 
                        g, erg->xc_ref_center[XX], erg->xc_ref_center[YY], erg->xc_ref_center[ZZ]);
                break;
            default:
                break;
            }
        }
        
        fprintf(fp, "#\n# Legend for the following data columns:\n");
        fprintf(fp, "#     ");
        print_aligned_short(fp, "t");
        nsets = 0;
        snew(setname, 4*rot->ngrp);
        
        for (g=0; g<rot->ngrp; g++)
        {
            rotg = &rot->grp[g];
            sprintf(buf, "theta_ref%d (degree)", g);
            print_aligned_group(fp, "theta_ref", g);
            setname[nsets] = strdup(buf);
            nsets++;
        }
        for (g=0; g<rot->ngrp; g++)
        {
            rotg = &rot->grp[g];
            if (rotg->eType==erotgFIXED || rotg->eType==erotgFIXED_PLANE)
            {
                sprintf(buf, "theta-av%d (degree)", g);
                print_aligned_group(fp, "theta_av", g);
                setname[nsets] = strdup(buf);
                nsets++;
                sprintf(buf, "tau%d (kJ/mol)", g);
                print_aligned_group(fp, "tau", g);
                setname[nsets] = strdup(buf);
                nsets++;
            }
            sprintf(buf, "energy%d (kJ/mol)", g);
            print_aligned_group(fp, "energy", g);
            setname[nsets] = strdup(buf);
            nsets++;
        }
        fprintf(fp, "\n#\n");
        
        if (nsets > 1)
            xvgr_legend(fp, nsets, setname, oenv);
        for(g=0; g<nsets; g++)
            sfree(setname[g]);
        sfree(setname);
        
        fflush(fp);
    }
    
    return fp;
}


/* Call on master only */
static FILE *open_angles_out(t_rot *rot, const char *fn)
{
    int      g;
    FILE     *fp=NULL;
    t_rotgrp *rotg;


    /* Open output file and write some information about it's structure: */
    fp = open_output_file(fn, rot->nstrout);
    fprintf(fp, "# All angles given in degrees, time in ps\n");
    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (rotg->eType == erotgFLEX1 || rotg->eType == erotgFLEX2)
            fprintf(fp, "# Rotation group %d (%s), slab distance %f nm, fit type %s\n", 
                    g, erotg_names[rotg->eType], rotg->slab_dist, erotg_fitnames[rotg->eFittype]);
    }
    fprintf(fp, "# The following columns will have the syntax:\n");
    fprintf(fp, "#     ");
    print_aligned_short(fp, "t");
    print_aligned_short(fp, "grp");
    print_aligned(fp, "theta_ref");
    print_aligned(fp, "theta_fit");
    print_aligned_short(fp, "slab");
    print_aligned_short(fp, "atoms");
    print_aligned(fp, "theta_fit");
    print_aligned_short(fp, "slab");
    print_aligned_short(fp, "atoms");
    print_aligned(fp, "theta_fit");
    fprintf(fp, " ...\n");
    fflush(fp);
    return fp;
}


/* Open torque output file and write some information about it's structure.
 * Call on master only */
static FILE *open_torque_out(t_rot *rot, const char *fn)
{
    FILE      *fp;
    int       g;
    t_rotgrp  *rotg;


    fp = open_output_file(fn, rot->nsttout);

    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (rotg->eType == erotgFLEX1 || rotg->eType == erotgFLEX2)
        {
            fprintf(fp, "# Rotation group %d (%s), slab distance %f nm\n", g, erotg_names[rotg->eType], rotg->slab_dist);
            fprintf(fp, "# The scalar tau is the torque [kJ/mol] in the direction of the rotation vector.\n");
            fprintf(fp, "# To obtain the vectorial torque, multiply tau with\n");
            fprintf(fp, "# rot_vec%d            %10.3e %10.3e %10.3e\n", g, rotg->vec[XX], rotg->vec[YY], rotg->vec[ZZ]);
            fprintf(fp, "#\n");
        }
    }
    fprintf(fp, "# The following columns will have the syntax (tau=torque for that slab):\n");
    fprintf(fp, "#     ");
    print_aligned_short(fp, "t");
    print_aligned_short(fp, "grp");
    print_aligned_short(fp, "slab");
    print_aligned(fp, "tau");
    print_aligned_short(fp, "slab");
    print_aligned(fp, "tau");
    fprintf(fp, " ...\n");
    fflush(fp);

    return fp;
}


static void swap_val(double* vec, int i, int j)
{
    double tmp = vec[j];
    
    
    vec[j]=vec[i];
    vec[i]=tmp;
}


static void swap_col(double **mat, int i, int j)
{
    double tmp[3] = {mat[0][j], mat[1][j], mat[2][j]};
    
    
    mat[0][j]=mat[0][i];
    mat[1][j]=mat[1][i];
    mat[2][j]=mat[2][i];
    
    mat[0][i]=tmp[0];
    mat[1][i]=tmp[1];
    mat[2][i]=tmp[2];
} 


/* Eigenvectors are stored in columns of eigen_vec */
static void diagonalize_symmetric(
        double **matrix,
        double **eigen_vec,
        double eigenval[3])
{
    int n_rot;
    
    
    jacobi(matrix,3,eigenval,eigen_vec,&n_rot);
    
    /* sort in ascending order */
    if (eigenval[0] > eigenval[1])
    {
        swap_val(eigenval, 0, 1);
        swap_col(eigen_vec, 0, 1);
    } 
    if (eigenval[1] > eigenval[2])
    {
        swap_val(eigenval, 1, 2);
        swap_col(eigen_vec, 1, 2);
    }
    if (eigenval[0] > eigenval[1])
    {
        swap_val(eigenval, 0, 1);
        swap_col(eigen_vec, 0, 1);
    }
}


static void align_with_z(
        rvec* s,           /* Structure to align */
        int natoms,
        rvec axis)
{
    int    i, j, k;
    rvec   zet = {0.0, 0.0, 1.0};
    rvec   rot_axis={0.0, 0.0, 0.0};
    rvec   *rotated_str=NULL;
    real   ooanorm;
    real   angle;
    matrix rotmat;
    
    
    snew(rotated_str, natoms);

    /* Normalize the axis */
    ooanorm = 1.0/norm(axis);
    svmul(ooanorm, axis, axis);
    
    /* Calculate the angle for the fitting procedure */
    cprod(axis, zet, rot_axis);
    angle = acos(axis[2]);
    if (angle < 0.0)
        angle += M_PI;
    
    /* Calculate the rotation matrix */
    calc_rotmat(rot_axis, angle*180.0/M_PI, rotmat);
    
    /* Apply the rotation matrix to s */
    for (i=0; i<natoms; i++)
    {    
        for(j=0; j<3; j++)
        {
            for(k=0; k<3; k++)
            {
                rotated_str[i][j] += rotmat[j][k]*s[i][k];
            }
        }
    }
    
    /* Rewrite the rotated structure to s */
    for(i=0; i<natoms; i++)
    {
        for(j=0; j<3; j++)
        {
            s[i][j]=rotated_str[i][j];
        }
    }
    
    sfree(rotated_str);
} 


static void calc_correl_matrix(rvec* Xstr, rvec* Ystr, double** Rmat, int natoms)
{    
    int i, j, k;
 
    
    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            Rmat[i][j] = 0.0;
    
    for (i=0; i<3; i++) 
        for (j=0; j<3; j++) 
            for (k=0; k<natoms; k++) 
                Rmat[i][j] += Ystr[k][i] * Xstr[k][j];
}


static void weigh_coords(rvec* str, real* weight, int natoms)
{
    int i, j;
    
    
    for(i=0; i<natoms; i++)
    {
        for(j=0; j<3; j++)
            str[i][j] *= sqrt(weight[i]);
    }  
}


static double opt_angle_analytic(
        rvec* ref_s,
        rvec* act_s,
        real* weight, 
        int natoms,
        rvec ref_com,
        rvec act_com,
        rvec axis)
{    
    int    i, j, k;
    rvec   *ref_s_1=NULL;
    rvec   *act_s_1=NULL;
    rvec   shift;
    double **Rmat, **RtR, **eigvec;
    double eigval[3];
    double V[3][3], WS[3][3];
    double rot_matrix[3][3];
    double opt_angle;
    
    
    /* Do not change the original coordinates */ 
    snew(ref_s_1, natoms);
    snew(act_s_1, natoms);
    for(i=0; i<natoms; i++)
    {
        copy_rvec(ref_s[i], ref_s_1[i]);
        copy_rvec(act_s[i], act_s_1[i]);
    }
    
    /* Translate the structures to the origin */
    shift[XX] = -ref_com[XX];
    shift[YY] = -ref_com[YY];
    shift[ZZ] = -ref_com[ZZ];
    translate_x(ref_s_1, natoms, shift);
    
    shift[XX] = -act_com[XX];
    shift[YY] = -act_com[YY];
    shift[ZZ] = -act_com[ZZ];
    translate_x(act_s_1, natoms, shift);
    
    /* Align rotation axis with z */
    align_with_z(ref_s_1, natoms, axis);
    align_with_z(act_s_1, natoms, axis);
    
    /* Correlation matrix */
    Rmat = allocate_square_matrix(3);
    
    for (i=0; i<natoms; i++)
    {
        ref_s_1[i][2]=0.0;
        act_s_1[i][2]=0.0;
    }
    
    /* Weight positions with sqrt(weight) */
    if (weight)
    {
        weigh_coords(ref_s_1, weight, natoms);
        weigh_coords(act_s_1, weight, natoms);
    }
    
    /* Calculate correlation matrices R=YXt (X=ref_s; Y=act_s) */
    calc_correl_matrix(ref_s_1, act_s_1, Rmat, natoms);
    
    /* Calculate RtR */
    RtR = allocate_square_matrix(3);
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            for (k=0; k<3; k++)
            {
                RtR[i][j] += Rmat[k][i] * Rmat[k][j];
            }
        }
    }
    /* Diagonalize RtR */
    snew(eigvec,3);
    for (i=0; i<3; i++)
        snew(eigvec[i],3);
    
    diagonalize_symmetric(RtR, eigvec, eigval);
    swap_col(eigvec,0,1);
    swap_col(eigvec,1,2);
    swap_val(eigval,0,1);
    swap_val(eigval,1,2);
    
    /* Calculate V */
    for(i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        {
            V[i][j]  = 0.0;
            WS[i][j] = 0.0;
        }
    }
    
    for (i=0; i<2; i++)
        for (j=0; j<2; j++)
            WS[i][j] = eigvec[i][j] / sqrt(eigval[j]);
    
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            for (k=0; k<3; k++)
            {
                V[i][j] += Rmat[i][k]*WS[k][j];
            }
        }
    }
    free_square_matrix(Rmat, 3);
    
    /* Calculate optimal rotation matrix */
    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            rot_matrix[i][j] = 0.0;
    
    for (i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        {
            for(k=0; k<3; k++){
                rot_matrix[i][j] += eigvec[i][k]*V[j][k];
            }
        }
    }
    rot_matrix[2][2] = 1.0;
        
    /* Determine the optimal rotation angle: */
    opt_angle = (-1.0)*acos(rot_matrix[0][0])*180.0/M_PI;
    if (rot_matrix[0][1] < 0.0)
        opt_angle = (-1.0)*opt_angle;
        
    /* Give back some memory */
    free_square_matrix(RtR, 3);
    sfree(ref_s_1);
    sfree(act_s_1);
    for (i=0; i<3; i++)
        sfree(eigvec[i]);
    sfree(eigvec);
    
    return opt_angle;
}


/* Determine actual angle of this slab by RMSD fit */
/* Not parallelized, call this routine only on the master */
static void flex_fit_angle(
        int  g,
        t_rotgrp *rotg,
        double t,
        real degangle,
        FILE *fp)
{
    int         i,l,n,islab,ind;
    rvec        curr_x, ref_x;
    rvec        *fitcoords=NULL;
    rvec        act_center;  /* Center of actual positions that are passed to the fit routine */
    rvec        ref_center;  /* Same for the reference positions */
    double      fitangle;    /* This will be the actual angle of the rotation group derived
                              * from an RMSD fit to the reference structure at t=0 */
    t_gmx_slabdata *sd;
    rvec        coord;
    real        scal;
    gmx_enfrotgrp_t erg;     /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;

    /**********************************/
    /* First collect the data we need */
    /**********************************/

    /* Loop over slabs */
    for (n = erg->slab_first; n <= erg->slab_last; n++)
    {
        islab = n - erg->slab_first; /* slab index */
        sd = &(rotg->enfrotgrp->slab_data[islab]);
        sd->nat = erg->lastatom[n]-erg->firstatom[n]+1;
        ind = 0;
        
        /* Loop over the relevant atoms in the slab */
        for (l=erg->firstatom[n]; l<=erg->lastatom[n]; l++)
        {            
            /* Current position of this atom: x[ii][XX/YY/ZZ] */
            copy_rvec(erg->xc[l], curr_x);
            
            /* The (unrotated) reference position of this atom is copied to ref_x.
             * Beware, the xc coords have been sorted in do_flex! */
            copy_rvec(erg->xc_ref_sorted[l], ref_x);
            
            /* Save data for doing angular RMSD fit later */
            /* Save the current atom position */
            copy_rvec(curr_x, sd->x[ind]);
            /* Save the corresponding reference position */
            copy_rvec(ref_x , sd->ref[ind]);
            /* Save the weight for this atom in this slab */
            sd->weight[ind] = gaussian_weight(curr_x, rotg, n);
            
            /* Next atom in this slab */
            ind++;
        }
    }

    /* Get the geometrical center of the whole rotation group: */
    get_center(erg->xc, NULL, rotg->nat, act_center);


    /******************************/
    /* Now do the fit calculation */
    /******************************/

    /* === Determine the optimal fit angle for the whole rotation group === */
    if (rotg->eFittype == erotgFitNORM)
    {
        /* Normalize every position to it's reference length 
         * prior to performing the fit */
        for (i=0; i<rotg->nat; i++)
        {
            /* First put geometrical center of positions into origin */
            rvec_sub(erg->xc[i], act_center, coord);
            /* Determine the scaling factor for the length: */
            scal = erg->xc_ref_length[erg->xc_sortind[i]] / norm(coord);
            /* Get position, multiply with the scaling factor and save in buf[i] */
            svmul(scal, coord, erg->xc_norm[i]);
        }
        fitcoords = erg->xc_norm;
    }
    else
    {
        fitcoords = erg->xc;
    }
    /* Note that from the point of view of the current positions, the reference has rotated backwards,
     * but we want to output the angle relative to the fixed reference, therefore the minus sign. */
    fitangle = -opt_angle_analytic(erg->xc_ref_sorted, fitcoords, NULL, rotg->nat, erg->xc_ref_center, act_center, rotg->vec);    
    fprintf(fp, "%12.3e%6d%12.3f%12.3lf", t, g, degangle, fitangle);


    /* === Now do RMSD fitting for each slab === */
    /* We require at least SLAB_MIN_ATOMS in a slab, such that the fit makes sense. */
#define SLAB_MIN_ATOMS 4

    for (n = erg->slab_first; n <= erg->slab_last; n++)
    {
        islab = n - erg->slab_first; /* slab index */
        sd = &(rotg->enfrotgrp->slab_data[islab]);
        if (sd->nat >= SLAB_MIN_ATOMS)
        {
            /* Get the center of the slabs reference and current positions */
            get_center(sd->ref, sd->weight, sd->nat, ref_center);
            get_center(sd->x  , sd->weight, sd->nat, act_center);
            if (rotg->eFittype == erotgFitNORM)
            {
                /* Normalize every position to it's reference length 
                 * prior to performing the fit */
                for (i=0; i<sd->nat;i++) /* Center */
                {
                    rvec_dec(sd->ref[i], ref_center);
                    rvec_dec(sd->x[i]  , act_center);
                    /* Normalize x_i such that it gets the same length as ref_i */
                    svmul( norm(sd->ref[i])/norm(sd->x[i]), sd->x[i], sd->x[i] );
                }
                /* We already subtracted the centers */
                clear_rvec(ref_center);
                clear_rvec(act_center);
            }
            fitangle = -opt_angle_analytic(sd->ref, sd->x, sd->weight, sd->nat, ref_center, act_center, rotg->vec);
            fprintf(fp, "%6d%6d%12.3f", n, sd->nat, fitangle);
        }
    }
    fprintf(fp     , "\n");

#undef SLAB_MIN_ATOMS
}


/* Shift x with is */
static inline void shift_single_coord(matrix box, rvec x, const ivec is)
{
    int tx,ty,tz;


    tx=is[XX];
    ty=is[YY];
    tz=is[ZZ];

    if(TRICLINIC(box))
    {
        x[XX] += tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
        x[YY] += ty*box[YY][YY]+tz*box[ZZ][YY];
        x[ZZ] += tz*box[ZZ][ZZ];
    } else
    {
        x[XX] += tx*box[XX][XX];
        x[YY] += ty*box[YY][YY];
        x[ZZ] += tz*box[ZZ][ZZ];
    }
}


/* Determine the 'home' slab of this atom which is the
 * slab with the highest Gaussian weight of all */
#define round(a) (int)(a+0.5)
static inline int get_homeslab(
        rvec curr_x,   /* The position for which the home slab shall be determined */ 
        rvec rotvec,   /* The rotation vector */
        real slabdist) /* The slab distance */
{
    real dist;
    
    
    /* The distance of the atom to the coordinate center (where the
     * slab with index 0) is */
    dist = iprod(rotvec, curr_x);
    
    return round(dist / slabdist); 
}


/* For a local atom determine the relevant slabs, i.e. slabs in
 * which the gaussian is larger than min_gaussian
 */
static int get_single_atom_gaussians(
        rvec      curr_x,
        t_commrec *cr,
        t_rotgrp  *rotg)
{
   int slab, homeslab;
   real g;
   int count = 0;
   gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */

   
   erg=rotg->enfrotgrp;
   
   /* Determine the 'home' slab of this atom: */
   homeslab = get_homeslab(curr_x, rotg->vec, rotg->slab_dist);

   /* First determine the weight in the atoms home slab: */
   g = gaussian_weight(curr_x, rotg, homeslab);
   
   erg->gn_atom[count] = g;
   erg->gn_slabind[count] = homeslab;
   count++;
   
   
   /* Determine the max slab */
   slab = homeslab;
   while (g > rotg->min_gaussian)
   {
       slab++;
       g = gaussian_weight(curr_x, rotg, slab);
       erg->gn_slabind[count]=slab;
       erg->gn_atom[count]=g;
       count++;
   }
   count--;
   
   /* Determine the max slab */
   slab = homeslab;
   do
   {
       slab--;
       g = gaussian_weight(curr_x, rotg, slab);       
       erg->gn_slabind[count]=slab;
       erg->gn_atom[count]=g;
       count++;
   }
   while (g > rotg->min_gaussian);
   count--;
   
   return count;
}


static real do_flex2_lowlevel(
        t_rotgrp  *rotg,
        real      sigma,    /* The Gaussian width sigma */
        rvec      x[],
        bool      bCalcTorque,
        matrix    box,
        t_commrec *cr)
{
    int  count,i,ic,ii,l,m,n,islab,iigrp;
    rvec dr;                /* difference vector between actual and reference position */
    real V;                 /* The rotation potential energy */
    real gaussian;          /* Gaussian weight */
    real gaussian_xi;
    real beta;
    rvec curr_x;            /* particle position */
    rvec xi;
    rvec curr_x_rel;        /* particle position relative to COG */
    rvec curr_COG;          /* the current slab's center of geometry (COG) */
    rvec ref_x, ref_x_cpy;  /* the reference position */
    rvec ref_COG;           /* the reference slab's COG */
    rvec r, yi, tmp;
    rvec force_n;           /* Single force from slab n on one atom */
    rvec s_ii;
    real inv_norm_ii;
    real sdotr_ii;          /* s_ii.r_ii */
    real tmp_s;
    rvec tmp_v, tmp_n1_v, tmp_n2_v, tmp_n3_v, tmp_n4_v;
    rvec sum_i1;
    real sum_i2;
    rvec s_i;
    real inv_norm_i;
    real sdotr_i;           /* s_i.r_i */
    real gauss_ratio;       /* g_n(x_ii)/Sum_i(g_n(x_i)*/
    real beta_sigma;        /* beta/sigma^2 */
    rvec sum_f_ii;          /* force on one atom summed over all slabs */
    gmx_enfrotgrp_t erg;    /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;
    
    /********************************************************/
    /* Main loop over all local atoms of the rotation group */
    /********************************************************/
    V = 0.0;
    for (l=0; l<erg->nat_loc; l++)
    {
        /* Local index of a rotation group atom  */
        ii = erg->ind_loc[l];
        /* Index of this rotation group atom with respect to the whole rotation group */
        iigrp = erg->xc_ref_ind[l];
        
        /* Current position of this atom: x[ii][XX/YY/ZZ] */
        rvec_sub(x[ii], erg->center, curr_x);

        /* Shift this atom such that it is near its reference */
        shift_single_coord(box, curr_x, erg->xc_shifts[iigrp]);

        /* Determine the slabs to loop over, i.e. the ones with contributions
         * larger than min_gaussian */
        count = get_single_atom_gaussians(curr_x, cr, rotg);
        
        clear_rvec(sum_f_ii);

        /* Loop over the relevant slabs for this atom */
        for (ic=0; ic < count; ic++)  
        {
            n = erg->gn_slabind[ic];
            
            /* Get the precomputed Gaussian value of curr_slab for curr_x */
            gaussian = erg->gn_atom[ic];

            islab = n - erg->slab_first; /* slab index */
            
            /* The (unrotated) reference position of this atom is copied to ref_x: */
            copy_rvec(erg->xc_ref[iigrp], ref_x);
            beta = calc_beta(curr_x, rotg,n);
            /* The center of geometry (COG) of this slab is copied to curr_COG: */
            copy_rvec(erg->slab_center[islab], curr_COG);
            /* The reference COG of this slab is copied to ref_COG: */
            copy_rvec(erg->slab_center_ref[islab+erg->slab_buffer], ref_COG);
            
            /* Rotate the reference position around the rotation axis through this slab's reference COG */
            /* 1. Subtract the reference slab COG from the reference position i */
            rvec_sub(ref_x, ref_COG, ref_x); /* ref_x =y_ii-y_0 */
            /* 2. Rotate reference position around origin: */
            copy_rvec(ref_x, ref_x_cpy);
            mvmul(erg->rotmat, ref_x_cpy, ref_x); /* ref_x = r_ii = Omega.(y_ii-y_0) */
            
            /* Now subtract the slab COG from current position i */
            rvec_sub(curr_x, curr_COG, curr_x_rel); /* curr_x_rel = x_ii-x_0 */
            
            /* Force on atom i is based on difference vector between current position 
             * and rotated reference position */
            /* Difference vector between current and reference position: */
            rvec_sub(curr_x_rel, ref_x, dr); /* dr=(x_ii-x_0)-Omega.(y_ii-y_0) */
            
            cprod(rotg->vec, curr_x_rel, s_ii); /* s_ii = v x (x_ii-x_0) */
            inv_norm_ii=1.0/norm(s_ii);         /* inv_norm_ii = 1/|v x (x_ii-x_0)| */
            unitv(s_ii, s_ii);                  /* s_ii = v x (x_ii-x_0)/|v x (x_ii-x_0)| */
            sdotr_ii=iprod(s_ii,ref_x);         /* sdotr_ii = ((v x (x_ii-x_0)/|v x (x_ii-x_0)|).Omega.(y_ii-y_0) */
            
            /*********************************/
            /* Add to the rotation potential */
            /*********************************/
            V += 0.5*rotg->k*gaussian*sqr(sdotr_ii);
            
            tmp_s=gaussian*sdotr_ii*inv_norm_ii;
            svmul(sdotr_ii, s_ii, tmp_n1_v);
            rvec_sub(ref_x, tmp_n1_v, tmp_n1_v);
            svmul(tmp_s, tmp_n1_v, tmp_n1_v);
            
            clear_rvec(sum_i1);
            sum_i2=0.0;
            
            for (i=erg->firstatom[n]; i<=erg->lastatom[n]; i++)
            {
                /* Coordinate xi of this atom */
                copy_rvec(erg->xc[i],xi);
                
                gaussian_xi = gaussian_weight(xi,rotg,n);

                /* Calculate r_i for atom i and slab n: */
                /* Unrotated reference y_i: */
                copy_rvec(erg->xc_ref_sorted[i],yi);
                
                /* COG y0 for this slab: */
                /* The reference COG of this slab is still in ref_COG */
                rvec_sub(yi, ref_COG, tmp);   /* tmp = y_i - y_0                       */
                mvmul(erg->rotmat, tmp, r);   /* r   = Omega*(y_i - y_0)               */
                rvec_sub(xi, curr_COG, tmp);  /* tmp = (x_i - x_0)                     */
                cprod(rotg->vec, tmp, s_i);   /* s_i = v x (x_i - x_0)                 */
                inv_norm_i=1.0/norm(s_i);     /* 1/|v x (x_i - x_0)|                   */
                unitv(s_i, s_i);              /* s_i = (v x (x_i-x_0))/|v x (x_i-x_0)| */
                sdotr_i=iprod(s_i,r);         /* sdotr_i = (s_i.r_i)                   */
                
                tmp_s=gaussian_xi*sdotr_i*inv_norm_i;     /* tmp_s = g_n(xi)*(s_i.r_i)*1/norm */
                /* sum_11 */
                svmul(sdotr_i, s_i, tmp_v);
                rvec_sub(tmp_v, r, tmp_v);
                svmul(tmp_s, tmp_v, tmp_v); /* tmp_v = g_n(xi)*(s_i.r_i)*1/norm *(-r_i+(r_i.s_i)s_i) */
                rvec_add(sum_i1, tmp_v, sum_i1); /* n2 */
                sum_i2 += tmp_s*(iprod(r,s_ii)-(iprod(s_i,s_ii)*sdotr_i)); /* n3 */
            } /* now we have the sum-over-atoms (over i) for the ii-th atom in the n-th slab */
            
            /* We can safely divide by slab_weights since we check in get_slab_centers
             * that it is non-zero. */
            gauss_ratio=gaussian/erg->slab_weights[islab];
            
            beta_sigma=beta/sqr(sigma);
            svmul(gauss_ratio, sum_i1, tmp_n2_v);   /* tmp_n2_v = gauss_ratio_ii*Sum_i(g_n(xi)*(s_i.r_i)*1/norm *(-r_i+(r_i.s_i)s_i)) */
            
            rvec_add(tmp_n1_v, tmp_n2_v, force_n);  /* adding up the terms perpendicular to v and to x_ii-x0 */
            cprod(force_n, rotg->vec, tmp_v);       /* now these are indeed perpendicular */
            svmul((-1.0*rotg->k), tmp_v, force_n);  /* multiplying by -k* we've got the final tangent contribution */
            /* calculating the parallel contribution */
            svmul((-1.0)*rotg->k*gauss_ratio*beta_sigma*sum_i2,rotg->vec,tmp_n3_v);
            svmul(0.5*rotg->k*beta_sigma*gaussian*sqr(sdotr_ii),rotg->vec,tmp_n4_v);
            rvec_add(tmp_n3_v, tmp_n4_v, tmp_n3_v); /* tmp_n3_v contains the final parallel contribution */
            rvec_add(tmp_n3_v, force_n, force_n);   /* force_n is the total force from slab n */
            /* sum the forces over slabs */
            rvec_add(sum_f_ii,force_n,sum_f_ii);
            
            /* Calculate the torque: */
            if (bCalcTorque)
            {
                /* The force on atom ii from slab n only: */
                erg->slab_torque_v[islab] += torque(rotg->vec, force_n, curr_x, curr_COG);
            }
        } /* END of loop over slabs */

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        for (m=0; m<DIM; m++)
            erg->f_rot_loc[l][m] = sum_f_ii[m];

#ifdef INFOF
        fprintf(stderr," FORCE on ATOM %d/%d  = (%15.8f %15.8f %15.8f)  \n",
                l,ii,rotg->sum_f_ii[XX], rotg->sum_f_ii[YY], rotg->sum_f_ii[ZZ]5);
#endif
    } /* END of loop over local atoms */

    return V;
}


static real do_flex_lowlevel(
        t_rotgrp *rotg,
        real      sigma,     /* The Gaussian width sigma */
        rvec      x[],
        bool      bCalcTorque,
        matrix    box,
        t_commrec *cr)
{
    int   count,i,ic,ii,iii,l,m,n,islab,iigrp;
    rvec  direction;         /* the direction for the force on atom i */
    rvec  sum_n1,sum_n2;     /* Two contributions to the rotation force */
    rvec  sum_i;             /* Inner part of sum_n2 */
    rvec  dummy1,tmp_f,s,tmp2;
    real  u;                 /* direction*dr */
    real  V;                 /* The rotation potential energy */
    real  gaussian;          /* Gaussian weight */
    real  gaussian_xi;
    real  fac,beta;
    rvec  curr_x;            /* position */
    rvec  xi;
    rvec  curr_x_rel;        /* position relative to COG */
    rvec  curr_COG;          /* the current slab's center of geometry (COG) */
    rvec  ref_x, ref_x_cpy;  /* the reference position */
    rvec  ref_COG;           /* the reference slab's COG */
    rvec  r, yi, tmp;
    rvec  force_n;           /* Single force from slab n on one atom */
    gmx_enfrotgrp_t erg;     /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;

    /********************************************************/
    /* Main loop over all local atoms of the rotation group */
    /********************************************************/
    V = 0.0;
    for (l=0; l<erg->nat_loc; l++)
    {
        /* Local index of a rotation group atom  */
        ii = erg->ind_loc[l];
        /* Index of this rotation group atom with respect to the whole rotation group */
        iigrp = erg->xc_ref_ind[l];
        
        /* Current position of this atom: x[ii][XX/YY/ZZ] */
        rvec_sub(x[ii], erg->center, curr_x);
        
        /* Shift this atom such that it is near its reference */
        shift_single_coord(box, curr_x, erg->xc_shifts[iigrp]);

        /* Determine the slabs to loop over, i.e. the ones with contributions
         * larger than min_gaussian */
        count = get_single_atom_gaussians(curr_x, cr, rotg);

        clear_rvec(sum_n1);
        clear_rvec(sum_n2);

        /* Loop over the relevant slabs for this atom */
        for (ic=0; ic < count; ic++)  
        {
            n = erg->gn_slabind[ic];
                
            /* Get the precomputed Gaussian value of curr_slab for curr_x */
            gaussian = erg->gn_atom[ic];

            islab = n - erg->slab_first; /* slab index */
            
            /* The (unrotated) reference position of this atom is copied to ref_x: */
            copy_rvec(erg->xc_ref[iigrp], ref_x);
            beta = calc_beta(curr_x, rotg,n);
            /* The center of geometry (COG) of this slab is copied to curr_COG: */
            copy_rvec(erg->slab_center[islab], curr_COG);
            /* The reference COG of this slab is copied to ref_COG: */
            copy_rvec(erg->slab_center_ref[islab+erg->slab_buffer], ref_COG);
            
            /* Rotate the reference position around the rotation axis through this slab's reference COG */
            /* 1. Subtract the reference slab COG from the reference position i */
            rvec_sub(ref_x, ref_COG, ref_x); /* ref_x =y_ii-y_0 */
            /* 2. Rotate reference position around origin: */
            copy_rvec(ref_x, ref_x_cpy);
            mvmul(erg->rotmat, ref_x_cpy, ref_x); /* ref_x = r_ii = Omega.(y_ii-y_0) */
            
            /* Now subtract the slab COG from current position i */
            rvec_sub(curr_x, curr_COG, curr_x_rel); /* curr_x_rel = x_ii-x_0 */
            
            /* Calculate the direction of the actual force */
            cprod(rotg->vec, ref_x, direction);
            unitv(direction,direction);
            u = iprod(direction, curr_x_rel);
            
            /*********************************/
            /* Add to the rotation potential */
            /*********************************/
            V += 0.5*rotg->k*gaussian*sqr(u);
            
            /*************************************************/
            /* sum_n1 is the main contribution to the force: */
            /*************************************************/
            /* gn*u*(vec_s - 0.5*u*beta/(sigma^2)*vec_a) */
            svmul(0.5*u*beta/sqr(sigma),rotg->vec,dummy1);
            /* Typically 'direction' will be the largest part */
            rvec_sub(direction,dummy1,dummy1);
            svmul(gaussian*u,dummy1,dummy1);
            /* Sum over n: */
            rvec_add(sum_n1,dummy1,sum_n1);
            
            /*************************************************************/
            /* For the term sum_n2 we need to loop over all atoms again: */
            /*************************************************************/
            clear_rvec(sum_i);
            
            GMX_MPE_LOG(ev_inner_loop_start);
            
            for (i=erg->firstatom[n]; i<=erg->lastatom[n]; i++)
            {
                /* Index of a rotation group atom  */
                iii = rotg->ind[i];
                
                /* Coordinate xi of this atom */
                copy_rvec(erg->xc[i],xi);
                
                gaussian_xi = gaussian_weight(xi,rotg,n);

                /* Calculate r=Omega*(y_i-y_0) for atom i and slab n: */
                /* Unrotated reference position y_i: */
                copy_rvec(erg->xc_ref[i],yi);
                
                /* COG y0 for this slab: */
                /* The reference COG of this slab is still in ref_COG */
                rvec_sub(yi, ref_COG, tmp);   /* tmp = y_i - y_0             */
                mvmul(erg->rotmat, tmp, r);   /* r   = Omega*(y_i - y_0)     */
                cprod(rotg->vec, r, tmp);     /* tmp = v x Omega*(y_i - y_0) */
                unitv(tmp, s);                /* s   = v x Omega*(y_i - y_0) / |s x Omega*(y_i - y_0)| */
                /* Now that we have si, let's calculate the i-sum: */
                /* tmp = x_0 - x_l */
                rvec_sub(curr_COG, curr_x, tmp);
                /* tmp2 = beta/sigma^2 * (s*(x_0 - x_l)) * v   */
                svmul(beta*iprod(tmp, s)/sqr(sigma), rotg->vec, tmp2);
                /* tmp = s + tmp2 */
                rvec_add(tmp2, s, tmp);
                /* tmp2 = xi - x0 */
                rvec_sub(xi, curr_COG, tmp2);
                /* fac = gn * s*(xi - x0 - ri) */
                fac = gaussian_xi*iprod(s, tmp2);
                /* tmp2 = gn * s*(xi - x0) * [beta/sigma^2 * (s*(x_0 - x_l)) * v] */
                svmul(fac, tmp, tmp2);
                rvec_add(sum_i, tmp2, sum_i);
                
            } /* now we have the i-sum for this atom in this slab */
            
            /* We can safely divide by slab_weights since we check in get_slab_centers
             * that it is non-zero. */
            svmul(gaussian/erg->slab_weights[islab], sum_i, sum_i);

            rvec_add(sum_n2, sum_i, sum_n2);
            
            GMX_MPE_LOG(ev_inner_loop_finish);

            /* Calculate the torque: */
            if (bCalcTorque)
            {
                /* The force on atom ii from slab n only: */
                rvec_sub(sum_i, dummy1, force_n);
                svmul(rotg->k, force_n, force_n);
                erg->slab_torque_v[islab] += torque(rotg->vec, force_n, curr_x, curr_COG);
            }
        } /* END of loop over slabs */

        /* Put both contributions together: */
        rvec_sub(sum_n2,sum_n1,tmp_f); /* F = -grad V */

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        for(m=0; m<DIM; m++)
            erg->f_rot_loc[l][m] = rotg->k*tmp_f[m];

#ifdef INFOF
        static bool bFirst=1;
        char buf[255];
        static FILE *fp;
        if (bFirst)
        {
            sprintf(buf, "forces%d.txt", cr->nodeid);
            fp = fopen(buf, "w");
            bFirst = 0;
        }

        fprintf(fp," FORCE on ATOM %d/%d  = (%15.8f %15.8f %15.8f)  1: %15.8f %15.8f %15.8f   2: %15.8f %15.8f %15.8f\n",
                l,ii,rotg->k*tmp_f[XX], rotg->k*tmp_f[YY], rotg->k*tmp_f[ZZ],
                -rotg->k*sum_n1[XX], -rotg->k*sum_n1[YY], -rotg->k*sum_n1[ZZ],
                rotg->k*sum_n2[XX],  rotg->k*sum_n2[YY],  rotg->k*sum_n2[ZZ]);
#endif
    } /* END of loop over local atoms */

    return V;
}


#ifdef PRINT_COORDS
static void print_coordinates(t_commrec *cr, t_rotgrp *rotg, rvec x[], matrix box, int step)
{
    int i;
    static FILE *fp;
    static char buf[STRLEN];
    static bool bFirst=1;


    if (bFirst)
    {
        sprintf(buf, "coords%d.txt", cr->nodeid);
        fp = fopen(buf, "w");
        bFirst = 0;
    }

    fprintf(fp, "\nStep %d\n", step);
    fprintf(fp, "box: %f %f %f %f %f %f %f %f %f\n",
            box[XX][XX], box[XX][YY], box[XX][ZZ],
            box[YY][XX], box[YY][YY], box[YY][ZZ],
            box[ZZ][XX], box[ZZ][ZZ], box[ZZ][ZZ]);
    for (i=0; i<rotg->nat; i++)
    {
        fprintf(fp, "%4d  %f %f %f\n", i,
                erg->xc[i][XX], erg->xc[i][YY], erg->xc[i][ZZ]);
    }
    fflush(fp);

}
#endif


static int projection_compare(const void *a, const void *b)
{
    sort_along_vec_t *xca, *xcb;
    
    
    xca = (sort_along_vec_t *)a;
    xcb = (sort_along_vec_t *)b;
    
    if (xca->xcproj < xcb->xcproj)
        return -1;
    else if (xca->xcproj > xcb->xcproj)
        return 1;
    else
        return 0;
}


/* Sort the collective coordinates along the rotation vector */
static void sort_collective_coordinates(
        t_rotgrp *rotg,         /* Rotation group */
        sort_along_vec_t *data, /* Buffer for sorting the positions */
        rvec *buf)              /* Buffer for sorting the positions */
{
    int i;
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;
    
    for (i=0; i<rotg->nat; i++)
    {
        data[i].xcproj = iprod(erg->xc[i], rotg->vec);
        data[i].ind    = i;
    }
    qsort(data, rotg->nat, sizeof(sort_along_vec_t), projection_compare);
    
    for (i=0; i<rotg->nat; i++)
    {
        copy_rvec(erg->xc[data[i].ind], buf[i]);
        copy_rvec(erg->xc_ref[data[i].ind], erg->xc_ref_sorted[i]);
        erg->xc_sortind[i] = data[i].ind;
    }

    for (i=0; i<rotg->nat; i++)
    {
        copy_rvec(buf[i], erg->xc[i]);
    }
}


/* For each slab, get the first and the last index of the sorted atom
 * indices */
static void get_firstlast_atom_per_slab(t_rotgrp *rotg, t_commrec *cr)
{
    int i,n;
    real beta;
    gmx_enfrotgrp_t erg;     /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;

    GMX_MPE_LOG(ev_get_firstlast_start);
    
    /* Find the first atom that needs to enter the calculation for each slab */
    n = erg->slab_first;
    i = 0; /* start with the first atom */
    do
    {
        /* Find the first atom that significantly contributes to this slab */
        do /* move forward in position until a large enough beta is found */
        {
            beta = calc_beta(erg->xc[i], rotg, n);
            i++;
        } while ((beta < -erg->max_beta) && (i < rotg->nat));
        i--;
        erg->firstatom[n] = i;
        /* Proceed to the next slab */
        n++;
    } while (n <= erg->slab_last);
    
    /* Find the last atom for each slab */
     n = erg->slab_last; /* start with last slab */
     i = rotg->nat-1;  /* start with the last atom */
     do
     {
         do /* move backward in position until a large enough beta is found */
         {
             beta = calc_beta(erg->xc[i], rotg, n);
             i--;
         } while ((beta > erg->max_beta) && (i > -1));
         i++;
         erg->lastatom[n] = i;
         /* Proceed to the next slab */
         n--;
     } while (n >= erg->slab_first);
    
     GMX_MPE_LOG(ev_get_firstlast_finish);
}


/* Determine the very first and very last slab that needs to be considered 
 * For the first slab that needs to be considered, we have to find the smallest
 * n that obeys:
 * 
 * x_first * v - n*Delta_x <= beta_max
 * 
 * slab index n, slab distance Delta_x, rotation vector v. For the last slab we 
 * have to find the largest n that obeys
 * 
 * x_last * v - n*Delta_x >= -beta_max
 *  
 */
static inline int get_first_slab(
        t_rotgrp *rotg,     /* The rotation group (inputrec data) */
        real     max_beta,  /* The max_beta value, instead of min_gaussian */
        rvec     firstatom) /* First atom after sorting along the rotation vector v */
{
    /* Find the first slab for the first atom */   
    return ceil((iprod(firstatom, rotg->vec) - max_beta)/rotg->slab_dist);
}


static inline int get_last_slab(
        t_rotgrp *rotg,     /* The rotation group (inputrec data) */
        real     max_beta,  /* The max_beta value, instead of min_gaussian */
        rvec     lastatom)  /* Last atom along v */
{
    /* Find the last slab for the last atom */
    return floor((iprod(lastatom, rotg->vec) + max_beta)/rotg->slab_dist);    
}


static void get_firstlast_slab_check(
        t_rotgrp        *rotg,     /* The rotation group (inputrec data) */
        t_gmx_enfrotgrp *erg,      /* The rotation group (data only accessible in this file) */
        rvec            firstatom, /* First atom after sorting along the rotation vector v */
        rvec            lastatom,  /* Last atom along v */
        int             g,         /* The rotation group number */
        t_commrec       *cr)
{
    erg->slab_first = get_first_slab(rotg, erg->max_beta, firstatom);
    erg->slab_last  = get_last_slab(rotg, erg->max_beta, lastatom);

    /* Check whether we have reference data to compare against */
    if (erg->slab_first < erg->slab_first_ref)
        gmx_fatal(FARGS, "Enforced rotation: No reference data for first slab (n=%d), unable to proceed", 
                  erg->slab_first);
    
    /* Check whether we have reference data to compare against */
    if (erg->slab_last > erg->slab_last_ref)
        gmx_fatal(FARGS, "Enforced rotation: No reference data for last slab (n=%d), unable to proceed",
                  erg->slab_last);
}


/* Enforced rotation with a flexible axis */
static void do_flexible(
        t_commrec *cr,
        gmx_enfrot_t enfrot,  /* Other rotation data            */
        t_rotgrp  *rotg,      /* The rotation group             */
        int       g,          /* Group number                   */
        rvec      x[],        /* The local positions            */
        matrix    box,
        double    t,          /* Time in picoseconds            */
        int       step,       /* The time step                  */
        bool      bOutstep)
{
    int          l,nslabs;
    real         sigma;       /* The Gaussian width sigma */
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;

    /* Define the sigma value */
    sigma = 0.7*rotg->slab_dist;
    
    /* Sort the collective coordinates erg->xc along the rotation vector. This is
     * an optimization for the inner loop.
     */
    sort_collective_coordinates(rotg, enfrot->data, enfrot->buf);
    
    /* Determine the first relevant slab for the first atom and the last
     * relevant slab for the last atom */
    get_firstlast_slab_check(rotg, erg, erg->xc[0], erg->xc[rotg->nat-1], g, cr);
    
    /* Determine for each slab depending on the min_gaussian cutoff criterium,
     * a first and a last atom index inbetween stuff needs to be calculated */
    get_firstlast_atom_per_slab(rotg, cr);

    /* Determine the gaussian-weighted center of positions for all slabs */
    get_slab_centers(rotg,erg->xc,cr,g,t,enfrot->out_slabs,bOutstep,FALSE);
        
    /* Clear the torque per slab from last time step: */
    nslabs = erg->slab_last - erg->slab_first + 1;
    for (l=0; l<nslabs; l++)
        erg->slab_torque_v[l] = 0.0;
    
    /* Call the rotational forces kernel */
    GMX_MPE_LOG(ev_flexll_start);
    if (rotg->eType == erotgFLEX1)
        erg->V = do_flex_lowlevel(rotg, sigma, x, bOutstep, box, cr);
    else if (rotg->eType == erotgFLEX2)
        erg->V = do_flex2_lowlevel(rotg, sigma, x, bOutstep, box, cr);
    else
        gmx_fatal(FARGS, "Unknown flexible rotation type");
    GMX_MPE_LOG(ev_flexll_finish);
    
    /* Determine actual angle of this slab by RMSD fit and output to file - Let's hope */
    /* this only happens once in a while, since this is not parallelized! */
    if (bOutstep && MASTER(cr))
        flex_fit_angle(g, rotg, t, erg->degangle, enfrot->out_angles);
}


/* Calculate the angle between reference and actual rotation group atom: */
static int angle(t_rotgrp *rotg,
        rvec x_act,
        rvec x_ref,
        real *alpha,
        real *weight,  /* atoms near the rotation axis should count less than atoms far away */
        rvec offset)
{
    rvec x , xr ;  /* actual and reference positions in new coordinate system */
    rvec xp, xrp;  /* dito, but projected on a plane perpendicular to pg->vec */
    real normxp;   /* saved for efficiency reasons */
    rvec dum;
    real cosalpha; /* cos of angle between projected reference and actual position */
    int sign;
    real denominator;


    /* Move the center of positions to rot_offset: */
    rvec_sub(x_act, offset, x);
    rvec_sub(x_ref, offset, xr);

    /* Project xr and x into a plane through the origin perpendicular to rot_vec: */
    /* Project xr: xrp = xr - (vec * xr) * vec */
    svmul(iprod(rotg->vec, xr), rotg->vec, dum);
    rvec_sub(xr, dum, xrp);
    /* Project x: */
    svmul(iprod(rotg->vec, x), rotg->vec, dum);
    rvec_sub(x, dum, xp);

    /* Calculate the angle between the projected positions: */
    normxp = norm(xp); /* save for later use */
    denominator = norm(xrp)*normxp;
    
    /* Can we actually calculate this angle? */
    if (0.0 == denominator)
        return -1;

    cosalpha = iprod(xrp, xp) / denominator;
    if (cosalpha < -1.0) cosalpha = -1.0;
    if (cosalpha >  1.0) cosalpha =  1.0;

    /* Retrieve some information about which vector precedes */
    cprod(xp, xrp, dum); /* if reference precedes, this is pointing into the same direction as vec */

    if (iprod(rotg->vec, dum) >= 0)
        /* This will be the case when the reference group runs ahead. Thus the sign for the
         * angle of the actual group (which we are interested in) is negative = behind */
        sign = -1;
    else
        sign = 1;

    /* Return the angle in radians */
    *alpha = sign * acos(cosalpha);
    /* Also return the weight */
    *weight = normxp;
    
    /* Everything OK */
    return 0;
}


/* Project first vector onto a plane perpendicular to the second vector 
 * dr = dr - (dr.v)v
 */
static inline void project_onto_plane(rvec dr, const rvec v)
{
    rvec tmp;
    
    
    svmul(iprod(dr,v),v,tmp);  /* tmp = (dr.v)v */
    rvec_dec(dr, tmp);         /* dr = dr - (dr.v)v */
}


/* Fixed rotation: The rotation reference group rotates around an axis */
/* The atoms of the actual rotation group are attached with imaginary  */
/* springs to the reference atoms.                                     */
static void do_fixed(
        t_commrec *cr,
        t_rotgrp  *rotg,        /* The rotation group          */
        rvec      x[],          /* The positions               */
        t_pbc     *pbc,
        matrix    box,          /* The simulation box          */
        double    t,            /* Time in picoseconds         */
        int       step,         /* The time step               */
        bool      bTorque)
{
    int       i,m;
    rvec      dr;
    rvec      tmp_f;           /* Force */
    real      alpha;           /* a single angle between an actual and a reference position */
    real      weight;          /* single weight for a single angle */
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;

    /* Clear values from last time step */
    erg->V            = 0.0;
    erg->fix_torque_v = 0.0;
    erg->fix_angles_v = 0.0;
    erg->fix_weight_v = 0.0;
    
    /* Loop over all local atoms of the rotation group */
    for (i=0; i<erg->nat_loc; i++)
    {
        /* Subtract COG / COM */
        rvec_dec(erg->x_loc_pbc[i], erg->center);
        
        /* Difference vector between reference and actual position (both centered): */
        rvec_sub(erg->xr_loc[i], erg->x_loc_pbc[i], dr);
        
        if (rotg->eType==erotgFIXED_PLANE)
            project_onto_plane(dr, rotg->vec);
            
        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        for (m=0; m<DIM; m++)
        {
            tmp_f[m]             = rotg->k*dr[m];
            erg->f_rot_loc[i][m] = tmp_f[m];
            erg->V              += 0.5*rotg->k*sqr(dr[m]);
        }
        
        if (bTorque)
        {
            /* Add to the torque of this rotation group */
            erg->fix_torque_v += torque(rotg->vec, tmp_f, erg->x_loc_pbc[i], rotg->pivot);
            
            /* Calculate the angle between reference and actual rotation group atom.
             * If this routine returns something else than 0, the angle could not 
             * be determined (e.g. because actual or reference position is on
             * the rotation axis) */
            if (0 == angle(rotg, erg->x_loc_pbc[i], erg->xr_loc[i], &alpha, &weight, rotg->pivot))  /* angle in rad, weighted */
            {
                erg->fix_angles_v += alpha * weight;
                erg->fix_weight_v += weight;
                /* Use the next two lines instead if you don't want weighting: */
                /*
                angles_v[g] += alpha;
                weight_v[g] += 1;
                 */
            }
        }
        /* Probably one does not want enforced rotation to influence */
        /* the virial. But if so, activate the following lines */
        /*
            if (MASTER(cr))
            {
               Add the rotation contribution to the virial
              for(j=0; j<DIM; j++)
                for(m=0;m<DIM;m++)
                  vir[j][m] += 0.5*f[ii][j]*dr[m];
            }
         */
#ifdef INFOF
        fprintf(stderr," FORCE on ATOM %d = (%15.8f %15.8f %15.8f)  torque=%15.8f\n",
                ii,erg->f_rot_loc[i][XX], erg->f_rot_loc[i][YY], erg->f_rot_loc[i][ZZ],erg->fix_torque_v);
#endif
    } /* end of loop over local rotation group atoms */
}


/* Determine the smallest and largest position vector (with respect to the 
 * rotation vector) for the reference group */
static void get_firstlast_atom_ref(
        t_rotgrp  *rotg, 
        int       *firstindex, 
        int       *lastindex)
{
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    int i;
    real xcproj;               /* The projection of a reference position on the 
                                  rotation vector */
    real minproj, maxproj;     /* Smallest and largest projection on v */
    

    
    erg=rotg->enfrotgrp;

    /* Start with some value */
    minproj = iprod(erg->xc_ref[0], rotg->vec);
    maxproj = minproj;
    
    /* This is just to ensure that it still works if all the atoms of the 
     * reference structure are situated in a plane perpendicular to the rotation 
     * vector */
    *firstindex = 0;
    *lastindex  = rotg->nat-1;
    
    /* Loop over all atoms of the reference group, 
     * project them on the rotation vector to find the extremes */
    for (i=0; i<rotg->nat; i++)
    {
        xcproj = iprod(erg->xc_ref[i], rotg->vec);
        if (xcproj < minproj)
        {
            minproj = xcproj;
            *firstindex = i;
        }
        if (xcproj > maxproj)
        {
            maxproj = xcproj;
            *lastindex = i;
        }
    }
}


/* Allocate memory for the slabs */
static void allocate_slabs(
        t_rotgrp  *rotg, 
        FILE      *fplog, 
        int       g, 
        t_commrec *cr)
{
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */
    int i, nslabs;
    
    
    erg=rotg->enfrotgrp;
    
    /* More slabs than are defined for the reference are never needed */
    nslabs = erg->slab_last_ref - erg->slab_first_ref + 1;
    
    /* Remember how many we allocated */
    erg->nslabs_alloc = nslabs;

    if (MASTER(cr))
        fprintf(fplog, "Enforced rotation: allocating memory to store data for %d slabs (rotation group %d).\n",nslabs,g);
    snew(erg->slab_center    , nslabs);
    snew(erg->slab_center_ref, nslabs);
    snew(erg->slab_weights   , nslabs);
    snew(erg->slab_torque_v  , nslabs);
    snew(erg->slab_data      , nslabs);
    snew(erg->gn_atom        , nslabs);
    snew(erg->gn_slabind     , nslabs);
    for (i=0; i<nslabs; i++)
    {
        snew(erg->slab_data[i].x     , rotg->nat);
        snew(erg->slab_data[i].ref   , rotg->nat);
        snew(erg->slab_data[i].weight, rotg->nat);
    }
    snew(erg->xc_ref_sorted, rotg->nat);
    snew(erg->xc_sortind   , rotg->nat);
    snew(erg->firstatom    , nslabs);
    snew(erg->lastatom     , nslabs);
}


extern void init_rot_group(
        FILE *fplog,t_commrec *cr,
        int g,t_rotgrp *rotg,
        rvec *x,       /* the positions */
        gmx_mtop_t *mtop,
        FILE *out_slabs,
        matrix box)
{
    int i,ii;
    char        filename[255];/* File to save the reference positions in for enforced rotation */
    rvec        f_box[3];     /* Box from reference file */
    rvec        coord;
    rvec        center;       /* COG or COM of rotation group */
    t_trnheader header;       /* Header information of reference file */
    bool        bSame;        /* Used for a check if box sizes agree */
    bool        bFlex;        /* Flexible rotation? */
    t_atom      *atom;
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */
    int         ref_firstindex, ref_lastindex;
    rvec        *xdum;
    rvec        transvec;
    
    
    bFlex = (rotg->eType == erotgFLEX1 || rotg->eType == erotgFLEX2);
    
    erg=rotg->enfrotgrp;
    
    snew(erg->xc_ref    , rotg->nat);
    snew(erg->xc_ref_ind, rotg->nat);
    snew(erg->f_rot_loc , rotg->nat);
    
    /* Allocate space for the masses if needed */
    erg->mc=NULL;
    erg->m_loc=NULL;
    if (rotg->eOrigin == erotgOriginCOM)
    {
        snew(erg->mc, rotg->nat);
        if (!bFlex)
            snew(erg->m_loc, rotg->nat);
    }
    
    /* xc_ref_ind needs to be set to identity in the serial case */
    if (!PAR(cr))
        for (i=0; i<rotg->nat; i++)
            erg->xc_ref_ind[i] = i;

    /* Allocate space for collective coordinates if needed */
    if (bFlex)
    {
        snew(erg->xc        , rotg->nat);
        snew(erg->xc_old    , rotg->nat);
        snew(erg->xc_shifts , rotg->nat);
        snew(erg->xc_eshifts, rotg->nat);
        if (rotg->eFittype == erotgFitNORM)
        {
            snew(erg->xc_ref_length, rotg->nat); /* in case fit type NORM is chosen */
            snew(erg->xc_norm      , rotg->nat);
        }
    }
    else
    {
        snew(erg->xr_loc   , rotg->nat);
        snew(erg->x_loc_pbc, rotg->nat);
    }

    /* Read in rotation reference positions from file, if it exists.
     * If not, write out the initial rotation group positions as the reference set */
    if (MASTER(cr))
    {
        /* Save the reference positions to trr */
        /* Make a trr for each rotation group */
        sprintf(filename, "ref_%d_%s.trr", g, erotg_names[rotg->eType]);
        if (gmx_fexist(filename)) /* Read rotation reference positions from file */
        {
            fprintf(fplog, "Enforced rotation: found file %s containing reference positions.\n", filename);
            read_trnheader(filename, &header);
            if (rotg->nat != header.natoms)
                gmx_fatal(FARGS,"Number of atoms in file %s (%d) does not match the number of atoms in rotation group (%d)!\n",
                        filename, header.natoms, rotg->nat);
            read_trn(filename, &header.step, &header.t, &header.lambda, f_box, &header.natoms, erg->xc_ref, NULL, NULL);
            fprintf(fplog , "Enforced rotation: read reference positions for group %d from %s.\n", g, filename);
            /* Check if the box is unchanged and output a warning if not: */
            bSame = TRUE;
            for (i=0; i<DIM; i++)
                for (ii=0; ii<DIM; ii++)
                    if (f_box[i][ii] != box[i][ii]) bSame = FALSE;

            if (!bSame)
            {
                sprintf(warn_buf, "Enforced rotation: Box size in reference file %s differs from actual box size!", filename);
                warning(NULL);
                pr_rvecs(stderr,0,"Your box is:",box  ,3);
                pr_rvecs(fplog ,0,"Your box is:",box  ,3);
                pr_rvecs(stderr,0,"Box in file:",f_box,3);
                pr_rvecs(fplog ,0,"Box in file:",f_box,3);
            }

            if (g != header.step)
            {   /* We use step to indicate the number of the rotation group here */
                sprintf(warn_buf,"Coordinates from %s will be used for rotation group %d", filename, g);
                warning(NULL);
            }
        }
        else /* Save the initial positions of the rotation group as reference */
        {
            for(i=0; i<rotg->nat; i++)
            {
                ii = rotg->ind[i];
                copy_rvec(x[ii], erg->xc_ref[i]);
            }
            write_trn(filename,g,0.0,0.0,box,rotg->nat,erg->xc_ref,NULL,NULL);
            fprintf(fplog, "Enforced rotation: saved %d coordinates of group %d to %s.\n",
                    rotg->nat, g, filename);
        }
        
        /* Follow the center of mass? Then we need the masses! */
        if (rotg->eOrigin == erotgOriginCOM)
        {
            /* We need to copy the masses to be able to determine the COM */
            for (i=0; i<rotg->nat; i++)
            {
                gmx_mtop_atomnr_to_atom(mtop,rotg->ind[i],&atom);
                erg->mc[i] = atom->m;
            }
        }
        
        /* When COM or COG is subtracted from the positions, the centers of 
         * reference and initial groups need to be determined */
        if (rotg->eOrigin != erotgOriginBox)
        {
            /* Save the center of the REFERENCE structure: */
            get_center(erg->xc_ref, erg->mc, rotg->nat, center);
            fprintf(fplog, "Enforced rotation: group %d reference %s =%14.7e %14.7e %14.7e is subtracted from ref. positions\n", g,
                    rotg->eOrigin==erotgOriginCOG? "COG":"COM", center[XX], center[YY], center[ZZ]);
            /* Subtract the center from the reference positions */
            svmul(-1.0, center, transvec);
            translate_x(erg->xc_ref, rotg->nat, transvec);

            /* For fixed rotation we need the center of the initial positions */
            if (!bFlex)
            {            
                snew(xdum, rotg->nat);            
                for (i=0; i<rotg->nat; i++)
                {
                    ii = rotg->ind[i];
                    copy_rvec(x[ii], xdum[i]);
                }
                get_center(xdum, erg->mc, rotg->nat, erg->center);
                sfree(xdum);
                fprintf(fplog, "Enforced rotation: group %d initial   %s =%14.7e %14.7e %14.7e\n", g,
                        rotg->eOrigin==erotgOriginCOG? "COG":"COM", erg->center[XX], erg->center[YY], erg->center[ZZ]);
            }
        }
        
        if (bFlex)
        {
            /* Save the original (whole) set of positions such that later the
             * molecule can always be made whole again */
            for (i=0; i<rotg->nat; i++)
            {
                ii = rotg->ind[i];
                copy_rvec(x[ii], erg->xc_old[i]);
            }
        }
    }
#ifdef GMX_MPI
    /* Broadcast from master node to all PP nodes */
    if (PAR(cr))
    {
        /* Broadcast reference positions to all PP nodes */
        gmx_bcast(rotg->nat*sizeof(erg->xc_ref[0]), erg->xc_ref, cr);
        /* Broadcast other stuff */
        if (bFlex)
            gmx_bcast(rotg->nat*sizeof(erg->xc_old[0]),erg->xc_old, cr);
        else
            gmx_bcast(sizeof(erg->center), erg->center, cr);
        
        /* Broadcast masses */
        if (rotg->eOrigin == erotgOriginCOM)
            gmx_bcast(rotg->nat*sizeof(erg->mc[0]), erg->mc, cr);
    }
#endif
    /* Enforced rotation with flexible axis */
    if (bFlex)
    {
        /* Calculate maximum beta value from minimum gaussian (performance opt.) */
        erg->max_beta = calc_beta_max(rotg->min_gaussian, rotg->slab_dist);

        /* Determine the smallest and larges coordinate with respect to the rotation vector */
        get_firstlast_atom_ref(rotg, &ref_firstindex, &ref_lastindex);
        
        /* From the extreme coordinates of the reference group, determine the first and last slab */
        erg->slab_first_ref = get_first_slab(rotg, erg->max_beta, erg->xc_ref[ref_firstindex]);
        erg->slab_last_ref  = get_last_slab( rotg, erg->max_beta, erg->xc_ref[ref_lastindex ]);

        /* Add some slabs on both sides such that the flexible group can move a bit */
        erg->slab_buffer = (erg->slab_last_ref - erg->slab_first_ref + 1) / 2;
        erg->slab_first_ref = erg->slab_first_ref - erg->slab_buffer;
        erg->slab_last_ref  = erg->slab_last_ref  + erg->slab_buffer;
        
        /* Allocate memory for the slabs */
        allocate_slabs(rotg, fplog, g, cr);

        /* Flexible rotation: determine the reference COGs for the rest of the simulation */
        erg->slab_first = erg->slab_first_ref;
        erg->slab_last = erg->slab_last_ref;
        get_slab_centers(rotg,erg->xc_ref,cr,g,-1,out_slabs,TRUE,TRUE);

        /* Also save the center of geometry of the reference structure (only needed for fitting): */
        get_center(erg->xc_ref, NULL, rotg->nat, erg->xc_ref_center);

        /* Length of each x_rotref vector from center (needed if fit routine NORM is chosen): */
        if (rotg->eFittype == erotgFitNORM)
        {
            for (i=0; i<rotg->nat; i++)
            {
                rvec_sub(erg->xc_ref[i], erg->xc_ref_center, coord);
                erg->xc_ref_length[i] = norm(coord);
            }
        }
    }
}


static void make_local_rotation_group(gmx_ga2la_t ga2la,
        t_rotgrp *rotg,int start,int end)
{
    int i,ii;
    gmx_enfrotgrp_t erg;      /* Pointer to enforced rotation group data */

    
    erg=rotg->enfrotgrp;

    erg->nat_loc = 0;
    for(i=0; i<rotg->nat; i++) {
        ii = rotg->ind[i];
        if (!ga2la_home(ga2la,ii,&ii)) {
            ii = -1;
        }

        if (ii >= start && ii < end) {
            /* This is a home atom, add it to the local rotation group */
            if (erg->nat_loc >= erg->nalloc_loc)
            {
                erg->nalloc_loc = over_alloc_dd(erg->nat_loc+1);
                srenew(erg->ind_loc,erg->nalloc_loc);
            }
            erg->ind_loc[erg->nat_loc] = ii;
            /* Copy the reference positions */
            /* Remember which of the x_rotref atoms are local: */
            erg->xc_ref_ind[erg->nat_loc]=i;  /* i is the number of the atom with respect to the whole rotation group */
            /* erg->ind[i] is the number with respect to the whole system! */
            erg->nat_loc++;
        }
    }
}

void dd_make_local_rotation_groups(gmx_domdec_t *dd,t_rot *rot,t_mdatoms *md)
{
    gmx_ga2la_t ga2la;
    int g;

    
    ga2la = dd->ga2la;

    for(g=0; g<rot->ngrp; g++)
        make_local_rotation_group(ga2la,&rot->grp[g],md->start,md->start+md->homenr);
}


void init_rot(FILE *fplog,t_inputrec *ir,int nfile,const t_filenm fnm[],
        t_commrec *cr, matrix box, rvec *x, gmx_mtop_t *mtop,
        const output_env_t oenv, unsigned long Flags)
{
    t_rot    *rot;
    t_rotgrp *rotg;
    int      g;
    int      nat_max=0;     /* Size of biggest rotation group */
    bool     bRerun;
    bool     bFlex=FALSE;
    gmx_enfrot_t er;        /* Pointer to the enforced rotation buffer variables */    
    gmx_enfrotgrp_t erg;    /* Pointer to enforced rotation group data */


    if ( (PAR(cr)) && !DOMAINDECOMP(cr) )
        gmx_fatal(FARGS, "Enforced rotation is only implemented for domain decomposition!");

    rot = ir->rot;
    snew(rot->enfrot, 1);
    er=rot->enfrot;
    
    /* Output every step for reruns */
    bRerun = Flags & MD_RERUN;
    if (bRerun)
    {
        if (fplog)
            fprintf(fplog, "Enforced rotation: rerun - will write rotation output every available step.\n");
        rot->nstrout = 1;
        rot->nsttout = 1;
    }

    er->out_slabs = NULL;
    if (MASTER(cr))
        er->out_slabs = open_slab_out(rot, opt2fn("-rs",nfile,fnm));

    for(g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];

        if (fplog)
            fprintf(fplog,"Enforced rotation: group %d type '%s'\n", g, erotg_names[rotg->eType]);

        if (rotg->eType == erotgFLEX1 || rotg->eType == erotgFLEX2)
            bFlex = TRUE;
        
        if (rotg->nat > 0)
        {
            /* Allocate space for the rotation group's data: */
            snew(rotg->enfrotgrp, 1);
            erg  = rotg->enfrotgrp;

            nat_max=max(nat_max, rotg->nat);
            
            if (PAR(cr))
            {
                erg->nat_loc    = 0;
                erg->nalloc_loc = 0;
                erg->ind_loc    = NULL;
            }
            else
            {
                erg->nat_loc = rotg->nat;
                erg->ind_loc = rotg->ind;
            }
            init_rot_group(fplog,cr,g,rotg,x,mtop,er->out_slabs,box);
        }
    }
    
    /* Allocate space for enforced rotation buffer variables */
    er->bufsize = nat_max;
    snew(er->data, nat_max);
    snew(er->buf , nat_max);

    /* Buffers for MPI reducing torques, angles, weights (for each group), and V */
    er->mpi_bufsize = 4*rot->ngrp; /* To start with */
    snew(er->mpi_inbuf , er->mpi_bufsize);
    snew(er->mpi_outbuf, er->mpi_bufsize);

    /* Only do I/O on the MASTER */
    er->out_angles  = NULL;
    er->out_rot     = NULL;
    er->out_torque  = NULL;
    if (MASTER(cr))
    {
        er->out_rot = open_rot_out(opt2fn("-ro",nfile,fnm), rot, oenv, Flags);
        if (bFlex)
        {
            if (rot->nstrout > 0)
                er->out_angles  = open_angles_out(rot, opt2fn("-ra",nfile,fnm));
            if (rot->nsttout > 0)
                er->out_torque  = open_torque_out(rot, opt2fn("-rt",nfile,fnm));
        }
    }
}


void finish_rot(FILE *fplog,t_rot *rot)
{
    gmx_enfrot_t er;        /* Pointer to the enforced rotation buffer variables */    

    
    er=rot->enfrot;
    if (er->out_rot)
        gmx_fio_fclose(er->out_rot);
    if (er->out_slabs)
        gmx_fio_fclose(er->out_slabs);
    if (er->out_angles)
        gmx_fio_fclose(er->out_angles);
    if (er->out_torque)
        gmx_fio_fclose(er->out_torque);
}


/* Rotate the local reference positions and store them in erg->xr_loc[0...(nat_loc-1)]*/
static void rotate_local_reference(t_rotgrp *rotg)
{
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec      xr, xrcpy;       /* rotated (reference) position */
    int i,iigrp;

    
    erg=rotg->enfrotgrp;
    
    for (i=0; i<erg->nat_loc; i++)
    {
        /* Index of this rotation group atom with respect to the whole rotation group */
        iigrp = erg->xc_ref_ind[i];
        /* Rotate this atom around dislocated rotation axis: */
        /* Move rotation axis, so that it runs through the origin: */
        rvec_sub(erg->xc_ref[iigrp], rotg->pivot, xr);
        /* Rotate around the origin: */
        mvmul(erg->rotmat, xr, xrcpy);
        /* Move back and store for later usage */
        rvec_add(xrcpy, rotg->pivot, erg->xr_loc[i]);
    }
}


/* Select the PBC representation for each local x position and store that
 * for later usage. The right PBC image of an x is the one nearest to its 
 * rotated reference */
static void choose_pbc_image(rvec x[], t_rotgrp *rotg, matrix box, int npbcdim)
{
    int d,i,ii,m;
    gmx_enfrotgrp_t erg;       /* Pointer to enforced rotation group data */
    rvec xref,xcurr,dx;
    ivec shift;
    
    
    erg=rotg->enfrotgrp;
    
    for (i=0; i<erg->nat_loc; i++)
    {
        clear_ivec(shift);
        
        /* Index of a rotation group atom  */
        ii = erg->ind_loc[i];

        /* Get the already rotated reference position */
        copy_rvec(erg->xr_loc[i], xref);
        
        /* Subtract the (old) center from the current positions
         * (just to determine the shifts!) */
        rvec_sub(x[ii], erg->center, xcurr);
        
        /* Shortest PBC distance between the atom and its reference */
        rvec_sub(xref, xcurr, dx);
        
        /* Determine the shift for this atom */
        /* TODO: check whether it is faster to save the shifts and
         * only update them in DD/NS steps */
        for(m=npbcdim-1; m>=0; m--)
        {
            while (dx[m] < -0.5*box[m][m])
            {
                for(d=0; d<DIM; d++)
                    dx[d] += box[m][d];
                shift[m]++;
            }
            while (dx[m] >= 0.5*box[m][m])
            {
                for(d=0; d<DIM; d++)
                    dx[d] -= box[m][d];
                shift[m]--;
            }
        }
            
        /* Apply the shift to the current atom */
        copy_rvec(x[ii], erg->x_loc_pbc[i]);
        shift_single_coord(box, erg->x_loc_pbc[i], shift); 
    }
}


extern void do_rotation(
        t_commrec *cr,
        t_inputrec *ir,
        matrix box,
        rvec x[],
        real t,
        int step,
        gmx_wallcycle_t wcycle,
        bool bNS)
{
    int      g,i,ii;
    t_pbc    pbc;
    t_rot    *rot;
    t_rotgrp *rotg;
    bool     outstep_torque;
    float    cycles_rot;
    gmx_enfrot_t er;     /* Pointer to the enforced rotation buffer variables */
    gmx_enfrotgrp_t erg; /* Pointer to enforced rotation group data           */
    rvec     transvec;
#ifdef TAKETIME
    double t0;
#endif
    
    
    rot=ir->rot;
    er=rot->enfrot;
    
    /* At which time steps do we want to output the torque */
    outstep_torque = do_per_step(step, rot->nsttout);
    
    /* Output time into rotation output file */
    if (outstep_torque && MASTER(cr))
        fprintf(er->out_rot, "%12.3e",t);

    /* First do ALL the communication! */
    for(g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        erg=rotg->enfrotgrp;

        /* Calculate the rotation matrix for this angle: */
        erg->degangle = rotg->rate * t;
        calc_rotmat(rotg->vec,erg->degangle,erg->rotmat);

        switch (rotg->eType)
        {
        case erotgFIXED:
        case erotgFIXED_PLANE:
            /* Rotate the local reference positions */
            rotate_local_reference(rotg);
            /* Choose the nearest PBC images of the group atoms with respect 
             * to the rotated reference positions */
            set_pbc(&pbc,ir->ePBC,box);
            choose_pbc_image(x, rotg, pbc.box, 3);
            /* For some cases we need the center of the local positions, which
             * invokes communication */
            if (rotg->eOrigin == erotgOriginBox)
                clear_rvec(erg->center);
            else
            {
                /* Prepare a local masses array; this array 
                 * changes in DD/neighborsearching steps */
                if (bNS && erg->m_loc != NULL)
                {
                    for (i=0; i<erg->nat_loc; i++)
                    {
                        /* Index of local atom w.r.t. the collective rotation group */
                        ii = erg->xc_ref_ind[i];
                        erg->m_loc[i] = erg->mc[ii];
                    }
                }
                get_center_comm(cr, erg->x_loc_pbc, erg->m_loc, erg->nat_loc, rotg->nat, erg->center);
            }
            break;
        case erotgFLEX1:
        case erotgFLEX2:
            /* Transfer the rotation group's positions such that every node has 
             * all of them. Every node contributes its local positions x and stores 
             * it in the collective erg->xc array. */
            communicate_group_positions(cr,erg->xc, erg->xc_shifts, erg->xc_eshifts, bNS,
                    x, rotg->nat, erg->nat_loc, erg->ind_loc, erg->xc_ref_ind, erg->xc_old, box);
            break;
        default:
            gmx_fatal(FARGS, "Enforced rotation: no such rotation type");
            break;
        }
    }
    
    /* Done communicating, we can start to count cycles now ... */
    wallcycle_start(wcycle, ewcROT);
    GMX_MPE_LOG(ev_rotcycles_start);
    
#ifdef TAKETIME
    t0 = MPI_Wtime();
#endif
    
    for(g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        erg=rotg->enfrotgrp;

        if (outstep_torque && MASTER(cr))
            fprintf(er->out_rot, "%12.4f", erg->degangle);

        switch(rotg->eType)
        {
        case erotgFIXED:
        case erotgFIXED_PLANE:
            do_fixed(cr,rotg,x,&pbc,box,t,step,outstep_torque);
            break;
        case erotgFLEX1:
        case erotgFLEX2:
            if (rotg->eOrigin != erotgOriginBox)
            {                
                /* Subtract the center of the rotation group */
                get_center(erg->xc, erg->mc, rotg->nat, erg->center);
                /* TODO: output the center somewhere */
                svmul(-1.0, erg->center, transvec);
                translate_x(erg->xc, rotg->nat, transvec);
            }
            do_flexible(cr,er,rotg,g,x,box,t,step,outstep_torque);
            break;
        default:
            break;
        }
    }

#ifdef TAKETIME
    if (MASTER(cr))
        fprintf(stderr, "Enforced rotation calculation (step %d) took %g seconds.\n", step, MPI_Wtime()-t0);
#endif

    /* Stop the cycle counter and add to the force cycles for load balancing */
    cycles_rot = wallcycle_stop(wcycle,ewcROT);
    if (DOMAINDECOMP(cr) && wcycle)
        dd_cycles_add(cr->dd,cycles_rot,ddCyclF);
    GMX_MPE_LOG(ev_rotcycles_finish);
}
