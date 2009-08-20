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
#include "edsam.h"


/* Enforce rotation / flexible: determine the angle of each slab */
typedef struct gmx_slabdata
{
    int  nat;     /* The number of coordinates belonging to this slab */
    rvec *x;      /* The coordinates belonging to this slab. In general, this should be all
                   * rotation group coordinates, but we can leave a few of them away if they have
                   * small enough weights. */
    rvec *ref;    /* Same for reference */
    real *weight; /* The weight for each atom */
} t_gmx_slabdata;


static double**  allocate_square_matrix(int dim)
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
    int      g,i,islab,k,nslabs=0;
    int      count; /* MPI element counter */
    t_rotgrp *rotg;
    
    
    /* Fill the MPI buffer with stuff to reduce: */
    if (PAR(cr))
    {
        count=0;
        for (g=0; g < rot->ngrp; g++)
        {
            rotg = &rot->grp[g];
            nslabs = 2*rotg->slab_max_nr+1;
            rot->inbuf[count++] = rotg->V;
            switch (rotg->eType)
            {
            case erotgFIXED:
                rot->inbuf[count++] = rotg->fix_torque_v;
                rot->inbuf[count++] = rotg->fix_angles_v;
                rot->inbuf[count++] = rotg->fix_weight_v;
                break;
            case erotgFLEX:
            case erotgFLEX2:
                /* (Re-)allocate memory for MPI buffer: */
                if (rot->bufsize < count+nslabs)
                {
                    rot->bufsize = count+nslabs;
                    srenew(rot->inbuf , rot->bufsize);
                    srenew(rot->outbuf, rot->bufsize);
                }
                for (i=0; i<nslabs; i++)
                    rot->inbuf[count++] = rotg->slab_torque_v[i];
                break;
            default:
                break;
            }
        }
#ifdef GMX_MPI
        MPI_Reduce(rot->inbuf, rot->outbuf, count, GMX_MPI_REAL, MPI_SUM, MASTERRANK(cr), cr->mpi_comm_mygroup);
#endif
        /* Copy back the reduced data from the buffer on the master */
        if (MASTER(cr))
        {
            count=0;
            for (g=0; g < rot->ngrp; g++)
            {
                rotg = &rot->grp[g];
                nslabs = 2*rotg->slab_max_nr+1;
                rotg->V = rot->outbuf[count++];
                switch (rotg->eType)
                {
                case erotgFIXED:
                    rotg->fix_torque_v = rot->outbuf[count++];
                    rotg->fix_angles_v = rot->outbuf[count++];
                    rotg->fix_weight_v = rot->outbuf[count++];
                    break;
                case erotgFLEX:
                case erotgFLEX2:
                    for (i=0; i<nslabs; i++)
                        rotg->slab_torque_v[i] = rot->outbuf[count++];
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
            
            /* Output to main rotation log file: */
            if (rotg->eType == erotgFIXED)
            {
                fprintf(rot->out_rot, "%12.4f%12.3e", 
                        (rotg->fix_angles_v/rotg->fix_weight_v)*180.0*M_1_PI,
                        rotg->fix_torque_v);
            }
            fprintf(rot->out_rot, "%12.3e", rotg->V);
                        
            /* Output to torque log file: */
            if (rotg->eType == erotgFLEX || rotg->eType == erotgFLEX2)
            {
                fprintf(rot->out_torque, "%12.3e%6d", t, g);
                k = rotg->slab_max_nr;
                for (i=-k; i <= k; i++)
                {
                    islab = i + rotg->slab_max_nr;  /* slab index */
                    /* Only output if enough weight is in slab */
                    if (rotg->slab_weights[islab] > rotg->min_gaussian)
                        fprintf(rot->out_torque, "%6d%12.3e", i, rotg->slab_torque_v[islab]);
                }
                fprintf(rot->out_torque , "\n");
            }
        }
        fprintf(rot->out_rot, "\n");
    }
}


/* Add the forces from enforced rotation potential to the local forces.
 * Should be called after the SR forces have been evaluated */
extern real add_rot_forces(t_rot *rot, rvec f[], t_commrec *cr, int step, real t)
{
    int g,l,ii;
    t_rotgrp *rotg;
    

    GMX_MPE_LOG(ev_add_rot_forces_start);
    
    /* Reduce energy,torque, angles etc. to get the sum values (per rotation group) 
     * on the master and output these values to file. */
    if (do_per_step(step, rot->nsttout))
        reduce_output(cr, rot, t);

    /* Total rotation potential is the sum over all rotation groups */
    rot->Vrot = 0.0; 
        
    /* Loop over enforced rotation groups (usually 1, though)
     * Apply the forces from rotation potentials */
    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        rot->Vrot += rotg->V;
        for (l=0; l<rotg->nat_loc; l++)
        {
            /* Get the right index of the local force */
            ii = rotg->ind_loc[l];
            /* Add */
            rvec_inc(f[ii],rotg->f_rot_loc[l]);
        }
    }
    
    GMX_MPE_LOG(ev_add_rot_forces_finish);

    return (MASTER(cr)? rot->Vrot : 0.0);
}


/* Calculate the box diagonal length */
static real diagonal_length(matrix box)
{
    rvec diag;

    
    copy_rvec(box[XX],diag);
    rvec_inc(diag,box[YY]);
    rvec_inc(diag,box[ZZ]);

    return norm(diag);
}


static inline real calc_beta(rvec curr_x, t_rotgrp *rotg, int n)
{
    return iprod(curr_x, rotg->vec) - rotg->slab_dist * n;
}


static inline real gaussian_weight(rvec curr_x, t_rotgrp *rotg, int n)
{
    /* norm is chosen such that the sum of the gaussians
     * over the slabs is approximately 1.0 everywhere */
    const real norm =  0.5698457353514458216;  /* = 1/1.7548609 */
    real        sigma;

    
    /* Define the sigma value */
    sigma = 0.7*rotg->slab_dist;
    /* Calculate the Gaussian value of slab n for coordinate curr_x */
    return norm * exp( -0.5 * sqr( calc_beta(curr_x, rotg, n)/sigma ) );
}


static void get_slab_centers(
        t_rotgrp *rotg,       /* The rotation group information */
        rvec      *xc,        /* The rotation group coordinates; will typically be pgrp->xc,
                               * but at first call it is pgrp->xc_ref */
        matrix    box,        /* The box coordinates        */
        t_commrec *cr,        /* Communication record       */
        int       g,          /* The rotation group number  */
        bool      bDynBox,    /* Is the box dynamic?        */
        real      time,       /* Used for output only       */
        FILE      *out_slabs, /* For outputting COG per slab information */
        bool      bOutStep,   /* Is this an output step?    */
        bool      bReference) /* If this routine is called from
                               * init_rotation_group_rot we need to store
                               * the reference slab COGs */
{
    rvec curr_x;          /* The coordinate of an atom        */
    rvec curr_x_weighted; /* The gaussian-weighted coordinate */
    real gaussian;        /* The gaussian                     */
    int i,j,k,islab,nslabs;
    real box_d;           /* The box diagonal                 */
    bool bFirstSet;


    /* If the box grows during the simulation, we might need more memory */
    if (bDynBox)
    {
        box_d = diagonal_length(box);

        /* The slab indices run from [-pgrp->slab_max_nr, -1, 0, +1, ..., +pgrp->slab_max_nr] */
        nslabs = 2*rotg->slab_max_nr + 1;

        /* The box diagonal divided by the slab distance gives the maximum number of slabs in positive direction: */
        if ( (int)ceil(box_d/rotg->slab_dist) > rotg->slab_max_nr )
        {
            while ( (int)ceil(box_d/rotg->slab_dist) > rotg->slab_max_nr )
                rotg->slab_max_nr++;
            /* TODO: it could still be that the rotation group diffuses out of the
             * box. Then we would have to allocate more slabs than fit in a box!
             */

            nslabs = 2*rotg->slab_max_nr + 1;
            fprintf(stdout, "Node %d reallocates memory to hold data for %d slabs (rotation group %d).\n", cr->nodeid,nslabs,g);
            srenew(rotg->slab_center  , nslabs);
            srenew(rotg->slab_weights , nslabs);
            srenew(rotg->slab_torque_v, nslabs);
            srenew(rotg->slab_data    , nslabs);
            for (i=0; i<nslabs; i++)
            {
                srenew(rotg->slab_data[i].x     , rotg->nat);
                srenew(rotg->slab_data[i].ref   , rotg->nat);
                srenew(rotg->slab_data[i].weight, rotg->nat);
            }
        }
    }

    /* Loop over slabs */
    bFirstSet = FALSE;
    k = rotg->slab_max_nr;
    rotg->slab_first = 1;
    rotg->slab_last  = 0;
    for (j = -k; j <= k; j++)
    {
        islab = j+rotg->slab_max_nr; /* slab index */
        /* Initialize data for this slab: */
        clear_rvec(rotg->slab_center[islab]);
        rotg->slab_weights[islab] = 0.0;

        /* loop over all atoms in the rotation group */
        for(i=0; i<rotg->nat;i++)
        {
            copy_rvec(xc[i], curr_x);
            gaussian = gaussian_weight(curr_x, rotg, j);
            svmul(gaussian, curr_x, curr_x_weighted);
            rvec_add(rotg->slab_center[islab], curr_x_weighted, rotg->slab_center[islab]);
            rotg->slab_weights[islab] += gaussian;
        } /* END of loop over rotation group atoms */

        /* Do the calculations ONLY if there is enough weight in the slab! */
        if (rotg->slab_weights[islab] > rotg->min_gaussian)
        {
            svmul(1.0/rotg->slab_weights[islab], rotg->slab_center[islab], rotg->slab_center[islab]);
            /* Remember which slabs to calculate for the low-level routines */
            if (!bFirstSet)
            {
                rotg->slab_first = j;
                bFirstSet = TRUE;
            }
            rotg->slab_last = j;
        }
        /* At first time step: save the COGs of the reference structure */
        if(bReference)
            copy_rvec(rotg->slab_center[islab], rotg->slab_center_ref[islab]);
    } /* END of loop over slabs */
    
    /* Output on the master */
    if (MASTER(cr) && bOutStep)
    {
        fprintf(out_slabs, "%12.3e", time);
        for (j = rotg->slab_first; j <= rotg->slab_last; j++)
        {
            islab = j+rotg->slab_max_nr; /* slab index */
            fprintf(out_slabs, "%6d%12.3e%12.3e%12.3e",
                    j,rotg->slab_center[islab][XX],rotg->slab_center[islab][YY],rotg->slab_center[islab][ZZ]);
        }
        if (!bFirstSet)
            fprintf(out_slabs, "WARNING: no weight in any of the slabs - nothing to calculate!");
        fprintf(out_slabs, "\n");
    }
}


static void calc_rotmat(
        rvec vec,
        real degangle,  /* Angle alpha of rotation at time t in degrees */
        matrix rotmat)  /* Rotation matrix */
{
    real radangle;            /* Rotation angle in radians */
    real cosa;                /* cosine alpha   */
    real sina;                /* sine alpha     */
    real OMcosa;              /* 1 - cos(alpha) */
    real dumxy, dumxz, dumyz; /* save computations */
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


/* Calculates torque m = r x f_perp */
static inline real torque(rvec v,  /* IN:  axis of rotation */
        rvec f,  /* IN:  force */
        rvec x,  /* IN:  coordinate of atom on which the force acts */
        rvec c)  /* IN:  center of mass of the slab (through this point
                                                     the axis of rotation goes */
{
    rvec f_perp,r,e;
    real alpha;
    rvec m, vectmp;


    /* Calculate alpha = v(x-c)/(v*v) */
    rvec_sub(x, c, vectmp);
    alpha = iprod(v, vectmp)/norm2(v); /* Actually v should be a unit vector anyway */

    /* Calculate r = x-a = x-(c+alpha*v) */
    svmul(alpha, v, vectmp);  /*  vectmp = alpha*v      */
    rvec_inc(vectmp, c);      /*  vectmp = alpha*v + c  */
    rvec_sub(x, vectmp, r);   /*  r      = x - a        */

    /* Calculate e, which is a unit vector in the direction of f_perp: */
    cprod(v, r, vectmp);
    unitv(vectmp, e);

    /* Calculate f_perp: */
    svmul(iprod(f,e),e,f_perp);

    /* Calculate torque m: */
    cprod(r, f_perp, m);

    /* Return the scalar value of m which is in the direction of the rotation axis v: */
    return iprod(m, v);
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


/* Open output file for slab COG data. Call on master only */
static FILE *open_slab_out(t_rot *rot)
{
    FILE      *fp=NULL;
    int       g;
    t_rotgrp  *rotg;


    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (rotg->eType == erotgFLEX || rotg->eType == erotgFLEX2)
        {
            if (NULL == fp)
                fp = ffopen("slabCOGs.log", "w");                
            fprintf(fp, "%% Rotation group %d (%s), slab distance %f nm\n", g, erotg_names[rotg->eType], rotg->slab_dist);
        }
    }
    
    if (fp != NULL)
    {
        fprintf(fp, "%% The following columns will have the syntax: (COG = center of geometry, gaussian weighted)\n");
        fprintf(fp, "%%     ");
        print_aligned_short(fp, "t");
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
static FILE *open_rot_out(t_rot *rot)
{
    FILE      *fp;
    int       g;
    t_rotgrp  *rotg;


    fp = ffopen("rotation.log", "w");
    fprintf(fp, "%% Output is written every %d time steps.\n", rot->nsttout);
    fprintf(fp, "%%\n");
    fprintf(fp, "%% The scalar tau is the torque in the direction of the rotation vector v.\n");
    fprintf(fp, "%% To obtain the vectorial torque, multiply tau with the group's rot_vec.\n");
    fprintf(fp, "%% Torques are given in [kJ/mol], anlges in degrees, time in ps.\n");

    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        fprintf(fp, "%% Rotation group %d (%s):\n", g, erotg_names[rotg->eType]);
        fprintf(fp, "%% rot_vec%d            %10.3e %10.3e %10.3e\n", g, rotg->vec[XX], rotg->vec[YY], rotg->vec[ZZ]);
        fprintf(fp, "%% rot_rate%d           %10.3e degree/ps\n",     g, rotg->rate);
        fprintf(fp, "%% rot_k%d              %10.3e kJ/(mol*nm^2)\n", g, rotg->k);

        switch (rotg->eType)
        {
        case erotgFIXED:
            fprintf(fp, "%% rot_offset%d         %10.3e %10.3e %10.3e\n", g, rotg->offset[XX], rotg->offset[YY], rotg->offset[ZZ]);
            break;
        case erotgFLEX:
        case erotgFLEX2:
            fprintf(fp, "%% rot_slab_distance%d   %f nm\n", g, rotg->slab_dist);
            fprintf(fp, "%% rot_min_gaussian%d   %10.3e\n", g, rotg->min_gaussian);
            break;
        default:
            break;
        }
    }
    
    fprintf(fp, "%%     ");
    print_aligned_short(fp, "t");
    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        print_aligned_group(fp, "theta_ref", g);
    }
    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (rotg->eType == erotgFIXED)
        {
            print_aligned_group(fp, "theta_av", g);
            print_aligned_group(fp, "tau", g);
        }
        print_aligned_group(fp, "energy", g);
    }
    fprintf(fp, "\n");
    fflush(fp);

    return fp;
}


/* Call on master only */
static FILE *open_angles_out(t_rot *rot, char filename[])
{
    int      g;
    FILE     *fp=NULL;
    t_rotgrp *rotg;


    /* Open output file and write some information about it's structure: */
    fp = ffopen(filename, "w");

    fprintf(fp, "%% Output will be written every %d steps\n", rot->nstrout);
    fprintf(fp, "%% All angles given in degrees, time in ps\n");
    fprintf(fp, "%% If more than one rotation group is present, each will appear in a separate line.\n");
    fprintf(fp, "%% n rotation groups will result in n lines of output all with the same time.\n");
    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (rotg->eType == erotgFLEX || rotg->eType == erotgFLEX2)
            fprintf(fp, "%% Rotation group %d (%s), slab distance %f nm\n", g, erotg_names[rotg->eType], rotg->slab_dist);
    }
    fprintf(fp, "%% The following columns will have the syntax:\n");
    fprintf(fp, "%%     ");
    print_aligned_short(fp, "t");
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
static FILE *open_torque_out(t_rot *rot)
{
    FILE      *fp;
    int       g;
    t_rotgrp  *rotg;


    fp = ffopen("torque.log", "w");
    fprintf(fp, "%% Output will be written every %d steps\n", rot->nsttout);

    for (g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        if (rotg->eType == erotgFLEX || rotg->eType == erotgFLEX2)
        {
            fprintf(fp, "%% Rotation group %d (%s), slab distance %f nm\n", g, erotg_names[rotg->eType], rotg->slab_dist);
            fprintf(fp, "%% The scalar tau is the torque [kJ/mol] in the direction of the rotation vector.\n");
            fprintf(fp, "%% To obtain the vectorial torque, multiply tau with\n");
            fprintf(fp, "%% rot_vec%d            %10.3e %10.3e %10.3e\n", g, rotg->vec[XX], rotg->vec[YY], rotg->vec[ZZ]);
            fprintf(fp, "%%\n");
        }
    }
    fprintf(fp, "%% The following columns will have the syntax (tau=torque for that slab):\n");
    fprintf(fp, "%%     ");
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


/* Determine geometrical center of structure with coordinates x */
static void get_center(rvec x[], real weight[], int nat, rvec center)
{
    int i;
    rvec coord;
    real weight_sum = 0.0;


    /* Zero out the center of mass */
    clear_rvec(center);

    /* Loop over all atoms and add their weighted position vectors */
    if (weight)
    {
        for (i=0; i<nat; i++)
        {
            weight_sum += weight[i];
            svmul(weight[i], x[i], coord);
            rvec_inc(center, coord);
        }

        /* Divide by the sum of weight */
        svmul(1.0/weight_sum, center, center);
    }
    else
    {
        for (i=0; i<nat; i++)
            rvec_inc(center, x[i]);

        /* Divide by the number of atoms */
        svmul(1.0/nat, center, center);
    }
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


static void weight_coords(rvec* str, real* weight, int natoms)
{
    int i, j;
    
    
    for(i=0; i<natoms; i++)
    {
        for(j=0; j<3; j++)
            str[i][j] *= sqrt(weight[i]);
    }  
}


static void trans(rvec* x, int natoms, double* vec)
{
    int i;
    
    
    for(i=0; i<natoms; i++)
    {
        x[i][0] += vec[0];
        x[i][1] += vec[1];
        x[i][2] += vec[2];
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
    double shift[3];
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
    for (i=0; i<3; i++)
        shift[i] = (-1.0)*ref_com[i];   
    trans(ref_s_1, natoms, shift);
    
    for (i=0; i<3; i++)
        shift[i] = (-1.0)*act_com[i];
    trans(act_s_1, natoms, shift);
    
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
    
    /* Weight coordinates with sqrt(weight) */
    if (weight)
    {
        weight_coords(ref_s_1, weight, natoms);
        weight_coords(act_s_1, weight, natoms);
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
        t_rotgrp *rotg,
        double t,
        real degangle,
        FILE *fp,
        FILE *fp_ref)
{
    int         i,l,n,islab,ind;
    rvec        curr_x, ref_x;
    real        gaussian;
    rvec        act_center;  /* Center of actual coordinates that are passed to the fit routine */
    rvec        ref_center;  /* Same for the reference coordinates */
    double      fitangle;    /* This will be the actual angle of the rotation group derived
                              * from an RMSD fit to the reference structure at t=0 */
    gmx_slabdata_t sd;
    rvec        coord;
    real        scal;


    /**********************************/
    /* First collect the data we need */
    /**********************************/

    /* Clear number of relevant atoms in all slabs */
    for (l=0; l<2*rotg->slab_max_nr+1; l++)
        rotg->slab_data[l].nat = 0;

    /* Loop over ALL rotation group atoms in all slabs */
    for(l=0; l<rotg->nat; l++)
    {
        for(n = -rotg->slab_first; n <= rotg->slab_last; n++)
        {
            islab = n+rotg->slab_max_nr; /* slab index */
            /* Current coordinate of this atom: x[ii][XX/YY/ZZ] */
            copy_rvec(rotg->xc[l], curr_x);
            /* Calculate the Gaussian value of curr_slab for curr_x */
            gaussian = gaussian_weight(curr_x, rotg, n);
            /* Only do the calculation for this slab if the gaussian is large enough: */
            if(gaussian > rotg->min_gaussian)
            {
                /* The (unrotated) reference coordinate of this atom is copied to ref_x: */
                copy_rvec(rotg->xc_ref[l], ref_x);
                /* Save data for doing angular RMSD fit later */
                sd = &(rotg->slab_data[islab]);
                ind = sd->nat;
                /* Save the current atom coordinate */
                copy_rvec(curr_x, sd->x[ind]);
                /* Save the corresponding reference coordinate */
                copy_rvec(ref_x , sd->ref[ind]);
                /* Save the weight for this atom in this slab */
                sd->weight[ind] = gaussian;
                /* Don't forget to keep track of the number of relevant atoms we found in this slab */
                sd->nat++;
            }
        }
    }

    /* Get the center of all rotation group coordinates: */
    get_center(rotg->xc, NULL, rotg->nat, act_center);


    /**************************/
    /* Now do the calculation */
    /**************************/

    /* === Determine the optimal fit angle for the whole rotation group === */

    /* METHOD 1 - coordinates stay as they are, 'normal' RMSD fit */
    /* Note that from the point of view of the current coordinates, the reference has rotated backwards,
     * but we want to output the angle relative to the fixed reference, therefore the minus sign. */
    fitangle = -opt_angle_analytic(rotg->xc_ref, rotg->xc, NULL, rotg->nat, rotg->xc_ref_center, act_center, rotg->vec);
    fprintf(fp, "%12.3e%12.3f%12.3lf", t, degangle, fitangle);

    /* METHOD 2 - now with normalization of every coordinate to it's reference coordinate length */
    for (i=0; i<rotg->nat; i++)
    {
        /* First put center of coordinates into origin */
        rvec_sub(rotg->xc[i], act_center, coord);
        /* Determine the scaling factor for the coordinate: */
        scal = rotg->xc_ref_length[i] / norm(coord);
        /* Get coordinate, multiply with the scaling factor and save in buf[i] */
        svmul(scal, coord, rotg->xc_norm[i]);
    }
    fitangle = -opt_angle_analytic(rotg->xc_ref, rotg->xc_norm, NULL, rotg->nat, rotg->xc_ref_center, act_center, rotg->vec);
    fprintf(fp_ref, "%12.3e%12.3f%12.3lf", t, degangle, fitangle);


    /* === Now do RMSD fitting with the 2 methods from above for each slab === */
    /* We require at least SLAB_MIN_ATOMS in a slab, such that the fit makes sense. */
#define SLAB_MIN_ATOMS 9

    /* METHOD 1 for each slab */
    for(n = -rotg->slab_first; n <= rotg->slab_last; n++)
    {
        islab = n+rotg->slab_max_nr; /* slab index */
        sd = &(rotg->slab_data[islab]);
        if (sd->nat >= SLAB_MIN_ATOMS)
        {
            /* Get the center of the slabs reference and current coordinates */
            get_center(sd->ref, sd->weight, sd->nat, ref_center);
            get_center(sd->x  , sd->weight, sd->nat, act_center);
            fitangle = -opt_angle_analytic(sd->ref, sd->x, sd->weight, sd->nat, ref_center, act_center, rotg->vec);
            fprintf(fp, "%6d%6d%12.3f", n, sd->nat, fitangle);
        }
    }

    /* METHOD 2 for each slab */
    for(n = -rotg->slab_first; n <= rotg->slab_last; n++)
    {
        islab = n+rotg->slab_max_nr; /* slab index */
        sd = &(rotg->slab_data[islab]);
        if (sd->nat >= SLAB_MIN_ATOMS)
        {
            /* Get the center of the slabs reference and current coordinates */
            get_center(sd->ref, sd->weight, sd->nat, ref_center);
            get_center(sd->x  , sd->weight, sd->nat, act_center);
            /* Center: */
            for (i=0; i<sd->nat;i++)
            {
                rvec_dec(sd->ref[i], ref_center);
                rvec_dec(sd->x[i]  , act_center);
                /* Normalize x_i such that it gets the same length as ref_i */
                svmul( norm(sd->ref[i])/norm(sd->x[i]), sd->x[i], sd->x[i] );
            }
            /* We already subtracted the centers */
            clear_rvec(ref_center);
            clear_rvec(act_center);
            fitangle = -opt_angle_analytic(sd->ref, sd->x, sd->weight, sd->nat, ref_center, act_center, rotg->vec);
            fprintf(fp_ref, "%6d%6d%12.3f", n, sd->nat, fitangle);
        }
    }

    fprintf(fp     , "\n");
    fprintf(fp_ref , "\n");

#undef SLAB_MIN_ATOMS
}


/* Get the shifts such that each atom is within closest
 * distance to its position at the last NS time step after shifting.
 * If we start with a whole structure, and always keep track of
 * shift changes, the structure will stay whole this way */
static void get_extra_shifts(
        t_rotgrp  *rotg,
        int       npbcdim,
        matrix    box,
        rvec      *xc)
{
    int  i,m,d;
    rvec dx;


    for (i=0; i < rotg->nat; i++)
        clear_ivec(rotg->xc_eshifts[i]);

    for (i=0; i<rotg->nat; i++)
    {
        /* The distance this atom moved since the last time step */
        /* If this is more than just a bit, it has changed its home pbc box */
        rvec_sub(xc[i],rotg->xc_old[i],dx);

        for(m=npbcdim-1; m>=0; m--)
        {
            while (dx[m] < -0.5*box[m][m])
            {
                for(d=0; d<DIM; d++)
                    dx[d] += box[m][d];
                rotg->xc_eshifts[i][m]++;
            }
            while (dx[m] >= 0.5*box[m][m])
            {
                for(d=0; d<DIM; d++)
                    dx[d] -= box[m][d];
                rotg->xc_eshifts[i][m]--;
            }
        }
    }
}


/* Assemble the coordinates such that every node has all of them */
static void get_coordinates(
        t_commrec *cr,
        t_rotgrp *rotg,
        rvec      x[],
        bool      bNeedShiftsUpdate, /* NS step, the shifts have changed */
        matrix    box)
{
    int i, ii, j, l;


    GMX_MPE_LOG(ev_get_coords_start);

    /* Zero out the collective coordinate array */
    clear_rvecs(rotg->nat, rotg->xc);

    /* Put the local coordinates that this node has into the right place of
     * the collective array. */
    for (l=0; l<rotg->nat_loc; l++)
    {
        /* Local index of a rotation group atom  */
        ii = rotg->ind_loc[l];
        /* Index of this rotation group atom w.r.t. the whole rotation group */
        j = rotg->xc_ref_ind[l];
        /* Sort the current x-coordinates into the right place of the array: */
        copy_rvec(x[ii], rotg->xc[j]);
    }
    if (PAR(cr))
    {
        /* Add the arrays from all nodes together */
        gmx_sum(DIM*rotg->nat, rotg->xc[0], cr);

        /* To make the molecule whole, start with a whole structure and each
         * step move the assembled coordinates at closest distance to the positions
         * from the last step. First shift the coordinates with the saved shift
         * vectors (these are 0 when this routine is called for the first time!) */
        ed_shift_coords(box, rotg->xc, rotg->xc_shifts, rotg->nat);

        /* Now check if some shifts changed since the last step.
         * This only needs to be done when the shifts are expected to have changed,
         * i.e. after neighboursearching */
        if (bNeedShiftsUpdate)
        {
            get_extra_shifts(rotg, 3, box, rotg->xc);

            /* Shift with the additional shifts such that we get a whole molecule now */
            ed_shift_coords(box, rotg->xc, rotg->xc_eshifts, rotg->nat);

            /* Add the shift vectors together for the next time step */
            for (i=0; i<rotg->nat; i++)
            {
                rotg->xc_shifts[i][XX] += rotg->xc_eshifts[i][XX];
                rotg->xc_shifts[i][YY] += rotg->xc_eshifts[i][YY];
                rotg->xc_shifts[i][ZZ] += rotg->xc_eshifts[i][ZZ];
            }

            /* Store current correctly-shifted coordinates for comparison in the next NS time step */
            for (i=0; i<rotg->nat; i++)
                copy_rvec(rotg->xc[i],rotg->xc_old[i]);
        }
    }

    GMX_MPE_LOG(ev_get_coords_finish);
}


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


static real do_flex2_lowlevel(
        t_rotgrp  *rotg,
        matrix    rotmat,
        real      sigma,    /* The Gaussian width sigma */
        rvec      x[],
        bool      bCalcTorque,
        matrix    box,
        t_commrec *cr)
{
    int  i,ii,l,m,n,islab,ipgrp;
    rvec dr;                /* difference vector between actual and reference coordinate */
    real V;                 /* The rotation potential energy */
    real gaussian;          /* Gaussian weight */
    real gaussian_xi;
    real beta;
    rvec curr_x;            /* particle coordinate */
    rvec xi;
    rvec curr_x_rel;        /* particle coordinate relative to COG */
    rvec curr_COG;          /* the current slab's center of geometry (COG) */
    rvec ref_x, ref_x_cpy;  /* the reference particle coordinate */
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

    
    /********************************************************/
    /* Main loop over all local atoms of the rotation group */
    /********************************************************/
    V = 0.0;
    for (l=0; l<rotg->nat_loc; l++)
    {
        /* Local index of a rotation group atom  */
        ii = rotg->ind_loc[l];
        /* Index of this rotation group atom with respect to the whole rotation group */
        ipgrp = rotg->xc_ref_ind[l];
        /* For each atom, loop over all slabs. We could have contributions from any slab */
        clear_rvec(sum_f_ii);

        for (n = -rotg->slab_first; n <= rotg->slab_last; n++)
        {
            clear_rvec(force_n);

            islab = n+rotg->slab_max_nr; /* slab index */
            /* Current coordinate of this atom: x[ii][XX/YY/ZZ] */
            copy_rvec(x[ii], curr_x);

            /* Shift this atom such that it is near its reference */
            shift_single_coord(box, curr_x, rotg->xc_shifts[ipgrp]);

            /* Calculate the Gaussian value of curr_slab for curr_x */
            gaussian = gaussian_weight(curr_x, rotg, n);

            /* Only do the calculation for this slab if the Gaussian is large enough: */
            if (gaussian > rotg->min_gaussian)
            {
                /* The (unrotated) reference coordinate of this atom is copied to ref_x: */
                copy_rvec(rotg->xc_ref[ipgrp], ref_x);
                beta = calc_beta(curr_x, rotg,n);
                /* The center of geometry (COG) of this slab is copied to curr_COG: */
                copy_rvec(rotg->slab_center[islab], curr_COG);
                /* The reference COG of this slab is copied to ref_COG: */
                copy_rvec(rotg->slab_center_ref[islab], ref_COG);

                /* Rotate the reference coordinate around the rotation axis through this slab's reference COG */
                /* 1. Subtract the reference slab COG from the reference coordinate i */
                rvec_sub(ref_x, ref_COG, ref_x); /* ref_x =y_ii-y_0 */
                /* 2. Rotate reference coordinate around origin: */
                copy_rvec(ref_x, ref_x_cpy);
                mvmul(rotmat, ref_x_cpy, ref_x); /* ref_x = r_ii = Omega.(y_ii-y_0) */

                /* Now subtract the slab COG from current coordinate i */
                rvec_sub(curr_x, curr_COG, curr_x_rel); /* curr_x_rel = x_ii-x_0 */

                /* Force on atom i is based on difference vector between current coordinate and rotated reference coordinate */
                /* Difference vector between current and reference coordinate: */
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

                for(i=0; i<rotg->nat; i++)
                {
                    /* Coordinate xi of this atom */
                    copy_rvec(rotg->xc[i],xi);

                    gaussian_xi = gaussian_weight(xi,rotg,n); /* g_n(xi)*/
                    if (gaussian_xi > rotg->min_gaussian)
                    {
                        /* Calculate r_i for atom i and slab n: */
                        /* Unrotated reference coordinate y_i: */
                        copy_rvec(rotg->xc_ref[i],yi);
                        
                        /* COG y0 for this slab: */
                        /* The reference COG of this slab is still in ref_COG */
                        rvec_sub(yi, ref_COG, tmp);   /* tmp = y_i - y_0                       */
                        mvmul(rotmat, tmp, r);        /* r   = Omega*(y_i - y_0)               */
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
                    }
                } /* now we have the sum-over-atoms (over i) for the ii-th atom in the n-th slab */

                gauss_ratio=gaussian/rotg->slab_weights[islab];
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
                    rotg->slab_torque_v[islab] += torque(rotg->vec, force_n, curr_x, curr_COG);
                }
            } /* END of "if gaussian > pg->min_gaussian" */
        } /* END of loop over slabs */

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        for (m=0; m<DIM; m++)
            rotg->f_rot_loc[l][m] = sum_f_ii[m];

#ifdef INFOF
        fprintf(stderr," FORCE on ATOM %d/%d  = (%15.8f %15.8f %15.8f)  \n",
                l,ii,rotg->sum_f_ii[XX], rotg->sum_f_ii[YY], rotg->sum_f_ii[ZZ]5);
#endif
    } /* END of loop over local atoms */

    return V;
}


static real do_flex_lowlevel(
        t_rotgrp *rotg,
        matrix    rotmat,
        real      sigma,     /* The Gaussian width sigma */
        rvec      x[],
        bool      bCalcTorque,
        matrix    box,
        t_commrec *cr)
{
    int   i,ii,iii,l,m,n,islab,ipgrp;
    rvec  direction;         /* the direction for the force on atom i */
    rvec  sum_n1,sum_n2;     /* Two contributions to the rotation force */
    rvec  sum_i;             /* Inner part of sum_n2 */
    rvec  dummy1,tmp_f,s,tmp2;
    real  u;                 /* direction*dr */
    rvec  dr;                /* difference vector between actual and reference coordinate */
    real  V;                 /* The rotation potential energy */
    real  gaussian;          /* Gaussian weight */
    real  gaussian_xi;
    real  fac,beta;
    rvec  curr_x;            /* particle coordinate */
    rvec  xi;
    rvec  curr_x_rel;        /* particle coordinate relative to COG */
    rvec  curr_COG;          /* the current slab's center of geometry (COG) */
    rvec  ref_x, ref_x_cpy;  /* the reference particle coordinate */
    rvec  ref_COG;           /* the reference slab's COG */
    rvec  r, yi, tmp;
    rvec  force_n;           /* Single force from slab n on one atom */


    /********************************************************/
    /* Main loop over all local atoms of the rotation group */
    /********************************************************/
    V = 0.0;
    for (l=0; l<rotg->nat_loc; l++)
    {
        /* Local index of a rotation group atom  */
        ii = rotg->ind_loc[l];
        /* Index of this rotation group atom with respect to the whole rotation group */
        ipgrp = rotg->xc_ref_ind[l];
        /* For each atom, loop over all slabs. We could have contributions from any slab */
        clear_rvec(sum_n1);
        clear_rvec(sum_n2);
        for (n = -rotg->slab_first; n <= rotg->slab_last; n++)
        {
            islab = n+rotg->slab_max_nr; /* slab index */
            /* Current coordinate of this atom: x[ii][XX/YY/ZZ] */
            copy_rvec(x[ii], curr_x);

            /* Shift this atom such that it is near its reference */
            shift_single_coord(box, curr_x, rotg->xc_shifts[ipgrp]);

            /* Calculate the Gaussian value of curr_slab for curr_x */
            gaussian = gaussian_weight(curr_x, rotg, n);

            /* Only do the calculation for this slab if the Gaussian is large enough: */
            if (gaussian > rotg->min_gaussian)
            {
                /* The (unrotated) reference coordinate of this atom is copied to ref_x: */
                copy_rvec(rotg->xc_ref[ipgrp], ref_x);
                beta = calc_beta(curr_x, rotg,n);
                /* The center of geometry (COG) of this slab is copied to curr_COG: */
                copy_rvec(rotg->slab_center[islab], curr_COG);
                /* The reference COG of this slab is copied to ref_COG: */
                copy_rvec(rotg->slab_center_ref[islab], ref_COG);

                /* Rotate the reference coordinate around the rotation axis through this slab's reference COG */
                /* 1. Subtract the reference slab COG from the reference coordinate i */
                rvec_sub(ref_x, ref_COG, ref_x); /* ref_x =y_ii-y_0 */
                /* 2. Rotate reference coordinate around origin: */
                copy_rvec(ref_x, ref_x_cpy);
                mvmul(rotmat, ref_x_cpy, ref_x); /* ref_x = r_ii = Omega.(y_ii-y_0) */

                /* Now subtract the slab COG from current coordinate i */
                rvec_sub(curr_x, curr_COG, curr_x_rel); /* curr_x_rel = x_ii-x_0 */

                /* Force on atom i is based on difference vector between current coordinate and rotated reference coordinate */
                /* Difference vector between current and reference coordinate: */
                rvec_sub(curr_x_rel, ref_x, dr); /* dr=(x_ii-x_0)-Omega.(y_ii-y_0) */

                /* Calculate the direction of the actual force */
                cprod(rotg->vec, ref_x, direction);
                unitv(direction,direction);
                u = iprod(direction,dr);

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
                for (i=0; i<rotg->nat; i++)
                {
                    /* Index of a rotation group atom  */
                    iii = rotg->ind[i];

                    /* Coordinate xi of this atom */
                    copy_rvec(rotg->xc[i],xi);

                    gaussian_xi = gaussian_weight(xi,rotg,n);
                    
                    if (gaussian_xi > rotg->min_gaussian)
                    {
                        /* Calculate r_i for atom i and slab n: */
                        /* Unrotated reference coordinate y_i: */
                        copy_rvec(rotg->xc_ref[i],yi);
                        
                        /* COG y0 for this slab: */
                        /* The reference COG of this slab is still in ref_COG */
                        rvec_sub(yi, ref_COG, tmp);   /* tmp = y_i - y_0             */
                        mvmul(rotmat, tmp, r);        /* r   = Omega*(y_i - y_0)     */
                        cprod(rotg->vec, r, tmp);     /* tmp = a x Omega*(y_i - y_0) */
                        unitv(tmp, s);                /* s   = a x Omega*(y_i - y_0) / |a x Omega*(y_i - y_0)| */
                        /* Now that we have ri and si, let's calculate the i-sum: */
                        /* tmp = x_0 - x_l */
                        rvec_sub(curr_COG, curr_x, tmp);
                        /* tmp2 = beta/sigma^2 * (s*(x_0 - x_l)) * a   */
                        svmul(beta*iprod(tmp, s)/sqr(sigma), rotg->vec, tmp2);
                        /* tmp = s + tmp2 */
                        rvec_add(tmp2, s, tmp);
                        /* tmp2 = xi - x0 */
                        rvec_sub(xi, curr_COG, tmp2);
                        /* tmp2 = xi - x0 - ri */
                        rvec_sub(tmp2, r, tmp2);
                        /* fac = gn * s*(xi - x0 - ri) */
                        fac = gaussian_xi*iprod(s, tmp2);
                        /* tmp2 = gn * s*(xi - x0 - ri) * [beta/sigma^2 * (s*(x_0 - x_l)) * a] */
                        svmul(fac, tmp, tmp2);
                        rvec_add(sum_i, tmp2, sum_i);
                    }
                } /* now we have the i-sum for this atom in this slab */
                svmul(gaussian/rotg->slab_weights[islab], sum_i, sum_i);
                rvec_add(sum_n2, sum_i, sum_n2);

                /* Calculate the torque: */
                if (bCalcTorque)
                {
                    /* The force on atom ii from slab n only: */
                    rvec_sub(sum_i, dummy1, force_n);
                    svmul(rotg->k, force_n, force_n);
                    rotg->slab_torque_v[islab] += torque(rotg->vec, force_n, curr_x, curr_COG);
                }
            } /* END of "if gaussian > pg->min_gaussian" */
        } /* END of loop over slabs */

        /* Put both contributions together: */
        rvec_sub(sum_n2,sum_n1,tmp_f); /* F = -grad V */

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        for(m=0; m<DIM; m++)
            rotg->f_rot_loc[l][m] = rotg->k*tmp_f[m];


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
                rotg->xc[i][XX], rotg->xc[i][YY], rotg->xc[i][ZZ]);
    }
    fflush(fp);

}
#endif


/* Enforced rotation with a flexible axis */
static void do_flexible(
        t_commrec *cr,
        t_rotgrp  *rotg,      /* The rotation group             */
        int       g,          /* Group number                   */
        real      degangle,   /* Angle theta_ref [degrees]      */
        matrix    rotmat,     /* Rotation matrix for this angle */
        rvec      x[],        /* The local coordinates          */
        matrix    box,
        double    t,          /* Time in picoseconds            */
        bool      bDynBox,    /* Is the box dynamic?            */
        int       step,       /* The time step                  */
        bool      bOutstep,
        FILE      *fp_slabs,
        FILE      *fp_torque,
        FILE      *fp_angles,
        FILE      *fp_nangles)
{
    int          l;
    real         sigma;             /* The Gaussian width sigma */


    /* Define the sigma value */
    sigma = 0.7*rotg->slab_dist;
    
    /* Determine the gaussian-weighted center of coordinates for all slabs,
     * also reallocate memory if the number of slabs grows (i.e. box expands) */
    get_slab_centers(rotg,rotg->xc,box,cr,g,bDynBox,t,fp_slabs,bOutstep,FALSE);
        
    /* Clear the torque per slab from last time step: */
    for (l=0; l<2*rotg->slab_max_nr+1; l++)
        rotg->slab_torque_v[l] = 0.0;
    
    /* Call the rotational forces kernel */
    GMX_MPE_LOG(ev_flexll_start);
    if (rotg->eType == erotgFLEX)
        rotg->V = do_flex_lowlevel(rotg, rotmat, sigma, x, bOutstep, box, cr);
    else if (rotg->eType == erotgFLEX2)
        rotg->V = do_flex2_lowlevel(rotg, rotmat, sigma, x, bOutstep, box, cr);
    else
        gmx_fatal(FARGS, "Unknown flexible rotation type");
    GMX_MPE_LOG(ev_flexll_finish);
    
    /* Determine actual angle of this slab by RMSD fit and output to file - Let's hope */
    /* this only happens once in a while, since this is not parallelized! */
    if (bOutstep && MASTER(cr))
        flex_fit_angle(rotg, t, degangle, fp_angles, fp_nangles);
}


/* Calculate the angle between reference and actual rotation group atom: */
/* This routine assumes that rot_vec is a unit vector! */
static void angle(t_rotgrp *rotg,
        rvec x_act,
        rvec x_ref,
        real *alpha,
        real *weight)  /* atoms near the rotation axis should count less than atoms far away */
{
    rvec x , xr ;  /* actual and reference coordinate in new coordinate system */
    rvec xp, xrp;  /* dito, but projected on a plane perpendicular to pg->vec */
    real normxp;   /* saved for efficiency reasons */
    rvec dum;
    real cosalpha; /* cos of angle between projected reference and actual coordinate */
    int sign;


    /* Move the center of coordinates to rot_offset: */
    rvec_sub(x_act, rotg->offset, x);
    rvec_sub(x_ref, rotg->offset, xr);

    /* Project xr and x into a plane through the origin perpendicular to rot_vec: */
    /* Project xr: xrp = xr - (rot_vec * xr) * rot_vec */
    svmul(iprod(rotg->vec, xr), rotg->vec, dum);    /* only works if pg->vec is a unit vector! */
    rvec_sub(xr, dum, xrp);
    /* Project x: */
    svmul(iprod(rotg->vec, x), rotg->vec, dum);    /* only works if pg->vec is a unit vector! */
    rvec_sub(x, dum, xp);

    /* Calculate the angle between the projected coordinates: */
    normxp = norm(xp); /* save for later use */
    cosalpha = iprod(xrp, xp) / (norm(xrp)*normxp);
    if (cosalpha < -1.0) cosalpha = -1.0;
    if (cosalpha >  1.0) cosalpha =  1.0;

    /* Retrieve some information about which vector precedes */
    cprod(xp, xrp, dum); /* if reference precedes, this is pointing into the same direction as pg->vec */

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
}


/* Fixed rotation: The rotation reference group rotates around an axis */
/* The atoms of the actual rotation group are attached with imaginary  */
/* springs to the reference atoms.                                     */
static void do_fixed(
        t_commrec *cr,
        t_rotgrp *rotg,         /* The rotation group       */
        matrix    rotmat,       /* rotary matrix            */
        rvec     x[],           /* The coordinates (natoms) */
        t_pbc    *pbc,
        double   t,             /* Time in picoseconds      */
        int      step,          /* The time step            */
        bool     bTorque)
{
    int       i,ii,m,iigrp;
    rvec      dr;
    rvec      x1;              /* particle coordinate */
    rvec      tmp_f;           /* Force */
    rvec      xr, xrcpy;       /* rotated (reference) particle coordinate */
    real      alpha;           /* a single angle between an actual and a reference coordinate */
    real      weight;          /* single weight for a single angle */

    
    /* Clear values from last time step */
    rotg->V            = 0.0;
    rotg->fix_torque_v = 0.0;
    rotg->fix_angles_v = 0.0;
    rotg->fix_weight_v = 0.0;
    
    /* Loop over all local atoms of the rotation group */
    for (i=0; i<rotg->nat_loc; i++)
    {
        /* Index of a rotation group atom  */
        ii = rotg->ind_loc[i];
        /* Actual coordinate of this atom: x[ii][XX/YY/ZZ] */
        copy_rvec(x[ii], x1);
        
        /* Index of this rotation group atom with respect to the whole rotation group */
        iigrp = rotg->xc_ref_ind[i];
        
        /* Copy the (unrotated) reference coordinate of this atom: */
        copy_rvec(rotg->xc_ref[iigrp], xr);
        /* Rotate this atom around dislocated rotation axis: */
        /* Move rotation axis, so that it runs through the origin: */
        rvec_sub(xr, rotg->offset, xr);
        /* Rotate around the origin: */
        copy_rvec(xr, xrcpy);
        mvmul(rotmat, xrcpy, xr);
        /* And move back: */
        rvec_add(xr, rotg->offset, xr);
        /* Difference vector between reference and actual coordinate: */
        pbc_dx(pbc,xr,x1, dr);
        
        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        for (m=0; m<DIM; m++)
        {
            tmp_f[m]              = rotg->k*dr[m];
            rotg->f_rot_loc[i][m] = tmp_f[m];
            rotg->V              +=  0.5*rotg->k*sqr(dr[m]);
        }
        
        if (bTorque)
        {
            /* Add to the torque of this rotation group */
            rotg->fix_torque_v += torque(rotg->vec, tmp_f, x1, rotg->offset);
            
            /* Calculate the angle between reference and actual rotation group atom: */
            angle(rotg, x1, xr, &alpha, &weight);  /* angle in rad, weighted */
            rotg->fix_angles_v += alpha * weight;
            rotg->fix_weight_v += weight;
            /* Use the next two lines instead if you don't want weighting: */
            /*
                angles_v[g] += alpha;
                weight_v[g] += 1;
             */
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
    } /* end of loop over local rotation group atoms */
}


extern void init_rot_group(FILE *fplog,t_commrec *cr,
        int g,t_rotgrp *rotg,
        rvec *x,       /* the coordinates */
        int rot_type,
        FILE *out_slabs,
        matrix box)
{
    int i,ii;
    real        box_d;        /* The box diagonal (needed for maximum slab count) */
    char        filename[255];/* File to save the reference coordinates in for enforced rotation */
    rvec        f_box[3];     /* Box from reference file */
    rvec        coord;
    t_trnheader header;       /* Header information of reference file */
    bool        bSame;        /* Used for a check if box sizes agree */
    int         nslabs;       /* Maximum number of slabs that fit into simulation box */
    bool        bFlex;

    /* Enforced rotation with fixed/flexible axis */    
    bFlex = (rot_type == erotgFLEX || rot_type == erotgFLEX2);
    
    snew(rotg->xc_ref       , rotg->nat);
    snew(rotg->xc_ref_length, rotg->nat);
    snew(rotg->xc_ref_ind   , rotg->nat);
    snew(rotg->f_rot_loc    , rotg->nat);

    /* xc_ref_ind needs to be set to identity in the serial case */
    if (!PAR(cr))
        for (i=0; i<rotg->nat; i++)
            rotg->xc_ref_ind[i] = i;

    /* Enforced rotation with flexible axis */
    if (bFlex)
    {
        snew(rotg->xc        , rotg->nat);
        snew(rotg->xc_norm   , rotg->nat);
        snew(rotg->xc_old    , rotg->nat);
        snew(rotg->xc_shifts , rotg->nat);
        snew(rotg->xc_eshifts, rotg->nat);

        /* A maximum of (box diagonal)/(slab distance) slabs are possible */
        box_d = diagonal_length(box);
        rotg->slab_max_nr = (int) ceil(box_d/rotg->slab_dist);
        nslabs = 2*rotg->slab_max_nr + 1;
        if (MASTER(cr))
            fprintf(stdout, "Enforced rotation: allocating memory to store data for %d slabs (rotation group %d).\n",nslabs,g);
        snew(rotg->slab_center    , nslabs);
        snew(rotg->slab_center_ref, nslabs);
        snew(rotg->slab_weights   , nslabs);
        snew(rotg->slab_torque_v  , nslabs);
        snew(rotg->slab_data      , nslabs);
        for (i=0; i<nslabs; i++)
        {
            snew(rotg->slab_data[i].x     , rotg->nat);
            snew(rotg->slab_data[i].ref   , rotg->nat);
            snew(rotg->slab_data[i].weight, rotg->nat);
        }
    }

    /* Read in rotation reference coordinates from file, if it exists.
     * If not, write out the initial rotation group coordinates as reference coordinates */
    if (MASTER(cr))
    {
        /* Save the reference coordinates to trr */
        /* Make a trr for each rotation group */
        sprintf(filename, "ref_%d_%s.trr", g, erotg_names[rotg->eType]);
        if (gmx_fexist(filename)) /* Read rotation reference coordinates from file */
        {
            fprintf(fplog, "Enforced rotation: found reference coordinate file %s.\n", filename);
            read_trnheader(filename, &header);
            if (rotg->nat != header.natoms)
                gmx_fatal(FARGS,"Number of atoms in coordinate file %s (%d) does not match the number of atoms in rotation group (%d)!\n",
                        filename, header.natoms, rotg->nat);
            read_trn(filename, &header.step, &header.t, &header.lambda, f_box, &header.natoms, rotg->xc_ref, NULL, NULL);
            fprintf(fplog , "Enforced rotation: read reference coordinates for group %d from %s.\n", g, filename);
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
        else /* Save the initial coordinates of the rotation group as reference */
        {
            for(i=0; i<rotg->nat; i++)
            {
                ii = rotg->ind[i];
                copy_rvec(x[ii], rotg->xc_ref[i]);
            }
            write_trn(filename,g,0.0,0.0,box,rotg->nat,rotg->xc_ref,NULL,NULL);
            fprintf(fplog, "Enforced rotation: saved %d coordinates of group %d to %s.\n",
                    rotg->nat, g, filename);
        }
        if (bFlex)
        {
            /* Save the original (whole) coordinates such that in the next time step
             * the molecule can be made whole again (in the parallel case) */
            for (i=0; i<rotg->nat; i++)
            {
                ii = rotg->ind[i];
                copy_rvec(x[ii], rotg->xc_old[i]);
            }
        }
    }
#ifdef GMX_MPI
    /* Copy reference coordinates to all PP nodes */
    if (PAR(cr))
    {
        gmx_bcast(rotg->nat*sizeof(rotg->xc_ref[0]), rotg->xc_ref, cr);
        if (bFlex)
            gmx_bcast(rotg->nat*sizeof(rotg->xc_old[0]),rotg->xc_old, cr);
    }
#endif

    if (bFlex)
    {
        /* Flexible rotation: determine the reference COGs for the rest of the simulation */
        get_slab_centers(rotg,rotg->xc_ref,box,cr,g,TRUE,-1,out_slabs,1,TRUE);

        /* Also save the center of coordinates for the reference structure: */
        get_center(rotg->xc_ref, NULL, rotg->nat, rotg->xc_ref_center);
        if (MASTER(cr))
            fprintf(fplog, "Enforced rotation: center of reference coordinates of group %d is at %f %f %f", g,
                    rotg->xc_ref_center[XX], rotg->xc_ref_center[YY], rotg->xc_ref_center[ZZ]);

        /* Length of each x_rotref vector from center (needed for later RMSD fit routines): */
        for (i=0; i<rotg->nat; i++)
        {
            rvec_sub(rotg->xc_ref[i], rotg->xc_ref_center, coord);
            rotg->xc_ref_length[i] = norm(coord);
        }
    }
}


static void make_local_rotation_group(gmx_ga2la_t ga2la,
        t_rotgrp *rotg,int start,int end)
{
    int i,ii;


    rotg->nat_loc = 0;
    for(i=0; i<rotg->nat; i++) {
        ii = rotg->ind[i];
        if (!ga2la_home(ga2la,ii,&ii)) {
            ii = -1;
        }

        if (ii >= start && ii < end) {
            /* This is a home atom, add it to the local rotation group */
            if (rotg->nat_loc >= rotg->nalloc_loc) {
                rotg->nalloc_loc = over_alloc_dd(rotg->nat_loc+1);
                srenew(rotg->ind_loc,rotg->nalloc_loc);
            }
            rotg->ind_loc[rotg->nat_loc] = ii;
            /* Copy the reference coordinates */
            if (rotg->xc_ref)
            {
                /* Remember which of the x_rotref coordinates are local: */
                rotg->xc_ref_ind[rotg->nat_loc]=i;  /* i is the number of the atom with respect to the whole rotation group */
                /* pg->ind[i] would be the number with respect to the whole system! */
            }
            rotg->nat_loc++;
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

    /* Indicate that the shift vectors for this structure need to be updated
     * at the next call to get_coordinates, since obviously we are in a NS step */
    rot->bUpdateShifts = TRUE;

}


void init_rot(FILE *fplog,t_inputrec *ir,
        t_commrec *cr, matrix box, rvec *x, unsigned long Flags)
{
    t_rot    *rot;
    t_rotgrp *rotg;
    int      g;
    bool     bRerun;
    bool     bFlex=FALSE;


    rot = ir->rot;

    /* Output every step for reruns */
    bRerun = Flags & MD_RERUN;
    if (bRerun)
    {
        if (fplog)
            fprintf(fplog, "Enforced rotation: rerun - will write rotation output every available step.\n");
        rot->nstrout = 1;
        rot->nsttout = 1;
    }

    rot->out_slabs = NULL;
    if (MASTER(cr))
        rot->out_slabs = open_slab_out(rot);

    for(g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        
        if (fplog)
            fprintf(fplog,"Enforced rotation: group %d type '%s'\n", g, erotg_names[rotg->eType]);

        if (rotg->eType == erotgFLEX || rotg->eType == erotgFLEX2)
            bFlex = TRUE;
        
        if (rotg->nat > 0)
        {
            if (PAR(cr))
            {
                rotg->nat_loc    = 0;
                rotg->nalloc_loc = 0;
                rotg->ind_loc    = NULL;
            }
            else
            {
                rotg->nat_loc = rotg->nat;
                rotg->ind_loc = rotg->ind;
            }
            init_rot_group(fplog,cr,g,rotg,x,rotg->eType,rot->out_slabs,box);
        }
    }
    
    /* Buffers for MPI reducing torques, angles, weights (for each group), and V */
    rot->bufsize = 4*rot->ngrp; /* To start with */
    snew(rot->inbuf , rot->bufsize);
    snew(rot->outbuf, rot->bufsize);

    /* Only do I/O on the MASTER */
    rot->out_angles  = NULL;
    rot->out_nangles = NULL;
    rot->out_rot     = NULL;
    rot->out_torque  = NULL;
    if (MASTER(cr))
    {
        rot->out_rot = open_rot_out(rot);
        if (bFlex)
        {
            if (rot->nstrout > 0)
            {
                rot->out_angles  = open_angles_out(rot, "angles.log");
                rot->out_nangles = open_angles_out(rot, "angles_n.log");
            }
            if (rot->nsttout > 0)
                rot->out_torque  = open_torque_out(rot);
        }
    }
}


extern void do_rotation(
        t_commrec *cr,
        t_inputrec *ir,
        matrix box,
        rvec x[],
        real t,
        int step,
        gmx_wallcycle_t wcycle)
{
    int      g;
    t_pbc    pbc;
    t_rot    *rot;
    t_rotgrp *rotg;
    bool     outstep_torque;
    float    cycles_rot;
    real     degangle;
    matrix   rotmat;

    
    if ( (PAR(cr)) && !DOMAINDECOMP(cr) )
        gmx_fatal(FARGS, "Enforced rotation is only implemented for domain decomposition!");

    rot=ir->rot;
    
    /* At which time steps do we want to output the torque */
    outstep_torque = do_per_step(step, rot->nsttout);

    /* Output time into rotation output file */
    if (outstep_torque && MASTER(cr))
        fprintf(rot->out_rot, "%12.3e",t);
    
    /* Let's first do ALL the communication! */
    for(g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        /* Transfer the rotation group's coordinates such that every node has all of them.
         * Every node contributes its local coordinates x and stores it in
         * the collective pg->xc array. */
        if (rotg->eType != erotgFIXED)
            get_coordinates(cr, rotg, x, rot->bUpdateShifts, box);
    }
    
    /* Done communicating, we can start to count cycles now ... */
    wallcycle_start(wcycle, ewcROT);
    GMX_MPE_LOG(ev_rotcycles_start);
    
    for(g=0; g<rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        degangle = rotg->rate * t; /* angle of rotation for this group: */
        if (outstep_torque && MASTER(cr))
            fprintf(rot->out_rot, "%12.4f", degangle);
        /* Calculate the rotation matrix for this angle: */
        calc_rotmat(rotg->vec,degangle,rotmat);

        if (rotg->eType == erotgFIXED)
        {
            set_pbc(&pbc,ir->ePBC,box);
            do_fixed(cr,rotg,rotmat,x,&pbc,t,step,outstep_torque);
        }
        else
        {
            do_flexible(cr,rotg,g,degangle,rotmat,x,box,t,DYNAMIC_BOX(*ir),step,outstep_torque,
                    rot->out_slabs,rot->out_torque,rot->out_angles,rot->out_nangles);
        }
    }
    /* If bUpdateShifts was TRUE then the shifts have just been updated in get_coordinates.
     * We do not need to update the shifts until the next NS step */
    rot->bUpdateShifts = FALSE;

    /* Stop the cycle counter and add to the force cycles for load balancing */
    cycles_rot = wallcycle_stop(wcycle,ewcROT);
    if (DOMAINDECOMP(cr) && wcycle)
        dd_cycles_add(cr->dd,cycles_rot,ddCyclF);
    GMX_MPE_LOG(ev_rotcycles_finish);
}
