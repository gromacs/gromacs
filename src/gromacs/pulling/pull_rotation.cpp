/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "gromacs/pulling/pull_rotation.h"

#include "config.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/nrjac.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

struct gmx_output_env_t;

static const std::string RotStr = { "Enforced rotation:" };

/* Set the minimum weight for the determination of the slab centers */
#define WEIGHT_MIN (10 * GMX_FLOAT_MIN)

//! Helper structure for sorting positions along rotation vector
struct sort_along_vec_t
{
    //! Projection of xc on the rotation vector
    real xcproj;
    //! Index of xc
    int ind;
    //! Mass
    real m;
    //! Position
    rvec x;
    //! Reference position
    rvec x_ref;
};


//! Enforced rotation / flexible: determine the angle of each slab
struct gmx_slabdata
{
    //! Number of atoms belonging to this slab
    int nat;
    /*! \brief The positions belonging to this slab.
     *
     * In general, this should be all positions of the whole
     * rotation group, but we leave those away that have a small
     * enough weight. */
    rvec* x;
    //! Same for reference
    rvec* ref;
    //! The weight for each atom
    real* weight;
};


//! Helper structure for potential fitting
struct gmx_potfit
{
    /*! \brief Set of angles for which the potential is calculated.
     *
     * The optimum fit is determined as the angle for with the
     * potential is minimal. */
    real* degangle;
    //! Potential for the different angles
    real* V;
    //! Rotation matrix corresponding to the angles
    matrix* rotmat;
};


//! Enforced rotation data for a single rotation group
struct gmx_enfrotgrp
{
    //! Input parameters for this group
    const t_rotgrp* rotg = nullptr;
    //! Index of this group within the set of groups
    int groupIndex;
    //! Rotation angle in degrees
    real degangle;
    //! Rotation matrix
    matrix rotmat;
    //! The atoms subject to enforced rotation
    std::unique_ptr<gmx::LocalAtomSet> atomSet;

    //! The normalized rotation vector
    rvec vec;
    //! Rotation potential for this rotation group
    real V;
    //! Array to store the forces on the local atoms resulting from enforced rotation potential
    rvec* f_rot_loc;

    /* Collective coordinates for the whole rotation group */
    //! Length of each x_rotref vector after x_rotref has been put into origin
    real* xc_ref_length;
    //! Center of the rotation group positions, may be mass weighted
    rvec xc_center;
    //! Center of the rotation group reference positions
    rvec xc_ref_center;
    //! Rotation group reference positions, perhaps centered
    std::vector<gmx::RVec> referencePositions;
    //! Current (collective) positions
    rvec* xc;
    //! Current (collective) shifts
    ivec* xc_shifts;
    //! Extra shifts since last DD step
    ivec* xc_eshifts;
    //! Old (collective) positions
    rvec* xc_old;
    //! Normalized form of the current positions
    rvec* xc_norm;
    //! Reference positions (sorted in the same order as xc when sorted)
    rvec* xc_ref_sorted;
    //! Where is a position found after sorting?
    int* xc_sortind;
    //! Collective masses
    real* mc;
    //! Collective masses sorted
    real* mc_sorted;
    //! one over the total mass of the rotation group
    real invmass;

    //! Torque in the direction of rotation vector
    real torque_v;
    //! Actual angle of the whole rotation group
    real angle_v;
    /* Fixed rotation only */
    //! Weights for angle determination
    real weight_v;
    //! Local reference coords, correctly rotated
    rvec* xr_loc;
    //! Local current coords, correct PBC image
    rvec* x_loc_pbc;
    //! Masses of the current local atoms
    real* m_loc;

    /* Flexible rotation only */
    //! Lowermost slab for that the calculation needs to be performed at a given time step
    int slab_first;
    //! Uppermost slab ...
    int slab_last;
    //! First slab for which ref. center is stored
    int slab_first_ref;
    //! Last ...
    int slab_last_ref;
    //! Slab buffer region around reference slabs
    int slab_buffer;
    //! First relevant atom for a slab
    int* firstatom;
    //! Last relevant atom for a slab
    int* lastatom;
    //! Gaussian-weighted slab center
    rvec* slab_center;
    //! Gaussian-weighted slab center for the reference positions
    rvec* slab_center_ref;
    //! Sum of gaussian weights in a slab
    real* slab_weights;
    //! Torque T = r x f for each slab. torque_v = m.v = angular momentum in the direction of v
    real* slab_torque_v;
    //! min_gaussian from t_rotgrp is the minimum value the gaussian must have so that the force is actually evaluated. max_beta is just another way to put it
    real max_beta;
    //! Precalculated gaussians for a single atom
    real* gn_atom;
    //! Tells to which slab each precalculated gaussian belongs
    int* gn_slabind;
    //! Inner sum of the flexible2 potential per slab; this is precalculated for optimization reasons
    rvec* slab_innersumvec;
    //! Holds atom positions and gaussian weights of atoms belonging to a slab
    gmx_slabdata* slab_data;

    /* For potential fits with varying angle: */
    //! Used for fit type 'potential'
    gmx_potfit* PotAngleFit;
};


//! Enforced rotation data for all groups
struct gmx_enfrot
{
    //! Input parameters.
    const t_rot* rot = nullptr;
    //! Output period for main rotation outfile
    int nstrout;
    //! Output period for per-slab data
    int nstsout;
    //! Output file for rotation data
    FILE* out_rot = nullptr;
    //! Output file for torque data
    FILE* out_torque = nullptr;
    //! Output file for slab angles for flexible type
    FILE* out_angles = nullptr;
    //! Output file for slab centers
    FILE* out_slabs = nullptr;
    //! Allocation size of buf
    int bufsize = 0;
    //! Coordinate buffer variable for sorting
    rvec* xbuf = nullptr;
    //! Masses buffer variable for sorting
    real* mbuf = nullptr;
    //! Buffer variable needed for position sorting
    sort_along_vec_t* data = nullptr;
    //! MPI buffer
    real* mpi_inbuf = nullptr;
    //! MPI buffer
    real* mpi_outbuf = nullptr;
    //! Allocation size of in & outbuf
    int mpi_bufsize = 0;
    //! If true, append output files
    gmx_bool restartWithAppending = false;
    //! Used to skip first output when appending to avoid duplicate entries in rotation outfiles
    gmx_bool bOut = false;
    //! Stores working data per group
    std::vector<gmx_enfrotgrp> enfrotgrp;
    ~gmx_enfrot();
};

gmx_enfrot::~gmx_enfrot()
{
    if (out_rot)
    {
        gmx_fio_fclose(out_rot);
    }
    if (out_slabs)
    {
        gmx_fio_fclose(out_slabs);
    }
    if (out_angles)
    {
        gmx_fio_fclose(out_angles);
    }
    if (out_torque)
    {
        gmx_fio_fclose(out_torque);
    }
}

namespace gmx
{

extern template LocalAtomSet LocalAtomSetManager::add<void, void>(ArrayRef<const int> globalAtomIndex);

class EnforcedRotation::Impl
{
public:
    gmx_enfrot enforcedRotation_;
};

EnforcedRotation::EnforcedRotation() : impl_(new Impl) {}

EnforcedRotation::~EnforcedRotation() = default;

gmx_enfrot* EnforcedRotation::getLegacyEnfrot()
{
    return &impl_->enforcedRotation_;
}

} // namespace gmx

/* Activate output of forces for correctness checks */
/* #define PRINT_FORCES */
#ifdef PRINT_FORCES
#    define PRINT_FORCE_J                       \
        fprintf(stderr,                         \
                "f%d = %15.8f %15.8f %15.8f\n", \
                erg->xc_ref_ind[j],             \
                erg->f_rot_loc[j][XX],          \
                erg->f_rot_loc[j][YY],          \
                erg->f_rot_loc[j][ZZ]);
#    define PRINT_POT_TAU                   \
        if (MAIN(cr))                       \
        {                                   \
            fprintf(stderr,                 \
                    "potential = %15.8f\n"  \
                    "torque    = %15.8f\n", \
                    erg->V,                 \
                    erg->torque_v);         \
        }
#else
#    define PRINT_FORCE_J
#    define PRINT_POT_TAU
#endif

/* Shortcuts for often used queries */
#define ISFLEX(rg)                                                                                         \
    (((rg)->eType == EnforcedRotationGroupType::Flex) || ((rg)->eType == EnforcedRotationGroupType::Flext) \
     || ((rg)->eType == EnforcedRotationGroupType::Flex2)                                                  \
     || ((rg)->eType == EnforcedRotationGroupType::Flex2t))
#define ISCOLL(rg)                                                                                         \
    (((rg)->eType == EnforcedRotationGroupType::Flex) || ((rg)->eType == EnforcedRotationGroupType::Flext) \
     || ((rg)->eType == EnforcedRotationGroupType::Flex2)                                                  \
     || ((rg)->eType == EnforcedRotationGroupType::Flex2t)                                                 \
     || ((rg)->eType == EnforcedRotationGroupType::Rmpf)                                                   \
     || ((rg)->eType == EnforcedRotationGroupType::Rm2pf))


/* Does any of the rotation groups use slab decomposition? */
static gmx_bool HaveFlexibleGroups(const t_rot* rot)
{
    for (const auto& grp : rot->grp)
    {
        if (ISFLEX(&grp))
        {
            return TRUE;
        }
    }

    return FALSE;
}


/* Is for any group the fit angle determined by finding the minimum of the
 * rotation potential? */
static gmx_bool HavePotFitGroups(const t_rot* rot)
{
    auto foundIt = std::find_if(rot->grp.begin(), rot->grp.end(), [](const auto& grp) {
        return RotationGroupFitting::Pot == grp.eFittype;
    });
    return foundIt != rot->grp.end();
}


static double** allocate_square_matrix(int dim)
{
    int      i;
    double** mat = nullptr;


    snew(mat, dim);
    for (i = 0; i < dim; i++)
    {
        snew(mat[i], dim);
    }

    return mat;
}


static void free_square_matrix(double** mat, int dim)
{
    int i;


    for (i = 0; i < dim; i++)
    {
        sfree(mat[i]);
    }
    sfree(mat);
}


/* Return the angle for which the potential is minimal */
static real get_fitangle(const gmx_enfrotgrp* erg)
{
    int         i;
    real        fitangle = -999.9;
    real        pot_min  = GMX_FLOAT_MAX;
    gmx_potfit* fit;


    fit = erg->PotAngleFit;

    for (i = 0; i < erg->rotg->PotAngle_nstep; i++)
    {
        if (fit->V[i] < pot_min)
        {
            pot_min  = fit->V[i];
            fitangle = fit->degangle[i];
        }
    }

    return fitangle;
}


/* Reduce potential angle fit data for this group at this time step? */
static inline gmx_bool bPotAngle(const gmx_enfrot* er, const t_rotgrp* rotg, int64_t step)
{
    return ((RotationGroupFitting::Pot == rotg->eFittype)
            && (do_per_step(step, er->nstsout) || do_per_step(step, er->nstrout)));
}

/* Reduce slab torqe data for this group at this time step? */
static inline gmx_bool bSlabTau(const gmx_enfrot* er, const t_rotgrp* rotg, int64_t step)
{
    return ((ISFLEX(rotg)) && do_per_step(step, er->nstsout));
}

/* Output rotation energy, torques, etc. for each rotation group */
static void reduce_output(const t_commrec* cr, gmx_enfrot* er, real t, int64_t step)
{
    int      i, islab, nslabs = 0;
    int      count; /* MPI element counter                               */
    real     fitangle;
    gmx_bool bFlex;

    /* Fill the MPI buffer with stuff to reduce. If items are added for reduction
     * here, the MPI buffer size has to be enlarged also in calc_mpi_bufsize() */
    if (PAR(cr))
    {
        count = 0;
        for (auto& ergRef : er->enfrotgrp)
        {
            gmx_enfrotgrp*  erg    = &ergRef;
            const t_rotgrp* rotg   = erg->rotg;
            nslabs                 = erg->slab_last - erg->slab_first + 1;
            er->mpi_inbuf[count++] = erg->V;
            er->mpi_inbuf[count++] = erg->torque_v;
            er->mpi_inbuf[count++] = erg->angle_v;
            er->mpi_inbuf[count++] =
                    erg->weight_v; /* weights are not needed for flex types, but this is just a single value */

            if (bPotAngle(er, rotg, step))
            {
                for (i = 0; i < rotg->PotAngle_nstep; i++)
                {
                    er->mpi_inbuf[count++] = erg->PotAngleFit->V[i];
                }
            }
            if (bSlabTau(er, rotg, step))
            {
                for (i = 0; i < nslabs; i++)
                {
                    er->mpi_inbuf[count++] = erg->slab_torque_v[i];
                }
            }
        }
        if (count > er->mpi_bufsize)
        {
            gmx_fatal(FARGS, "%s MPI buffer overflow, please report this error.", RotStr.c_str());
        }

#if GMX_MPI
        MPI_Reduce(er->mpi_inbuf, er->mpi_outbuf, count, GMX_MPI_REAL, MPI_SUM, MAINRANK(cr), cr->mpi_comm_mygroup);
#endif

        /* Copy back the reduced data from the buffer on the main */
        if (MAIN(cr))
        {
            count = 0;
            for (auto& ergRef : er->enfrotgrp)
            {
                gmx_enfrotgrp*  erg  = &ergRef;
                const t_rotgrp* rotg = erg->rotg;
                nslabs               = erg->slab_last - erg->slab_first + 1;
                erg->V               = er->mpi_outbuf[count++];
                erg->torque_v        = er->mpi_outbuf[count++];
                erg->angle_v         = er->mpi_outbuf[count++];
                erg->weight_v        = er->mpi_outbuf[count++];

                if (bPotAngle(er, rotg, step))
                {
                    for (int i = 0; i < rotg->PotAngle_nstep; i++)
                    {
                        erg->PotAngleFit->V[i] = er->mpi_outbuf[count++];
                    }
                }
                if (bSlabTau(er, rotg, step))
                {
                    for (int i = 0; i < nslabs; i++)
                    {
                        erg->slab_torque_v[i] = er->mpi_outbuf[count++];
                    }
                }
            }
        }
    }

    /* Output */
    if (MAIN(cr))
    {
        /* Angle and torque for each rotation group */
        for (auto& ergRef : er->enfrotgrp)
        {
            gmx_enfrotgrp*  erg  = &ergRef;
            const t_rotgrp* rotg = erg->rotg;
            bFlex                = ISFLEX(rotg);

            /* Output to main rotation output file: */
            if (do_per_step(step, er->nstrout))
            {
                if (RotationGroupFitting::Pot == rotg->eFittype)
                {
                    fitangle = get_fitangle(erg);
                }
                else
                {
                    if (bFlex)
                    {
                        fitangle = erg->angle_v; /* RMSD fit angle */
                    }
                    else
                    {
                        fitangle = (erg->angle_v / erg->weight_v) * 180.0 * M_1_PI;
                    }
                }
                fprintf(er->out_rot, "%12.4f", fitangle);
                fprintf(er->out_rot, "%12.3e", erg->torque_v);
                fprintf(er->out_rot, "%12.3e", erg->V);
            }

            if (do_per_step(step, er->nstsout))
            {
                /* Output to torque log file: */
                if (bFlex)
                {
                    fprintf(er->out_torque, "%12.3e%6d", t, erg->groupIndex);
                    for (int i = erg->slab_first; i <= erg->slab_last; i++)
                    {
                        islab = i - erg->slab_first; /* slab index */
                        /* Only output if enough weight is in slab */
                        if (erg->slab_weights[islab] > rotg->min_gaussian)
                        {
                            fprintf(er->out_torque, "%6d%12.3e", i, erg->slab_torque_v[islab]);
                        }
                    }
                    fprintf(er->out_torque, "\n");
                }

                /* Output to angles log file: */
                if (RotationGroupFitting::Pot == rotg->eFittype)
                {
                    fprintf(er->out_angles, "%12.3e%6d%12.4f", t, erg->groupIndex, erg->degangle);
                    /* Output energies at a set of angles around the reference angle */
                    for (int i = 0; i < rotg->PotAngle_nstep; i++)
                    {
                        fprintf(er->out_angles, "%12.3e", erg->PotAngleFit->V[i]);
                    }
                    fprintf(er->out_angles, "\n");
                }
            }
        }
        if (do_per_step(step, er->nstrout))
        {
            fprintf(er->out_rot, "\n");
        }
    }
}


/* Add the forces from enforced rotation potential to the local forces.
 * Should be called after the SR forces have been evaluated */
real add_rot_forces(gmx_enfrot* er, gmx::ArrayRef<gmx::RVec> force, const t_commrec* cr, int64_t step, real t)
{
    real Vrot = 0.0; /* If more than one rotation group is present, Vrot
                        assembles the local parts from all groups         */

    /* Loop over enforced rotation groups (usually 1, though)
     * Apply the forces from rotation potentials */
    for (auto& ergRef : er->enfrotgrp)
    {
        gmx_enfrotgrp* erg = &ergRef;
        Vrot += erg->V; /* add the local parts from the nodes */
        const auto& localRotationGroupIndex = erg->atomSet->localIndex();
        for (gmx::Index l = 0; l < localRotationGroupIndex.ssize(); l++)
        {
            /* Get the right index of the local force */
            int ii = localRotationGroupIndex[l];
            /* Add */
            rvec_inc(force[ii], erg->f_rot_loc[l]);
        }
    }

    /* Reduce energy,torque, angles etc. to get the sum values (per rotation group)
     * on the main and output these values to file. */
    if ((do_per_step(step, er->nstrout) || do_per_step(step, er->nstsout)) && er->bOut)
    {
        reduce_output(cr, er, t, step);
    }

    /* When appending, er->bOut is FALSE the first time to avoid duplicate entries */
    er->bOut = TRUE;

    PRINT_POT_TAU

    return Vrot;
}


/* The Gaussian norm is chosen such that the sum of the gaussian functions
 * over the slabs is approximately 1.0 everywhere */
#define GAUSS_NORM 0.569917543430618


/* Calculate the maximum beta that leads to a gaussian larger min_gaussian,
 * also does some checks
 */
static double calc_beta_max(real min_gaussian, real slab_dist)
{
    double sigma;
    double arg;


    /* Actually the next two checks are already made in grompp */
    if (slab_dist <= 0)
    {
        gmx_fatal(FARGS, "Slab distance of flexible rotation groups must be >=0 !");
    }
    if (min_gaussian <= 0)
    {
        gmx_fatal(FARGS, "Cutoff value for Gaussian must be > 0. (You requested %f)", min_gaussian);
    }

    /* Define the sigma value */
    sigma = 0.7 * slab_dist;

    /* Calculate the argument for the logarithm and check that the log() result is negative or 0 */
    arg = min_gaussian / GAUSS_NORM;
    if (arg > 1.0)
    {
        gmx_fatal(FARGS, "min_gaussian of flexible rotation groups must be <%g", GAUSS_NORM);
    }

    return std::sqrt(-2.0 * sigma * sigma * std::log(min_gaussian / GAUSS_NORM));
}


static inline real calc_beta(rvec curr_x, const gmx_enfrotgrp* erg, int n)
{
    return iprod(curr_x, erg->vec) - erg->rotg->slab_dist * n;
}


static inline real gaussian_weight(rvec curr_x, const gmx_enfrotgrp* erg, int n)
{
    const real norm = GAUSS_NORM;
    real       sigma;


    /* Define the sigma value */
    sigma = 0.7 * erg->rotg->slab_dist;
    /* Calculate the Gaussian value of slab n for position curr_x */
    return norm * std::exp(-0.5 * gmx::square(calc_beta(curr_x, erg, n) / sigma));
}


/* Returns the weight in a single slab, also calculates the Gaussian- and mass-
 * weighted sum of positions for that slab */
static real get_slab_weight(int                            j,
                            const gmx_enfrotgrp*           erg,
                            gmx::ArrayRef<const gmx::RVec> xc,
                            const real                     mc[],
                            rvec*                          x_weighted_sum)
{
    rvec curr_x;           /* The position of an atom                      */
    rvec curr_x_weighted;  /* The gaussian-weighted position               */
    real gaussian;         /* A single gaussian weight                     */
    real wgauss;           /* gaussian times current mass                  */
    real slabweight = 0.0; /* The sum of weights in the slab               */

    clear_rvec(*x_weighted_sum);

    /* Loop over all atoms in the rotation group */
    for (int i = 0; i < erg->rotg->nat; i++)
    {
        copy_rvec(xc[i], curr_x);
        gaussian = gaussian_weight(curr_x, erg, j);
        wgauss   = gaussian * mc[i];
        svmul(wgauss, curr_x, curr_x_weighted);
        rvec_add(*x_weighted_sum, curr_x_weighted, *x_weighted_sum);
        slabweight += wgauss;
    } /* END of loop over rotation group atoms */

    return slabweight;
}


static void get_slab_centers(gmx_enfrotgrp* erg, /* Enforced rotation group working data */
                             gmx::ArrayRef<const gmx::RVec> xc, /* The rotation group positions;
                                                   will typically be enfrotgrp->xc, but at first
                                                   call it is enfrotgrp->xc_ref */
                             real*    mc,         /* The masses of the rotation group atoms       */
                             real     time,       /* Used for output only                         */
                             FILE*    out_slabs,  /* For outputting center per slab information   */
                             gmx_bool bOutStep,   /* Is this an output step?                      */
                             gmx_bool bReference) /* If this routine is called from
                                                     init_rot_group we need to store
                                                     the reference slab centers                   */
{
    /* Loop over slabs */
    for (int j = erg->slab_first; j <= erg->slab_last; j++)
    {
        int slabIndex                = j - erg->slab_first;
        erg->slab_weights[slabIndex] = get_slab_weight(j, erg, xc, mc, &erg->slab_center[slabIndex]);

        /* We can do the calculations ONLY if there is weight in the slab! */
        if (erg->slab_weights[slabIndex] > WEIGHT_MIN)
        {
            svmul(1.0 / erg->slab_weights[slabIndex], erg->slab_center[slabIndex], erg->slab_center[slabIndex]);
        }
        else
        {
            /* We need to check this here, since we divide through slab_weights
             * in the flexible low-level routines! */
            gmx_fatal(FARGS, "Not enough weight in slab %d. Slab center cannot be determined!", j);
        }

        /* At first time step: save the centers of the reference structure */
        if (bReference)
        {
            copy_rvec(erg->slab_center[slabIndex], erg->slab_center_ref[slabIndex]);
        }
    } /* END of loop over slabs */

    /* Output on the main */
    if ((nullptr != out_slabs) && bOutStep)
    {
        fprintf(out_slabs, "%12.3e%6d", time, erg->groupIndex);
        for (int j = erg->slab_first; j <= erg->slab_last; j++)
        {
            int slabIndex = j - erg->slab_first;
            fprintf(out_slabs,
                    "%6d%12.3e%12.3e%12.3e",
                    j,
                    erg->slab_center[slabIndex][XX],
                    erg->slab_center[slabIndex][YY],
                    erg->slab_center[slabIndex][ZZ]);
        }
        fprintf(out_slabs, "\n");
    }
}


static void calc_rotmat(const rvec vec,
                        real   degangle, /* Angle alpha of rotation at time t in degrees       */
                        matrix rotmat)   /* Rotation matrix                                    */
{
    real radangle;            /* Rotation angle in radians */
    real cosa;                /* cosine alpha              */
    real sina;                /* sine alpha                */
    real OMcosa;              /* 1 - cos(alpha)            */
    real dumxy, dumxz, dumyz; /* save computations         */
    rvec rot_vec;             /* Rotate around rot_vec ... */


    radangle = degangle * M_PI / 180.0;
    copy_rvec(vec, rot_vec);

    /* Precompute some variables: */
    cosa   = std::cos(radangle);
    sina   = std::sin(radangle);
    OMcosa = 1.0 - cosa;
    dumxy  = rot_vec[XX] * rot_vec[YY] * OMcosa;
    dumxz  = rot_vec[XX] * rot_vec[ZZ] * OMcosa;
    dumyz  = rot_vec[YY] * rot_vec[ZZ] * OMcosa;

    /* Construct the rotation matrix for this rotation group: */
    /* 1st column: */
    rotmat[XX][XX] = cosa + rot_vec[XX] * rot_vec[XX] * OMcosa;
    rotmat[YY][XX] = dumxy + rot_vec[ZZ] * sina;
    rotmat[ZZ][XX] = dumxz - rot_vec[YY] * sina;
    /* 2nd column: */
    rotmat[XX][YY] = dumxy - rot_vec[ZZ] * sina;
    rotmat[YY][YY] = cosa + rot_vec[YY] * rot_vec[YY] * OMcosa;
    rotmat[ZZ][YY] = dumyz + rot_vec[XX] * sina;
    /* 3rd column: */
    rotmat[XX][ZZ] = dumxz + rot_vec[YY] * sina;
    rotmat[YY][ZZ] = dumyz - rot_vec[XX] * sina;
    rotmat[ZZ][ZZ] = cosa + rot_vec[ZZ] * rot_vec[ZZ] * OMcosa;

#ifdef PRINTMATRIX
    int iii, jjj;

    for (iii = 0; iii < 3; iii++)
    {
        for (jjj = 0; jjj < 3; jjj++)
        {
            fprintf(stderr, " %10.8f ", rotmat[iii][jjj]);
        }
        fprintf(stderr, "\n");
    }
#endif
}


/* Calculates torque on the rotation axis tau = position x force */
static inline real torque(const rvec rotvec, /* rotation vector; MUST be normalized!          */
                          rvec force, /* force                                                */
                          rvec x,     /* position of atom on which the force acts             */
                          rvec pivot) /* pivot point of rotation axis                         */
{
    rvec vectmp, tau;


    /* Subtract offset */
    rvec_sub(x, pivot, vectmp);

    /* position x force */
    cprod(vectmp, force, tau);

    /* Return the part of the torque which is parallel to the rotation vector */
    return iprod(tau, rotvec);
}


/* Right-aligned output of value with standard width */
static void print_aligned(FILE* fp, char const* str)
{
    fprintf(fp, "%12s", str);
}


/* Right-aligned output of value with standard short width */
static void print_aligned_short(FILE* fp, char const* str)
{
    fprintf(fp, "%6s", str);
}


static FILE* open_output_file(const char* fn, int steps, const char what[])
{
    FILE* fp;


    fp = gmx_ffopen(fn, "w");

    fprintf(fp, "# Output of %s is written in intervals of %d time step%s.\n#\n", what, steps, steps > 1 ? "s" : "");

    return fp;
}


/* Open output file for slab center data. Call on main only */
static FILE* open_slab_out(const char* fn, gmx_enfrot* er)
{
    FILE* fp;

    if (er->restartWithAppending)
    {
        fp = gmx_fio_fopen(fn, "a");
    }
    else
    {
        fp = open_output_file(fn, er->nstsout, "gaussian weighted slab centers");

        for (auto& ergRef : er->enfrotgrp)
        {
            gmx_enfrotgrp* erg = &ergRef;
            if (ISFLEX(erg->rotg))
            {
                fprintf(fp,
                        "# Rotation group %d (%s), slab distance %f nm, %s.\n",
                        erg->groupIndex,
                        enumValueToString(erg->rotg->eType),
                        erg->rotg->slab_dist,
                        erg->rotg->bMassW ? "centers of mass" : "geometrical centers");
            }
        }

        fprintf(fp, "# Reference centers are listed first (t=-1).\n");
        fprintf(fp, "# The following columns have the syntax:\n");
        fprintf(fp, "#     ");
        print_aligned_short(fp, "t");
        print_aligned_short(fp, "grp");
        /* Print legend for the first two entries only ... */
        for (int i = 0; i < 2; i++)
        {
            print_aligned_short(fp, "slab");
            print_aligned(fp, "X center");
            print_aligned(fp, "Y center");
            print_aligned(fp, "Z center");
        }
        fprintf(fp, " ...\n");
        fflush(fp);
    }

    return fp;
}


/* Adds 'buf' to 'str' */
static void add_to_string(char** str, char* buf)
{
    int len;


    len = strlen(*str) + strlen(buf) + 1;
    srenew(*str, len);
    strcat(*str, buf);
}


static void add_to_string_aligned(char** str, char* buf)
{
    char buf_aligned[STRLEN];

    sprintf(buf_aligned, "%12s", buf);
    add_to_string(str, buf_aligned);
}


/* Open output file and print some general information about the rotation groups.
 * Call on main only */
static FILE* open_rot_out(const char* fn, const gmx_output_env_t* oenv, gmx_enfrot* er)
{
    FILE*                    fp;
    std::vector<std::string> setname;
    char                     buf[50];
    gmx_bool                 bFlex;
    char*                    LegendStr = nullptr;
    const t_rot*             rot       = er->rot;

    if (er->restartWithAppending)
    {
        fp = gmx_fio_fopen(fn, "a");
    }
    else
    {
        fp = xvgropen(fn,
                      "Rotation angles and energy",
                      "Time (ps)",
                      "angles (degrees) and energies (kJ/mol)",
                      oenv);
        fprintf(fp,
                "# Output of enforced rotation data is written in intervals of %d time "
                "step%s.\n#\n",
                er->nstrout,
                er->nstrout > 1 ? "s" : "");
        fprintf(fp,
                "# The scalar tau is the torque (kJ/mol) in the direction of the rotation vector "
                "v.\n");
        fprintf(fp, "# To obtain the vectorial torque, multiply tau with the group's rot-vec.\n");
        fprintf(fp,
                "# For flexible groups, tau(t,n) from all slabs n have been summed in a single "
                "value tau(t) here.\n");
        fprintf(fp, "# The torques tau(t,n) are found in the rottorque.log (-rt) output file\n");

        for (int g = 0; g < gmx::ssize(rot->grp); g++)
        {
            const t_rotgrp*      rotg = &rot->grp[g];
            const gmx_enfrotgrp* erg  = &er->enfrotgrp[g];
            bFlex                     = ISFLEX(rotg);

            fprintf(fp, "#\n");
            fprintf(fp, "# ROTATION GROUP %d, potential type '%s':\n", g, enumValueToString(rotg->eType));
            fprintf(fp, "# rot-massw%d          %s\n", g, booleanValueToString(rotg->bMassW));
            fprintf(fp,
                    "# rot-vec%d            %12.5e %12.5e %12.5e\n",
                    g,
                    erg->vec[XX],
                    erg->vec[YY],
                    erg->vec[ZZ]);
            fprintf(fp, "# rot-rate%d           %12.5e degrees/ps\n", g, rotg->rate);
            fprintf(fp, "# rot-k%d              %12.5e kJ/(mol*nm^2)\n", g, rotg->k);
            if (rotg->eType == EnforcedRotationGroupType::Iso || rotg->eType == EnforcedRotationGroupType::Pm
                || rotg->eType == EnforcedRotationGroupType::Rm
                || rotg->eType == EnforcedRotationGroupType::Rm2)
            {
                fprintf(fp,
                        "# rot-pivot%d          %12.5e %12.5e %12.5e  nm\n",
                        g,
                        rotg->pivot[XX],
                        rotg->pivot[YY],
                        rotg->pivot[ZZ]);
            }

            if (bFlex)
            {
                fprintf(fp, "# rot-slab-distance%d   %f nm\n", g, rotg->slab_dist);
                fprintf(fp, "# rot-min-gaussian%d   %12.5e\n", g, rotg->min_gaussian);
            }

            /* Output the centers of the rotation groups for the pivot-free potentials */
            if ((rotg->eType == EnforcedRotationGroupType::Isopf)
                || (rotg->eType == EnforcedRotationGroupType::Pmpf)
                || (rotg->eType == EnforcedRotationGroupType::Rmpf)
                || (rotg->eType == EnforcedRotationGroupType::Rm2pf
                    || (rotg->eType == EnforcedRotationGroupType::Flext)
                    || (rotg->eType == EnforcedRotationGroupType::Flex2t)))
            {
                fprintf(fp,
                        "# ref. grp. %d center  %12.5e %12.5e %12.5e\n",
                        g,
                        erg->xc_ref_center[XX],
                        erg->xc_ref_center[YY],
                        erg->xc_ref_center[ZZ]);

                fprintf(fp,
                        "# grp. %d init.center  %12.5e %12.5e %12.5e\n",
                        g,
                        erg->xc_center[XX],
                        erg->xc_center[YY],
                        erg->xc_center[ZZ]);
            }

            if ((rotg->eType == EnforcedRotationGroupType::Rm2)
                || (rotg->eType == EnforcedRotationGroupType::Flex2)
                || (rotg->eType == EnforcedRotationGroupType::Flex2t))
            {
                fprintf(fp, "# rot-eps%d            %12.5e nm^2\n", g, rotg->eps);
            }
            if (RotationGroupFitting::Pot == rotg->eFittype)
            {
                fprintf(fp, "#\n");
                fprintf(fp,
                        "# theta_fit%d is determined by first evaluating the potential for %d "
                        "angles around theta_ref%d.\n",
                        g,
                        rotg->PotAngle_nstep,
                        g);
                fprintf(fp,
                        "# The fit angle is the one with the smallest potential. It is given as "
                        "the deviation\n");
                fprintf(fp,
                        "# from the reference angle, i.e. if theta_ref=X and theta_fit=Y, then the "
                        "angle with\n");
                fprintf(fp,
                        "# minimal value of the potential is X+Y. Angular resolution is %g "
                        "degrees.\n",
                        rotg->PotAngle_step);
            }
        }

        /* Print a nice legend */
        snew(LegendStr, 1);
        LegendStr[0] = '\0';
        sprintf(buf, "#     %6s", "time");
        add_to_string_aligned(&LegendStr, buf);

        for (int g = 0; g < gmx::ssize(rot->grp); g++)
        {
            sprintf(buf, "theta_ref%d", g);
            add_to_string_aligned(&LegendStr, buf);

            setname.emplace_back(gmx::formatString("%s (degrees)", buf));
        }
        for (int g = 0; g < gmx::ssize(rot->grp); g++)
        {
            const t_rotgrp* rotg = &rot->grp[g];
            bFlex                = ISFLEX(rotg);

            /* For flexible axis rotation we use RMSD fitting to determine the
             * actual angle of the rotation group */
            if (bFlex || RotationGroupFitting::Pot == rotg->eFittype)
            {
                sprintf(buf, "theta_fit%d", g);
            }
            else
            {
                sprintf(buf, "theta_av%d", g);
            }
            add_to_string_aligned(&LegendStr, buf);
            setname.emplace_back(gmx::formatString("%s (degrees)", buf));

            sprintf(buf, "tau%d", g);
            add_to_string_aligned(&LegendStr, buf);
            setname.emplace_back(gmx::formatString("%s (kJ/mol)", buf));

            sprintf(buf, "energy%d", g);
            add_to_string_aligned(&LegendStr, buf);
            setname.emplace_back(gmx::formatString("%s (kJ/mol)", buf));
        }
        fprintf(fp, "#\n");

        if (setname.size() > 1)
        {
            xvgrLegend(fp, setname, oenv);
        }
        fprintf(fp, "#\n# Legend for the following data columns:\n");
        fprintf(fp, "%s\n", LegendStr);
        sfree(LegendStr);

        fflush(fp);
    }

    return fp;
}


/* Call on main only */
static FILE* open_angles_out(const char* fn, gmx_enfrot* er)
{
    FILE*        fp;
    char         buf[100];
    const t_rot* rot = er->rot;

    if (er->restartWithAppending)
    {
        fp = gmx_fio_fopen(fn, "a");
    }
    else
    {
        /* Open output file and write some information about it's structure: */
        fp = open_output_file(fn, er->nstsout, "rotation group angles");
        fprintf(fp, "# All angles given in degrees, time in ps.\n");
        for (int g = 0; g < gmx::ssize(rot->grp); g++)
        {
            const t_rotgrp*      rotg = &rot->grp[g];
            const gmx_enfrotgrp* erg  = &er->enfrotgrp[g];

            /* Output for this group happens only if potential type is flexible or
             * if fit type is potential! */
            if (ISFLEX(rotg) || (RotationGroupFitting::Pot == rotg->eFittype))
            {
                if (ISFLEX(rotg))
                {
                    sprintf(buf, " slab distance %f nm, ", rotg->slab_dist);
                }
                else
                {
                    buf[0] = '\0';
                }

                fprintf(fp,
                        "#\n# ROTATION GROUP %d '%s',%s fit type '%s'.\n",
                        g,
                        enumValueToString(rotg->eType),
                        buf,
                        enumValueToString(rotg->eFittype));

                /* Special type of fitting using the potential minimum. This is
                 * done for the whole group only, not for the individual slabs. */
                if (RotationGroupFitting::Pot == rotg->eFittype)
                {
                    fprintf(fp,
                            "#    To obtain theta_fit%d, the potential is evaluated for %d angles "
                            "around theta_ref%d\n",
                            g,
                            rotg->PotAngle_nstep,
                            g);
                    fprintf(fp,
                            "#    The fit angle in the rotation standard outfile is the one with "
                            "minimal energy E(theta_fit) [kJ/mol].\n");
                    fprintf(fp, "#\n");
                }

                fprintf(fp, "# Legend for the group %d data columns:\n", g);
                fprintf(fp, "#     ");
                print_aligned_short(fp, "time");
                print_aligned_short(fp, "grp");
                print_aligned(fp, "theta_ref");

                if (RotationGroupFitting::Pot == rotg->eFittype)
                {
                    /* Output the set of angles around the reference angle */
                    for (int i = 0; i < rotg->PotAngle_nstep; i++)
                    {
                        sprintf(buf, "E(%g)", erg->PotAngleFit->degangle[i]);
                        print_aligned(fp, buf);
                    }
                }
                else
                {
                    /* Output fit angle for each slab */
                    print_aligned_short(fp, "slab");
                    print_aligned_short(fp, "atoms");
                    print_aligned(fp, "theta_fit");
                    print_aligned_short(fp, "slab");
                    print_aligned_short(fp, "atoms");
                    print_aligned(fp, "theta_fit");
                    fprintf(fp, " ...");
                }
                fprintf(fp, "\n");
            }
        }
        fflush(fp);
    }

    return fp;
}


/* Open torque output file and write some information about it's structure.
 * Call on main only */
static FILE* open_torque_out(const char* fn, gmx_enfrot* er)
{
    FILE*        fp;
    const t_rot* rot = er->rot;

    if (er->restartWithAppending)
    {
        fp = gmx_fio_fopen(fn, "a");
    }
    else
    {
        fp = open_output_file(fn, er->nstsout, "torques");

        for (int g = 0; g < gmx::ssize(rot->grp); g++)
        {
            const t_rotgrp*      rotg = &rot->grp[g];
            const gmx_enfrotgrp* erg  = &er->enfrotgrp[g];
            if (ISFLEX(rotg))
            {
                fprintf(fp,
                        "# Rotation group %d (%s), slab distance %f nm.\n",
                        g,
                        enumValueToString(rotg->eType),
                        rotg->slab_dist);
                fprintf(fp,
                        "# The scalar tau is the torque (kJ/mol) in the direction of the rotation "
                        "vector.\n");
                fprintf(fp, "# To obtain the vectorial torque, multiply tau with\n");
                fprintf(fp,
                        "# rot-vec%d            %10.3e %10.3e %10.3e\n",
                        g,
                        erg->vec[XX],
                        erg->vec[YY],
                        erg->vec[ZZ]);
                fprintf(fp, "#\n");
            }
        }
        fprintf(fp, "# Legend for the following data columns: (tau=torque for that slab):\n");
        fprintf(fp, "#     ");
        print_aligned_short(fp, "t");
        print_aligned_short(fp, "grp");
        print_aligned_short(fp, "slab");
        print_aligned(fp, "tau");
        print_aligned_short(fp, "slab");
        print_aligned(fp, "tau");
        fprintf(fp, " ...\n");
        fflush(fp);
    }

    return fp;
}


static void swap_val(double* vec, int i, int j)
{
    double tmp = vec[j];


    vec[j] = vec[i];
    vec[i] = tmp;
}


static void swap_col(double** mat, int i, int j)
{
    double tmp[3] = { mat[0][j], mat[1][j], mat[2][j] };


    mat[0][j] = mat[0][i];
    mat[1][j] = mat[1][i];
    mat[2][j] = mat[2][i];

    mat[0][i] = tmp[0];
    mat[1][i] = tmp[1];
    mat[2][i] = tmp[2];
}


/* Eigenvectors are stored in columns of eigen_vec */
static void diagonalize_symmetric(double** mat, double** eigen_vec, double eigenval[3])
{
    int n_rot;


    jacobi(mat, 3, eigenval, eigen_vec, &n_rot);

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


static void align_with_z(rvec* s, /* Structure to align */
                         int   natoms,
                         rvec  axis)
{
    int    i, j, k;
    rvec   zet         = { 0.0, 0.0, 1.0 };
    rvec   rot_axis    = { 0.0, 0.0, 0.0 };
    rvec*  rotated_str = nullptr;
    real   ooanorm;
    real   angle;
    matrix rotmat;


    snew(rotated_str, natoms);

    /* Normalize the axis */
    ooanorm = 1.0 / norm(axis);
    svmul(ooanorm, axis, axis);

    /* Calculate the angle for the fitting procedure */
    cprod(axis, zet, rot_axis);
    angle = std::acos(axis[2]);
    if (angle < 0.0)
    {
        angle += M_PI;
    }

    /* Calculate the rotation matrix */
    calc_rotmat(rot_axis, angle * 180.0 / M_PI, rotmat);

    /* Apply the rotation matrix to s */
    for (i = 0; i < natoms; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
            {
                rotated_str[i][j] += rotmat[j][k] * s[i][k];
            }
        }
    }

    /* Rewrite the rotated structure to s */
    for (i = 0; i < natoms; i++)
    {
        for (j = 0; j < 3; j++)
        {
            s[i][j] = rotated_str[i][j];
        }
    }

    sfree(rotated_str);
}


static void calc_correl_matrix(rvec* Xstr, rvec* Ystr, double** Rmat, int natoms)
{
    int i, j, k;


    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            Rmat[i][j] = 0.0;
        }
    }

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < natoms; k++)
            {
                Rmat[i][j] += Ystr[k][i] * Xstr[k][j];
            }
        }
    }
}


static void weigh_coords(rvec* str, real* weight, int natoms)
{
    int i, j;


    for (i = 0; i < natoms; i++)
    {
        for (j = 0; j < 3; j++)
        {
            str[i][j] *= std::sqrt(weight[i]);
        }
    }
}


static real opt_angle_analytic(rvec*      ref_s,
                               rvec*      act_s,
                               real*      weight,
                               int        natoms,
                               const rvec ref_com,
                               const rvec act_com,
                               rvec       axis)
{
    int      i, j, k;
    rvec*    ref_s_1 = nullptr;
    rvec*    act_s_1 = nullptr;
    rvec     shift;
    double **Rmat, **RtR, **eigvec;
    double   eigval[3];
    double   V[3][3], WS[3][3];
    double   rot_matrix[3][3];
    double   opt_angle;


    /* Do not change the original coordinates */
    snew(ref_s_1, natoms);
    snew(act_s_1, natoms);
    for (i = 0; i < natoms; i++)
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

    for (i = 0; i < natoms; i++)
    {
        ref_s_1[i][2] = 0.0;
        act_s_1[i][2] = 0.0;
    }

    /* Weight positions with sqrt(weight) */
    if (nullptr != weight)
    {
        weigh_coords(ref_s_1, weight, natoms);
        weigh_coords(act_s_1, weight, natoms);
    }

    /* Calculate correlation matrices R=YXt (X=ref_s; Y=act_s) */
    calc_correl_matrix(ref_s_1, act_s_1, Rmat, natoms);

    /* Calculate RtR */
    RtR = allocate_square_matrix(3);
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
            {
                RtR[i][j] += Rmat[k][i] * Rmat[k][j];
            }
        }
    }
    /* Diagonalize RtR */
    snew(eigvec, 3);
    for (i = 0; i < 3; i++)
    {
        snew(eigvec[i], 3);
    }

    diagonalize_symmetric(RtR, eigvec, eigval);
    swap_col(eigvec, 0, 1);
    swap_col(eigvec, 1, 2);
    swap_val(eigval, 0, 1);
    swap_val(eigval, 1, 2);

    /* Calculate V */
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            V[i][j]  = 0.0;
            WS[i][j] = 0.0;
        }
    }

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            WS[i][j] = eigvec[i][j] / std::sqrt(eigval[j]);
        }
    }

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
            {
                V[i][j] += Rmat[i][k] * WS[k][j];
            }
        }
    }
    free_square_matrix(Rmat, 3);

    /* Calculate optimal rotation matrix */
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            rot_matrix[i][j] = 0.0;
        }
    }

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
            {
                rot_matrix[i][j] += eigvec[i][k] * V[j][k];
            }
        }
    }
    rot_matrix[2][2] = 1.0;

    /* In some cases abs(rot_matrix[0][0]) can be slighly larger
     * than unity due to numerical inacurracies. To be able to calculate
     * the acos function, we put these values back in range. */
    if (rot_matrix[0][0] > 1.0)
    {
        rot_matrix[0][0] = 1.0;
    }
    else if (rot_matrix[0][0] < -1.0)
    {
        rot_matrix[0][0] = -1.0;
    }

    /* Determine the optimal rotation angle: */
    opt_angle = (-1.0) * std::acos(rot_matrix[0][0]) * 180.0 / M_PI;
    if (rot_matrix[0][1] < 0.0)
    {
        opt_angle = (-1.0) * opt_angle;
    }

    /* Give back some memory */
    free_square_matrix(RtR, 3);
    sfree(ref_s_1);
    sfree(act_s_1);
    for (i = 0; i < 3; i++)
    {
        sfree(eigvec[i]);
    }
    sfree(eigvec);

    return static_cast<real>(opt_angle);
}


/* Determine angle of the group by RMSD fit to the reference */
/* Not parallelized, call this routine only on the main */
static real flex_fit_angle(gmx_enfrotgrp* erg)
{
    rvec* fitcoords = nullptr;
    rvec  center;   /* Center of positions passed to the fit routine */
    real  fitangle; /* Angle of the rotation group derived by fitting */
    rvec  coord;
    real  scal;

    /* Get the center of the rotation group.
     * Note, again, erg->xc has been sorted in do_flexible */
    get_center(erg->xc, erg->mc_sorted, erg->rotg->nat, center);

    /* === Determine the optimal fit angle for the rotation group === */
    if (erg->rotg->eFittype == RotationGroupFitting::Norm)
    {
        /* Normalize every position to it's reference length */
        for (int i = 0; i < erg->rotg->nat; i++)
        {
            /* Put the center of the positions into the origin */
            rvec_sub(erg->xc[i], center, coord);
            /* Determine the scaling factor for the length: */
            scal = erg->xc_ref_length[erg->xc_sortind[i]] / norm(coord);
            /* Get position, multiply with the scaling factor and save  */
            svmul(scal, coord, erg->xc_norm[i]);
        }
        fitcoords = erg->xc_norm;
    }
    else
    {
        fitcoords = erg->xc;
    }
    /* From the point of view of the current positions, the reference has rotated
     * backwards. Since we output the angle relative to the fixed reference,
     * we need the minus sign. */
    fitangle = -opt_angle_analytic(
            erg->xc_ref_sorted, fitcoords, erg->mc_sorted, erg->rotg->nat, erg->xc_ref_center, center, erg->vec);

    return fitangle;
}


/* Determine actual angle of each slab by RMSD fit to the reference */
/* Not parallelized, call this routine only on the main */
static void flex_fit_angle_perslab(gmx_enfrotgrp* erg, double t, real degangle, FILE* fp)
{
    rvec curr_x, ref_x;
    rvec act_center; /* Center of actual positions that are passed to the fit routine */
    rvec ref_center; /* Same for the reference positions */
    real fitangle;   /* Angle of a slab derived from an RMSD fit to
                      * the reference structure at t=0  */
    gmx_slabdata* sd;
    real          OOm_av; /* 1/average_mass of a rotation group atom */
    real          m_rel;  /* Relative mass of a rotation group atom  */


    /* Average mass of a rotation group atom: */
    OOm_av = erg->invmass * erg->rotg->nat;

    /**********************************/
    /* First collect the data we need */
    /**********************************/

    /* Collect the data for the individual slabs */
    for (int n = erg->slab_first; n <= erg->slab_last; n++)
    {
        int slabIndex = n - erg->slab_first; /* slab index */
        sd            = &(erg->slab_data[slabIndex]);
        sd->nat       = erg->lastatom[slabIndex] - erg->firstatom[slabIndex] + 1;
        int ind       = 0;

        /* Loop over the relevant atoms in the slab */
        for (int l = erg->firstatom[slabIndex]; l <= erg->lastatom[slabIndex]; l++)
        {
            /* Current position of this atom: x[ii][XX/YY/ZZ] */
            copy_rvec(erg->xc[l], curr_x);

            /* The (unrotated) reference position of this atom is copied to ref_x.
             * Beware, the xc coords have been sorted in do_flexible */
            copy_rvec(erg->xc_ref_sorted[l], ref_x);

            /* Save data for doing angular RMSD fit later */
            /* Save the current atom position */
            copy_rvec(curr_x, sd->x[ind]);
            /* Save the corresponding reference position */
            copy_rvec(ref_x, sd->ref[ind]);

            /* Maybe also mass-weighting was requested. If yes, additionally
             * multiply the weights with the relative mass of the atom. If not,
             * multiply with unity. */
            m_rel = erg->mc_sorted[l] * OOm_av;

            /* Save the weight for this atom in this slab */
            sd->weight[ind] = gaussian_weight(curr_x, erg, n) * m_rel;

            /* Next atom in this slab */
            ind++;
        }
    }

    /******************************/
    /* Now do the fit calculation */
    /******************************/

    fprintf(fp, "%12.3e%6d%12.3f", t, erg->groupIndex, degangle);

    /* === Now do RMSD fitting for each slab === */
    /* We require at least SLAB_MIN_ATOMS in a slab, such that the fit makes sense. */
#define SLAB_MIN_ATOMS 4

    for (int n = erg->slab_first; n <= erg->slab_last; n++)
    {
        int slabIndex = n - erg->slab_first; /* slab index */
        sd            = &(erg->slab_data[slabIndex]);
        if (sd->nat >= SLAB_MIN_ATOMS)
        {
            /* Get the center of the slabs reference and current positions */
            get_center(sd->ref, sd->weight, sd->nat, ref_center);
            get_center(sd->x, sd->weight, sd->nat, act_center);
            if (erg->rotg->eFittype == RotationGroupFitting::Norm)
            {
                /* Normalize every position to it's reference length
                 * prior to performing the fit */
                for (int i = 0; i < sd->nat; i++) /* Center */
                {
                    rvec_dec(sd->ref[i], ref_center);
                    rvec_dec(sd->x[i], act_center);
                    /* Normalize x_i such that it gets the same length as ref_i */
                    svmul(norm(sd->ref[i]) / norm(sd->x[i]), sd->x[i], sd->x[i]);
                }
                /* We already subtracted the centers */
                clear_rvec(ref_center);
                clear_rvec(act_center);
            }
            fitangle = -opt_angle_analytic(
                    sd->ref, sd->x, sd->weight, sd->nat, ref_center, act_center, erg->vec);
            fprintf(fp, "%6d%6d%12.3f", n, sd->nat, fitangle);
        }
    }
    fprintf(fp, "\n");

#undef SLAB_MIN_ATOMS
}


/* Shift x with is */
static inline void shift_single_coord(const matrix box, rvec x, const ivec is)
{
    int tx, ty, tz;


    tx = is[XX];
    ty = is[YY];
    tz = is[ZZ];

    if (TRICLINIC(box))
    {
        x[XX] += tx * box[XX][XX] + ty * box[YY][XX] + tz * box[ZZ][XX];
        x[YY] += ty * box[YY][YY] + tz * box[ZZ][YY];
        x[ZZ] += tz * box[ZZ][ZZ];
    }
    else
    {
        x[XX] += tx * box[XX][XX];
        x[YY] += ty * box[YY][YY];
        x[ZZ] += tz * box[ZZ][ZZ];
    }
}


/* Determine the 'home' slab of this atom which is the
 * slab with the highest Gaussian weight of all */
static inline int get_homeslab(rvec curr_x, /* The position for which the home slab shall be determined */
                               const rvec rotvec, /* The rotation vector */
                               real       slabdist)     /* The slab distance */
{
    real dist;


    /* The distance of the atom to the coordinate center (where the
     * slab with index 0) is */
    dist = iprod(rotvec, curr_x);

    return gmx::roundToInt(dist / slabdist);
}


/* For a local atom determine the relevant slabs, i.e. slabs in
 * which the gaussian is larger than min_gaussian
 */
static int get_single_atom_gaussians(rvec curr_x, gmx_enfrotgrp* erg)
{

    /* Determine the 'home' slab of this atom: */
    int homeslab = get_homeslab(curr_x, erg->vec, erg->rotg->slab_dist);

    /* First determine the weight in the atoms home slab: */
    real g                 = gaussian_weight(curr_x, erg, homeslab);
    int  count             = 0;
    erg->gn_atom[count]    = g;
    erg->gn_slabind[count] = homeslab;
    count++;


    /* Determine the max slab */
    int slab = homeslab;
    while (g > erg->rotg->min_gaussian)
    {
        slab++;
        g                      = gaussian_weight(curr_x, erg, slab);
        erg->gn_slabind[count] = slab;
        erg->gn_atom[count]    = g;
        count++;
    }
    count--;

    /* Determine the min slab */
    slab = homeslab;
    do
    {
        slab--;
        g                      = gaussian_weight(curr_x, erg, slab);
        erg->gn_slabind[count] = slab;
        erg->gn_atom[count]    = g;
        count++;
    } while (g > erg->rotg->min_gaussian);
    count--;

    return count;
}


static void flex2_precalc_inner_sum(const gmx_enfrotgrp* erg)
{
    rvec xi;       /* positions in the i-sum                        */
    rvec xcn, ycn; /* the current and the reference slab centers    */
    real gaussian_xi;
    rvec yi0;
    rvec rin; /* Helper variables                              */
    real fac, fac2;
    rvec innersumvec;
    real OOpsii, OOpsiistar;
    real sin_rin; /* s_ii.r_ii */
    rvec s_in, tmpvec, tmpvec2;
    real mi, wi; /* Mass-weighting of the positions                 */
    real N_M;    /* N/M                                             */


    N_M = erg->rotg->nat * erg->invmass;

    /* Loop over all slabs that contain something */
    for (int n = erg->slab_first; n <= erg->slab_last; n++)
    {
        int slabIndex = n - erg->slab_first; /* slab index */

        /* The current center of this slab is saved in xcn: */
        copy_rvec(erg->slab_center[slabIndex], xcn);
        /* ... and the reference center in ycn: */
        copy_rvec(erg->slab_center_ref[slabIndex + erg->slab_buffer], ycn);

        /*** D. Calculate the whole inner sum used for second and third sum */
        /* For slab n, we need to loop over all atoms i again. Since we sorted
         * the atoms with respect to the rotation vector, we know that it is sufficient
         * to calculate from firstatom to lastatom only. All other contributions will
         * be very small. */
        clear_rvec(innersumvec);
        for (int i = erg->firstatom[slabIndex]; i <= erg->lastatom[slabIndex]; i++)
        {
            /* Coordinate xi of this atom */
            copy_rvec(erg->xc[i], xi);

            /* The i-weights */
            gaussian_xi = gaussian_weight(xi, erg, n);
            mi          = erg->mc_sorted[i]; /* need the sorted mass here */
            wi          = N_M * mi;

            /* Calculate rin */
            copy_rvec(erg->xc_ref_sorted[i], yi0); /* Reference position yi0   */
            rvec_sub(yi0, ycn, tmpvec2);           /* tmpvec2 = yi0 - ycn      */
            mvmul(erg->rotmat, tmpvec2, rin);      /* rin = Omega.(yi0 - ycn)  */

            /* Calculate psi_i* and sin */
            rvec_sub(xi, xcn, tmpvec2); /* tmpvec2 = xi - xcn       */

            /* In rare cases, when an atom position coincides with a slab center
             * (tmpvec2 == 0) we cannot compute the vector product for s_in.
             * However, since the atom is located directly on the pivot, this
             * slab's contribution to the force on that atom will be zero
             * anyway. Therefore, we continue with the next atom. */
            if (gmx_numzero(norm(tmpvec2))) /* 0 == norm(xi - xcn) */
            {
                continue;
            }

            cprod(erg->vec, tmpvec2, tmpvec);            /* tmpvec = v x (xi - xcn)  */
            OOpsiistar = norm2(tmpvec) + erg->rotg->eps; /* OOpsii* = 1/psii* = |v x (xi-xcn)|^2 + eps */
            OOpsii = norm(tmpvec);                       /* OOpsii = 1 / psii = |v x (xi - xcn)| */

            /*                           *         v x (xi - xcn)          */
            unitv(tmpvec, s_in); /*  sin = ----------------         */
                                 /*        |v x (xi - xcn)|         */

            sin_rin = iprod(s_in, rin); /* sin_rin = sin . rin             */

            /* Now the whole sum */
            fac = OOpsii / OOpsiistar;
            svmul(fac, rin, tmpvec);
            fac2 = fac * fac * OOpsii;
            svmul(fac2 * sin_rin, s_in, tmpvec2);
            rvec_dec(tmpvec, tmpvec2);

            svmul(wi * gaussian_xi * sin_rin, tmpvec, tmpvec2);

            rvec_inc(innersumvec, tmpvec2);
        } /* now we have the inner sum, used both for sum2 and sum3 */

        /* Save it to be used in do_flex2_lowlevel */
        copy_rvec(innersumvec, erg->slab_innersumvec[slabIndex]);
    } /* END of loop over slabs */
}


static void flex_precalc_inner_sum(const gmx_enfrotgrp* erg)
{
    rvec xi;       /* position                                      */
    rvec xcn, ycn; /* the current and the reference slab centers    */
    rvec qin, rin; /* q_i^n and r_i^n                               */
    real bin;
    rvec tmpvec;
    rvec innersumvec; /* Inner part of sum_n2                          */
    real gaussian_xi; /* Gaussian weight gn(xi)                        */
    real mi, wi;      /* Mass-weighting of the positions               */
    real N_M;         /* N/M                                           */

    N_M = erg->rotg->nat * erg->invmass;

    /* Loop over all slabs that contain something */
    for (int n = erg->slab_first; n <= erg->slab_last; n++)
    {
        int slabIndex = n - erg->slab_first; /* slab index */

        /* The current center of this slab is saved in xcn: */
        copy_rvec(erg->slab_center[slabIndex], xcn);
        /* ... and the reference center in ycn: */
        copy_rvec(erg->slab_center_ref[slabIndex + erg->slab_buffer], ycn);

        /* For slab n, we need to loop over all atoms i again. Since we sorted
         * the atoms with respect to the rotation vector, we know that it is sufficient
         * to calculate from firstatom to lastatom only. All other contributions will
         * be very small. */
        clear_rvec(innersumvec);
        for (int i = erg->firstatom[slabIndex]; i <= erg->lastatom[slabIndex]; i++)
        {
            /* Coordinate xi of this atom */
            copy_rvec(erg->xc[i], xi);

            /* The i-weights */
            gaussian_xi = gaussian_weight(xi, erg, n);
            mi          = erg->mc_sorted[i]; /* need the sorted mass here */
            wi          = N_M * mi;

            /* Calculate rin and qin */
            rvec_sub(erg->xc_ref_sorted[i], ycn, tmpvec); /* tmpvec = yi0-ycn */

            /* In rare cases, when an atom position coincides with a slab center
             * (tmpvec == 0) we cannot compute the vector product for qin.
             * However, since the atom is located directly on the pivot, this
             * slab's contribution to the force on that atom will be zero
             * anyway. Therefore, we continue with the next atom. */
            if (gmx_numzero(norm(tmpvec))) /* 0 == norm(yi0 - ycn) */
            {
                continue;
            }

            mvmul(erg->rotmat, tmpvec, rin); /* rin = Omega.(yi0 - ycn)  */
            cprod(erg->vec, rin, tmpvec);    /* tmpvec = v x Omega*(yi0-ycn) */

            /*                                *        v x Omega*(yi0-ycn)    */
            unitv(tmpvec, qin); /* qin = ---------------------   */
                                /*       |v x Omega*(yi0-ycn)|   */

            /* Calculate bin */
            rvec_sub(xi, xcn, tmpvec); /* tmpvec = xi-xcn          */
            bin = iprod(qin, tmpvec);  /* bin  = qin*(xi-xcn)      */

            svmul(wi * gaussian_xi * bin, qin, tmpvec);

            /* Add this contribution to the inner sum: */
            rvec_add(innersumvec, tmpvec, innersumvec);
        } /* now we have the inner sum vector S^n for this slab */
          /* Save it to be used in do_flex_lowlevel */
        copy_rvec(innersumvec, erg->slab_innersumvec[slabIndex]);
    }
}


static real do_flex2_lowlevel(gmx_enfrotgrp*                 erg,
                              real                           sigma, /* The Gaussian width sigma */
                              gmx::ArrayRef<const gmx::RVec> coords,
                              gmx_bool                       bOutstepRot,
                              gmx_bool                       bOutstepSlab,
                              const matrix                   box)
{
    int  count, ii, iigrp;
    rvec xj;          /* position in the i-sum                         */
    rvec yj0;         /* the reference position in the j-sum           */
    rvec xcn, ycn;    /* the current and the reference slab centers    */
    real V;           /* This node's part of the rotation pot. energy  */
    real gaussian_xj; /* Gaussian weight                               */
    real beta;

    real numerator, fit_numerator;
    rvec rjn, fit_rjn; /* Helper variables                              */
    real fac, fac2;

    real     OOpsij, OOpsijstar;
    real     OOsigma2; /* 1/(sigma^2)                                   */
    real     sjn_rjn;
    real     betasigpsi;
    rvec     sjn, tmpvec, tmpvec2, yj0_ycn;
    rvec     sum1vec_part, sum1vec, sum2vec_part, sum2vec, sum3vec, sum4vec, innersumvec;
    real     sum3, sum4;
    real     mj, wj; /* Mass-weighting of the positions               */
    real     N_M;    /* N/M                                           */
    real     Wjn;    /* g_n(x_j) m_j / Mjn                            */
    gmx_bool bCalcPotFit;

    /* To calculate the torque per slab */
    rvec slab_force; /* Single force from slab n on one atom          */
    rvec slab_sum1vec_part;
    real slab_sum3part, slab_sum4part;
    rvec slab_sum1vec, slab_sum2vec, slab_sum3vec, slab_sum4vec;

    /* Pre-calculate the inner sums, so that we do not have to calculate
     * them again for every atom */
    flex2_precalc_inner_sum(erg);

    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (RotationGroupFitting::Pot == erg->rotg->eFittype);

    /********************************************************/
    /* Main loop over all local atoms of the rotation group */
    /********************************************************/
    N_M                                      = erg->rotg->nat * erg->invmass;
    V                                        = 0.0;
    OOsigma2                                 = 1.0 / (sigma * sigma);
    const auto& localRotationGroupIndex      = erg->atomSet->localIndex();
    const auto& collectiveRotationGroupIndex = erg->atomSet->collectiveIndex();

    for (gmx::Index j = 0; j < localRotationGroupIndex.ssize(); j++)
    {
        /* Local index of a rotation group atom  */
        ii = localRotationGroupIndex[j];
        /* Position of this atom in the collective array */
        iigrp = collectiveRotationGroupIndex[j];
        /* Mass-weighting */
        mj = erg->mc[iigrp]; /* need the unsorted mass here */
        wj = N_M * mj;

        /* Current position of this atom: x[ii][XX/YY/ZZ]
         * Note that erg->xc_center contains the center of mass in case the flex2-t
         * potential was chosen. For the flex2 potential erg->xc_center must be
         * zero. */
        rvec_sub(coords[ii], erg->xc_center, xj);

        /* Shift this atom such that it is near its reference */
        shift_single_coord(box, xj, erg->xc_shifts[iigrp]);

        /* Determine the slabs to loop over, i.e. the ones with contributions
         * larger than min_gaussian */
        count = get_single_atom_gaussians(xj, erg);

        clear_rvec(sum1vec_part);
        clear_rvec(sum2vec_part);
        sum3 = 0.0;
        sum4 = 0.0;
        /* Loop over the relevant slabs for this atom */
        for (int ic = 0; ic < count; ic++)
        {
            int n = erg->gn_slabind[ic];

            /* Get the precomputed Gaussian value of curr_slab for curr_x */
            gaussian_xj = erg->gn_atom[ic];

            int slabIndex = n - erg->slab_first; /* slab index */

            /* The (unrotated) reference position of this atom is copied to yj0: */
            copy_rvec(erg->referencePositions[iigrp], yj0);

            beta = calc_beta(xj, erg, n);

            /* The current center of this slab is saved in xcn: */
            copy_rvec(erg->slab_center[slabIndex], xcn);
            /* ... and the reference center in ycn: */
            copy_rvec(erg->slab_center_ref[slabIndex + erg->slab_buffer], ycn);

            rvec_sub(yj0, ycn, yj0_ycn); /* yj0_ycn = yj0 - ycn      */

            /* Rotate: */
            mvmul(erg->rotmat, yj0_ycn, rjn); /* rjn = Omega.(yj0 - ycn)  */

            /* Subtract the slab center from xj */
            rvec_sub(xj, xcn, tmpvec2); /* tmpvec2 = xj - xcn       */

            /* In rare cases, when an atom position coincides with a slab center
             * (tmpvec2 == 0) we cannot compute the vector product for sjn.
             * However, since the atom is located directly on the pivot, this
             * slab's contribution to the force on that atom will be zero
             * anyway. Therefore, we directly move on to the next slab.       */
            if (gmx_numzero(norm(tmpvec2))) /* 0 == norm(xj - xcn) */
            {
                continue;
            }

            /* Calculate sjn */
            cprod(erg->vec, tmpvec2, tmpvec); /* tmpvec = v x (xj - xcn)  */

            OOpsijstar = norm2(tmpvec) + erg->rotg->eps; /* OOpsij* = 1/psij* = |v x (xj-xcn)|^2 + eps */

            numerator = gmx::square(iprod(tmpvec, rjn));

            /*********************************/
            /* Add to the rotation potential */
            /*********************************/
            V += 0.5 * erg->rotg->k * wj * gaussian_xj * numerator / OOpsijstar;

            /* If requested, also calculate the potential for a set of angles
             * near the current reference angle */
            if (bCalcPotFit)
            {
                for (int ifit = 0; ifit < erg->rotg->PotAngle_nstep; ifit++)
                {
                    mvmul(erg->PotAngleFit->rotmat[ifit], yj0_ycn, fit_rjn);
                    fit_numerator = gmx::square(iprod(tmpvec, fit_rjn));
                    erg->PotAngleFit->V[ifit] +=
                            0.5 * erg->rotg->k * wj * gaussian_xj * fit_numerator / OOpsijstar;
                }
            }

            /*************************************/
            /* Now calculate the force on atom j */
            /*************************************/

            OOpsij = norm(tmpvec); /* OOpsij = 1 / psij = |v x (xj - xcn)| */

            /*                              *         v x (xj - xcn)          */
            unitv(tmpvec, sjn); /*  sjn = ----------------         */
                                /*        |v x (xj - xcn)|         */

            sjn_rjn = iprod(sjn, rjn); /* sjn_rjn = sjn . rjn             */


            /*** A. Calculate the first of the four sum terms: ****************/
            fac = OOpsij / OOpsijstar;
            svmul(fac, rjn, tmpvec);
            fac2 = fac * fac * OOpsij;
            svmul(fac2 * sjn_rjn, sjn, tmpvec2);
            rvec_dec(tmpvec, tmpvec2);
            fac2 = wj * gaussian_xj; /* also needed for sum4 */
            svmul(fac2 * sjn_rjn, tmpvec, slab_sum1vec_part);
            /********************/
            /*** Add to sum1: ***/
            /********************/
            rvec_inc(sum1vec_part, slab_sum1vec_part); /* sum1 still needs to vector multiplied with v */

            /*** B. Calculate the forth of the four sum terms: ****************/
            betasigpsi = beta * OOsigma2 * OOpsij; /* this is also needed for sum3 */
            /********************/
            /*** Add to sum4: ***/
            /********************/
            slab_sum4part = fac2 * betasigpsi * fac * sjn_rjn
                            * sjn_rjn; /* Note that fac is still valid from above */
            sum4 += slab_sum4part;

            /*** C. Calculate Wjn for second and third sum */
            /* Note that we can safely divide by slab_weights since we check in
             * get_slab_centers that it is non-zero. */
            Wjn = gaussian_xj * mj / erg->slab_weights[slabIndex];

            /* We already have precalculated the inner sum for slab n */
            copy_rvec(erg->slab_innersumvec[slabIndex], innersumvec);

            /* Weigh the inner sum vector with Wjn */
            svmul(Wjn, innersumvec, innersumvec);

            /*** E. Calculate the second of the four sum terms: */
            /********************/
            /*** Add to sum2: ***/
            /********************/
            rvec_inc(sum2vec_part, innersumvec); /* sum2 still needs to be vector crossproduct'ed with v */

            /*** F. Calculate the third of the four sum terms: */
            slab_sum3part = betasigpsi * iprod(sjn, innersumvec);
            sum3 += slab_sum3part; /* still needs to be multiplied with v */

            /*** G. Calculate the torque on the local slab's axis: */
            if (bOutstepRot)
            {
                /* Sum1 */
                cprod(slab_sum1vec_part, erg->vec, slab_sum1vec);
                /* Sum2 */
                cprod(innersumvec, erg->vec, slab_sum2vec);
                /* Sum3 */
                svmul(slab_sum3part, erg->vec, slab_sum3vec);
                /* Sum4 */
                svmul(slab_sum4part, erg->vec, slab_sum4vec);

                /* The force on atom ii from slab n only: */
                for (int m = 0; m < DIM; m++)
                {
                    slab_force[m] = erg->rotg->k
                                    * (-slab_sum1vec[m] + slab_sum2vec[m] - slab_sum3vec[m]
                                       + 0.5 * slab_sum4vec[m]);
                }

                erg->slab_torque_v[slabIndex] += torque(erg->vec, slab_force, xj, xcn);
            }
        } /* END of loop over slabs */

        /* Construct the four individual parts of the vector sum: */
        cprod(sum1vec_part, erg->vec, sum1vec); /* sum1vec =   { } x v  */
        cprod(sum2vec_part, erg->vec, sum2vec); /* sum2vec =   { } x v  */
        svmul(sum3, erg->vec, sum3vec);         /* sum3vec =   { } . v  */
        svmul(sum4, erg->vec, sum4vec);         /* sum4vec =   { } . v  */

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        for (int m = 0; m < DIM; m++)
        {
            erg->f_rot_loc[j][m] =
                    erg->rotg->k * (-sum1vec[m] + sum2vec[m] - sum3vec[m] + 0.5 * sum4vec[m]);
        }

#ifdef SUM_PARTS
        fprintf(stderr,
                "sum1: %15.8f %15.8f %15.8f\n",
                -erg->rotg->k * sum1vec[XX],
                -erg->rotg->k * sum1vec[YY],
                -erg->rotg->k * sum1vec[ZZ]);
        fprintf(stderr,
                "sum2: %15.8f %15.8f %15.8f\n",
                erg->rotg->k * sum2vec[XX],
                erg->rotg->k * sum2vec[YY],
                erg->rotg->k * sum2vec[ZZ]);
        fprintf(stderr,
                "sum3: %15.8f %15.8f %15.8f\n",
                -erg->rotg->k * sum3vec[XX],
                -erg->rotg->k * sum3vec[YY],
                -erg->rotg->k * sum3vec[ZZ]);
        fprintf(stderr,
                "sum4: %15.8f %15.8f %15.8f\n",
                0.5 * erg->rotg->k * sum4vec[XX],
                0.5 * erg->rotg->k * sum4vec[YY],
                0.5 * erg->rotg->k * sum4vec[ZZ]);
#endif

        PRINT_FORCE_J

    } /* END of loop over local atoms */

    return V;
}


static real do_flex_lowlevel(gmx_enfrotgrp* erg,
                             real sigma, /* The Gaussian width sigma                      */
                             gmx::ArrayRef<const gmx::RVec> coords,
                             gmx_bool                       bOutstepRot,
                             gmx_bool                       bOutstepSlab,
                             const matrix                   box)
{
    int      count, iigrp;
    rvec     xj, yj0;        /* current and reference position                */
    rvec     xcn, ycn;       /* the current and the reference slab centers    */
    rvec     yj0_ycn;        /* yj0 - ycn                                     */
    rvec     xj_xcn;         /* xj - xcn                                      */
    rvec     qjn, fit_qjn;   /* q_i^n                                         */
    rvec     sum_n1, sum_n2; /* Two contributions to the rotation force       */
    rvec     innersumvec;    /* Inner part of sum_n2                          */
    rvec     s_n;
    rvec     force_n;                /* Single force from slab n on one atom          */
    rvec     force_n1, force_n2;     /* First and second part of force_n              */
    rvec     tmpvec, tmpvec2, tmp_f; /* Helper variables                              */
    real     V;                      /* The rotation potential energy                 */
    real     OOsigma2;               /* 1/(sigma^2)                                   */
    real     beta;                   /* beta_n(xj)                                    */
    real     bjn, fit_bjn;           /* b_j^n                                         */
    real     gaussian_xj;            /* Gaussian weight gn(xj)                        */
    real     betan_xj_sigma2;
    real     mj, wj; /* Mass-weighting of the positions               */
    real     N_M;    /* N/M                                           */
    gmx_bool bCalcPotFit;

    /* Pre-calculate the inner sums, so that we do not have to calculate
     * them again for every atom */
    flex_precalc_inner_sum(erg);

    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (RotationGroupFitting::Pot == erg->rotg->eFittype);

    /********************************************************/
    /* Main loop over all local atoms of the rotation group */
    /********************************************************/
    OOsigma2                                 = 1.0 / (sigma * sigma);
    N_M                                      = erg->rotg->nat * erg->invmass;
    V                                        = 0.0;
    const auto& localRotationGroupIndex      = erg->atomSet->localIndex();
    const auto& collectiveRotationGroupIndex = erg->atomSet->collectiveIndex();

    for (gmx::Index j = 0; j < localRotationGroupIndex.ssize(); j++)
    {
        /* Local index of a rotation group atom  */
        int ii = localRotationGroupIndex[j];
        /* Position of this atom in the collective array */
        iigrp = collectiveRotationGroupIndex[j];
        /* Mass-weighting */
        mj = erg->mc[iigrp]; /* need the unsorted mass here */
        wj = N_M * mj;

        /* Current position of this atom: x[ii][XX/YY/ZZ]
         * Note that erg->xc_center contains the center of mass in case the flex-t
         * potential was chosen. For the flex potential erg->xc_center must be
         * zero. */
        rvec_sub(coords[ii], erg->xc_center, xj);

        /* Shift this atom such that it is near its reference */
        shift_single_coord(box, xj, erg->xc_shifts[iigrp]);

        /* Determine the slabs to loop over, i.e. the ones with contributions
         * larger than min_gaussian */
        count = get_single_atom_gaussians(xj, erg);

        clear_rvec(sum_n1);
        clear_rvec(sum_n2);

        /* Loop over the relevant slabs for this atom */
        for (int ic = 0; ic < count; ic++)
        {
            int n = erg->gn_slabind[ic];

            /* Get the precomputed Gaussian for xj in slab n */
            gaussian_xj = erg->gn_atom[ic];

            int slabIndex = n - erg->slab_first; /* slab index */

            /* The (unrotated) reference position of this atom is saved in yj0: */
            copy_rvec(erg->referencePositions[iigrp], yj0);

            beta = calc_beta(xj, erg, n);

            /* The current center of this slab is saved in xcn: */
            copy_rvec(erg->slab_center[slabIndex], xcn);
            /* ... and the reference center in ycn: */
            copy_rvec(erg->slab_center_ref[slabIndex + erg->slab_buffer], ycn);

            rvec_sub(yj0, ycn, yj0_ycn); /* yj0_ycn = yj0 - ycn */

            /* In rare cases, when an atom position coincides with a reference slab
             * center (yj0_ycn == 0) we cannot compute the normal vector qjn.
             * However, since the atom is located directly on the pivot, this
             * slab's contribution to the force on that atom will be zero
             * anyway. Therefore, we directly move on to the next slab.       */
            if (gmx_numzero(norm(yj0_ycn))) /* 0 == norm(yj0 - ycn) */
            {
                continue;
            }

            /* Rotate: */
            mvmul(erg->rotmat, yj0_ycn, tmpvec2); /* tmpvec2= Omega.(yj0-ycn) */

            /* Subtract the slab center from xj */
            rvec_sub(xj, xcn, xj_xcn); /* xj_xcn = xj - xcn         */

            /* Calculate qjn */
            cprod(erg->vec, tmpvec2, tmpvec); /* tmpvec= v x Omega.(yj0-ycn) */

            /*                         *         v x Omega.(yj0-ycn)    */
            unitv(tmpvec, qjn); /*  qjn = ---------------------   */
                                /*        |v x Omega.(yj0-ycn)|   */

            bjn = iprod(qjn, xj_xcn); /* bjn = qjn * (xj - xcn) */

            /*********************************/
            /* Add to the rotation potential */
            /*********************************/
            V += 0.5 * erg->rotg->k * wj * gaussian_xj * gmx::square(bjn);

            /* If requested, also calculate the potential for a set of angles
             * near the current reference angle */
            if (bCalcPotFit)
            {
                for (int ifit = 0; ifit < erg->rotg->PotAngle_nstep; ifit++)
                {
                    /* As above calculate Omega.(yj0-ycn), now for the other angles */
                    mvmul(erg->PotAngleFit->rotmat[ifit], yj0_ycn, tmpvec2); /* tmpvec2= Omega.(yj0-ycn) */
                    /* As above calculate qjn */
                    cprod(erg->vec, tmpvec2, tmpvec); /* tmpvec= v x Omega.(yj0-ycn) */
                    /*                                                        *             v x Omega.(yj0-ycn) */
                    unitv(tmpvec, fit_qjn);           /*  fit_qjn = ---------------------   */
                                                      /*            |v x Omega.(yj0-ycn)|   */
                    fit_bjn = iprod(fit_qjn, xj_xcn); /* fit_bjn = fit_qjn * (xj - xcn) */
                    /* Add to the rotation potential for this angle */
                    erg->PotAngleFit->V[ifit] +=
                            0.5 * erg->rotg->k * wj * gaussian_xj * gmx::square(fit_bjn);
                }
            }

            /****************************************************************/
            /* sum_n1 will typically be the main contribution to the force: */
            /****************************************************************/
            betan_xj_sigma2 = beta * OOsigma2; /*  beta_n(xj)/sigma^2  */

            /* The next lines calculate
             *  qjn - (bjn*beta(xj)/(2sigma^2))v  */
            svmul(bjn * 0.5 * betan_xj_sigma2, erg->vec, tmpvec2);
            rvec_sub(qjn, tmpvec2, tmpvec);

            /* Multiply with gn(xj)*bjn: */
            svmul(gaussian_xj * bjn, tmpvec, tmpvec2);

            /* Sum over n: */
            rvec_inc(sum_n1, tmpvec2);

            /* We already have precalculated the Sn term for slab n */
            copy_rvec(erg->slab_innersumvec[slabIndex], s_n);
            /*                                                             *          beta_n(xj) */
            svmul(betan_xj_sigma2 * iprod(s_n, xj_xcn), erg->vec, tmpvec); /* tmpvec = ---------- s_n (xj-xcn) */
            /*            sigma^2               */

            rvec_sub(s_n, tmpvec, innersumvec);

            /* We can safely divide by slab_weights since we check in get_slab_centers
             * that it is non-zero. */
            svmul(gaussian_xj / erg->slab_weights[slabIndex], innersumvec, innersumvec);

            rvec_add(sum_n2, innersumvec, sum_n2);

            /* Calculate the torque: */
            if (bOutstepRot)
            {
                /* The force on atom ii from slab n only: */
                svmul(-erg->rotg->k * wj, tmpvec2, force_n1);    /* part 1 */
                svmul(erg->rotg->k * mj, innersumvec, force_n2); /* part 2 */
                rvec_add(force_n1, force_n2, force_n);
                erg->slab_torque_v[slabIndex] += torque(erg->vec, force_n, xj, xcn);
            }
        } /* END of loop over slabs */

        /* Put both contributions together: */
        svmul(wj, sum_n1, sum_n1);
        svmul(mj, sum_n2, sum_n2);
        rvec_sub(sum_n2, sum_n1, tmp_f); /* F = -grad V */

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        for (int m = 0; m < DIM; m++)
        {
            erg->f_rot_loc[j][m] = erg->rotg->k * tmp_f[m];
        }

        PRINT_FORCE_J

    } /* END of loop over local atoms */

    return V;
}

static void sort_collective_coordinates(gmx_enfrotgrp*    erg,
                                        sort_along_vec_t* data) /* Buffer for sorting the positions */
{
    /* The projection of the position vector on the rotation vector is
     * the relevant value for sorting. Fill the 'data' structure */
    for (int i = 0; i < erg->rotg->nat; i++)
    {
        data[i].xcproj = iprod(erg->xc[i], erg->vec); /* sort criterium */
        data[i].m      = erg->mc[i];
        data[i].ind    = i;
        copy_rvec(erg->xc[i], data[i].x);
        copy_rvec(erg->referencePositions[i], data[i].x_ref);
    }
    /* Sort the 'data' structure */
    std::sort(data, data + erg->rotg->nat, [](const sort_along_vec_t& a, const sort_along_vec_t& b) {
        return a.xcproj < b.xcproj;
    });

    /* Copy back the sorted values */
    for (int i = 0; i < erg->rotg->nat; i++)
    {
        copy_rvec(data[i].x, erg->xc[i]);
        copy_rvec(data[i].x_ref, erg->xc_ref_sorted[i]);
        erg->mc_sorted[i]  = data[i].m;
        erg->xc_sortind[i] = data[i].ind;
    }
}


/* For each slab, get the first and the last index of the sorted atom
 * indices */
static void get_firstlast_atom_per_slab(const gmx_enfrotgrp* erg)
{
    real beta;

    /* Find the first atom that needs to enter the calculation for each slab */
    int n = erg->slab_first; /* slab */
    int i = 0;               /* start with the first atom */
    do
    {
        /* Find the first atom that significantly contributes to this slab */
        do /* move forward in position until a large enough beta is found */
        {
            beta = calc_beta(erg->xc[i], erg, n);
            i++;
        } while ((beta < -erg->max_beta) && (i < erg->rotg->nat));
        i--;
        int slabIndex             = n - erg->slab_first; /* slab index */
        erg->firstatom[slabIndex] = i;
        /* Proceed to the next slab */
        n++;
    } while (n <= erg->slab_last);

    /* Find the last atom for each slab */
    n = erg->slab_last;     /* start with last slab */
    i = erg->rotg->nat - 1; /* start with the last atom */
    do
    {
        do /* move backward in position until a large enough beta is found */
        {
            beta = calc_beta(erg->xc[i], erg, n);
            i--;
        } while ((beta > erg->max_beta) && (i > -1));
        i++;
        int slabIndex            = n - erg->slab_first; /* slab index */
        erg->lastatom[slabIndex] = i;
        /* Proceed to the next slab */
        n--;
    } while (n >= erg->slab_first);
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
static inline int get_first_slab(const gmx_enfrotgrp* erg,
                                 const gmx::RVec& firstatom) /* First atom after sorting along the rotation vector v */
{
    /* Find the first slab for the first atom */
    return static_cast<int>(std::ceil(
            static_cast<double>((iprod(firstatom, erg->vec) - erg->max_beta) / erg->rotg->slab_dist)));
}


static inline int get_last_slab(const gmx_enfrotgrp* erg, const gmx::RVec& lastatom) /* Last atom along v */
{
    /* Find the last slab for the last atom */
    return static_cast<int>(std::floor(
            static_cast<double>((iprod(lastatom, erg->vec) + erg->max_beta) / erg->rotg->slab_dist)));
}


static void get_firstlast_slab_check(
        gmx_enfrotgrp*   erg,       /* The rotation group (data only accessible in this file) */
        const gmx::RVec& firstatom, /* First atom after sorting along the rotation vector v */
        const gmx::RVec& lastatom)  /* Last atom along v */
{
    erg->slab_first = get_first_slab(erg, firstatom);
    erg->slab_last  = get_last_slab(erg, lastatom);

    /* Calculate the slab buffer size, which changes when slab_first changes */
    erg->slab_buffer = erg->slab_first - erg->slab_first_ref;

    /* Check whether we have reference data to compare against */
    if (erg->slab_first < erg->slab_first_ref)
    {
        gmx_fatal(FARGS,
                  "%s No reference data for first slab (n=%d), unable to proceed.",
                  RotStr.c_str(),
                  erg->slab_first);
    }

    /* Check whether we have reference data to compare against */
    if (erg->slab_last > erg->slab_last_ref)
    {
        gmx_fatal(FARGS,
                  "%s No reference data for last slab (n=%d), unable to proceed.",
                  RotStr.c_str(),
                  erg->slab_last);
    }
}


/* Enforced rotation with a flexible axis */
static void do_flexible(gmx_bool       bMain,
                        gmx_enfrot*    enfrot, /* Other rotation data                        */
                        gmx_enfrotgrp* erg,
                        gmx::ArrayRef<const gmx::RVec> coords, /* The local positions */
                        const matrix                   box,
                        double   t,            /* Time in picoseconds                        */
                        gmx_bool bOutstepRot,  /* Output to main rotation output file        */
                        gmx_bool bOutstepSlab) /* Output per-slab data                       */
{
    int  nslabs;
    real sigma; /* The Gaussian width sigma */

    /* Define the sigma value */
    sigma = 0.7 * erg->rotg->slab_dist;

    /* Sort the collective coordinates erg->xc along the rotation vector. This is
     * an optimization for the inner loop. */
    sort_collective_coordinates(erg, enfrot->data);

    /* Determine the first relevant slab for the first atom and the last
     * relevant slab for the last atom */
    get_firstlast_slab_check(erg, erg->xc[0], erg->xc[erg->rotg->nat - 1]);

    /* Determine for each slab depending on the min_gaussian cutoff criterium,
     * a first and a last atom index inbetween stuff needs to be calculated */
    get_firstlast_atom_per_slab(erg);

    /* Determine the gaussian-weighted center of positions for all slabs */
    get_slab_centers(erg,
                     gmx::arrayRefFromArray(reinterpret_cast<gmx::RVec*>(erg->xc), erg->rotg->nat),
                     erg->mc_sorted,
                     t,
                     enfrot->out_slabs,
                     bOutstepSlab,
                     FALSE);

    /* Clear the torque per slab from last time step: */
    nslabs = erg->slab_last - erg->slab_first + 1;
    for (int l = 0; l < nslabs; l++)
    {
        erg->slab_torque_v[l] = 0.0;
    }

    /* Call the rotational forces kernel */
    if (erg->rotg->eType == EnforcedRotationGroupType::Flex
        || erg->rotg->eType == EnforcedRotationGroupType::Flext)
    {
        erg->V = do_flex_lowlevel(erg, sigma, coords, bOutstepRot, bOutstepSlab, box);
    }
    else if (erg->rotg->eType == EnforcedRotationGroupType::Flex2
             || erg->rotg->eType == EnforcedRotationGroupType::Flex2t)
    {
        erg->V = do_flex2_lowlevel(erg, sigma, coords, bOutstepRot, bOutstepSlab, box);
    }
    else
    {
        gmx_fatal(FARGS, "Unknown flexible rotation type");
    }

    /* Determine angle by RMSD fit to the reference - Let's hope this */
    /* only happens once in a while, since this is not parallelized! */
    if (bMain && (RotationGroupFitting::Pot != erg->rotg->eFittype))
    {
        if (bOutstepRot)
        {
            /* Fit angle of the whole rotation group */
            erg->angle_v = flex_fit_angle(erg);
        }
        if (bOutstepSlab)
        {
            /* Fit angle of each slab */
            flex_fit_angle_perslab(erg, t, erg->degangle, enfrot->out_angles);
        }
    }

    /* Lump together the torques from all slabs: */
    erg->torque_v = 0.0;
    for (int l = 0; l < nslabs; l++)
    {
        erg->torque_v += erg->slab_torque_v[l];
    }
}


/* Calculate the angle between reference and actual rotation group atom,
 * both projected into a plane perpendicular to the rotation vector: */
static void angle(const gmx_enfrotgrp* erg,
                  rvec                 x_act,
                  rvec                 x_ref,
                  real*                alpha,
                  real* weight) /* atoms near the rotation axis should count less than atoms far away */
{
    rvec xp, xrp; /* current and reference positions projected on a plane perpendicular to pg->vec */
    rvec dum;


    /* Project x_ref and x into a plane through the origin perpendicular to rot_vec: */
    /* Project x_ref: xrp = x_ref - (vec * x_ref) * vec */
    svmul(iprod(erg->vec, x_ref), erg->vec, dum);
    rvec_sub(x_ref, dum, xrp);
    /* Project x_act: */
    svmul(iprod(erg->vec, x_act), erg->vec, dum);
    rvec_sub(x_act, dum, xp);

    /* Retrieve information about which vector precedes. gmx_angle always
     * returns a positive angle. */
    cprod(xp, xrp, dum); /* if reference precedes, this is pointing into the same direction as vec */

    if (iprod(erg->vec, dum) >= 0)
    {
        *alpha = -gmx_angle(xrp, xp);
    }
    else
    {
        *alpha = +gmx_angle(xrp, xp);
    }

    /* Also return the weight */
    *weight = norm(xp);
}


/* Project first vector onto a plane perpendicular to the second vector
 * dr = dr - (dr.v)v
 * Note that v must be of unit length.
 */
static inline void project_onto_plane(rvec dr, const rvec v)
{
    rvec tmp;


    svmul(iprod(dr, v), v, tmp); /* tmp = (dr.v)v */
    rvec_dec(dr, tmp);           /* dr = dr - (dr.v)v */
}


/* Fixed rotation: The rotation reference group rotates around the v axis. */
/* The atoms of the actual rotation group are attached with imaginary  */
/* springs to the reference atoms.                                     */
static void do_fixed(gmx_enfrotgrp* erg,
                     gmx_bool       bOutstepRot, /* Output to main rotation output file        */
                     gmx_bool       bOutstepSlab)      /* Output per-slab data                       */
{
    rvec     dr;
    rvec     tmp_f;  /* Force */
    real     alpha;  /* a single angle between an actual and a reference position */
    real     weight; /* single weight for a single angle */
    rvec     xi_xc;  /* xi - xc */
    gmx_bool bCalcPotFit;
    rvec     fit_xr_loc;

    /* for mass weighting: */
    real wi;   /* Mass-weighting of the positions */
    real N_M;  /* N/M */
    real k_wi; /* k times wi */

    gmx_bool bProject;

    bProject = (erg->rotg->eType == EnforcedRotationGroupType::Pm)
               || (erg->rotg->eType == EnforcedRotationGroupType::Pmpf);
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (RotationGroupFitting::Pot == erg->rotg->eFittype);

    N_M                                      = erg->rotg->nat * erg->invmass;
    const auto& collectiveRotationGroupIndex = erg->atomSet->collectiveIndex();
    /* Each process calculates the forces on its local atoms */
    for (size_t j = 0; j < erg->atomSet->numAtomsLocal(); j++)
    {
        /* Calculate (x_i-x_c) resp. (x_i-u) */
        rvec_sub(erg->x_loc_pbc[j], erg->xc_center, xi_xc);

        /* Calculate Omega*(y_i-y_c)-(x_i-x_c) */
        rvec_sub(erg->xr_loc[j], xi_xc, dr);

        if (bProject)
        {
            project_onto_plane(dr, erg->vec);
        }

        /* Mass-weighting */
        wi = N_M * erg->m_loc[j];

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        k_wi = erg->rotg->k * wi;
        for (int m = 0; m < DIM; m++)
        {
            tmp_f[m]             = k_wi * dr[m];
            erg->f_rot_loc[j][m] = tmp_f[m];
            erg->V += 0.5 * k_wi * gmx::square(dr[m]);
        }

        /* If requested, also calculate the potential for a set of angles
         * near the current reference angle */
        if (bCalcPotFit)
        {
            for (int ifit = 0; ifit < erg->rotg->PotAngle_nstep; ifit++)
            {
                /* Index of this rotation group atom with respect to the whole rotation group */
                int jj = collectiveRotationGroupIndex[j];

                /* Rotate with the alternative angle. Like rotate_local_reference(),
                 * just for a single local atom */
                mvmul(erg->PotAngleFit->rotmat[ifit],
                      erg->referencePositions[jj],
                      fit_xr_loc); /* fit_xr_loc = Omega*(y_i-y_c) */

                /* Calculate Omega*(y_i-y_c)-(x_i-x_c) */
                rvec_sub(fit_xr_loc, xi_xc, dr);

                if (bProject)
                {
                    project_onto_plane(dr, erg->vec);
                }

                /* Add to the rotation potential for this angle: */
                erg->PotAngleFit->V[ifit] += 0.5 * k_wi * norm2(dr);
            }
        }

        if (bOutstepRot)
        {
            /* Add to the torque of this rotation group */
            erg->torque_v += torque(erg->vec, tmp_f, erg->x_loc_pbc[j], erg->xc_center);

            /* Calculate the angle between reference and actual rotation group atom. */
            angle(erg, xi_xc, erg->xr_loc[j], &alpha, &weight); /* angle in rad, weighted */
            erg->angle_v += alpha * weight;
            erg->weight_v += weight;
        }
        /* If you want enforced rotation to contribute to the virial,
         * activate the following lines:
            if (MAIN(cr))
            {
               Add the rotation contribution to the virial
              for(j=0; j<DIM; j++)
                for(m=0;m<DIM;m++)
                  vir[j][m] += 0.5*f[ii][j]*dr[m];
            }
         */

        PRINT_FORCE_J

    } /* end of loop over local rotation group atoms */
}


/* Calculate the radial motion potential and forces */
static void do_radial_motion(gmx_enfrotgrp* erg,
                             gmx_bool bOutstepRot,  /* Output to main rotation output file        */
                             gmx_bool bOutstepSlab) /* Output per-slab data                       */
{
    rvec     tmp_f;  /* Force */
    real     alpha;  /* a single angle between an actual and a reference position */
    real     weight; /* single weight for a single angle */
    rvec     xj_u;   /* xj - u */
    rvec     tmpvec, fit_tmpvec;
    real     fac, fac2, sum = 0.0;
    rvec     pj;
    gmx_bool bCalcPotFit;

    /* For mass weighting: */
    real wj;  /* Mass-weighting of the positions */
    real N_M; /* N/M */

    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (RotationGroupFitting::Pot == erg->rotg->eFittype);

    N_M                                      = erg->rotg->nat * erg->invmass;
    const auto& collectiveRotationGroupIndex = erg->atomSet->collectiveIndex();
    /* Each process calculates the forces on its local atoms */
    for (size_t j = 0; j < erg->atomSet->numAtomsLocal(); j++)
    {
        /* Calculate (xj-u) */
        rvec_sub(erg->x_loc_pbc[j], erg->xc_center, xj_u); /* xj_u = xj-u */

        /* Calculate Omega.(yj0-u) */
        cprod(erg->vec, erg->xr_loc[j], tmpvec); /* tmpvec = v x Omega.(yj0-u) */

        /*                       *         v x Omega.(yj0-u)     */
        unitv(tmpvec, pj); /*  pj = ---------------------   */
                           /*       | v x Omega.(yj0-u) |   */

        fac  = iprod(pj, xj_u); /* fac = pj.(xj-u) */
        fac2 = fac * fac;

        /* Mass-weighting */
        wj = N_M * erg->m_loc[j];

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        svmul(-erg->rotg->k * wj * fac, pj, tmp_f);
        copy_rvec(tmp_f, erg->f_rot_loc[j]);
        sum += wj * fac2;

        /* If requested, also calculate the potential for a set of angles
         * near the current reference angle */
        if (bCalcPotFit)
        {
            for (int ifit = 0; ifit < erg->rotg->PotAngle_nstep; ifit++)
            {
                /* Index of this rotation group atom with respect to the whole rotation group */
                int jj = collectiveRotationGroupIndex[j];

                /* Rotate with the alternative angle. Like rotate_local_reference(),
                 * just for a single local atom */
                mvmul(erg->PotAngleFit->rotmat[ifit],
                      erg->referencePositions[jj],
                      fit_tmpvec); /* fit_tmpvec = Omega*(yj0-u) */

                /* Calculate Omega.(yj0-u) */
                cprod(erg->vec, fit_tmpvec, tmpvec); /* tmpvec = v x Omega.(yj0-u) */
                /*                                     *         v x Omega.(yj0-u)     */
                unitv(tmpvec, pj); /*  pj = ---------------------   */
                                   /*       | v x Omega.(yj0-u) |   */

                fac  = iprod(pj, xj_u); /* fac = pj.(xj-u) */
                fac2 = fac * fac;

                /* Add to the rotation potential for this angle: */
                erg->PotAngleFit->V[ifit] += 0.5 * erg->rotg->k * wj * fac2;
            }
        }

        if (bOutstepRot)
        {
            /* Add to the torque of this rotation group */
            erg->torque_v += torque(erg->vec, tmp_f, erg->x_loc_pbc[j], erg->xc_center);

            /* Calculate the angle between reference and actual rotation group atom. */
            angle(erg, xj_u, erg->xr_loc[j], &alpha, &weight); /* angle in rad, weighted */
            erg->angle_v += alpha * weight;
            erg->weight_v += weight;
        }

        PRINT_FORCE_J

    } /* end of loop over local rotation group atoms */
    erg->V = 0.5 * erg->rotg->k * sum;
}


/* Calculate the radial motion pivot-free potential and forces */
static void do_radial_motion_pf(gmx_enfrotgrp*                 erg,
                                gmx::ArrayRef<const gmx::RVec> coords, /* The positions */
                                const matrix box, /* The simulation box                         */
                                gmx_bool     bOutstepRot, /* Output to main rotation output file  */
                                gmx_bool     bOutstepSlab)    /* Output per-slab data */
{
    rvec     xj;      /* Current position */
    rvec     xj_xc;   /* xj  - xc  */
    rvec     yj0_yc0; /* yj0 - yc0 */
    rvec     tmp_f;   /* Force */
    real     alpha;   /* a single angle between an actual and a reference position */
    real     weight;  /* single weight for a single angle */
    rvec     tmpvec, tmpvec2;
    rvec     innersumvec; /* Precalculation of the inner sum */
    rvec     innersumveckM;
    real     fac, fac2, V = 0.0;
    rvec     qi, qj;
    gmx_bool bCalcPotFit;

    /* For mass weighting: */
    real mj, wi, wj; /* Mass-weighting of the positions */
    real N_M;        /* N/M */

    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (RotationGroupFitting::Pot == erg->rotg->eFittype);

    N_M = erg->rotg->nat * erg->invmass;

    /* Get the current center of the rotation group: */
    get_center(erg->xc, erg->mc, erg->rotg->nat, erg->xc_center);

    /* Precalculate Sum_i [ wi qi.(xi-xc) qi ] which is needed for every single j */
    clear_rvec(innersumvec);
    for (int i = 0; i < erg->rotg->nat; i++)
    {
        /* Mass-weighting */
        wi = N_M * erg->mc[i];

        /* Calculate qi. Note that xc_ref_center has already been subtracted from
         * x_ref_original in init_rot_group.*/
        mvmul(erg->rotmat, erg->referencePositions[i], tmpvec); /* tmpvec  = Omega.(yi0-yc0) */

        cprod(erg->vec, tmpvec, tmpvec2); /* tmpvec2 = v x Omega.(yi0-yc0) */

        /*                                             *         v x Omega.(yi0-yc0)     */
        unitv(tmpvec2, qi); /*  qi = -----------------------   */
                            /*       | v x Omega.(yi0-yc0) |   */

        rvec_sub(erg->xc[i], erg->xc_center, tmpvec); /* tmpvec = xi-xc */

        svmul(wi * iprod(qi, tmpvec), qi, tmpvec2);

        rvec_inc(innersumvec, tmpvec2);
    }
    svmul(erg->rotg->k * erg->invmass, innersumvec, innersumveckM);

    /* Each process calculates the forces on its local atoms */
    const auto& localRotationGroupIndex      = erg->atomSet->localIndex();
    const auto& collectiveRotationGroupIndex = erg->atomSet->collectiveIndex();
    for (gmx::Index j = 0; j < localRotationGroupIndex.ssize(); j++)
    {
        /* Local index of a rotation group atom  */
        int ii = localRotationGroupIndex[j];
        /* Position of this atom in the collective array */
        int iigrp = collectiveRotationGroupIndex[j];
        /* Mass-weighting */
        mj = erg->mc[iigrp]; /* need the unsorted mass here */
        wj = N_M * mj;

        /* Current position of this atom: x[ii][XX/YY/ZZ] */
        copy_rvec(coords[ii], xj);

        /* Shift this atom such that it is near its reference */
        shift_single_coord(box, xj, erg->xc_shifts[iigrp]);

        /* The (unrotated) reference position is yj0. yc0 has already
         * been subtracted in init_rot_group */
        copy_rvec(erg->referencePositions[iigrp], yj0_yc0); /* yj0_yc0 = yj0 - yc0      */

        /* Calculate Omega.(yj0-yc0) */
        mvmul(erg->rotmat, yj0_yc0, tmpvec2); /* tmpvec2 = Omega.(yj0 - yc0)  */

        cprod(erg->vec, tmpvec2, tmpvec); /* tmpvec = v x Omega.(yj0-yc0) */

        /*                     *         v x Omega.(yj0-yc0)     */
        unitv(tmpvec, qj); /*  qj = -----------------------   */
                           /*       | v x Omega.(yj0-yc0) |   */

        /* Calculate (xj-xc) */
        rvec_sub(xj, erg->xc_center, xj_xc); /* xj_xc = xj-xc */

        fac  = iprod(qj, xj_xc); /* fac = qj.(xj-xc) */
        fac2 = fac * fac;

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        svmul(-erg->rotg->k * wj * fac, qj, tmp_f); /* part 1 of force */
        svmul(mj, innersumveckM, tmpvec);           /* part 2 of force */
        rvec_inc(tmp_f, tmpvec);
        copy_rvec(tmp_f, erg->f_rot_loc[j]);
        V += wj * fac2;

        /* If requested, also calculate the potential for a set of angles
         * near the current reference angle */
        if (bCalcPotFit)
        {
            for (int ifit = 0; ifit < erg->rotg->PotAngle_nstep; ifit++)
            {
                /* Rotate with the alternative angle. Like rotate_local_reference(),
                 * just for a single local atom */
                mvmul(erg->PotAngleFit->rotmat[ifit], yj0_yc0, tmpvec2); /* tmpvec2 = Omega*(yj0-yc0) */

                /* Calculate Omega.(yj0-u) */
                cprod(erg->vec, tmpvec2, tmpvec); /* tmpvec = v x Omega.(yj0-yc0) */
                /*                                  *         v x Omega.(yj0-yc0)     */
                unitv(tmpvec, qj); /*  qj = -----------------------   */
                                   /*       | v x Omega.(yj0-yc0) |   */

                fac  = iprod(qj, xj_xc); /* fac = qj.(xj-xc) */
                fac2 = fac * fac;

                /* Add to the rotation potential for this angle: */
                erg->PotAngleFit->V[ifit] += 0.5 * erg->rotg->k * wj * fac2;
            }
        }

        if (bOutstepRot)
        {
            /* Add to the torque of this rotation group */
            erg->torque_v += torque(erg->vec, tmp_f, xj, erg->xc_center);

            /* Calculate the angle between reference and actual rotation group atom. */
            angle(erg, xj_xc, yj0_yc0, &alpha, &weight); /* angle in rad, weighted */
            erg->angle_v += alpha * weight;
            erg->weight_v += weight;
        }

        PRINT_FORCE_J

    } /* end of loop over local rotation group atoms */
    erg->V = 0.5 * erg->rotg->k * V;
}


/* Precalculate the inner sum for the radial motion 2 forces */
static void radial_motion2_precalc_inner_sum(const gmx_enfrotgrp* erg, rvec innersumvec)
{
    int  i;
    rvec xi_xc; /* xj - xc */
    rvec tmpvec, tmpvec2;
    real fac;
    rvec ri, si;
    real siri;
    rvec v_xi_xc; /* v x (xj - u) */
    real psii, psiistar;
    real wi;  /* Mass-weighting of the positions */
    real N_M; /* N/M */
    rvec sumvec;

    N_M = erg->rotg->nat * erg->invmass;

    /* Loop over the collective set of positions */
    clear_rvec(sumvec);
    for (i = 0; i < erg->rotg->nat; i++)
    {
        /* Mass-weighting */
        wi = N_M * erg->mc[i];

        rvec_sub(erg->xc[i], erg->xc_center, xi_xc); /* xi_xc = xi-xc         */

        /* Calculate ri. Note that xc_ref_center has already been subtracted from
         * x_ref_original in init_rot_group.*/
        mvmul(erg->rotmat, erg->referencePositions[i], ri); /* ri  = Omega.(yi0-yc0) */

        cprod(erg->vec, xi_xc, v_xi_xc); /* v_xi_xc = v x (xi-u)  */

        fac = norm2(v_xi_xc);
        /*                                 *                      1           */
        psiistar = 1.0 / (fac + erg->rotg->eps); /* psiistar = --------------------- */
        /*            |v x (xi-xc)|^2 + eps */

        psii = gmx::invsqrt(fac); /*                 1                */
                                  /*  psii    = -------------         */
                                  /*            |v x (xi-xc)|         */

        svmul(psii, v_xi_xc, si); /*  si = psii * (v x (xi-xc) )     */

        siri = iprod(si, ri); /* siri = si.ri           */

        svmul(psiistar / psii, ri, tmpvec);
        svmul(psiistar * psiistar / (psii * psii * psii) * siri, si, tmpvec2);
        rvec_dec(tmpvec, tmpvec2);
        cprod(tmpvec, erg->vec, tmpvec2);

        svmul(wi * siri, tmpvec2, tmpvec);

        rvec_inc(sumvec, tmpvec);
    }
    svmul(erg->rotg->k * erg->invmass, sumvec, innersumvec);
}


/* Calculate the radial motion 2 potential and forces */
static void do_radial_motion2(gmx_enfrotgrp*                 erg,
                              gmx::ArrayRef<const gmx::RVec> coords, /* The positions */
                              const matrix box,     /* The simulation box                         */
                              gmx_bool bOutstepRot, /* Output to main rotation output file        */
                              gmx_bool bOutstepSlab) /* Output per-slab data */
{
    rvec     xj;      /* Position */
    real     alpha;   /* a single angle between an actual and a reference position */
    real     weight;  /* single weight for a single angle */
    rvec     xj_u;    /* xj - u */
    rvec     yj0_yc0; /* yj0 -yc0 */
    rvec     tmpvec, tmpvec2;
    real     fac, fit_fac, fac2, Vpart = 0.0;
    rvec     rj, fit_rj, sj;
    real     sjrj;
    rvec     v_xj_u; /* v x (xj - u) */
    real     psij, psijstar;
    real     mj, wj; /* For mass-weighting of the positions */
    real     N_M;    /* N/M */
    gmx_bool bPF;
    rvec     innersumvec;
    gmx_bool bCalcPotFit;

    bPF         = erg->rotg->eType == EnforcedRotationGroupType::Rm2pf;
    bCalcPotFit = (bOutstepRot || bOutstepSlab) && (RotationGroupFitting::Pot == erg->rotg->eFittype);

    clear_rvec(yj0_yc0); /* Make the compiler happy */

    clear_rvec(innersumvec);
    if (bPF)
    {
        /* For the pivot-free variant we have to use the current center of
         * mass of the rotation group instead of the pivot u */
        get_center(erg->xc, erg->mc, erg->rotg->nat, erg->xc_center);

        /* Also, we precalculate the second term of the forces that is identical
         * (up to the weight factor mj) for all forces */
        radial_motion2_precalc_inner_sum(erg, innersumvec);
    }

    N_M = erg->rotg->nat * erg->invmass;

    /* Each process calculates the forces on its local atoms */
    const auto& localRotationGroupIndex      = erg->atomSet->localIndex();
    const auto& collectiveRotationGroupIndex = erg->atomSet->collectiveIndex();
    for (gmx::Index j = 0; j < localRotationGroupIndex.ssize(); j++)
    {
        if (bPF)
        {
            /* Local index of a rotation group atom  */
            int ii = localRotationGroupIndex[j];
            /* Position of this atom in the collective array */
            int iigrp = collectiveRotationGroupIndex[j];
            /* Mass-weighting */
            mj = erg->mc[iigrp];

            /* Current position of this atom: x[ii] */
            copy_rvec(coords[ii], xj);

            /* Shift this atom such that it is near its reference */
            shift_single_coord(box, xj, erg->xc_shifts[iigrp]);

            /* The (unrotated) reference position is yj0. yc0 has already
             * been subtracted in init_rot_group */
            copy_rvec(erg->referencePositions[iigrp], yj0_yc0); /* yj0_yc0 = yj0 - yc0  */

            /* Calculate Omega.(yj0-yc0) */
            mvmul(erg->rotmat, yj0_yc0, rj); /* rj = Omega.(yj0-yc0)  */
        }
        else
        {
            mj = erg->m_loc[j];
            copy_rvec(erg->x_loc_pbc[j], xj);
            copy_rvec(erg->xr_loc[j], rj); /* rj = Omega.(yj0-u)    */
        }
        /* Mass-weighting */
        wj = N_M * mj;

        /* Calculate (xj-u) resp. (xj-xc) */
        rvec_sub(xj, erg->xc_center, xj_u); /* xj_u = xj-u           */

        cprod(erg->vec, xj_u, v_xj_u); /* v_xj_u = v x (xj-u)   */

        fac = norm2(v_xj_u);
        /*                                      *                      1           */
        psijstar = 1.0 / (fac + erg->rotg->eps); /*  psistar = --------------------  */
        /*                                      *            |v x (xj-u)|^2 + eps  */

        psij = gmx::invsqrt(fac); /*                 1                */
                                  /*  psij    = ------------          */
                                  /*            |v x (xj-u)|          */

        svmul(psij, v_xj_u, sj); /*  sj = psij * (v x (xj-u) )       */

        fac  = iprod(v_xj_u, rj); /* fac = (v x (xj-u)).rj */
        fac2 = fac * fac;

        sjrj = iprod(sj, rj); /* sjrj = sj.rj          */

        svmul(psijstar / psij, rj, tmpvec);
        svmul(psijstar * psijstar / (psij * psij * psij) * sjrj, sj, tmpvec2);
        rvec_dec(tmpvec, tmpvec2);
        cprod(tmpvec, erg->vec, tmpvec2);

        /* Store the additional force so that it can be added to the force
         * array after the normal forces have been evaluated */
        svmul(-erg->rotg->k * wj * sjrj, tmpvec2, tmpvec);
        svmul(mj, innersumvec, tmpvec2); /* This is != 0 only for the pivot-free variant */

        rvec_add(tmpvec2, tmpvec, erg->f_rot_loc[j]);
        Vpart += wj * psijstar * fac2;

        /* If requested, also calculate the potential for a set of angles
         * near the current reference angle */
        if (bCalcPotFit)
        {
            for (int ifit = 0; ifit < erg->rotg->PotAngle_nstep; ifit++)
            {
                if (bPF)
                {
                    mvmul(erg->PotAngleFit->rotmat[ifit], yj0_yc0, fit_rj); /* fit_rj = Omega.(yj0-yc0) */
                }
                else
                {
                    /* Position of this atom in the collective array */
                    int iigrp = collectiveRotationGroupIndex[j];
                    /* Rotate with the alternative angle. Like rotate_local_reference(),
                     * just for a single local atom */
                    mvmul(erg->PotAngleFit->rotmat[ifit], erg->referencePositions[iigrp], fit_rj); /* fit_rj = Omega*(yj0-u) */
                }
                fit_fac = iprod(v_xj_u, fit_rj); /* fac = (v x (xj-u)).fit_rj */
                /* Add to the rotation potential for this angle: */
                erg->PotAngleFit->V[ifit] += 0.5 * erg->rotg->k * wj * psijstar * fit_fac * fit_fac;
            }
        }

        if (bOutstepRot)
        {
            /* Add to the torque of this rotation group */
            erg->torque_v += torque(erg->vec, erg->f_rot_loc[j], xj, erg->xc_center);

            /* Calculate the angle between reference and actual rotation group atom. */
            angle(erg, xj_u, rj, &alpha, &weight); /* angle in rad, weighted */
            erg->angle_v += alpha * weight;
            erg->weight_v += weight;
        }

        PRINT_FORCE_J

    } /* end of loop over local rotation group atoms */
    erg->V = 0.5 * erg->rotg->k * Vpart;
}


/* Determine the smallest and largest position vector (with respect to the
 * rotation vector) for the reference group */
static void get_firstlast_atom_ref(const gmx_enfrotgrp* erg, int* firstindex, int* lastindex)
{
    int  i;
    real xcproj;           /* The projection of a reference position on the
                              rotation vector */
    real minproj, maxproj; /* Smallest and largest projection on v */

    /* Start with some value */
    minproj = iprod(erg->referencePositions[0], erg->vec);
    maxproj = minproj;

    /* This is just to ensure that it still works if all the atoms of the
     * reference structure are situated in a plane perpendicular to the rotation
     * vector */
    *firstindex = 0;
    *lastindex  = erg->rotg->nat - 1;

    /* Loop over all atoms of the reference group,
     * project them on the rotation vector to find the extremes */
    for (i = 0; i < erg->rotg->nat; i++)
    {
        xcproj = iprod(erg->referencePositions[i], erg->vec);
        if (xcproj < minproj)
        {
            minproj     = xcproj;
            *firstindex = i;
        }
        if (xcproj > maxproj)
        {
            maxproj    = xcproj;
            *lastindex = i;
        }
    }
}


/* Allocate memory for the slabs */
static void allocate_slabs(gmx_enfrotgrp* erg, FILE* fplog, gmx_bool bVerbose)
{
    /* More slabs than are defined for the reference are never needed */
    int nslabs = erg->slab_last_ref - erg->slab_first_ref + 1;

    if ((nullptr != fplog) && bVerbose)
    {
        fprintf(fplog,
                "%s allocating memory to store data for %d slabs (rotation group %d).\n",
                RotStr.c_str(),
                nslabs,
                erg->groupIndex);
    }
    snew(erg->slab_center, nslabs);
    snew(erg->slab_center_ref, nslabs);
    snew(erg->slab_weights, nslabs);
    snew(erg->slab_torque_v, nslabs);
    snew(erg->slab_data, nslabs);
    snew(erg->gn_atom, nslabs);
    snew(erg->gn_slabind, nslabs);
    snew(erg->slab_innersumvec, nslabs);
    for (int i = 0; i < nslabs; i++)
    {
        snew(erg->slab_data[i].x, erg->rotg->nat);
        snew(erg->slab_data[i].ref, erg->rotg->nat);
        snew(erg->slab_data[i].weight, erg->rotg->nat);
    }
    snew(erg->xc_ref_sorted, erg->rotg->nat);
    snew(erg->xc_sortind, erg->rotg->nat);
    snew(erg->firstatom, nslabs);
    snew(erg->lastatom, nslabs);
}


/* From the extreme positions of the reference group, determine the first
 * and last slab of the reference. We can never have more slabs in the real
 * simulation than calculated here for the reference.
 */
static void get_firstlast_slab_ref(gmx_enfrotgrp* erg, real mc[], int ref_firstindex, int ref_lastindex)
{
    rvec dummy;

    int first = get_first_slab(erg, erg->referencePositions[ref_firstindex]);
    int last  = get_last_slab(erg, erg->referencePositions[ref_lastindex]);

    while (get_slab_weight(first, erg, erg->referencePositions, mc, &dummy) > WEIGHT_MIN)
    {
        first--;
    }
    erg->slab_first_ref = first + 1;
    while (get_slab_weight(last, erg, erg->referencePositions, mc, &dummy) > WEIGHT_MIN)
    {
        last++;
    }
    erg->slab_last_ref = last - 1;
}


/* Special version of copy_rvec:
 * During the copy procedure of xcurr to b, the correct PBC image is chosen
 * such that the copied vector ends up near its reference position xref */
static inline void copy_correct_pbc_image(const rvec xcurr, /* copy vector xcurr ... */
                                          rvec       b, /* ... to b ...                         */
                                          const rvec xref, /* choosing the PBC image such that b ends up near xref */
                                          const matrix box,
                                          int          npbcdim)
{
    rvec dx;
    int  d, m;
    ivec shift;


    /* Shortest PBC distance between the atom and its reference */
    rvec_sub(xcurr, xref, dx);

    /* Determine the shift for this atom */
    clear_ivec(shift);
    for (m = npbcdim - 1; m >= 0; m--)
    {
        while (dx[m] < -0.5 * box[m][m])
        {
            for (d = 0; d < DIM; d++)
            {
                dx[d] += box[m][d];
            }
            shift[m]++;
        }
        while (dx[m] >= 0.5 * box[m][m])
        {
            for (d = 0; d < DIM; d++)
            {
                dx[d] -= box[m][d];
            }
            shift[m]--;
        }
    }

    /* Apply the shift to the position */
    copy_rvec(xcurr, b);
    shift_single_coord(box, b, shift);
}


static void init_rot_group(FILE*             fplog,
                           const t_commrec*  cr,
                           gmx_enfrotgrp*    erg,
                           rvec*             x,
                           const gmx_mtop_t& mtop,
                           gmx_bool          bVerbose,
                           FILE*             out_slabs,
                           const matrix      box,
                           t_inputrec*       ir,
                           gmx_bool          bOutputCenters)
{
    rvec            coord, xref, *xdum;
    gmx_bool        bFlex, bColl;
    int             ref_firstindex, ref_lastindex;
    real            mass, totalmass;
    real            start = 0.0;
    double          t_start;
    const t_rotgrp* rotg = erg->rotg;


    /* Do we have a flexible axis? */
    bFlex = ISFLEX(rotg);
    /* Do we use a global set of coordinates? */
    bColl = ISCOLL(rotg);

    /* Allocate space for collective coordinates if needed */
    if (bColl)
    {
        snew(erg->xc, erg->rotg->nat);
        snew(erg->xc_shifts, erg->rotg->nat);
        snew(erg->xc_eshifts, erg->rotg->nat);
        snew(erg->xc_old, erg->rotg->nat);

        if (erg->rotg->eFittype == RotationGroupFitting::Norm)
        {
            snew(erg->xc_ref_length, erg->rotg->nat); /* in case fit type NORM is chosen */
            snew(erg->xc_norm, erg->rotg->nat);
        }
    }
    else
    {
        snew(erg->xr_loc, erg->rotg->nat);
        snew(erg->x_loc_pbc, erg->rotg->nat);
    }

    copy_rvec(erg->rotg->inputVec, erg->vec);
    snew(erg->f_rot_loc, erg->rotg->nat);

    /* Make space for the calculation of the potential at other angles (used
     * for fitting only) */
    if (RotationGroupFitting::Pot == erg->rotg->eFittype)
    {
        snew(erg->PotAngleFit, 1);
        snew(erg->PotAngleFit->degangle, erg->rotg->PotAngle_nstep);
        snew(erg->PotAngleFit->V, erg->rotg->PotAngle_nstep);
        snew(erg->PotAngleFit->rotmat, erg->rotg->PotAngle_nstep);

        /* Get the set of angles around the reference angle */
        start = -0.5 * (erg->rotg->PotAngle_nstep - 1) * erg->rotg->PotAngle_step;
        for (int i = 0; i < erg->rotg->PotAngle_nstep; i++)
        {
            erg->PotAngleFit->degangle[i] = start + i * erg->rotg->PotAngle_step;
        }
    }
    else
    {
        erg->PotAngleFit = nullptr;
    }

    /* Copy the masses so that the center can be determined. For all types of
     * enforced rotation, we store the masses in the erg->mc array. */
    snew(erg->mc, erg->rotg->nat);
    if (bFlex)
    {
        snew(erg->mc_sorted, erg->rotg->nat);
    }
    if (!bColl)
    {
        snew(erg->m_loc, erg->rotg->nat);
    }
    totalmass = 0.0;
    int molb  = 0;
    for (int i = 0; i < erg->rotg->nat; i++)
    {
        if (erg->rotg->bMassW)
        {
            mass = mtopGetAtomMass(mtop, erg->rotg->ind[i], &molb);
        }
        else
        {
            mass = 1.0;
        }
        erg->mc[i] = mass;
        totalmass += mass;
    }
    erg->invmass = 1.0 / totalmass;

    /* Set xc_ref_center for any rotation potential */
    if ((erg->rotg->eType == EnforcedRotationGroupType::Iso)
        || (erg->rotg->eType == EnforcedRotationGroupType::Pm)
        || (erg->rotg->eType == EnforcedRotationGroupType::Rm)
        || (erg->rotg->eType == EnforcedRotationGroupType::Rm2))
    {
        /* Set the pivot point for the fixed, stationary-axis potentials. This
         * won't change during the simulation */
        copy_rvec(erg->rotg->pivot, erg->xc_ref_center);
        copy_rvec(erg->rotg->pivot, erg->xc_center);
    }
    else
    {
        /* Center of the reference positions */
        get_center(as_rvec_array(erg->rotg->x_ref_original.data()), erg->mc, erg->rotg->nat, erg->xc_ref_center);

        /* Center of the actual positions */
        if (MAIN(cr))
        {
            snew(xdum, erg->rotg->nat);
            for (int i = 0; i < erg->rotg->nat; i++)
            {
                int ii = erg->rotg->ind[i];
                copy_rvec(x[ii], xdum[i]);
            }
            get_center(xdum, erg->mc, erg->rotg->nat, erg->xc_center);
            sfree(xdum);
        }
#if GMX_MPI
        if (PAR(cr))
        {
            gmx_bcast(sizeof(erg->xc_center), erg->xc_center, cr->mpi_comm_mygroup);
        }
#endif
    }

    erg->referencePositions = erg->rotg->x_ref_original;
    if (bColl)
    {
        /* Save the original (whole) set of positions in xc_old such that at later
         * steps the rotation group can always be made whole again. If the simulation is
         * restarted, we compute the starting reference positions (given the time)
         * and assume that the correct PBC image of each position is the one nearest
         * to the current reference */
        if (MAIN(cr))
        {
            /* Calculate the rotation matrix for this angle: */
            t_start       = ir->init_t + ir->init_step * ir->delta_t;
            erg->degangle = erg->rotg->rate * t_start;
            calc_rotmat(erg->vec, erg->degangle, erg->rotmat);

            for (int i = 0; i < erg->rotg->nat; i++)
            {
                int ii = erg->rotg->ind[i];

                /* Subtract pivot, rotate, and add pivot again. This will yield the
                 * reference position for time t */
                rvec_sub(erg->referencePositions[i], erg->xc_ref_center, coord);
                mvmul(erg->rotmat, coord, xref);
                rvec_inc(xref, erg->xc_ref_center);

                copy_correct_pbc_image(x[ii], erg->xc_old[i], xref, box, 3);
            }
        }
#if GMX_MPI
        if (PAR(cr))
        {
            gmx_bcast(erg->rotg->nat * sizeof(erg->xc_old[0]), erg->xc_old, cr->mpi_comm_mygroup);
        }
#endif
    }

    if ((erg->rotg->eType != EnforcedRotationGroupType::Flex)
        && (erg->rotg->eType != EnforcedRotationGroupType::Flex2))
    {
        /* Put the reference positions into origin: */
        for (int i = 0; i < erg->rotg->nat; i++)
        {
            erg->referencePositions[i] -= erg->xc_ref_center;
        }
    }

    /* Enforced rotation with flexible axis */
    if (bFlex)
    {
        /* Calculate maximum beta value from minimum gaussian (performance opt.) */
        erg->max_beta = calc_beta_max(erg->rotg->min_gaussian, erg->rotg->slab_dist);

        /* Determine the smallest and largest coordinate with respect to the rotation vector */
        get_firstlast_atom_ref(erg, &ref_firstindex, &ref_lastindex);

        /* From the extreme positions of the reference group, determine the first
         * and last slab of the reference. */
        get_firstlast_slab_ref(erg, erg->mc, ref_firstindex, ref_lastindex);

        /* Allocate memory for the slabs */
        allocate_slabs(erg, fplog, bVerbose);

        /* Flexible rotation: determine the reference centers for the rest of the simulation */
        erg->slab_first = erg->slab_first_ref;
        erg->slab_last  = erg->slab_last_ref;
        get_slab_centers(erg, erg->referencePositions, erg->mc, -1, out_slabs, bOutputCenters, TRUE);

        /* Length of each x_rotref vector from center (needed if fit routine NORM is chosen): */
        if (erg->rotg->eFittype == RotationGroupFitting::Norm)
        {
            for (int i = 0; i < erg->rotg->nat; i++)
            {
                rvec_sub(erg->referencePositions[i], erg->xc_ref_center, coord);
                erg->xc_ref_length[i] = norm(coord);
            }
        }
    }
}

/* Calculate the size of the MPI buffer needed in reduce_output() */
static int calc_mpi_bufsize(const gmx_enfrot* er)

{
    int count_total = 0;
    for (int g = 0; g < gmx::ssize(er->rot->grp); g++)
    {
        const t_rotgrp*      rotg = &er->rot->grp[g];
        const gmx_enfrotgrp* erg  = &er->enfrotgrp[g];

        /* Count the items that are transferred for this group: */
        int count_group = 4; /* V, torque, angle, weight */

        /* Add the maximum number of slabs for flexible groups */
        if (ISFLEX(rotg))
        {
            count_group += erg->slab_last_ref - erg->slab_first_ref + 1;
        }

        /* Add space for the potentials at different angles: */
        if (RotationGroupFitting::Pot == erg->rotg->eFittype)
        {
            count_group += erg->rotg->PotAngle_nstep;
        }

        /* Add to the total number: */
        count_total += count_group;
    }

    return count_total;
}


std::unique_ptr<gmx::EnforcedRotation> init_rot(FILE*                       fplog,
                                                t_inputrec*                 ir,
                                                int                         nfile,
                                                const t_filenm              fnm[],
                                                const t_commrec*            cr,
                                                gmx::LocalAtomSetManager*   atomSets,
                                                const t_state*              globalState,
                                                const gmx_mtop_t&           mtop,
                                                const gmx_output_env_t*     oenv,
                                                const gmx::MdrunOptions&    mdrunOptions,
                                                const gmx::StartingBehavior startingBehavior)
{
    int nat_max = 0; /* Size of biggest rotation group */

    if (MAIN(cr) && mdrunOptions.verbose)
    {
        fprintf(stdout, "%s Initializing ...\n", RotStr.c_str());
    }

    auto        enforcedRotation = std::make_unique<gmx::EnforcedRotation>();
    gmx_enfrot* er               = enforcedRotation->getLegacyEnfrot();
    // TODO When this module implements IMdpOptions, the ownership will become more clear.
    er->rot                  = ir->rot.get();
    er->restartWithAppending = (startingBehavior == gmx::StartingBehavior::RestartWithAppending);

    /* When appending, skip first output to avoid duplicate entries in the data files */
    er->bOut = !er->restartWithAppending;

    if (MAIN(cr) && er->bOut)
    {
        please_cite(fplog, "Kutzner2011");
    }

    /* Output every step for reruns */
    if (mdrunOptions.rerun)
    {
        if (nullptr != fplog)
        {
            fprintf(fplog, "%s rerun - will write rotation output every available step.\n", RotStr.c_str());
        }
        er->nstrout = 1;
        er->nstsout = 1;
    }
    else
    {
        er->nstrout = er->rot->nstrout;
        er->nstsout = er->rot->nstsout;
    }

    er->out_slabs = nullptr;
    if (MAIN(cr) && HaveFlexibleGroups(er->rot))
    {
        er->out_slabs = open_slab_out(opt2fn("-rs", nfile, fnm), er);
    }

    /* Space for the pbc-correct atom positions */
    std::vector<gmx::RVec> x_pbc;

    if (MAIN(cr))
    {
        /* Remove pbc, make molecule whole.
         * When ir->bContinuation=TRUE this has already been done, but ok. */
        x_pbc.resize(mtop.natoms);
        std::copy(globalState->x.begin(), globalState->x.end(), x_pbc.begin());
        do_pbc_first_mtop(nullptr, ir->pbcType, false, nullptr, globalState->box, &mtop, x_pbc, {});
        /* All molecules will be whole now, but not necessarily in the home box.
         * Additionally, if a rotation group consists of more than one molecule
         * (e.g. two strands of DNA), each one of them can end up in a different
         * periodic box. This is taken care of in init_rot_group.  */
    }

    /* Allocate space for the per-rotation-group data: */
    er->enfrotgrp.resize(er->rot->grp.size());
    int groupIndex = 0;
    for (auto& ergRef : er->enfrotgrp)
    {
        gmx_enfrotgrp* erg = &ergRef;
        erg->rotg          = &er->rot->grp[groupIndex];
        erg->atomSet       = std::make_unique<gmx::LocalAtomSet>(
                atomSets->add({ erg->rotg->ind, erg->rotg->ind + erg->rotg->nat }));
        erg->groupIndex = groupIndex;

        if (nullptr != fplog)
        {
            fprintf(fplog,
                    "%s group %d type '%s'\n",
                    RotStr.c_str(),
                    groupIndex,
                    enumValueToString(erg->rotg->eType));
        }

        if (erg->rotg->nat > 0)
        {
            nat_max = std::max(nat_max, erg->rotg->nat);

            init_rot_group(fplog,
                           cr,
                           erg,
                           as_rvec_array(x_pbc.data()),
                           mtop,
                           mdrunOptions.verbose,
                           er->out_slabs,
                           MAIN(cr) ? globalState->box : nullptr,
                           ir,
                           !er->restartWithAppending); /* Do not output the reference centers
                                                        * again if we are appending */
        }
        ++groupIndex;
    }

    /* Allocate space for enforced rotation buffer variables */
    er->bufsize = nat_max;
    snew(er->data, nat_max);
    snew(er->xbuf, nat_max);
    snew(er->mbuf, nat_max);

    /* Buffers for MPI reducing torques, angles, weights (for each group), and V */
    if (PAR(cr))
    {
        er->mpi_bufsize = calc_mpi_bufsize(er) + 100; /* larger to catch errors */
        snew(er->mpi_inbuf, er->mpi_bufsize);
        snew(er->mpi_outbuf, er->mpi_bufsize);
    }
    else
    {
        er->mpi_bufsize = 0;
        er->mpi_inbuf   = nullptr;
        er->mpi_outbuf  = nullptr;
    }

    /* Only do I/O on the MAIN */
    er->out_angles = nullptr;
    er->out_rot    = nullptr;
    er->out_torque = nullptr;
    if (MAIN(cr))
    {
        er->out_rot = open_rot_out(opt2fn("-ro", nfile, fnm), oenv, er);

        if (er->nstsout > 0)
        {
            if (HaveFlexibleGroups(er->rot) || HavePotFitGroups(er->rot))
            {
                er->out_angles = open_angles_out(opt2fn("-ra", nfile, fnm), er);
            }
            if (HaveFlexibleGroups(er->rot))
            {
                er->out_torque = open_torque_out(opt2fn("-rt", nfile, fnm), er);
            }
        }
    }
    return enforcedRotation;
}

/* Rotate the local reference positions and store them in
 * erg->xr_loc[0...(nat_loc-1)]
 *
 * Note that we already subtracted u or y_c from the reference positions
 * in init_rot_group().
 */
static void rotate_local_reference(gmx_enfrotgrp* erg)
{
    const auto& collectiveRotationGroupIndex = erg->atomSet->collectiveIndex();
    for (size_t i = 0; i < erg->atomSet->numAtomsLocal(); i++)
    {
        /* Index of this rotation group atom with respect to the whole rotation group */
        int ii = collectiveRotationGroupIndex[i];
        /* Rotate */
        mvmul(erg->rotmat, erg->referencePositions[ii], erg->xr_loc[i]);
    }
}


/* Select the PBC representation for each local x position and store that
 * for later usage. We assume the right PBC image of an x is the one nearest to
 * its rotated reference */
static void choose_pbc_image(gmx::ArrayRef<const gmx::RVec> coords, gmx_enfrotgrp* erg, const matrix box, int npbcdim)
{
    const auto& localRotationGroupIndex = erg->atomSet->localIndex();
    for (gmx::Index i = 0; i < localRotationGroupIndex.ssize(); i++)
    {
        /* Index of a rotation group atom  */
        int ii = localRotationGroupIndex[i];

        /* Get the correctly rotated reference position. The pivot was already
         * subtracted in init_rot_group() from the reference positions. Also,
         * the reference positions have already been rotated in
         * rotate_local_reference(). For the current reference position we thus
         * only need to add the pivot again. */
        rvec xref;
        copy_rvec(erg->xr_loc[i], xref);
        rvec_inc(xref, erg->xc_ref_center);

        copy_correct_pbc_image(coords[ii], erg->x_loc_pbc[i], xref, box, npbcdim);
    }
}


void do_rotation(const t_commrec*               cr,
                 gmx_enfrot*                    er,
                 const matrix                   box,
                 gmx::ArrayRef<const gmx::RVec> coords,
                 real                           t,
                 int64_t                        step,
                 gmx_bool                       bNS)
{
    gmx_bool    outstep_slab, outstep_rot;
    gmx_bool    bColl;
    rvec        transvec;
    gmx_potfit* fit = nullptr; /* For fit type 'potential' determine the fit
                                    angle via the potential minimum            */

    GMX_ASSERT(er, "Enforced rotation needs a valid pointer to its data object");

#ifdef TAKETIME
    double t0;
#endif

    /* When to output in main rotation output file */
    outstep_rot = do_per_step(step, er->nstrout) && er->bOut;
    /* When to output per-slab data */
    outstep_slab = do_per_step(step, er->nstsout) && er->bOut;

    /* Output time into rotation output file */
    if (outstep_rot && MAIN(cr))
    {
        fprintf(er->out_rot, "%12.3e", t);
    }

    /**************************************************************************/
    /* First do ALL the communication! */
    for (auto& ergRef : er->enfrotgrp)
    {
        gmx_enfrotgrp*  erg  = &ergRef;
        const t_rotgrp* rotg = erg->rotg;

        /* Do we use a collective (global) set of coordinates? */
        bColl = ISCOLL(rotg);

        /* Calculate the rotation matrix for this angle: */
        erg->degangle = rotg->rate * t;
        calc_rotmat(erg->vec, erg->degangle, erg->rotmat);

        if (bColl)
        {
            /* Transfer the rotation group's positions such that every node has
             * all of them. Every node contributes its local positions x and stores
             * it in the collective erg->xc array. */
            communicate_group_positions(cr,
                                        erg->xc,
                                        erg->xc_shifts,
                                        erg->xc_eshifts,
                                        bNS,
                                        as_rvec_array(coords.data()),
                                        rotg->nat,
                                        erg->atomSet->numAtomsLocal(),
                                        erg->atomSet->localIndex().data(),
                                        erg->atomSet->collectiveIndex().data(),
                                        erg->xc_old,
                                        box);
        }
        else
        {
            /* Fill the local masses array;
             * this array changes in DD/neighborsearching steps */
            if (bNS)
            {
                const auto& collectiveRotationGroupIndex = erg->atomSet->collectiveIndex();
                for (gmx::Index i = 0; i < collectiveRotationGroupIndex.ssize(); i++)
                {
                    /* Index of local atom w.r.t. the collective rotation group */
                    int ii        = collectiveRotationGroupIndex[i];
                    erg->m_loc[i] = erg->mc[ii];
                }
            }

            /* Calculate Omega*(y_i-y_c) for the local positions */
            rotate_local_reference(erg);

            /* Choose the nearest PBC images of the group atoms with respect
             * to the rotated reference positions */
            choose_pbc_image(coords, erg, box, 3);

            /* Get the center of the rotation group */
            if ((rotg->eType == EnforcedRotationGroupType::Isopf)
                || (rotg->eType == EnforcedRotationGroupType::Pmpf))
            {
                get_center_comm(
                        cr, erg->x_loc_pbc, erg->m_loc, erg->atomSet->numAtomsLocal(), rotg->nat, erg->xc_center);
            }
        }

    } /* End of loop over rotation groups */

    /**************************************************************************/
    /* Done communicating, we can start to count cycles for the load balancing now ... */
    if (haveDDAtomOrdering(*cr))
    {
        ddReopenBalanceRegionCpu(cr->dd);
    }

#ifdef TAKETIME
    t0 = MPI_Wtime();
#endif

    for (auto& ergRef : er->enfrotgrp)
    {
        gmx_enfrotgrp*  erg  = &ergRef;
        const t_rotgrp* rotg = erg->rotg;

        if (outstep_rot && MAIN(cr))
        {
            fprintf(er->out_rot, "%12.4f", erg->degangle);
        }

        /* Calculate angles and rotation matrices for potential fitting: */
        if ((outstep_rot || outstep_slab) && (RotationGroupFitting::Pot == rotg->eFittype))
        {
            fit = erg->PotAngleFit;
            for (int i = 0; i < rotg->PotAngle_nstep; i++)
            {
                calc_rotmat(erg->vec, erg->degangle + fit->degangle[i], fit->rotmat[i]);

                /* Clear value from last step */
                erg->PotAngleFit->V[i] = 0.0;
            }
        }

        /* Clear values from last time step */
        erg->V        = 0.0;
        erg->torque_v = 0.0;
        erg->angle_v  = 0.0;
        erg->weight_v = 0.0;

        switch (rotg->eType)
        {
            case EnforcedRotationGroupType::Iso:
            case EnforcedRotationGroupType::Isopf:
            case EnforcedRotationGroupType::Pm:
            case EnforcedRotationGroupType::Pmpf: do_fixed(erg, outstep_rot, outstep_slab); break;
            case EnforcedRotationGroupType::Rm:
                do_radial_motion(erg, outstep_rot, outstep_slab);
                break;
            case EnforcedRotationGroupType::Rmpf:
                do_radial_motion_pf(erg, coords, box, outstep_rot, outstep_slab);
                break;
            case EnforcedRotationGroupType::Rm2:
            case EnforcedRotationGroupType::Rm2pf:
                do_radial_motion2(erg, coords, box, outstep_rot, outstep_slab);
                break;
            case EnforcedRotationGroupType::Flext:
            case EnforcedRotationGroupType::Flex2t:
                /* Subtract the center of the rotation group from the collective positions array
                 * Also store the center in erg->xc_center since it needs to be subtracted
                 * in the low level routines from the local coordinates as well */
                get_center(erg->xc, erg->mc, rotg->nat, erg->xc_center);
                svmul(-1.0, erg->xc_center, transvec);
                translate_x(erg->xc, rotg->nat, transvec);
                do_flexible(MAIN(cr), er, erg, coords, box, t, outstep_rot, outstep_slab);
                break;
            case EnforcedRotationGroupType::Flex:
            case EnforcedRotationGroupType::Flex2:
                /* Do NOT subtract the center of mass in the low level routines! */
                clear_rvec(erg->xc_center);
                do_flexible(MAIN(cr), er, erg, coords, box, t, outstep_rot, outstep_slab);
                break;
            default: gmx_fatal(FARGS, "No such rotation potential.");
        }
    }

#ifdef TAKETIME
    if (MAIN(cr))
    {
        fprintf(stderr,
                "%s calculation (step %" PRId64 ") took %g seconds.\n",
                RotStr.c_str(),
                step,
                MPI_Wtime() - t0);
    }
#endif
}
