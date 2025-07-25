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

#include "edsam.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/nrjac.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/edsamhistory.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

namespace
{
/*! \brief Identifies the type of ED: none, normal ED, flooding. */
enum class EssentialDynamicsType
{
    //! No essential dynamics
    None,
    //! Essential dynamics sampling
    EDSampling,
    //! Essential dynamics flooding
    Flooding
};

/*! \brief Identify on which structure to operate. */
enum class EssentialDynamicsStructure
{
    //! Reference structure
    Reference,
    //! Average Structure
    Average,
    //! Origin Structure
    Origin,
    //! Target Structure
    Target
};

/*! \brief Essential dynamics vector.
 * TODO: split into working data and input data
 * NOTE: called eigvec, because traditionally eigenvectors from PCA analysis
 * were used as essential dynamics vectors, however, vectors used for ED need
 * not have any special properties
 */
struct t_eigvec
{
    //! nr of eigenvectors
    int neig = 0;
    //! index nrs of eigenvectors
    int* ieig = nullptr;
    //! stepsizes (per eigenvector)
    real* stpsz = nullptr;
    //! eigenvector components
    rvec** vec = nullptr;
    //! instantaneous x projections
    real* xproj = nullptr;
    //! instantaneous f projections
    real* fproj = nullptr;
    //! instantaneous radius
    real radius = 0.;
    //! starting or target projections
    real* refproj = nullptr;
};

/*! \brief Essential dynamics vectors per method implementation.
 */
struct t_edvecs
{
    //! only monitored, no constraints
    t_eigvec mon = {};
    //! fixed linear constraints
    t_eigvec linfix = {};
    //! acceptance linear constraints
    t_eigvec linacc = {};
    //! fixed radial constraints (exp)
    t_eigvec radfix = {};
    //! acceptance radial constraints (exp)
    t_eigvec radacc = {};
    //! acceptance rad. contraction constr.
    t_eigvec radcon = {};
};

/*! \brief Essential dynamics flooding parameters and work data.
 * TODO: split into working data and input parameters
 * NOTE: The implementation here follows:
 * O.E. Lange, L.V. Schafer, and H. Grubmuller,
 * “Flooding in GROMACS: Accelerated barrier crossings in molecular dynamics,”
 *  J. Comp. Chem., 27 1693–1702 (2006)
 */
struct t_edflood
{
    /*! \brief Target destabilisation free energy.
     * Controls the flooding potential strength.
     * Linked to the expected speed-up of mean first passage time out of flooded minimum */
    real deltaF0 = 0;
    //! Do not calculate a flooding potential, instead flood with a constant force
    bool bConstForce = false;
    //! Relaxation time scale for slowest flooded degree of freedom
    real tau = 0;
    //! Current estimated destabilisation free energy
    real deltaF = 0;
    //! Flooding energy, acting as a proportionality factor for the flooding potential
    real Efl = 0;
    //! Boltzmann constant times temperature, provided by user
    real kT = 0;
    //! The flooding potential
    real Vfl = 0;
    //! Integration time step
    real dt = 0;
    //! Inital flooding strenth
    real constEfl = 0;
    //! Empirical global scaling parameter of the essential dynamics vectors.
    real alpha2 = 0;
    //! The forces from flooding in atom coordinate space (in contrast to projections onto essential dynamics vectors)
    rvec* forces_cartesian = nullptr;
    //! The vectors along which to flood
    t_eigvec vecs = {};
    //! Use flooding for harmonic restraint on the eigenvector
    bool bHarmonic = false;
    //! The initial reference projection of the flooding vectors. Only with harmonic restraint.
    real* initialReferenceProjection = nullptr;
    //! The current reference projection is the initialReferenceProjection + step * slope. Only with harmonic restraint.
    real* referenceProjectionSlope = nullptr;
};
} // namespace


/* This type is for the average, reference, target, and origin structure    */
struct gmx_edx
{
    int  nr         = 0;       /* number of atoms this structure contains  */
    int  nr_loc     = 0;       /* number of atoms on local node            */
    int* anrs       = nullptr; /* atom index numbers                       */
    int* anrs_loc   = nullptr; /* local atom index numbers                 */
    int  nalloc_loc = 0;       /* allocation size of anrs_loc              */
    int* c_ind      = nullptr; /* at which position of the whole anrs
                                * array is a local atom?, i.e.
                                * c_ind[0...nr_loc-1] gives the atom index
                                * with respect to the collective
                                * anrs[0...nr-1] array                     */
    rvec* x     = nullptr;     /* positions for this structure             */
    rvec* x_old = nullptr;     /* Last positions which have the correct PBC
                                  representation of the ED group. In
                                  combination with keeping track of the
                                  shift vectors, the ED group can always
                                  be made whole                            */
    real* m     = nullptr;     /* masses                                   */
    real  mtot  = 0.;          /* total mass (only used in sref)           */
    real* sqrtm = nullptr;     /* sqrt of the masses used for mass-
                                * weighting of analysis (only used in sav) */
};


typedef struct edpar
{
    int      nini     = 0;     /* total Nr of atoms                    */
    gmx_bool fitmas   = false; /* true if trans fit with cm            */
    gmx_bool pcamas   = false; /* true if mass-weighted PCA            */
    int      presteps = 0;     /* number of steps to run without any
                                *    perturbations ... just monitoring */
    int outfrq     = 0;        /* freq (in steps) of writing to edo    */
    int maxedsteps = 0;        /* max nr of steps per cycle            */

    /* all gmx_edx datasets are copied to all nodes in the parallel case   */
    gmx_edx sref = {};         /* reference positions, to these fitting
                                * will be done                         */
    gmx_bool bRefEqAv = false; /* If true, reference & average indices
                                * are the same. Used for optimization  */
    gmx_edx sav  = {};         /* average positions                    */
    gmx_edx star = {};         /* target positions                     */
    gmx_edx sori = {};         /* origin positions                     */

    t_edvecs vecs  = {}; /* eigenvectors                         */
    real     slope = 0;  /* minimal slope in acceptance radexp   */

    t_edflood           flood = {};      /* parameters especially for flooding   */
    struct t_ed_buffer* buf   = nullptr; /* handle to local buffers              */
} t_edpar;


struct gmx_edsam
{
    ~gmx_edsam();
    //! The type of ED
    EssentialDynamicsType eEDtype = EssentialDynamicsType::None;
    //! output file pointer
    FILE*                edo = nullptr;
    std::vector<t_edpar> edpar;
    gmx_bool             bFirst = false;
};
gmx_edsam::~gmx_edsam()
{
    if (edo)
    {
        /* edo is opened sometimes with xvgropen, sometimes with
         * gmx_fio_fopen, so we use the least common denominator for
         * closing. */
        gmx_fio_fclose(edo);
    }
}

struct t_do_edsam
{
    real  oldrad;
    rvec* xcoll;                  /* Positions from all nodes, this is the
                                     collective set we work on.
                                     These are the positions of atoms with
                                     average structure indices */
    rvec*    xc_ref;              /* same but with reference structure indices */
    ivec*    shifts_xcoll;        /* Shifts for xcoll  */
    ivec*    extra_shifts_xcoll;  /* xcoll shift changes since last NS step */
    ivec*    shifts_xc_ref;       /* Shifts for xc_ref */
    ivec*    extra_shifts_xc_ref; /* xc_ref shift changes since last NS step */
    gmx_bool bUpdateShifts;       /* TRUE in NS steps to indicate that the
                                     ED shifts for this ED group need to
                                     be updated */
};


/* definition of ED buffer structure */
struct t_ed_buffer
{
    struct t_fit_to_ref* fit_to_ref;
    struct t_do_edfit*   do_edfit;
    struct t_do_edsam*   do_edsam;
    struct t_do_radcon*  do_radcon;
};

namespace gmx
{
class EssentialDynamics::Impl
{
public:
    // TODO: move all fields from gmx_edsam here and remove gmx_edsam
    gmx_edsam essentialDynamics_;
};
EssentialDynamics::EssentialDynamics() : impl_(new Impl) {}
EssentialDynamics::~EssentialDynamics() = default;

gmx_edsam* EssentialDynamics::getLegacyED()
{
    return &impl_->essentialDynamics_;
}
} // namespace gmx

/* Function declarations */
static void fit_to_reference(rvec* xcoll, rvec transvec, matrix rotmat, t_edpar* edi);
static void translate_and_rotate(rvec* x, int nat, rvec transvec, matrix rotmat);
static real rmsd_from_structure(rvec* x, struct gmx_edx* s);
namespace
{
/*! \brief Read in the essential dynamics input file and return its contents.
 * \param[in] fn the file name of the edi file to be read
 * \param[in] nr_mdatoms the number of atoms in the simulation
 * \returns A vector of containing the essentiail dyanmics parameters.
 * NOTE: edi files may that it may contain several ED data sets from concatenated edi files.
 * The standard case would be a single ED data set, though. */
std::vector<t_edpar> read_edi_file(const char* fn, int nr_mdatoms);
} // namespace
static void crosscheck_edi_file_vs_checkpoint(const gmx_edsam& ed, edsamhistory_t* EDstate);
static void init_edsamstate(const gmx_edsam& ed, edsamhistory_t* EDstate);
static void write_edo_legend(gmx_edsam* ed, int nED, const gmx_output_env_t* oenv);
/* End function declarations */

/* Define formats for the column width in the output file */
const char EDcol_sfmt[] = "%17s";
const char EDcol_efmt[] = "%17.5e";
const char EDcol_ffmt[] = "%17f";

/* Define a formats for reading, otherwise cppcheck complains for scanf without width limits */
const char max_ev_fmt_d[]       = "%7d"; /* Number of eigenvectors. 9,999,999 should be enough */
const char max_ev_fmt_lf[]      = "%12lf";
const char max_ev_fmt_dlf[]     = "%7d%12lf";
const char max_ev_fmt_dlflflf[] = "%7d%12lf%12lf%12lf";
const char max_ev_fmt_lelele[]  = "%12le%12le%12le";

namespace
{
/*!\brief Determines whether to perform essential dynamics constraints other than flooding.
 * \param[in] edi the essential dynamics parameters
 * \returns true if essential dyanmics constraints need to be performed
 */
bool bNeedDoEdsam(const t_edpar& edi)
{
    return (edi.vecs.mon.neig != 0) || (edi.vecs.linfix.neig != 0) || (edi.vecs.linacc.neig != 0)
           || (edi.vecs.radfix.neig != 0) || (edi.vecs.radacc.neig != 0) || (edi.vecs.radcon.neig != 0);
}
} // namespace


/* Multiple ED groups will be labeled with letters instead of numbers
 * to avoid confusion with eigenvector indices */
static char get_EDgroupChar(int nr_edi, int nED)
{
    if (nED == 1)
    {
        return ' ';
    }

    /* nr_edi = 1 -> A
     * nr_edi = 2 -> B ...
     */
    return 'A' + nr_edi - 1;
}

namespace
{
/*! \brief The mass-weighted inner product of two coordinate vectors.
 * Does not subtract average positions, projection on single eigenvector is returned
 * used by: do_linfix, do_linacc, do_radfix, do_radacc, do_radcon
 * Average position is subtracted in ed_apply_constraints prior to calling projectx
 * \param[in] edi Essential dynamics parameters
 * \param[in] xcoll vector of atom coordinates
 * \param[in] vec vector of coordinates to project onto
 * \return mass-weighted projection.
 */
real projectx(const t_edpar& edi, rvec* xcoll, rvec* vec)
{
    int  i;
    real proj = 0.0;


    for (i = 0; i < edi.sav.nr; i++)
    {
        proj += edi.sav.sqrtm[i] * iprod(vec[i], xcoll[i]);
    }

    return proj;
}
/*!\brief Project coordinates onto vector after substracting average position.
 * projection is stored in vec->refproj which is used for radacc, radfix,
 * radcon and center of flooding potential.
 * Average positions are first substracted from x, then added back again.
 * \param[in] edi essential dynamics parameters with average position
 * \param[in] x Coordinates to be projected
 * \param[out] vec eigenvector, radius and refproj are overwritten here
 */
void rad_project(const t_edpar& edi, rvec* x, t_eigvec* vec)
{
    int  i;
    real rad = 0.0;

    /* Subtract average positions */
    for (i = 0; i < edi.sav.nr; i++)
    {
        rvec_dec(x[i], edi.sav.x[i]);
    }

    for (i = 0; i < vec->neig; i++)
    {
        vec->refproj[i] = projectx(edi, x, vec->vec[i]);
        rad += gmx::square((vec->refproj[i] - vec->xproj[i]));
    }
    vec->radius = std::sqrt(rad);

    /* Add average positions */
    for (i = 0; i < edi.sav.nr; i++)
    {
        rvec_inc(x[i], edi.sav.x[i]);
    }
}

/*!\brief Projects coordinates onto eigenvectors and stores result in vec->xproj.
 * Mass-weighting is applied. Subtracts average positions prior to projection and add
 * them afterwards to retain the unchanged vector.
 * \param[in] x The coordinates to project to an eigenvector
 * \param[in,out] vec The eigenvectors
 * \param[in] edi essential dynamics parameters holding average structure and masses
 */
void project_to_eigvectors(rvec* x, t_eigvec* vec, const t_edpar& edi)
{
    if (!vec->neig)
    {
        return;
    }

    /* Subtract average positions */
    for (int i = 0; i < edi.sav.nr; i++)
    {
        rvec_dec(x[i], edi.sav.x[i]);
    }

    for (int i = 0; i < vec->neig; i++)
    {
        vec->xproj[i] = projectx(edi, x, vec->vec[i]);
    }

    /* Add average positions */
    for (int i = 0; i < edi.sav.nr; i++)
    {
        rvec_inc(x[i], edi.sav.x[i]);
    }
}
} // namespace

/* Project vector x onto all edi->vecs (mon, linfix,...) */
static void project(rvec*    x,   /* positions to project */
                    t_edpar* edi) /* edi data set */
{
    /* It is not more work to subtract the average position in every
     * subroutine again, because these routines are rarely used simultaneously */
    project_to_eigvectors(x, &edi->vecs.mon, *edi);
    project_to_eigvectors(x, &edi->vecs.linfix, *edi);
    project_to_eigvectors(x, &edi->vecs.linacc, *edi);
    project_to_eigvectors(x, &edi->vecs.radfix, *edi);
    project_to_eigvectors(x, &edi->vecs.radacc, *edi);
    project_to_eigvectors(x, &edi->vecs.radcon, *edi);
}

namespace
{
/*!\brief Evaluates the distance from reference to current eigenvector projection.
 * \param[in] vec eigenvector
 * \returns distance
 */
real calc_radius(const t_eigvec& vec)
{
    real rad = 0.0;

    for (int i = 0; i < vec.neig; i++)
    {
        rad += gmx::square((vec.refproj[i] - vec.xproj[i]));
    }

    return std::sqrt(rad);
}
} // namespace

struct t_do_edfit
{
    double** omega;
    double** om;
};

static void do_edfit(int natoms, rvec* xp, rvec* x, matrix R, t_edpar* edi)
{
    /* this is a copy of do_fit with some modifications */
    int    c, r, n, j, i, irot;
    double d[6], xnr, xpc;
    matrix vh, vk, u;
    int    index;
    real   max_d;

    struct t_do_edfit* loc;
    gmx_bool           bFirst;

    if (edi->buf->do_edfit != nullptr)
    {
        bFirst = FALSE;
    }
    else
    {
        bFirst = TRUE;
        snew(edi->buf->do_edfit, 1);
    }
    loc = edi->buf->do_edfit;

    if (bFirst)
    {
        snew(loc->omega, 2 * DIM);
        snew(loc->om, 2 * DIM);
        for (i = 0; i < 2 * DIM; i++)
        {
            snew(loc->omega[i], 2 * DIM);
            snew(loc->om[i], 2 * DIM);
        }
    }

    for (i = 0; (i < 6); i++)
    {
        d[i] = 0;
        for (j = 0; (j < 6); j++)
        {
            loc->omega[i][j] = 0;
            loc->om[i][j]    = 0;
        }
    }

    /* calculate the matrix U */
    clear_mat(u);
    for (n = 0; (n < natoms); n++)
    {
        for (c = 0; (c < DIM); c++)
        {
            xpc = xp[n][c];
            for (r = 0; (r < DIM); r++)
            {
                xnr = x[n][r];
                u[c][r] += xnr * xpc;
            }
        }
    }

    /* construct loc->omega */
    /* loc->omega is symmetric -> loc->omega==loc->omega' */
    for (r = 0; (r < 6); r++)
    {
        for (c = 0; (c <= r); c++)
        {
            if ((r >= 3) && (c < 3))
            {
                loc->omega[r][c] = u[r - 3][c];
                loc->omega[c][r] = u[r - 3][c];
            }
            else
            {
                loc->omega[r][c] = 0;
                loc->omega[c][r] = 0;
            }
        }
    }

    /* determine h and k */
    jacobi(loc->omega, 6, d, loc->om, &irot);

    if (irot == 0)
    {
        fprintf(stderr, "IROT=0\n");
    }

    index = 0; /* For the compiler only */

    for (j = 0; (j < 3); j++)
    {
        max_d = -1000;
        for (i = 0; (i < 6); i++)
        {
            if (d[i] > max_d)
            {
                max_d = d[i];
                index = i;
            }
        }
        d[index] = -10000;
        for (i = 0; (i < 3); i++)
        {
            vh[j][i] = M_SQRT2 * loc->om[i][index];
            vk[j][i] = M_SQRT2 * loc->om[i + DIM][index];
        }
    }

    /* determine R */
    for (c = 0; (c < 3); c++)
    {
        for (r = 0; (r < 3); r++)
        {
            R[c][r] = vk[0][r] * vh[0][c] + vk[1][r] * vh[1][c] + vk[2][r] * vh[2][c];
        }
    }
    if (det(R) < 0)
    {
        for (c = 0; (c < 3); c++)
        {
            for (r = 0; (r < 3); r++)
            {
                R[c][r] = vk[0][r] * vh[0][c] + vk[1][r] * vh[1][c] - vk[2][r] * vh[2][c];
            }
        }
    }
}


static void rmfit(int nat, rvec* xcoll, const rvec transvec, matrix rotmat)
{
    rvec   vec;
    matrix tmat;


    /* Remove rotation.
     * The inverse rotation is described by the transposed rotation matrix */
    transpose(rotmat, tmat);
    rotate_x(xcoll, nat, tmat);

    /* Remove translation */
    vec[XX] = -transvec[XX];
    vec[YY] = -transvec[YY];
    vec[ZZ] = -transvec[ZZ];
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

// TODO split this into multiple files
namespace
{
/*!\brief Output flooding simulation settings and results to file.
 * \param[in] edi Essential dynamics input parameters
 * \param[in] fp output file
 * \param[in] rmsd rmsd to reference structure
 */
void write_edo_flood(const t_edpar& edi, FILE* fp, real rmsd)
{
    /* Output how well we fit to the reference structure */
    fprintf(fp, EDcol_ffmt, rmsd);

    for (int i = 0; i < edi.flood.vecs.neig; i++)
    {
        fprintf(fp, EDcol_efmt, edi.flood.vecs.xproj[i]);

        /* Check whether the reference projection changes with time (this can happen
         * in case flooding is used as harmonic restraint). If so, output the
         * current reference projection */
        if (edi.flood.bHarmonic && edi.flood.referenceProjectionSlope[i] != 0.0)
        {
            fprintf(fp, EDcol_efmt, edi.flood.vecs.refproj[i]);
        }

        /* Output Efl if we are doing adaptive flooding */
        if (0 != edi.flood.tau)
        {
            fprintf(fp, EDcol_efmt, edi.flood.Efl);
        }
        fprintf(fp, EDcol_efmt, edi.flood.Vfl);

        /* Output deltaF if we are doing adaptive flooding */
        if (0 != edi.flood.tau)
        {
            fprintf(fp, EDcol_efmt, edi.flood.deltaF);
        }
        fprintf(fp, EDcol_efmt, edi.flood.vecs.fproj[i]);
    }
}
} // namespace


/* From flood.xproj compute the Vfl(x) at this point */
static real flood_energy(t_edpar* edi, int64_t step)
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
    int  i;


    /* Each time this routine is called (i.e. each time step), we add a small
     * value to the reference projection. This way a harmonic restraint towards
     * a moving reference is realized. If no value for the additive constant
     * is provided in the edi file, the reference will not change. */
    if (edi->flood.bHarmonic)
    {
        for (i = 0; i < edi->flood.vecs.neig; i++)
        {
            edi->flood.vecs.refproj[i] = edi->flood.initialReferenceProjection[i]
                                         + step * edi->flood.referenceProjectionSlope[i];
        }
    }

    sum = 0.0;
    /* Compute sum which will be the exponent of the exponential */
    for (i = 0; i < edi->flood.vecs.neig; i++)
    {
        /* stpsz stores the reciprocal eigenvalue 1/sigma_i */
        sum += edi->flood.vecs.stpsz[i] * (edi->flood.vecs.xproj[i] - edi->flood.vecs.refproj[i])
               * (edi->flood.vecs.xproj[i] - edi->flood.vecs.refproj[i]);
    }

    /* Compute the Gauss function*/
    if (edi->flood.bHarmonic)
    {
        Vfl = -0.5 * edi->flood.Efl * sum; /* minus sign because Efl is negative, if restrain is on. */
    }
    else
    {
        Vfl = edi->flood.Efl != 0
                      ? edi->flood.Efl
                                * std::exp(-edi->flood.kT / 2 / edi->flood.Efl / edi->flood.alpha2 * sum)
                      : 0;
    }

    return Vfl;
}


/* From the position and from Vfl compute forces in subspace -> store in edi->vec.flood.fproj */
static void flood_forces(t_edpar* edi)
{
    /* compute the forces in the subspace of the flooding eigenvectors
     * by the formula F_i= V_{fl}(c) * ( \frac {kT} {E_{fl}} \lambda_i c_i */

    int  i;
    real energy = edi->flood.Vfl;


    if (edi->flood.bHarmonic)
    {
        for (i = 0; i < edi->flood.vecs.neig; i++)
        {
            edi->flood.vecs.fproj[i] = edi->flood.Efl * edi->flood.vecs.stpsz[i]
                                       * (edi->flood.vecs.xproj[i] - edi->flood.vecs.refproj[i]);
        }
    }
    else
    {
        for (i = 0; i < edi->flood.vecs.neig; i++)
        {
            /* if Efl is zero the forces are zero if not use the formula */
            edi->flood.vecs.fproj[i] =
                    edi->flood.Efl != 0
                            ? edi->flood.kT / edi->flood.Efl / edi->flood.alpha2 * energy
                                      * edi->flood.vecs.stpsz[i]
                                      * (edi->flood.vecs.xproj[i] - edi->flood.vecs.refproj[i])
                            : 0;
        }
    }
}

namespace
{
/*!\brief Raise forces from subspace into cartesian space.
 * This function lifts the forces from the subspace to the cartesian space
 * all the values not contained in the subspace are assumed to be zero and then
 * a coordinate transformation from eigenvector to cartesian vectors is performed
 * The nonexistent values don't have to be set to zero explicitly, they would occur
 * as zero valued summands, hence we just stop to compute this part of the sum.
 * For every atom we add all the contributions to this atom from all the different eigenvectors.
 * NOTE: one could add directly to the forcefield forces, would mean we wouldn't have to clear the
 * field forces_cart prior the computation, but we compute the forces separately
 * to have them accessible for diagnostics
 *
 * \param[in] edi Essential dynamics input parameters
 * \param[out] forces_cart The cartesian forces
 */

void flood_blowup(const t_edpar& edi, rvec* forces_cart)
{
    const real* forces_sub = edi.flood.vecs.fproj;
    /* Calculate the cartesian forces for the local atoms */

    /* Clear forces first */
    for (int j = 0; j < edi.sav.nr_loc; j++)
    {
        clear_rvec(forces_cart[j]);
    }

    /* Now compute atomwise */
    for (int j = 0; j < edi.sav.nr_loc; j++)
    {
        /* Compute forces_cart[edi.sav.anrs[j]] */
        for (int eig = 0; eig < edi.flood.vecs.neig; eig++)
        {
            rvec addedForce;
            /* Force vector is force * eigenvector (compute only atom j) */
            svmul(forces_sub[eig], edi.flood.vecs.vec[eig][edi.sav.c_ind[j]], addedForce);
            /* Add this vector to the cartesian forces */
            rvec_inc(forces_cart[j], addedForce);
        }
    }
}

} // namespace


/* Update the values of Efl, deltaF depending on tau and Vfl */
static void update_adaption(t_edpar* edi)
{
    /* this function updates the parameter Efl and deltaF according to the rules given in
     * 'predicting unimolecular chemical reactions: chemical flooding' M Mueller et al,
     * J. chem Phys. */

    if ((edi->flood.tau < 0 ? -edi->flood.tau : edi->flood.tau) > 0.00000001)
    {
        edi->flood.Efl = edi->flood.Efl
                         + edi->flood.dt / edi->flood.tau * (edi->flood.deltaF0 - edi->flood.deltaF);
        /* check if restrain (inverted flooding) -> don't let EFL become positive */
        if (edi->flood.alpha2 < 0 && edi->flood.Efl > -0.00000001)
        {
            edi->flood.Efl = 0;
        }

        edi->flood.deltaF = (1 - edi->flood.dt / edi->flood.tau) * edi->flood.deltaF
                            + edi->flood.dt / edi->flood.tau * edi->flood.Vfl;
    }
}


static void do_single_flood(FILE*                          edo,
                            gmx::ArrayRef<const gmx::RVec> coords,
                            gmx::ArrayRef<gmx::RVec>       force,
                            t_edpar*                       edi,
                            int64_t                        step,
                            const matrix                   box,
                            const gmx::MpiComm&            mpiComm,
                            gmx_bool bNS) /* Are we in a neighbor searching step? */
{
    int                i;
    matrix             rotmat;   /* rotation matrix */
    matrix             tmat;     /* inverse rotation */
    rvec               transvec; /* translation vector */
    real               rmsdev;
    struct t_do_edsam* buf;


    buf = edi->buf->do_edsam;


    /* Broadcast the positions of the AVERAGE structure such that they are known on
     * every processor. Each node contributes its local positions x and stores them in
     * the collective ED array buf->xcoll */
    communicate_group_positions(mpiComm,
                                buf->xcoll,
                                buf->shifts_xcoll,
                                buf->extra_shifts_xcoll,
                                bNS,
                                as_rvec_array(coords.data()),
                                edi->sav.nr,
                                edi->sav.nr_loc,
                                edi->sav.anrs_loc,
                                edi->sav.c_ind,
                                edi->sav.x_old,
                                box);

    /* Only assembly REFERENCE positions if their indices differ from the average ones */
    if (!edi->bRefEqAv)
    {
        communicate_group_positions(mpiComm,
                                    buf->xc_ref,
                                    buf->shifts_xc_ref,
                                    buf->extra_shifts_xc_ref,
                                    bNS,
                                    as_rvec_array(coords.data()),
                                    edi->sref.nr,
                                    edi->sref.nr_loc,
                                    edi->sref.anrs_loc,
                                    edi->sref.c_ind,
                                    edi->sref.x_old,
                                    box);
    }

    /* If bUpdateShifts was TRUE, the shifts have just been updated in get_positions.
     * We do not need to update the shifts until the next NS step */
    buf->bUpdateShifts = FALSE;

    /* Now all nodes have all of the ED/flooding positions in edi->sav->xcoll,
     * as well as the indices in edi->sav.anrs */

    /* Fit the reference indices to the reference structure */
    if (edi->bRefEqAv)
    {
        fit_to_reference(buf->xcoll, transvec, rotmat, edi);
    }
    else
    {
        fit_to_reference(buf->xc_ref, transvec, rotmat, edi);
    }

    /* Now apply the translation and rotation to the ED structure */
    translate_and_rotate(buf->xcoll, edi->sav.nr, transvec, rotmat);

    /* Project fitted structure onto supbspace -> store in edi->flood.vecs.xproj */
    project_to_eigvectors(buf->xcoll, &edi->flood.vecs, *edi);

    if (!edi->flood.bConstForce)
    {
        /* Compute Vfl(x) from flood.xproj */
        edi->flood.Vfl = flood_energy(edi, step);

        update_adaption(edi);

        /* Compute the flooding forces */
        flood_forces(edi);
    }

    /* Translate them into cartesian positions */
    flood_blowup(*edi, edi->flood.forces_cartesian);

    /* Rotate forces back so that they correspond to the given structure and not to the fitted one */
    /* Each node rotates back its local forces */
    transpose(rotmat, tmat);
    rotate_x(edi->flood.forces_cartesian, edi->sav.nr_loc, tmat);

    /* Finally add forces to the main force variable */
    for (i = 0; i < edi->sav.nr_loc; i++)
    {
        rvec_inc(force[edi->sav.anrs_loc[i]], edi->flood.forces_cartesian[i]);
    }

    /* Output is written by the main process */
    if (do_per_step(step, edi->outfrq) && mpiComm.isMainRank())
    {
        /* Output how well we fit to the reference */
        if (edi->bRefEqAv)
        {
            /* Indices of reference and average structures are identical,
             * thus we can calculate the rmsd to SREF using xcoll */
            rmsdev = rmsd_from_structure(buf->xcoll, &edi->sref);
        }
        else
        {
            /* We have to translate & rotate the reference atoms first */
            translate_and_rotate(buf->xc_ref, edi->sref.nr, transvec, rotmat);
            rmsdev = rmsd_from_structure(buf->xc_ref, &edi->sref);
        }

        write_edo_flood(*edi, edo, rmsdev);
    }
}


/* Main flooding routine, called from do_force */
void do_flood(const gmx::MpiComm&            mpiComm,
              const t_inputrec&              ir,
              gmx::ArrayRef<const gmx::RVec> coords,
              gmx::ArrayRef<gmx::RVec>       force,
              gmx_edsam*                     ed,
              const matrix                   box,
              int64_t                        step,
              bool                           bNS)
{
    /* Write time to edo, when required. Output the time anyhow since we need
     * it in the output file for ED constraints. */
    if (mpiComm.isMainRank() && do_per_step(step, ed->edpar.begin()->outfrq))
    {
        fprintf(ed->edo, "\n%12f", ir.init_t + step * ir.delta_t);
    }

    if (ed->eEDtype != EssentialDynamicsType::Flooding)
    {
        return;
    }

    for (auto& edi : ed->edpar)
    {
        /* Call flooding for one matrix */
        if (edi.flood.vecs.neig)
        {
            do_single_flood(ed->edo, coords, force, &edi, step, box, mpiComm, bNS);
        }
    }
}


/* Called by init_edi, configure some flooding related variables and structures,
 * print headers to output files */
static void init_flood(t_edpar* edi, gmx_edsam* ed, real dt)
{
    int i;


    edi->flood.Efl = edi->flood.constEfl;
    edi->flood.Vfl = 0;
    edi->flood.dt  = dt;

    if (edi->flood.vecs.neig)
    {
        /* If in any of the ED groups we find a flooding vector, flooding is turned on */
        ed->eEDtype = EssentialDynamicsType::Flooding;

        fprintf(stderr,
                "ED: Flooding %d eigenvector%s.\n",
                edi->flood.vecs.neig,
                edi->flood.vecs.neig > 1 ? "s" : "");

        if (edi->flood.bConstForce)
        {
            /* We have used stpsz as a vehicle to carry the fproj values for constant
             * force flooding. Now we copy that to flood.vecs.fproj. Note that
             * in const force flooding, fproj is never changed. */
            for (i = 0; i < edi->flood.vecs.neig; i++)
            {
                edi->flood.vecs.fproj[i] = edi->flood.vecs.stpsz[i];

                fprintf(stderr,
                        "ED: applying on eigenvector %d a constant force of %g\n",
                        edi->flood.vecs.ieig[i],
                        edi->flood.vecs.fproj[i]);
            }
        }
    }
}


#ifdef DEBUGHELPERS
/*********** Energy book keeping ******/
static void get_flood_enx_names(t_edpar* edi, char** names, int* nnames) /* get header of energies */
{
    t_edpar* actual;
    int      count;
    char     buf[STRLEN];
    actual = edi;
    count  = 1;
    while (actual)
    {
        srenew(names, count);
        sprintf(buf, "Vfl_%d", count);
        names[count - 1] = gmx_strdup(buf);
        actual           = actual->next_edi;
        count++;
    }
    *nnames = count - 1;
}


static void get_flood_energies(t_edpar* edi, real Vfl[], int nnames)
{
    /*fl has to be big enough to capture nnames-many entries*/
    t_edpar* actual;
    int      count;


    actual = edi;
    count  = 1;
    while (actual)
    {
        Vfl[count - 1] = actual->flood.Vfl;
        actual         = actual->next_edi;
        count++;
    }
    if (nnames != count - 1)
    {
        gmx_fatal(FARGS, "Number of energies is not consistent with t_edi structure");
    }
}
/************* END of FLOODING IMPLEMENTATION ****************************/
#endif


/* This function opens the ED input and output files, reads in all datasets it finds
 * in the input file, and cross-checks whether the .edi file information is consistent
 * with the essential dynamics data found in the checkpoint file (if present).
 * gmx make_edi can be used to create an .edi input file.
 */
static std::unique_ptr<gmx::EssentialDynamics> ed_open(int                         natoms,
                                                       ObservablesHistory*         oh,
                                                       const char*                 ediFileName,
                                                       const char*                 edoFileName,
                                                       const gmx::StartingBehavior startingBehavior,
                                                       const gmx_output_env_t*     oenv,
                                                       const gmx::MpiComm&         mpiComm)
{
    auto  edHandle = std::make_unique<gmx::EssentialDynamics>();
    auto* ed       = edHandle->getLegacyED();
    /* We want to perform ED (this switch might later be upgraded to EssentialDynamicsType::Flooding) */
    ed->eEDtype = EssentialDynamicsType::EDSampling;

    if (mpiComm.isMainRank())
    {
        // If we start from a checkpoint file, we already have an edsamHistory struct
        if (oh->edsamHistory == nullptr)
        {
            oh->edsamHistory = std::make_unique<edsamhistory_t>(edsamhistory_t{});
        }
        edsamhistory_t* EDstate = oh->edsamHistory.get();

        /* Read the edi input file: */
        ed->edpar = read_edi_file(ediFileName, natoms);

        /* Make sure the checkpoint was produced in a run using this .edi file */
        if (EDstate->bFromCpt)
        {
            crosscheck_edi_file_vs_checkpoint(*ed, EDstate);
        }
        else
        {
            EDstate->nED = ed->edpar.size();
        }
        init_edsamstate(*ed, EDstate);

        /* The main opens the ED output file */
        if (startingBehavior == gmx::StartingBehavior::RestartWithAppending)
        {
            ed->edo = gmx_fio_fopen(edoFileName, "a+");
        }
        else
        {
            ed->edo = xvgropen(edoFileName,
                               "Essential dynamics / flooding output",
                               "Time (ps)",
                               "RMSDs (nm), projections on EVs (nm), ...",
                               oenv);

            /* Make a descriptive legend */
            write_edo_legend(ed, EDstate->nED, oenv);
        }
    }
    return edHandle;
}


/* Broadcasts the structure data */
static void bc_ed_positions(const gmx::MpiComm& mpiComm, struct gmx_edx* s, EssentialDynamicsStructure stype)
{
    snew_bc(mpiComm.isMainRank(), s->anrs, s->nr); /* Index numbers     */
    snew_bc(mpiComm.isMainRank(), s->x, s->nr);    /* Positions         */
    nblock_bc(mpiComm.comm(), s->nr, s->anrs);
    nblock_bc(mpiComm.comm(), s->nr, s->x);

    /* For the average & reference structures we need an array for the collective indices,
     * and we need to broadcast the masses as well */
    if (stype == EssentialDynamicsStructure::Average || stype == EssentialDynamicsStructure::Reference)
    {
        /* We need these additional variables in the parallel case: */
        snew(s->c_ind, s->nr); /* Collective indices */
        /* Local atom indices get assigned in dd_make_local_group_indices.
         * There, also memory is allocated */
        s->nalloc_loc = 0;                              /* allocation size of s->anrs_loc */
        snew_bc(mpiComm.isMainRank(), s->x_old, s->nr); /* To be able to always make the ED molecule whole, ... */
        nblock_bc(mpiComm.comm(), s->nr, s->x_old); /* ... keep track of shift changes with the help of old coords */
    }

    /* broadcast masses for the reference structure (for mass-weighted fitting) */
    if (stype == EssentialDynamicsStructure::Reference)
    {
        snew_bc(mpiComm.isMainRank(), s->m, s->nr);
        nblock_bc(mpiComm.comm(), s->nr, s->m);
    }

    /* For the average structure we might need the masses for mass-weighting */
    if (stype == EssentialDynamicsStructure::Average)
    {
        snew_bc(mpiComm.isMainRank(), s->sqrtm, s->nr);
        nblock_bc(mpiComm.comm(), s->nr, s->sqrtm);
        snew_bc(mpiComm.isMainRank(), s->m, s->nr);
        nblock_bc(mpiComm.comm(), s->nr, s->m);
    }
}


/* Broadcasts the eigenvector data */
static void bc_ed_vecs(const gmx::MpiComm& mpiComm, t_eigvec* ev, int length)
{
    int i;

    snew_bc(mpiComm.isMainRank(), ev->ieig, ev->neig);    /* index numbers of eigenvector  */
    snew_bc(mpiComm.isMainRank(), ev->stpsz, ev->neig);   /* stepsizes per eigenvector     */
    snew_bc(mpiComm.isMainRank(), ev->xproj, ev->neig);   /* instantaneous x projection    */
    snew_bc(mpiComm.isMainRank(), ev->fproj, ev->neig);   /* instantaneous f projection    */
    snew_bc(mpiComm.isMainRank(), ev->refproj, ev->neig); /* starting or target projection */

    nblock_bc(mpiComm.comm(), ev->neig, ev->ieig);
    nblock_bc(mpiComm.comm(), ev->neig, ev->stpsz);
    nblock_bc(mpiComm.comm(), ev->neig, ev->xproj);
    nblock_bc(mpiComm.comm(), ev->neig, ev->fproj);
    nblock_bc(mpiComm.comm(), ev->neig, ev->refproj);

    snew_bc(mpiComm.isMainRank(), ev->vec, ev->neig); /* Eigenvector components        */
    for (i = 0; i < ev->neig; i++)
    {
        snew_bc(mpiComm.isMainRank(), ev->vec[i], length);
        nblock_bc(mpiComm.comm(), length, ev->vec[i]);
    }
}


/* Broadcasts the ED / flooding data to other nodes
 * and allocates memory where needed */
static void broadcast_ed_data(const gmx::MpiComm& mpiComm, gmx_edsam* ed)
{
    /* Main lets the other nodes know if its ED only or also flooding */
    gmx_bcast(sizeof(ed->eEDtype), &(ed->eEDtype), mpiComm.comm());

    int numedis = ed->edpar.size();
    /* First let everybody know how many ED data sets to expect */
    gmx_bcast(sizeof(numedis), &numedis, mpiComm.comm());
    nblock_abc(mpiComm.isMainRank(), mpiComm.comm(), numedis, &(ed->edpar));
    /* Now transfer the ED data set(s) */
    for (auto& edi : ed->edpar)
    {
        /* Broadcast a single ED data set */
        block_bc(mpiComm.comm(), edi);

        /* Broadcast positions */
        bc_ed_positions(mpiComm, &(edi.sref), EssentialDynamicsStructure::Reference); /* reference positions (don't broadcast masses)    */
        bc_ed_positions(mpiComm, &(edi.sav), EssentialDynamicsStructure::Average); /* average positions (do broadcast masses as well) */
        bc_ed_positions(mpiComm, &(edi.star), EssentialDynamicsStructure::Target); /* target positions */
        bc_ed_positions(mpiComm, &(edi.sori), EssentialDynamicsStructure::Origin); /* origin positions */

        /* Broadcast eigenvectors */
        bc_ed_vecs(mpiComm, &edi.vecs.mon, edi.sav.nr);
        bc_ed_vecs(mpiComm, &edi.vecs.linfix, edi.sav.nr);
        bc_ed_vecs(mpiComm, &edi.vecs.linacc, edi.sav.nr);
        bc_ed_vecs(mpiComm, &edi.vecs.radfix, edi.sav.nr);
        bc_ed_vecs(mpiComm, &edi.vecs.radacc, edi.sav.nr);
        bc_ed_vecs(mpiComm, &edi.vecs.radcon, edi.sav.nr);
        /* Broadcast flooding eigenvectors and, if needed, values for the moving reference */
        bc_ed_vecs(mpiComm, &edi.flood.vecs, edi.sav.nr);

        /* For harmonic restraints the reference projections can change with time */
        if (edi.flood.bHarmonic)
        {
            snew_bc(mpiComm.isMainRank(), edi.flood.initialReferenceProjection, edi.flood.vecs.neig);
            snew_bc(mpiComm.isMainRank(), edi.flood.referenceProjectionSlope, edi.flood.vecs.neig);
            nblock_bc(mpiComm.comm(), edi.flood.vecs.neig, edi.flood.initialReferenceProjection);
            nblock_bc(mpiComm.comm(), edi.flood.vecs.neig, edi.flood.referenceProjectionSlope);
        }
    }
}


/* init-routine called for every *.edi-cycle, initialises t_edpar structure */
static void init_edi(const gmx_mtop_t& mtop, t_edpar* edi)
{
    int  i;
    real totalmass = 0.0;
    rvec com;

    /* NOTE Init_edi is executed on the main process only
     * The initialized data sets are then transmitted to the
     * other nodes in broadcast_ed_data */

    /* evaluate masses (reference structure) */
    snew(edi->sref.m, edi->sref.nr);
    int molb = 0;
    for (i = 0; i < edi->sref.nr; i++)
    {
        if (edi->fitmas)
        {
            edi->sref.m[i] = mtopGetAtomMass(mtop, edi->sref.anrs[i], &molb);
        }
        else
        {
            edi->sref.m[i] = 1.0;
        }

        /* Check that every m > 0. Bad things will happen otherwise. */
        if (edi->sref.m[i] <= 0.0)
        {
            gmx_fatal(FARGS,
                      "Reference structure atom %d (sam.edi index %d) has a mass of %g.\n"
                      "For a mass-weighted fit, all reference structure atoms need to have a mass "
                      ">0.\n"
                      "Either make the covariance analysis non-mass-weighted, or exclude massless\n"
                      "atoms from the reference structure by creating a proper index group.\n",
                      i,
                      edi->sref.anrs[i] + 1,
                      edi->sref.m[i]);
        }

        totalmass += edi->sref.m[i];
    }
    edi->sref.mtot = totalmass;

    /* Masses m and sqrt(m) for the average structure. Note that m
     * is needed if forces have to be evaluated in do_edsam */
    snew(edi->sav.sqrtm, edi->sav.nr);
    snew(edi->sav.m, edi->sav.nr);
    for (i = 0; i < edi->sav.nr; i++)
    {
        edi->sav.m[i] = mtopGetAtomMass(mtop, edi->sav.anrs[i], &molb);
        if (edi->pcamas)
        {
            edi->sav.sqrtm[i] = std::sqrt(edi->sav.m[i]);
        }
        else
        {
            edi->sav.sqrtm[i] = 1.0;
        }

        /* Check that every m > 0. Bad things will happen otherwise. */
        if (edi->sav.sqrtm[i] <= 0.0)
        {
            gmx_fatal(FARGS,
                      "Average structure atom %d (sam.edi index %d) has a mass of %g.\n"
                      "For ED with mass-weighting, all average structure atoms need to have a mass "
                      ">0.\n"
                      "Either make the covariance analysis non-mass-weighted, or exclude massless\n"
                      "atoms from the average structure by creating a proper index group.\n",
                      i,
                      edi->sav.anrs[i] + 1,
                      edi->sav.m[i]);
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


static void check(const char* line, const char* label)
{
    if (!std::strstr(line, label))
    {
        gmx_fatal(FARGS,
                  "Could not find input parameter %s at expected position in edsam input-file "
                  "(.edi)\nline read instead is %s",
                  label,
                  line);
    }
}


static int read_checked_edint(FILE* file, const char* label)
{
    char line[STRLEN + 1];
    int  idum;

    fgets2(line, STRLEN, file);
    check(line, label);
    fgets2(line, STRLEN, file);
    sscanf(line, max_ev_fmt_d, &idum);
    return idum;
}

static bool read_checked_edbool(FILE* file, const char* label)
{
    return read_checked_edint(file, label) != 0;
}

static int read_edint(FILE* file, bool* bEOF)
{
    char  line[STRLEN + 1];
    int   idum;
    char* eof;


    eof = fgets2(line, STRLEN, file);
    if (eof == nullptr)
    {
        *bEOF = TRUE;
        return -1;
    }
    eof = fgets2(line, STRLEN, file);
    if (eof == nullptr)
    {
        *bEOF = TRUE;
        return -1;
    }
    sscanf(line, max_ev_fmt_d, &idum);
    *bEOF = FALSE;
    return idum;
}


static real read_checked_edreal(FILE* file, const char* label)
{
    char   line[STRLEN + 1];
    double rdum;


    fgets2(line, STRLEN, file);
    check(line, label);
    fgets2(line, STRLEN, file);
    sscanf(line, max_ev_fmt_lf, &rdum);
    return static_cast<real>(rdum); /* always read as double and convert to single */
}


static void read_edx(FILE* file, int number, int* anrs, rvec* x)
{
    int    i, j;
    char   line[STRLEN + 1];
    double d[3];


    for (i = 0; i < number; i++)
    {
        fgets2(line, STRLEN, file);
        sscanf(line, max_ev_fmt_dlflflf, &anrs[i], &d[0], &d[1], &d[2]);
        anrs[i]--; /* we are reading FORTRAN indices */
        for (j = 0; j < 3; j++)
        {
            x[i][j] = d[j]; /* always read as double and convert to single */
        }
    }
}

namespace
{
/*!\brief Read eigenvector coordinates for multiple eigenvectors from a file.
 * \param[in] in the file to read from
 * \param[in] numAtoms the number of atoms/coordinates for the eigenvector
 * \param[out] vec the eigen-vectors
 * \param[in] nEig the number of eigenvectors
 */
void scan_edvec(FILE* in, int numAtoms, rvec*** vec, int nEig)
{
    snew(*vec, nEig);
    for (int iEigenvector = 0; iEigenvector < nEig; iEigenvector++)
    {
        snew((*vec)[iEigenvector], numAtoms);
        for (int iAtom = 0; iAtom < numAtoms; iAtom++)
        {
            char   line[STRLEN + 1];
            double x, y, z;
            fgets2(line, STRLEN, in);
            sscanf(line, max_ev_fmt_lelele, &x, &y, &z);
            (*vec)[iEigenvector][iAtom][XX] = x;
            (*vec)[iEigenvector][iAtom][YY] = y;
            (*vec)[iEigenvector][iAtom][ZZ] = z;
        }
    }
}
/*!\brief Read an essential dynamics eigenvector and corresponding step size.
 * \param[in] in input file
 * \param[in] nrAtoms number of atoms/coordinates
 * \param[out] tvec the eigenvector
 */
void read_edvec(FILE* in, int nrAtoms, t_eigvec* tvec)
{
    tvec->neig = read_checked_edint(in, "NUMBER OF EIGENVECTORS");
    if (tvec->neig <= 0)
    {
        return;
    }

    snew(tvec->ieig, tvec->neig);
    snew(tvec->stpsz, tvec->neig);
    for (int i = 0; i < tvec->neig; i++)
    {
        char line[STRLEN + 1];
        fgets2(line, STRLEN, in);
        int       numEigenVectors;
        double    stepSize;
        const int nscan = sscanf(line, max_ev_fmt_dlf, &numEigenVectors, &stepSize);
        if (nscan != 2)
        {
            gmx_fatal(FARGS, "Expected 2 values for flooding vec: <nr> <stpsz>\n");
        }
        tvec->ieig[i]  = numEigenVectors;
        tvec->stpsz[i] = stepSize;
    } /* end of loop over eigenvectors */

    scan_edvec(in, nrAtoms, &tvec->vec, tvec->neig);
}
/*!\brief Read an essential dynamics eigenvector and initial reference projections and slopes if available.
 *
 * Eigenvector and its initial reference and slope are stored on the
 * same line in the .edi format. To avoid re-winding the .edi file,
 * the reference eigen-vector and reference data are read in one go.
 *
 * \param[in] in input file
 * \param[in] nrAtoms number of atoms/coordinates
 * \param[out] tvec the eigenvector
 * \param[out] initialReferenceProjection The projection onto eigenvectors as read from file.
 * \param[out] referenceProjectionSlope The slope of the reference projections.
 */
bool readEdVecWithReferenceProjection(FILE*     in,
                                      int       nrAtoms,
                                      t_eigvec* tvec,
                                      real**    initialReferenceProjection,
                                      real**    referenceProjectionSlope)
{
    bool providesReference = false;
    tvec->neig             = read_checked_edint(in, "NUMBER OF EIGENVECTORS");
    if (tvec->neig <= 0)
    {
        return false;
    }

    snew(tvec->ieig, tvec->neig);
    snew(tvec->stpsz, tvec->neig);
    snew(*initialReferenceProjection, tvec->neig);
    snew(*referenceProjectionSlope, tvec->neig);
    for (int i = 0; i < tvec->neig; i++)
    {
        char line[STRLEN + 1];
        fgets2(line, STRLEN, in);
        int       numEigenVectors;
        double    stepSize            = 0.;
        double    referenceProjection = 0.;
        double    referenceSlope      = 0.;
        const int nscan               = sscanf(
                line, max_ev_fmt_dlflflf, &numEigenVectors, &stepSize, &referenceProjection, &referenceSlope);
        /* Zero out values which were not scanned */
        switch (nscan)
        {
            case 4:
                /* Every 4 values read, including reference position */
                providesReference = true;
                break;
            case 3:
                /* A reference position is provided */
                providesReference = true;
                /* No value for slope, set to 0 */
                referenceSlope = 0.0;
                break;
            case 2:
                /* No values for reference projection and slope, set to 0 */
                referenceProjection = 0.0;
                referenceSlope      = 0.0;
                break;
            default:
                gmx_fatal(FARGS,
                          "Expected 2 - 4 (not %d) values for flooding vec: <nr> <spring const> "
                          "<refproj> <refproj-slope>\n",
                          nscan);
        }
        (*initialReferenceProjection)[i] = referenceProjection;
        (*referenceProjectionSlope)[i]   = referenceSlope;

        tvec->ieig[i]  = numEigenVectors;
        tvec->stpsz[i] = stepSize;
    } /* end of loop over eigenvectors */

    scan_edvec(in, nrAtoms, &tvec->vec, tvec->neig);
    return providesReference;
}

/*!\brief Allocate working memory for the eigenvectors.
 * \param[in,out] tvec eigenvector for which memory will be allocated
 */
void setup_edvec(t_eigvec* tvec)
{
    snew(tvec->xproj, tvec->neig);
    snew(tvec->fproj, tvec->neig);
    snew(tvec->refproj, tvec->neig);
}
} // namespace


/* Check if the same atom indices are used for reference and average positions */
static gmx_bool check_if_same(struct gmx_edx sref, struct gmx_edx sav)
{
    int i;


    /* If the number of atoms differs between the two structures,
     * they cannot be identical */
    if (sref.nr != sav.nr)
    {
        return FALSE;
    }

    /* Now that we know that both structures have the same number of atoms,
     * check if also the indices are identical */
    for (i = 0; i < sav.nr; i++)
    {
        if (sref.anrs[i] != sav.anrs[i])
        {
            return FALSE;
        }
    }
    fprintf(stderr,
            "ED: Note: Reference and average structure are composed of the same atom indices.\n");

    return TRUE;
}

namespace
{

//! Oldest (lowest) magic number indicating unsupported essential dynamics input
constexpr int c_oldestUnsupportedMagicNumber = 666;
//! Oldest (lowest) magic number indicating supported essential dynamics input
constexpr int c_oldestSupportedMagicNumber = 669;
//! Latest (highest) magic number indicating supported essential dynamics input
constexpr int c_latestSupportedMagicNumber = 670;

/*!\brief Set up essential dynamics work parameters.
 * \param[in] edi Essential dynamics working structure.
 */
void setup_edi(t_edpar* edi)
{
    snew(edi->sref.x_old, edi->sref.nr);
    edi->sref.sqrtm = nullptr;
    snew(edi->sav.x_old, edi->sav.nr);
    if (edi->star.nr > 0)
    {
        edi->star.sqrtm = nullptr;
    }
    if (edi->sori.nr > 0)
    {
        edi->sori.sqrtm = nullptr;
    }
    setup_edvec(&edi->vecs.linacc);
    setup_edvec(&edi->vecs.mon);
    setup_edvec(&edi->vecs.linfix);
    setup_edvec(&edi->vecs.linacc);
    setup_edvec(&edi->vecs.radfix);
    setup_edvec(&edi->vecs.radacc);
    setup_edvec(&edi->vecs.radcon);
    setup_edvec(&edi->flood.vecs);
}

/*!\brief Check if edi file version supports CONST_FORCE_FLOODING.
 * \param[in] magicNumber indicates the essential dynamics input file version
 * \returns true if CONST_FORCE_FLOODING is supported
 */
bool ediFileHasConstForceFlooding(int magicNumber)
{
    return magicNumber > c_oldestSupportedMagicNumber;
};

/*!\brief Reports reading success of the essential dynamics input file magic number.
 * \param[in] in pointer to input file
 * \param[out] magicNumber determines the edi file version
 * \returns true if setting file version from input file was successful.
 */
bool ediFileHasValidDataSet(FILE* in, int* magicNumber)
{
    /* the edi file is not free format, so expect problems if the input is corrupt. */
    bool reachedEndOfFile = true;
    /* check the magic number */
    *magicNumber = read_edint(in, &reachedEndOfFile);
    /* Check whether we have reached the end of the input file */
    return !reachedEndOfFile;
}

/*!\brief Trigger fatal error if magic number is unsupported.
 * \param[in] magicNumber A magic number read from the edi file
 * \param[in] fn name of input file for error reporting
 */
void exitWhenMagicNumberIsInvalid(int magicNumber, const char* fn)
{

    if (magicNumber < c_oldestSupportedMagicNumber || magicNumber > c_latestSupportedMagicNumber)
    {
        if (magicNumber >= c_oldestUnsupportedMagicNumber && magicNumber < c_oldestSupportedMagicNumber)
        {
            gmx_fatal(FARGS,
                      "Wrong magic number: Use newest version of make_edi to produce edi file");
        }
        else
        {
            gmx_fatal(FARGS, "Wrong magic number %d in %s", magicNumber, fn);
        }
    }
}

/*!\brief reads an essential dynamics input file into a essential dynamics data structure.
 *
 * \param[in] in input file
 * \param[in] nr_mdatoms the number of md atoms, used for consistency checking
 * \param[in] hasConstForceFlooding determines if field "CONST_FORCE_FLOODING" is available in input file
 * \param[in] fn the file name of the input file for error reporting
 * \returns edi essential dynamics parameters
 */
t_edpar read_edi(FILE* in, int nr_mdatoms, bool hasConstForceFlooding, const char* fn)
{
    t_edpar edi;
    bool    bEOF;
    /* check the number of atoms */
    edi.nini = read_edint(in, &bEOF);
    if (edi.nini != nr_mdatoms)
    {
        gmx_fatal(FARGS, "Nr of atoms in %s (%d) does not match nr of md atoms (%d)", fn, edi.nini, nr_mdatoms);
    }

    /* Done checking. For the rest we blindly trust the input */
    edi.fitmas     = read_checked_edbool(in, "FITMAS");
    edi.pcamas     = read_checked_edbool(in, "ANALYSIS_MAS");
    edi.outfrq     = read_checked_edint(in, "OUTFRQ");
    edi.maxedsteps = read_checked_edint(in, "MAXLEN");
    edi.slope      = read_checked_edreal(in, "SLOPECRIT");

    edi.presteps        = read_checked_edint(in, "PRESTEPS");
    edi.flood.deltaF0   = read_checked_edreal(in, "DELTA_F0");
    edi.flood.deltaF    = read_checked_edreal(in, "INIT_DELTA_F");
    edi.flood.tau       = read_checked_edreal(in, "TAU");
    edi.flood.constEfl  = read_checked_edreal(in, "EFL_NULL");
    edi.flood.alpha2    = read_checked_edreal(in, "ALPHA2");
    edi.flood.kT        = read_checked_edreal(in, "KT");
    edi.flood.bHarmonic = read_checked_edbool(in, "HARMONIC");
    if (hasConstForceFlooding)
    {
        edi.flood.bConstForce = read_checked_edbool(in, "CONST_FORCE_FLOODING");
    }
    else
    {
        edi.flood.bConstForce = FALSE;
    }
    edi.sref.nr = read_checked_edint(in, "NREF");

    /* allocate space for reference positions and read them */
    snew(edi.sref.anrs, edi.sref.nr);
    snew(edi.sref.x, edi.sref.nr);
    read_edx(in, edi.sref.nr, edi.sref.anrs, edi.sref.x);

    /* average positions. they define which atoms will be used for ED sampling */
    edi.sav.nr = read_checked_edint(in, "NAV");
    snew(edi.sav.anrs, edi.sav.nr);
    snew(edi.sav.x, edi.sav.nr);
    read_edx(in, edi.sav.nr, edi.sav.anrs, edi.sav.x);

    /* Check if the same atom indices are used for reference and average positions */
    edi.bRefEqAv = check_if_same(edi.sref, edi.sav);

    /* eigenvectors */

    read_edvec(in, edi.sav.nr, &edi.vecs.mon);
    read_edvec(in, edi.sav.nr, &edi.vecs.linfix);
    read_edvec(in, edi.sav.nr, &edi.vecs.linacc);
    read_edvec(in, edi.sav.nr, &edi.vecs.radfix);
    read_edvec(in, edi.sav.nr, &edi.vecs.radacc);
    read_edvec(in, edi.sav.nr, &edi.vecs.radcon);

    /* Was a specific reference point for the flooding/umbrella potential provided in the edi file? */
    bool bHaveReference = false;
    if (edi.flood.bHarmonic)
    {
        bHaveReference = readEdVecWithReferenceProjection(in,
                                                          edi.sav.nr,
                                                          &edi.flood.vecs,
                                                          &edi.flood.initialReferenceProjection,
                                                          &edi.flood.referenceProjectionSlope);
    }
    else
    {
        read_edvec(in, edi.sav.nr, &edi.flood.vecs);
    }

    /* target positions */
    edi.star.nr = read_edint(in, &bEOF);
    if (edi.star.nr > 0)
    {
        snew(edi.star.anrs, edi.star.nr);
        snew(edi.star.x, edi.star.nr);
        read_edx(in, edi.star.nr, edi.star.anrs, edi.star.x);
    }

    /* positions defining origin of expansion circle */
    edi.sori.nr = read_edint(in, &bEOF);
    if (edi.sori.nr > 0)
    {
        if (bHaveReference)
        {
            /* Both an -ori structure and a at least one manual reference point have been
             * specified. That's ambiguous and probably not intentional. */
            gmx_fatal(FARGS,
                      "ED: An origin structure has been provided and a at least one (moving) "
                      "reference\n"
                      "    point was manually specified in the edi file. That is ambiguous. "
                      "Aborting.\n");
        }
        snew(edi.sori.anrs, edi.sori.nr);
        snew(edi.sori.x, edi.sori.nr);
        read_edx(in, edi.sori.nr, edi.sori.anrs, edi.sori.x);
    }
    return edi;
}

std::vector<t_edpar> read_edi_file(const char* fn, int nr_mdatoms)
{
    std::vector<t_edpar> essentialDynamicsParameters;
    FILE*                in;
    /* This routine is executed on the main only */

    /* Open the .edi parameter input file */
    in = gmx_fio_fopen(fn, "r");
    fprintf(stderr, "ED: Reading edi file %s\n", fn);

    /* Now read a sequence of ED input parameter sets from the edi file */
    int ediFileMagicNumber;
    while (ediFileHasValidDataSet(in, &ediFileMagicNumber))
    {
        exitWhenMagicNumberIsInvalid(ediFileMagicNumber, fn);
        bool hasConstForceFlooding = ediFileHasConstForceFlooding(ediFileMagicNumber);
        auto edi                   = read_edi(in, nr_mdatoms, hasConstForceFlooding, fn);
        setup_edi(&edi);
        essentialDynamicsParameters.emplace_back(edi);

        /* Since we arrived within this while loop we know that there is still another data set to be read in */
        /* We need to allocate space for the data: */
    }
    if (essentialDynamicsParameters.empty())
    {
        gmx_fatal(FARGS, "No complete ED data set found in edi file %s.", fn);
    }

    fprintf(stderr,
            "ED: Found %zu ED group%s.\n",
            essentialDynamicsParameters.size(),
            essentialDynamicsParameters.size() > 1 ? "s" : "");

    /* Close the .edi file again */
    gmx_fio_fclose(in);

    return essentialDynamicsParameters;
}
} // anonymous namespace


struct t_fit_to_ref
{
    rvec* xcopy; /* Working copy of the positions in fit_to_reference */
};

/* Fit the current positions to the reference positions
 * Do not actually do the fit, just return rotation and translation.
 * Note that the COM of the reference structure was already put into
 * the origin by init_edi. */
static void fit_to_reference(rvec*    xcoll,    /* The positions to be fitted */
                             rvec     transvec, /* The translation vector */
                             matrix   rotmat,   /* The rotation matrix */
                             t_edpar* edi)      /* Just needed for do_edfit */
{
    rvec                 com; /* center of mass */
    int                  i;
    struct t_fit_to_ref* loc;


    /* Allocate memory the first time this routine is called for each edi group */
    if (nullptr == edi->buf->fit_to_ref)
    {
        snew(edi->buf->fit_to_ref, 1);
        snew(edi->buf->fit_to_ref->xcopy, edi->sref.nr);
    }
    loc = edi->buf->fit_to_ref;

    /* We do not touch the original positions but work on a copy. */
    for (i = 0; i < edi->sref.nr; i++)
    {
        copy_rvec(xcoll[i], loc->xcopy[i]);
    }

    /* Calculate the center of mass */
    get_center(loc->xcopy, edi->sref.m, edi->sref.nr, com);

    transvec[XX] = -com[XX];
    transvec[YY] = -com[YY];
    transvec[ZZ] = -com[ZZ];

    /* Subtract the center of mass from the copy */
    translate_x(loc->xcopy, edi->sref.nr, transvec);

    /* Determine the rotation matrix */
    do_edfit(edi->sref.nr, edi->sref.x, loc->xcopy, rotmat, edi);
}


static void translate_and_rotate(rvec*  x,        /* The positions to be translated and rotated */
                                 int    nat,      /* How many positions are there? */
                                 rvec   transvec, /* The translation vector */
                                 matrix rotmat)   /* The rotation matrix */
{
    /* Translation */
    translate_x(x, nat, transvec);

    /* Rotation */
    rotate_x(x, nat, rotmat);
}


/* Gets the rms deviation of the positions to the structure s */
/* fit_to_structure has to be called before calling this routine! */
static real rmsd_from_structure(rvec* x,           /* The positions under consideration */
                                struct gmx_edx* s) /* The structure from which the rmsd shall be computed */
{
    real rmsd = 0.0;
    int  i;


    for (i = 0; i < s->nr; i++)
    {
        rmsd += distance2(s->x[i], x[i]);
    }

    rmsd /= static_cast<real>(s->nr);
    rmsd = std::sqrt(rmsd);

    return rmsd;
}


void dd_make_local_ed_indices(gmx_domdec_t* dd, struct gmx_edsam* ed)
{
    if (ed->eEDtype != EssentialDynamicsType::None)
    {
        /* Loop over ED groups */

        for (auto& edi : ed->edpar)
        {
            /* Local atoms of the reference structure (for fitting), need only be assembled
             * if their indices differ from the average ones */
            if (!edi.bRefEqAv)
            {
                dd_make_local_group_indices(dd->ga2la.get(),
                                            edi.sref.nr,
                                            edi.sref.anrs,
                                            &edi.sref.nr_loc,
                                            &edi.sref.anrs_loc,
                                            &edi.sref.nalloc_loc,
                                            edi.sref.c_ind);
            }

            /* Local atoms of the average structure (on these ED will be performed) */
            dd_make_local_group_indices(dd->ga2la.get(),
                                        edi.sav.nr,
                                        edi.sav.anrs,
                                        &edi.sav.nr_loc,
                                        &edi.sav.anrs_loc,
                                        &edi.sav.nalloc_loc,
                                        edi.sav.c_ind);

            /* Indicate that the ED shift vectors for this structure need to be updated
             * at the next call to communicate_group_positions, since obviously we are in a NS step */
            edi.buf->do_edsam->bUpdateShifts = TRUE;
        }
    }
}


static inline void ed_unshift_single_coord(const matrix box, const rvec x, const ivec is, rvec xu)
{
    int tx, ty, tz;


    tx = is[XX];
    ty = is[YY];
    tz = is[ZZ];

    if (TRICLINIC(box))
    {
        xu[XX] = x[XX] - tx * box[XX][XX] - ty * box[YY][XX] - tz * box[ZZ][XX];
        xu[YY] = x[YY] - ty * box[YY][YY] - tz * box[ZZ][YY];
        xu[ZZ] = x[ZZ] - tz * box[ZZ][ZZ];
    }
    else
    {
        xu[XX] = x[XX] - tx * box[XX][XX];
        xu[YY] = x[YY] - ty * box[YY][YY];
        xu[ZZ] = x[ZZ] - tz * box[ZZ][ZZ];
    }
}

namespace
{
/*!\brief Apply fixed linear constraints to essential dynamics variable.
 * \param[in,out] xcoll The collected coordinates.
 * \param[in] edi the essential dynamics parameters
 * \param[in] step the current simulation step
 */
void do_linfix(rvec* xcoll, const t_edpar& edi, int64_t step)
{
    /* loop over linfix vectors */
    for (int i = 0; i < edi.vecs.linfix.neig; i++)
    {
        /* calculate the projection */
        real proj = projectx(edi, xcoll, edi.vecs.linfix.vec[i]);

        /* calculate the correction */
        real preFactor = edi.vecs.linfix.refproj[i] + step * edi.vecs.linfix.stpsz[i] - proj;

        /* apply the correction */
        preFactor /= edi.sav.sqrtm[i];
        for (int j = 0; j < edi.sav.nr; j++)
        {
            rvec differenceVector;
            svmul(preFactor, edi.vecs.linfix.vec[i][j], differenceVector);
            rvec_inc(xcoll[j], differenceVector);
        }
    }
}

/*!\brief Apply acceptance linear constraints to essential dynamics variable.
 * \param[in,out] xcoll The collected coordinates.
 * \param[in] edi the essential dynamics parameters
 */
void do_linacc(rvec* xcoll, t_edpar* edi)
{
    /* loop over linacc vectors */
    for (int i = 0; i < edi->vecs.linacc.neig; i++)
    {
        /* calculate the projection */
        real proj = projectx(*edi, xcoll, edi->vecs.linacc.vec[i]);

        /* calculate the correction */
        real preFactor = 0.0;
        if (edi->vecs.linacc.stpsz[i] > 0.0)
        {
            if ((proj - edi->vecs.linacc.refproj[i]) < 0.0)
            {
                preFactor = edi->vecs.linacc.refproj[i] - proj;
            }
        }
        if (edi->vecs.linacc.stpsz[i] < 0.0)
        {
            if ((proj - edi->vecs.linacc.refproj[i]) > 0.0)
            {
                preFactor = edi->vecs.linacc.refproj[i] - proj;
            }
        }

        /* apply the correction */
        preFactor /= edi->sav.sqrtm[i];
        for (int j = 0; j < edi->sav.nr; j++)
        {
            rvec differenceVector;
            svmul(preFactor, edi->vecs.linacc.vec[i][j], differenceVector);
            rvec_inc(xcoll[j], differenceVector);
        }

        /* new positions will act as reference */
        edi->vecs.linacc.refproj[i] = proj + preFactor;
    }
}
} // namespace


static void do_radfix(rvec* xcoll, t_edpar* edi)
{
    int   i, j;
    real *proj, rad = 0.0, ratio;
    rvec  vec_dum;


    if (edi->vecs.radfix.neig == 0)
    {
        return;
    }

    snew(proj, edi->vecs.radfix.neig);

    /* loop over radfix vectors */
    for (i = 0; i < edi->vecs.radfix.neig; i++)
    {
        /* calculate the projections, radius */
        proj[i] = projectx(*edi, xcoll, edi->vecs.radfix.vec[i]);
        rad += gmx::square(proj[i] - edi->vecs.radfix.refproj[i]);
    }

    rad   = std::sqrt(rad);
    ratio = (edi->vecs.radfix.stpsz[0] + edi->vecs.radfix.radius) / rad - 1.0;
    edi->vecs.radfix.radius += edi->vecs.radfix.stpsz[0];

    /* loop over radfix vectors */
    for (i = 0; i < edi->vecs.radfix.neig; i++)
    {
        proj[i] -= edi->vecs.radfix.refproj[i];

        /* apply the correction */
        proj[i] /= edi->sav.sqrtm[i];
        proj[i] *= ratio;
        for (j = 0; j < edi->sav.nr; j++)
        {
            svmul(proj[i], edi->vecs.radfix.vec[i][j], vec_dum);
            rvec_inc(xcoll[j], vec_dum);
        }
    }

    sfree(proj);
}


static void do_radacc(rvec* xcoll, t_edpar* edi)
{
    int   i, j;
    real *proj, rad = 0.0, ratio = 0.0;
    rvec  vec_dum;


    if (edi->vecs.radacc.neig == 0)
    {
        return;
    }

    snew(proj, edi->vecs.radacc.neig);

    /* loop over radacc vectors */
    for (i = 0; i < edi->vecs.radacc.neig; i++)
    {
        /* calculate the projections, radius */
        proj[i] = projectx(*edi, xcoll, edi->vecs.radacc.vec[i]);
        rad += gmx::square(proj[i] - edi->vecs.radacc.refproj[i]);
    }
    rad = std::sqrt(rad);

    /* only correct when radius decreased */
    if (rad < edi->vecs.radacc.radius)
    {
        ratio = edi->vecs.radacc.radius / rad - 1.0;
    }
    else
    {
        edi->vecs.radacc.radius = rad;
    }

    /* loop over radacc vectors */
    for (i = 0; i < edi->vecs.radacc.neig; i++)
    {
        proj[i] -= edi->vecs.radacc.refproj[i];

        /* apply the correction */
        proj[i] /= edi->sav.sqrtm[i];
        proj[i] *= ratio;
        for (j = 0; j < edi->sav.nr; j++)
        {
            svmul(proj[i], edi->vecs.radacc.vec[i][j], vec_dum);
            rvec_inc(xcoll[j], vec_dum);
        }
    }
    sfree(proj);
}


struct t_do_radcon
{
    real* proj;
};

static void do_radcon(rvec* xcoll, t_edpar* edi)
{
    int                 i, j;
    real                rad = 0.0, ratio = 0.0;
    struct t_do_radcon* loc;
    gmx_bool            bFirst;
    rvec                vec_dum;


    if (edi->buf->do_radcon != nullptr)
    {
        bFirst = FALSE;
    }
    else
    {
        bFirst = TRUE;
        snew(edi->buf->do_radcon, 1);
    }
    loc = edi->buf->do_radcon;

    if (edi->vecs.radcon.neig == 0)
    {
        return;
    }

    if (bFirst)
    {
        snew(loc->proj, edi->vecs.radcon.neig);
    }

    /* loop over radcon vectors */
    for (i = 0; i < edi->vecs.radcon.neig; i++)
    {
        /* calculate the projections, radius */
        loc->proj[i] = projectx(*edi, xcoll, edi->vecs.radcon.vec[i]);
        rad += gmx::square(loc->proj[i] - edi->vecs.radcon.refproj[i]);
    }
    rad = std::sqrt(rad);
    /* only correct when radius increased */
    if (rad > edi->vecs.radcon.radius)
    {
        ratio = edi->vecs.radcon.radius / rad - 1.0;

        /* loop over radcon vectors */
        for (i = 0; i < edi->vecs.radcon.neig; i++)
        {
            /* apply the correction */
            loc->proj[i] -= edi->vecs.radcon.refproj[i];
            loc->proj[i] /= edi->sav.sqrtm[i];
            loc->proj[i] *= ratio;

            for (j = 0; j < edi->sav.nr; j++)
            {
                svmul(loc->proj[i], edi->vecs.radcon.vec[i][j], vec_dum);
                rvec_inc(xcoll[j], vec_dum);
            }
        }
    }
    else
    {
        edi->vecs.radcon.radius = rad;
    }
}


static void ed_apply_constraints(rvec* xcoll, t_edpar* edi, int64_t step)
{
    int i;


    /* subtract the average positions */
    for (i = 0; i < edi->sav.nr; i++)
    {
        rvec_dec(xcoll[i], edi->sav.x[i]);
    }

    /* apply the constraints */
    if (step >= 0)
    {
        do_linfix(xcoll, *edi, step);
    }
    do_linacc(xcoll, edi);
    if (step >= 0)
    {
        do_radfix(xcoll, edi);
    }
    do_radacc(xcoll, edi);
    do_radcon(xcoll, edi);

    /* add back the average positions */
    for (i = 0; i < edi->sav.nr; i++)
    {
        rvec_inc(xcoll[i], edi->sav.x[i]);
    }
}


namespace
{
/*!\brief Write out the projections onto the eigenvectors.
 * The order of output corresponds to ed_output_legend().
 * \param[in] edi The essential dyanmics parameters
 * \param[in] fp The output file
 * \param[in] rmsd the rmsd to the reference structure
 */
void write_edo(const t_edpar& edi, FILE* fp, real rmsd)
{
    /* Output how well we fit to the reference structure */
    fprintf(fp, EDcol_ffmt, rmsd);

    for (int i = 0; i < edi.vecs.mon.neig; i++)
    {
        fprintf(fp, EDcol_efmt, edi.vecs.mon.xproj[i]);
    }

    for (int i = 0; i < edi.vecs.linfix.neig; i++)
    {
        fprintf(fp, EDcol_efmt, edi.vecs.linfix.xproj[i]);
    }

    for (int i = 0; i < edi.vecs.linacc.neig; i++)
    {
        fprintf(fp, EDcol_efmt, edi.vecs.linacc.xproj[i]);
    }

    for (int i = 0; i < edi.vecs.radfix.neig; i++)
    {
        fprintf(fp, EDcol_efmt, edi.vecs.radfix.xproj[i]);
    }
    if (edi.vecs.radfix.neig)
    {
        fprintf(fp, EDcol_ffmt, calc_radius(edi.vecs.radfix)); /* fixed increment radius */
    }

    for (int i = 0; i < edi.vecs.radacc.neig; i++)
    {
        fprintf(fp, EDcol_efmt, edi.vecs.radacc.xproj[i]);
    }
    if (edi.vecs.radacc.neig)
    {
        fprintf(fp, EDcol_ffmt, calc_radius(edi.vecs.radacc)); /* acceptance radius */
    }

    for (int i = 0; i < edi.vecs.radcon.neig; i++)
    {
        fprintf(fp, EDcol_efmt, edi.vecs.radcon.xproj[i]);
    }
    if (edi.vecs.radcon.neig)
    {
        fprintf(fp, EDcol_ffmt, calc_radius(edi.vecs.radcon)); /* contracting radius */
    }
}
} // namespace


/* Returns if any constraints are switched on */
static bool ed_constraints(EssentialDynamicsType edtype, const t_edpar& edi)
{
    if (edtype == EssentialDynamicsType::EDSampling || edtype == EssentialDynamicsType::Flooding)
    {
        return ((edi.vecs.linfix.neig != 0) || (edi.vecs.linacc.neig != 0) || (edi.vecs.radfix.neig != 0)
                || (edi.vecs.radacc.neig != 0) || (edi.vecs.radcon.neig != 0));
    }
    return false;
}


/* Copies reference projection 'refproj' to fixed 'initialReferenceProjection' variable for
 * flooding/ umbrella sampling simulations. */
static void copyEvecReference(t_eigvec* floodvecs, real* initialReferenceProjection)
{
    int i;


    if (nullptr == initialReferenceProjection)
    {
        snew(initialReferenceProjection, floodvecs->neig);
    }

    for (i = 0; i < floodvecs->neig; i++)
    {
        initialReferenceProjection[i] = floodvecs->refproj[i];
    }
}


/* Call on MAIN only. Check whether the essential dynamics / flooding
 * groups of the checkpoint file are consistent with the provided .edi file. */
static void crosscheck_edi_file_vs_checkpoint(const gmx_edsam& ed, edsamhistory_t* EDstate)
{
    if (nullptr == EDstate->nref || nullptr == EDstate->nav)
    {
        gmx_fatal(FARGS,
                  "Essential dynamics and flooding can only be switched on (or off) at the\n"
                  "start of a new simulation. If a simulation runs with/without ED constraints,\n"
                  "it must also continue with/without ED constraints when checkpointing.\n"
                  "To switch on (or off) ED constraints, please prepare a new .tpr to start\n"
                  "from without a checkpoint.\n");
    }

    for (size_t edinum = 0; edinum < ed.edpar.size(); ++edinum)
    {
        /* Check number of atoms in the reference and average structures */
        if (EDstate->nref[edinum] != ed.edpar[edinum].sref.nr)
        {
            gmx_fatal(FARGS,
                      "The number of reference structure atoms in ED group %c is\n"
                      "not the same in .cpt (NREF=%d) and .edi (NREF=%d) files!\n",
                      get_EDgroupChar(edinum + 1, 0),
                      EDstate->nref[edinum],
                      ed.edpar[edinum].sref.nr);
        }
        if (EDstate->nav[edinum] != ed.edpar[edinum].sav.nr)
        {
            gmx_fatal(FARGS,
                      "The number of average structure atoms in ED group %c is\n"
                      "not the same in .cpt (NREF=%d) and .edi (NREF=%d) files!\n",
                      get_EDgroupChar(edinum + 1, 0),
                      EDstate->nav[edinum],
                      ed.edpar[edinum].sav.nr);
        }
    }

    if (gmx::ssize(ed.edpar) != EDstate->nED)
    {
        gmx_fatal(FARGS,
                  "The number of essential dynamics / flooding groups is not consistent.\n"
                  "There are %d ED groups in the .cpt file, but %zu in the .edi file!\n"
                  "Are you sure this is the correct .edi file?\n",
                  EDstate->nED,
                  ed.edpar.size());
    }
}


/* The edsamstate struct stores the information we need to make the ED group
 * whole again after restarts from a checkpoint file. Here we do the following:
 * a) If we did not start from .cpt, we prepare the struct for proper .cpt writing,
 * b) if we did start from .cpt, we copy over the last whole structures from .cpt,
 * c) in any case, for subsequent checkpoint writing, we set the pointers in
 * edsamstate to the x_old arrays, which contain the correct PBC representation of
 * all ED structures at the last time step. */
static void init_edsamstate(const gmx_edsam& ed, edsamhistory_t* EDstate)
{
    snew(EDstate->old_sref_p, EDstate->nED);
    snew(EDstate->old_sav_p, EDstate->nED);

    /* If we did not read in a .cpt file, these arrays are not yet allocated */
    if (!EDstate->bFromCpt)
    {
        snew(EDstate->nref, EDstate->nED);
        snew(EDstate->nav, EDstate->nED);
    }

    /* Loop over all ED/flooding data sets (usually only one, though) */
    for (int nr_edi = 0; nr_edi < EDstate->nED; nr_edi++)
    {
        const auto& edi = ed.edpar[nr_edi];
        /* We always need the last reference and average positions such that
         * in the next time step we can make the ED group whole again
         * if the atoms do not have the correct PBC representation */
        if (EDstate->bFromCpt)
        {
            /* Copy the last whole positions of reference and average group from .cpt */
            for (int i = 0; i < edi.sref.nr; i++)
            {
                copy_rvec(EDstate->old_sref[nr_edi][i], edi.sref.x_old[i]);
            }
            for (int i = 0; i < edi.sav.nr; i++)
            {
                copy_rvec(EDstate->old_sav[nr_edi][i], edi.sav.x_old[i]);
            }
        }
        else
        {
            EDstate->nref[nr_edi] = edi.sref.nr;
            EDstate->nav[nr_edi]  = edi.sav.nr;
        }

        /* For subsequent checkpoint writing, set the edsamstate pointers to the edi arrays: */
        EDstate->old_sref_p[nr_edi] = edi.sref.x_old;
        EDstate->old_sav_p[nr_edi]  = edi.sav.x_old;
    }
}

static void nice_legend(std::vector<std::string>* setname,
                        std::string*              LegendStr,
                        const std::string&        value,
                        const std::string&        unit,
                        char                      EDgroupchar)
{
    auto tmp = gmx::formatString("%c %s", EDgroupchar, value.c_str());
    LegendStr->append(gmx::formatString(EDcol_sfmt, tmp.c_str()));
    tmp += gmx::formatString(" (%s)", unit.c_str());
    setname->emplace_back(tmp);
}


static void nice_legend_evec(std::vector<std::string>* setname,
                             std::string*              LegendStr,
                             t_eigvec*                 evec,
                             char                      EDgroupChar,
                             const std::string&        EDtype)
{
    for (int i = 0; i < evec->neig; i++)
    {
        auto tmp = gmx::formatString("EV%dprj%s", evec->ieig[i], EDtype.c_str());
        nice_legend(setname, LegendStr, tmp, "nm", EDgroupChar);
    }
}


/* Makes a legend for the xvg output file. Call on MAIN only! */
static void write_edo_legend(gmx_edsam* ed, int nED, const gmx_output_env_t* oenv)
{
    int                      n_flood, n_edsam;
    std::vector<std::string> setname;
    std::string              LegendStr;


    auto edi = ed->edpar.begin();

    fprintf(ed->edo, "# Output will be written every %d step%s\n", edi->outfrq, edi->outfrq != 1 ? "s" : "");

    for (int nr_edi = 1; nr_edi <= nED; nr_edi++)
    {
        fprintf(ed->edo, "#\n");
        fprintf(ed->edo,
                "# Summary of applied con/restraints for the ED group %c\n",
                get_EDgroupChar(nr_edi, nED));
        fprintf(ed->edo, "# Atoms in average structure: %d\n", edi->sav.nr);
        fprintf(ed->edo, "#    monitor  : %d vec%s\n", edi->vecs.mon.neig, edi->vecs.mon.neig != 1 ? "s" : "");
        fprintf(ed->edo,
                "#    LINFIX   : %d vec%s\n",
                edi->vecs.linfix.neig,
                edi->vecs.linfix.neig != 1 ? "s" : "");
        fprintf(ed->edo,
                "#    LINACC   : %d vec%s\n",
                edi->vecs.linacc.neig,
                edi->vecs.linacc.neig != 1 ? "s" : "");
        fprintf(ed->edo,
                "#    RADFIX   : %d vec%s\n",
                edi->vecs.radfix.neig,
                edi->vecs.radfix.neig != 1 ? "s" : "");
        fprintf(ed->edo,
                "#    RADACC   : %d vec%s\n",
                edi->vecs.radacc.neig,
                edi->vecs.radacc.neig != 1 ? "s" : "");
        fprintf(ed->edo,
                "#    RADCON   : %d vec%s\n",
                edi->vecs.radcon.neig,
                edi->vecs.radcon.neig != 1 ? "s" : "");
        fprintf(ed->edo, "#    FLOODING : %d vec%s  ", edi->flood.vecs.neig, edi->flood.vecs.neig != 1 ? "s" : "");

        if (edi->flood.vecs.neig)
        {
            /* If in any of the groups we find a flooding vector, flooding is turned on */
            ed->eEDtype = EssentialDynamicsType::Flooding;

            /* Print what flavor of flooding we will do */
            if (0 == edi->flood.tau) /* constant flooding strength */
            {
                fprintf(ed->edo, "Efl_null = %g", edi->flood.constEfl);
                if (edi->flood.bHarmonic)
                {
                    fprintf(ed->edo, ", harmonic");
                }
            }
            else /* adaptive flooding */
            {
                fprintf(ed->edo, ", adaptive");
            }
        }
        fprintf(ed->edo, "\n");

        ++edi;
    }

    LegendStr.append(gmx::formatString("#     %6s", "time"));

    /* In the mdrun time step in a first function call (do_flood()) the flooding
     * forces are calculated and in a second function call (do_edsam()) the
     * ED constraints. To get a corresponding legend, we need to loop twice
     * over the edi groups and output first the flooding, then the ED part */

    /* The flooding-related legend entries, if flooding is done */
    if (EssentialDynamicsType::Flooding == ed->eEDtype)
    {
        auto edi = ed->edpar.begin();
        for (int nr_edi = 1; nr_edi <= nED; nr_edi++)
        {
            /* Always write out the projection on the flooding EVs. Of course, this can also
             * be achieved with the monitoring option in do_edsam() (if switched on by the
             * user), but in that case the positions need to be communicated in do_edsam(),
             * which is not necessary when doing flooding only. */
            nice_legend(&setname, &LegendStr, "RMSD to ref", "nm", get_EDgroupChar(nr_edi, nED));

            for (int i = 0; i < edi->flood.vecs.neig; i++)
            {
                auto buf = gmx::formatString("EV%dprjFLOOD", edi->flood.vecs.ieig[i]);
                nice_legend(&setname, &LegendStr, buf, "nm", get_EDgroupChar(nr_edi, nED));

                /* Output the current reference projection if it changes with time;
                 * this can happen when flooding is used as harmonic restraint */
                if (edi->flood.bHarmonic && edi->flood.referenceProjectionSlope[i] != 0.0)
                {
                    auto buf = gmx::formatString("EV%d ref.prj.", edi->flood.vecs.ieig[i]);
                    nice_legend(&setname, &LegendStr, buf, "nm", get_EDgroupChar(nr_edi, nED));
                }

                /* For flooding we also output Efl, Vfl, deltaF, and the flooding forces */
                if (0 != edi->flood.tau) /* only output Efl for adaptive flooding (constant otherwise) */
                {
                    auto buf = gmx::formatString("EV%d-Efl", edi->flood.vecs.ieig[i]);
                    nice_legend(&setname, &LegendStr, buf, "kJ/mol", get_EDgroupChar(nr_edi, nED));
                }

                buf = gmx::formatString("EV%d-Vfl", edi->flood.vecs.ieig[i]);
                nice_legend(&setname, &LegendStr, buf, "kJ/mol", get_EDgroupChar(nr_edi, nED));

                if (0 != edi->flood.tau) /* only output deltaF for adaptive flooding (zero otherwise) */
                {
                    auto buf = gmx::formatString("EV%d-deltaF", edi->flood.vecs.ieig[i]);
                    nice_legend(&setname, &LegendStr, buf, "kJ/mol", get_EDgroupChar(nr_edi, nED));
                }

                buf = gmx::formatString("EV%d-FLforces", edi->flood.vecs.ieig[i]);
                nice_legend(&setname, &LegendStr, buf, "kJ/mol/nm", get_EDgroupChar(nr_edi, nED));
            }

            ++edi;
        } /* End of flooding-related legend entries */
    }
    n_flood = setname.size();

    /* Now the ED-related entries, if essential dynamics is done */
    edi = ed->edpar.begin();
    for (int nr_edi = 1; nr_edi <= nED; nr_edi++)
    {
        if (bNeedDoEdsam(*edi)) /* Only print ED legend if at least one ED option is on */
        {
            nice_legend(&setname, &LegendStr, "RMSD to ref", "nm", get_EDgroupChar(nr_edi, nED));

            /* Essential dynamics, projections on eigenvectors */
            nice_legend_evec(
                    &setname, &LegendStr, &edi->vecs.mon, get_EDgroupChar(nr_edi, nED), "MON");
            nice_legend_evec(
                    &setname, &LegendStr, &edi->vecs.linfix, get_EDgroupChar(nr_edi, nED), "LINFIX");
            nice_legend_evec(
                    &setname, &LegendStr, &edi->vecs.linacc, get_EDgroupChar(nr_edi, nED), "LINACC");
            nice_legend_evec(
                    &setname, &LegendStr, &edi->vecs.radfix, get_EDgroupChar(nr_edi, nED), "RADFIX");
            if (edi->vecs.radfix.neig)
            {
                nice_legend(&setname, &LegendStr, "RADFIX radius", "nm", get_EDgroupChar(nr_edi, nED));
            }
            nice_legend_evec(
                    &setname, &LegendStr, &edi->vecs.radacc, get_EDgroupChar(nr_edi, nED), "RADACC");
            if (edi->vecs.radacc.neig)
            {
                nice_legend(&setname, &LegendStr, "RADACC radius", "nm", get_EDgroupChar(nr_edi, nED));
            }
            nice_legend_evec(
                    &setname, &LegendStr, &edi->vecs.radcon, get_EDgroupChar(nr_edi, nED), "RADCON");
            if (edi->vecs.radcon.neig)
            {
                nice_legend(&setname, &LegendStr, "RADCON radius", "nm", get_EDgroupChar(nr_edi, nED));
            }
        }
        ++edi;
    } /* end of 'pure' essential dynamics legend entries */
    n_edsam = setname.size() - n_flood;

    xvgrLegend(ed->edo, setname, oenv);

    fprintf(ed->edo,
            "#\n"
            "# Legend for %d column%s of flooding plus %d column%s of essential dynamics data:\n",
            n_flood,
            1 == n_flood ? "" : "s",
            n_edsam,
            1 == n_edsam ? "" : "s");
    fprintf(ed->edo, "%s", LegendStr.c_str());

    std::fflush(ed->edo);
}


/* Init routine for ED and flooding. Calls init_edi in a loop for every .edi-cycle
 * contained in the input file, creates a NULL terminated list of t_edpar structures */
std::unique_ptr<gmx::EssentialDynamics> init_edsam(const gmx::MDLogger&        mdlog,
                                                   const char*                 ediFileName,
                                                   const char*                 edoFileName,
                                                   const gmx_mtop_t&           mtop,
                                                   const t_inputrec&           ir,
                                                   const gmx::MpiComm&         mpiComm,
                                                   const gmx_domdec_t*         dd,
                                                   gmx::Constraints*           constr,
                                                   const t_state*              globalState,
                                                   ObservablesHistory*         oh,
                                                   const gmx_output_env_t*     oenv,
                                                   const gmx::StartingBehavior startingBehavior)
{
    int    i, avindex;
    rvec * xfit = nullptr, *xstart = nullptr; /* dummy arrays to determine initial RMSDs  */
    rvec   fit_transvec;                      /* translation ... */
    matrix fit_rotmat;                        /* ... and rotation from fit to reference structure */
    rvec*  ref_x_old = nullptr;               /* helper pointer */


    if (mpiComm.isMainRank())
    {
        fprintf(stderr, "ED: Initializing essential dynamics constraints.\n");

        if (oh->edsamHistory != nullptr && oh->edsamHistory->bFromCpt && !gmx_fexist(ediFileName))
        {
            gmx_fatal(FARGS,
                      "The checkpoint file you provided is from an essential dynamics or flooding\n"
                      "simulation. Please also set the .edi file on the command line with -ei.\n");
        }
    }
    GMX_LOG(mdlog.info)
            .asParagraph()
            .appendText(
                    "Activating essential dynamics via files passed to "
                    "gmx mdrun -edi will change in future to be files passed to "
                    "gmx grompp and the related .mdp options may change also.");

    /* Open input and output files, allocate space for ED data structure */
    auto edHandle = ed_open(mtop.natoms, oh, ediFileName, edoFileName, startingBehavior, oenv, mpiComm);
    auto* ed = edHandle->getLegacyED();
    GMX_RELEASE_ASSERT(constr != nullptr, "Must have valid constraints object");
    constr->saveEdsamPointer(ed);

    /* Needed for initializing radacc radius in do_edsam */
    ed->bFirst = TRUE;

    /* The input file is read by the main and the edi structures are
     * initialized here. Input is stored in ed->edpar. Then the edi
     * structures are transferred to the other nodes */
    if (mpiComm.isMainRank())
    {
        /* Initialization for every ED/flooding group. Flooding uses one edi group per
         * flooding vector, Essential dynamics can be applied to more than one structure
         * as well, but will be done in the order given in the edi file, so
         * expect different results for different order of edi file concatenation! */
        for (auto& edi : ed->edpar)
        {
            init_edi(mtop, &edi);
            init_flood(&edi, ed, ir.delta_t);
        }
    }

    std::vector<gmx::RVec> x_pbc; /* positions of the whole MD system with pbc removed  */

    /* The main does the work here. The other nodes get the positions
     * not before dd_partition_system which is called after init_edsam */
    if (mpiComm.isMainRank())
    {
        const edsamhistory_t* EDstate = oh->edsamHistory.get();

        if (!EDstate->bFromCpt)
        {
            /* Remove PBC, make molecule(s) subject to ED whole. */
            x_pbc.resize(mtop.natoms);
            std::copy(globalState->x.begin(), globalState->x.end(), x_pbc.begin());
            do_pbc_first_mtop(nullptr, ir.pbcType, false, nullptr, globalState->box, &mtop, x_pbc, {});
        }
        /* Reset pointer to first ED data set which contains the actual ED data */
        auto edi = ed->edpar.begin();
        GMX_ASSERT(!ed->edpar.empty(), "Essential dynamics parameters is undefined");

        /* Loop over all ED/flooding data sets (usually only one, though) */
        for (int nr_edi = 1; nr_edi <= EDstate->nED; nr_edi++)
        {
            /* For multiple ED groups we use the output frequency that was specified
             * in the first set */
            if (nr_edi > 1)
            {
                edi->outfrq = ed->edpar.begin()->outfrq;
            }

            /* Extract the initial reference and average positions. When starting
             * from .cpt, these have already been read into sref.x_old
             * in init_edsamstate() */
            if (!EDstate->bFromCpt)
            {
                /* If this is the first run (i.e. no checkpoint present) we assume
                 * that the starting positions give us the correct PBC representation */
                for (i = 0; i < edi->sref.nr; i++)
                {
                    copy_rvec(x_pbc[edi->sref.anrs[i]], edi->sref.x_old[i]);
                }

                for (i = 0; i < edi->sav.nr; i++)
                {
                    copy_rvec(x_pbc[edi->sav.anrs[i]], edi->sav.x_old[i]);
                }
            }

            /* Now we have the PBC-correct start positions of the reference and
               average structure. We copy that over to dummy arrays on which we
               can apply fitting to print out the RMSD. We srenew the memory since
               the size of the buffers is likely different for every ED group */
            srenew(xfit, edi->sref.nr);
            srenew(xstart, edi->sav.nr);
            if (edi->bRefEqAv)
            {
                /* Reference indices are the same as average indices */
                ref_x_old = edi->sav.x_old;
            }
            else
            {
                ref_x_old = edi->sref.x_old;
            }
            copy_rvecn(ref_x_old, xfit, 0, edi->sref.nr);
            copy_rvecn(edi->sav.x_old, xstart, 0, edi->sav.nr);

            /* Make the fit to the REFERENCE structure, get translation and rotation */
            fit_to_reference(xfit, fit_transvec, fit_rotmat, &(*edi));

            /* Output how well we fit to the reference at the start */
            translate_and_rotate(xfit, edi->sref.nr, fit_transvec, fit_rotmat);
            fprintf(stderr,
                    "ED: Initial RMSD from reference after fit = %f nm",
                    rmsd_from_structure(xfit, &edi->sref));
            if (EDstate->nED > 1)
            {
                fprintf(stderr, " (ED group %c)", get_EDgroupChar(nr_edi, EDstate->nED));
            }
            fprintf(stderr, "\n");

            /* Now apply the translation and rotation to the atoms on which ED sampling will be performed */
            translate_and_rotate(xstart, edi->sav.nr, fit_transvec, fit_rotmat);

            /* calculate initial projections */
            project(xstart, &(*edi));

            /* For the target and origin structure both a reference (fit) and an
             * average structure can be provided in make_edi. If both structures
             * are the same, make_edi only stores one of them in the .edi file.
             * If they differ, first the fit and then the average structure is stored
             * in star (or sor), thus the number of entries in star/sor is
             * (n_fit + n_av) with n_fit the size of the fitting group and n_av
             * the size of the average group. */

            /* process target structure, if required */
            if (edi->star.nr > 0)
            {
                fprintf(stderr, "ED: Fitting target structure to reference structure\n");

                /* get translation & rotation for fit of target structure to reference structure */
                fit_to_reference(edi->star.x, fit_transvec, fit_rotmat, &(*edi));
                /* do the fit */
                translate_and_rotate(edi->star.x, edi->star.nr, fit_transvec, fit_rotmat);
                if (edi->star.nr == edi->sav.nr)
                {
                    avindex = 0;
                }
                else /* edi->star.nr = edi->sref.nr + edi->sav.nr */
                {
                    /* The last sav.nr indices of the target structure correspond to
                     * the average structure, which must be projected */
                    avindex = edi->star.nr - edi->sav.nr;
                }
                rad_project(*edi, &edi->star.x[avindex], &edi->vecs.radcon);
            }
            else
            {
                rad_project(*edi, xstart, &edi->vecs.radcon);
            }

            /* process structure that will serve as origin of expansion circle */
            if ((EssentialDynamicsType::Flooding == ed->eEDtype) && (!edi->flood.bConstForce))
            {
                fprintf(stderr,
                        "ED: Setting center of flooding potential (0 = average structure)\n");
            }

            if (edi->sori.nr > 0)
            {
                fprintf(stderr, "ED: Fitting origin structure to reference structure\n");

                /* fit this structure to reference structure */
                fit_to_reference(edi->sori.x, fit_transvec, fit_rotmat, &(*edi));
                /* do the fit */
                translate_and_rotate(edi->sori.x, edi->sori.nr, fit_transvec, fit_rotmat);
                if (edi->sori.nr == edi->sav.nr)
                {
                    avindex = 0;
                }
                else /* edi->sori.nr = edi->sref.nr + edi->sav.nr */
                {
                    /* For the projection, we need the last sav.nr indices of sori */
                    avindex = edi->sori.nr - edi->sav.nr;
                }

                rad_project(*edi, &edi->sori.x[avindex], &edi->vecs.radacc);
                rad_project(*edi, &edi->sori.x[avindex], &edi->vecs.radfix);
                if ((EssentialDynamicsType::Flooding == ed->eEDtype) && (!edi->flood.bConstForce))
                {
                    fprintf(stderr,
                            "ED: The ORIGIN structure will define the flooding potential "
                            "center.\n");
                    /* Set center of flooding potential to the ORIGIN structure */
                    rad_project(*edi, &edi->sori.x[avindex], &edi->flood.vecs);
                    /* We already know that no (moving) reference position was provided,
                     * therefore we can overwrite refproj[0]*/
                    copyEvecReference(&edi->flood.vecs, edi->flood.initialReferenceProjection);
                }
            }
            else /* No origin structure given */
            {
                rad_project(*edi, xstart, &edi->vecs.radacc);
                rad_project(*edi, xstart, &edi->vecs.radfix);
                if ((EssentialDynamicsType::Flooding == ed->eEDtype) && (!edi->flood.bConstForce))
                {
                    if (edi->flood.bHarmonic)
                    {
                        fprintf(stderr,
                                "ED: A (possibly changing) ref. projection will define the "
                                "flooding potential center.\n");
                        for (i = 0; i < edi->flood.vecs.neig; i++)
                        {
                            edi->flood.vecs.refproj[i] = edi->flood.initialReferenceProjection[i];
                        }
                    }
                    else
                    {
                        fprintf(stderr,
                                "ED: The AVERAGE structure will define the flooding potential "
                                "center.\n");
                        /* Set center of flooding potential to the center of the covariance matrix,
                         * i.e. the average structure, i.e. zero in the projected system */
                        for (i = 0; i < edi->flood.vecs.neig; i++)
                        {
                            edi->flood.vecs.refproj[i] = 0.0;
                        }
                    }
                }
            }
            /* For convenience, output the center of the flooding potential for the eigenvectors */
            if ((EssentialDynamicsType::Flooding == ed->eEDtype) && (!edi->flood.bConstForce))
            {
                for (i = 0; i < edi->flood.vecs.neig; i++)
                {
                    fprintf(stdout,
                            "ED: EV %d flooding potential center: %11.4e",
                            edi->flood.vecs.ieig[i],
                            edi->flood.vecs.refproj[i]);
                    if (edi->flood.bHarmonic)
                    {
                        fprintf(stdout, " (adding %11.4e/timestep)", edi->flood.referenceProjectionSlope[i]);
                    }
                    fprintf(stdout, "\n");
                }
            }

            /* set starting projections for linsam */
            rad_project(*edi, xstart, &edi->vecs.linacc);
            rad_project(*edi, xstart, &edi->vecs.linfix);

            /* Prepare for the next edi data set: */
            ++edi;
        }
        /* Cleaning up on the main node: */
        sfree(xfit);
        sfree(xstart);

    } /* end of MAIN only section */

    if (dd)
    {
        /* Broadcast the essential dynamics / flooding data to all nodes.
         * In a single-rank case, only the necessary memory allocation is done. */
        broadcast_ed_data(mpiComm, ed);
    }
    else
    {
        /* In the non-DD case, point the local atom numbers pointers to the global
         * one, so that we can use the same notation in serial and parallel case: */
        /* Loop over all ED data sets (usually only one, though) */
        for (auto edi = ed->edpar.begin(); edi != ed->edpar.end(); ++edi)
        {
            edi->sref.anrs_loc = edi->sref.anrs;
            edi->sav.anrs_loc  = edi->sav.anrs;
            edi->star.anrs_loc = edi->star.anrs;
            edi->sori.anrs_loc = edi->sori.anrs;
            /* For the same reason as above, make a dummy c_ind array: */
            snew(edi->sav.c_ind, edi->sav.nr);
            /* Initialize the array */
            for (i = 0; i < edi->sav.nr; i++)
            {
                edi->sav.c_ind[i] = i;
            }
            /* In the general case we will need a different-sized array for the reference indices: */
            if (!edi->bRefEqAv)
            {
                snew(edi->sref.c_ind, edi->sref.nr);
                for (i = 0; i < edi->sref.nr; i++)
                {
                    edi->sref.c_ind[i] = i;
                }
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
    }

    /* Allocate space for ED buffer variables */
    /* Again, loop over ED data sets */
    for (auto edi = ed->edpar.begin(); edi != ed->edpar.end(); ++edi)
    {
        /* Allocate space for ED buffer variables */
        snew_bc(mpiComm.isMainRank(), edi->buf, 1); /* MAIN has already allocated edi->buf in init_edi() */
        snew(edi->buf->do_edsam, 1);

        /* Space for collective ED buffer variables */

        /* Collective positions of atoms with the average indices */
        snew(edi->buf->do_edsam->xcoll, edi->sav.nr);
        snew(edi->buf->do_edsam->shifts_xcoll, edi->sav.nr); /* buffer for xcoll shifts */
        snew(edi->buf->do_edsam->extra_shifts_xcoll, edi->sav.nr);
        /* Collective positions of atoms with the reference indices */
        if (!edi->bRefEqAv)
        {
            snew(edi->buf->do_edsam->xc_ref, edi->sref.nr);
            snew(edi->buf->do_edsam->shifts_xc_ref, edi->sref.nr); /* To store the shifts in */
            snew(edi->buf->do_edsam->extra_shifts_xc_ref, edi->sref.nr);
        }

        /* Get memory for flooding forces */
        snew(edi->flood.forces_cartesian, edi->sav.nr);
    }

    /* Flush the edo file so that the user can check some things
     * when the simulation has started */
    if (ed->edo)
    {
        std::fflush(ed->edo);
    }

    return edHandle;
}


void do_edsam(const t_inputrec*        ir,
              int64_t                  step,
              const gmx::MpiComm&      mpiComm,
              gmx::ArrayRef<gmx::RVec> coords,
              gmx::ArrayRef<gmx::RVec> velocities,
              const matrix             box,
              gmx_edsam*               ed)
{
    int    i, iupdate = 500;
    matrix rotmat;         /* rotation matrix */
    rvec   transvec;       /* translation vector */
    rvec   dv, dx, x_unsh; /* tmp vectors for velocity, distance, unshifted x coordinate */
    real   dt_1;           /* 1/dt */
    struct t_do_edsam* buf;
    real     rmsdev    = -1; /* RMSD from reference structure prior to applying the constraints */
    gmx_bool bSuppress = FALSE; /* Write .xvg output file on main? */


    /* Check if ED sampling has to be performed */
    if (ed->eEDtype == EssentialDynamicsType::None)
    {
        return;
    }

    dt_1 = 1.0 / ir->delta_t;

    /* Loop over all ED groups (usually one) */
    for (auto& edi : ed->edpar)
    {
        if (bNeedDoEdsam(edi))
        {

            buf = edi.buf->do_edsam;

            if (ed->bFirst)
            {
                /* initialize radacc radius for slope criterion */
                buf->oldrad = calc_radius(edi.vecs.radacc);
            }

            /* Copy the positions into buf->xc* arrays and after ED
             * feed back corrections to the official positions */

            /* Broadcast the ED positions such that every node has all of them
             * Every node contributes its local positions xs and stores it in
             * the collective buf->xcoll array. Note that for edinr > 1
             * xs could already have been modified by an earlier ED */

            communicate_group_positions(mpiComm,
                                        buf->xcoll,
                                        buf->shifts_xcoll,
                                        buf->extra_shifts_xcoll,
                                        mpiComm.size() > 1 ? buf->bUpdateShifts : TRUE,
                                        as_rvec_array(coords.data()),
                                        edi.sav.nr,
                                        edi.sav.nr_loc,
                                        edi.sav.anrs_loc,
                                        edi.sav.c_ind,
                                        edi.sav.x_old,
                                        box);

            /* Only assembly reference positions if their indices differ from the average ones */
            if (!edi.bRefEqAv)
            {
                communicate_group_positions(mpiComm,
                                            buf->xc_ref,
                                            buf->shifts_xc_ref,
                                            buf->extra_shifts_xc_ref,
                                            mpiComm.size() > 1 ? buf->bUpdateShifts : TRUE,
                                            as_rvec_array(coords.data()),
                                            edi.sref.nr,
                                            edi.sref.nr_loc,
                                            edi.sref.anrs_loc,
                                            edi.sref.c_ind,
                                            edi.sref.x_old,
                                            box);
            }

            /* If bUpdateShifts was TRUE then the shifts have just been updated in communicate_group_positions.
             * We do not need to update the shifts until the next NS step. Note that dd_make_local_ed_indices
             * set bUpdateShifts=TRUE in the parallel case. */
            buf->bUpdateShifts = FALSE;

            /* Now all nodes have all of the ED positions in edi.sav->xcoll,
             * as well as the indices in edi.sav.anrs */

            /* Fit the reference indices to the reference structure */
            if (edi.bRefEqAv)
            {
                fit_to_reference(buf->xcoll, transvec, rotmat, &edi);
            }
            else
            {
                fit_to_reference(buf->xc_ref, transvec, rotmat, &edi);
            }

            /* Now apply the translation and rotation to the ED structure */
            translate_and_rotate(buf->xcoll, edi.sav.nr, transvec, rotmat);

            /* Find out how well we fit to the reference (just for output steps) */
            if (do_per_step(step, edi.outfrq) && mpiComm.isMainRank())
            {
                if (edi.bRefEqAv)
                {
                    /* Indices of reference and average structures are identical,
                     * thus we can calculate the rmsd to SREF using xcoll */
                    rmsdev = rmsd_from_structure(buf->xcoll, &edi.sref);
                }
                else
                {
                    /* We have to translate & rotate the reference atoms first */
                    translate_and_rotate(buf->xc_ref, edi.sref.nr, transvec, rotmat);
                    rmsdev = rmsd_from_structure(buf->xc_ref, &edi.sref);
                }
            }

            /* update radsam references, when required */
            if (do_per_step(step, edi.maxedsteps) && step >= edi.presteps)
            {
                project(buf->xcoll, &edi);
                rad_project(edi, buf->xcoll, &edi.vecs.radacc);
                rad_project(edi, buf->xcoll, &edi.vecs.radfix);
                buf->oldrad = -1.e5;
            }

            /* update radacc references, when required */
            if (do_per_step(step, iupdate) && step >= edi.presteps)
            {
                edi.vecs.radacc.radius = calc_radius(edi.vecs.radacc);
                if (edi.vecs.radacc.radius - buf->oldrad < edi.slope)
                {
                    project(buf->xcoll, &edi);
                    rad_project(edi, buf->xcoll, &edi.vecs.radacc);
                    buf->oldrad = 0.0;
                }
                else
                {
                    buf->oldrad = edi.vecs.radacc.radius;
                }
            }

            /* apply the constraints */
            if (step >= edi.presteps && ed_constraints(ed->eEDtype, edi))
            {
                /* ED constraints should be applied already in the first MD step
                 * (which is step 0), therefore we pass step+1 to the routine */
                ed_apply_constraints(buf->xcoll, &edi, step + 1 - ir->init_step);
            }

            /* write to edo, when required */
            if (do_per_step(step, edi.outfrq))
            {
                project(buf->xcoll, &edi);
                if (mpiComm.isMainRank() && !bSuppress)
                {
                    write_edo(edi, ed->edo, rmsdev);
                }
            }

            /* Copy back the positions unless monitoring only */
            if (ed_constraints(ed->eEDtype, edi))
            {
                /* remove fitting */
                rmfit(edi.sav.nr, buf->xcoll, transvec, rotmat);

                /* Copy the ED corrected positions into the coordinate array */
                /* Each node copies its local part. In the serial case, nat_loc is the
                 * total number of ED atoms */
                for (i = 0; i < edi.sav.nr_loc; i++)
                {
                    /* Unshift local ED coordinate and store in x_unsh */
                    ed_unshift_single_coord(
                            box, buf->xcoll[edi.sav.c_ind[i]], buf->shifts_xcoll[edi.sav.c_ind[i]], x_unsh);

                    /* dx is the ED correction to the positions: */
                    rvec_sub(x_unsh, coords[edi.sav.anrs_loc[i]], dx);

                    if (!velocities.empty())
                    {
                        /* dv is the ED correction to the velocity: */
                        svmul(dt_1, dx, dv);
                        /* apply the velocity correction: */
                        rvec_inc(velocities[edi.sav.anrs_loc[i]], dv);
                    }
                    /* Finally apply the position correction due to ED: */
                    copy_rvec(x_unsh, coords[edi.sav.anrs_loc[i]]);
                }
            }
        } /* END of if ( bNeedDoEdsam(edi) ) */

        /* Prepare for the next ED group */

    } /* END of loop over ED groups */

    ed->bFirst = FALSE;
}
