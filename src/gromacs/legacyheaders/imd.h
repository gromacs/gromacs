#ifndef IMD_H__
#define IMD_H__

#ifdef GMX_IMD

/*how long shall we wait until we check for a connection again?*/
#define IMDLOOPWAIT 0.1

/*how long shall we check for the IMD_GO?.*/
#define IMDCONNECTWAIT 2

/*IMD default port*/
#define IMDDEFPORT 8888

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include <Windows.h>
#define NOFLAGS 0
#endif

/* Put this before all Gromacs/IMD-related output */
const static char IMDstr[] = "IMD:";

#include <limits.h>

/*GROMACS Includes*/
#include "typedefs.h"
#include "imdsocket.h"
#include "readinp.h"
#include "types/commrec.h"
#include "oenv.h"
#include "filenm.h"

/*We define int32 which is the 32bit integer used in the IMD functions.*/
#if ( INT_MAX == 2147483647 )
typedef int int32;
#else
typedef short int32;
#endif

/*IMD Protocol Version*/
#define HEADERSIZE 8
#define IMDVERSION 2

/*We use the same records as the NAMD/VMD IMD implementation.*/
typedef enum IMDType_t
{
    IMD_DISCONNECT, /*client disconnect*/
    IMD_ENERGIES, /*energy data*/
    IMD_FCOORDS, /*atomic coordinates*/
    IMD_GO, /*start command for the simulation*/
    IMD_HANDSHAKE, /*handshake to determine little/big endian*/
    IMD_KILL, /*terminates the simulation*/
    IMD_MDCOMM, /*force data*/
    IMD_PAUSE, /*pause the simulation*/
    IMD_TRATE, /*sets the IMD transmission, and processing rate*/
    IMD_IOERROR, /*I/O error*/
    IMD_NR
/*<-GROMACS specific extension to access message names*/
} IMDMessageType;

/*Macros to access names for the IMDMessageType*/
#define UNDEFINED       "UNDEFINED"
#define ENUM_NAME(e,max,names)  ((((e)<0)||((e)>=(max)))?UNDEFINED:(names)[e])

/*Energy record as in the original IMD implementation, energys in Kcal/mol*/
/*NOTE: We return the energies in GROMACS / SI units, so they also show up as SI in VMD.*/
typedef struct
{
    int32 tstep; /*Timestep*/
    float T_abs; /*Absolute Temperature*/
    float E_tot; /*Total Energy*/
    float E_pot; /*Potential Energy*/
    float E_vdw; /*Van der Waals Energy*/
    float E_coul; /*Coulomb Interaction Energy*/
    float E_bond; /*Bond Energy*/
    float E_angle; /*Angle Energy*/
    float E_dihe; /*Dihedral Energy*/
    float E_impr; /*Improper Dihedral Energy*/
} IMDEnergyBlock;

/* definitions for the IMD header & protocol version */
typedef struct
{
    int32 type;
    int32 length;
} IMDHeader;

/* Public functions, see imd.c */

void write_imdatoms(t_inputrec *ir, t_state *state, gmx_mtop_t *sys, const char *fn);

void dd_make_local_IMD_atoms(gmx_domdec_t *dd, t_IMD *imd);

extern void init_imd(
        t_inputrec *ir,
        t_commrec *cr,
        gmx_mtop_t *top_global,
        FILE *fplog,
        int defnstimd,
        int nat_total,
        rvec x[],
        t_mdatoms *md,
        int nfile,
        const t_filenm fnm[],
        const output_env_t oenv,
        int imdport,
        int imdfreq,
        unsigned long Flags);

void imd_send_x_E(gmx_mtop_t *mtop, rvec *x, matrix box, t_inputrec *ir);

void do_imd_prepare_energies(t_gmx_IMD IMDSetup, gmx_enerdata_t *enerd, int step, gmx_bool bHaveNewEnergies);

/*returns if we should do IMD communication in this step. Also checks for new IMD connection and syncs the nodes.*/
extern gmx_bool do_IMD(
        int          step,
        t_commrec   *cr,
        gmx_bool     bNS,       /* Is this a ns step?          */
        matrix       box,
        rvec         x[],       /* The atomic positions local to this node    */
        t_inputrec  *ir,
        double       t);

int imd_get_step(t_gmx_IMD IMDSetup);

extern void do_imd_send_positions(
        t_IMD *imd);


void imd_apply_forces(t_inputrec *ir, t_commrec *cr, rvec *f);
void imd_finalize(t_inputrec *ir);

#endif
#endif
