#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_IMD

#include "imd.h"
#include "imdsocket.h"
#include "smalloc.h"
#include "network.h"
#include "mdrun.h"
#include "sighandler.h"
#include "gmx_ga2la.h"
#include "xvgr.h"
#include "partdec.h"
#include "groupcoord.h"
#include "mtop_util.h"
#include "confio.h"
#include "vec.h"

#include <string.h>

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include <Windows.h>
#else
#include <unistd.h>
#endif

#define  block_bc(cr,   d) gmx_bcast(     sizeof(d),     &(d),(cr))
#define  nblock_bc(cr,nr,d) gmx_bcast((nr)*sizeof((d)[0]), (d),(cr))

/* Local IMD data */
typedef struct gmx_IMD
{
    FILE      *outf;            /* Output file for IMD data, mainly forces    */

    int        nat;             /* Number of atoms that can be pulled via IMD */
    int        nat_loc;         /* Part of the atoms that are local           */
    atom_id   *ind;             /* Global indices of the IMD atoms            */
    atom_id   *ind_loc;         /* Local indices of the IMD atoms             */
    int        nalloc_loc;      /* Allocation size for ind_loc                */
    rvec      *xa;              /* Positions for all IMD atoms assembled on
                                   master node                                */
    ivec      *xa_shifts;       /* Shifts for all IMD atoms, to make
                                   molecule(s) whole                          */
    ivec      *xa_eshifts;      /* Extra shifts since last DD step            */
    rvec      *xa_old;          /* Old positions for all IMD atoms, master    */
    int       *xa_ind;          /* Position of each local atom in the
                                 * collective array                           */

    int        nstimd;          /* global frequency, known to all nodes       */
    int        nstimd_new;      /* new frequency from IMD client, on master   */
    int        nstimd_def;      /* def. frequency to set to when disconnected */

    int        port;            /* Port to use for network socket             */
    IMDSocket *socket;          /* The IMDsocket on the master node           */
    IMDSocket *clientsocket;    /* ... and client                             */
    int        length;          /* length we got with last header             */

    gmx_bool   bWConnect;       /* shall we block and wait for connection?    */
    gmx_bool   bTerminated;     /* set if md is terminated                    */
    gmx_bool   bTerminatable;   /* set if md is terminatable                  */
    gmx_bool   bConnected;      /* set if connection present                  */
    gmx_bool   bNewForces;      /* set if we received new forces              */
    gmx_bool   bForceActivated; /* set if pulling from vmd is allowed         */

    gmx_bool   bCommSinceNS;    /* flag to see if we communicated coordinates */

    IMDEnergyBlock *energies;   /* pointer to energies we send back           */

    int32      vmd_nforces;     /* number of vmd forces                       */
    int32     *vmd_f_ind;       /* vmd indices                                */
    float     *vmd_forces;      /* the vmd forces flat in memory              */
    int        nforces;         /* same thing for the actual md forces,
                                 * this gets communicated to the clients.     */
    atom_id   *f_ind;
    rvec      *f;               /* The IMD pulling forces                     */

    int       npdlocalf;        /* Number of local particle decomp forces     */
    int      *pdlocalf;         /* buffer for local particle decomp forces    */

    char     *forcesendbuf;     /* buffer for force sending                   */
    char     *coordsendbuf;     /* buffer for coordinate sending              */
    char     *energysendbuf;    /* same for energies                          */
    rvec     *sendxbuf;         /* buf. to make mols. whole before sending    */

    t_block  mols;              /* molecules block in IMD group           */

    /* The next block is used on the master node only to reduce the output
     * without sacrificing information. If any of these values changes,
     * we need to write output */
    int       old_nforces;
    atom_id  *old_f_ind;
    rvec     *old_forces;

} t_gmx_IMD_setup;

/*Names of the IMDType for error messages*/
const char *eIMDType_names[IMD_NR + 1] =
{ "IMD_DISCONNECT", "IMD_ENERGIES", "IMD_FCOORDS", "IMD_GO", "IMD_HANDSHAKE", "IMD_KILL", "IMD_MDCOMM", "IMD_PAUSE", "IMD_TRATE", "IMD_IOERROR", NULL };


/*Re-implementation of basic IMD functions from NAMD/VMD*/

/*Fills the header with message and the length argument*/
static void fill_header(IMDHeader *header, IMDMessageType type, int32 length)
{
    /*we (ab-)use htonl network function for the correct endianness*/
    header->type = htonl((int32) type);
    header->length = htonl(length);
}

/*Swaps the endianess of the header*/
static void swap_header(IMDHeader *header)
{
    /*and vice versa...*/
    header->type = ntohl(header->type);
    header->length = ntohl(header->length);
}

/*Reads multiple bytes from socket*/
static int32 imd_read_multiple(IMDSocket *socket, char *datptr, int32 toread)
{
    int32 leftcount, countread;
    leftcount = toread;
    /*Try to read while we haven't reached toread*/
    while (leftcount != 0)
    {
        if ((countread = imdsock_read(socket, datptr, leftcount)) < 0)
        {
            /*interrupted function call, try again...*/
            if (errno == EINTR)
            {
                countread = 0;
            }
            /*this is an unexpected error, return what we got*/
            else
            {
                return toread - leftcount;
            }

            /*nothing read, finished*/
        }
        else if (countread == 0)
        {
            break;
        }
        leftcount -= countread;
        datptr += countread;
    } /* end while */
    /*return nr of bytes read.*/
    return toread - leftcount;
}

/*Writes multiple bytes to socket in analogy to imd_read_multiple*/
static int32 imd_write_multiple(IMDSocket *socket, const char *datptr, int32 towrite)
{
    int32 leftcount, countwritten;
    leftcount = towrite;
    while (leftcount != 0)
    {
        if ((countwritten = imdsock_write(socket, datptr, leftcount)) <= 0)
        {
            if (errno == EINTR)
            {
                countwritten = 0;
            }
            else
            {
                return towrite - leftcount;
            }
        }
        leftcount -= countwritten;
        datptr += countwritten;
    } /* end while */
    return towrite - leftcount;
}

/*handshake with IMD client*/
static int imd_handshake(IMDSocket *socket)
{
    IMDHeader header;
    fill_header(&header, IMD_HANDSHAKE, 1);
    header.length = IMDVERSION; /*client wants unswapped version*/
    return (imd_write_multiple(socket, (char *) &header, HEADERSIZE) != HEADERSIZE);
}

/*send energies using the energy block and the send buffer*/
static int imd_send_energies(IMDSocket *socket, const IMDEnergyBlock *energies, char *buffer)
{
    int32 recsize = HEADERSIZE + sizeof(IMDEnergyBlock);
    fill_header((IMDHeader *) buffer, IMD_ENERGIES, 1);
    memcpy(buffer + HEADERSIZE, energies, sizeof(IMDEnergyBlock));
    return (imd_write_multiple(socket, buffer, recsize) != recsize);
}

/*receive IMD header from socket, sets the length and returns the IMD message*/
static IMDMessageType imd_recv_header(IMDSocket *socket, int32 *length)
{
    IMDHeader header;
    if (imd_read_multiple(socket, (char *) &header, HEADERSIZE) != HEADERSIZE)
    {
        return IMD_IOERROR;
    }
    swap_header(&header);
    *length = header.length;
    return (IMDMessageType) header.type;
}

/*receive force indices and forces itself, the number of forces was previously communicated via the header.*/
static int imd_recv_mdcomm(IMDSocket *socket, int32 nforces, int32 *forcendx, float *forces)
{
    int retsize, retbytes;
    /*reading indices*/
    retsize = sizeof(int32) * nforces;
    retbytes = imd_read_multiple(socket, (char *) forcendx, retsize);
    if (retbytes != retsize)
    {
        return FALSE;
    }

    /*reading forces as float array*/
    retsize = 3 * sizeof(float) * nforces;
    retbytes = imd_read_multiple(socket, (char *) forces, retsize);
    if (retbytes != retsize)
    {
        return FALSE;
    }
    return TRUE;
}

/* GROMACS specific functions for the IMD implementation */

/* Called from grompp.
 *
 * Write out the group of atoms selected for interactive manipulation.
 */
extern void write_imdatoms(t_inputrec *ir, t_state *state, gmx_mtop_t *sys, const char *fn)
{
    t_atoms IMDatoms;


    IMDatoms = gmx_mtop_global_atoms(sys);
    write_sto_conf_indexed(fn,"IMDgroup",&IMDatoms,
    state->x,state->v,ir->ePBC,state->box,ir->imd->nat, ir->imd->ind);
}


/* Make a selection of the local atoms to be communicated via the IMD protocol.
 *
 * Each node checks which of the atoms from "ind" are local and puts its local
 * atom numbers into the "ind_local" array. Furthermore, in "xa_ind" it is
 * stored at which position each local atom belongs in the assembled/collective
 * array, so that on the master node all positions can be merged into the
 * assembled array correctly.
 */
extern void dd_make_local_IMD_atoms(gmx_domdec_t *dd, t_IMD *imd)
{
    gmx_ga2la_t        ga2la;
    t_gmx_IMD_setup    *IMDsetup;


    IMDsetup = imd->setup;
    ga2la = dd->ga2la;

    dd_make_local_group_indices(
            ga2la, IMDsetup->nat, IMDsetup->ind, &IMDsetup->nat_loc,
            &IMDsetup->ind_loc, &IMDsetup->nalloc_loc, IMDsetup->xa_ind);
}



/* Send positions from rvec, therefore we need a separate send buffer and conversion to Angstroem*/
static int imd_send_rvecs(IMDSocket *socket, int nat, rvec *x, char *buffer, atom_id *ind)
{
    int32 size;
    int i;
    float sendx[3];
    int tuplesize = 3 * sizeof(float);


    /* Required size for the send buffer */
    size = HEADERSIZE + 3 * sizeof(float) * nat;

    /* Prepare header */
    fill_header((IMDHeader *) buffer, IMD_FCOORDS, (int32) nat);
    for (i = 0; i < nat; i++)
    {
        sendx[0] = (float) x[i][0] * NM2A;
        sendx[1] = (float) x[i][1] * NM2A;
        sendx[2] = (float) x[i][2] * NM2A;
        memcpy(buffer + HEADERSIZE + i * tuplesize, sendx, tuplesize);
    }

    return (imd_write_multiple(socket, buffer, size) != size);
}


/* Initializes the IMD private data */
static t_gmx_IMD_setup* imd_create(int imdatoms)
{
    t_gmx_IMD_setup *IMDsetup = NULL;

    snew(IMDsetup, 1);
    IMDsetup->nat             = imdatoms;
    IMDsetup->bTerminated     = FALSE;
    IMDsetup->bTerminatable   = FALSE;
    IMDsetup->bWConnect       = FALSE;
    IMDsetup->bConnected      = FALSE;
    IMDsetup->bForceActivated = FALSE;
    IMDsetup->bNewForces      = FALSE;
    IMDsetup->bForceActivated = FALSE;
    IMDsetup->bCommSinceNS    = FALSE;
    IMDsetup->port            = IMDDEFPORT;
    IMDsetup->nstimd          = 1;
    IMDsetup->nstimd_new      = 1;

    return IMDsetup;
}

/*prepares the socket on the MASTER process/thread*/
static void imd_prepare_master_socket(t_gmx_IMD_setup *IMDsetup)
{
    int ret;

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    /*Winsock requires separate initialization*/
    fprintf(stderr,"%s Initializing winsock.\n",IMDstr);
    if (!imdsock_winsockinit())
    {
        gmx_fatal(FARGS,"%s Failed to initialize winsock.\n",IMDstr);
    }
#endif

    /*The rest is identical, first create and bind a socket and set to listen then.*/
    fprintf(stderr, "%s Setting up incoming socket.\n", IMDstr);
    IMDsetup->socket = imdsock_create();
    if (!IMDsetup->socket)
    {
        gmx_fatal(FARGS, "%s Failed to create socket.", IMDstr);
    }

    ret = imdsock_bind(IMDsetup->socket, IMDsetup->port);
    if (ret != 0)
    {
        gmx_fatal(FARGS, "%s binding socket failed with error %d.\n", IMDstr, ret);
    }

    fprintf(stderr, "%s Listening for IMD connection on port %d.\n", IMDstr, IMDsetup->port);
    ret = imd_sock_listen(IMDsetup->socket);
    if (ret != 0)
    {
        gmx_fatal(FARGS, "%s socket listen failed with error %d.\n", IMDstr, ret);
    }
}


/* Disconnect client */
static void imd_disconnect(t_gmx_IMD_setup *IMDsetup)
{
    /* Write out any buffered pulling data */
    fflush(IMDsetup->outf);

    /*we first try to shut down the clientsocket*/
    imdsock_shutdown(IMDsetup->clientsocket);
    if (!imdsock_destroy(IMDsetup->clientsocket))
    {
        fprintf(stderr, "%s Failed to destroy socket.\n", IMDstr);
    }

    /*then we reset the IMD step to its default, and reset the connection boolean*/
    IMDsetup->nstimd_new = IMDsetup->nstimd_def;
    IMDsetup->clientsocket = NULL;
    IMDsetup->bConnected = FALSE;
}

/*this is the imd_fatal which does not terminate gromacs (e.g. mdrun) but prints message and disconnects the client.*/
static void imd_fatal(t_gmx_IMD_setup *IMDsetup, const char *msg)
{
    fprintf(stderr, "%s %s", IMDstr, msg);
    imd_disconnect(IMDsetup);
    fprintf(stderr, "%s disconnected.\n", IMDstr);
}

/*check if we got an incoming connection*/
static gmx_bool imd_tryconnect(t_gmx_IMD_setup *IMDsetup)
{
    if (imdsock_tryread(IMDsetup->socket, 0, 0) > 0)
    {
        /*yea we got sth, accept on clientsocket*/
        IMDsetup->clientsocket = imdsock_accept(IMDsetup->socket);
        if (!IMDsetup->clientsocket)
        {
            fprintf(stderr, "%s Accepting the connection on the socket failed.\n", IMDstr);
            return FALSE;
        }

        /*handshake with client.*/
        if (imd_handshake(IMDsetup->clientsocket))
        {
            imd_fatal(IMDsetup, "Connection failed.\n");
            return FALSE;
        }

        fprintf(stderr, "%s Connection established, checking if I got IMD_GO orders.\n", IMDstr);

        /*Check if we get the proper "GO" command from client.*/
        if (imdsock_tryread(IMDsetup->clientsocket, IMDCONNECTWAIT, 0) != 1 || imd_recv_header(IMDsetup->clientsocket, &(IMDsetup->length)) != IMD_GO)
        {
            imd_fatal(IMDsetup, "No IMD_GO order received. IMD connection failed.\n");
        }

        /*IMD connected*/
        IMDsetup->bConnected = TRUE;
        return TRUE;
    }
    return FALSE;
}

/*this wraps tryconnect in order to make it "blocking". Used when the simulation should wait for an incoming connection (GMX_IMDWAIT).*/
static void imd_blockconnect(t_gmx_IMD_setup *IMDsetup)
{
    /*do not wait for connection, when e.g. ctrl+c is pressed and we will terminate anyways.*/
    if (!(int) gmx_get_stop_condition() == gmx_stop_cond_none)
    {
        return;
    }

    fprintf(stderr, "%s Will wait until I have a connection and IMD_GO orders.\n", IMDstr);

    /*while we have no clientsocket... 2nd part:we should still react on ctrl+c*/
    while ((!IMDsetup->clientsocket) && ((int) gmx_get_stop_condition() == gmx_stop_cond_none))
    {
        imd_tryconnect(IMDsetup);
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
        /*for whatever reason, it is called Sleep on windows*/
        Sleep(IMDLOOPWAIT);
#else
        sleep(IMDLOOPWAIT);
#endif
    }
}


/*make sure that our array holding the forces received via IMD is large enough.*/
static void imd_prepare_vmd_Forces(t_gmx_IMD_setup *IMDsetup)
{
    srenew((IMDsetup->vmd_f_ind), IMDsetup->vmd_nforces);
    srenew((IMDsetup->vmd_forces), 3*IMDsetup->vmd_nforces);
}

/*reads forces received via IMD*/
static void imd_read_vmd_Forces(t_gmx_IMD_setup *IMDsetup)
{
    /*the length of the previously received header tells us the nr of forces we will receive.*/
    IMDsetup->vmd_nforces = IMDsetup->length;
    /*prepare the arrays*/
    imd_prepare_vmd_Forces(IMDsetup);
    /*Now we read the forces...*/
    if (!(imd_recv_mdcomm(IMDsetup->clientsocket, IMDsetup->vmd_nforces, IMDsetup->vmd_f_ind, IMDsetup->vmd_forces)))
    {
        imd_fatal(IMDsetup, "Error while reading forces from remote. Disconnecting\n");
    }
}

/*prepares the md force arrays*/
static void imd_prepare_MD_Forces(t_gmx_IMD_setup *IMDsetup)
{
    srenew((IMDsetup->f_ind), IMDsetup->nforces);
    srenew((IMDsetup->f    ), IMDsetup->nforces);
}

/*prepares the index array of local forces for particle decomposition.*/
static void imd_prepare_PartDecomp(t_gmx_IMD_setup *IMDsetup)
{
    srenew((IMDsetup->pdlocalf), IMDsetup->nforces);
}

/* Copy IMD forces to MD forces, i.e. do conversion from Cal->Joule and from Angstrom -> nm and from a pointer array to arrays to 3*N array */
static void imd_copyto_MD_Forces(t_gmx_IMD_setup *IMDsetup)
{
    int i;
    real conversion = CAL2JOULE * NM2A;


    for (i = 0; i < IMDsetup->nforces; i++)
    {
        /* Copy the indices, a copy is important because we may update the incoming forces
         * whenever we receive new forces while the MD forces are only communicated upon
         * IMD communication.*/
        IMDsetup->f_ind[i] = IMDsetup->vmd_f_ind[i];

        /* Convert to rvecs and do a proper unit conversion */
        IMDsetup->f[i][0] = IMDsetup->vmd_forces[3*i    ] * conversion;
        IMDsetup->f[i][1] = IMDsetup->vmd_forces[3*i + 1] * conversion;
        IMDsetup->f[i][2] = IMDsetup->vmd_forces[3*i + 2] * conversion;
    }
}


/* Return TRUE if any of the forces or indices changed */
static gmx_bool bForcesChanged(t_gmx_IMD_setup *IMDsetup)
{
    int i;


    /* First, check whether the number of pulled atoms changed */
    if (IMDsetup->nforces != IMDsetup->old_nforces)
        return TRUE;

    /* Second, check whether any of the involved atoms changed */
    for (i=0; i<IMDsetup->nforces; i++)
    {
        if (IMDsetup->f_ind[i] != IMDsetup->old_f_ind[i])
            return TRUE;
    }

    /* Third, check whether all forces are the same */
    for (i=0; i<IMDsetup->nforces; i++)
    {
        if (IMDsetup->f[i][XX] != IMDsetup->old_forces[i][XX]) return TRUE;
        if (IMDsetup->f[i][YY] != IMDsetup->old_forces[i][YY]) return TRUE;
        if (IMDsetup->f[i][ZZ] != IMDsetup->old_forces[i][ZZ]) return TRUE;
    }

    /* All old and new forces are identical! */
    return FALSE;
}


/* Fill the
 *    old_f_ind and
 *    old_forces
 * arrays with the new, old values */
static void keep_old_values(t_gmx_IMD_setup *IMDsetup, double time)
{
    int i;


    IMDsetup->old_nforces = IMDsetup->nforces;

    for (i = 0; i < IMDsetup->nforces; i++)
    {
        IMDsetup->old_f_ind[i] = IMDsetup->f_ind[i];
        copy_rvec(IMDsetup->f[i], IMDsetup->old_forces[i]);
    }
}


/* Returns true if any component of the two rvecs differs */
static inline gmx_bool rvecs_differ(const rvec v1, const rvec v2)
{
    int i;


    for (i=0; i<DIM; i++)
    {
        if (v1[i] != v2[i])
        {
            return TRUE;
        }
    }

    return FALSE;
}

/* Write the applied pull forces to logfile. Call on master only! */
static void output_imd_forces(t_inputrec *ir, double time)
{
    t_gmx_IMD_setup *IMDsetup;
    int i;


    IMDsetup = ir->imd->setup;

    if (bForcesChanged(IMDsetup))
    {
        /* Write time and total number of applied IMD forces */
        fprintf(IMDsetup->outf, "%14.6e%6d", time, IMDsetup->nforces);

        /* Write out the global atom indices of the pulled atoms and the forces itself,
         * write out a force only if it has changed since the last output */
        for (i = 0; i < IMDsetup->nforces; i++)
        {
            if (rvecs_differ(IMDsetup->f[i], IMDsetup->old_forces[i]))
            {
                fprintf(IMDsetup->outf, "%9d", IMDsetup->ind[IMDsetup->f_ind[i]] + 1);
                fprintf(IMDsetup->outf, "%12.4e%12.4e%12.4e", IMDsetup->f[i][0], IMDsetup->f[i][1], IMDsetup->f[i][2]);
            }
        }
        fprintf(IMDsetup->outf, "\n");

        keep_old_values(IMDsetup, time);
    }
}

/*sync the different nodes*/
static void imd_sync_nodes(t_inputrec *ir, t_commrec *cr, double t)
{
    int new_nforces = 0;
    t_gmx_IMD_setup *IMDsetup;
    int start, end, i;


    /*just a local pointer...*/
    IMDsetup = ir->imd->setup;

    /*notify the non master nodes if we are still connected.*/
    if (PAR(cr))
    {
        block_bc(cr, IMDsetup->bConnected);
    }

    /*...if not connected, the job is done here.*/
    if (!IMDsetup->bConnected)
    {
        return;
    }

    /*let the non master nodes know if we got a new imd synchronization frequency.*/
    if (PAR(cr))
    {
        block_bc(cr, IMDsetup->nstimd_new);
    }

    /*Now we all set the (new) nstimd communication time step*/
    IMDsetup->nstimd = IMDsetup->nstimd_new;

    /*We're done if we don't allow pulling at all*/
    if (!(IMDsetup->bForceActivated))
    {
        return;
    }

    /*OK, lets check if we have received forces which we need to communicate to the non master nodes*/
    if (MASTER(cr))
    {
        if (IMDsetup->bNewForces)
        {
            new_nforces = IMDsetup->vmd_nforces;
        }
        /*make the "new_forces" negative, when we did not receive new ones.*/
        else
        {
            new_nforces = IMDsetup->vmd_nforces * -1;
        }
    }

    /*make new_forces known to the clients*/
    if (PAR(cr))
    {
        block_bc(cr, new_nforces);
    }

    /*When new_natoms < 0 then we know that these are still the same forces so we don't communicate them, otherwise...*/
    if (new_nforces >= 0)
    {
        /*set local vmd and nforces*/
        IMDsetup->vmd_nforces = new_nforces;
        IMDsetup->nforces = IMDsetup->vmd_nforces;

        /*now everybody knows the number of forces in f_ind, so we can prepare the target arrays for indices and forces*/
        imd_prepare_MD_Forces(IMDsetup);

        /*we first update the md forces on the master by converting the VMD forces*/
        if (MASTER(cr))
        {
            imd_copyto_MD_Forces(IMDsetup);
            /* We also write out forces on every update, so that we know which forces are applied for every step */
            if (IMDsetup->outf)
            {
                output_imd_forces(ir, t);
            }
        }

        /* In parallel mode we communicate the to-be-applied forces to the other nodes */
        if (PAR(cr))
        {
            nblock_bc(cr, IMDsetup->nforces, IMDsetup->f_ind);
            nblock_bc(cr, IMDsetup->nforces, IMDsetup->f    );
        }

        /*done communicating the forces, reset bNewForces*/
        IMDsetup->bNewForces = FALSE;

        /*for particle decomposition, we fill an array which holds the indices of the local forces*/
        if (PARTDECOMP(cr))
        {
            /*first make sure that the array is large enough*/
            imd_prepare_PartDecomp(IMDsetup);
            /*find out start and end atom of this pd node*/
            pd_at_range(cr, &start, &end);
            /*reset number of local forces*/
            IMDsetup->npdlocalf = 0;
            /*go through all forces...*/
            for (i = 0; i < IMDsetup->nforces; i++)
            {
                /*check if they are local on this pd node*/
                if (IMDsetup->f_ind[i] >= start && IMDsetup->f_ind[i] < end)
                {
                    /*if so add them to the local pd force indices and increment the counter*/
                    IMDsetup->pdlocalf[IMDsetup->npdlocalf] = i;
                    IMDsetup->npdlocalf++;
                }
            }
        }
        /*particle decomposition end*/
    }
}

/*reads header from the client and decides what to do*/
static void imd_readcommand(t_gmx_IMD_setup *IMDsetup)
{
    gmx_bool IMDpaused = FALSE;
    IMDMessageType itype;

    while (IMDsetup->clientsocket && (imdsock_tryread(IMDsetup->clientsocket, 0, 0) > 0 || IMDpaused))
    {
        itype = imd_recv_header(IMDsetup->clientsocket, &(IMDsetup->length));
        /*lets see what we got:*/
        switch (itype)
        {
            /*IMD asks us to terminate the simulation, check if the user allowed this via the GMX_IMDTERM variable.*/
            case IMD_KILL:
                if (IMDsetup->bTerminatable)
                {
                    fprintf(stderr, " %s Terminating connection and running simulation (if supported by integrator).\n", IMDstr);
                    IMDsetup->bTerminated = TRUE;
                    IMDsetup->bWConnect = FALSE;
                    gmx_set_stop_condition(gmx_stop_cond_next);
                }
                else
                {
                    fprintf(stderr, " %s Set -imdterm command line switch to allow mdrun termination from within IMD.\n", IMDstr);
                }

                break;

                /*the client doen't want to talk to us anymore*/
            case IMD_DISCONNECT:
                fprintf(stderr, " %s Disconnecting client.\n", IMDstr);
                imd_disconnect(IMDsetup);
                break;

                /*we got new forces, read them and set bNewForces flag*/
            case IMD_MDCOMM:
                imd_read_vmd_Forces(IMDsetup);
                IMDsetup->bNewForces = TRUE;
                break;

                /*the client asks us to (un)pause the simulation. So we toggle the IMDpaused state.*/
            case IMD_PAUSE:
                if (IMDpaused)
                {
                    fprintf(stderr, " %s Un-pause command received.\n", IMDstr);
                    IMDpaused = FALSE;
                }
                else
                {
                    fprintf(stderr, " %s Pause command received.\n", IMDstr);
                    IMDpaused = TRUE;
                }

                break;

                /*the client sets a new transfer rate, if we get 0, we reset the rate to the default. VMD filters 0 however.*/
            case IMD_TRATE:
                IMDsetup->nstimd_new = (IMDsetup->length > 0) ? IMDsetup->length : IMDsetup->nstimd_def;
                fprintf(stderr, " %s Update frequency will be set to %d.\n", IMDstr, IMDsetup->nstimd_new);
                break;

                /*Catch all rule for the remaining IMD types which we don't expect*/
            default:
                fprintf(stderr, " %s Received unexpected %s.\n", IMDstr, ENUM_NAME(itype,IMD_NR,eIMDType_names));
                imd_fatal(IMDsetup, "Terminating connection\n");
                break;
        } /*end switch */
    } /* end while */
}


/* Open IMD output file and write header information. Call on master only. */
static FILE *open_imd_out(
        int                  nfile,
        const t_filenm       fnm[],
        t_gmx_IMD_setup      *IMDsetup,
        int                  nat_total,      /* Number of atoms in the system */
        const output_env_t   oenv,
        unsigned long        Flags)
{
    FILE *fp;
    const char *fn;


    /* Check whether a log file of applied IMD forces was requested */
    if (opt2bSet("-if", nfile, fnm))
    {
        fn = opt2fn("-if", nfile, fnm);
        /* If we append to an existing file, all the header information is already there */
        if (Flags & MD_APPENDFILES)
        {
            fp = gmx_fio_fopen(fn, "a+");
        }
        else
        {
            fp = gmx_fio_fopen(fn, "w+");
            if (IMDsetup->nat == nat_total)
                fprintf(fp, "# Note that you can select an IMD index group in the .mdp file if a subset of the atoms suffices.\n");

            xvgr_header(fp, "IMD Pull Forces", "Time (ps)", "# of Forces / Atom IDs / Forces (kJ/mol)", exvggtNONE, oenv);

            fprintf(fp, "# Can display and manipulate %d (of a total of %d) atoms via IMD.\n", IMDsetup->nat, nat_total);
            fprintf(fp, "# column 1    : time (ps)\n");
            fprintf(fp, "# column 2    : total number of atoms feeling an IMD pulling force at that time\n");
            fprintf(fp, "# cols. 3.-6  : global atom number of pulled atom, x-force, y-force, z-force (kJ/mol)\n");
            fprintf(fp, "# then follow : atom-ID, f[x], f[y], f[z] for more atoms in case the force on multiple atoms is changed simultaneously.\n");
            fprintf(fp, "# Note that the force on any atom is always equal to the last value for that atom-ID found in the data.\n");
            fflush(fp);
        }

        /* To reduce the output file size we remember the old values and output only
         * when something changed */
        snew(IMDsetup->old_f_ind , IMDsetup->nat); /* One can never pull on more atoms */
        snew(IMDsetup->old_forces, IMDsetup->nat);

        return fp;
    } /*end logfile check*/

    fprintf(stdout, "%s For a log of the IMD pull forces explicitly specify '-if' on the command line.\n", IMDstr);

    return NULL;
}


/*imd finalize currently only closes the force output file*/
extern void imd_finalize(t_inputrec *ir)
{
    if (ir->bIMD)
    {
        if (ir->imd->setup->outf)
        {
            gmx_fio_fclose(ir->imd->setup->outf);
        }
    }
}

/* creates the molecule start-end position array only of molecules in the IMD group*/
static void init_imd_prepare_mols_in_imdgroup(t_gmx_IMD_setup *IMDsetup,gmx_mtop_t *top_global){
    int i,ii;
    int gstart,gend,count;
    t_block gmols,lmols;
    int nat;
    atom_id *ind;

    gmols = top_global->mols;
    nat = IMDsetup->nat;
    ind = IMDsetup->ind;

    lmols.nr=0;

    //check if index sorted
    for (i=0;i<nat-1;i++){
        if (ind[i]>ind[i+1]){
            gmx_fatal(FARGS, "%s IMD index is not sorted. This is currently not supported.\n", IMDstr);
        }
    }

    snew(lmols.index,gmols.nr+1);
    lmols.index[0]=0;

    for (i=0;i<gmols.nr;i++){
        gstart = gmols.index[i];
        gend   = gmols.index[i+1];
        count=0;
        for (ii=0;ii<nat;ii++)
        {
            if ((ind[ii]>=gstart)&&(ind[ii]<gend)){
                count+=1;
            }
        }
        if (count>0){
            lmols.index[lmols.nr+1]=lmols.index[lmols.nr]+count;
            lmols.nr+=1;
        }
    }
    srenew(lmols.index,lmols.nr+1);
    lmols.nalloc_index=lmols.nr+1;
    IMDsetup->mols = lmols;
}

/* copied and modified from groupcoord.c shift_positions_group*/

static void shift_positions(
        matrix box,
        rvec   x[],      /* The positions [0..nr] */
        ivec   is,       /* The shift [0..nr] */
        int    nr)       /* The number of positions */
{
    int      i,tx,ty,tz;

    /* Loop over the group's atoms */
    if(TRICLINIC(box))
    {
        for (i=0; i < nr; i++)
        {
            tx=is[XX];
            ty=is[YY];
            tz=is[ZZ];

            x[i][XX]=x[i][XX]-tx*box[XX][XX]-ty*box[YY][XX]-tz*box[ZZ][XX];
            x[i][YY]=x[i][YY]-ty*box[YY][YY]-tz*box[ZZ][YY];
            x[i][ZZ]=x[i][ZZ]-tz*box[ZZ][ZZ];
        }
    }
    else
    {
        for (i=0; i < nr; i++)
        {
            tx=is[XX];
            ty=is[YY];
            tz=is[ZZ];

            x[i][XX]=x[i][XX]-tx*box[XX][XX];
            x[i][YY]=x[i][YY]-ty*box[YY][YY];
            x[i][ZZ]=x[i][ZZ]-tz*box[ZZ][ZZ];
        }
    }
}


/* removes shifts of molecules diffused outside of the box*/
static void imd_remove_molshifts(t_gmx_IMD_setup *IMDsetup,matrix box){
    int i,ii,molsize;
    ivec largest,smallest,shift;
    t_block mols;
    mols = IMDsetup->mols;

    /*for each molecule also present in IMD group*/
    for (i=0;i<mols.nr;i++){
        /* first we determine the minimum and maximum shifts for each molecule*/

        clear_ivec(largest);
        clear_ivec(smallest);
        clear_ivec(shift);

        copy_ivec(IMDsetup->xa_shifts[mols.index[i]],largest);
        copy_ivec(IMDsetup->xa_shifts[mols.index[i]],smallest);

        for (ii=mols.index[i]+1;ii<mols.index[i+1];ii++){
            if (IMDsetup->xa_shifts[ii][XX]>largest[XX])  largest[XX]  = IMDsetup->xa_shifts[ii][XX];
            if (IMDsetup->xa_shifts[ii][XX]<smallest[XX]) smallest[XX] = IMDsetup->xa_shifts[ii][XX];

            if (IMDsetup->xa_shifts[ii][YY]>largest[YY])  largest[YY]  = IMDsetup->xa_shifts[ii][YY];
            if (IMDsetup->xa_shifts[ii][YY]<smallest[YY]) smallest[YY] = IMDsetup->xa_shifts[ii][YY];

            if (IMDsetup->xa_shifts[ii][ZZ]>largest[ZZ])  largest[ZZ]  = IMDsetup->xa_shifts[ii][ZZ];
            if (IMDsetup->xa_shifts[ii][ZZ]<smallest[ZZ]) smallest[ZZ] = IMDsetup->xa_shifts[ii][ZZ];

        }

        /*check if we what we can subtract/add to the coordinates to put them back to the box*/
        if (smallest[XX] > 0) shift[XX] = smallest[XX];
        if (smallest[YY] > 0) shift[YY] = smallest[YY];
        if (smallest[ZZ] > 0) shift[ZZ] = smallest[ZZ];

        if (largest[XX] < 0) shift[XX] = largest[XX];
        if (largest[YY] < 0) shift[YY] = largest[YY];
        if (largest[ZZ] < 0) shift[ZZ] = largest[ZZ];

        /*is there a shift at all?*/
        if ((shift[XX])||(shift[YY])||(shift[ZZ]))
        {
            molsize=mols.index[i+1]-mols.index[i];
            /*shift the coordinates*/
            shift_positions(box,&(IMDsetup->xa[mols.index[i]]),shift,molsize);
        }

    }
}


/* Initialize arrays used to assemble the positions from the other nodes */
static void init_imd_prepare_for_x_assembly(t_commrec *cr, rvec x[], t_mdatoms *md, t_gmx_IMD_setup *IMDsetup)
{
    int i,ii;

    snew(IMDsetup->xa        , IMDsetup->nat);
    snew(IMDsetup->xa_ind    , IMDsetup->nat);
    snew(IMDsetup->xa_shifts , IMDsetup->nat);
    snew(IMDsetup->xa_eshifts, IMDsetup->nat);
    snew(IMDsetup->xa_old    , IMDsetup->nat);

    /* Save the original (whole) set of positions such that later the
     * molecule can always be made whole again */
    if (MASTER(cr))
    {
        for (i=0; i<IMDsetup->nat; i++)
        {
            ii = IMDsetup->ind[i];
            copy_rvec(x[ii], IMDsetup->xa_old[i]);
        }
    }

    if (!PAR(cr))
    {
        IMDsetup->nat_loc = IMDsetup->nat;
        IMDsetup->ind_loc = IMDsetup->ind;
    }

    if (!DOMAINDECOMP(cr))
    {
        /* xa_ind needs to be set to identity for
         * serial runs or particle decomposition */
        for (i=0; i<IMDsetup->nat; i++)
            IMDsetup->xa_ind[i] = i;
    }

    /* Store local atom indices and # of local atoms for PD*/
    if (PARTDECOMP(cr))
    {
        snew(IMDsetup->ind_loc, IMDsetup->nat);
        IMDsetup->nat_loc=0;
        for (i=0; i<IMDsetup->nat; i++)
        {
            ii = IMDsetup->ind[i];
            if (ii >= md->start && ii < md->start+md->homenr)
            {
                IMDsetup->ind_loc[IMDsetup->nat_loc++] = ii;
            }
        }
    }

    /*Communicate initial coordinates xa_old to all processes*/
#ifdef GMX_MPI
    if (PAR(cr))
        gmx_bcast(IMDsetup->nat * sizeof(IMDsetup->xa_old[0]), IMDsetup->xa_old, cr);
#endif
}

/* Check for non-working integrator / parallel options */
static void imd_check_integrator_parallel(t_inputrec *ir,t_commrec *cr)
{
    if (PARTDECOMP(cr))
    {
        gmx_fatal(FARGS, "%s Particle decomposition is currently not supported by IMD.\n", IMDstr);
        return;
    }

    if (DOMAINDECOMP(cr))
    {
        if (((ir->eI) == eiSteep)||((ir->eI) == eiCG)||((ir->eI) == eiLBFGS)||((ir->eI) == eiNM))
        {
            gmx_fatal(FARGS, "%s Energy, minimization via steep, CG, lbfgs and nm in parallel is currently not supported by IMD.\n", IMDstr);
            return;
        }
    }
}

/* This function initializes or disables IMD, it is called before the main loop,
 * and it must be called prior to any call to dd_partition_system if in parallel */
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
        unsigned long Flags)
{
    int i;
    t_gmx_IMD_setup *IMDsetup;
    gmx_bool bIMDallowed, bIMDenabled;
    int32 bufxsize;
    gmx_bool bIMD = FALSE;


    if (MASTER(cr))
    {
        /* The default is disabled IMD, we first check if tpr allows IMD, if so we check if it's actually
         * enabled by an environment variable.*/
        bIMDenabled = FALSE;

        /*Check the command line switches if imd was enabled (by setting at least one.).*/
        if ((imdport) || (imdfreq) || (Flags & MD_IMDWAIT) || (Flags & MD_IMDTERM) || (Flags & MD_IMDPULL))
        {
            bIMDenabled=TRUE;
        }

        /*Write out reason why IMD is not activated.*/
        /*IMD not allowed?*/
        if (!ir->bIMD)
        {
            fprintf(stderr, "%s IMD disabled in TPR.\n"
                    "%s This run will not accept IMD connections\n", IMDstr, IMDstr);
        }
        /*IMD is allowed by TPR*/
        else
        {
            /*IMD not enabled by one of the -imd switches?*/
            if (!bIMDenabled)
            {
                fprintf(stderr, "%s None of the -imd switches was used.\n"
                        "%s This run will not accept incoming IMD connections\n", IMDstr,IMDstr);
            }
            /*one of the IMD switches was used...*/
            else
            {
                /*Multiple simulations or replica exchange -> Doh!*/
                if (MULTISIM(cr))
                {
                    fprintf(stderr, "%s IMD for multiple simulations or replica exchange is currently not supported.\n", IMDstr);
                    bIMD = FALSE;
                }
                /*OK, IMD seems to be allowed and turned on...*/
                else
                {
                    fprintf(stderr,"%s Enabled. This simulation will accept incoming IMD connections.\n",IMDstr);
                    bIMD = TRUE;
                }

            }
        }
    } /*end master only*/

    /* Let the other nodes know if we want IMD */
    if (PAR(cr))
    {
        block_bc(cr, bIMD);
    }

    /* set our local inputrec bIMD to check if we want IMD during the simulation */
    ir->bIMD = bIMD;

    /*... if not we are done.*/
    if (!ir->bIMD)
    {
        return;
    }

    /*check if we're using a sane integrator / parallel combination*/
    imd_check_integrator_parallel(ir,cr);


    /*******************************************************/
    /** From here on we assume that IMD is turned on      **/
    /*******************************************************/

    /* Initialize IMD setup structure. If we read in a pre-IMD .tpr file, imd->nat
     * will be zero. For those cases we transfer _all_ atomic positions */
    ir->imd->setup = imd_create(ir->imd->nat > 0 ? ir->imd->nat : nat_total);
    IMDsetup = ir->imd->setup;

    /* We might need to open an output file for IMD forces data */
    if (MASTER(cr))
    {
        IMDsetup->outf = open_imd_out(nfile, fnm, ir->imd->setup, nat_total, oenv, Flags);
    }


    /* Make sure that we operate with a valid atom index array for the IMD atoms */
    if (ir->imd->nat > 0)
    {
        /* Point to the user-supplied array of atom numbers */
        IMDsetup->ind = ir->imd->ind;
    }
    else
    {
        /* Make a dummy (ind[i] = i) array of all atoms */
        snew(IMDsetup->ind, nat_total);
        for (i=0; i<nat_total; i++)
        {
            IMDsetup->ind[i] = i;
        }
    }

    /*read environment on master and prepare socket for incoming connections*/
    if (MASTER(cr))
    {
        /*we allocate memory for our IMD energy structure*/
        int32 recsize = HEADERSIZE + sizeof(IMDEnergyBlock);
        snew(IMDsetup->energysendbuf, recsize);

        /*Now we check if we should set the listening port to sth else than the default*/
        if (imdport)
        {
            IMDsetup->port=imdport;
            fprintf(stderr, "%s Setting non-default port for incoming connections to %d (-imdport).\n", IMDstr, IMDsetup->port);
        }

        /*Shall we use a special IMD update/communication frequency?*/
        if (imdfreq)
        {
            IMDsetup->nstimd=imdfreq;
            IMDsetup->nstimd_new=imdfreq;
            fprintf(stderr, "%s Setting IMD update frequency to %d steps. (-imdfreq).\n", IMDstr, IMDsetup->nstimd);
        }

        /*Shall we wait for a connection?*/
        if (Flags & MD_IMDWAIT)
        {
            IMDsetup->bWConnect = TRUE;
            fprintf(stderr, "%s Pausing simulation while no IMD connection present (-imdwait).\n", IMDstr);
        }

        /*Will the IMD clients be able to terminate the simulation?*/
        if (Flags & MD_IMDTERM)
        {
            IMDsetup->bTerminatable = TRUE;
            fprintf(stderr, "%s Allow termination of the simulation from IMD client (-imdterm).\n", IMDstr);
        }

        /* Is pulling from IMD client allowed? */
        if (Flags & MD_IMDPULL)
        {
            IMDsetup->bForceActivated = TRUE;
            fprintf(stderr, "%s Pulling from IMD remote is enabled (-imdpull).\n", IMDstr);
        }

        /* Initialize send buffers with constant size */
        snew(IMDsetup->sendxbuf, IMDsetup->nat);
        snew(IMDsetup->energies, 1);
        bufxsize = HEADERSIZE + 3 * sizeof(float) * IMDsetup->nat;
        snew(IMDsetup->coordsendbuf, bufxsize);

    }

    /*do we allow interactive pulling? If so let the other nodes know.*/
    if (PAR(cr))
    {
        block_bc(cr, IMDsetup->bForceActivated);
    }

    /*set the current, the new and the default update frequency.*/
    IMDsetup->nstimd     = defnstimd;
    IMDsetup->nstimd_def = defnstimd;
    IMDsetup->nstimd_new = defnstimd;

    /*everybody should know the IMD update frequency so that we can do communication together.*/
    if (PAR(cr))
    {
        block_bc(cr, IMDsetup->nstimd);
    }

    /* setup the listening socket on master process*/
    if (MASTER(cr))
    {
        fprintf(fplog, "%s Setting port for connection requests to %d.\n", IMDstr, IMDsetup->port);
        fprintf(stderr, "%s Turning on IMD - port for incoming requests is %d.\n", IMDstr, IMDsetup->port);
        imd_prepare_master_socket(IMDsetup);
        /*Wait until we have a connection if specified before*/
        if (IMDsetup->bWConnect)
        {
            imd_blockconnect(IMDsetup);
        }
        else
        {
            fprintf(stderr, "%s GMX_IMDWAIT is not set, starting simulation.\n", IMDstr);
        }
    }
    /* Let the other nodes know whether we are connected */
    imd_sync_nodes(ir, cr, 0);

    /* Initialize arrays used to assemble the positions from the other nodes */
    init_imd_prepare_for_x_assembly(cr, x, md,IMDsetup);

    /* Initialize molecule blocks to make them whole later...*/
    if (MASTER(cr))
    {
        init_imd_prepare_mols_in_imdgroup(IMDsetup,top_global);
    }
}


/* Returns if we should do IMD communication in this step.
 * Also checks for new IMD connection and syncs the nodes. */
extern gmx_bool do_IMD(
        int          step,
        t_commrec   *cr,
        gmx_bool     bNS,       /* Is this a neighborsearching step?          */
        matrix       box,
        rvec         x[],       /* The atomic positions local to this node    */
        t_inputrec  *ir,
        double       t)
{
    gmx_bool imdstep = FALSE;
    t_gmx_IMD_setup *IMDsetup;

    /*IMD at all?*/
    if (!ir->bIMD)
    {
        return FALSE;
    }

    /*local pointer*/
    IMDsetup = ir->imd->setup;

    /*read command from client and check if new incoming connection*/
    if (MASTER(cr))
    {
        /*If not already connected, check for new connections*/
        if (!IMDsetup->clientsocket)
        {
            if (IMDsetup->bWConnect)
            {
                imd_blockconnect(IMDsetup);
            }
            else
            {
                imd_tryconnect(IMDsetup);
            }
        }

        /*Lets see if we have new IMD messages for us*/
        if (IMDsetup->clientsocket)
        {
            imd_readcommand(IMDsetup);
        }
    }

    /*is this an IMD communication step?*/
    imdstep = do_per_step(step, IMDsetup->nstimd);

    /*OK so this is an IMD step...*/
    if (imdstep)
    {
        /*First we sync all nodes that the others also know if we are connected*/
        imd_sync_nodes(ir, cr, t);

        /*If a client is connected, we collect the coordinates and put molecules back in the box if necessary*/
        if (IMDsetup->bConnected)
        {

            /* Transfer the IMD positions to the master node. Every node contributes
               its local positions x and stores it in the assembled "xa" array. */
            communicate_group_positions(cr,
                    IMDsetup->xa, IMDsetup->xa_shifts, IMDsetup->xa_eshifts, bNS,
                    x, IMDsetup->nat, IMDsetup->nat_loc, IMDsetup->ind_loc, IMDsetup->xa_ind, IMDsetup->xa_old, box);

            /*Avoid extra-communication at next NS step.*/
            IMDsetup->bCommSinceNS = TRUE;

            /*If connected and master remove shifts*/
            if (MASTER(cr))
            {
                imd_remove_molshifts(IMDsetup,box);
            }
        }

    } /*end if imdstep*/

    else
    {
        if (bNS)
        {
            /*We communicated the coordinates already since last ns step*/
            if (IMDsetup->bCommSinceNS){
                IMDsetup->bCommSinceNS=FALSE;
            }
            else
            {
               /*independently on imdstep, we communicate positions at each NS step*/
               /* Transfer the IMD positions to the master node. Every node contributes
                * its local positions x and stores it in the assembled "xa" array. */
               communicate_group_positions(cr,
                   IMDsetup->xa, IMDsetup->xa_shifts, IMDsetup->xa_eshifts, bNS,
                   x, IMDsetup->nat, IMDsetup->nat_loc, IMDsetup->ind_loc, IMDsetup->xa_ind, IMDsetup->xa_old, box);
            }
        } /* end ns step */
    } /*end not imdstep*/

    return imdstep;
}



/* Copy energies and convert to float from enerdata to the IMD energy record.
 * We do no conversion so units in client are SI */
extern void do_imd_prepare_energies(t_gmx_IMD IMDsetup, gmx_enerdata_t *enerd, int step, gmx_bool bHaveNewEnergies)
{
    IMDEnergyBlock *ene;


    if (IMDsetup->clientsocket)
    {
        ene = IMDsetup->energies;

        ene->tstep   = step;

        /* In parallel simulations the energies are not accessible at every time step.
         * We update them only if we have new values. */
        if (bHaveNewEnergies)
        {
            ene->T_abs   = (float)  enerd->term[F_TEMP   ];
            ene->E_tot   = (float)  enerd->term[F_ETOT   ];
            ene->E_bond  = (float)  enerd->term[F_BONDS  ];
            ene->E_angle = (float)  enerd->term[F_ANGLES ];
            ene->E_dihe  = (float)  enerd->term[F_PDIHS  ];
            ene->E_impr  = (float)  enerd->term[F_IDIHS  ];
            ene->E_vdw   = (float)  enerd->term[F_LJ     ];
            ene->E_coul  = (float) (enerd->term[F_COUL_SR] + enerd->term[F_COUL_LR]);
        }
    }
}


/* Send positions and energy to the client */
extern void do_imd_send_positions(t_IMD *imd)
{
    t_gmx_IMD IMDsetup;


    IMDsetup = imd->setup;

    if (IMDsetup->clientsocket)
    {

        if (imd_send_energies(IMDsetup->clientsocket, IMDsetup->energies, IMDsetup->energysendbuf))
        {
            imd_fatal(IMDsetup, "Error sending updated energies. Disconnecting client.\n");
        }

        if (imd_send_rvecs(IMDsetup->clientsocket, IMDsetup->nat, IMDsetup->xa, IMDsetup->coordsendbuf, IMDsetup->ind))
        {
            imd_fatal(IMDsetup, "Error sending updated positions. Disconnecting client.\n");
        }
    }
}

/*access nstimd, the IMD update/communication frequency*/
extern int imd_get_step(t_gmx_IMD IMDsetup)
{
    return IMDsetup->nstimd;
}

/* Apply the interactive pulling forces */
extern void imd_apply_forces(t_inputrec *ir, t_commrec *cr, rvec *f)
{
    int i, j;
    int locndx;
    t_gmx_IMD_setup *IMDsetup;


    IMDsetup = ir->imd->setup;

    /* Are forces allowed at all? If not we're done */
    if (!IMDsetup->bForceActivated)
    {
        return;
    }

    if (PARTDECOMP(cr))
    {
        /* Parallel run with particle decomposition*/
        for (i = 0; i < IMDsetup->npdlocalf; i++)
        {
            /*j are the indices in the "System group"
             * IMDsetup->pdlocalf[i] is the index of the i-th local force in the global IMD force array*/
            j = IMDsetup->ind[IMDsetup->f_ind[IMDsetup->pdlocalf[i]]];
            rvec_inc(f[j], IMDsetup->f[IMDsetup->pdlocalf[i]]);
        }
    }
    else
    {
        /* Serial or domain decomposition */
        for (i = 0; i < IMDsetup->nforces; i++)
        {
            /* j are the indices in the "System group".*/
            j = IMDsetup->ind[IMDsetup->f_ind[i]];

            /* check if this is a local atom and find out locndx */
            if (PAR(cr) && ga2la_get_home(cr->dd->ga2la, j, &locndx))
            {
                j = locndx;
            }

            rvec_inc(f[j], IMDsetup->f[i]);
        }
    }
}

#endif
