/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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

/*! \internal \file
 *
 * \brief
 * Implements Interactive Molecular Dynamics
 *
 * Re-implementation of basic IMD functions to work with VMD,
 * see imdsocket.h for references to the IMD API.
 *
 * \author Martin Hoefling, Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \ingroup module_imd
 */
#include "gmxpre.h"

#include "imd.h"

#include "config.h"

#include <cerrno>
#include <cstdio>
#include <cstring>

#include <array>
#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/imd/imdsocket.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
class ForceProviders;
class IMDOutputProvider;
class IMdpOptionProvider;
struct IMDSocket;
struct MDModulesNotifiers;

/*! \brief How long shall we wait in seconds until we check for a connection again? */
constexpr int c_loopWait = 1;

/*! \brief How long shall we check for the IMD_GO? */
constexpr int c_connectWait = 1;

/*! \brief IMD Header Size. */
constexpr int c_headerSize = 8;

/*! \brief IMD Protocol Version. */
constexpr int c_protocolVersion = 2;

enum class IMDMessageType : int;

/*! \internal
 * \brief
 * IMD (interactive molecular dynamics) energy record.
 *
 * As in the original IMD implementation. Energies in kcal/mol.
 * NOTE: We return the energies in GROMACS / SI units,
 * so they also show up as SI in VMD.
 *
 */
typedef struct
{
    int32_t tstep;   /**< time step                                     */
    float   T_abs;   /**< absolute temperature                          */
    float   E_tot;   /**< total energy                                  */
    float   E_pot;   /**< potential energy                              */
    float   E_vdw;   /**< van der Waals energy                          */
    float   E_coul;  /**< Coulomb interaction energy                    */
    float   E_bond;  /**< bonds energy                                  */
    float   E_angle; /**< angles energy                                 */
    float   E_dihe;  /**< dihedrals energy                              */
    float   E_impr;  /**< improper dihedrals energy                     */
} IMDEnergyBlock;


/*! \internal
 * \brief IMD (interactive molecular dynamics) communication structure.
 *
 * This structure defines the IMD communication message header & protocol version.
 */
typedef struct
{
    int32_t type;   /**< Type of IMD message, see IMDType_t above      */
    int32_t length; /**< Length                                        */
} IMDHeader;


/*! \internal
 * \brief Implementation type for the IMD session
 *
 * \todo Make this class (or one extracted from it) model
 * IForceProvider.
 * \todo Use RAII for files and allocations
 */
class ImdSession::Impl
{
public:
    //! Constructor
    Impl(const MDLogger& mdlog);
    ~Impl();

    /*! \brief Prepare the socket on the MAIN. */
    void prepareMainSocket();
    /*! \brief Disconnect the client. */
    void disconnectClient();
    /*! \brief Prints an error message and disconnects the client.
     *
     *  Does not terminate mdrun!
     */
    void issueFatalError(const char* msg);
    /*! \brief Check whether we got an incoming connection. */
    bool tryConnect();
    /*! \brief Wrap imd_tryconnect in order to make it blocking.
     *
     * Used when the simulation should wait for an incoming connection.
     */
    void blockConnect();
    /*! \brief Make sure that our array holding the forces received via IMD is large enough. */
    void prepareVmdForces();
    /*! \brief Reads forces received via IMD. */
    void readVmdForces();
    /*! \brief Prepares the MD force arrays. */
    void prepareMDForces();
    /*! \brief Copy IMD forces to MD forces.
     *
     * Do conversion from Cal->Joule and from
     * Angstrom -> nm and from a pointer array to arrays to 3*N array.
     */
    void copyToMDForces() const;
    /*! \brief Return true if any of the forces or indices changed. */
    bool bForcesChanged() const;
    /*! \brief Update the old_f_ind and old_forces arrays to contain the current values. */
    void keepOldValues();
    /*! \brief Write the applied pull forces to logfile.
     *
     * Call on main only!
     */
    void outputForces(double time);
    /*! \brief Synchronize the nodes. */
    void syncNodes(const t_commrec* cr, double t);
    /*! \brief Reads header from the client and decides what to do. */
    void readCommand();
    /*! \brief Open IMD output file and write header information.
     *
     * Call on main only.
     */
    void openOutputFile(const char* fn, int nat_total, const gmx_output_env_t* oenv, StartingBehavior startingBehavior);
    /*! \brief Creates the molecule start-end position array of molecules in the IMD group. */
    void prepareMoleculesInImdGroup(const gmx_mtop_t& top_global);
    /*! \brief Removes shifts of molecules diffused outside of the box. */
    void removeMolecularShifts(const matrix box) const;
    /*! \brief Initialize arrays used to assemble the positions from the other nodes. */
    void prepareForPositionAssembly(const t_commrec* cr, gmx::ArrayRef<const gmx::RVec> coords);
    /*! \brief Interact with any connected VMD session */
    bool run(int64_t step, bool bNS, const matrix box, gmx::ArrayRef<const gmx::RVec> coords, double t);

    // TODO rename all the data members to have underscore suffixes

    //! True if tpr and mdrun input combine to permit IMD sessions
    bool sessionPossible = false;
    //! Output file for IMD data, mainly forces.
    FILE* outf = nullptr;

    //! Number of atoms that can be pulled via IMD.
    int nat = 0;
    //! Part of the atoms that are local.
    int nat_loc = 0;
    //! Global indices of the IMD atoms.
    int* ind = nullptr;
    //! Local indices of the IMD atoms.
    int* ind_loc = nullptr;
    //! Allocation size for ind_loc.
    int nalloc_loc = 0;
    //! Positions for all IMD atoms assembled on the main node.
    rvec* xa = nullptr;
    //! Shifts for all IMD atoms, to make molecule(s) whole.
    ivec* xa_shifts = nullptr;
    //! Extra shifts since last DD step.
    ivec* xa_eshifts = nullptr;
    //! Old positions for all IMD atoms on main.
    rvec* xa_old = nullptr;
    //! Position of each local atom in the collective array.
    int* xa_ind = nullptr;

    //! Global IMD frequency, known to all ranks.
    int nstimd = 1;
    //! New frequency from IMD client, main only.
    int nstimd_new = 1;
    //! Default IMD frequency when disconnected.
    int defaultNstImd = -1;

    //! Port to use for network socket.
    int port = 0;
    //! The IMD socket on the main node.
    IMDSocket* socket = nullptr;
    //! The IMD socket on the client.
    IMDSocket* clientsocket = nullptr;
    //! Length we got with last header.
    int length = 0;

    //! Shall we block and wait for connection?
    bool bWConnect = false;
    //! Set if MD can be terminated.
    bool bTerminatable = false;
    //! Set if connection is present.
    bool bConnected = false;
    //! Set if we received new forces.
    bool bNewForces = false;
    //! Set if pulling from VMD is allowed.
    bool bForceActivated = false;

    //! Pointer to energies we send back.
    IMDEnergyBlock* energies = nullptr;

    //! Number of VMD forces.
    int32_t vmd_nforces = 0;
    //! VMD forces indices.
    int32_t* vmd_f_ind = nullptr;
    //! The VMD forces flat in memory.
    float* vmd_forces = nullptr;
    //! Number of actual MD forces; this gets communicated to the clients.
    int nforces = 0;
    //! Force indices.
    int* f_ind = nullptr;
    //! The IMD pulling forces.
    rvec* f = nullptr;

    //! Buffer for coordinate sending.
    char* coordsendbuf = nullptr;
    //! Send buffer for energies.
    char* energysendbuf = nullptr;
    //! Buffer to make molecules whole before sending.
    rvec* sendxbuf = nullptr;

    //! Molecules block in IMD group.
    t_block mols;

    /* The next block is used on the main node only to reduce the output
     * without sacrificing information. If any of these values changes,
     * we need to write output */
    //! Old value for nforces.
    int old_nforces = 0;
    //! Old values for force indices.
    int* old_f_ind = nullptr;
    //! Old values for IMD pulling forces.
    rvec* old_forces = nullptr;

    //! Logger
    const MDLogger& mdLog_;
    //! Commmunication object
    const t_commrec* cr_ = nullptr;
    //! Wallcycle counting manager.
    gmx_wallcycle* wcycle = nullptr;
    //! Energy output handler
    gmx_enerdata_t* enerd = nullptr;
};

/*! \internal
 * \brief Implement interactive molecular dynamics.
 *
 * \todo Some aspects of this module provides forces (when the user
 * pulls on things in VMD), so in future it should have a class that
 * models IForceProvider and is contributed to the collection of such
 * things.
 */
class InteractiveMolecularDynamics final : public IMDModule
{
    // From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return nullptr; }
    IMDOutputProvider*  outputProvider() override { return nullptr; }
    void                initForceProviders(ForceProviders* /* forceProviders */) override {}
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}
};

std::unique_ptr<IMDModule> createInteractiveMolecularDynamicsModule()
{
    return std::make_unique<InteractiveMolecularDynamics>();
}

/*! \brief Enum for types of IMD messages.
 *
 * We use the same records as the NAMD/VMD IMD implementation.
 */
enum class IMDMessageType : int
{
    Disconnect, /**< client disconnect                               */
    Energies,   /**< energy data                                     */
    FCoords,    /**< atomic coordinates                              */
    Go,         /**< start command for the simulation                */
    Handshake,  /**< handshake to determine little/big endianness    */
    Kill,       /**< terminates the simulation                       */
    Mdcomm,     /**< force data                                      */
    Pause,      /**< pauses the simulation                           */
    TRate,      /**< sets the IMD transmission and processing rate   */
    IOerror,    /**< I/O error                                       */
    Count       /**< number of entries                               */
};


/*! \brief Names of the IMDType for error messages.
 */
static const char* enumValueToString(IMDMessageType enumValue)
{
    constexpr gmx::EnumerationArray<IMDMessageType, const char*> imdMessageTypeNames = {
        "IMD_DISCONNECT", "IMD_ENERGIES", "IMD_FCOORDS", "IMD_GO",    "IMD_HANDSHAKE",
        "IMD_KILL",       "IMD_MDCOMM",   "IMD_PAUSE",   "IMD_TRATE", "IMD_IOERROR"
    };
    return imdMessageTypeNames[enumValue];
}


/*! \brief Fills the header with message and the length argument. */
static void fill_header(IMDHeader* header, IMDMessageType type, int32_t length)
{
    /* We (ab-)use htonl network function for the correct endianness */
    header->type   = imd_htonl(static_cast<int>(type));
    header->length = imd_htonl(length);
}


/*! \brief Swaps the endianess of the header. */
static void swap_header(IMDHeader* header)
{
    /* and vice versa... */
    header->type   = imd_ntohl(static_cast<int>(header->type));
    header->length = imd_ntohl(header->length);
}


/*! \brief Reads multiple bytes from socket. */
static int32_t imd_read_multiple(IMDSocket* socket, char* datptr, int32_t toread)
{
    int32_t leftcount, countread;


    leftcount = toread;
    /* Try to read while we haven't reached toread */
    while (leftcount != 0)
    {
        if ((countread = imdsock_read(socket, datptr, leftcount)) < 0)
        {
            /* interrupted function call, try again... */
            if (errno == EINTR)
            {
                countread = 0;
            }
            /* this is an unexpected error, return what we got */
            else
            {
                return toread - leftcount;
            }

            /* nothing read, finished */
        }
        else if (countread == 0)
        {
            break;
        }
        leftcount -= countread;
        datptr += countread;
    } /* end while */

    /* return nr of bytes read */
    return toread - leftcount;
}


/*! \brief Writes multiple bytes to socket in analogy to imd_read_multiple. */
static int32_t imd_write_multiple(IMDSocket* socket, const char* datptr, int32_t towrite)
{
    int32_t leftcount, countwritten;


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


/*! \brief Handshake with IMD client. */
static int imd_handshake(IMDSocket* socket)
{
    IMDHeader header;


    fill_header(&header, IMDMessageType::Handshake, 1);
    header.length = c_protocolVersion; /* client wants unswapped version */

    return static_cast<int>(imd_write_multiple(socket, reinterpret_cast<char*>(&header), c_headerSize)
                            != c_headerSize);
}


/*! \brief Send energies using the energy block and the send buffer. */
static int imd_send_energies(IMDSocket* socket, const IMDEnergyBlock* energies, char* buffer)
{
    int32_t recsize;


    recsize = c_headerSize + sizeof(IMDEnergyBlock);
    fill_header(reinterpret_cast<IMDHeader*>(buffer), IMDMessageType::Energies, 1);
    memcpy(buffer + c_headerSize, energies, sizeof(IMDEnergyBlock));

    return static_cast<int>(imd_write_multiple(socket, buffer, recsize) != recsize);
}


/*! \brief Receive IMD header from socket, sets the length and returns the IMD message. */
static IMDMessageType imd_recv_header(IMDSocket* socket, int32_t* length)
{
    IMDHeader header;


    if (imd_read_multiple(socket, reinterpret_cast<char*>(&header), c_headerSize) != c_headerSize)
    {
        return IMDMessageType::IOerror;
    }
    swap_header(&header);
    *length = header.length;

    return static_cast<IMDMessageType>(header.type);
}


/*! \brief Receive force indices and forces.
 *
 * The number of forces was previously communicated via the header.
 */
static bool imd_recv_mdcomm(IMDSocket* socket, int32_t nforces, int32_t* forcendx, float* forces)
{
    /* reading indices */
    int retsize  = sizeof(int32_t) * nforces;
    int retbytes = imd_read_multiple(socket, reinterpret_cast<char*>(forcendx), retsize);
    if (retbytes != retsize)
    {
        return false;
    }

    /* reading forces as float array */
    retsize  = 3 * sizeof(float) * nforces;
    retbytes = imd_read_multiple(socket, reinterpret_cast<char*>(forces), retsize);
    return (retbytes == retsize);
}

/* GROMACS specific functions for the IMD implementation */
void write_IMDgroup_to_file(bool              bIMD,
                            t_inputrec*       ir,
                            const t_state*    state,
                            const gmx_mtop_t& sys,
                            int               nfile,
                            const t_filenm    fnm[])
{
    t_atoms IMDatoms;


    if (bIMD)
    {
        IMDatoms = gmx_mtop_global_atoms(sys);
        write_sto_conf_indexed(opt2fn("-imd", nfile, fnm),
                               "IMDgroup",
                               &IMDatoms,
                               state->x.rvec_array(),
                               state->v.rvec_array(),
                               ir->pbcType,
                               state->box,
                               ir->imd->nat,
                               ir->imd->ind);
    }
}


void ImdSession::dd_make_local_IMD_atoms(const gmx_domdec_t* dd)
{
    if (!impl_->sessionPossible)
    {
        return;
    }

    dd_make_local_group_indices(
            dd->ga2la.get(), impl_->nat, impl_->ind, &impl_->nat_loc, &impl_->ind_loc, &impl_->nalloc_loc, impl_->xa_ind);
}


/*! \brief Send positions from rvec.
 *
 * We need a separate send buffer and conversion to Angstrom.
 */
static int imd_send_rvecs(IMDSocket* socket, int nat, rvec* x, char* buffer)
{
    int32_t size;
    int     i;
    float   sendx[3];
    int     tuplesize = 3 * sizeof(float);


    /* Required size for the send buffer */
    size = c_headerSize + 3 * sizeof(float) * nat;

    /* Prepare header */
    fill_header(reinterpret_cast<IMDHeader*>(buffer), IMDMessageType::FCoords, static_cast<int32_t>(nat));
    for (i = 0; i < nat; i++)
    {
        sendx[0] = static_cast<float>(x[i][0]) * gmx::c_nm2A;
        sendx[1] = static_cast<float>(x[i][1]) * gmx::c_nm2A;
        sendx[2] = static_cast<float>(x[i][2]) * gmx::c_nm2A;
        memcpy(buffer + c_headerSize + i * tuplesize, sendx, tuplesize);
    }

    return static_cast<int>(imd_write_multiple(socket, buffer, size) != size);
}


void ImdSession::Impl::prepareMainSocket()
{
    if (imdsock_winsockinit() == -1)
    {
        gmx_fatal(FARGS, "%s Failed to initialize winsock.\n", IMDstr);
    }

    /* The rest is identical, first create and bind a socket and set to listen then. */
    GMX_LOG(mdLog_.warning).appendTextFormatted("%s Setting up incoming socket.", IMDstr);
    socket = imdsock_create();
    if (!socket)
    {
        gmx_fatal(FARGS, "%s Failed to create socket.", IMDstr);
    }

    /* Bind to port */
    int ret = imdsock_bind(socket, port);
    if (ret)
    {
        gmx_fatal(FARGS, "%s binding socket to port %d failed with error %d.\n", IMDstr, port, ret);
    }

    if (imd_sock_listen(socket))
    {
        gmx_fatal(FARGS, "%s socket listen failed with error %d.\n", IMDstr, ret);
    }

    if (imdsock_getport(socket, &port))
    {
        gmx_fatal(FARGS, "%s Could not determine port number.\n", IMDstr);
    }

    GMX_LOG(mdLog_.warning).appendTextFormatted("%s Listening for IMD connection on port %d.", IMDstr, port);
}


void ImdSession::Impl::disconnectClient()
{
    /* Write out any buffered pulling data */
    fflush(outf);

    /* we first try to shut down the clientsocket */
    imdsock_shutdown(clientsocket);
    if (!imdsock_destroy(clientsocket))
    {
        GMX_LOG(mdLog_.warning).appendTextFormatted("%s Failed to destroy socket.", IMDstr);
    }

    /* then we reset the IMD step to its default, and reset the connection boolean */
    nstimd_new   = defaultNstImd;
    clientsocket = nullptr;
    bConnected   = false;
}


void ImdSession::Impl::issueFatalError(const char* msg)
{
    GMX_LOG(mdLog_.warning).appendTextFormatted("%s %s", IMDstr, msg);
    disconnectClient();
    GMX_LOG(mdLog_.warning).appendTextFormatted("%s disconnected.", IMDstr);
}


bool ImdSession::Impl::tryConnect()
{
    if (imdsock_tryread(socket, 0, 0) > 0)
    {
        /* yes, we got something, accept on clientsocket */
        clientsocket = imdsock_accept(socket);
        if (!clientsocket)
        {
            GMX_LOG(mdLog_.warning)
                    .appendTextFormatted("%s Accepting the connection on the socket failed.", IMDstr);
            return false;
        }

        /* handshake with client */
        if (imd_handshake(clientsocket))
        {
            issueFatalError("Connection failed.");
            return false;
        }

        GMX_LOG(mdLog_.warning)
                .appendTextFormatted("%s Connection established, checking if I got IMD_GO orders.", IMDstr);

        /* Check if we get the proper "GO" command from client. */
        if (imdsock_tryread(clientsocket, c_connectWait, 0) != 1
            || imd_recv_header(clientsocket, &(length)) != IMDMessageType::Go)
        {
            issueFatalError("No IMD_GO order received. IMD connection failed.");
        }

        /* IMD connected */
        bConnected = true;

        return true;
    }

    return false;
}


void ImdSession::Impl::blockConnect()
{
    /* do not wait for connection, when e.g. ctrl+c is pressed and we will terminate anyways. */
    if (!(gmx_get_stop_condition() == StopCondition::None))
    {
        return;
    }

    GMX_LOG(mdLog_.warning)
            .appendTextFormatted("%s Will wait until I have a connection and IMD_GO orders.", IMDstr);

    /* while we have no clientsocket... 2nd part: we should still react on ctrl+c */
    while ((!clientsocket) && (gmx_get_stop_condition() == StopCondition::None))
    {
        tryConnect();
        imd_sleep(c_loopWait);
    }
}


void ImdSession::Impl::prepareVmdForces()
{
    srenew((vmd_f_ind), vmd_nforces);
    srenew((vmd_forces), 3 * vmd_nforces);
}


void ImdSession::Impl::readVmdForces()
{
    /* the length of the previously received header tells us the nr of forces we will receive */
    vmd_nforces = length;
    /* prepare the arrays */
    prepareVmdForces();
    /* Now we read the forces... */
    if (!(imd_recv_mdcomm(clientsocket, vmd_nforces, vmd_f_ind, vmd_forces)))
    {
        issueFatalError("Error while reading forces from remote. Disconnecting");
    }
}


void ImdSession::Impl::prepareMDForces()
{
    srenew((f_ind), nforces);
    srenew((f), nforces);
}


void ImdSession::Impl::copyToMDForces() const
{
    int  i;
    real conversion = gmx::c_cal2Joule * gmx::c_nm2A;


    for (i = 0; i < nforces; i++)
    {
        /* Copy the indices, a copy is important because we may update the incoming forces
         * whenever we receive new forces while the MD forces are only communicated upon
         * IMD communication.*/
        f_ind[i] = vmd_f_ind[i];

        /* Convert to rvecs and do a proper unit conversion */
        f[i][0] = vmd_forces[3 * i] * conversion;
        f[i][1] = vmd_forces[3 * i + 1] * conversion;
        f[i][2] = vmd_forces[3 * i + 2] * conversion;
    }
}


/*! \brief Returns true if any component of the two rvecs differs. */
static inline bool rvecs_differ(const rvec v1, const rvec v2)
{
    for (int i = 0; i < DIM; i++)
    {
        if (v1[i] != v2[i])
        {
            return true;
        }
    }

    return false;
}

bool ImdSession::Impl::bForcesChanged() const
{
    /* First, check whether the number of pulled atoms changed */
    if (nforces != old_nforces)
    {
        return true;
    }

    /* Second, check whether any of the involved atoms changed */
    for (int i = 0; i < nforces; i++)
    {
        if (f_ind[i] != old_f_ind[i])
        {
            return true;
        }
    }

    /* Third, check whether all forces are the same */
    for (int i = 0; i < nforces; i++)
    {
        if (rvecs_differ(f[i], old_forces[i]))
        {
            return true;
        }
    }

    /* All old and new forces are identical! */
    return false;
}


void ImdSession::Impl::keepOldValues()
{
    old_nforces = nforces;

    for (int i = 0; i < nforces; i++)
    {
        old_f_ind[i] = f_ind[i];
        copy_rvec(f[i], old_forces[i]);
    }
}


void ImdSession::Impl::outputForces(double time)
{
    if (!bForcesChanged())
    {
        return;
    }

    /* Write time and total number of applied IMD forces */
    fprintf(outf, "%14.6e%6d", time, nforces);

    /* Write out the global atom indices of the pulled atoms and the forces itself,
     * write out a force only if it has changed since the last output */
    for (int i = 0; i < nforces; i++)
    {
        if (rvecs_differ(f[i], old_forces[i]))
        {
            fprintf(outf, "%9d", ind[f_ind[i]] + 1);
            fprintf(outf, "%12.4e%12.4e%12.4e", f[i][0], f[i][1], f[i][2]);
        }
    }
    fprintf(outf, "\n");

    keepOldValues();
}


void ImdSession::Impl::syncNodes(const t_commrec* cr, double t)
{
    /* Notify the other nodes whether we are still connected. */
    if (PAR(cr))
    {
        block_bc(cr->mpi_comm_mygroup, bConnected);
    }

    /* ...if not connected, the job is done here. */
    if (!bConnected)
    {
        return;
    }

    /* Let the other nodes know whether we got a new IMD synchronization frequency. */
    if (PAR(cr))
    {
        block_bc(cr->mpi_comm_mygroup, nstimd_new);
    }

    /* Now we all set the (new) nstimd communication time step */
    nstimd = nstimd_new;

    /* We're done if we don't allow pulling at all */
    if (!(bForceActivated))
    {
        return;
    }

    /* OK, let's check if we have received forces which we need to communicate
     * to the other nodes */
    int new_nforces = 0;
    if (MAIN(cr))
    {
        if (bNewForces)
        {
            new_nforces = vmd_nforces;
        }
        /* make the "new_forces" negative, when we did not receive new ones */
        else
        {
            new_nforces = vmd_nforces * -1;
        }
    }

    /* make new_forces known to the clients */
    if (PAR(cr))
    {
        block_bc(cr->mpi_comm_mygroup, new_nforces);
    }

    /* When new_natoms < 0 then we know that these are still the same forces
     * so we don't communicate them, otherwise... */
    if (new_nforces < 0)
    {
        return;
    }

    /* set local VMD and nforces */
    vmd_nforces = new_nforces;
    nforces     = vmd_nforces;

    /* now everybody knows the number of forces in f_ind, so we can prepare
     * the target arrays for indices and forces */
    prepareMDForces();

    /* we first update the MD forces on the main by converting the VMD forces */
    if (MAIN(cr))
    {
        copyToMDForces();
        /* We also write out forces on every update, so that we know which
         * forces are applied for every step */
        if (outf)
        {
            outputForces(t);
        }
    }

    /* In parallel mode we communicate the to-be-applied forces to the other nodes */
    if (PAR(cr))
    {
        nblock_bc(cr->mpi_comm_mygroup, nforces, f_ind);
        nblock_bc(cr->mpi_comm_mygroup, nforces, f);
    }

    /* done communicating the forces, reset bNewForces */
    bNewForces = false;
}


void ImdSession::Impl::readCommand()
{
    bool IMDpaused = false;

    while (clientsocket && (imdsock_tryread(clientsocket, 0, 0) > 0 || IMDpaused))
    {
        IMDMessageType itype = imd_recv_header(clientsocket, &(length));
        /* let's see what we got: */
        switch (itype)
        {
            /* IMD asks us to terminate the simulation, check if the user allowed this */
            case IMDMessageType::Kill:
                if (bTerminatable)
                {
                    GMX_LOG(mdLog_.warning)
                            .appendTextFormatted(
                                    " %s Terminating connection and running simulation (if "
                                    "supported by integrator).",
                                    IMDstr);
                    bWConnect = false;
                    gmx_set_stop_condition(StopCondition::Next);
                }
                else
                {
                    GMX_LOG(mdLog_.warning)
                            .appendTextFormatted(
                                    " %s Set -imdterm command line switch to allow mdrun "
                                    "termination from within IMD.",
                                    IMDstr);
                }

                break;

            /* the client doen't want to talk to us anymore */
            case IMDMessageType::Disconnect:
                GMX_LOG(mdLog_.warning).appendTextFormatted(" %s Disconnecting client.", IMDstr);
                disconnectClient();
                break;

            /* we got new forces, read them and set bNewForces flag */
            case IMDMessageType::Mdcomm:
                readVmdForces();
                bNewForces = true;
                break;

            /* the client asks us to (un)pause the simulation. So we toggle the IMDpaused state */
            case IMDMessageType::Pause:
                if (IMDpaused)
                {
                    GMX_LOG(mdLog_.warning).appendTextFormatted(" %s Un-pause command received.", IMDstr);
                    IMDpaused = false;
                }
                else
                {
                    GMX_LOG(mdLog_.warning).appendTextFormatted(" %s Pause command received.", IMDstr);
                    IMDpaused = true;
                }

                break;

            /* the client sets a new transfer rate, if we get 0, we reset the rate
             * to the default. VMD filters 0 however */
            case IMDMessageType::TRate:
                nstimd_new = (length > 0) ? length : defaultNstImd;
                GMX_LOG(mdLog_.warning)
                        .appendTextFormatted(" %s Update frequency will be set to %d.", IMDstr, nstimd_new);
                break;

            /* Catch all rule for the remaining IMD types which we don't expect */
            default:
                GMX_LOG(mdLog_.warning)
                        .appendTextFormatted(
                                " %s Received unexpected %s.", IMDstr, enumValueToString(itype));
                issueFatalError("Terminating connection");
                break;
        } /* end switch */
    }     /* end while  */
}


void ImdSession::Impl::openOutputFile(const char*                 fn,
                                      int                         nat_total,
                                      const gmx_output_env_t*     oenv,
                                      const gmx::StartingBehavior startingBehavior)
{
    /* Open log file of applied IMD forces if requested */
    if (!fn || !oenv)
    {
        fprintf(stdout,
                "%s For a log of the IMD pull forces explicitly specify '-if' on the command "
                "line.\n"
                "%s (Not possible with energy minimization.)\n",
                IMDstr,
                IMDstr);
        return;
    }

    /* If we append to an existing file, all the header information is already there */
    if (startingBehavior == StartingBehavior::RestartWithAppending)
    {
        outf = gmx_fio_fopen(fn, "a+");
    }
    else
    {
        outf = gmx_fio_fopen(fn, "w+");
        if (nat == nat_total)
        {
            fprintf(outf,
                    "# Note that you can select an IMD index group in the .mdp file if a subset of "
                    "the atoms suffices.\n");
        }

        xvgr_header(outf, "IMD Pull Forces", "Time (ps)", "# of Forces / Atom IDs / Forces (kJ/mol)", exvggtNONE, oenv);

        fprintf(outf, "# Can display and manipulate %d (of a total of %d) atoms via IMD.\n", nat, nat_total);
        fprintf(outf, "# column 1    : time (ps)\n");
        fprintf(outf,
                "# column 2    : total number of atoms feeling an IMD pulling force at that "
                "time\n");
        fprintf(outf,
                "# cols. 3.-6  : global atom number of pulled atom, x-force, y-force, z-force "
                "(kJ/mol)\n");
        fprintf(outf,
                "# then follow : atom-ID, f[x], f[y], f[z] for more atoms in case the force on "
                "multiple atoms is changed simultaneously.\n");
        fprintf(outf,
                "# Note that the force on any atom is always equal to the last value for that "
                "atom-ID found in the data.\n");
        fflush(outf);
    }

    /* To reduce the output file size we remember the old values and output only
     * when something changed */
    snew(old_f_ind, nat); /* One can never pull on more atoms */
    snew(old_forces, nat);
}


ImdSession::Impl::Impl(const MDLogger& mdlog) : mdLog_(mdlog)
{
    init_block(&mols);
}

ImdSession::Impl::~Impl()
{
    if (outf)
    {
        gmx_fio_fclose(outf);
    }
    done_block(&mols);
}


void ImdSession::Impl::prepareMoleculesInImdGroup(const gmx_mtop_t& top_global)
{
    /* check whether index is sorted */
    for (int i = 0; i < nat - 1; i++)
    {
        if (ind[i] > ind[i + 1])
        {
            gmx_fatal(FARGS, "%s IMD index is not sorted. This is currently not supported.\n", IMDstr);
        }
    }

    RangePartitioning gmols = gmx_mtop_molecules(top_global);
    t_block           lmols;
    lmols.nr = 0;
    snew(lmols.index, gmols.numBlocks() + 1);
    lmols.index[0] = 0;

    for (int i = 0; i < gmols.numBlocks(); i++)
    {
        auto mol   = gmols.block(i);
        int  count = 0;
        for (int ii = 0; ii < nat; ii++)
        {
            if (mol.isInRange(ind[ii]))
            {
                count += 1;
            }
        }
        if (count > 0)
        {
            lmols.index[lmols.nr + 1] = lmols.index[lmols.nr] + count;
            lmols.nr += 1;
        }
    }
    srenew(lmols.index, lmols.nr + 1);
    lmols.nalloc_index = lmols.nr + 1;
    mols               = lmols;
}


/*! \brief Copied and modified from groupcoord.c shift_positions_group(). */
static void shift_positions(const matrix box,
                            rvec         x[], /* The positions [0..nr] */
                            const ivec   is,  /* The shift [0..nr] */
                            int          nr)           /* The number of positions */
{
    int i, tx, ty, tz;

    /* Loop over the group's atoms */
    if (TRICLINIC(box))
    {
        for (i = 0; i < nr; i++)
        {
            tx = is[XX];
            ty = is[YY];
            tz = is[ZZ];

            x[i][XX] = x[i][XX] - tx * box[XX][XX] - ty * box[YY][XX] - tz * box[ZZ][XX];
            x[i][YY] = x[i][YY] - ty * box[YY][YY] - tz * box[ZZ][YY];
            x[i][ZZ] = x[i][ZZ] - tz * box[ZZ][ZZ];
        }
    }
    else
    {
        for (i = 0; i < nr; i++)
        {
            tx = is[XX];
            ty = is[YY];
            tz = is[ZZ];

            x[i][XX] = x[i][XX] - tx * box[XX][XX];
            x[i][YY] = x[i][YY] - ty * box[YY][YY];
            x[i][ZZ] = x[i][ZZ] - tz * box[ZZ][ZZ];
        }
    }
}


void ImdSession::Impl::removeMolecularShifts(const matrix box) const
{
    /* for each molecule also present in IMD group */
    for (int i = 0; i < mols.nr; i++)
    {
        /* first we determine the minimum and maximum shifts for each molecule */

        ivec largest, smallest, shift;
        clear_ivec(largest);
        clear_ivec(smallest);
        clear_ivec(shift);

        copy_ivec(xa_shifts[mols.index[i]], largest);
        copy_ivec(xa_shifts[mols.index[i]], smallest);

        for (int ii = mols.index[i] + 1; ii < mols.index[i + 1]; ii++)
        {
            if (xa_shifts[ii][XX] > largest[XX])
            {
                largest[XX] = xa_shifts[ii][XX];
            }
            if (xa_shifts[ii][XX] < smallest[XX])
            {
                smallest[XX] = xa_shifts[ii][XX];
            }

            if (xa_shifts[ii][YY] > largest[YY])
            {
                largest[YY] = xa_shifts[ii][YY];
            }
            if (xa_shifts[ii][YY] < smallest[YY])
            {
                smallest[YY] = xa_shifts[ii][YY];
            }

            if (xa_shifts[ii][ZZ] > largest[ZZ])
            {
                largest[ZZ] = xa_shifts[ii][ZZ];
            }
            if (xa_shifts[ii][ZZ] < smallest[ZZ])
            {
                smallest[ZZ] = xa_shifts[ii][ZZ];
            }
        }

        /* check if we what we can subtract/add to the positions
         * to put them back into the central box */
        if (smallest[XX] > 0)
        {
            shift[XX] = smallest[XX];
        }
        if (smallest[YY] > 0)
        {
            shift[YY] = smallest[YY];
        }
        if (smallest[ZZ] > 0)
        {
            shift[ZZ] = smallest[ZZ];
        }

        if (largest[XX] < 0)
        {
            shift[XX] = largest[XX];
        }
        if (largest[YY] < 0)
        {
            shift[YY] = largest[YY];
        }
        if (largest[ZZ] < 0)
        {
            shift[ZZ] = largest[ZZ];
        }

        /* is there a shift at all? */
        if ((shift[XX]) || (shift[YY]) || (shift[ZZ]))
        {
            int molsize = mols.index[i + 1] - mols.index[i];
            /* shift the positions */
            shift_positions(box, &(xa[mols.index[i]]), shift, molsize);
        }
    }
}


void ImdSession::Impl::prepareForPositionAssembly(const t_commrec* cr, gmx::ArrayRef<const gmx::RVec> coords)
{
    snew(xa, nat);
    snew(xa_ind, nat);
    snew(xa_shifts, nat);
    snew(xa_eshifts, nat);
    snew(xa_old, nat);

    /* Save the original (whole) set of positions such that later the
     * molecule can always be made whole again */
    if (MAIN(cr))
    {
        for (int i = 0; i < nat; i++)
        {
            int ii = ind[i];
            copy_rvec(coords[ii], xa_old[i]);
        }
    }

    if (!haveDDAtomOrdering(*cr))
    {
        nat_loc = nat;
        ind_loc = ind;

        /* xa_ind[i] needs to be set to i for serial runs */
        for (int i = 0; i < nat; i++)
        {
            xa_ind[i] = i;
        }
    }

    /* Communicate initial coordinates xa_old to all processes */
    if (cr && havePPDomainDecomposition(cr))
    {
        gmx_bcast(nat * sizeof(xa_old[0]), xa_old, cr->mpi_comm_mygroup);
    }
}


/*! \brief Check for non-working integrator / parallel options. */
static void imd_check_integrator_parallel(const t_inputrec* ir, const t_commrec* cr)
{
    if (PAR(cr))
    {
        if (((ir->eI) == IntegrationAlgorithm::Steep) || ((ir->eI) == IntegrationAlgorithm::CG)
            || ((ir->eI) == IntegrationAlgorithm::LBFGS) || ((ir->eI) == IntegrationAlgorithm::NM))
        {
            gmx_fatal(FARGS,
                      "%s Energy minimization via steep, CG, lbfgs and nm in parallel is currently "
                      "not supported by IMD.\n",
                      IMDstr);
        }
    }
}

std::unique_ptr<ImdSession> makeImdSession(const t_inputrec*              ir,
                                           const t_commrec*               cr,
                                           gmx_wallcycle*                 wcycle,
                                           gmx_enerdata_t*                enerd,
                                           const gmx_multisim_t*          ms,
                                           const gmx_mtop_t&              top_global,
                                           const MDLogger&                mdlog,
                                           gmx::ArrayRef<const gmx::RVec> coords,
                                           int                            nfile,
                                           const t_filenm                 fnm[],
                                           const gmx_output_env_t*        oenv,
                                           const ImdOptions&              options,
                                           const gmx::StartingBehavior    startingBehavior)
{
    std::unique_ptr<ImdSession> session(new ImdSession(mdlog));
    auto*                       impl = session->impl_.get();

    /* We will allow IMD sessions only if supported by the binary and
       explicitly enabled in the .tpr file */
    if (!GMX_IMD || !ir->bIMD)
    {
        return session;
    }

    // TODO As IMD is intended for interactivity, and the .tpr file
    // opted in for that, it is acceptable to write more terminal
    // output than in a typical simulation. However, all the GMX_LOG
    // statements below should go to both the log file and to the
    // terminal. This is probably be implemented by adding a logging
    // stream named like ImdInfo, to separate it from warning and to
    // send it to both destinations.

    if (EI_DYNAMICS(ir->eI))
    {
        impl->defaultNstImd = ir->nstcalcenergy;
    }
    else if (EI_ENERGY_MINIMIZATION(ir->eI))
    {
        impl->defaultNstImd = 1;
    }
    else
    {
        GMX_LOG(mdlog.warning)
                .appendTextFormatted(
                        "%s Integrator '%s' is not supported for Interactive Molecular Dynamics, "
                        "running normally instead",
                        IMDstr,
                        enumValueToString(ir->eI));
        return session;
    }
    if (isMultiSim(ms))
    {
        GMX_LOG(mdlog.warning)
                .appendTextFormatted(
                        "%s Cannot use IMD for multiple simulations or replica exchange, running "
                        "normally instead",
                        IMDstr);
        return session;
    }

    bool createSession = false;
    /* It seems we have a .tpr file that defines an IMD group and thus allows IMD connections.
     * Check whether we can actually provide the IMD functionality for this setting: */
    if (MAIN(cr))
    {
        /* Check whether IMD was enabled by one of the command line switches: */
        if (options.wait || options.terminatable || options.pull)
        {
            GMX_LOG(mdlog.warning)
                    .appendTextFormatted(
                            "%s Enabled. This simulation will accept incoming IMD connections.", IMDstr);
            createSession = true;
        }
        else
        {
            GMX_LOG(mdlog.warning)
                    .appendTextFormatted(
                            "%s None of the -imd switches was used.\n"
                            "%s This run will not accept incoming IMD connections",
                            IMDstr,
                            IMDstr);
        }
    } /* end main only */

    /* Let the other nodes know whether we want IMD */
    if (PAR(cr))
    {
        block_bc(cr->mpi_comm_mygroup, createSession);
    }

    /*... if not we are done.*/
    if (!createSession)
    {
        return session;
    }


    /* check if we're using a sane integrator / parallel combination */
    imd_check_integrator_parallel(ir, cr);


    /*
     *****************************************************
     * From here on we assume that IMD is turned on      *
     *****************************************************
     */

    int nat_total = top_global.natoms;

    /* Initialize IMD session. If we read in a pre-IMD .tpr file, ir->imd->nat
     * will be zero. For those cases we transfer _all_ atomic positions */
    impl->sessionPossible = true;
    impl->nat             = ir->imd->nat > 0 ? ir->imd->nat : nat_total;
    if (options.port >= 1)
    {
        impl->port = options.port;
    }
    impl->cr_    = cr;
    impl->wcycle = wcycle;
    impl->enerd  = enerd;

    /* We might need to open an output file for IMD forces data */
    if (MAIN(cr))
    {
        impl->openOutputFile(opt2fn("-if", nfile, fnm), nat_total, oenv, startingBehavior);
    }

    /* Make sure that we operate with a valid atom index array for the IMD atoms */
    if (ir->imd->nat > 0)
    {
        /* Point to the user-supplied array of atom numbers */
        impl->ind = ir->imd->ind;
    }
    else
    {
        /* Make a dummy (ind[i] = i) array of all atoms */
        snew(impl->ind, nat_total);
        for (int i = 0; i < nat_total; i++)
        {
            impl->ind[i] = i;
        }
    }

    /* read environment on main and prepare socket for incoming connections */
    if (MAIN(cr))
    {
        /* we allocate memory for our IMD energy structure */
        int32_t recsize = c_headerSize + sizeof(IMDEnergyBlock);
        snew(impl->energysendbuf, recsize);

        /* Shall we wait for a connection? */
        if (options.wait)
        {
            impl->bWConnect = true;
            GMX_LOG(mdlog.warning)
                    .appendTextFormatted(
                            "%s Pausing simulation while no IMD connection present (-imdwait).", IMDstr);
        }

        /* Will the IMD clients be able to terminate the simulation? */
        if (options.terminatable)
        {
            impl->bTerminatable = true;
            GMX_LOG(mdlog.warning)
                    .appendTextFormatted(
                            "%s Allow termination of the simulation from IMD client (-imdterm).", IMDstr);
        }

        /* Is pulling from IMD client allowed? */
        if (options.pull)
        {
            impl->bForceActivated = true;
            GMX_LOG(mdlog.warning)
                    .appendTextFormatted("%s Pulling from IMD remote is enabled (-imdpull).", IMDstr);
        }

        /* Initialize send buffers with constant size */
        snew(impl->sendxbuf, impl->nat);
        snew(impl->energies, 1);
        int32_t bufxsize = c_headerSize + 3 * sizeof(float) * impl->nat;
        snew(impl->coordsendbuf, bufxsize);
    }

    /* do we allow interactive pulling? If so let the other nodes know. */
    if (PAR(cr))
    {
        block_bc(cr->mpi_comm_mygroup, impl->bForceActivated);
    }

    /* setup the listening socket on main process */
    if (MAIN(cr))
    {
        GMX_LOG(mdlog.warning).appendTextFormatted("%s Setting port for connection requests to %d.", IMDstr, impl->port);
        impl->prepareMainSocket();
        /* Wait until we have a connection if specified before */
        if (impl->bWConnect)
        {
            impl->blockConnect();
        }
        else
        {
            GMX_LOG(mdlog.warning).appendTextFormatted("%s -imdwait not set, starting simulation.", IMDstr);
        }
    }
    /* Let the other nodes know whether we are connected */
    impl->syncNodes(cr, 0);

    /* Initialize arrays used to assemble the positions from the other nodes */
    impl->prepareForPositionAssembly(cr, coords);

    /* Initialize molecule blocks to make them whole later...*/
    if (MAIN(cr))
    {
        impl->prepareMoleculesInImdGroup(top_global);
    }

    return session;
}


bool ImdSession::Impl::run(int64_t step, bool bNS, const matrix box, gmx::ArrayRef<const gmx::RVec> coords, double t)
{
    /* IMD at all? */
    if (!sessionPossible)
    {
        return false;
    }

    wallcycle_start(wcycle, WallCycleCounter::Imd);

    /* read command from client and check if new incoming connection */
    if (MAIN(cr_))
    {
        /* If not already connected, check for new connections */
        if (!clientsocket)
        {
            if (bWConnect)
            {
                blockConnect();
            }
            else
            {
                tryConnect();
            }
        }

        /* Let's see if we have new IMD messages for us */
        if (clientsocket)
        {
            readCommand();
        }
    }

    /* is this an IMD communication step? */
    bool imdstep = do_per_step(step, nstimd);

    /* OK so this is an IMD step ... */
    if (imdstep)
    {
        /* First we sync all nodes to let everybody know whether we are connected to VMD */
        syncNodes(cr_, t);
    }

    /* If a client is connected, we collect the positions
     * and put molecules back into the box before transfer */
    if ((imdstep && bConnected) || bNS) /* independent of imdstep, we communicate positions at each NS step */
    {
        /* Transfer the IMD positions to the main node. Every node contributes
         * its local positions x and stores them in the assembled xa array. */
        communicate_group_positions(
                cr_, xa, xa_shifts, xa_eshifts, true, as_rvec_array(coords.data()), nat, nat_loc, ind_loc, xa_ind, xa_old, box);

        /* If connected and main -> remove shifts */
        if ((imdstep && bConnected) && MAIN(cr_))
        {
            removeMolecularShifts(box);
        }
    }

    wallcycle_stop(wcycle, WallCycleCounter::Imd);

    return imdstep;
}

bool ImdSession::run(int64_t step, bool bNS, const matrix box, gmx::ArrayRef<const gmx::RVec> coords, double t)
{
    return impl_->run(step, bNS, box, coords, t);
}

void ImdSession::fillEnergyRecord(int64_t step, bool bHaveNewEnergies)
{
    if (!impl_->sessionPossible || !impl_->clientsocket)
    {
        return;
    }

    IMDEnergyBlock* ene = impl_->energies;

    ene->tstep = step;

    /* In MPI-parallel simulations the energies are not accessible a at every time step.
     * We update them if we have new values, otherwise, the energy values from the
     * last global communication step are still on display in the viewer. */
    if (bHaveNewEnergies)
    {
        ene->T_abs   = static_cast<float>(impl_->enerd->term[F_TEMP]);
        ene->E_pot   = static_cast<float>(impl_->enerd->term[F_EPOT]);
        ene->E_tot   = static_cast<float>(impl_->enerd->term[F_ETOT]);
        ene->E_bond  = static_cast<float>(impl_->enerd->term[F_BONDS]);
        ene->E_angle = static_cast<float>(impl_->enerd->term[F_ANGLES]);
        ene->E_dihe  = static_cast<float>(impl_->enerd->term[F_PDIHS]);
        ene->E_impr  = static_cast<float>(impl_->enerd->term[F_IDIHS]);
        ene->E_vdw   = static_cast<float>(impl_->enerd->term[F_LJ]);
        ene->E_coul  = static_cast<float>(impl_->enerd->term[F_COUL_SR]);
    }
}


void ImdSession::sendPositionsAndEnergies()
{
    if (!impl_->sessionPossible || !impl_->clientsocket)
    {
        return;
    }

    if (imd_send_energies(impl_->clientsocket, impl_->energies, impl_->energysendbuf))
    {
        impl_->issueFatalError("Error sending updated energies. Disconnecting client.");
    }

    if (imd_send_rvecs(impl_->clientsocket, impl_->nat, impl_->xa, impl_->coordsendbuf))
    {
        impl_->issueFatalError("Error sending updated positions. Disconnecting client.");
    }
}


void ImdSession::updateEnergyRecordAndSendPositionsAndEnergies(bool bIMDstep, int64_t step, bool bHaveNewEnergies)
{
    if (!impl_->sessionPossible)
    {
        return;
    }

    wallcycle_start(impl_->wcycle, WallCycleCounter::Imd);

    /* Update time step for IMD and prepare IMD energy record if we have new energies. */
    fillEnergyRecord(step, bHaveNewEnergies);

    if (bIMDstep)
    {
        /* Send positions and energies to VMD client via IMD */
        sendPositionsAndEnergies();
    }

    wallcycle_stop(impl_->wcycle, WallCycleCounter::Imd);
}

void ImdSession::applyForces(gmx::ArrayRef<gmx::RVec> force)
{
    if (!impl_->sessionPossible || !impl_->bForceActivated)
    {
        return;
    }

    wallcycle_start(impl_->wcycle, WallCycleCounter::Imd);

    for (int i = 0; i < impl_->nforces; i++)
    {
        /* j are the indices in the "System group".*/
        int j = impl_->ind[impl_->f_ind[i]];

        /* check if this is a local atom and find out locndx */
        const int*       locndx;
        const t_commrec* cr = impl_->cr_;
        if (PAR(cr) && (locndx = cr->dd->ga2la->findHome(j)))
        {
            j = *locndx;
        }

        rvec_inc(force[j], impl_->f[i]);
    }

    wallcycle_stop(impl_->wcycle, WallCycleCounter::Imd);
}

ImdSession::ImdSession(const MDLogger& mdlog) : impl_(new Impl(mdlog)) {}

ImdSession::~ImdSession() = default;

} // namespace gmx
