#include "computeGlobalsElement.h"

#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/utility/basedefinitions.h"

#include "programs/mdrun/customintegrator/statemanager.h"
#include "programs/mdrun/customintegrator/helperFunctionsTemporary.h"


ComputeGlobals::ComputeGlobals(gmx_bool bEnergy,
                                gmx_bool bReadEkin,
                                gmx_bool bScaleEkin,
                                gmx_bool bEkinAveVel,
                                gmx_bool bStopCM,
                                gmx_bool bTemp,
                                gmx_bool bPres,
                                gmx_bool bConstraint,
                                bool bInterSimSignal,
                                bool bIntraSimSignal): bIntraSimSignal(bIntraSimSignal),
                                                       bInterSimSignal(bInterSimSignal)
{
    /*
    bEnergy     : Calculate dispersion correction and add it to the potential energy
    bReadEkin   : If TRUE will make bEkinAveVel=TRUE
    bScaleEkin  : When bEkinAveVel=True, The Kinetic Energies are scaled only when bScaleEkin=False ?? 
    bStopCM     : Remove com motion from velocities, depending on the vcm option
    bTemp       : Calculate the kinetic energy and temperature
    bPres       : Calculate dispersion correction and pressure
    bConstraint : Does the same thing as bPres (May not be neccessary)
    bEkinAveVel : Determines the way the Kinetic Energy is calculated
        * if bEkinAveVel is FALSE (may used if only halfstep velocities are present)
            the Kinetic energy is calculated via the average of the halfstep kinetic energies
        * if bEkinAveVel is TRUE (may used if the full step velocity is present)
            the Kinetic energy is calculated from the given velocities
    */
    cglo_flags =  (bEnergy ? CGLO_ENERGY : 0)
                | (bReadEkin ? CGLO_READEKIN : 0)
                | (bScaleEkin ? CGLO_SCALEEKIN : 0)
                | (bStopCM ? CGLO_STOPCM : 0)
                | (bTemp ? CGLO_TEMPERATURE : 0)
                | (bPres ? CGLO_PRESSURE : 0)
                | (bConstraint ? CGLO_CONSTRAINT : 0)
                | (bEkinAveVel ? CGLO_EKINAVEVEL : 0);
}

void ComputeGlobals::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void ComputeGlobals::run()
{
    cglo_flags = cglo_flags
                 | (_data->bGStat ? CGLO_GSTAT : 0);

    if (cglo_flags & CGLO_GSTAT)
    {
        cglo_flags = cglo_flags 
                     | (_data->shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0);
        // set signaller if inter or intra simulation signalling is required
        gmx::SimulationSignaller signaller(nullptr, nullptr, false, false);
        
        if (bInterSimSignal || bIntraSimSignal)
        { 
            gmx::SimulationSignaller signaller(&_data->signals, _data->cr, bInterSimSignal, bIntraSimSignal);
        }

        compute_globals(_data->fplog, _data->gstat, _data->cr, _data->ir, _data->fr, _data->ekind, _data->state, _data->mdatoms, _data->nrnb, _data->vcm,
                        _data->wcycle, _data->enerd, _data->force_vir, _data->shake_vir, _data->total_vir, _data->pres, _data->mu_tot,
                        _data->constr, &signaller,
                        _data->lastbox,
                        &_data->totalNumberOfBondedInteractions, &_data->bSumEkinhOld,
                        cglo_flags);

        checkNumberOfBondedInteractions(_data->fplog, _data->cr, _data->totalNumberOfBondedInteractions,
                                        _data->top_global, _data->top, _data->state,
                                        &_data->shouldCheckNumberOfBondedInteractions);

        _data->enerd->term[F_ETOT] = _data->enerd->term[F_EPOT] + _data->enerd->term[F_EKIN];
    }
}