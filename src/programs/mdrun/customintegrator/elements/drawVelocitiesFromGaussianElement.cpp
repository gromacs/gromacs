#include "drawVelocitiesFromGaussianElement.h"

#include "gromacs/math/units.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/math/vec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"

#include "programs/mdrun/customintegrator/statemanager.h"

void DrawVelocitiesFromGaussian::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void DrawVelocitiesFromGaussian::run()
{
    /*
    v ~ N(0, sqrt(kT/m))
    */
    int * gatindex = DOMAINDECOMP(_data->cr) ? _data->cr->dd->gatindex : nullptr;

    /* Note: the temperature information is obtained from mdp
    * It will take the first temperature group's temperature to 
    * draw velocities for the whole system.
    * It might be worth re-thinking what the code should do if 
    * there are multiple temperature groups.
    * Or introducing another temperature field in the mdp for HMC integrators
    * may be more useful*/
    
    real kT = real(BOLTZ*_data->ir->opts.ref_t[0]);
    int  gf = 0;
    int  start  = 0;
    int  homenr = _data->mdatoms->homenr;
    int  nrend  = start+homenr;
    
    int nth = gmx_omp_nthreads_get(emntUpdate);
    // loop over the threads
#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th, end_th;

            start_th = start + ((nrend-start)* th   )/nth;
            end_th   = start + ((nrend-start)*(th+1))/nth;
            
            rvec             *v = as_rvec_array(_data->state->v.data());
            
            // Even 0 bits internal counter gives 2x64 ints (more than enough for three table lookups)
            gmx::ThreeFry2x64<0> rng(unsigned(_data->ir->ld_seed),
                                     gmx::RandomDomain::UpdateCoordinates);
            gmx::TabulatedNormalDistribution<real, 14> gauss;

            for (int n = start_th; n < end_th; n++)
            {
                int  ng = gatindex ? gatindex[n] : n;

                rng.restart(_data->step, ng);
                gauss.reset();
                if (_data->mdatoms->cFREEZE)
                {    
                    gf = _data->mdatoms->cFREEZE[n];
                }
                real sigma = std::sqrt(_data->mdatoms->invmass[n]*kT);
                for (int d = 0; d < DIM; d++)
                {
                    if ((_data->mdatoms->ptype[n] != eptVSite) &&
                        (_data->mdatoms->ptype[n] != eptShell) &&
                        !_data->ir->opts.nFreeze[gf][d])
                    {
                        v[n][d] = gauss(rng)*sigma;
                    }
                    else
                    {
                        v[n][d] = 0.0;
                    }
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }  

}