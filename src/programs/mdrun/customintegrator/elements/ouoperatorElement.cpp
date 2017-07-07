#include "ouoperatorElement.h"

#include "tgmath.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/math/units.h"
#include "programs/mdrun/customintegrator/statemanager.h"

void OUOperator::initialize(StateManager& dataRef)
{
    // casts made explicit to ease debugging

    _data = &dataRef;
    unsigned int ngtc = unsigned(_data->ir->opts.ngtc);
    _a.resize(ngtc);
    _sqrt_one_a2_beta.resize(ngtc);
    real dt = real(_data->ir->delta_t);
    for (unsigned int gt = 0; gt < ngtc; ++gt)
    {
        real a = real(exp(-dt * _step_fraction / _data->ir->opts.tau_t[gt]));
        real kT = real(BOLTZ*_data->ir->opts.ref_t[gt]);
        _a.at(gt) = a;
        _sqrt_one_a2_beta.at(gt) = std::sqrt((1 - a*a) * kT);
    }
}

void OUOperator::run()
{
    int * gatindex = DOMAINDECOMP(_data->cr) ? _data->cr->dd->gatindex : nullptr;

    int  start  = 0;
    int  homenr = _data->mdatoms->homenr;
    int  nrend  = start+homenr;

    int nth = gmx_omp_nthreads_get(emntUpdate);

#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th, end_th;

            start_th = start + ((nrend-start)* th   )/nth;
            end_th   = start + ((nrend-start)*(th+1))/nth;

            rvec         *v = as_rvec_array(_data->state->v.data());

            int             gf = 0, gt = 0;
            real            ism;
            int             n, d;

            // Even 0 bits internal counter gives 2x64 ints (more than enough for three table lookups)
            gmx::ThreeFry2x64<0> rng(unsigned(_data->ir->ld_seed),
                                     gmx::RandomDomain::UpdateCoordinates);
            gmx::TabulatedNormalDistribution<real, 14> dist;

            for (n = start_th; n < end_th; n++)
            {
                int  ng = gatindex ? gatindex[n] : n;

                rng.restart(_data->step, ng);
                dist.reset();

                ism = std::sqrt(_data->mdatoms->invmass[n]);

                if (_data->mdatoms->cFREEZE)
                {
                    gf  = _data->mdatoms->cFREEZE[n];
                }
                if (_data->mdatoms->cTC)
                {
                    gt  = _data->mdatoms->cTC[n];
                }

                for (d = 0; d < DIM; d++)
                {
                    if ((_data->mdatoms->ptype[n] != eptVSite) &&
                        (_data->mdatoms->ptype[n] != eptShell) &&
                        !_data->ir->opts.nFreeze[gf][d])
                    {
                        v[n][d] = v[n][d]*_a[gt] +
                                  ism*_sqrt_one_a2_beta[gt]*dist(rng);
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