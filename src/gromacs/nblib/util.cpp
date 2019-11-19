//
// Created by sebkelle on 19.11.19.
//

#include "util.h"

#include "gromacs/math/units.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/utility/fatalerror.h"

static std::vector<gmx::RVec> low_mspeed(real tempi,  std::vector<real> const& masses, gmx::ThreeFry2x64<>* rng)
{
    int                                    nrdf;
    real                                   boltz;
    real                                   ekin, temp;
    gmx::TabulatedNormalDistribution<real> normalDist;

    std::vector<gmx::RVec> velocities(masses.size());

    boltz = BOLTZ * tempi;
    ekin  = 0.0;
    nrdf  = 0;
    for(size_t i = 0; i < masses.size(); i++) {
        real mass  = masses[i];
        if (mass > 0)
        {
            rng->restart(i, 0);
            real sd = std::sqrt(boltz / mass);
            for (int m = 0; (m < DIM); m++)
            {
                velocities[i][m] = sd * normalDist(*rng);
                ekin += 0.5 * mass * velocities[i][m] * velocities[i][m];
            }
            nrdf += DIM;
        }

    }
    temp = (2.0 * ekin) / (nrdf * BOLTZ);
    if (temp > 0)
    {
        real scal = std::sqrt(tempi / temp);
        for(auto &vel : velocities)
        {
            for (int m = 0; (m < DIM); m++)
            {
                vel[m] *= scal;
            }
        }
    }
    fprintf(stderr, "Velocities were taken from a Maxwell distribution at %g K\n", tempi);
    if (debug)
    {
        fprintf(debug,
                "Velocities were taken from a Maxwell distribution\n"
                "Initial generated temperature: %12.5e (scaled to: %12.5e)\n",
                temp, tempi);
    }

    return velocities;
}

//! generate Velocites from a Maxwell Boltzmann distro, masses should be the
//! same as the ones specified for the Topology object
std::vector<gmx::RVec> generateVelocity(real tempi, unsigned int seed, std::vector<real> const& masses)
{

    if (seed == 0)
    {
        seed = static_cast<int>(gmx::makeRandomSeed());
        fprintf(stderr, "Using random seed %u for generating velocities\n", seed);
    }
    gmx::ThreeFry2x64<> rng(seed, gmx::RandomDomain::MaxwellVelocities);

    return low_mspeed(tempi, masses, &rng);
}

bool checkNumericValues(const std::vector<gmx::RVec> &values)
{
    for (auto val : values)
    {
        for (int m = 0; (m < DIM); m++)
        {
            if(std::isnan(val[m]) or std::isinf(val[m]))
            {
                return false;
            }
        }
    }
    return true;
}
