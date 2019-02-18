#include "gmxpre.h"

#include "gromacs/mdlib/lincs.h"

#include "gromacs/mdlib/lincs_cuda.h"


LincsCuda::LincsCuda(int nAtom,
                     int nIter,
                     int nOrder)
{
    GMX_UNUSED_VALUE(nAtom);
    GMX_UNUSED_VALUE(nIter);
    GMX_UNUSED_VALUE(nOrder);
    GMX_ASSERT(false, "Should not be here.");
}

LincsCuda::~LincsCuda()
{
}

/* externally visible function to perform lincs on GPU */
void LincsCuda::apply(const bool       updateVelocities,
                      const real       invdt,
                      const gmx_bool   bCalcVir,
                      tensor           virialScaled)
{
    GMX_UNUSED_VALUE(updateVelocities);
    GMX_UNUSED_VALUE(invdt);
    GMX_UNUSED_VALUE(bCalcVir);
    GMX_UNUSED_VALUE(virialScaled);
    GMX_ASSERT(false, "Should not be here.");
}



void LincsCuda::setPBC(t_pbc *pbc)
{
    GMX_UNUSED_VALUE(pbc);
    GMX_ASSERT(false, "Should not be here.");
}

void LincsCuda::set(const t_idef     &idef,
                    const t_mdatoms  &md)
{
    GMX_UNUSED_VALUE(idef);
    GMX_UNUSED_VALUE(md);
    GMX_ASSERT(false, "Should not be here.");
}

/*! \brief 
 * Copy coordinates and velocities from provided CPU location to GPU.
 *
 * Copies the coordinates before the integration step (x), coordinates 
 * after the integration step (xp) and velocities (v) from the provided 
 * CPU location to GPU. The data are assumed to be in float3/fvec format 
 * (single precision).
 *
 * \param[in] *x  CPU pointer where coordinates should be copied from. 
 * \param[in] *xp CPU pointer where coordinates should be copied from. 
 * \param[in] *v  CPU pointer where velocities should be copied from. 
 */
void LincsCuda::copyCoordinatesToGpu(const rvec * x, const rvec * xp, const rvec * v)
{
    GMX_UNUSED_VALUE(x);
    GMX_UNUSED_VALUE(xp);
    GMX_UNUSED_VALUE(v);
    GMX_ASSERT(false, "Should not be here.");
}

/*! \brief 
 * Copy coordinates from GPU to provided CPU location.
 *
 * Copies the constrained coordinates to the provided location. The coordinates 
 * are assumed to be in float3/fvec format (single precision).
 *
 * \param[out] *xp CPU pointer where coordinates should be copied to. 
 */
void LincsCuda::copyCoordinatesFromGpu(rvec * xp)
{
    GMX_UNUSED_VALUE(xp);
    GMX_ASSERT(false, "Should not be here.");
}

/*! \brief 
 * Copy velocities from GPU to provided CPU location.
 *
 * The velocities are assumed to be in float3/fvec format (single precision).
 *
 et the internal GPU-memory x, xprime and v pointers.
 640  *
 641  * Data is not copied. The data are assumed to be in float3/fvec format
 642  * (float3 is used internaly, but the data layout should be identical).
 643  *
 644  * \param[in] 
 */
void LincsCuda::copyVelocitiesFromGpu(rvec * v)
{
    GMX_UNUSED_VALUE(v);
    GMX_ASSERT(false, "Should not be here.");
}

/*! \brief 
 * Set the internal GPU-memory x, xprime and v pointers.
 *
 * \param[in] *xpDevice Pointer to the coordinates after integrator update, before update (on GPU)   
 * \param[in] *vDevice  Pointer to the velocities before integrator update (on GPU)   
 */
void LincsCuda::setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice)
{
    GMX_UNUSED_VALUE(xDevice);
    GMX_UNUSED_VALUE(xpDevice);
    GMX_UNUSED_VALUE(vDevice);
    GMX_ASSERT(false, "Should not be here.");
}
