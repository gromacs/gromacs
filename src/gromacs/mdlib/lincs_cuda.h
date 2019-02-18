#ifndef GMX_MDLIB_LINCS_CUDA_H
#define GMX_MDLIB_LINCS_CUDA_H


#include "gromacs/mdlib/nbnxn_gpu.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdlib/nbnxn_gpu_types.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"

#if defined(__CUDACC__)
#include "cuda.h"
#include "cuda_runtime.h"
#include "gromacs/gpu_utils/gpu_vec.cuh"
#include "gromacs/pbcutil/gpu_pbc.cuh"
#else
typedef struct {
    float x, y, z;
} float3;

struct PbcAiuc
{
    float invBoxDiagZ;
    float boxZX;
    float boxZY;
    float boxZZ;
    float invBoxDiagY;
    float boxYX;
    float boxYY;
    float invBoxDiagX;
    float boxXX;
};
#endif

/*! \brief Class with interfaces and data for CUDA version of LINCS.
 *
 * The class provides major interfaces to constrain bonds using LINCS on GPU.
 * Current implementation is developed for H_Bond constraints. Cant handle constraints triangles. 
 *
 *
 */
class LincsCuda {

    public:
        /*! \brief Constructor.
         *  Initializes objects 
         */
        LincsCuda(int N,
                  int nOrder,
                  int nIter);
        ~LincsCuda();

        void apply(bool       updateVelocities,
                   real       invdt,
                   gmx_bool   bCalcVir,
                   tensor     virialScaled);
        void setPBC(t_pbc *pbc);
        void set(const t_idef &idef,
                    const t_mdatoms      &md);
        /*void set(std::vector<gmx::AtomPair>                        atoms,
                    std::vector<int>                                  blnr,
                    std::vector<int>                                  blbnb,
                    std::vector<real>                                 blmf,
                    std::vector<real, gmx::AlignedAllocator<real>>    bllen,
                    const real                                       *invmass);*/
        void copyCoordinatesToGpu(const rvec * x, const rvec * xp, const rvec * v);
        void copyCoordinatesFromGpu(rvec * xp);
        void copyVelocitiesFromGpu(rvec * v);
        void setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice);
        
    private:
        
        PbcAiuc pbcAiuc;
        
        int    nAtom;                                   //!< Number of atoms
        
        int    nOrder;                                  //!< Order of expansion when inversing the matrix
        int    nIter;                                   //!< Number of iterations used to correct the projection
        
        float3  *xDevice;                                 //!< Coordinates before the timestep (in the GPU memory)
        float3  *xpDevice;                                //!< Coordinates after the timestep, before constraining (on GPU).
                                                        //   These will be changed upon constraining.
        float3  *vDevice;                                 //!< Velocities of atoms (on GPU)
        
        real  *invmassDevice;                          //!< 1/mass for all atoms (GPU)
        float *mlambdaDevice;                          //!< Scaled lagrange multipliers (GPU)

        real  *virialScaledHost;                       //!< Scaled virial tensor (9 reals, GPU)
        real  *virialScaledDevice;                     //!< Scaled virial tensor (9 reals, GPU)
        
        int maxConstraintsNumberSoFar;                  //!< Maximum total of constraints assigned to this class so far.
                                                        //   If the new assigned number is larger, the GPU data arrays are reallocated.
        int nConstraintsThreads;                        //!< Total number of threads, which covers all constraints and gaps in the 
                                                        //   ends of the thread blocks that are nessesary to avoid inter-block 
                                                        //   syncronizations.
    
        std::vector<gmx::AtomPair> constraintsHost;     //!< List of constrained atoms (CPU memory)
        gmx::AtomPair *constraintsDevice;               //!< List of constrained atoms (GPU memory)

        std::vector<float> constraintsR0Host;            //!< Equilibrium distances for the constraints (CPU)
        float *constraintsR0Device;                      //!< Equilibrium distances for the constraints (GPU)
        
        std::vector<int> coupledConstraintsCountsHost;  //!< Number of constraints, coupled with the current one (CPU)
        int *coupledConstraintsCountsDevice;            //!< Number of constraints, coupled with the current one (GPU)
        
        std::vector<int> coupledConstraintsIdxesHost;   //!< List of coupled with the current one (CPU)
        int *coupledConstraintsIdxesDevice;             //!< List of coupled with the current one (GPU)
        
        float *matrixADevice;                            //!< Elements of the coupling matrix. 
    
        std::vector<float> massFactorsHost;              //!< Mass factors (CPU)
        float* massFactorsDevice;                        //!< Mass factors (GPU)
        
};

#endif
