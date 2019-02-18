#include "gmxpre.h"
#include <assert.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>

#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/pbcutil/gpu_pbc.cuh"
#include "gromacs/pbcutil/pbc.h"


#include "gromacs/mdlib/lincs_cuda.h"

#if defined(_MSVC)
#include <limits>
#endif

#define BLOCK_SIZE 256


/*! \brief Main kernel for lincs constraints.
 *
 * \todo Combine arguments
 * \todo Move everythong to local/shared memory, try to get rid of atomics.
 */
__global__ void lincs_combined_kernel(int         ncons,
                                          const float3 *x,
                                          float3 *xp,
            
                                          int nOrder, int nIter,
                                          const gmx::AtomPair  *constraints,
                                          const real  *constraintsR0,

                                          const int  *coupledConstraintsCounts,
                                          const int  *coupledConstraintsIdxes,
                                          const real *massFactors,
                                          real       *matrixA,

                                          int width, //Make it into gridDim*blockSize??


                                          const PbcAiuc    pbcAiuc,

                                          real       *mlambda,
                                          const real *invmass,
                                          bool updateVelocities, float3* v, real invdt)
{

    int b = blockIdx.x*blockDim.x+threadIdx.x;
    int bs = threadIdx.x;
    int blockStarts = blockIdx.x*blockDim.x;

    extern __shared__ float3 r[];
    extern __shared__ float  rhs[];

    if (b < ncons)
    {

        float mvb;
        
        gmx::AtomPair pair = constraints[b];
        int i = pair.index1; 
        int j = pair.index2;
        
        real ml = 0.0;
        
        if(i != -1){
            float bllen = constraintsR0[b];
           
            float im1      = invmass[i];
            float im2      = invmass[j];
            
            float blc = rsqrt(im1 + im2);
            
            
            float3 xi = x[i];
            float3 xj = x[j];
            
            float3 dx = pbcDxAiucFloat3(pbcAiuc, xi, xj);
            float rlen = rsqrtf(dx.x*dx.x + dx.y*dx.y + dx.z*dx.z);
            
            float3 rb = rlen*dx;
            r[bs] = rb;

            xi = xp[i];
            xj = xp[j];
            dx = pbcDxAiucFloat3(pbcAiuc, xi, xj);
            
            mvb = blc*((rb.x*dx.x + rb.y*dx.y + rb.z*dx.z) - bllen);

            float sol  = mvb; 

            __syncthreads();
            
            int coupledConstraintsCount = coupledConstraintsCounts[b];

            for (int n = 0; n < coupledConstraintsCount; n++)
            {
                int index = n*width + b;
                int b1 = coupledConstraintsIdxes[index];  //Can be moved to local/shared memory
                
                float3 rb1 = r[b1-blockStarts];
                matrixA[index] = massFactors[index]*(rb.x*rb1.x + rb.y*rb1.y + rb.z*rb1.z);

            }
            
            rhs[bs] = mvb;  
            __syncthreads();
            
            for (int rec = 0; rec < nOrder; rec++)
            {
                mvb = 0;
                 
                for (int n = 0; n < coupledConstraintsCount; n++)
                {
                    int index = n*width + b;
                    int b1 = coupledConstraintsIdxes[index];
                    
                    mvb = mvb + matrixA[index]*rhs[b1-blockStarts + blockDim.x*(rec % 2)];

                }
                rhs[bs + blockDim.x*((rec + 1) % 2)] = mvb;
                sol  = sol + mvb;
                __syncthreads();
                
            }

            ml = blc*sol;

            mvb      = ml;

            float3 tmp     = rb*mvb;

            atomicAdd(&xp[i], -tmp*im1);
            atomicAdd(&xp[j], tmp*im2);

            __syncthreads();

            for (int iter = 0; iter < nIter; iter++)
            {
                float len2, dlen2;

                xi = xp[i];
                xj = xp[j];

               
                dx = pbcDxAiucFloat3(pbcAiuc, xi, xj);
               
                len2  = bllen*bllen;
                dlen2 = 2.0f*len2 - norm2(dx);

                if (dlen2 > 0)
                {
                    mvb = blc*(bllen - dlen2*rsqrt(dlen2));
                }
                else
                {
                    mvb = blc*bllen;
                }
                
                rhs[bs]  = mvb;
                sol  = mvb;
                __syncthreads();
                
                for (int rec = 0; rec < nOrder; rec++)
                {
                    mvb = 0;

                    for (int n = 0; n < coupledConstraintsCount; n++)
                    {
                        int index = n*width + b;
                        int b1 = coupledConstraintsIdxes[index];
                        
                        mvb = mvb + matrixA[index]*rhs[b1-blockStarts + blockDim.x*(rec % 2)];

                    }
                    rhs[bs + blockDim.x*((rec + 1) % 2)] = mvb;
                    sol  = sol + mvb;
                    __syncthreads();

                }

                mvb         = blc*sol;
                float blc_sol  = mvb;
                ml += mvb;

                mvb      = blc_sol;
                
                float3 tmp = rb*mvb;
                
                atomicAdd(&xp[i], -tmp*im1);
                atomicAdd(&xp[j], tmp*im2);
                __syncthreads();
            }
            
            if(updateVelocities)
            {
                mvb      = invdt*ml;
                
                float3 tmp     = rb*mvb;
                
                atomicAdd(&v[i], -tmp*im1);
                atomicAdd(&v[j], tmp*im2);
            }
        
        }
        mlambda[b] = ml; // Needed for virial
    }

    return;
}


__global__ void lincs_update_virial_kernel_new(int         ncons,
                                           const float3 *x,
                                           real      *virialScaled,
                                           const real *mlambda,
                                           const gmx::AtomPair  *constraints,
                                           const real *constraintsR0,
                                           const PbcAiuc    pbcAiuc)
{


    int b = blockIdx.x*blockDim.x+threadIdx.x;

    // Constraint virial
    if (b < ncons)
    {
        float tmp0, tmp1;
        
        gmx::AtomPair pair = constraints[b];
        int i = pair.index1; 
        int j = pair.index2;

        if(i != -1){
            
            float3 xi = x[i];
            float3 xj = x[j];
            
            float3 dx = pbcDxAiucFloat3(pbcAiuc, xi, xj);
            float rlen = rsqrtf(dx.x*dx.x + dx.y*dx.y + dx.z*dx.z);
                        
            rvec rb;

            rb[0] = rlen*dx.x;
            rb[1] = rlen*dx.y;
            rb[2] = rlen*dx.z;

            tmp0 = -constraintsR0[b]*mlambda[b];
            for (int d1 = 0; d1 < DIM; d1++)
            {
                tmp1 = tmp0*rb[d1];
                for (int d2 = 0; d2 < DIM; d2++)
                {
                    atomicAdd(&(virialScaled[d1*DIM+d2]), -tmp1*rb[d2]);
                }
            }
        }
    }
    return;
}


struct gpuUpdateConstraintsData{

    cudaStream_t     stream;
};

//TEMPORARY Global variable to store a copy of above struct within this module.
//TODO move this to a central location (e.g. gpu_nbv) and pass through fn args.
static gpuUpdateConstraintsData gpuUCDmod;

#define cudaCheckError() {                                          \
        cudaError_t e = cudaGetLastError();                                 \
        if (e != cudaSuccess) {                                              \
            printf("Cuda failure %s:%d: '%s'\n", __FILE__, __LINE__, cudaGetErrorString(e));           \
            exit(0); \
        }                                                                 \
}

/*! \brief Create LINCS object
 *
 * \param [in] nAtom  Number of atoms that will be handles by LINCS.
 *                    Used to compute the memory size in allocations and copy.
 * \param [in] nOrder Number of iterations used to compute inverse matrix.
 * \param [in] nOrder LINCS projection order for correcting the dirrection of comstraint.
 */
LincsCuda::LincsCuda(int nAtom,
                     int nIter,
                     int nOrder)
                     : nAtom(nAtom), nIter(nIter), nOrder(nOrder)

{
    GMX_ASSERT(sizeof(real) == sizeof(float), "Real numbers should be in single precision in GPU code.");
    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;
    cudaMalloc(&xDevice, nAtom*DIM*sizeof(float));    
    cudaMalloc(&xpDevice, nAtom*DIM*sizeof(float));    
    cudaMalloc(&vDevice, nAtom*DIM*sizeof(float));

    cudaMalloc(&virialScaledDevice, DIM*DIM*sizeof(float));
    maxConstraintsNumberSoFar = 0;
    cudaStreamCreate(&gpuUCD->stream);
    cudaCheckError();
}

LincsCuda::~LincsCuda()
{
}


/*! \brief Apply LINCS.
 *
 * Applies LINCS to coordinates and velocities, stored on GPU.
 * Data at pointers xPrime and v (class fields) change in the GPU
 * memory. The results are not automatically copied back to the CPU 
 * memory. Method uses this class data structures which should be 
 * updated when needed using update method. 
 *
 * \param[in] updateVelocities  If the velocities should be constrained.
 * \param[in] invdt             Inversed timestep (to scale lagrange 
 *                              multipliers when velocities are updated)
 * \param[in] bCalcVir          If virial should be updated.
 * \param[in] scaleLambda       If the Lagrange multipliers should be scaled
 *                              before virial is computed.
 * \param[in,out] virialScaled  Scaled virial tensor to be updated.
 */
void LincsCuda::apply(const bool       updateVelocities,
                      const real       invdt,
                      const gmx_bool   bCalcVir,
                      tensor           virialScaled)
{

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;
    cudaStream_t stream = gpuUCD->stream;

    cudaCheckError();

    int blockSize = BLOCK_SIZE;
    int blockCount = (nConstraintsThreads + blockSize - 1)/blockSize;
    int width = coupledConstraintsCountsHost.size();
    
 

    lincs_combined_kernel
            <<< blockCount, blockSize, blockSize*DIM*sizeof(float), stream>>>
            (nConstraintsThreads, xDevice, xpDevice,
            nOrder, nIter,
            constraintsDevice, constraintsR0Device,
            coupledConstraintsCountsDevice, coupledConstraintsIdxesDevice, 
            massFactorsDevice, matrixADevice,            
            width,
            pbcAiuc,
            mlambdaDevice,
            invmassDevice,
            updateVelocities, vDevice, invdt);

    cudaCheckError();


    if (bCalcVir)
    {
       cudaMemcpy(virialScaledDevice, virialScaled, DIM*DIM*sizeof(real), cudaMemcpyHostToDevice);

        lincs_update_virial_kernel_new
        <<< blockCount, blockSize, 0, stream>>>
        (nConstraintsThreads, xDevice, virialScaledDevice, mlambdaDevice, constraintsDevice, constraintsR0Device, pbcAiuc);

        cudaCheckError();

        cudaMemcpy(virialScaled, virialScaledDevice, DIM*DIM*sizeof(real), cudaMemcpyDeviceToHost);
        
        cudaCheckError();
    }

    cudaStreamSynchronize(stream);
    cudaCheckError();

    return;
}



/*! \brief 
 * Update PBC data.
 *
 * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
 *
 * \param[in] *pbc The PBC data in t_pbc format.
 */
void LincsCuda::setPBC(t_pbc *pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &pbcAiuc);
}


/*! \brief 
 * Update data-structures (e.g. after NB search step).
 *
 * Updates the constraints data and copies it to the GPU. Should be 
 * called if the particles were sorted, redistributed betwen domains, etc.
 * This version uses common data fromats so it can be called from anywhere 
 * in the code. Does not recycle the data preparation routines from the CPU 
 * version. Works only with simple case when all the constraints in idef are 
 * are handled by a single GPU. Triangles are not handled as special case. 
 *
 * Information about constraints is taken from:
 *     idef.il[F_CONSTR].iatoms  --- type (T) of constraint and two atom indeces (i1, i2)
 *     idef.iparams[T].constr.dA --- target length for constraint of type T
 * From t_mdatom, the code takes:
 *     md.invmass  --- array of inverse square root of masses for each atokm in the system.
 * 
 * \param[in] idef  Local topology data to get information on constraints from.
 * \param[in] md    Atoms data to get atom masses from.
 */
 void LincsCuda::set(const t_idef &idef,
                        const t_mdatoms      &md)
{

    int blockSize = BLOCK_SIZE;
    
    //t_idef idef = top.idef;
    t_iatom *iatoms = idef.il[F_CONSTR].iatoms;
    const int atomsCount = md.homenr;
    const int constraintsCount = idef.il[F_CONSTR].nr/3;    
    // Constructing adjacency list.
    std::vector<std::vector<std::tuple<int, int, int>>> atomsAdjacencyList(atomsCount);
    for(int c = 0; c < constraintsCount; c++)
    {
        //int type   = iatoms[3*c];
        int a1     = iatoms[3*c + 1];
        int a2     = iatoms[3*c + 2];
        
        // Each constraint will be represented as a tuple, containing index of the second constrained atom,
        // index of the constraint and a sign that indicates the order of atoms in which they are listed.
        // Sign is needed to compute the mass factors.
        atomsAdjacencyList.at(a1).push_back(std::make_tuple(a2, c, +1));
        atomsAdjacencyList.at(a2).push_back(std::make_tuple(a1, c, -1));
    }
    
    // Compute, how many coupled constraints are in front of each constraint.
    // Needed to introduce splits in data so that all coupled constraints will be computed in a single GPU block.
    std::vector<int> spaceNeeded;
    spaceNeeded.resize(constraintsCount, -1);
    int maxSpaceNeeded = 0;
    for (int c = 0; c < constraintsCount; c++)
    {
        int a1     = iatoms[3*c + 1];
        int a2     = iatoms[3*c + 2];
        
        // Constraint 'c' is counted twice, but it should be excluded altogether. Hence '-2'.
        spaceNeeded.at(c) = atomsAdjacencyList.at(a1).size() + atomsAdjacencyList.at(a2).size() - 2; 
        
        if(spaceNeeded.at(c) > maxSpaceNeeded)
        {
            maxSpaceNeeded = spaceNeeded.at(c);
        }
    }
    
    // Map of splits in the constraints data. For each 'old' constraint index gives 'new' which 
    // takes into account the empty spaces which might be needed in the end of each thread block.
    std::vector<int> splitMap;
    splitMap.resize(constraintsCount, -1);
    int currentMapIndex = 0;
    for (int c = 0; c < constraintsCount; c++)
    {
        if(currentMapIndex / blockSize != (currentMapIndex + spaceNeeded.at(c)) / blockSize){
            currentMapIndex = ((currentMapIndex/blockSize) + 1) * blockSize;
        }
        splitMap.at(c) = currentMapIndex;
        currentMapIndex++;
    }

    // Initialize constraints and their target indexes taking into account the splits in the
    // data arrays.
    gmx::AtomPair pair;
    pair.index1 = -1;
    pair.index2 = -1;
    constraintsHost.resize(currentMapIndex + blockSize - currentMapIndex % blockSize, pair);
    std::fill(constraintsHost.begin(), constraintsHost.end(), pair);
    constraintsR0Host.resize(currentMapIndex + blockSize - currentMapIndex % blockSize, 0.0);
    std::fill(constraintsR0Host.begin(), constraintsR0Host.end(), 0.0);
    for (int c = 0; c < constraintsCount; c++)
    {
        int a1     = iatoms[3*c + 1];
        int a2     = iatoms[3*c + 2];
        int type   = iatoms[3*c];
        
        gmx::AtomPair pair;
        pair.index1 = a1;
        pair.index2 = a2;
        constraintsHost.at(splitMap.at(c)) = pair;
        
        constraintsR0Host.at(splitMap.at(c)) = idef.iparams[type].constr.dA;

    }
    
    // The adjacency list of coupled constraints
    // We map a sigle thread to a single constraint, hence each thread 'c' will be using one element from 
    // coupledConstraintsCountsHost array, which is the number of constraints coupled to the constraint 'c'.
    // The coupled constraints indexes are placed into the coupledConstraintsIdxesHost array. Latter is organized 
    // as a one-dimentional array to ensure good memory alignment. It is adressed as [c + i*width], where i goes
    // from zero to the number of constraints coupled to 'c' minus 1.
    nConstraintsThreads = constraintsHost.size();
    coupledConstraintsCountsHost.resize(constraintsHost.size(), 0);
    coupledConstraintsIdxesHost.resize(maxSpaceNeeded*constraintsHost.size(), -1);
    massFactorsHost.resize(maxSpaceNeeded*constraintsHost.size(), -1);
   
    for (int c1 = 0; c1 < constraintsCount; c1++)
    {
        int c1a1     = iatoms[3*c1 + 1];
        int c1a2     = iatoms[3*c1 + 2];
        int c2;
        int c2a1;
        int c2a2;
        
        int sign;
        
        c2a1 = c1a1;
        for(unsigned j = 0; j < atomsAdjacencyList.at(c1a1).size(); j++)
        {
             
             std::tie(c2a2, c2, sign) = atomsAdjacencyList.at(c1a1).at(j);
             
             if(c1 != c2)
             {
                int index = constraintsHost.size()*coupledConstraintsCountsHost.at(splitMap.at(c1)) + splitMap.at(c1);
                
                coupledConstraintsIdxesHost.at(index) = c2;
                
                int center = c1a1;
                
                real blc1 = 1.0/sqrt(md.invmass[c1a1] + md.invmass[c1a2]);
                real blc2 = 1.0/sqrt(md.invmass[c2a1] + md.invmass[c2a2]);

                massFactorsHost.at(index) = -sign*md.invmass[center]*blc1*blc2;                
                
                coupledConstraintsCountsHost.at(splitMap.at(c1)) ++;
                
             }
        }
        
        c2a1 = c1a2;
        for(unsigned j = 0; j < atomsAdjacencyList.at(c1a2).size(); j++)
        {
             
             std::tie(c2a2, c2, sign) = atomsAdjacencyList.at(c1a2).at(j);
             
             if(c1 != c2)
             {
                int index = constraintsHost.size()*coupledConstraintsCountsHost.at(splitMap.at(c1)) + splitMap.at(c1);
                
                coupledConstraintsIdxesHost.at(index) = c2;
                
                int center = c1a2;
                
                real blc1 = 1.0/sqrt(md.invmass[c1a1] + md.invmass[c1a2]);
                real blc2 = 1.0/sqrt(md.invmass[c2a1] + md.invmass[c2a2]);

                massFactorsHost.at(index) = sign*md.invmass[center]*blc1*blc2;                
                
                coupledConstraintsCountsHost.at(splitMap.at(c1)) ++;
                
             }
        }
    }
    
    if(constraintsCount > maxConstraintsNumberSoFar)
    {
        
        if (maxConstraintsNumberSoFar > 0)
        {
            cudaFree(invmassDevice);
            cudaFree(mlambdaDevice);
            
            cudaFree(constraintsDevice);
            cudaFree(constraintsR0Device);
            
            cudaFree(coupledConstraintsCountsDevice);
            cudaFree(coupledConstraintsIdxesDevice);
            cudaFree(massFactorsDevice);
            cudaFree(matrixADevice);
            
        }
        maxConstraintsNumberSoFar = constraintsCount;
        
        cudaMalloc(&invmassDevice, nAtom*sizeof(real));
        cudaMalloc(&mlambdaDevice, constraintsHost.size()*sizeof(float));
        
        cudaMalloc(&constraintsDevice, constraintsHost.size()*sizeof(gmx::AtomPair));
        cudaMalloc(&constraintsR0Device, constraintsHost.size()*sizeof(float));
        
        cudaMalloc(&coupledConstraintsCountsDevice, coupledConstraintsCountsHost.size()*sizeof(int));
        cudaMalloc(&coupledConstraintsIdxesDevice, coupledConstraintsIdxesHost.size()*sizeof(int));
        cudaMalloc(&massFactorsDevice, massFactorsHost.size()*sizeof(float));
        cudaMalloc(&matrixADevice, coupledConstraintsIdxesHost.size()*sizeof(float));
                
    }
        
    cudaMemcpy(constraintsDevice, constraintsHost.data(), constraintsHost.size()*sizeof(gmx::AtomPair), cudaMemcpyHostToDevice);
    cudaMemcpy(constraintsR0Device, constraintsR0Host.data(), constraintsHost.size()*sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(coupledConstraintsCountsDevice, coupledConstraintsCountsHost.data(), 
                coupledConstraintsCountsHost.size()*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(coupledConstraintsIdxesDevice, coupledConstraintsIdxesHost.data(), 
                coupledConstraintsIdxesHost.size()*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(massFactorsDevice, massFactorsHost.data(), massFactorsHost.size()*sizeof(float), cudaMemcpyHostToDevice);

    cudaCheckError();
    
    GMX_ASSERT(md.invmass != nullptr, "Masses of attoms should be specified.\n");
    cudaMemcpy(invmassDevice, md.invmass, nAtom*sizeof(real), cudaMemcpyHostToDevice);

    cudaCheckError();
    
}


/*! \brief 
 * Update data-structures (e.g. after NB search step).
 *
 * Updates the constraints data and copies it to the GPU. Should be 
 * called after the particles were sorted, redistributed betwen domains, etc.
 * This is the version of the initialization that uses the data structures, prepared by 
 * CPU version of the algorithm. This is not efficient, since some of the structures are 
 * re-arranged twice. However, it can be usefull to sccess the performance of the rewritten 
 * preparation above.
 *
 * \param [in] atoms   Constraints in AtomPair format.
 * \param [in] blnr    Map of the coupled constraints (blnr[c] to blnr[c+1] are indexes in blbnb 
 *                     of constraints, coupled with constraint c.
 * \param [in] blbnb   List of coupled constraints (addressed using blnr).
 * \param [in] blmf    Mass-factors ( (+/-) * (1/sqrt(1/m1 + 1/m2)) * (1/m2) * 1/sqrt(1/m2 + 1/m3)),
 *                     where m1 and m3 are coupled through m2 and sign + or - indicates the order, 
 *                     in which they are arranged in atoms array.
 * \param [in] bllen   Target length for each constraint.
 * \param [in] invmass Inverse mass for each atom.
 * 
 */
/*void LincsCuda::set(std::vector<gmx::AtomPair> atoms,
                                std::vector<int>     blnr,
                                std::vector<int> blbnb,
                                std::vector<real>    blmf,
                                std::vector<real, gmx::AlignedAllocator < real>> bllen,
                                const real            *invmass)
{

    int blockSize = BLOCK_SIZE;
    
    int nConstraints = atoms.size();
    
    std::vector<int> atomsToUpdateGroups;
    std::vector<int> constraintsToUpdateGroups;
    
    atomsToUpdateGroups.resize(nAtom, -1);
    std::fill(atomsToUpdateGroups.begin(), atomsToUpdateGroups.end(), -1);
    constraintsToUpdateGroups.resize(nConstraints, -1);
    int currentGroup = 0;
    int updateGroupSet = -1;
    for (int c = 0; c < nConstraints; c++)
    {
        int i = atoms.at(c).index1;
        int j = atoms.at(c).index2;
        
        GMX_ASSERT((i < nAtom) && (j < nAtom), "Atom index is out of range.");
        if(atomsToUpdateGroups.at(i) != -1 && atomsToUpdateGroups.at(j) != -1)
        {
            if(atomsToUpdateGroups.at(i) == atomsToUpdateGroups.at(j))
            {
                updateGroupSet = atomsToUpdateGroups.at(i);
                // Found triangle ??
            } 
            else
            {
                // ERROR?
            }
        }
        else if(atomsToUpdateGroups.at(i) != -1 && atomsToUpdateGroups.at(j) == -1)
        {
            atomsToUpdateGroups.at(j) = atomsToUpdateGroups.at(i);
            updateGroupSet = atomsToUpdateGroups.at(i);
        }
        else if(atomsToUpdateGroups.at(j) != -1 && atomsToUpdateGroups.at(i) == -1)
        {
            atomsToUpdateGroups.at(i) = atomsToUpdateGroups.at(j);
            updateGroupSet = atomsToUpdateGroups.at(j);
        }
        else if(atomsToUpdateGroups.at(i) == -1 && atomsToUpdateGroups.at(j) == -1)
        {
            atomsToUpdateGroups.at(j) = currentGroup;
            atomsToUpdateGroups.at(i) = currentGroup;
            updateGroupSet = currentGroup;
            currentGroup ++;
        }
        else
        {
            updateGroupSet = -1;
        }
        constraintsToUpdateGroups.at(c) = updateGroupSet;
    }
    int nUpdateGroup = currentGroup;
        
    std::vector<int> groupSizes;
    groupSizes.resize(nUpdateGroup, 0);
    for (int c = 0; c < nConstraints; c++)
    {
        groupSizes.at(constraintsToUpdateGroups.at(c))++;
    }
    
    std::vector<int> spaceNeeded;
    spaceNeeded.resize(nConstraints, -1);
    int maxSpaceNeeded = 0;
    for (int c = 0; c < nConstraints; c++)
    {
        spaceNeeded.at(c) = 0;
        int c1 = 1;
        while(c + c1 < nConstraints && constraintsToUpdateGroups.at(c) == constraintsToUpdateGroups.at(c+c1))
        {
            spaceNeeded.at(c) ++;
            c1++;
            GMX_ASSERT(c1 < blockSize, "Group is larger then the block size.\n");
        }
        if(spaceNeeded.at(c) > maxSpaceNeeded)
        {
            maxSpaceNeeded = spaceNeeded.at(c);
        }
    }    
    
    std::vector<int> splitMap;
    splitMap.resize(nConstraints, -1);
    int currentMapIndex = 0;
    for (int c = 0; c < nConstraints; c++)
    {
        if(currentMapIndex / blockSize != (currentMapIndex + spaceNeeded.at(c)) / blockSize){
            //currentMapIndex += blockSize - currentMapIndex % blockSize;
            currentMapIndex = ((currentMapIndex/blockSize) + 1) * blockSize;
        }
        splitMap.at(c) = currentMapIndex;
        currentMapIndex++;
    }
    
    
    gmx::AtomPair pair;
    pair.index1 = -1;
    pair.index2 = -1;
    constraintsHost.resize(currentMapIndex + blockSize - currentMapIndex % blockSize, pair);
    std::fill(constraintsHost.begin(), constraintsHost.end(), pair);
    constraintsR0Host.resize(currentMapIndex + blockSize - currentMapIndex % blockSize, 0.0);
    std::fill(constraintsR0Host.begin(), constraintsR0Host.end(), 0.0);
    for (int c = 0; c < nConstraints; c++)
    {
        pair = atoms.at(c);
        constraintsHost.at(splitMap.at(c)) = pair;
        constraintsR0Host.at(splitMap.at(c)) = bllen.at(c);
    }
    
    nConstraintsThreads = constraintsHost.size();
    coupledConstraintsCountsHost.resize(constraintsHost.size(), 0);
    coupledConstraintsIdxesHost.resize(maxSpaceNeeded*constraintsHost.size(), -1);
    massFactorsHost.resize(maxSpaceNeeded*constraintsHost.size(), -1);
    for (int c = 0; c < nConstraints; c++)
    {
        coupledConstraintsCountsHost.at(splitMap.at(c)) = blnr.at(c+1)-blnr.at(c);
        for (int n = blnr.at(c); n < blnr.at(c+1); n++)
        {
            int index = (n - blnr.at(c))*constraintsHost.size() + splitMap.at(c);
            coupledConstraintsIdxesHost.at(index) = splitMap.at(blbnb.at(n));
            massFactorsHost.at(index) = blmf.at(n);
        }
    }
    
    if(nConstraints > maxConstraintsNumberSoFar)
    {
        
        if (maxConstraintsNumberSoFar > 0)
        {
            cudaFree(invmassDevice);
            cudaFree(mlambdaDevice);
            
            cudaFree(constraintsDevice);
            cudaFree(constraintsR0Device);
            
            cudaFree(coupledConstraintsCountsDevice);
            cudaFree(coupledConstraintsIdxesDevice);
            cudaFree(massFactorsDevice);
            cudaFree(matrixADevice);
            
        }
        maxConstraintsNumberSoFar = nConstraints;
        
        cudaMalloc(&invmassDevice, nAtom*sizeof(real));
        cudaMalloc(&mlambdaDevice, constraintsHost.size()*sizeof(float));
        
        cudaMalloc(&constraintsDevice, constraintsHost.size()*sizeof(gmx::AtomPair));
        cudaMalloc(&constraintsR0Device, constraintsHost.size()*sizeof(float));
        
        cudaMalloc(&coupledConstraintsCountsDevice, coupledConstraintsCountsHost.size()*sizeof(int));
        cudaMalloc(&coupledConstraintsIdxesDevice, coupledConstraintsIdxesHost.size()*sizeof(int));
        cudaMalloc(&massFactorsDevice, massFactorsHost.size()*sizeof(float));
        cudaMalloc(&matrixADevice, coupledConstraintsIdxesHost.size()*sizeof(float));
                
    }
        
    cudaMemcpy(constraintsDevice, constraintsHost.data(), constraintsHost.size()*sizeof(gmx::AtomPair), cudaMemcpyHostToDevice);
    cudaMemcpy(constraintsR0Device, constraintsR0Host.data(), constraintsHost.size()*sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(coupledConstraintsCountsDevice, coupledConstraintsCountsHost.data(), 
                coupledConstraintsCountsHost.size()*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(coupledConstraintsIdxesDevice, coupledConstraintsIdxesHost.data(), 
                coupledConstraintsIdxesHost.size()*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(massFactorsDevice, massFactorsHost.data(), massFactorsHost.size()*sizeof(float), cudaMemcpyHostToDevice);

    cudaCheckError();
    
    if (invmass != nullptr)
    {
        cudaMemcpy(invmassDevice, invmass, nAtom*sizeof(real), cudaMemcpyHostToDevice);
    }

    cudaCheckError();

}*/

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
    cudaMemcpy(xDevice, x, nAtom*sizeof(float3), cudaMemcpyHostToDevice);
    cudaCheckError();
    cudaMemcpy(xpDevice, xp, nAtom*sizeof(float3), cudaMemcpyHostToDevice);
    cudaCheckError();
    if(v != nullptr)
    {
        cudaMemcpy(vDevice, v, nAtom*sizeof(float3), cudaMemcpyHostToDevice);
    }
    
    cudaCheckError();
    cudaDeviceSynchronize();
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
    cudaMemcpy(xp, xpDevice, nAtom*sizeof(float3), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaCheckError();
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
    cudaMemcpy(v, vDevice, nAtom*sizeof(float3), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaCheckError();
}

/*! \brief 
 * Set the internal GPU-memory x, xprime and v pointers.
 *
 * Data is not copied. The data are assumed to be in float3/fvec format
 * (float3 is used internaly, but the data layout should be identical).
 *
 * \param[in] *xDevice  Pointer to the coordinates before integrator update (on GPU)   
 * \param[in] *xpDevice Pointer to the coordinates after integrator update, before update (on GPU)   
 * \param[in] *vDevice  Pointer to the velocities before integrator update (on GPU)   
 */
void LincsCuda::setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice)
{
    this->xDevice = (float3*)xDevice;
    this->xpDevice = (float3*)xpDevice;
    this->vDevice = (float3*)vDevice;
}
