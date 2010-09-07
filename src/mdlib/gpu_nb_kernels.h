
inline __host__ __device__ float3 make_float3(float x)
{
    return make_float3(x, x, x);
}
__global__ void k_calc_nb(float3 *f, float3 *x, int natoms)
{
    unsigned int    tnb = blockDim.x * gridDim.x;
    unsigned int    tid = blockIdx.x * blockDim.x + threadIdx.x;

    /* accessing float3-s is not really a good idea, but ATM that's what we'll do*/
    for (int i = tid; i < natoms; i += tnb)
    {
        f[i] = make_float3(0.0);
    }
}
