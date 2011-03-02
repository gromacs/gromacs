#ifndef _NB_KERNEL_UTILS_CUH_
#define _NB_KERNEL_UTILS_CUH_

#define CELL_SIZE_2         (CELL_SIZE * CELL_SIZE)
#define STRIDE_DIM          (CELL_SIZE_2)
#define STRIDE_SI           (3*STRIDE_DIM)
#define GPU_FACEL           (138.935485f)

/* texture reference bound to the cudata.nbfp array */
texture<float, 1, cudaReadModeElementType> tex_nbfp;

/* texture reference bound to the cudata.coulomb_tab array */
texture<float, 1, cudaReadModeElementType> tex_coulomb_tab;

/* Original idea: OpenMM */
static inline __device__ float interpolate_coulomb_force_r(float r, float scale)
{  
    float   normalized = scale * r;
    int     index = (int) normalized;
    float   fract2 = normalized - index;
    float   fract1 = 1.0f - fract2;

    return  fract1 * tex1Dfetch(tex_coulomb_tab, index) 
            + fract2 * tex1Dfetch(tex_coulomb_tab, index + 1);
}

static inline __device__ 
void reduce_force_j_generic(float *fbuf, float4 *fout, 
                            int tidxi, int tidxj, int aidx)
{
    if (tidxi == 0)
    {
        float4 f = make_float4(0.0f);
        for (int j = tidxj * CELL_SIZE; j < (tidxj + 1) * CELL_SIZE; j++)
        {
            f.x += fbuf[                 j];
            f.y += fbuf[    STRIDE_DIM + j];
            f.z += fbuf[2 * STRIDE_DIM + j];
        }

        atomicAdd(&fout[aidx].x, f.x);
        atomicAdd(&fout[aidx].y, f.y);
        atomicAdd(&fout[aidx].z, f.z);
    }
}

static inline __device__ 
void reduce_force_i_generic(float *fbuf, float4 *fout, 
                            float3 *fbuf_shift, gmx_bool calc_fshift, 
                            int tidxi, int tidxj, int aidx)
{
    if (tidxj == 0)
    {
        float4 f = make_float4(0.0f);
        for (int j = tidxi; j < CELL_SIZE_2; j += CELL_SIZE)
        {
            f.x += fbuf[                 j];
            f.y += fbuf[    STRIDE_DIM + j];
            f.z += fbuf[2 * STRIDE_DIM + j];
        }

        atomicAdd(&fout[aidx].x, f.x);
        atomicAdd(&fout[aidx].y, f.y);
        atomicAdd(&fout[aidx].z, f.z);

        if (calc_fshift)
        {
            fbuf_shift->x += f.x;
            fbuf_shift->y += f.y;
            fbuf_shift->z += f.z;
        }
    }
}

static inline __device__ 
void reduce_force_i_pow2(volatile float *fbuf, float4 *fout, 
                         float3 *fbuf_shift, gmx_bool calc_fshift, 
                         int tidxi, int tidxj, int aidx)
{
    int     i, j; 
    float4  f = make_float4(0.0f);

    /* Reduce the initial CELL_SIZE values for each i atom to half
       every step by using CELL_SIZE * i threads. */
    i = CELL_SIZE/2;
    # pragma unroll 5
    for (j = CELL_SIZE_POW2_EXPONENT - 1; j > 0; j--)
    {
        if (tidxj < i)
        {

            fbuf[                 tidxj * CELL_SIZE + tidxi] += fbuf[                 (tidxj + i) * CELL_SIZE + tidxi];
            fbuf[    STRIDE_DIM + tidxj * CELL_SIZE + tidxi] += fbuf[    STRIDE_DIM + (tidxj + i) * CELL_SIZE + tidxi];
            fbuf[2 * STRIDE_DIM + tidxj * CELL_SIZE + tidxi] += fbuf[2 * STRIDE_DIM + (tidxj + i) * CELL_SIZE + tidxi];
        }
        i >>= 1;
    }

    /* i == 1, last reduction step, writing to global mem */
    if (tidxj == 0)
    {
        f.x = fbuf[                 tidxj * CELL_SIZE + tidxi] + fbuf[                 (tidxj + i) * CELL_SIZE + tidxi];
        f.y = fbuf[    STRIDE_DIM + tidxj * CELL_SIZE + tidxi] + fbuf[    STRIDE_DIM + (tidxj + i) * CELL_SIZE + tidxi]; 
        f.z = fbuf[2 * STRIDE_DIM + tidxj * CELL_SIZE + tidxi] + fbuf[2 * STRIDE_DIM + (tidxj + i) * CELL_SIZE + tidxi];

        atomicAdd(&fout[aidx].x, f.x);
        atomicAdd(&fout[aidx].y, f.y);
        atomicAdd(&fout[aidx].z, f.z);

        if (calc_fshift)
        {
            fbuf_shift->x += f.x;
            fbuf_shift->y += f.y;
            fbuf_shift->z += f.z;
        }
    }
}

static inline __device__ 
void reduce_force_i(float *fbuf, float4 *f,
                    float3 *fbuf_shift, gmx_bool calc_fshift, 
                    int tidxi, int tidxj, int ai)
{
    if ((CELL_SIZE & (CELL_SIZE - 1)))
    {
        reduce_force_i_generic(fbuf, f, fbuf_shift, calc_fshift, tidxi, tidxj, ai);
    }
    else
    {
        reduce_force_i_pow2(fbuf, f, fbuf_shift, calc_fshift, tidxi, tidxj, ai);
    }
}

static inline __device__ void reduce_energy_pow2(volatile float *buf,
                                          float *e_lj, float *e_el,
                                          unsigned int tidx)
{
    int     i, j; 
    float   e1, e2;

    i = CELL_SIZE_2/2;

# pragma unroll 10
    for (j = 2 * CELL_SIZE_POW2_EXPONENT - 1; j > 0; j--)
    {
        if (tidx < i)
        {
            buf[             tidx] += buf[             tidx + i];
            buf[STRIDE_DIM + tidx] += buf[STRIDE_DIM + tidx + i];
        }
        i >>= 1;
    }

    /* last reduction step, writing to global mem */
    if (tidx == 0)
    {
        e1 = buf[             tidx] + buf[             tidx + i];
        e2 = buf[STRIDE_DIM + tidx] + buf[STRIDE_DIM + tidx + i];

        atomicAdd(e_lj, e1);
        atomicAdd(e_el, e2); 
    }
}

/*********************************************************************************/
/* Old stuff  */
#if 0
inline __device__ float coulomb(float q1, 
                                float q2,
                                float r2, 
                                float inv_r, 
                                float inv_r2, 
                                float beta,
                                float erfc_tab_scale)
{
    float x      = r2 * inv_r * beta;
    float x2     = x * x; 
    // float inv_x2 = inv_r2 / (beta * beta); 
    float res    =
        q1 * q2 * (erfc(x) * inv_r + beta * exp(-x2)) * inv_r2;
    return res;
}
#endif 

#endif /*_NB_KERNEL_UTILS_CUH_*/ 
