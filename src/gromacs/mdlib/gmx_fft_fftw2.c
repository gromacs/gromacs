/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 *
 * Gromacs 4.0                         Copyright (c) 1991-2003
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_FFT_FFTW2

#include <errno.h>
#include <string.h>
#include <stdlib.h>


#include "gmx_fft.h"
#include "gmx_fatal.h"



#ifdef FFTW2_NAME_FFTW
#  include<fftw.h>
#  include<rfftw.h>
#elif defined FFTW2_NAME_SFFTW
#  include<sfftw.h>
#  include<srfftw.h>
#elif defined FFTW2_NAME_DFFTW
#  include<dfftw.h>
#  include<drfftw.h>
#else
#error No FFTW2 name defined - you must define one of:
#error FFTW2_NAME_FFTW, FFTW2_NAME_SFFTW, FFTW2_NAME_DFFTW
#endif


/* Contents of the FFTW2 setup */
struct gmx_fft 
{
    int               ndim;         /**< Number of dimensions in transform.   */
    int               nx;           /**< Data X dimension                     */                          
    int               ny;           /**< Data Y dimension                     */
    int               nz;           /**< Data Z dimension                     */
    /* Arrays with fftw2 plans. 
     * First index is 0 for out-of-place, 1 for in-place transform.
     * Second index is 0 for backward, 1 for forward.
     */
    fftw_plan         single[2][2]; /**< Plans for 1d transforms.             */
    fftwnd_plan       multi[2][2];  /**< Plans for n-d transforms.            */
    real *            work;         /**< Avoid overwriting input for c2r ffts */
};



int
gmx_fft_init_1d(gmx_fft_t *        pfft,
                int                nx,
                gmx_fft_flag       flags) 
{
    int i,j;
    gmx_fft_t           fft;
    int                    fftw_flags;

    /* FFTW2 is slow to measure, so we do not use it */
    
    /* If you change this, add an #ifndef for GMX_DISABLE_FFTW_MEASURE around it! */
    fftw_flags = FFTW_ESTIMATE;    

    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if( (fft = (gmx_fft_t)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    

    
    fft->single[0][0] = fftw_create_plan(nx,FFTW_BACKWARD,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->single[0][1] = fftw_create_plan(nx,FFTW_FORWARD,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->single[1][0] = fftw_create_plan(nx,FFTW_BACKWARD,FFTW_IN_PLACE|fftw_flags);
    fft->single[1][1] = fftw_create_plan(nx,FFTW_FORWARD,FFTW_IN_PLACE|fftw_flags);

    
    fft->multi[0][0] = NULL;
    fft->multi[0][1] = NULL;
    fft->multi[1][0] = NULL;
    fft->multi[1][1] = NULL;

    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            if(fft->single[i][j] == NULL)
            {
                gmx_fatal(FARGS,"Error initializing FFTW2 plan.");
                gmx_fft_destroy(fft);
                return -1;
            }        
        }
    }
    
    /* No workspace needed for complex-to-complex FFTs */
    fft->work = NULL;
    
    fft->ndim = 1;
    fft->nx   = nx;
    
    *pfft = fft;
    return 0;
}



int
gmx_fft_init_1d_real(gmx_fft_t *        pfft,
                     int                nx,
                     gmx_fft_flag       flags) 
{
    int i,j;
    gmx_fft_t             fft;
    int                    fftw_flags;
    
    /* FFTW2 is slow to measure, so we do not use it */
    /* If you change this, add an #ifndef for GMX_DISABLE_FFTW_MEASURE around it! */
    fftw_flags = FFTW_ESTIMATE;    
    
    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if( (fft = (gmx_fft_t)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
 
 
    fft->single[0][0] = rfftw_create_plan(nx,FFTW_COMPLEX_TO_REAL,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->single[0][1] = rfftw_create_plan(nx,FFTW_REAL_TO_COMPLEX,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->single[1][0] = rfftw_create_plan(nx,FFTW_COMPLEX_TO_REAL,FFTW_IN_PLACE|fftw_flags);
    fft->single[1][1] = rfftw_create_plan(nx,FFTW_REAL_TO_COMPLEX,FFTW_IN_PLACE|fftw_flags);

    
    fft->multi[0][0] = NULL;
    fft->multi[0][1] = NULL;
    fft->multi[1][0] = NULL;
    fft->multi[1][1] = NULL;
    
    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            if(fft->single[i][j] == NULL)
            {
                gmx_fatal(FARGS,"Error initializing FFTW2 plan.");
                gmx_fft_destroy(fft);
                return -1;
            }        
        }
    }
    
    /* FFTW2 overwrites the input when doing out-of-place complex-to-real FFTs.
     * This is not acceptable for the Gromacs interface, so we define a
     * work array and copy the data there before doing complex-to-real FFTs.
     */
    fft->work = (real *)malloc(sizeof(real)*( (nx/2 + 1)*2) );
    if(fft->work == NULL)
    {
        gmx_fatal(FARGS,"Cannot allocate complex-to-real FFT workspace.");
        gmx_fft_destroy(fft);
        return ENOMEM;
    }
    
    fft->ndim = 1;
    fft->nx   = nx;
  
    *pfft = fft;
    return 0;
}


	    
int
gmx_fft_init_2d(gmx_fft_t *        pfft,
                int                nx, 
                int                ny,
                gmx_fft_flag       flags) 
{
    int i,j;
    gmx_fft_t             fft;
    int                    fftw_flags;

    
    /* FFTW2 is slow to measure, so we do not use it */
    /* If you change this, add an #ifndef for GMX_DISABLE_FFTW_MEASURE around it! */
    fftw_flags = FFTW_ESTIMATE;    

    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if( (fft = (gmx_fft_t)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
 
    fft->single[0][0] = NULL;
    fft->single[0][1] = NULL;
    fft->single[1][0] = NULL;
    fft->single[1][1] = NULL;
        
    fft->multi[0][0] = fftw2d_create_plan(nx,ny,FFTW_BACKWARD,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->multi[0][1] = fftw2d_create_plan(nx,ny,FFTW_FORWARD,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->multi[1][0] = fftw2d_create_plan(nx,ny,FFTW_BACKWARD,FFTW_IN_PLACE|fftw_flags);
    fft->multi[1][1] = fftw2d_create_plan(nx,ny,FFTW_FORWARD,FFTW_IN_PLACE|fftw_flags);
    
    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            if(fft->multi[i][j] == NULL)
            {
                gmx_fatal(FARGS,"Error initializing FFTW2 plan.");
                gmx_fft_destroy(fft);
                return -1;
            }        
        }
    }

    /* No workspace needed for complex-to-complex FFTs */
    fft->work = NULL;
    
    fft->ndim = 2;
    fft->nx   = nx;
    fft->ny   = ny;

    *pfft = fft;
    return 0;
}




int
gmx_fft_init_2d_real(gmx_fft_t *        pfft,
                     int                nx, 
                     int                ny,
                     gmx_fft_flag       flags) 
{
    int i,j;
    gmx_fft_t             fft;
    int                    fftw_flags;

    
    /* FFTW2 is slow to measure, so we do not use it */
    /* If you change this, add an #ifndef for GMX_DISABLE_FFTW_MEASURE around it! */
    fftw_flags = FFTW_ESTIMATE;    

    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if( (fft = (gmx_fft_t)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    fft->single[0][0] = NULL;
    fft->single[0][1] = NULL;
    fft->single[1][0] = NULL;
    fft->single[1][1] = NULL;
    
    
    fft->multi[0][0] = rfftw2d_create_plan(nx,ny,FFTW_COMPLEX_TO_REAL,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->multi[0][1] = rfftw2d_create_plan(nx,ny,FFTW_REAL_TO_COMPLEX,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->multi[1][0] = rfftw2d_create_plan(nx,ny,FFTW_COMPLEX_TO_REAL,FFTW_IN_PLACE|fftw_flags);
    fft->multi[1][1] = rfftw2d_create_plan(nx,ny,FFTW_REAL_TO_COMPLEX,FFTW_IN_PLACE|fftw_flags);
    

    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            if(fft->multi[i][j] == NULL)
            {
                gmx_fatal(FARGS,"Error initializing FFTW2 plan.");
                gmx_fft_destroy(fft);
                return -1;
            }        
        }
    }
        
    /* FFTW2 overwrites the input when doing out-of-place complex-to-real FFTs.
     * This is not acceptable for the Gromacs interface, so we define a
     * work array and copy the data there before doing complex-to-real FFTs.
     */
    fft->work = (real *)malloc(sizeof(real)*( nx*(ny/2 + 1)*2) );
    if(fft->work == NULL)
    {
        gmx_fatal(FARGS,"Cannot allocate complex-to-real FFT workspace.");
        gmx_fft_destroy(fft);
        return ENOMEM;
    }
    

    fft->ndim = 2;
    fft->nx   = nx;
    fft->ny   = ny;

    *pfft = fft;
    return 0;
}


int
gmx_fft_init_3d(gmx_fft_t *        pfft,
                int                nx, 
                int                ny,
                int                nz,
                gmx_fft_flag       flags) 
{
    int i,j;
    gmx_fft_t             fft;
    int                    fftw_flags;

    
    /* FFTW2 is slow to measure, so we do not use it */
    /* If you change this, add an #ifndef for GMX_DISABLE_FFTW_MEASURE around it! */
    fftw_flags = FFTW_ESTIMATE;    

    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if( (fft = (gmx_fft_t)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    fft->single[0][0] = NULL;
    fft->single[0][1] = NULL;
    fft->single[1][0] = NULL;
    fft->single[1][1] = NULL;
    
    
    fft->multi[0][0] = fftw3d_create_plan(nx,ny,nz,FFTW_BACKWARD,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->multi[0][1] = fftw3d_create_plan(nx,ny,nz,FFTW_FORWARD,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->multi[1][0] = fftw3d_create_plan(nx,ny,nz,FFTW_BACKWARD,FFTW_IN_PLACE|fftw_flags);
    fft->multi[1][1] = fftw3d_create_plan(nx,ny,nz,FFTW_FORWARD,FFTW_IN_PLACE|fftw_flags);
    
    
    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            if(fft->multi[i][j] == NULL)
            {
                gmx_fatal(FARGS,"Error initializing FFTW2 plan.");
                gmx_fft_destroy(fft);
                return -1;
            }        
        }
    }
    
    /* No workspace needed for complex-to-complex FFTs */
    fft->work = NULL;
    fft->nx   = nx;
    fft->ny   = ny;
    fft->nz   = nz;
    
    fft->ndim = 3;
    
    *pfft = fft;
    return 0;
} 




int
gmx_fft_init_3d_real(gmx_fft_t *        pfft,
                     int                nx, 
                     int                ny,
                     int                nz,
                     gmx_fft_flag       flags) 
{
    int i,j;
    gmx_fft_t            fft;
    int                    fftw_flags;

    
    /* FFTW2 is slow to measure, so we do not use it */
    /* If you change this, add an #ifndef for GMX_DISABLE_FFTW_MEASURE around it! */
    fftw_flags = FFTW_ESTIMATE;    

    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if( (fft = (gmx_fft_t)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    fft->single[0][0] = NULL;
    fft->single[0][1] = NULL;
    fft->single[1][0] = NULL;
    fft->single[1][1] = NULL;
    
    
    fft->multi[0][0] = rfftw3d_create_plan(nx,ny,nz,FFTW_COMPLEX_TO_REAL,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->multi[0][1] = rfftw3d_create_plan(nx,ny,nz,FFTW_REAL_TO_COMPLEX,FFTW_OUT_OF_PLACE|fftw_flags);
    fft->multi[1][0] = rfftw3d_create_plan(nx,ny,nz,FFTW_COMPLEX_TO_REAL,FFTW_IN_PLACE|fftw_flags);
    fft->multi[1][1] = rfftw3d_create_plan(nx,ny,nz,FFTW_REAL_TO_COMPLEX,FFTW_IN_PLACE|fftw_flags);
    
    
    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            if(fft->multi[i][j] == NULL)
            {
                gmx_fatal(FARGS,"Error initializing FFTW2 plan.");
                gmx_fft_destroy(fft);
                return -1;
            }        
        }
    }
    
    /* FFTW2 overwrites the input when doing out-of-place complex-to-real FFTs.
     * This is not acceptable for the Gromacs interface, so we define a
     * work array and copy the data there before doing complex-to-real FFTs.
     */
    fft->work = (real *)malloc(sizeof(real)*( nx*ny*(nz/2 + 1)*2) );
    if(fft->work == NULL)
    {
        gmx_fatal(FARGS,"Cannot allocate complex-to-real FFT workspace.");
        gmx_fft_destroy(fft);
        return ENOMEM;
    }
    
    fft->ndim = 3;
    fft->nx   = nx;
    fft->ny   = ny;
    fft->nz   = nz;
    
    *pfft = fft;
    return 0;
} 


int 
gmx_fft_1d(gmx_fft_t                  fft,
           enum gmx_fft_direction     dir,
           void *                     in_data,
           void *                     out_data)
{
    int inplace   = (in_data == out_data);
    int isforward = (dir == GMX_FFT_FORWARD);

    if((fft->ndim != 1) ||
       ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)))
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }
    
    fftw_one(fft->single[inplace][isforward],(fftw_complex *)in_data,(fftw_complex *)out_data);
    
  return 0;
}

int 
gmx_fft_1d_real(gmx_fft_t                  fft,
                enum gmx_fft_direction     dir,
                void *                     in_data,
                void *                     out_data)
{
    /* FFTW2 1-dimensional real transforms are special.
     *
     * First, the complex data is stored in a special packed half-complex
     * fashion. To enable a standard common Gromacs interface this forces us
     * to always use out-of-place FFTs, and permute the data after 
     * real-to-complex FFTs or before complex-to-real FFTs.
     *
     * The input is also destroyed for out-of-place complex-to-real FFTs, but
     * this doesn't matter since we need to permute and copy the data into 
     * the work array first anyway.
     */
    real *     work = fft->work;
    t_complex *  data;
    int              n    = fft->nx;
    int              i;
    
    if((fft->ndim != 1) ||
       ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)))
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }
    
    if(dir==GMX_FFT_REAL_TO_COMPLEX) 
    {
        rfftw_one(fft->single[0][1],(fftw_real *)in_data,(fftw_real *)work);
        /* permute it back into data, in standard complex format 
         * instead of halfcomplex...
         */
        data = (t_complex *)out_data;
        
        data[0].re = work[0];
        data[0].im = 0;
        
        for(i=1;i<n/2;i++)
        {
            data[i].re = work[i];
            data[i].im = work[n-i];
        }

        data[i].re=work[i];
        
        if(2*i==n) 
        {
            data[i].im=0;
        }
        else
        {
            data[i].im=work[n-i];
        }
    }
    else
    {
        /* Complex-to-real. First permute standard format into halfcomplex */
        data = (t_complex *)in_data;
        
        work[0]=data[0].re;
        
        for(i=1;i<n/2;i++) 
        {
            work[i]  =data[i].re;
            work[n-i]=data[i].im;
        }      
        
        if(2*i!=n)
        {
            work[n-i]=data[i].im;
        }
        
        rfftw_one(fft->single[0][0],(fftw_real *)work,(fftw_real *)out_data);
    }

    return 0;
}


int 
gmx_fft_2d(gmx_fft_t                  fft,
           enum gmx_fft_direction     dir,
           void *                     in_data,
           void *                     out_data)
{
    int inplace   = (in_data == out_data);
    int isforward = (dir == GMX_FFT_FORWARD);
    
    if((fft->ndim != 2) ||
       ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)))
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }
    
    fftwnd_one(fft->multi[inplace][isforward],(fftw_complex *)in_data,(fftw_complex *)out_data);
    
    return 0;
}


int 
gmx_fft_2d_real(gmx_fft_t                  fft,
                enum gmx_fft_direction     dir,
                void *                     in_data,
                void *                     out_data)
{
    int inplace   = (in_data == out_data);
    int isforward = (dir == GMX_FFT_REAL_TO_COMPLEX);
    int sz;

    if((fft->ndim != 2) ||
       ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)))
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    if(inplace == 0)
    {
        /* Copy data to avoid overwriting input, and redirect input ptr to work array */
        sz = fft->nx*(fft->ny/2 + 1)*2;
        memcpy(fft->work,in_data,sz*sizeof(real));
        in_data = fft->work;
    }
    
    if(isforward)
    {
        rfftwnd_one_real_to_complex(fft->multi[inplace][isforward],(fftw_real *)in_data,(fftw_complex *)out_data);
    }
    else
    {
        rfftwnd_one_complex_to_real(fft->multi[inplace][isforward],(fftw_complex *)in_data,(fftw_real *)out_data);
    }
    
    return 0;
}


int 
gmx_fft_3d(gmx_fft_t                  fft,
           enum gmx_fft_direction     dir,
           void *                     in_data,
           void *                     out_data)
{
    int inplace   = (in_data == out_data);
    int isforward = (dir == GMX_FFT_FORWARD);
    
    if((fft->ndim != 3) ||
       ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)))
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }
    
    fftwnd_one(fft->multi[inplace][isforward],(fftw_complex *)in_data,(fftw_complex *)out_data);
    
    return 0;
}


int 
gmx_fft_3d_real(gmx_fft_t                  fft,
                enum gmx_fft_direction     dir,
                void *                     in_data,
                void *                     out_data)
{
    int inplace   = (in_data == out_data);
    int isforward = (dir == GMX_FFT_REAL_TO_COMPLEX);
    int sz;
    
    if((fft->ndim != 3) ||
       ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)))
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }
    
    if(inplace == 0)
    {
        /* Copy data to avoid overwriting input, and redirect input ptr to work array */
        sz = fft->nx*fft->ny*(fft->nz/2 + 1)*2;
        memcpy(fft->work,in_data,sz*sizeof(real));
        in_data = fft->work;
    }    

    if(isforward)
    {
        rfftwnd_one_real_to_complex(fft->multi[inplace][isforward],(fftw_real *)in_data,(fftw_complex *)out_data);
    }
    else
    {
        rfftwnd_one_complex_to_real(fft->multi[inplace][isforward],(fftw_complex *)in_data,(fftw_real *)out_data);
    }
    
    return 0;
}




void
gmx_fft_destroy(gmx_fft_t    fft)
{
    int i,j;
    
    if(fft != NULL)
    {
        for(i=0;i<2;i++)
        {
            for(j=0;j<2;j++)
            {
                if(fft->single[i][j] != NULL)
                {
                    rfftw_destroy_plan(fft->single[i][j]);
                    fft->single[i][j] = NULL;
                }
                if(fft->multi[i][j] != NULL)
                {
                    rfftwnd_destroy_plan(fft->multi[i][j]);
                    fft->multi[i][j] = NULL;
                }
            }
        }
        free(fft);
    }
}


#else
int
gmx_fft_fftw2_empty;
#endif /* GMX_FFT_FFTW2 */
