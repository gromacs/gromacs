/* ************************************************************************
 * Copyright 2013 Advanced Micro Devices, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ************************************************************************/


// clfft.accessors.cpp : Defines all the getters/setters for the Plan
//

#include "stdafx.h"
#include "private.h"
#include "repo.h"

using std::vector;

clfftStatus clfftGetPlanBatchSize( const clfftPlanHandle plHandle, size_t* batchsize )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanBatchSize" ) );

	*batchsize   = fftPlan->batchsize;
	return CLFFT_SUCCESS;
}

clfftStatus clfftSetPlanBatchSize( clfftPlanHandle plHandle, size_t batchsize )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetPlanBatchSize" ) );

	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked		= false;
	fftPlan->batchsize  = batchsize;
	return CLFFT_SUCCESS;
}

clfftStatus clfftGetPlanContext( const clfftPlanHandle plHandle, cl_context* context )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanContext" ) );

	*context = fftPlan->context;
	return CLFFT_SUCCESS;
}

clfftStatus clfftGetPlanPrecision( const clfftPlanHandle plHandle, clfftPrecision* precision )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanPrecision" ) );

	*precision	= fftPlan->precision;

	return	CLFFT_SUCCESS;
}


clfftStatus clfftSetPlanPrecision( clfftPlanHandle plHandle, clfftPrecision precision )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetPlanPrecision" ) );

	if( precision >= ENDPRECISION )
		return CLFFT_INVALID_ARG_VALUE;

	//	We do not support CLFFT_*_FAST currently
	if( precision == CLFFT_SINGLE_FAST || precision == CLFFT_DOUBLE_FAST )
		return CLFFT_NOTIMPLEMENTED;



	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked		= false;
	fftPlan->precision	= precision;

	return	CLFFT_SUCCESS;
}

clfftStatus clfftGetPlanScale( const clfftPlanHandle plHandle, clfftDirection dir, cl_float* scale )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanScale" ) );

	if( dir >= ENDDIRECTION )
		return CLFFT_INVALID_ARG_VALUE;

	if( dir == CLFFT_FORWARD || dir == CLFFT_MINUS )
		*scale = (cl_float)(fftPlan->forwardScale);
	else
		*scale = (cl_float)(fftPlan->backwardScale);

	return CLFFT_SUCCESS;
}

clfftStatus clfftSetPlanScale( clfftPlanHandle plHandle, clfftDirection dir, cl_float scale )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetPlanScale" ) );

	if( dir >= ENDDIRECTION )
		return CLFFT_INVALID_ARG_VALUE;

	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked		= false;

	if( dir == CLFFT_FORWARD || dir == CLFFT_MINUS )
		fftPlan->forwardScale = scale;
	else
		fftPlan->backwardScale = scale;

	return CLFFT_SUCCESS;
}

clfftStatus clfftGetPlanDim( const clfftPlanHandle plHandle, clfftDim* dim, cl_uint* size )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanDim" ) );

	*dim		= fftPlan->dim;

	switch( fftPlan->dim )
	{
		case CLFFT_1D:
		{
			*size = 1;
		}
			break;
		case CLFFT_2D:
		{
			*size = 2;
		}
			break;
		case CLFFT_3D:
		{
			*size = 3;
		}
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}

	return CLFFT_SUCCESS;
}

clfftStatus clfftSetPlanDim( clfftPlanHandle plHandle, const clfftDim dim )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanDim" ) );

	//	We resize the vectors in the plan to keep their sizes consistent with the value of the dimension
	switch( dim )
	{
		case CLFFT_1D:
		{
			fftPlan->length.resize( 1 );
			fftPlan->inStride.resize( 1 );
			fftPlan->outStride.resize( 1 );
		}
			break;
		case CLFFT_2D:
		{
			fftPlan->length.resize( 2 );
			fftPlan->inStride.resize( 2 );
			fftPlan->outStride.resize( 2 );
		}
			break;
		case CLFFT_3D:
		{
			fftPlan->length.resize( 3 );
			fftPlan->inStride.resize( 3 );
			fftPlan->outStride.resize( 3 );
		}
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}

	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked	= false;
	fftPlan->dim	= dim;

	return CLFFT_SUCCESS;
}

clfftStatus clfftGetPlanLength( const clfftPlanHandle plHandle, const clfftDim dim, size_t* clLengths )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanLength" ) );

	if( clLengths == NULL )
		return CLFFT_INVALID_HOST_PTR;

	if( fftPlan->length.empty( ) )
		return CLFFT_INVALID_ARG_INDEX;

	switch( dim )
	{
		case CLFFT_1D:
		{
			clLengths[ DimX ] = fftPlan->length[ DimX ];
		}
			break;
		case CLFFT_2D:
		{
			if( fftPlan->length.size( ) < 2 )
				return CLFFT_INVALID_ARG_INDEX;

			clLengths[ DimX ] = fftPlan->length[ DimX ];
			clLengths[ DimY ] = fftPlan->length[ DimY ];
		}
			break;
		case CLFFT_3D:
        {
			if( fftPlan->length.size( ) < 3 )
				return CLFFT_INVALID_ARG_INDEX;

			clLengths[ DimX ] = fftPlan->length[ DimX ];
			clLengths[ DimY ] = fftPlan->length[ DimY ];
			clLengths[ DimZ ] = fftPlan->length[ DimZ ];
		}
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}
	return	CLFFT_SUCCESS;
}

clfftStatus clfftSetPlanLength( clfftPlanHandle plHandle, const clfftDim dim, const size_t* clLengths )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetPlanLength" ) );

	if( clLengths == NULL )
		return CLFFT_INVALID_HOST_PTR;

	//	Simplest to clear any previous contents, because it's valid for user to shrink dimension
	fftPlan->length.clear( );
	switch( dim )
	{
		case CLFFT_1D:
		{
			//	Minimum length size is 1
			if( clLengths[ DimX ] == 0 )
				return CLFFT_INVALID_ARG_VALUE;

			if( !IsASupportedLength( clLengths[ DimX ] ) )
				return CLFFT_NOTIMPLEMENTED;

			fftPlan->length.push_back( clLengths[ DimX ] );
		}
			break;
		case CLFFT_2D:
		{
			//	Minimum length size is 1
			if( clLengths[ DimX ] == 0 || clLengths[ DimY ] == 0 )
				return CLFFT_INVALID_ARG_VALUE;

			if(!fftPlan->transflag)
			{
				if( !IsASupportedLength( clLengths[ DimX ] ) || !IsASupportedLength( clLengths[ DimY ] ) )
				{
					return CLFFT_NOTIMPLEMENTED;
				}
			}

			fftPlan->length.push_back( clLengths[ DimX ] );
			fftPlan->length.push_back( clLengths[ DimY ] );
		}
			break;
		case CLFFT_3D:
		{
			//	Minimum length size is 1
			if( clLengths[ DimX ] == 0 || clLengths[ DimY ] == 0 || clLengths[ DimZ ] == 0)
				return CLFFT_INVALID_ARG_VALUE;

			if( !IsASupportedLength( clLengths[ DimX ] ) || !IsASupportedLength( clLengths[ DimY ] ) ||
				!IsASupportedLength( clLengths[ DimZ ] ) )
			{
				return CLFFT_NOTIMPLEMENTED;
			}

			fftPlan->length.push_back( clLengths[ DimX ] );
			fftPlan->length.push_back( clLengths[ DimY ] );
			fftPlan->length.push_back( clLengths[ DimZ ] );
		}
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}

	fftPlan->dim = dim;

	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked	= false;

	return	CLFFT_SUCCESS;
}

clfftStatus clfftGetPlanInStride( const clfftPlanHandle plHandle, const clfftDim dim, size_t* clStrides )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanInStride" ) );

	if( clStrides == NULL )
		return CLFFT_INVALID_HOST_PTR;

	switch( dim )
	{
		case CLFFT_1D:
		{
			if( fftPlan->inStride.size( ) > 0 )
				clStrides[ DimX ] = fftPlan->inStride[ DimX ];
			else
				return CLFFT_INVALID_ARG_INDEX;
		}
			break;
		case CLFFT_2D:
		{
			if( fftPlan->inStride.size( ) > 1 )
			{
				clStrides[ DimX ] = fftPlan->inStride[ DimX ];
				clStrides[ DimY ] = fftPlan->inStride[ DimY ];
			}
			else
				return CLFFT_INVALID_ARG_INDEX;
		}
			break;
		case CLFFT_3D:
		{
			if( fftPlan->inStride.size( ) > 2 )
			{
				clStrides[ DimX ] = fftPlan->inStride[ DimX ];
				clStrides[ DimY ] = fftPlan->inStride[ DimY ];
				clStrides[ DimZ ] = fftPlan->inStride[ DimZ ];
			}
			else
				return CLFFT_INVALID_ARG_INDEX;
		}
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}

	return CLFFT_SUCCESS;
}

clfftStatus clfftSetPlanInStride( clfftPlanHandle plHandle, const clfftDim dim, size_t* clStrides )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetPlanInStride" ) );

	if( clStrides == NULL )
		return CLFFT_INVALID_HOST_PTR;

	//	Simplest to clear any previous contents, because it's valid for user to shrink dimension
	fftPlan->inStride.clear( );
	switch( dim )
	{
		case CLFFT_1D:
		{
			fftPlan->inStride.push_back( clStrides[ DimX ] );
		}
			break;
		case CLFFT_2D:
		{
			fftPlan->inStride.push_back( clStrides[ DimX ] );
			fftPlan->inStride.push_back( clStrides[ DimY ] );
		}
			break;
		case CLFFT_3D:
		{
			fftPlan->inStride.push_back( clStrides[ DimX ] );
			fftPlan->inStride.push_back( clStrides[ DimY ] );
			fftPlan->inStride.push_back( clStrides[ DimZ ] );
		}
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}

	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked	= false;

	return CLFFT_SUCCESS;
}

clfftStatus clfftGetPlanOutStride( const clfftPlanHandle plHandle, const clfftDim dim, size_t* clStrides )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanOutStride" ) );

	if( clStrides == NULL )
		return CLFFT_INVALID_HOST_PTR;

	switch( dim )
	{
		case CLFFT_1D:
		{
			if( fftPlan->outStride.size( ) > 0 )
				clStrides[ DimX ] = fftPlan->outStride[ DimX ];
			else
				return CLFFT_INVALID_ARG_INDEX;
		}
			break;
		case CLFFT_2D:
		{
			if( fftPlan->outStride.size( ) > 1 )
			{
				clStrides[ DimX ] = fftPlan->outStride[ DimX ];
				clStrides[ DimY ] = fftPlan->outStride[ DimY ];
			}
			else
				return CLFFT_INVALID_ARG_INDEX;
		}
			break;
		case CLFFT_3D:
		{
			if( fftPlan->outStride.size( ) > 2 )
			{
				clStrides[ DimX ] = fftPlan->outStride[ DimX ];
				clStrides[ DimY ] = fftPlan->outStride[ DimY ];
				clStrides[ DimZ ] = fftPlan->outStride[ DimZ ];
			}
			else
				return CLFFT_INVALID_ARG_INDEX;
		}
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}

	return CLFFT_SUCCESS;
}

clfftStatus clfftSetPlanOutStride( clfftPlanHandle plHandle, const clfftDim dim, size_t* clStrides )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetPlanOutStride" ) );

	if( clStrides == NULL )
		return CLFFT_INVALID_HOST_PTR;

	//	Simplest to clear any previous contents, because it's valid for user to shrink dimension
	fftPlan->outStride.clear( );
	switch( dim )
	{
		case CLFFT_1D:
		{
			fftPlan->outStride.push_back( clStrides[ DimX ] );
		}
			break;
		case CLFFT_2D:
		{
			fftPlan->outStride.push_back( clStrides[ DimX ] );
			fftPlan->outStride.push_back( clStrides[ DimY ] );
		}
			break;
		case CLFFT_3D:
		{
			fftPlan->outStride.push_back( clStrides[ DimX ] );
			fftPlan->outStride.push_back( clStrides[ DimY ] );
			fftPlan->outStride.push_back( clStrides[ DimZ ] );
		}
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}

	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked	= false;

	return CLFFT_SUCCESS;
}

clfftStatus clfftGetPlanDistance( const clfftPlanHandle plHandle, size_t* iDist, size_t* oDist )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanDistance" ) );

	*iDist				= fftPlan->iDist;
	*oDist				= fftPlan->oDist;

	return	CLFFT_SUCCESS;
}

clfftStatus clfftSetPlanDistance( clfftPlanHandle plHandle, size_t iDist, size_t oDist )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetPlanDistance" ) );

	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked	= false;
	fftPlan->iDist	= iDist;
	fftPlan->oDist	= oDist;

	return CLFFT_SUCCESS;
}

clfftStatus clfftGetLayout( const clfftPlanHandle plHandle, clfftLayout* iLayout, clfftLayout* oLayout )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetLayout" ) );

	*iLayout			= fftPlan->inputLayout;
	*oLayout			= fftPlan->outputLayout;

	return	CLFFT_SUCCESS;
}

clfftStatus clfftSetLayout( clfftPlanHandle plHandle, clfftLayout iLayout, clfftLayout oLayout )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetLayout" ) );

	//	Basic error checking on parameter
	if( ( iLayout >= ENDLAYOUT ) || ( oLayout >= ENDLAYOUT ) )
		return CLFFT_INVALID_ARG_VALUE;

	//	We currently only support a subset of formats
	switch( iLayout )
	{
		case CLFFT_COMPLEX_INTERLEAVED:
			{
				if( (oLayout == CLFFT_HERMITIAN_INTERLEAVED) || (oLayout == CLFFT_HERMITIAN_PLANAR) || (oLayout == CLFFT_REAL))
					return CLFFT_NOTIMPLEMENTED;
			}
			break;
		case CLFFT_COMPLEX_PLANAR:
			{
				if( (oLayout == CLFFT_HERMITIAN_INTERLEAVED) || (oLayout == CLFFT_HERMITIAN_PLANAR) || (oLayout == CLFFT_REAL))
					return CLFFT_NOTIMPLEMENTED;
			}
			break;
		case CLFFT_HERMITIAN_INTERLEAVED:
			{
				if(oLayout != CLFFT_REAL) return CLFFT_NOTIMPLEMENTED;
			}
			break;
		case CLFFT_HERMITIAN_PLANAR:
			{
				if(oLayout != CLFFT_REAL) return CLFFT_NOTIMPLEMENTED;
			}
			break;
		case CLFFT_REAL:
			{
				if((oLayout == CLFFT_REAL) || (oLayout == CLFFT_COMPLEX_INTERLEAVED) || (oLayout == CLFFT_COMPLEX_PLANAR))
					return CLFFT_NOTIMPLEMENTED;
			}
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}

	//	We currently only support a subset of formats
	switch( oLayout )
	{
		case CLFFT_COMPLEX_PLANAR:
		case CLFFT_COMPLEX_INTERLEAVED:
		case CLFFT_HERMITIAN_INTERLEAVED:
		case CLFFT_HERMITIAN_PLANAR:
		case CLFFT_REAL:
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}

	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked	= false;
	fftPlan->inputLayout	= iLayout;
	fftPlan->outputLayout	= oLayout;

	return	CLFFT_SUCCESS;
}

clfftStatus clfftGetResultLocation( const clfftPlanHandle plHandle, clfftResultLocation* placeness )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetResultLocation" ) );

	*placeness	= fftPlan->placeness;

	return	CLFFT_SUCCESS;
}

clfftStatus clfftSetResultLocation( clfftPlanHandle plHandle, clfftResultLocation placeness )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetResultLocation" ) );

	//	Basic error checking on parameter
	if( placeness >= ENDPLACE )
		return CLFFT_INVALID_ARG_VALUE;

	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked		= false;
	fftPlan->placeness	= placeness;

	return	CLFFT_SUCCESS;
}


clfftStatus clfftGetPlanTransposeResult( const clfftPlanHandle plHandle, clfftResultTransposed * transposed )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetResultLocation" ) );

	*transposed	= fftPlan->transposed;

	return	CLFFT_SUCCESS;
}

clfftStatus clfftSetPlanTransposeResult( clfftPlanHandle plHandle, clfftResultTransposed transposed )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetResultLocation" ) );

	//	Basic error checking on parameter
	if( transposed >= ENDTRANSPOSED )
		return CLFFT_INVALID_ARG_VALUE;

	//	If we modify the state of the plan, we assume that we can't trust any pre-calculated contents anymore
	fftPlan->baked		= false;
	fftPlan->transposed	= transposed;

	return	CLFFT_SUCCESS;
}

clfftStatus clfftGetTmpBufSize( const clfftPlanHandle plHandle, size_t* buffersize )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanBatchSize" ) );

	if (fftPlan->baked == true)
	{
		*buffersize   = fftPlan->tmpBufSize;
		return CLFFT_SUCCESS;
	}

	return CLFFT_INVALID_OPERATION;
}

clfftStatus clfftLocalMemSize( const clfftPlanHandle plHandle, cl_ulong* local_mem_size )
{
	FFTRepo& repo = FFTRepo::getInstance( );
	FFTPlan* plan = NULL;
	lockRAII* lock = NULL;

	OPENCL_V( repo.getPlan( plHandle, plan, lock ), _T( "repo.getPlan failed" ) );
	scopedLock sLock( *lock, _T( "clfftLocalMemSize" ) );

	*local_mem_size = plan->envelope.limit_LocalMemSize;
	return CLFFT_SUCCESS;
}

clfftStatus clfftSetPlanCallback(clfftPlanHandle plHandle, const char* funcName, 
								 const char* funcString, int localMemSize, 
								 clfftCallbackType callbackType, cl_mem *userdata, int numUserdataBuffers)
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftSetPlanCallback" ) );

	switch (callbackType)
	{
	case PRECALLBACK:
		{
			ARG_CHECK(funcName != NULL);
			ARG_CHECK(funcString != NULL);
			ARG_CHECK(numUserdataBuffers >= 0);

			//	We do not currently support multiple user data buffers
			if( numUserdataBuffers > 1 )
				return CLFFT_NOTIMPLEMENTED;

			fftPlan->hasPreCallback = true;

			fftPlan->preCallback.funcname = funcName;
			fftPlan->preCallback.funcstring = funcString;
			fftPlan->preCallback.localMemSize = (localMemSize > 0) ? localMemSize : 0;

			cl_mem userdataBuf = NULL;
			
			if (userdata)
				userdataBuf = userdata[0];

			fftPlan->precallUserData = userdataBuf;
		}

		break;
	case POSTCALLBACK:
		{
			ARG_CHECK(funcName != NULL);
			ARG_CHECK(funcString != NULL);
			ARG_CHECK(numUserdataBuffers >= 0);

			//	We do not currently support multiple user data buffers
			if( numUserdataBuffers > 1 )
		return CLFFT_NOTIMPLEMENTED;

			fftPlan->hasPostCallback = true;
			fftPlan->postCallbackParam.funcname = funcName;
			fftPlan->postCallbackParam.funcstring = funcString;
			fftPlan->postCallbackParam.localMemSize = (localMemSize > 0) ? localMemSize : 0;

			cl_mem userdataBuf = NULL;
			
			if (userdata)
				userdataBuf = userdata[0];

			fftPlan->postcallUserData = userdataBuf;
		}
		break;
	default:
		ARG_CHECK (false);
	}

	return	CLFFT_SUCCESS;
}