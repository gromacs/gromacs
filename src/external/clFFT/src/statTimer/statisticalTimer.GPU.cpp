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


// StatTimer.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <cassert>
#include <iomanip>
#include <functional>
#include <cmath>
#include "statisticalTimer.GPU.h"
#include "../library/private.h"

//	Format an unsigned number with comma thousands separator
//
template< typename T >		// T could be 32-bit or 64-bit
std::basic_string<TCHAR> commatize (T number)
{
	static TCHAR scratch [8*sizeof(T)];

	register TCHAR * ptr = scratch + countOf( scratch );
	*(--ptr) = 0;

	for (int digits = 3; ; )
	{
		*(--ptr) = '0' + int (number % 10);
		number /= 10;
		if (0 == number)
			break;
		if (--digits <= 0)
		{
			*(--ptr) = ',';
			digits = 3;
		}
	}

	return std::basic_string<TCHAR> (ptr);
}

//	Functor object to help with accumulating values in vectors
template< typename T >
struct Accumulator: public std::unary_function< T, void >
{
	T acc;

	Accumulator( ): acc( 0 ) {}
	void operator( )(T x) { acc += x; }
};

//	Functor object to help with accumulating values in vectors
template< >
struct Accumulator< StatData >
{
	StatData acc;

	Accumulator( ) {}
	void operator( )( const StatData& x )
	{
		acc.deltaNanoSec += x.deltaNanoSec;
	}
};


//	Unary predicate used for remove_if() algorithm
//	Currently, RangeType is expected to be a floating point type, and ValType an integer type
template< typename T, typename R >
struct PruneRange: public std::binary_function< T, R, bool >
{
	R lower, upper;

	PruneRange( R mean, R stdev ): lower( mean-stdev ), upper( mean+stdev ) {}

	bool operator( )( T val )
	{
		//	These comparisons can be susceptible to signed/unsigned casting problems
		//	This is why we cast ValType to RangeType, because RangeType should always be floating and signed
		if( static_cast< R >( val ) < lower )
			return true;
		else if( static_cast< R >( val ) > upper )
			return true;

		return false;
	}
};

//	Template specialization for StatData datatypes
template< >
struct PruneRange< StatData, cl_double >
{
	StatData mean;
	cl_double stdDev;

	PruneRange( StatData m, cl_double s ): mean( m ), stdDev( s ) {}

	bool operator( )( StatData val )
	{
		//	These comparisons can be susceptible to signed/unsigned casting problems
		//	This is why we cast ValType to RangeType, because RangeType should always be floating and signed
		if( val.doubleNanoSec < (mean.doubleNanoSec - stdDev) )
			return true;
		else if( val.doubleNanoSec > (mean.doubleNanoSec + stdDev) )
			return true;

		return false;
	}
};

//	Sorting operator for struct StatData, such that it can be used in a map
bool operator<( const StatData& lhs, const StatData& rhs)
{
	if( lhs.deltaNanoSec < rhs.deltaNanoSec )
		return true;
	else
		return false;
}

GpuStatTimer&
GpuStatTimer::getInstance( )
{
	static	GpuStatTimer	timer;
	return	timer;
}

GpuStatTimer::GpuStatTimer( ): nEvents( 0 ), nSamples( 0 ), currID( 0 ), currSample( 0 ), currRecord( 0 )
{}

GpuStatTimer::~GpuStatTimer( )
{}

void
GpuStatTimer::Clear( )
{
	labelID.clear( );
	timerData.clear( );

	nEvents = 0;
	nSamples = 0;
	currID = 0;
	currSample = 0;
	currRecord = 0;
}

//	The caller can pre-allocate memory, to improve performance.
//	nEvents is an approximate value for how many seperate events the caller will think
//	they will need, and nSamples is a hint on how many samples we think we will take
//	per event
void
GpuStatTimer::Reserve( size_t nE, size_t nS )
{
	Clear( );
	nEvents		= std::max< size_t >( 1, nE );
	nSamples	= std::max< size_t >( 1, nS );

	labelID.reserve( nEvents );
	timerData.resize( nEvents );
}

void
GpuStatTimer::Reset( )
{
	if( nEvents == 0 || nSamples == 0 )
		throw	std::runtime_error( "StatisticalTimer::Reserve( ) was not called before Reset( )" );

	ReleaseEvents();
	Reserve( nEvents, nSamples );

	return;
}

void
GpuStatTimer::setNormalize( bool norm )
{
}

void
GpuStatTimer::Start( size_t id )
{
	currID		= id;
	currSample	= 0;
}

void
GpuStatTimer::Stop( size_t id )
{
	++currRecord;
}

void
GpuStatTimer::AddSample( clfftPlanHandle plHandle, FFTPlan* plan, cl_kernel kern, cl_uint numEvents, cl_event* ev,
	const std::vector< size_t >& gWorkSize, const std::vector< size_t >& lWorkSize )
{
	if( (numEvents != 0) && (ev == NULL) )
		return;

	if( timerData.empty( ) )
		return;

	for( size_t i = 0; i < numEvents; ++i )
	{
		::clRetainEvent(ev[i]);
	}

	if( currRecord == 0 )
	{
		timerData.at( currID ).push_back( StatDataVec( ) );
		timerData.at( currID ).back( ).reserve( nSamples );
		timerData.at( currID ).back( ).push_back( StatData( plHandle, plan, kern, numEvents, ev, gWorkSize, lWorkSize) );
	}
	else
	{
		timerData.at( currID ).at( currSample )
			.push_back( StatData( plHandle, plan, kern, numEvents, ev, gWorkSize, lWorkSize ) );
		++currSample;
	}
}

//	This function's purpose is to provide a mapping from a 'friendly' human readable text string
//	to an index into internal data structures.
size_t
GpuStatTimer::getUniqueID( const std::string& label, cl_uint groupID )
{
	//	I expect labelID will hardly ever grow beyond 30, so it's not of any use
	//	to keep this sorted and do a binary search

	idPair	sItem	= std::make_pair( label, groupID );

	idVector::iterator	iter;
	iter	= std::find( labelID.begin(), labelID.end(), sItem );

	if( iter != labelID.end( ) )
		return	std::distance( labelID.begin( ), iter );

	labelID.push_back( sItem );

	return	labelID.size( ) - 1;

}

void GpuStatTimer::ReleaseEvents()
{
	for( cl_uint id = 0; id < labelID.size( ); ++id )
	{
		for( size_t s = 0; s < timerData.at( id ).size( ); ++s )
		{
			for( size_t n = 0; n < timerData.at( id ).at( s ).size( ); ++n )
			{
				StatData& sd = timerData[ id ][ s ][ n ];

				for( size_t i = 0; i < sd.outEvents.size( ); ++i )
				{
					::clReleaseEvent(sd.outEvents[ i ]);
				}

			}
		}
	}
}

void GpuStatTimer::queryOpenCL( size_t id )
{
	for( size_t s = 0; s < timerData.at( id ).size( ); ++s )
	{
		for( size_t n = 0; n < timerData.at( id ).at( s ).size( ); ++n )
		{
			StatData& sd = timerData[ id ][ s ][ n ];

			cl_ulong profStart, profEnd = 0;
			sd.deltaNanoSec = 0;

			for( size_t i = 0; i < sd.outEvents.size( ); ++i )
			{
				if( ::clGetEventProfilingInfo( sd.outEvents[ i ], CL_PROFILING_COMMAND_START, sizeof( cl_ulong ), &profStart, NULL ) != CL_SUCCESS )
				{
					profStart = 0;
				}

				if( ::clGetEventProfilingInfo( sd.outEvents[ i ], CL_PROFILING_COMMAND_END, sizeof( cl_ulong ), &profEnd, NULL ) != CL_SUCCESS )
				{
					profEnd = 0;
				}
				sd.deltaNanoSec += (profEnd - profStart);
			}

			sd.doubleNanoSec = static_cast< cl_double >( sd.deltaNanoSec );
		}
	}
}

std::vector< StatData >
GpuStatTimer::getMean( size_t id )
{
	//	Prep the data; query openCL for the timer information
	queryOpenCL( id );

	std::vector< StatData > meanVec;
	for( size_t s = 0; s < timerData.at( id ).size( ); ++s )
	{
		Accumulator< StatData > sum = std::for_each( timerData.at( id ).at( s ).begin( ), timerData.at( id ).at( s ).end( ),
			Accumulator< StatData >() );

		StatData tmp = timerData[ id ][ s ].front( );
		tmp.doubleNanoSec = static_cast< cl_double >( sum.acc.deltaNanoSec ) / timerData.at( id ).at( s ).size( );
		meanVec.push_back( tmp );
	}

	return meanVec;
}

std::vector< StatData >
GpuStatTimer::getVariance( size_t id )
{
	std::vector< StatData > variance = getMean( id );

	for( cl_uint v = 0; v < variance.size( ); ++v )
	{
		double sum = 0;
		for( cl_uint n = 0; n < timerData[ id ][ v ].size( ); ++n )
		{
			cl_double	diff	= static_cast< cl_double >( timerData[ id ][ v ][ n ].deltaNanoSec ) - variance[ v ].doubleNanoSec;
			diff	*= diff;
			sum		+= diff;
		}

		variance[ v ].doubleNanoSec = sum / timerData[ id ][ v ].size( );
	}

	return variance;
}

std::vector< StatData >
GpuStatTimer::getStdDev( size_t id )
{
	std::vector< StatData > stddev = getVariance( id );

	for( cl_uint v = 0; v < stddev.size( ); ++v )
	{
		stddev[ v ].doubleNanoSec = sqrt( stddev[ v ].doubleNanoSec );
	}

	return stddev;
}

std::vector< StatData >
GpuStatTimer::getAverageTime( size_t id )
{
	return getMean( id );
}

std::vector< StatData >
GpuStatTimer::getMinimumTime( size_t id )
{
	//	Prep the data; query openCL for the timer information
	queryOpenCL( id );

	std::vector< StatData > minTime;
	for( size_t s = 0; s < timerData.at( id ).size( ); ++s )
	{
		StatDataVec::iterator iter
			= std::min_element( timerData.at( id ).at( s ).begin( ), timerData.at( id ).at( s ).end( ) );

		if( iter != timerData.at( id ).at( s ).end( ) )
		{
			iter->doubleNanoSec = static_cast< cl_double >( iter->deltaNanoSec ) / timerData.at( id ).at( s ).size( );
			minTime.push_back( *iter );
		}
		else
			return std::vector< StatData >( );
	}

	return minTime;
}

std::vector< size_t >
GpuStatTimer::pruneOutliers( size_t id , cl_double multiple )
{
	std::vector< StatData > mean	= getMean( id );
	std::vector< StatData > stdDev	= getStdDev( id );

	std::vector< size_t > totalPrune;
	for( size_t s = 0; s < timerData.at( id ).size( ); ++s )
	{
		//	Look on p. 379, "The C++ Standard Library"
		//	std::remove_if does not actually erase, it only copies elements, it returns new 'logical' end
		StatDataVec::iterator newEnd	= std::remove_if( timerData.at( id ).at( s ).begin( ), timerData.at( id ).at( s ).end( ),
			PruneRange< StatData,cl_double >( mean[ s ], multiple * stdDev[ s ].doubleNanoSec ) );

		StatDataVec::difference_type dist	= std::distance( newEnd, timerData.at( id ).at( s ).end( ) );

		if( dist != 0 )
			timerData.at( id ).at( s ).erase( newEnd, timerData.at( id ).at( s ).end( ) );

		totalPrune.push_back( dist );
	}

	return totalPrune;
}

size_t
GpuStatTimer::pruneOutliers( cl_double multiple )
{
	const int tableWidth = 60;
	const int tableHalf = tableWidth / 2;
	const int tableThird = tableWidth / 3;
	const int tableFourth = tableWidth / 4;
	const int tableFifth = tableWidth / 5;

	//	Print label of timer, in a header
	std::string header( "StdDev" );
	size_t	sizeTitle = (header.size( ) + 6) /2;

	std::cout << std::endl;
	std::cout << std::setfill( '=' ) << std::setw( tableHalf ) << header << " ( " << multiple << " )"
			<< std::setw( tableHalf - sizeTitle ) << "=" << std::endl;
	tout << std::setfill( _T( ' ' ) );

	size_t tCount = 0;
	for( cl_uint l = 0; l < labelID.size( ); ++l )
	{
		std::vector< size_t > lCount = pruneOutliers( l , multiple );

		for( cl_uint c = 0; c < lCount.size( ); ++c )
		{
			std::cout << labelID[l].first << "[ " << c << " ]" << ": Pruning " << lCount[ c ] << " samples out of " << currRecord << std::endl;
			tCount += lCount[ c ];
		}
	}

	return tCount;
}

void
GpuStatTimer::Print( )
{
	const int tableWidth = 60;
	const int tableHalf = tableWidth / 2;
	const int tableThird = tableWidth / 3;
	const int tableFourth = tableWidth / 4;
	const int tableFifth = tableWidth / 5;

	for( cl_uint id = 0; id < labelID.size( ); ++id )
	{
		size_t	halfString = labelID[ id ].first.size( ) / 2;

		//	Print label of timer, in a header
		std::cout << std::endl << std::setw( tableHalf + halfString ) << std::setfill( '=' ) << labelID[ id ].first
				<< std::setw( tableHalf - halfString ) << "=" << std::endl;
		tout << std::setfill( _T( ' ' ) );

		std::vector< StatData > mean	= getMean( id );

		//	Print each individual dimension
		tstringstream catLengths;
		for( cl_uint t = 0; t < mean.size( ); ++t )
		{
			cl_double time		= 0;
			if( mean[ t ].kernel == NULL )
			{
				for( cl_uint m = 0; m < t; ++m )
				{
					if( mean[ m ].plHandle == mean[ t ].planX ||
						mean[ m ].plHandle == mean[ t ].planY ||
						mean[ m ].plHandle == mean[ t ].planZ ||
						mean[ m ].plHandle == mean[ t ].planTX ||
						mean[ m ].plHandle == mean[ t ].planTY ||
						mean[ m ].plHandle == mean[ t ].planTZ ||
						mean[ m ].plHandle == mean[ t ].planRCcopy ||
						mean[ m ].plHandle == mean[ t ].planCopy )
					{
						time	+= mean[ m ].doubleNanoSec;
					}
				}
				mean[ t ].doubleNanoSec = time;
			}
			else
			{
				time	= mean[ t ].doubleNanoSec;
			}
			double gFlops = mean[ t ].calcFlops( ) / time;

			tout << std::setw( tableFourth ) << _T( "Handle:" )
				<< std::setw( tableThird )  << mean[ t ].plHandle << std::endl;

			if( mean[ t ].kernel != 0 )
			{
				tout << std::setw( tableFourth ) << _T( "Kernel:" )
					<< std::setw( tableThird )  << reinterpret_cast<void*>( mean[ t ].kernel ) << std::endl;
			}

			if( ( mean[ t ].planX + mean[ t ].planY + mean[ t ].planZ ) > 0 ||
				( mean[ t ].planTX + mean[ t ].planTY + mean[ t ].planTZ ) > 0  || 
				( mean[ t ].planRCcopy + mean[ t ].planCopy ) > 0 )
			{
				tout << std::setw( tableFourth ) << _T( "Child Handles:" );
				catLengths.str( _T( "" ) );
				catLengths << _T( "(" );
				if( mean[ t ].planX != 0 )
					catLengths << mean[ t ].planX;
				if( mean[ t ].planTX != 0 )
				{
					catLengths << _T( "," );
					catLengths << mean[ t ].planTX;
				}
				if( mean[ t ].planY != 0 )
				{
					catLengths << _T( "," );
					catLengths << mean[ t ].planY;
				}
				if( mean[ t ].planTY != 0 )
				{
					catLengths << _T( "," );
					catLengths << mean[ t ].planTY;
				}
				if( mean[ t ].planZ != 0 )
				{
					catLengths << _T( "," );
					catLengths << mean[ t ].planZ;
				}
				if( mean[ t ].planTZ != 0 )
				{
					catLengths << _T( "," );
					catLengths << mean[ t ].planTZ;
				}
				if( mean[ t ].planRCcopy != 0 )
				{
					catLengths << _T( "," );
					catLengths << mean[ t ].planRCcopy;
				}
				if( mean[ t ].planCopy != 0 )
				{
					catLengths << _T( "," );
					catLengths << mean[ t ].planCopy;
				}
				catLengths << _T( ")" );
				tout << std::setw( tableThird ) << catLengths.str( ) << std::endl;
			}

			if( mean[ t ].outEvents.size( ) != 0 )
			{
				tout << std::setw( tableFourth ) << _T( "OutEvents:" ) << std::setw( tableThird );
				for( size_t i = 0; i < mean[ t ].outEvents.size( ); ++i )
				{
					tout << mean[ t ].outEvents[ i ];
					if( i < (mean[ t ].outEvents.size( )-1) )
					{
						tout << _T( "," ) << std::endl;
						tout << std::setw( tableFourth+tableThird );
					}
				}
				tout << std::endl;
			}


			tout << std::setw(tableFourth) << _T("Generator:");
			switch(mean[t].gen)
			{
			case Stockham:				tout << std::setw(tableThird) << _T("Stockham"); break;
			case Transpose_GCN:			tout << std::setw(tableThird) << _T("Transpose_GCN"); break;
			case Transpose_SQUARE:		tout << std::setw(tableThird) << _T("Transpose_SQUARE"); break;
			case Transpose_NONSQUARE:	tout << std::setw(tableThird) << _T("Transpose_NONSQUARE"); break;
			case Copy:					tout << std::setw(tableThird) << _T("Copy"); break;
			}
			tout << std::endl;


			tout << std::setw( tableFourth ) << _T( "Length:" );
			catLengths.str( _T( "" ) );
			catLengths << _T( "(" );
			for( size_t i = 0; i < mean[ t ].lengths.size( ); ++i )
			{
				catLengths << mean[ t ].lengths.at( i );
				if( i < (mean[ t ].lengths.size( )-1) )
					catLengths << _T( "," );
			}
			catLengths << _T( ")" );
			tout << std::setw( tableThird ) << catLengths.str( ) << std::endl;

			if( mean[ t ].batchSize > 1 )
			{
				tout << std::setw( tableFourth ) << _T( "Batch:" )
					<< std::setw( tableThird )  << mean[ t ].batchSize << std::endl;
			}

			tout << std::setw(tableFourth) << _T("Placeness:") << std::setw(tableThird)
				<< ( mean[t].placeness == CLFFT_INPLACE ? "InPlace": "OutOfPlace" ) << std::endl;

			tout << std::setw(tableFourth) << _T("Input Dist:") << std::setw(tableThird) << mean[t].iDist << std::endl;
			tout << std::setw(tableFourth) << _T("Output Dist:") << std::setw(tableThird) << mean[t].oDist << std::endl;

			tout << std::setw( tableFourth ) << _T( "Input Stride:" );

			catLengths.str( _T( "" ) );
			catLengths << _T( "(" );
			for( size_t i = 0; i < mean[ t ].inStride.size( ); ++i )
			{
				catLengths << mean[ t ].inStride.at( i );
				if( i < (mean[ t ].inStride.size( )-1) )
					catLengths << _T( "," );
			}
			catLengths << _T( ")" );
			tout << std::setw( tableThird ) << catLengths.str( ) << std::endl;

			tout << std::setw( tableFourth ) << _T( "Output Stride:" );

			catLengths.str( _T( "" ) );
			catLengths << _T( "(" );
			for( size_t i = 0; i < mean[ t ].outStride.size( ); ++i )
			{
				catLengths << mean[ t ].outStride.at( i );
				if( i < (mean[ t ].outStride.size( )-1) )
					catLengths << _T( "," );
			}
			catLengths << _T( ")" );
			tout << std::setw( tableThird ) << catLengths.str( ) << std::endl;

			if( mean[ t ].enqueueWorkSize.size( ) != 0 )
			{
				tout << std::setw( tableFourth ) << _T( "Global Work:" );
				catLengths.str( _T( "" ) );
				catLengths << _T( "(" );
				for( size_t i = 0; i < mean[ t ].enqueueWorkSize.size( ); ++i )
				{
					catLengths << mean[ t ].enqueueWorkSize.at( i );
					if( i < (mean[ t ].enqueueWorkSize.size( )-1) )
						catLengths << _T( "," );
				}
				catLengths << _T( ")" );
				tout << std::setw( tableThird ) << catLengths.str( ) << std::endl;

				tout << std::setw(tableFourth) << _T("Local Work:");
				catLengths.str(_T(""));
				catLengths << _T("(");
				for (size_t i = 0; i < mean[t].enqueueLocalWorkSize.size(); ++i)
				{
					catLengths << mean[t].enqueueLocalWorkSize.at(i);
					if (i < (mean[t].enqueueLocalWorkSize.size() - 1))
						catLengths << _T(",");
				}
				catLengths << _T(")");
				tout << std::setw(tableThird) << catLengths.str() << std::endl;
			}

			tout << std::setw( tableFourth ) << _T( "Gflops:" )
				<< std::setw( 2*tableFourth ) << gFlops << std::endl;
			tout << std::setw( tableFourth ) << _T( "Time (ns):" )
				<< std::setw( 3*tableFourth ) << commatize( static_cast< cl_ulong >( time ) ) << std::endl;
			tout << std::endl;
		}
	}
}

//	Defining an output print operator
std::ostream&
operator<<( std::ostream& os, const GpuStatTimer& st )
{
	//if( st.clkTicks.empty( ) )
	//	return	os;

	//std::ios::fmtflags bckup	= os.flags( );

	//for( cl_uint l = 0; l < st.labelID.size( ); ++l )
	//{
	//	cl_ulong min	= 0;
	//	clkVector::const_iterator iter	= std::min_element( st.clkTicks.at( l ).begin( ), st.clkTicks.at( l ).end( ) );

	//	if( iter != st.clkTicks.at( l ).end( ) )
	//		min		= *iter;

	//	os << st.labelID[l].first << ", " << st.labelID[l].second << std::fixed << std::endl;
	//	os << "Min:," << min << std::endl;
	//	os << "Mean:," << st.getMean( l ) << std::endl;
	//	os << "StdDev:," << st.getStdDev( l ) << std::endl;
	//	os << "AvgTime:," << st.getAverageTime( l ) << std::endl;
	//	os << "MinTime:," << st.getMinimumTime( l ) << std::endl;

	//	//for( cl_uint	t = 0; t < st.clkTicks[l].size( ); ++t )
	//	//{
	//	//	os << st.clkTicks[l][t]<< ",";
	//	//}
	//	os << "\n" << std::endl;

	//}

	//os.flags( bckup );

	return	os;
}
