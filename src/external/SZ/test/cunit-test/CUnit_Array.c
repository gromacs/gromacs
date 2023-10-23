#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include <stdio.h>

void CU_ASSERT_EQUAL_ARRAY_INT(int* actual, int* expected, int count)
{
	int result = 1, i;
	for(i=0;i<count;i++)
	{
		if(actual[i]!=expected[i])
		{
			result = 0;
			break;
		}
	}
	
	if(result==1)
	{	
		CU_ASSERT(CU_TRUE);
	}
	else
	{
		CU_ASSERT(CU_FALSE);
	}
}

void CU_ASSERT_EQUAL_ARRAY_BYTE(unsigned char* actual, unsigned char* expected, int count)
{
	int result = 1, i;
	for(i=0;i<count;i++)
	{
		if(actual[i]!=expected[i])
		{
			result = 0;
			break;
		}
	}
	
	if(result==1)
	{	
		CU_ASSERT(CU_TRUE);
	}
	else
	{
		CU_ASSERT(CU_FALSE);
	}
}

void CU_ASSERT_EQUAL_ARRAY_FLOAT(float* actual, float* expected, int count, double granularity)
{
	int result = 1, i;
	float value = 0;
	for(i=0;i<count;i++)
	{
		value = actual[i] - expected[i];
		if(value < -granularity || value > granularity)
		{
			result = 0;
			break;
		}
	}
	
	if(result==1)
	{	
		CU_ASSERT(CU_TRUE);
	}
	else
	{
		CU_ASSERT(CU_FALSE);
	}	
}

void CU_ASSERT_EQUAL_ARRAY_DOUBLE(double* actual, double* expected, int count, double granularity)
{
	int result = 1, i;
	double value = 0;
	for(i=0;i<count;i++)
	{
		value = actual[i] - expected[i];
		if(value < -granularity || value > granularity)
		{
			result = 0;
			break;
		}
	}
	
	if(result==1)
	{	
		CU_ASSERT(CU_TRUE);
	}
	else
	{
		CU_ASSERT(CU_FALSE);
	}
}

