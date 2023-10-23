#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"

#include "sz.h"

#include <stdio.h>  // for printf
#include <stdlib.h>

/* Test Suite setup and cleanup functions: */

int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }

/************* Test case functions ****************/

void test_encode_decode_uniform_distribution_256(void)
{
	//initialization
	maxRangeRadius = 128;
	stateNum = maxRangeRadius*2;
	allNodes = maxRangeRadius*4;
	
	intvCapacity = maxRangeRadius*2;
	intvRadius = maxRangeRadius;
	SZ_Reset();
	
	//testing
	int data[1000000];
	int i = 0;
	for(i=0;i<1000000;i++)
		data[i] = i%256;
	unsigned char* out;
	int outSize;
	encode_withTree(data, 1000000, &out, &outSize);
	printf("outSize=%d\n", outSize);	
	int actual[1000000];
	decode_withTree(out, 1000000, actual);
	
	CU_ASSERT_EQUAL_ARRAY_INT(data, actual, 1000000);
	//CU_ASSERT(outSize<1000000);
	
	free(out);
	SZ_ReleaseHuffman();		
}

void test_encode_decode_linear_distribution_256(void)
{
	//initialization
	maxRangeRadius = 128;
	stateNum = maxRangeRadius*2;
	allNodes = maxRangeRadius*4;
	
	intvCapacity = maxRangeRadius*2;
	intvRadius = maxRangeRadius;

	SZ_Reset();
	
	//testing
	int data[1011840]; 
	int a = 0, i = 0, j = 0,k = 0;
	for(a = 0;a<31;a++)
		for(i=0;i<256;i++) //(1+255)*255/2=32640
		{
			for(j=0;j<i;j++)
			{
				data[k++] = i;
			}
		}
	
	unsigned char* out;
	int outSize;
	encode_withTree(data, 1011840, &out, &outSize);
	printf("outSize=%d\n", outSize);	
	int actual[1011840];
	decode_withTree(out, 1011840, actual);
	
	CU_ASSERT_EQUAL_ARRAY_INT(data, actual, 1011840);
	//CU_ASSERT(outSize<1000000);
	
	free(out);
	SZ_ReleaseHuffman();	
}

void test_encode_decode_extreme_distribution_256(void)
{
	//initialization
	maxRangeRadius = 128;
	stateNum = maxRangeRadius*2;
	allNodes = maxRangeRadius*4;
	
	intvCapacity = maxRangeRadius*2;
	intvRadius = maxRangeRadius;

	SZ_Reset();
	
	//testing
	int data[1000000];
	int i = 0;
	for(i=0;i<1000000;i++)
		data[i] = 0;
	for(i=0;i<10;i++)
		data[i] = 1;
	for(i=10;i<100;i++)
		data[i] = 2;
	
	unsigned char* out;
	int outSize;
	encode_withTree(data, 1000000, &out, &outSize);
	printf("outSize=%d\n", outSize);	
	int actual[1000000];
	decode_withTree(out, 1000000, actual);
	
	CU_ASSERT_EQUAL_ARRAY_INT(data, actual, 1000000);
	//CU_ASSERT(outSize<1000000);
	
	free(out);
	SZ_ReleaseHuffman();	
}

void test_encode_decode_uniform_distribution_65536(void)
{
	//initialization
	maxRangeRadius = 32768;
	stateNum = maxRangeRadius*2;
	allNodes = maxRangeRadius*4;
	
	intvCapacity = maxRangeRadius*2;
	intvRadius = maxRangeRadius;
	SZ_Reset();
	//testing	
	int data[1000000];
	int i = 0;
	for(i=0;i<1000000;i++)
		data[i] = i%65536;
	unsigned char* out;
	int outSize;
	encode_withTree(data, 1000000, &out, &outSize);
	printf("outSize=%d\n", outSize);	
	int actual[1000000];
	decode_withTree(out, 1000000, actual);
	
	CU_ASSERT_EQUAL_ARRAY_INT(data, actual, 1000000);
	//CU_ASSERT(outSize<1000000*2);
	
	free(out);
	SZ_ReleaseHuffman();	
}

void test_encode_decode_extreme_distribution_65536(void)
{
	//initialization
	maxRangeRadius = 32768;
	stateNum = maxRangeRadius*2;
	allNodes = maxRangeRadius*4;
	
	intvCapacity = maxRangeRadius*2;
	intvRadius = maxRangeRadius;
	SZ_Reset();
	//testing	
	int data[1000000];
	int i = 0;
	for(i=0;i<1000000;i++)
		data[i] = 0;
	for(i=0;i<10;i++)
		data[i] = 1;
	for(i=10;i<100;i++)
		data[i] = 2;
					
	unsigned char* out;
	int outSize;
	encode_withTree(data, 1000000, &out, &outSize);
	printf("outSize=%d\n", outSize);	
	int actual[1000000];
	decode_withTree(out, 1000000, actual);
	
	CU_ASSERT_EQUAL_ARRAY_INT(data, actual, 1000000);
	//CU_ASSERT(outSize<1000000*2);
	
	free(out);
	SZ_ReleaseHuffman();	
}

/************* Test Runner Code goes here **************/

int main ( void )
{
   CU_pSuite pSuite = NULL;

   /* initialize the CUnit test registry */
   if ( CUE_SUCCESS != CU_initialize_registry() )
      return CU_get_error();

   /* add a suite to the registry */
   pSuite = CU_add_suite( "test_conf_suite", init_suite, clean_suite );
   if ( NULL == pSuite ) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* add the tests to the suite */
   if ( (NULL == CU_add_test(pSuite, "test_encode_decode_uniform_distribution_256", test_encode_decode_uniform_distribution_256)) ||
	    (NULL == CU_add_test(pSuite, "test_encode_decode_linear_distribution_256", test_encode_decode_linear_distribution_256)) ||
	    (NULL == CU_add_test(pSuite, "test_encode_decode_linear_distribution_256", test_encode_decode_extreme_distribution_256)) ||
        (NULL == CU_add_test(pSuite, "test_encode_decode_uniform_distribution_65536", test_encode_decode_uniform_distribution_65536)) ||
        (NULL == CU_add_test(pSuite, "test_encode_decode_extreme_distribution_65536", test_encode_decode_extreme_distribution_65536)) 
      )
   {
      CU_cleanup_registry();
      return CU_get_error();
   }

   // Run all tests using the basic interface
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   printf("\n");
   CU_basic_show_failures(CU_get_failure_list());
   printf("\n\n");
/*
   // Run all tests using the automated interface
   CU_automated_run_tests();
   CU_list_tests_to_file();

   // Run all tests using the console interface
   CU_console_run_tests();
*/
   /* Clean up registry and return */
   CU_cleanup_registry();
   return CU_get_error();
}
