#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"

#include "sz.h"

#include <stdio.h>  // for printf

/* Test Suite setup and cleanup functions: */

int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }

/************* Test case functions ****************/

void test_bytesToInt_bigEndian(void)
{
	unsigned char bytes[4] = {0,0,0,1};
	CU_ASSERT_EQUAL(bytesToInt_bigEndian(bytes), 1);
	unsigned char bytes2[4] = {1,2,3,4};
	CU_ASSERT_EQUAL(bytesToInt_bigEndian(bytes2), 16909060);
}

void test_intToBytes_bigEndian(void)
{
	int value[2] = {1,16909060};
	unsigned char bytes[4] = {0,0,0,1};
	unsigned char bytes2[4] = {1,2,3,4};
	unsigned char bytesBuf[4];
	intToBytes_bigEndian(bytesBuf, value[0]);
	CU_ASSERT_EQUAL_ARRAY_BYTE(bytesBuf, bytes, 4);
	intToBytes_bigEndian(bytesBuf, value[1]);
	CU_ASSERT_EQUAL_ARRAY_BYTE(bytesBuf, bytes2, 4);
}

void test_intToBytes_bytesToInt_bigEndian(void)
{
	int i = 0, value;
	unsigned char bytesBuf[4];

	for(i=0;i<1000000;i+=100000)
	{
		intToBytes_bigEndian(bytesBuf, i);
		value = bytesToInt_bigEndian(bytesBuf);	
		CU_ASSERT_EQUAL(value, i);
	}
}

void test_bytesToLong_bigEndian(void)
{
	unsigned char bytes[8] = {0, 0, 0, 0, 0, 0, 0, 1};
	CU_ASSERT_EQUAL(bytesToLong_bigEndian(bytes), 1);
	unsigned char bytes2[8] = {0, 0, 0, 0, 1, 2, 3, 4};
	CU_ASSERT_EQUAL(bytesToLong_bigEndian(bytes2), 16909060);
}

void test_longToBytes_bigEndian(void)
{
	long value[2] = {1,16909060};
	unsigned char bytes[8] = {0,0,0,0,0,0,0,1};
	unsigned char bytes2[8] = {0,0,0,0,1,2,3,4};
	unsigned char bytesBuf[8];
	longToBytes_bigEndian(bytesBuf, value[0]);
	CU_ASSERT_EQUAL_ARRAY_BYTE(bytesBuf, bytes, 8);
	longToBytes_bigEndian(bytesBuf, value[1]);
	CU_ASSERT_EQUAL_ARRAY_BYTE(bytesBuf, bytes2, 8);	
}

void test_longToBytes_bytesToLong_bigEndian(void)
{
	long i = 0, value;
	unsigned char bytesBuf[8];

	for(i=1;i<2000000000;i+=500000000)
	{
		longToBytes_bigEndian(bytesBuf, i);
		value = bytesToLong_bigEndian(bytesBuf);	
		CU_ASSERT_EQUAL(value, i);
	}
}

void test_doubleToOSEndianLong(void)
{
	//deprecated
}

void test_floatToOSEndianInt(void)
{
	//deprecated
}

void test_getExponent_float(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{
		CU_ASSERT_EQUAL(getExponent_float(123.123), 6);
		CU_ASSERT_EQUAL(getExponent_float(12.3123), 3);
		CU_ASSERT_EQUAL(getExponent_float(1.23123), 0);
		CU_ASSERT_EQUAL(getExponent_float(0.123123), -4);
		CU_ASSERT_EQUAL(getExponent_float(0.0123123), -7);						
	}
}

void test_getPrecisionReqLength_float(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{	
		CU_ASSERT_EQUAL(getPrecisionReqLength_float(10), 3);
		CU_ASSERT_EQUAL(getPrecisionReqLength_float(1), 0);
		CU_ASSERT_EQUAL(getPrecisionReqLength_float(0.1), -4);
		CU_ASSERT_EQUAL(getPrecisionReqLength_float(0.01), -7);
		CU_ASSERT_EQUAL(getPrecisionReqLength_float(0.001), -10);
	}
}

void test_getExponent_double(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{	
		CU_ASSERT_EQUAL(getExponent_double(123.123), 6);
		CU_ASSERT_EQUAL(getExponent_double(12.3123), 3);
		CU_ASSERT_EQUAL(getExponent_double(1.23123), 0);
		CU_ASSERT_EQUAL(getExponent_double(0.123123), -4);
		CU_ASSERT_EQUAL(getExponent_double(0.0123123), -7);		
	}
}

void test_getPrecisionReqLength_double(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{
		CU_ASSERT_EQUAL(getExponent_double(123.123), 6);
		CU_ASSERT_EQUAL(getExponent_double(12.3123), 3);
		CU_ASSERT_EQUAL(getExponent_double(1.23123), 0);
		CU_ASSERT_EQUAL(getExponent_double(0.123123), -4);
		CU_ASSERT_EQUAL(getExponent_double(0.0123123), -7);	
	}
}

void test_numberOfLeadingZeros_Int(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{	
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(1), 31);
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(10), 28);
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(50), 26);	
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(123), 25);
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(123123123), 5);
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(-123), 0);
	}
}

void test_numberOfLeadingZeros_Long(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{	
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(1), 31);
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(10), 28);
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(50), 26);	
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(123), 25);
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(123123123), 5);
		CU_ASSERT_EQUAL(numberOfLeadingZeros_Int(-123), 0);	
	}
}

void test_getLeadingNumbers_Int(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{	
		CU_ASSERT_EQUAL(getLeadingNumbers_Int(123123123, 123123456), 22);
		CU_ASSERT_EQUAL(getLeadingNumbers_Int(1234, 4567), 19);
	}
}

void test_getLeadingNumbers_Long(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{
		CU_ASSERT_EQUAL(getLeadingNumbers_Long(123123123, 123123456), 54);
		CU_ASSERT_EQUAL(getLeadingNumbers_Long(1234, 4567), 51);
	}
}

void test_bytesToShort(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{	
		unsigned char bytes[2] = {1,3};
		CU_ASSERT_EQUAL(bytesToShort(bytes), 769);
	}
}

void test_bytesToInt(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{	
		unsigned char bytes[4] = {1,2,3,4};
		CU_ASSERT_EQUAL(bytesToInt(bytes), 67305985);
		unsigned char bytes2[4] = {100,50,25,12};
		CU_ASSERT_EQUAL(bytesToInt(bytes2), 202977892);
	}
}

void test_bytesToLong(void)
{
	if(dataEndianType == LITTLE_ENDIAN_DATA)
	{
		long expected = 578437695752307201;
		unsigned char bytes[8] = {1,2,3,4,5,6,7,8};
		CU_ASSERT_EQUAL(bytesToLong(bytes), expected);		
	}
}

void test_bytesToFloat(void)
{
	//see test_floatToBytes_bytesToFloat
}

void test_floatToBytes(void)
{
	//see test_floatToBytes_bytesToFloat
}

void test_floatToBytes_bytesToFloat()
{
	float value = 123.456;
	unsigned char bytes[4];
	floatToBytes(bytes, value);
	float newValue = bytesToFloat(bytes);
	CU_ASSERT_DOUBLE_EQUAL(value, newValue, 1E-4);
}

void test_bytesToDouble(void)
{
	//see test_doubleToBytes_bytesToDouble
}

void test_doubleToBytes(void)
{
	//see test_doubleToBytes_bytesToDouble
}

void test_doubleToBytes_bytesToDouble()
{
	double value = 123.456;
	unsigned char bytes[sizeof(double)];
	doubleToBytes(bytes, value);
	double newValue = bytesToDouble(bytes);
	CU_ASSERT_DOUBLE_EQUAL(value, newValue, 1E-10);
}

void test_extractBytes(void)
{
	//TODO
}

void test_getMaskRightCode(void)
{
	//TODO
}

void test_getLeftMovingCode(void)
{
	//TODO
}

void test_getRightMovingSteps(void)
{
	//TODO
}

void test_getRightMovingCode(void)
{
	//TODO
}

void test_convertByteDataToShortArray(void)
{
	//TODO
}

void test_convertShortArrayToBytes(void)
{
	//TODO
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
   if ( (NULL == CU_add_test(pSuite, "test_bytesToInt_bigEndian", test_bytesToInt_bigEndian)) ||
        (NULL == CU_add_test(pSuite, "test_intToBytes_bigEndian", test_intToBytes_bigEndian)) ||
        (NULL == CU_add_test(pSuite, "test_intToBytes_bytesToInt_bigEndian", test_intToBytes_bytesToInt_bigEndian)) ||
		(NULL == CU_add_test(pSuite, "test_bytesToLong_bigEndian", test_bytesToLong_bigEndian)) ||
        (NULL == CU_add_test(pSuite, "test_longToBytes_bigEndian", test_longToBytes_bigEndian)) ||
        (NULL == CU_add_test(pSuite, "test_longToBytes_bytesToLong_bigEndian", test_longToBytes_bytesToLong_bigEndian)) ||
        (NULL == CU_add_test(pSuite, "test_getExponent_float", test_getExponent_float)) ||
        (NULL == CU_add_test(pSuite, "test_getPrecisionReqLength_float", test_getPrecisionReqLength_float)) ||
        (NULL == CU_add_test(pSuite, "test_getExponent_double", test_getExponent_double)) ||
        (NULL == CU_add_test(pSuite, "test_getPrecisionReqLength_double", test_getPrecisionReqLength_double)) ||
        (NULL == CU_add_test(pSuite, "test_numberOfLeadingZeros_Int", test_numberOfLeadingZeros_Int)) ||
        (NULL == CU_add_test(pSuite, "test_numberOfLeadingZeros_Long", test_numberOfLeadingZeros_Long)) ||        
        (NULL == CU_add_test(pSuite, "test_getLeadingNumbers_Int", test_getLeadingNumbers_Int)) ||
        (NULL == CU_add_test(pSuite, "test_getLeadingNumbers_Long", test_getLeadingNumbers_Long)) ||
        (NULL == CU_add_test(pSuite, "test_bytesToShort", test_bytesToShort)) ||
        (NULL == CU_add_test(pSuite, "test_bytesToInt", test_bytesToInt)) ||
        (NULL == CU_add_test(pSuite, "test_bytesToLong", test_bytesToLong)) ||
        (NULL == CU_add_test(pSuite, "test_floatToBytes_bytesToFloat", test_floatToBytes_bytesToFloat)) ||
        (NULL == CU_add_test(pSuite, "test_doubleToBytes_bytesToDouble", test_doubleToBytes_bytesToDouble))
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
	 unsigned int num_failures = CU_get_number_of_failures();
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
   return num_failures || CU_get_error();
}
