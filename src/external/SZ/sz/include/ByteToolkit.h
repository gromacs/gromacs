/**
 *  @file ByteToolkit.h
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief Header file for the ByteToolkit.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _ByteToolkit_H
#define _ByteToolkit_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

//ByteToolkit.c

extern unsigned short bytesToUInt16_bigEndian(unsigned char* bytes);
extern unsigned int bytesToUInt32_bigEndian(unsigned char* bytes);
extern unsigned long bytesToUInt64_bigEndian(unsigned char* b);

extern short bytesToInt16_bigEndian(unsigned char* bytes);
extern int bytesToInt32_bigEndian(unsigned char* bytes);
extern long bytesToInt64_bigEndian(unsigned char* b);
extern int bytesToInt_bigEndian(unsigned char* bytes);

extern void intToBytes_bigEndian(unsigned char *b, unsigned int num);

extern void int64ToBytes_bigEndian(unsigned char *b, uint64_t num);
extern void int32ToBytes_bigEndian(unsigned char *b, uint32_t num);
extern void int16ToBytes_bigEndian(unsigned char *b, uint16_t num);

extern long bytesToLong_bigEndian(unsigned char* b);
extern void longToBytes_bigEndian(unsigned char *b, unsigned long num);
long doubleToOSEndianLong(double value);
int floatToOSEndianInt(float value);
extern short getExponent_float(float value);
extern short getPrecisionReqLength_float(float precision);
extern short getExponent_double(double value);
extern short getPrecisionReqLength_double(double precision);
unsigned char numberOfLeadingZeros_Int(int i);
unsigned char numberOfLeadingZeros_Long(long i);
unsigned char getLeadingNumbers_Int(int v1, int v2);
unsigned char getLeadingNumbers_Long(long v1, long v2);
short bytesToShort(unsigned char* bytes);
void shortToBytes(unsigned char* b, short value);
int bytesToInt(unsigned char* bytes);
long bytesToLong(unsigned char* bytes);
extern float bytesToFloat(unsigned char* bytes);
extern void floatToBytes(unsigned char *b, float num);
extern double bytesToDouble(unsigned char* bytes);
extern void doubleToBytes(unsigned char *b, double num);
int extractBytes(unsigned char* byteArray, size_t k, int validLength);
int getMaskRightCode(int m);
extern int getLeftMovingCode(int kMod8);
extern int getRightMovingSteps(int kMod8, int resiBitLength);
extern int getRightMovingCode(int kMod8, int resiBitLength);
short* convertByteDataToShortArray(unsigned char* bytes, size_t byteLength);
unsigned short* convertByteDataToUShortArray(unsigned char* bytes, size_t byteLength);

void convertShortArrayToBytes(short* states, size_t stateLength, unsigned char* bytes);
void convertUShortArrayToBytes(unsigned short* states, size_t stateLength, unsigned char* bytes);
void convertIntArrayToBytes(int* states, size_t stateLength, unsigned char* bytes);
void convertUIntArrayToBytes(unsigned int* states, size_t stateLength, unsigned char* bytes);
void convertLongArrayToBytes(int64_t* states, size_t stateLength, unsigned char* bytes);
void convertULongArrayToBytes(uint64_t* states, size_t stateLength, unsigned char* bytes);

extern size_t bytesToSize(unsigned char* bytes);
extern void sizeToBytes(unsigned char* outBytes, size_t size);

void put_codes_to_output(unsigned int buf, int bitSize, unsigned char** p, int* lackBits, size_t *outSize);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _ByteToolkit_H  ----- */

