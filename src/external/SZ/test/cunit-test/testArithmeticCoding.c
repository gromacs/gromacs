#include <stdio.h>
#include <stdlib.h>
#include <sz.h>

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start()
{
        totalCost = 0;
        gettimeofday(&costStart, NULL);
}

void cost_end()
{
        double elapsed;
        struct timeval costEnd;
        gettimeofday(&costEnd, NULL);
        elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
        totalCost += elapsed;
}

void main(int argc, char* argv[])
{
	int status = 0;
	char inputFile[100];
	size_t byteLen = 0;
	sprintf(inputFile, "%s", argv[1]);
	unsigned char* bytes = readByteData(inputFile, &byteLen, &status);
	int* codes = (int*)malloc(sizeof(int)*byteLen);
	size_t i = 0;
	for(i=0;i<byteLen;i++)
		codes[i] = bytes[i];
	unsigned char* cmprData = (unsigned char*)malloc(sizeof(unsigned char)*byteLen*1.2);
	unsigned char* ariCoderBytes = NULL;
	
	//compression
	cost_start();
	AriCoder *ariCoder = createAriCoder(256, codes, byteLen);
	unsigned int padSize = pad_ariCoder(ariCoder, &ariCoderBytes);
	memcpy(cmprData, ariCoderBytes, padSize);
	unsigned char* p = cmprData + padSize;
	size_t encodeSize = 0;
	ari_encode(ariCoder, codes, byteLen, p, &encodeSize);
	cost_end();

	size_t totalCmprSize = padSize + encodeSize;
	
	char cmprFile[100];
	sprintf(cmprFile, "%s.ari", inputFile);
	writeByteData(cmprData, totalCmprSize, cmprFile, &status);

	printf("compressed data size is: %zu\n", totalCmprSize);
	printf("compression ratio is: %f\n", 1.0*byteLen/totalCmprSize);
	printf("compression time: %f\n", totalCost);

	//decompression
	cost_start();
	AriCoder *ariCoder2 = NULL;
	int offset = unpad_ariCoder(&ariCoder2, cmprData);

	for(i=0;i<ariCoder->numOfRealStates;i++)
	{
		if(ariCoder->cumulative_frequency[i].high != ariCoder2->cumulative_frequency[i].high)
		{
			printf("i=%zu, %d vs. %d\n", i, ariCoder->cumulative_frequency[i].high, ariCoder2->cumulative_frequency[i].high);
			break;
		}
	}
	printf("done checking\n");

	int* decData = (int*)malloc(sizeof(int)*byteLen);
	ari_decode(ariCoder2, cmprData+offset, totalCmprSize-offset, byteLen, decData);
	freeAriCoder(ariCoder2);
	cost_end();

	int same = 1;
	for(i=0;i<byteLen;i++)
	{
		if(codes[i]!=decData[i])
		{
			printf("Error: i = %zu, codes[i] = %d, decData[i] = %d\n", i, codes[i], decData[i]);
			same = 0;
			break;
		}
	}

	if(same)
		printf("decompression is correct!\n");
	else
		printf("decomipression is wrong!\n");


	printf("decompression time: %f\n", totalCost);

	freeAriCoder(ariCoder);
	free(ariCoderBytes);
	free(codes);
	free(cmprData);
	free(bytes);
	free(decData);
}
