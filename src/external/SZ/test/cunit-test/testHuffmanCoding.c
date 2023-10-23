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
	
	cost_start();
	unsigned char* out = NULL;
	size_t outSize = 0;
	HuffmanTree* huffmanTree = createHuffmanTree(256);
	encode_withTree(huffmanTree, codes, byteLen, &out, &outSize);
	SZ_ReleaseHuffman(huffmanTree);	
	cost_end();

	size_t totalCmprSize = outSize;
	
	printf("compressed data size is: %zu\n", totalCmprSize);
	printf("compression ratio is: %f\n", 1.0*byteLen/totalCmprSize);
	printf("compression time: %f\n", totalCost);
	free(codes);
	free(out);
	free(bytes);
}
