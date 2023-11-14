/**
 *  @file ArithmeticCoding.c
 *  @author Sheng Di, Mark Thomas Nelson
 *  @date April, 2016
 *  @brief Byte Toolkit
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *  (C) The MIT License (MIT), this code was modified from Mark's arithmetic coding code: http://www.drdobbs.com/cpp/data-compression-with-arithmetic-encodin/240169251?pgno=1
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include "SZ3/encoder/ArithmeticEncoder.hpp"
#include "SZ3/utils/FileUtil.hpp"

unsigned char *readByteData(char *srcFilePath, size_t *byteLength, int *status) {
    FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL) {
        printf("Failed to open input file. 1\n");
        return 0;
    }
    fseek(pFile, 0, SEEK_END);
    *byteLength = ftell(pFile);
    fclose(pFile);

    unsigned char *byteBuf = (unsigned char *) malloc((*byteLength) * sizeof(unsigned char)); //sizeof(char)==1

    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL) {
        printf("Failed to open input file. 2\n");
        return 0;
    }
    fread(byteBuf, 1, *byteLength, pFile);
    fclose(pFile);
    *status = 0;
    return byteBuf;
}

void writeByteData(unsigned char *bytes, size_t byteLength, char *tgtFilePath, int *status) {
    FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL) {
        printf("Failed to open input file. 3\n");
        *status = 0;
        return;
    }

    fwrite(bytes, 1, byteLength, pFile); //write outSize bytes
    fclose(pFile);
    *status = 0;
}

using namespace SZ3;

int main(int argc, char *argv[]) {
    int status = 0;
    char inputFile[100];
    size_t byteLen = 0;
    snprintf(inputFile, 100, "%s", argv[1]);
//    unsigned char *bytes = readByteData(inputFile, &byteLen, &status);
//    std::vector<int> codes(byteLen);
    size_t i = 0;
//    for (i = 0; i < byteLen; i++)
//        codes[i] = bytes[i];

    auto read = readfile<int>(argv[1], byteLen);
//    read[0] = 251;
//    byteLen = 4;
    std::vector<int> codes(read.get(), read.get() + byteLen);

    unsigned char *cmprData = (unsigned char *) malloc(sizeof(unsigned char) * byteLen * 1.2);
    unsigned char *ariCoderBytes = cmprData;

    //compression
//    cost_start();
    ArithmeticEncoder<int> encoder;

    encoder.preprocess_encode(codes, 4096);
    encoder.save(ariCoderBytes);
    encoder.encode(codes, ariCoderBytes);
    encoder.postprocess_encode();
//    cost_end();

    size_t totalCmprSize = ariCoderBytes - cmprData;

    char cmprFile[100];
    snprintf(cmprFile, 100, "%s.ari", inputFile);
    writefile(cmprFile, cmprData, totalCmprSize);

    printf("compressed data size is: %zu\n", totalCmprSize);
    printf("compression ratio is: %f\n", 1.0 * byteLen / totalCmprSize);
//    printf("compression time: %f\n", totalCost);
    const uchar *cmprData2;
    cmprData2 = readfile<uchar>(cmprFile, totalCmprSize).get();

    //decompression
//    ArithmeticEncoder<int> decoder;
    encoder.load(cmprData2, totalCmprSize);

//    for (i = 0; i < encoder.ariCoder.numOfRealStates; i++) {
//        if (encoder.ariCoder.cumulative_frequency[i].high != decoder.ariCoder.cumulative_frequency[i].high) {
//            printf("i=%zu, %ld vs. %ld\n", i, encoder.ariCoder.cumulative_frequency[i].high,
//                   decoder.ariCoder.cumulative_frequency[i].high);
//            break;
//        }
//    }
    printf("done checking\n");
//    cost_start();
    auto decData = encoder.decode(cmprData2, byteLen);
    encoder.postprocess_decode();
//    cost_end();

    int same = 1;
    for (i = 0; i < byteLen; i++) {
        if (codes[i] != decData[i]) {
            printf("Error: i = %zu, codes[i] = %d, decData[i] = %d\n", i, codes[i], decData[i]);
            same = 0;
//            break;
        }
    }

    if (same)
        printf("decompression is correct!\n");
    else
        printf("decomipression is wrong!\n");


//    printf("decompression time: %f\n", totalCost);

    free(cmprData);
//    free(bytes);
}



