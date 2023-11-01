//
// Created by Kai Zhao on 10/27/22.
//

#include "sz3c.h"
#include "SZ3/api/sz.hpp"


using namespace SZ3;

unsigned char *SZ_compress_args(int dataType, void *data, size_t *outSize,
                                int errBoundMode, double absErrBound, double relBoundRatio, double pwrBoundRatio,
                                size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {

    SZ3::Config conf;
    if (r2 == 0) {
        conf = SZ3::Config(r1);
    } else if (r3 == 0) {
        conf = SZ3::Config(r2, r1);
    } else if (r4 == 0) {
        conf = SZ3::Config(r3, r2, r1);
    } else if (r5 == 0) {
        conf = SZ3::Config(r4, r3, r2, r1);
    } else {
        conf = SZ3::Config(r5 * r4, r3, r2, r1);
    }
//    conf.loadcfg(conPath);
    conf.absErrorBound = absErrBound;
    conf.relErrorBound = relBoundRatio;
//    conf.pwrErrorBound = pwrBoundRatio;
    if (errBoundMode == ABS) {
        conf.errorBoundMode = EB_ABS;
    } else if (errBoundMode == REL) {
        conf.errorBoundMode = EB_REL;
    } else if (errBoundMode == ABS_AND_REL) {
        conf.errorBoundMode = EB_ABS_AND_REL;
    } else if (errBoundMode == ABS_OR_REL) {
        conf.errorBoundMode = EB_ABS_OR_REL;
    } else {
        printf("errBoundMode %d not support\n ", errBoundMode);
        exit(0);
    }

    unsigned char *cmpr_data = NULL;
    if (dataType == SZ_FLOAT) {
        cmpr_data = (unsigned char *) SZ_compress<float>(conf, (float *) data, *outSize);
    } else if (dataType == SZ_DOUBLE) {
        cmpr_data = (unsigned char *) SZ_compress<double>(conf, (double *) data, *outSize);
    } else {
        printf("dataType %d not support\n", dataType);
        exit(0);
    }

    //convert c++ memory (by 'new' operator) to c memory (by malloc)
    auto *cmpr = (unsigned char *) malloc(*outSize);
    memcpy(cmpr, cmpr_data, *outSize);
    delete[]cmpr_data;

    return cmpr;

}

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength,
                    size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {
    size_t n = 0;
    if (r2 == 0) {
        n = r1;
    } else if (r3 == 0) {
        n = r1 * r2;
    } else if (r4 == 0) {
        n = r1 * r2 * r3;
    } else if (r5 == 0) {
        n = r1 * r2 * r3 * r4;
    } else {
        n = r1 * r2 * r3 * r4 * r5;
    }

    SZ3::Config conf;
    if (dataType == SZ_FLOAT) {
        auto dec_data = (float *) malloc(n * sizeof(float));
        SZ_decompress<float>(conf, (char *) bytes, byteLength, dec_data);
        return dec_data;
    } else if (dataType == SZ_DOUBLE) {
        auto dec_data = (double *) malloc(n * sizeof(double));
        SZ_decompress<double>(conf, (char *) bytes, byteLength, dec_data);
        return dec_data;
    } else {
        printf("dataType %d not support\n", dataType);
        exit(0);
    }
}