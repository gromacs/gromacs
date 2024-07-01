//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_Config_HPP
#define SZ_Config_HPP

#include <iostream>
#include <vector>
#include <numeric>
#include <cstdint>
#include "SZ3/def.hpp"
#include "MemoryUtil.hpp"
#include "SZ3/utils/inih/INIReader.h"

#define SZ_FLOAT 0
#define SZ_DOUBLE 1
#define SZ_UINT8 2
#define SZ_INT8 3
#define SZ_UINT16 4
#define SZ_INT16 5
#define SZ_UINT32 6
#define SZ_INT32 7
#define SZ_UINT64 8
#define SZ_INT64 9

namespace SZ3 {
    enum EB {
        EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL
    };
    
    constexpr const char *EB_STR[] = {"ABS", "REL", "PSNR", "NORM", "ABS_AND_REL", "ABS_OR_REL"};
    constexpr EB EB_OPTIONS[] = {EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL};
    
    enum ALGO {
        ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP, ALGO_BIOMD, ALGO_BIOMDXTC
    };
    constexpr const char *ALGO_STR[] = {
        "ALGO_LORENZO_REG", "ALGO_INTERP_LORENZO", "ALGO_INTERP", "ALGO_BIOMD", "ALGO_BIOMDXTC"
    };
    constexpr const ALGO ALGO_OPTIONS[] = {
        ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP, ALGO_BIOMD, ALGO_BIOMDXTC
    };
    
    enum INTERP_ALGO {
        INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC
    };
    
    constexpr const char *INTERP_ALGO_STR[] = {"INTERP_ALGO_LINEAR", "INTERP_ALGO_CUBIC"};
    constexpr INTERP_ALGO INTERP_ALGO_OPTIONS[] = {INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC};
    
    template<class T>
    const char *enum2Str(T e) {
        if (std::is_same<T, ALGO>::value) {
            return ALGO_STR[e];
        } else if (std::is_same<T, INTERP_ALGO>::value) {
            return INTERP_ALGO_STR[e];
        } else if (std::is_same<T, EB>::value) {
            return EB_STR[e];
        } else {
            printf("invalid enum type for enum2Str()\n ");
            exit(0);
        }
    }
    
    class Config {
     public:
        template<class... Dims>
        Config(Dims... args) {
            dims = std::vector<size_t>{static_cast<size_t>(std::forward<Dims>(args))...};
            setDims(dims.begin(), dims.end());
        }
        
        template<class Iter>
        size_t setDims(Iter begin, Iter end) {
            auto dims_ = std::vector<size_t>(begin, end);
            dims.clear();
            for (auto dim : dims_) {
                if (dim > 1) {
                    dims.push_back(dim);
                }
            }
            if (dims.empty()) {
                dims = {1};
            }
            N = dims.size();
            num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());
            pred_dim = N;
            blockSize = (N == 1 ? 128 : (N == 2 ? 16 : 6));
            stride = blockSize;
            return num;
        }
        
        void loadcfg(const std::string &cfgpath) {
            INIReader cfg(cfgpath);
            
            if (cfg.ParseError() != 0) {
                std::cout << "Can't load cfg file " << cfgpath << std::endl;
                exit(0);
            }
            
            auto cmprAlgoStr = cfg.Get("GlobalSettings", "CmprAlgo", "");
            if (cmprAlgoStr == ALGO_STR[ALGO_LORENZO_REG]) {
                cmprAlgo = ALGO_LORENZO_REG;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_INTERP_LORENZO]) {
                cmprAlgo = ALGO_INTERP_LORENZO;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_INTERP]) {
                cmprAlgo = ALGO_INTERP;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_BIOMD]) {
                cmprAlgo = ALGO_BIOMD;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_BIOMDXTC]) {
                cmprAlgo = ALGO_BIOMDXTC;
            }
            auto ebModeStr = cfg.Get("GlobalSettings", "ErrorBoundMode", "");
            if (ebModeStr == EB_STR[EB_ABS]) {
                errorBoundMode = EB_ABS;
            } else if (ebModeStr == EB_STR[EB_REL]) {
                errorBoundMode = EB_REL;
            } else if (ebModeStr == EB_STR[EB_PSNR]) {
                errorBoundMode = EB_PSNR;
            } else if (ebModeStr == EB_STR[EB_L2NORM]) {
                errorBoundMode = EB_L2NORM;
            } else if (ebModeStr == EB_STR[EB_ABS_AND_REL]) {
                errorBoundMode = EB_ABS_AND_REL;
            } else if (ebModeStr == EB_STR[EB_ABS_OR_REL]) {
                errorBoundMode = EB_ABS_OR_REL;
            }
            absErrorBound = cfg.GetReal("GlobalSettings", "AbsErrorBound", absErrorBound);
            relErrorBound = cfg.GetReal("GlobalSettings", "RelErrorBound", relErrorBound);
            psnrErrorBound = cfg.GetReal("GlobalSettings", "PSNRErrorBound", psnrErrorBound);
            l2normErrorBound = cfg.GetReal("GlobalSettings", "L2NormErrorBound", l2normErrorBound);
            
            openmp = cfg.GetBoolean("GlobalSettings", "OpenMP", openmp);
            lorenzo = cfg.GetBoolean("AlgoSettings", "Lorenzo", lorenzo);
            lorenzo2 = cfg.GetBoolean("AlgoSettings", "Lorenzo2ndOrder", lorenzo2);
            regression = cfg.GetBoolean("AlgoSettings", "Regression", regression);
            regression2 = cfg.GetBoolean("AlgoSettings", "Regression2ndOrder", regression2);
            
            auto interpAlgoStr = cfg.Get("AlgoSettings", "InterpolationAlgo", "");
            if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_LINEAR]) {
                interpAlgo = INTERP_ALGO_LINEAR;
            } else if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_CUBIC]) {
                interpAlgo = INTERP_ALGO_CUBIC;
            }
            interpDirection = cfg.GetInteger("AlgoSettings", "InterpolationDirection", interpDirection);
            interpBlockSize = cfg.GetInteger("AlgoSettings", "InterpolationBlockSize", interpBlockSize);
            blockSize = cfg.GetInteger("AlgoSettings", "BlockSize", blockSize);
            quantbinCnt = cfg.GetInteger("AlgoSettings", "QuantizationBinTotal", quantbinCnt);
        }
        
        void save(unsigned char *&c) {
            write(N, c);
            write(dims.data(), dims.size(), c);
            write(num, c);
            write(cmprAlgo, c);
            write(errorBoundMode, c);
            write(absErrorBound, c);
            write(relErrorBound, c);
            write(lorenzo, c);
            write(lorenzo2, c);
            write(regression, c);
            write(regression2, c);
            write(interpAlgo, c);
            write(interpDirection, c);
            write(interpBlockSize, c);
            write(lossless, c);
            write(encoder, c);
            write(quantbinCnt, c);
            write(blockSize, c);
            write(stride, c);
            write(pred_dim, c);
            write(openmp, c);
            write(dataType, c);
        };
        
        void load(const unsigned char *&c) {
            read(N, c);
            dims.resize(N);
            read(dims.data(), N, c);
            read(num, c);
            read(cmprAlgo, c);
            read(errorBoundMode, c);
            read(absErrorBound, c);
            read(relErrorBound, c);
            read(lorenzo, c);
            read(lorenzo2, c);
            read(regression, c);
            read(regression2, c);
            read(interpAlgo, c);
            read(interpDirection, c);
            read(interpBlockSize, c);
            read(lossless, c);
            read(encoder, c);
            read(quantbinCnt, c);
            read(blockSize, c);
            read(stride, c);
            read(pred_dim, c);
            read(openmp, c);
            read(dataType, c);
        }
        
        void print() {
            printf("===================== Begin SZ3 Configuration =====================\n");
            printf("N = %d\n", N);
            printf("dims = ");
            for (auto dim : dims) {
                printf("%zu ", dim);
            }
            printf("\nnum = %zu\n", num);
            printf("CmprAlgo = %s\n", enum2Str((ALGO) cmprAlgo));
            printf("ErrorBoundMode = %s\n", enum2Str((EB) errorBoundMode));
            printf("AbsErrorBound = %f\n", absErrorBound);
            printf("RelErrorBound = %f\n", relErrorBound);
            printf("PSNRErrorBound = %f\n", psnrErrorBound);
            printf("L2NormErrorBound = %f\n", l2normErrorBound);
            printf("Lorenzo = %d\n", lorenzo);
            printf("Lorenzo2ndOrder = %d\n", lorenzo2);
            printf("Regression = %d\n", regression);
            printf("Regression2ndOrder = %d\n", regression2);
            printf("OpenMP = %d\n", openmp);
            printf("DataType = %d\n", dataType);
            printf("Lossless = %d\n", lossless);
            printf("Encoder = %d\n", encoder);
            printf("InterpolationAlgo = %s\n", enum2Str((INTERP_ALGO) interpAlgo));
            printf("InterpolationDirection = %d\n", interpDirection);
            printf("InterpolationBlockSize = %d\n", interpBlockSize);
            printf("QuantizationBinTotal = %d\n", quantbinCnt);
            printf("BlockSize = %d\n", blockSize);
            printf("Stride = %d\n", stride);
            printf("PredDim = %d\n", pred_dim);
            printf("===================== End SZ3 Configuration =====================\n");
        }
        
        static size_t size_est() {
            return sizeof(Config) + sizeof(size_t) * 5 + 32; //sizeof(size_t) * 5 is for dims vector, 32 is for redundancy
            return sizeof(size_t) * 5 + sizeof(double) * 4 + sizeof(bool) * 5 + sizeof(uint8_t) * 7 + sizeof(int) * 5 +
                50; //50 is for redundancy
        }
        
        char N;
        std::vector<size_t> dims;
        size_t num;
        uint8_t cmprAlgo = ALGO_INTERP_LORENZO;
        uint8_t errorBoundMode = EB_ABS;
        double absErrorBound = 1e-3;
        double relErrorBound = 0;
        double psnrErrorBound = 0;
        double l2normErrorBound = 0;
        bool lorenzo = true;
        bool lorenzo2 = false;
        bool regression = true;
        bool regression2 = false;
        bool openmp = false;
        uint8_t dataType = SZ_FLOAT; // dataType is only used in HDF5 filter
        uint8_t lossless = 1; // 0-> skip lossless(use lossless_bypass); 1-> zstd
        uint8_t encoder = 1; // 0-> skip encoder; 1->HuffmanEncoder; 2->ArithmeticEncoder
        uint8_t interpAlgo = INTERP_ALGO_CUBIC;
        uint8_t interpDirection = 0;
        int interpBlockSize = 32;
        int quantbinCnt = 65536;
        int blockSize;
        int stride; //not used now
        int pred_dim; // not used now
    };
}

#endif //SZ_CONFIG_HPP
