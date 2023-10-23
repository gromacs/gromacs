/**
 *
 * A set of python bindings for SZ 
 * 
 * Developed by Robert Underwood while he was at Clemson University
 * This material is based upon work supported by the National Science 
 * Foundation under Grant No. 1633608.
 * 
 * Copyright © 2019 Robert Underwood
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY Robert Underwood ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL Robert Underwood BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation
 * are those of the authors and should not be interpreted as representing
 * official policies, either expressed or implied, of Robert Underwood.
 */
#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include "sz.h"
#include "exafelSZ.h"
#include "pysz_private.h"


struct Config {
  sz_params params;
  std::string app;
};


class ConfigBuilder {
  public:
  ConfigBuilder();
  ConfigBuilder(std::string const& config_file_path);
  Config build();

	ConfigBuilder& app(std::string app) noexcept;
	ConfigBuilder& absErrBound(double value) noexcept;
	ConfigBuilder& dataType(int value) noexcept;
	ConfigBuilder& errorBoundMode(int  value) noexcept;
	ConfigBuilder& gzipMode(int value) noexcept;
	ConfigBuilder& losslessCompressor(int value) noexcept;
	ConfigBuilder& maxRangeRadius(unsigned int value) noexcept;
	ConfigBuilder& max_quant_intervals(unsigned int value) noexcept;
	ConfigBuilder& predThreshold(float value) noexcept;
	ConfigBuilder& predictionMode(int value) noexcept;
	ConfigBuilder& psnr(double value) noexcept;
	ConfigBuilder& pw_relBoundRatio(double value) noexcept;
	ConfigBuilder& pwr_type(int value) noexcept;
	ConfigBuilder& quantization_intervals(unsigned int value) noexcept;
	ConfigBuilder& randomAccess(int value) noexcept;
	ConfigBuilder& relBoundRatio(double value) noexcept;
	ConfigBuilder& sampleDistance(int value) noexcept;
	ConfigBuilder& segment_size(int value) noexcept;
	ConfigBuilder& snapshotCmprStep(int value) noexcept;
	ConfigBuilder& sol_ID(int value) noexcept;
	ConfigBuilder& szMode(int value) noexcept;
  ConfigBuilder& accelerate_pw_rel_compression(int value) noexcept;
  ConfigBuilder& plus_bits(int value) noexcept;
  private:
  Config building;
};


class ExaFELConfigBuilder {
  public:
  exafelSZ_params* build() { 
    exafelSZ_params* params = new exafelSZ_params;
    params->peaksSegs = new uint16_t [peaks.size()];
    params->peaksRows = new uint16_t [peaks.size()];
    params->peaksCols = new uint16_t [peaks.size()];
    params->calibPanel = new uint8_t [calibPanel.size()];
    std::copy(std::begin(peaks), std::end(peaks), params->peaksSegs);
    std::copy(std::begin(peakRows), std::end(peakRows), params->peaksRows);
    std::copy(std::begin(peakCols), std::end(peakCols), params->peaksCols);
    std::copy(std::begin(calibPanel), std::end(calibPanel), params->calibPanel);
    params->binSize = binSize;
    params->tolerance = tolerance;
    params->szDim = szDim;
    params->peakSize = peakSize;
    return params;
  }

  int binSize; //Binning: (pr->binSize x pr->binSize) to (1 x 1)
  double tolerance; //SZ pr->tolerance
  int szDim; //1D/2D/3D compression/decompression
  int peakSize; //MUST BE ODD AND NOT EVEN! Each peak will have size of: (peakSize x peakSize)


  static void free(exafelSZ_params* params) {
    delete[] params->peaksSegs;
    delete[] params->peaksRows;
    delete[] params->peaksCols;
    delete[] params->calibPanel;
    delete params;
  }

  void setPeaks(std::vector<uint16_t> const& peaks) {
    this->peaks = peaks;
  }
  void setPeaksRows(std::vector<uint16_t> const& rows) {
    this->peakRows = rows;
  }
  void setPeaksCols(std::vector<uint16_t> const& cols) {
    this->peakCols = cols;
  }
  void setCalibPanel(std::vector<uint8_t> const& calibPanel) {
    this->calibPanel = calibPanel;
  }

  private:
  std::vector<uint16_t> peaks;
  std::vector<uint16_t> peakRows;
  std::vector<uint16_t> peakCols;
  std::vector<uint8_t> calibPanel;
};

class Compressor {
  public:
  /*
   * Curently SZ has global state that needs to initialized only once; Thus
   * compressor should be a move only type for now and use RAII to manage the
   * global state
   */
  Compressor(Config config);
  ~Compressor();
  
  Compressor(Compressor const&)=delete;
  Compressor& operator=(Compressor&)=delete;

  //
  //Yes, this could be done using varatic arguments, but due to limitations in how SWIG handles templates
  //We need to explicitly list all of these out.
  //
  //even though sz supports 5d arrays, they aren't supported well by numpy, omit it for now
  //

  template<class T>
  std::string Compress1(T* data, size_t r1, void* params=nullptr) {
    size_t outsize = 0;
    int status;
    char* tmp = reinterpret_cast<char*>(SZ_compress_customize(app.c_str(), params, SZTypeToTypeID<T>::value, data, 0, 0, 0, 0, r1, &outsize, &status));
    return std::string(tmp, outsize);
  }
  template<class T>
  std::string Compress2(T* data, size_t r1, size_t r2, void* params=nullptr) {
    int status;
    size_t outsize = 0;
    auto tmp = reinterpret_cast<char*>(SZ_compress_customize(app.c_str(), params, SZTypeToTypeID<T>::value, data, 0, 0, 0, r2, r1, &outsize, &status));
    return std::string(tmp, outsize);
  }
  template<class T>
  std::string Compress3(T* data, size_t r1, size_t r2, size_t r3, void* params=nullptr) {
    int status;
    size_t outsize = 0;
    auto tmp = reinterpret_cast<char*>(SZ_compress_customize(app.c_str(), params, SZTypeToTypeID<T>::value, data, 0, 0, r3, r2, r1, &outsize, &status));
    return std::string(tmp, outsize);
  }
  template<class T>
  std::string Compress4(T* data, size_t r1, size_t r2, size_t r3, size_t r4, void* params=nullptr) {
    int status;
    size_t outsize = 0;
    auto tmp = reinterpret_cast<char*>(SZ_compress_customize(app.c_str(), params, SZTypeToTypeID<T>::value, data, 0, r4, r3, r2, r1, &outsize, &status));
    return std::string(tmp, outsize);

  }

  template<class T>
  std::vector<T> Decompress(std::string data, std::vector<int> r, void* params=nullptr) {
    T* decompressed = nullptr;
    size_t len = 0;
    int status;
    switch(r.size()) {
      case 1:

        decompressed = (T*)SZ_decompress_customize(app.c_str(), params, SZTypeToTypeID<T>::value, (unsigned char*)(data.c_str()), data.length(), 0,0,0,0,r[0], &status);
        len = r[0];
        break;
      case 2:
        decompressed = (T*)SZ_decompress_customize(app.c_str(), params, SZTypeToTypeID<T>::value, (unsigned char*)(data.c_str()), data.length(), 0,0,0,r[1],r[0], &status);
        len = r[0]*r[1];
        break;
      case 3:
        decompressed = (T*)SZ_decompress_customize(app.c_str(), params, SZTypeToTypeID<T>::value, (unsigned char*)(data.c_str()), data.length(), 0,0,r[2],r[1],r[0], &status);
        len = r[0]*r[1]*r[2];
        break;
      case 4:
        decompressed = (T*)SZ_decompress_customize(app.c_str(), params, SZTypeToTypeID<T>::value, (unsigned char*)(data.c_str()), data.length(), 0,r[3],r[2],r[1],r[0], &status);
        len = r[0]*r[1]*r[2]*r[3];
        break;
      default:
        printf("%zu dimensional arrays not supported\n", r.size());
    }
    return std::vector<T>(decompressed, decompressed+len);
  }

  private:
  std::string app;
};
