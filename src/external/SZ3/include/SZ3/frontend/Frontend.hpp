#ifndef SZ3_FRONTEND_INTERFACE
#define SZ3_FRONTEND_INTERFACE
/**
 * Frontend is the combination of Predictor and Quantizer
 * For compression, it takes the original data as input, and outputs integer values
 * which will be used for lossless compression by the Encoder and Lossless modules.
 */

#include "SZ3/def.hpp"
#include <vector>

namespace SZ {


    namespace concepts {

        template<class T, uint N>
        class FrontendInterface {
        public:

            virtual ~FrontendInterface() = default;

            virtual std::vector<int> compress(T *data) = 0;

            virtual T *decompress(std::vector<int> &quant_inds, T *dec_data) = 0;

            virtual void save(uchar *&c) = 0;

            virtual void load(const uchar *&c, size_t &remaining_length) = 0;

            virtual size_t size_est() = 0;

            virtual int get_radius() const = 0;

            virtual size_t get_num_elements() const = 0;

            virtual void print() = 0;

            virtual void clear() = 0;
        };

    }

}

#endif
