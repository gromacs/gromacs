#ifndef SZ3_FRONTEND_INTERFACE
#define SZ3_FRONTEND_INTERFACE
/**
 * Frontend is the combination of Predictor and Quantizer
 * For compression, it takes the original data as input, and outputs integer values
 * which will be used for lossless compression by the Encoder and Lossless modules.
 */

#include "SZ3/def.hpp"
#include <vector>

namespace SZ3 {


    namespace concepts {

        /**
         * Frontend is the combination of predictor and quantizer
         * @tparam T
         * @tparam N
         */
        template<class T, uint N>
        class FrontendInterface {
        public:

            virtual ~FrontendInterface() = default;

            /**
             * predict the data and quantize the error
             * @param data original input
             * @return quantized prediction error
             */
            virtual std::vector<int> compress(T *data) = 0;

            /**
             * reverse of compress(), reconstruct the data
             * @param quant_inds quantized prediction error
             * @param dec_data place to write the reconstructed data
             * @return same value with dec_data
             */
            virtual T *decompress(std::vector<int> &quant_inds, T *dec_data) = 0;

            /**
             * serialize the frontend and store it to a buffer
             * @param c One large buffer is pre-allocated, and the start location of the serialized frontend in the buffer is indicated by c.
             *          After saving the frontend to the buffer, this function should change c to indicate the next empty location in the buffer
             */
            virtual void save(uchar *&c) = 0;

            /**
             * deserialize the frontend from a buffer
             * @param c start location of the frontend in the buffer
             * @param remaining_length the remaining length of the buffer
             */
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
