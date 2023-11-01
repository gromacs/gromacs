#ifndef _SZ_ENCODER_HPP
#define _SZ_ENCODER_HPP

#include <vector>
namespace SZ3 {
    namespace concepts {

        template<class T>
        class EncoderInterface {
        public:

            virtual ~EncoderInterface() = default;

            /**
             * init the encoder
             * E.g., Huffman will build tree in this step
             * @param bins to-be-encoded integers
             * @param stateNum stateNum > 0 indicates the bins has a range of [0, stateNum). stateNum == 0 means no such guarantee
             */
            virtual void preprocess_encode(const std::vector<T> &bins, int stateNum) = 0;

            /**
             * encode the input (in vector<T> format) to a more compact representative(in byte stream format)
             * @param bins input in vector
             * @param bytes output in byte stream
             * @return size of output (# of bytes)
             */
            virtual size_t encode(const std::vector<T> &bins, uchar *&bytes) = 0;

            virtual void postprocess_encode() = 0;

            virtual void preprocess_decode() = 0;

            /**
             * reverse of encode()
             * @param bytes input in byte stream
             * @param targetLength size of the output vector
             * @return output in vector
             */
            virtual std::vector<T> decode(const uchar *&bytes, size_t targetLength) = 0;

            virtual void postprocess_decode() = 0;

            /**
             * serialize the encoder and store it to a buffer
             * @param c One large buffer is pre-allocated, and the start location of the serialized encoder in the buffer is indicated by c.
             *          After saving the encoder to the buffer, this function should change c to indicate the next empty location in the buffer
             */
            virtual void save(uchar *&c) = 0;

            /**
             * deserialize the encoder from a buffer
             * @param c start location of the encoder in the buffer
             * @param remaining_length the remaining length of the buffer
             */
            virtual void load(const uchar *&c, size_t &remaining_length) = 0;

            // return the size of the encoder itself (such as the tree size of the huffman encoder)
            virtual size_t size_est() {
                return 0;
            }

        };
    }
}
#endif
