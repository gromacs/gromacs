#ifndef _SZ_ENCODER_HPP
#define _SZ_ENCODER_HPP

namespace SZ {
    namespace concepts {

        template<class T>
        class EncoderInterface {
        public:

            virtual ~EncoderInterface() = default;

            virtual void preprocess_encode(const std::vector<T> &bins, int stateNum) = 0;

            virtual size_t encode(const std::vector<T> &bins, uchar *&bytes) = 0;

            virtual void postprocess_encode() = 0;

            virtual void preprocess_decode() = 0;

            virtual std::vector<T> decode(const uchar *&bytes, size_t targetLength) = 0;

            virtual void postprocess_decode() = 0;

            virtual uint save(uchar *&c) = 0;

            virtual void load(const uchar *&c, size_t &remaining_length) = 0;

            // return the size of the encoder itself (such as the tree size of the huffman encoder)
            virtual size_t size_est() {
                return 0;
            }

        };
    }
}
#endif
