#ifndef _SZ_QUANTIZER_HPP
#define _SZ_QUANTIZER_HPP


namespace SZ {
    namespace concepts {

        template<class T>
        class QuantizerInterface {
        public:

            virtual ~QuantizerInterface() = default;

            virtual int quantize(T data, T pred) = 0;

            virtual int quantize_and_overwrite(T &data, T pred) = 0;

            virtual T recover(T pred, int quant_index) = 0;

            /**
             * reset quantizer to initial state
             */
            virtual void clear() = 0;

            virtual void save(uchar *&c) const = 0;

            virtual void load(const uchar *&c, size_t &remaining_length) = 0;

            virtual void precompress_data() = 0;

            virtual void predecompress_data() = 0;

            /**
             * this function is always executed before save()
             * DO NOT reset non-temporary variables (such as unpredictable_data) in this function.
             */
            virtual void postcompress_data() = 0;

            /**
            * DO NOT reset non-temporary variables (such as unpredictable_data) in this function.
            */
            virtual void postdecompress_data() = 0;

        };
    }
}

#endif
