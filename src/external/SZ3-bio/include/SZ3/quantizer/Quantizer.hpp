#ifndef _SZ_QUANTIZER_HPP
#define _SZ_QUANTIZER_HPP


namespace SZ3 {
    namespace concepts {

        template<class T>
        class QuantizerInterface {
        public:

            virtual ~QuantizerInterface() = default;

            /**
             * quantize the error (error=data-pred) based on error bound
             * @param data single data point
             * @param pred predicted value for the data point
             * @return quantized error
             */
            virtual int quantize(T data, T pred) = 0;

            /**
             * quantize the error (error=data-pred) based on error bound, and overwrite the data with reconstructed value
             * @param data single data point
             * @param pred predicted value for this data point
             * @return quantized error
             */
            virtual int quantize_and_overwrite(T &data, T pred) = 0;

            /**
             * reconstructed the data point
             * @param pred predicted value for the data point
             * @param quant_index quantized error
             * @return reconstructed value of the data point
             */
            virtual T recover(T pred, int quant_index) = 0;

            /**
             * reset quantizer to initial state
             */
            virtual void clear() = 0;

            /**
             ** serialize the quantizer and store it to a buffer
             * @param c One large buffer is pre-allocated, and the start location of the serialized quantizer in the buffer is indicated by c.
             *          After saving the quantizer to the buffer, this function should change c to indicate the next empty location in the buffer
             */
            virtual void save(uchar *&c) const = 0;

            /**
             * deserialize the quantizer from a buffer
             * @param c start location of the quantizer in the buffer
             * @param remaining_length the remaining length of the buffer
             */
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
