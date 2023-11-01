#ifndef SZ_ArithmeticEncoder_HPP
#define SZ_ArithmeticEncoder_HPP

#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include <cassert>
#include <iostream>

namespace SZ3 {
#define ONE_FOURTH (0x40000000000) //44 bits are absolutely enough to deal with a large dataset (support at most 16TB per process)
#define ONE_HALF (0x80000000000)
#define THREE_FOURTHS (0xC0000000000)
#define MAX_CODE (0xFFFFFFFFFFF)
#define MAX_INTERVALS 1048576 //the limit to the arithmetic coding (at most 2^(20) intervals)


    template<class T>
    class ArithmeticEncoder : public concepts::EncoderInterface<T> {
    public:
        struct Prob {
            size_t low;
            size_t high;
            int state;
        };

        struct AriCoder {
            int numOfRealStates; //the # real states menas the number of states after the optimization of # intervals
            int numOfValidStates; //the # valid states means the number of non-zero frequency cells (some states/codes actually didn't appear)
            size_t total_frequency;
            Prob *cumulative_frequency; //used to encode data more efficiencly
        };


//        ArithmeticEncoder(int stateNum, bool transform = false) {
        ArithmeticEncoder(bool transform = false) {
//            assert(stateNum <= 4096 && "StateNum of Arithmetic Encoder should be <= 4096");
//            ariCoder.numOfRealStates = stateNum;
            ariCoder.numOfRealStates = 0;
            ariCoder.numOfValidStates = 0;
            ariCoder.total_frequency = 0;
            ariCoder.cumulative_frequency = NULL;
            this->transform = transform;
        }

        ~ArithmeticEncoder() {
            if (ariCoder.cumulative_frequency != NULL) {
                free(ariCoder.cumulative_frequency);
            }
        }

        void postprocess_encode() {
            if (transform) {
                bins_transform.clear();
            }
        };

        void preprocess_decode() {};

        void postprocess_decode() {};


        void preprocess_encode(const std::vector<T> &bins, int stateNum) {
            assert(stateNum <= 4096 && "StateNum of Arithmetic Encoder should be <= 4096");
            ariCoder.numOfRealStates = stateNum;
            const T *s = bins.data();
            size_t length = bins.size();
            if (transform) {
                bins_transform = bins;
                for (size_t i = 0; i < bins_transform.size(); i++) {
                    T x = bins_transform[i];
                    bins_transform[i] = fabs(x - ariCoder.numOfRealStates / 2) * 2;
                    if (x - ariCoder.numOfRealStates / 2 < 0) {
                        bins_transform[i] -= 1;
                    }
                }

                s = bins_transform.data();
            }


            size_t i; //# states is in the range of integer.
            int index = 0;
            size_t *freq = (size_t *) malloc(ariCoder.numOfRealStates * sizeof(size_t));
            memset(freq, 0, ariCoder.numOfRealStates * sizeof(size_t));
            for (i = 0; i < length; i++) {
                index = s[i];
                freq[index]++;
            }

            int counter = 0;
            size_t _sum = 0, sum = 0, freqDiv = 0;
            ariCoder.cumulative_frequency = (Prob *) malloc(ariCoder.numOfRealStates * sizeof(Prob));

            memset(ariCoder.cumulative_frequency, 0, ariCoder.numOfRealStates * sizeof(Prob));

            if (length <= MAX_INTERVALS) {
                for (index = 0; index < ariCoder.numOfRealStates; index++) {
                    if (freq[index]) {
                        sum += freq[index];
                        (ariCoder.cumulative_frequency[index]).low = _sum;
                        (ariCoder.cumulative_frequency[index]).high = sum;
                        (ariCoder.cumulative_frequency[index]).state = index;
                        _sum = sum;
                        counter++;
                    }
                }
                ariCoder.numOfValidStates = counter;
                ariCoder.total_frequency = sum;
            } else {
                int intvSize = length % MAX_INTERVALS == 0 ? length / MAX_INTERVALS : length / MAX_INTERVALS + 1;
                for (index = 0; index < ariCoder.numOfRealStates; index++) {
                    if (freq[index]) {
                        freqDiv = freq[index] /
                                  intvSize; //control the sum of frequency to be no greater than MAX_INTERVALS
                        if (freqDiv == 0)
                            freqDiv = 1;
                        sum += freqDiv;
                        (ariCoder.cumulative_frequency[index]).low = _sum;
                        (ariCoder.cumulative_frequency[index]).high = sum;
                        (ariCoder.cumulative_frequency[index]).state = index;
                        _sum = sum;
                        counter++;
                    }
                }
                ariCoder.numOfValidStates = counter;
                ariCoder.total_frequency = sum;
            }

            free(freq);
        }


        uint save(uchar *&p) {

            int numOfRealStates = ariCoder.numOfRealStates;
            int numOfValidStates = ariCoder.numOfValidStates;
            uint64_t total_frequency = ariCoder.total_frequency;
            Prob *cumulative_frequency = ariCoder.cumulative_frequency;

            unsigned int outSize = 0;

            intToBytes_bigEndian(p, numOfRealStates);
            p += sizeof(int);
            intToBytes_bigEndian(p, numOfValidStates);
            p += sizeof(int);
            int64ToBytes_bigEndian(p, total_frequency);
            p += sizeof(uint64_t);
            size_t i = 0;
            if (total_frequency <= 65536) {
                uint16_t low, high;
                if (numOfRealStates <= 256) {
                    for (i = 0; i < numOfRealStates; i++) {
                        high = (uint16_t) (cumulative_frequency[i].high);
                        if (high != 0) //if this state cell is not null
                        {
                            low = (uint16_t) (cumulative_frequency[i].low);
                            int16ToBytes_bigEndian(p, low);
                            p += sizeof(uint16_t);
                            int16ToBytes_bigEndian(p, high);
                            p += sizeof(uint16_t);
                            *(p++) = (unsigned char) cumulative_frequency[i].state;
                            //if(((unsigned char)cumulative_frequency[i].state)==129)
                        }
                    }
                    outSize =
                            2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 5; //2*sizeof(uint16_t)+1
                } else if (numOfRealStates <= 65536) {
                    for (i = 0; i < numOfRealStates; i++) {
                        high = (uint16_t) (cumulative_frequency[i].high);
                        if (high != 0) {
                            low = (uint16_t) (cumulative_frequency[i].low);
                            int16ToBytes_bigEndian(p, low);
                            p += sizeof(uint16_t);
                            int16ToBytes_bigEndian(p, high);
                            p += sizeof(uint16_t);
                            uint16_t state = (uint16_t) cumulative_frequency[i].state;
                            int16ToBytes_bigEndian(p, state);
                            p += sizeof(uint16_t);
                        }
                    }
                    outSize = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 6;
                } else {
                    for (i = 0; i < numOfRealStates; i++) {
                        high = (uint16_t) (cumulative_frequency[i].high);
                        if (high != 0) {
                            low = (uint16_t) (cumulative_frequency[i].low);
                            int16ToBytes_bigEndian(p, low);
                            p += sizeof(uint16_t);
                            int16ToBytes_bigEndian(p, high);
                            p += sizeof(uint16_t);
                            int32ToBytes_bigEndian(p, cumulative_frequency[i].state);
                            p += sizeof(uint32_t);
                        }
                    }
                    outSize = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 8;
                }
            } else if (total_frequency <= 4294967296) {
                uint32_t low, high;
                if (numOfRealStates <= 256) {
                    for (i = 0; i < numOfRealStates; i++) {
                        high = (uint32_t) (cumulative_frequency[i].high);
                        if (high != 0) {
                            low = (uint32_t) (cumulative_frequency[i].low);
                            int32ToBytes_bigEndian(p, low);
                            p += sizeof(uint32_t);
                            int32ToBytes_bigEndian(p, high);
                            p += sizeof(uint32_t);
                            *(p++) = (unsigned char) cumulative_frequency[i].state;
                        }
                    }
                    outSize = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 9;
                } else if (numOfRealStates <= 65536) {
                    for (i = 0; i < numOfRealStates; i++) {
                        high = (uint32_t) (cumulative_frequency[i].high);
                        if (high != 0) {
                            low = (uint32_t) (cumulative_frequency[i].low);
                            int32ToBytes_bigEndian(p, low);
                            p += sizeof(uint32_t);
                            int32ToBytes_bigEndian(p, high);
                            p += sizeof(uint32_t);
                            uint16_t state = (uint16_t) cumulative_frequency[i].state;
                            int16ToBytes_bigEndian(p, state);
                            p += sizeof(uint16_t);

                        }
                    }
                    outSize = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 10;
                } else {
                    for (i = 0; i < numOfRealStates; i++) {
                        high = (uint32_t) (cumulative_frequency[i].high);
                        if (high != 0) {
                            low = (uint32_t) (cumulative_frequency[i].low);
                            int32ToBytes_bigEndian(p, low);
                            p += sizeof(uint32_t);
                            int32ToBytes_bigEndian(p, high);
                            p += sizeof(uint32_t);
                            int32ToBytes_bigEndian(p, cumulative_frequency[i].state);
                            p += sizeof(uint32_t);
                        }
                    }
                    outSize = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 12;
                }
            } else {
                uint64_t low, high;
                if (numOfRealStates <= 256) {
                    for (i = 0; i < numOfRealStates; i++) {
                        high = (uint64_t) (cumulative_frequency[i].high);
                        if (high != 0) {
                            low = (uint64_t) (cumulative_frequency[i].low);
                            int64ToBytes_bigEndian(p, low);
                            p += sizeof(uint64_t);
                            int64ToBytes_bigEndian(p, high);
                            p += sizeof(uint64_t);
                            *(p++) = (unsigned char) cumulative_frequency[i].state;
                        }
                    }
                    outSize = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 17;
                } else if (numOfRealStates <= 65536) {
                    for (i = 0; i < numOfRealStates; i++) {
                        high = (uint64_t) (cumulative_frequency[i].high);
                        if (high != 0) {
                            low = (uint64_t) (cumulative_frequency[i].low);
                            int64ToBytes_bigEndian(p, low);
                            p += sizeof(uint64_t);
                            int64ToBytes_bigEndian(p, high);
                            p += sizeof(uint64_t);
                            uint16_t state = (uint16_t) cumulative_frequency[i].state;
                            int16ToBytes_bigEndian(p, state);
                            p += sizeof(uint16_t);
                        }
                    }
                    outSize = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 18;
                } else {
                    for (i = 0; i < numOfRealStates; i++) {
                        high = (uint64_t) (cumulative_frequency[i].high);
                        if (high != 0) {
                            low = (uint64_t) (cumulative_frequency[i].low);
                            int64ToBytes_bigEndian(p, low);
                            p += sizeof(uint64_t);
                            int64ToBytes_bigEndian(p, high);
                            p += sizeof(uint64_t);
                            int32ToBytes_bigEndian(p, cumulative_frequency[i].state);
                            p += sizeof(uint32_t);
                        }
                    }
                    outSize = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 20;
                }
            }
            return outSize;
        }

/**
 * Reconstruct AriCoder based on the bytes loaded from compressed data
 * @param AriCoder** ariCoder (ourput)
 * @param unsigned char* bytes (input)
 *
 * @return offset
 * */
        void load(const uchar *&p, size_t &remaining_length) {

//        int unpad_ariCoder(AriCoder **ariCoder, unsigned char *bytes) {
            int offset = 0;

            int numOfRealStates = ariCoder.numOfRealStates = bytesToInt_bigEndian(p);
            p += sizeof(int);
            int numOfValidStates = ariCoder.numOfValidStates = bytesToInt_bigEndian(p);
            p += sizeof(int);
            size_t total_frequency = ariCoder.total_frequency = bytesToInt64_bigEndian(p);
            p += sizeof(uint64_t);

            ariCoder.cumulative_frequency = (Prob *) malloc(ariCoder.numOfRealStates * sizeof(Prob));
            memset(ariCoder.cumulative_frequency, 0, ariCoder.numOfRealStates * sizeof(Prob));

            size_t i = 0;
            const uchar *low_p = NULL, *high_p = NULL, *state_p = NULL;
            int state = 0;
            if (total_frequency <= 65536) {
                if (numOfRealStates <= 256) {
                    for (i = 0; i < numOfValidStates; i++) {
                        low_p = p;
                        high_p = low_p + sizeof(uint16_t);
                        state_p = high_p + sizeof(uint16_t);
                        state = *state_p;
                        ariCoder.cumulative_frequency[state].low = bytesToUInt16_bigEndian(low_p);
                        ariCoder.cumulative_frequency[state].high = bytesToUInt16_bigEndian(high_p);
                        ariCoder.cumulative_frequency[state].state = state;

                        p = state_p + 1;
                    }
                    offset = 2 * sizeof(int) + sizeof(uint64_t) +
                             ariCoder.numOfValidStates * 5; //2*sizeof(uint16_t)+1
                } else if (numOfRealStates <= 65536) {
                    for (i = 0; i < numOfValidStates; i++) {
                        low_p = p;
                        high_p = low_p + sizeof(uint16_t);
                        state_p = high_p + sizeof(uint16_t);
                        state = bytesToUInt16_bigEndian(state_p);

                        ariCoder.cumulative_frequency[state].low = bytesToUInt16_bigEndian(low_p);
                        ariCoder.cumulative_frequency[state].high = bytesToUInt16_bigEndian(high_p);
                        ariCoder.cumulative_frequency[state].state = state;

                        p = state_p + sizeof(uint16_t);
                    }
                    offset = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 6;
                } else {
                    for (i = 0; i < numOfValidStates; i++) {
                        low_p = p;
                        high_p = low_p + sizeof(uint16_t);
                        state_p = high_p + sizeof(uint16_t);
                        state = bytesToUInt32_bigEndian(state_p);

                        ariCoder.cumulative_frequency[state].low = bytesToUInt16_bigEndian(low_p);
                        ariCoder.cumulative_frequency[state].high = bytesToUInt16_bigEndian(high_p);
                        ariCoder.cumulative_frequency[state].state = state;

                        p = state_p + sizeof(uint32_t);
                    }
                    offset = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 8;
                }
            } else if (total_frequency <= 4294967296) {
                if (numOfRealStates <= 256) {
                    for (i = 0; i < numOfValidStates; i++) {
                        low_p = p;
                        high_p = low_p + sizeof(uint32_t);
                        state_p = high_p + sizeof(uint32_t);
                        state = *state_p;

                        ariCoder.cumulative_frequency[state].low = bytesToUInt32_bigEndian(low_p);
                        ariCoder.cumulative_frequency[state].high = bytesToUInt32_bigEndian(high_p);
                        ariCoder.cumulative_frequency[state].state = state;

                        p = state_p + 1;
                    }
                    offset = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 9;
                } else if (numOfRealStates <= 65536) {
                    for (i = 0; i < numOfValidStates; i++) {
                        low_p = p;
                        high_p = low_p + sizeof(uint32_t);
                        state_p = high_p + sizeof(uint32_t);
                        state = bytesToUInt16_bigEndian(state_p);

                        ariCoder.cumulative_frequency[state].low = bytesToUInt32_bigEndian(low_p);
                        ariCoder.cumulative_frequency[state].high = bytesToUInt32_bigEndian(high_p);
                        ariCoder.cumulative_frequency[state].state = state;

                        p = state_p + sizeof(uint16_t);
                    }
                    offset = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 10;
                } else {
                    for (i = 0; i < numOfValidStates; i++) {
                        low_p = p;
                        high_p = low_p + sizeof(uint32_t);
                        state_p = high_p + sizeof(uint32_t);
                        state = bytesToUInt32_bigEndian(state_p);

                        ariCoder.cumulative_frequency[state].low = bytesToUInt32_bigEndian(low_p);
                        ariCoder.cumulative_frequency[state].high = bytesToUInt32_bigEndian(high_p);
                        ariCoder.cumulative_frequency[state].state = state;

                        p = state_p + sizeof(uint32_t);
                    }
                    offset = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 12;
                }
            } else {
                if (numOfRealStates <= 256) {
                    for (i = 0; i < numOfValidStates; i++) {
                        low_p = p;
                        high_p = low_p + sizeof(uint64_t);
                        state_p = high_p + sizeof(uint64_t);
                        state = *state_p;

                        ariCoder.cumulative_frequency[state].low = bytesToUInt64_bigEndian(low_p);
                        ariCoder.cumulative_frequency[state].high = bytesToUInt64_bigEndian(high_p);
                        ariCoder.cumulative_frequency[state].state = state;

                        p = state_p + 1;
                    }
                    offset = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 17;
                } else if (numOfRealStates <= 65536) {
                    for (i = 0; i < numOfValidStates; i++) {
                        low_p = p;
                        high_p = low_p + sizeof(uint64_t);
                        state_p = high_p + sizeof(uint64_t);
                        state = bytesToUInt16_bigEndian(state_p);

                        ariCoder.cumulative_frequency[state].low = bytesToUInt64_bigEndian(low_p);
                        ariCoder.cumulative_frequency[state].high = bytesToUInt64_bigEndian(high_p);
                        ariCoder.cumulative_frequency[state].state = state;

                        p = state_p + sizeof(uint16_t);
                    }
                    offset = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 18;
                } else {
                    for (i = 0; i < numOfValidStates; i++) {
                        low_p = p;
                        high_p = low_p + sizeof(uint64_t);
                        state_p = high_p + sizeof(uint64_t);
                        state = bytesToUInt32_bigEndian(state_p);

                        ariCoder.cumulative_frequency[state].low = bytesToUInt64_bigEndian(low_p);
                        ariCoder.cumulative_frequency[state].high = bytesToUInt64_bigEndian(high_p);
                        ariCoder.cumulative_frequency[state].state = state;

                        p = state_p + sizeof(uint32_t);
                    }
                    offset = 2 * sizeof(int) + sizeof(uint64_t) + ariCoder.numOfValidStates * 20;
                }
            }
            remaining_length -= offset;
        }

/**
 * Arithmetic Encoding
 * @param AriCoder *ariCoder (input)
 * @param int *s (input)
 * @param size_t length (input)
 * @param unsigned char *out (output)
 * @param size_t *outSize (output)
 *
 * */
        //        void ari_encode(AriCoder *ariCoder, int *s, size_t length, unsigned char *out, size_t *outSize) {
        size_t encode(const std::vector<T> &bins, uchar *&bytes) {
            const T *s = transform ? bins_transform.data() : bins.data();
            size_t length = transform ? bins_transform.size() : bins.size();
//            unsigned char *bytes = out;
            size_t outSize = 0;

            int pending_bits = 0;
            size_t low = 0;
            size_t high = MAX_CODE;
            size_t i = 0, range = 0;
            size_t count = ariCoder.total_frequency;
            int c = 0, lackBits = 0;


            Prob *cumulative_frequency = ariCoder.cumulative_frequency;
            unsigned int buf = 0;

            for (i = 0; i < length; i++) {
                c = s[i];
                Prob p = cumulative_frequency[c];
                range = high - low + 1;
                high = low + (range * p.high / count) - 1;
                low = low + (range * p.low / count);
                for (;;) {
                    if (high < ONE_HALF) {
                        buf = output_bit_0_plus_pending(pending_bits);
                        put_codes_to_output(buf, pending_bits + 1, &bytes, &lackBits, &outSize);
                        pending_bits = 0;
                    } else if (low >= ONE_HALF) {
                        buf = output_bit_1_plus_pending(pending_bits);
                        put_codes_to_output(buf, pending_bits + 1, &bytes, &lackBits, &outSize);
                        pending_bits = 0;
                    } else if (low >= ONE_FOURTH && high < THREE_FOURTHS) {
                        pending_bits++;
                        low -= ONE_FOURTH;
                        high -= ONE_FOURTH;
                    } else
                        break;
                    high <<= 1;
                    high++;
                    low <<= 1;
                    high &= MAX_CODE;
                    low &= MAX_CODE;
                }
            }
            pending_bits++;
            if (low < ONE_FOURTH) {
                buf = output_bit_0_plus_pending(pending_bits);
                put_codes_to_output(buf, pending_bits + 1, &bytes, &lackBits, &outSize);
            } else {
                buf = output_bit_1_plus_pending(pending_bits);
                put_codes_to_output(buf, pending_bits + 1, &bytes, &lackBits, &outSize);
            }
            bytes += 1;
            return outSize;
        }

        /**
 * Arithmetic Decoding algorithm
 * @param AriCoder *ariCoder (input): the encoder with the constructed frequency information
 * @param unsigned char *s (input): the compressed stream of bytes
 * @param size_t s_len (input): the number of bytes in the 'unsigned char *s'
 * @param size_t targetLength (input): the target number of elements in the type array
 * @param int *out (output) : the result (type array decompressed from the stream 's')
 *
 * */
        std::vector<T> decode(const uchar *&bytes, size_t targetLength) {
            std::vector<T> out(targetLength);

//        void ari_decode(AriCoder *ariCoder, unsigned char *s, size_t s_len, size_t targetLength, int *out) {
            size_t high = MAX_CODE;
            size_t low = 0, i = 0;
            size_t range = 0, scaled_value = 0;
            size_t total_frequency = ariCoder.total_frequency;
            const uchar *sp = bytes + 5;
            unsigned int offset = 4;
            size_t value = (bytesToUInt64_bigEndian(bytes) >> 20); //alignment with the MAX_CODE
            size_t s_counter = sizeof(int);

            for (i = 0; i < targetLength; i++) {
                range = high - low + 1;
                scaled_value = ((value - low + 1) * ariCoder.total_frequency - 1) / range;
                Prob *p = getCode(scaled_value);
//                out[i] = p->state;  //output the state to the 'out' array
                if (transform) {
                    T x = p->state;  //output the state to the 'out' array
                    if (x % 2 == 0) {
                        out[i] = ariCoder.numOfRealStates / 2 + std::ceil(x / 2.0);
                    } else {
                        out[i] = ariCoder.numOfRealStates / 2 - std::ceil(x / 2.0);
                    }
                } else {
                    out[i] = p->state;  //output the state to the 'out' array
                }

                if (i == targetLength - 1) {
                    break;
                }
                high = low + (range * p->high) / total_frequency - 1;
                low = low + (range * p->low) / total_frequency;

                for (;;) {
                    if (high < ONE_HALF) {
                        //do nothing, bit is a zero
                    } else if (low >= ONE_HALF) {
                        value -= ONE_HALF;  //subtract one half from all three code values
                        low -= ONE_HALF;
                        high -= ONE_HALF;
                    } else if (low >= ONE_FOURTH && high < THREE_FOURTHS) {
                        value -= ONE_FOURTH;
                        low -= ONE_FOURTH;
                        high -= ONE_FOURTH;
                    } else
                        break;
                    low <<= 1;
                    high <<= 1;
                    high++;
                    value <<= 1;
                    //load one bit from the input byte stream
//                    if (s_counter < s_len) {
                    value += get_bit(sp, offset++);
                    if (offset == 8) {
                        sp++;
                        s_counter++;
                        offset = 0;
                    }
//                    }
                }
            }
            bytes += s_counter;
            return out;
        }

        AriCoder ariCoder;

    private:
        bool transform;
        std::vector<T> bins_transform;

        inline void output_bit_1(unsigned int *buf) {
            (*buf) = (*buf) << 1;
            (*buf) |= 1;
        }

        inline void output_bit_0(unsigned int *buf) {
            (*buf) = (*buf) << 1;
            //(*byte) |= 0; //actually doesn't have to set the bit to 0
        }

        //TODO: problematic
        inline unsigned int output_bit_1_plus_pending(int pending_bits) {
            unsigned int buf = 0, pbits = pending_bits;
            output_bit_1(&buf);
            while (pbits--)
                output_bit_0(&buf);
            buf = buf << (32 - (pending_bits +
                                1)); //alignment to the left leading bit, which would be easier for the final output
            return buf;
        }

        inline unsigned int output_bit_0_plus_pending(int pending_bits) {
            unsigned int buf = 0, pbits = pending_bits;
            //output_bit_0(&buf);
            while (pbits--)
                output_bit_1(&buf);
            buf = buf << (32 - (pending_bits + 1)); //alignment to the left leading bit
            return buf;
        }

/**
 * Get the integer code based on Arithmetic Coding Value
 * @param AriCoder *ariCoder (input)
 * @param size_t scaled_value (input)
 *
 * @return Prob* (output)
 *
 * */
        Prob *getCode(size_t scaled_value) {
            int numOfRealStates = ariCoder.numOfRealStates;
            int i = 0;
            Prob *p = ariCoder.cumulative_frequency;
            for (i = 0; i < numOfRealStates; i++, p++) {
                if (scaled_value < p->high)
                    break;
            }
            return p;
        }

/**
 * Get one bit from the input stream of bytes
 * @param unsigned char* p (input): the current location to be read (byte) of the byte stream
 * @param int offset (input): the offset of the specified byte in the byte stream
 *
 * @return unsigned char (output) : 1 or 0
 * */
        inline unsigned char get_bit(const uchar *p, int offset) {
            return ((*p) >> (7 - offset)) & 0x01;
        }

        /**
 * put 'buf_nbBits' bits represented by buf into a long byte stream (the current output byte pointer is p, where offset is the number of bits already filled out for this byte so far)
 * */
        void put_codes_to_output(unsigned int buf, int bitSize, unsigned char **p, int *lackBits, size_t *outSize) {
            int byteSize, byteSizep;
            if (*lackBits == 0) {
                byteSize = bitSize % 8 == 0 ? bitSize / 8 : bitSize / 8 +
                                                            1; //it's equal to the number of bytes involved (for *outSize)
                byteSizep = bitSize >> 3; //it's used to move the pointer p for next data
                intToBytes_bigEndian(*p, buf);
                (*p) += byteSizep;
                *outSize += byteSize;
                (*lackBits) = bitSize % 8 == 0 ? 0 : 8 - bitSize % 8;
            } else {
                **p = (**p) | (unsigned char) (buf >> (32 - *lackBits));
                if ((*lackBits) < bitSize) {
                    (*p)++;
                    int newCode = buf << (*lackBits);
                    intToBytes_bigEndian(*p, newCode);
                    bitSize -= *lackBits;
                    byteSizep = bitSize >> 3; // =bitSize/8
                    byteSize = bitSize % 8 == 0 ? byteSizep : byteSizep + 1;
                    *p += byteSizep;
                    (*outSize) += byteSize;
                    (*lackBits) = bitSize % 8 == 0 ? 0 : 8 - bitSize % 8;
                } else {
                    (*lackBits) -= bitSize;
                    if (*lackBits == 0)
                        (*p)++;
                }
            }
        }

    };

}
#endif /* ----- #ifndef _ArithmeticEncoder_H  ----- */

