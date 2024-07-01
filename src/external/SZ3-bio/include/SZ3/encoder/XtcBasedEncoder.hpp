/*! \file XTC compression encoder
 * This is based on libxdrf.cpp from GROMACS (2024-05-18)
 * License: LGPL 2.1 or later
 * \author The GROMACS Authors
 * \author Magnus Lundborg: Modifications to fit as SZ3 encoder
 */

#ifndef _SZ_XTC3_ENCODER_HPP
#define _SZ_XTC3_ENCODER_HPP

#include <climits>

#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"

//#define DEBUG_OUTPUT

/* What follows are the C routine to read/write compressed coordinates together
 * with some routines to assist in this task (those are marked
 * static and cannot be called from user programs)
 */

// Integers above 2^24 do not have unique representations in
// 32-bit floats ie with 24 bits of precision.  We use maxAbsoluteInt
// to check that float values can be transformed into an in-range
// 32-bit integer. There is no need to ensure we are within the range
// of ints with exact floating-point representations. However, we should
// reject all floats above that which converts to an in-range 32-bit integer.
const float maxAbsoluteInt = std::nextafterf(float(INT_MAX), 0.F); // NOLINT(cert-err58-cpp)

#ifndef SQR
#    define SQR(x) ((x) * (x))
#endif
static const int magicInts[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8,
    10, 12, 16, 20, 25, 32, 40, 50, 64, 80,
    101, 128, 161, 203, 256, 322, 406, 512, 645, 812,
    1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192,
    10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570,
    104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, 832255,
    1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607,
    10568983, 13316085, 16777216
};

#define FIRSTIDX 9
/* note that magicInts[FIRSTIDX-1] == 0 */
#define LASTIDX static_cast<int>((sizeof(magicInts) / sizeof(*magicInts)))

namespace SZ3 {
    
    struct DataBuffer {
        std::size_t index;
        int lastbits;
        unsigned int lastbyte;
        unsigned char *data;
    };

/*! \brief encode num into buf using the specified number of bits
 *
 * This routines appends the value of num to the bits already present in
 * the databuffer. You need to give it the number of bits to use and you
 * better make sure that this number of bits is enough to hold the value
 * Also num must be positive.
 *
 */
    
    static void sendbits(struct DataBuffer *buffer, int num_of_bits, int num) {
        
        unsigned int lastbyte;
        int lastbits;
        
        lastbits = buffer->lastbits;
        lastbyte = buffer->lastbyte;
        while (num_of_bits >= CHAR_BIT) {
            lastbyte = (lastbyte << CHAR_BIT) | ((num >> (num_of_bits - CHAR_BIT)) /* & 0xff*/);
            buffer->data[buffer->index++] = lastbyte >> lastbits;
            num_of_bits -= CHAR_BIT;
        }
        if (num_of_bits > 0) {
            lastbyte = (lastbyte << num_of_bits) | num;
            lastbits += num_of_bits;
            if (lastbits >= CHAR_BIT) {
                lastbits -= CHAR_BIT;
                buffer->data[buffer->index++] = lastbyte >> lastbits;
            }
        }
        buffer->lastbits = lastbits;
        buffer->lastbyte = lastbyte;
        if (lastbits > 0) {
            buffer->data[buffer->index] = lastbyte << (CHAR_BIT - lastbits);
        }
    }

/*! \brief calculate bitSize of an integer
 *
 * return the number of bits needed to store an integer with given max size
 *
 */
    
    static int sizeofint(const int size) {
        int num = 1;
        int num_of_bits = 0;
        
        while (size >= num && num_of_bits < 32) {
            num_of_bits++;
            num <<= 1;
        }
        return num_of_bits;
    }

/*! \brief calculate 'bitSize' of compressed ints
 *
 * given the number of small unsigned integers and the maximum value
 * return the number of bits needed to read or write them with the
 * routines receiveints and sendints. You need this parameter when
 * calling these routines. Note that for many calls I can use
 * the variable 'smallIdx' which is exactly the number of bits, and
 * So I don't need to call 'sizeofints for those calls.
 */
    
    static int sizeofints(const int num_of_ints, const unsigned int sizes[]) {
        int i, num;
        int bytes[32];
        unsigned int num_of_bytes, num_of_bits, bytecnt, tmp;
        num_of_bytes = 1;
        bytes[0] = 1;
        num_of_bits = 0;
        for (i = 0; i < num_of_ints; i++) {
            tmp = 0;
            for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
                tmp = bytes[bytecnt] * sizes[i] + tmp;
                bytes[bytecnt] = tmp & 0xff;
                tmp >>= CHAR_BIT;
            }
            while (tmp != 0) {
                bytes[bytecnt++] = tmp & 0xff;
                tmp >>= CHAR_BIT;
            }
            num_of_bytes = bytecnt;
        }
        num = 1;
        num_of_bytes--;
        while (bytes[num_of_bytes] >= num) {
            num_of_bits++;
            num *= 2;
        }
        return num_of_bits + num_of_bytes * CHAR_BIT;
    }

/*! \brief send a small set of small integers in compressed format
 *
 * this routine is used internally by xdr3dfcoord, to send a set of
 * small integers to the buffer.
 * Multiplication with fixed (specified maximum ) sizes is used to get
 * to one big, multibyte integer. Allthough the routine could be
 * modified to handle sizes bigger than 16777216, or more than just
 * a few integers, this is not done, because the gain in compression
 * isn't worth the effort. Note that overflowing the multiplication
 * or the byte buffer (32 bytes) is unchecked and causes bad results.
 *
 */
    
    static void sendints(struct DataBuffer *buffer,
                         const int num_of_ints,
                         const int num_of_bits,
                         unsigned int sizes[],
                         unsigned int nums[]) {
        
        int i, num_of_bytes, bytecnt;
        unsigned int bytes[32], tmp;
        
        tmp = nums[0];
        num_of_bytes = 0;
        do {
            bytes[num_of_bytes++] = tmp & 0xff;
            tmp >>= CHAR_BIT;
        } while (tmp != 0);
        
        for (i = 1; i < num_of_ints; i++) {
            if (nums[i] >= sizes[i]) {
                fprintf(stderr,
                        "major breakdown in sendints num %u doesn't "
                        "match size %u\n",
                        nums[i],
                        sizes[i]);
                exit(1);
            }
            /* use one step multiply */
            tmp = nums[i];
            for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
                tmp = bytes[bytecnt] * sizes[i] + tmp;
                bytes[bytecnt] = tmp & 0xff;
                tmp >>= CHAR_BIT;
            }
            while (tmp != 0) {
                bytes[bytecnt++] = tmp & 0xff;
                tmp >>= CHAR_BIT;
            }
            num_of_bytes = bytecnt;
        }
        if (num_of_bits >= num_of_bytes * CHAR_BIT) {
            for (i = 0; i < num_of_bytes; i++) {
                sendbits(buffer, CHAR_BIT, bytes[i]);
            }
            sendbits(buffer, num_of_bits - num_of_bytes * CHAR_BIT, 0);
        } else {
            for (i = 0; i < num_of_bytes - 1; i++) {
                sendbits(buffer, CHAR_BIT, bytes[i]);
            }
            sendbits(buffer, num_of_bits - (num_of_bytes - 1) * CHAR_BIT, bytes[i]);
        }
    }

/*! \brief decode number from buffer using specified number of bits
 *
 * extract the number of bits from the data array in buffer and construct an integer
 * from it. Return that value.
 *
 */
    
    static int receivebits(struct DataBuffer *buffer, int num_of_bits) {
        
        int num, lastbits;
        unsigned int lastbyte;
        int mask = (1 << num_of_bits) - 1;
        
        lastbits = buffer->lastbits;
        lastbyte = buffer->lastbyte;
        
        num = 0;
        while (num_of_bits >= CHAR_BIT) {
            lastbyte = (lastbyte << CHAR_BIT) | buffer->data[buffer->index++];
            num |= (lastbyte >> lastbits) << (num_of_bits - CHAR_BIT);
            num_of_bits -= CHAR_BIT;
        }
        if (num_of_bits > 0) {
            if (lastbits < num_of_bits) {
                lastbits += CHAR_BIT;
                lastbyte = (lastbyte << CHAR_BIT) | buffer->data[buffer->index++];
            }
            lastbits -= num_of_bits;
            num |= (lastbyte >> lastbits) & ((1 << num_of_bits) - 1);
        }
        num &= mask;
        buffer->lastbits = lastbits;
        buffer->lastbyte = lastbyte;
        return num;
    }

/*! \brief  decode 'small' integers from the buf array
 *
 * this routine is the inverse from sendints() and decodes the small integers
 * written to buf by calculating the remainder and doing divisions with
 * the given sizes[]. You need to specify the total number of bits to be
 * used from buf in num_of_bits.
 *
 */
    
    static void receiveints(struct DataBuffer *buffer,
                            const int num_of_ints,
                            int num_of_bits,
                            const unsigned int sizes[],
                            int nums[]) {
        int bytes[32];
        int i, j, num_of_bytes, p, num;
        
        bytes[0] = bytes[1] = bytes[2] = bytes[3] = 0;
        num_of_bytes = 0;
        while (num_of_bits > CHAR_BIT) {
            bytes[num_of_bytes++] = receivebits(buffer, CHAR_BIT);
            num_of_bits -= CHAR_BIT;
        }
        if (num_of_bits > 0) {
            bytes[num_of_bytes++] = receivebits(buffer, num_of_bits);
        }
        for (i = num_of_ints - 1; i > 0; i--) {
            num = 0;
            for (j = num_of_bytes - 1; j >= 0; j--) {
                num = (num << CHAR_BIT) | bytes[j];
                p = num / sizes[i];
                bytes[j] = p;
                num = num - p * sizes[i];
            }
            nums[i] = num;
        }
        nums[0] = bytes[0] | (bytes[1] << CHAR_BIT) | (bytes[2] << 16) | (bytes[3] << 24);
    }
    
    template<class T>
    class XtcBasedEncoder : public concepts::EncoderInterface<T> {
     
     private:
        Config conf_;
     
     public:
        void preprocess_encode(const std::vector<T> &quantData, int stateNum) {}
        
        /*! \brief Compress 3d coordinates to memory.
         *
         * this routine writes a large number of, already quantized (as integers with a
         * given precision) 3d coordinates.
         * The minimum and maximum value are calculated to determine the range.
         * The limited range of integers so found, is used to compress the coordinates.
         * In addition the differences between succesive coordinates is calculated.
         * If the difference happens to be 'small' then only the difference is saved,
         * compressing the data even more. The notion of 'small' is changed dynamically
         * and is enlarged or reduced whenever needed or possible.
         * Extra compression is achieved in the case of commonly used water models.
         * In those the Oxygen position is followed by  the two hydrogens. In order to
         * make the differences smaller (and thereby compression the data better) the
         * order is changed into first one hydrogen then the oxygen, followed by the
         * other hydrogen. This is rather special, but it shouldn't harm in the general case.
         *
         */
        size_t encode(const std::vector<T> &quantData, unsigned char *&bytes) {
            size_t size3 = quantData.size();
            
            size_t bufferSize = size3 * 1.2;
            struct DataBuffer buffer;
            int *intBufferPoiner = reinterpret_cast<int *>(malloc(size3 * sizeof(*intBufferPoiner)));
            buffer.data = reinterpret_cast<unsigned char *>(malloc(bufferSize * sizeof(int)));
            if (buffer.data == nullptr) {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
            buffer.index = 0;
            buffer.lastbits = 0;
            buffer.lastbyte = 0;
            
            unsigned char *charOutputPtr = bytes;
            unsigned int *intOutputPtr = reinterpret_cast<unsigned int *>(charOutputPtr);
            uint64_t numTriplets = size3 / 3;

#ifdef DEBUG_OUTPUT
            printf("Encoding: %llu triplets to write.\n", numTriplets);
            for (size_t i = 0; i < numTriplets; i++)
            {
                printf("Triplet %ld: %d %d %d\n", i, quantData[i * 3], quantData[i * 3 + 1], quantData[i * 3 + 2]);
            }
#endif
            
            size_t coordDataOffset = reinterpret_cast<unsigned char *>(intOutputPtr) - charOutputPtr;
            int *localIntBufferPointer = intBufferPoiner;
            int minInt[3] = {INT_MAX, INT_MAX, INT_MAX};
            int maxInt[3] = {INT_MIN, INT_MIN, INT_MIN};
            int minDiff = INT_MAX;
            int oldLocalValue1 = 0;
            int oldLocalValue2 = 0;
            int oldLocalValue3 = 0;
            const int *inputDataPtr = quantData.data();
            while (inputDataPtr < quantData.data() + size3) {
                int localValue1 = *inputDataPtr++;
                if (localValue1 < minInt[0]) {
                    minInt[0] = localValue1;
                }
                if (localValue1 > maxInt[0]) {
                    maxInt[0] = localValue1;
                }
                *localIntBufferPointer++ = localValue1;
                
                int localValue2 = *inputDataPtr++;
                if (localValue2 < minInt[1]) {
                    minInt[1] = localValue2;
                }
                if (localValue2 > maxInt[1]) {
                    maxInt[1] = localValue2;
                }
                *localIntBufferPointer++ = localValue2;
                
                int localValue3 = *inputDataPtr++;
                if (localValue3 < minInt[2]) {
                    minInt[2] = localValue3;
                }
                if (localValue3 > maxInt[2]) {
                    maxInt[2] = localValue3;
                }
                *localIntBufferPointer++ = localValue3;
                int diff = std::abs(oldLocalValue1 - localValue1) + std::abs(oldLocalValue2 - localValue2)
                    + std::abs(oldLocalValue3 - localValue3);
                if (diff < minDiff && inputDataPtr > quantData.data() + 3) {
                    minDiff = diff;
                }
                oldLocalValue1 = localValue1;
                oldLocalValue2 = localValue2;
                oldLocalValue3 = localValue3;
            }
            
            for (int i = 0; i < 3; i++) {
                *intOutputPtr++ = minInt[i];
            }
            for (int i = 0; i < 3; i++) {
                *intOutputPtr++ = maxInt[i];
            }
            
            if (static_cast<float>(maxInt[0]) - static_cast<float>(minInt[0]) >= maxAbsoluteInt
                || static_cast<float>(maxInt[1]) - static_cast<float>(minInt[1]) >= maxAbsoluteInt
                || static_cast<float>(maxInt[2]) - static_cast<float>(minInt[2]) >= maxAbsoluteInt) {
                /* turning value in unsigned by subtracting minInt
                 * would cause overflow
                 */
                fprintf(stderr,
                        "Error. Turning value in unsigned by subtracting minInt would cause "
                        "overflow.\n");
            }
            unsigned int sizeInt[3];
            sizeInt[0] = maxInt[0] - minInt[0] + 1;
            sizeInt[1] = maxInt[1] - minInt[1] + 1;
            sizeInt[2] = maxInt[2] - minInt[2] + 1;
            unsigned int bitSizeInt[3];
            int bitSize;

#ifdef DEBUG_OUTPUT
            printf("    minInt %d %d %d, maxInt %d %d %d\n",
                   minInt[0],
                   minInt[1],
                   minInt[2],
                   maxInt[0],
                   maxInt[1],
                   maxInt[2]);
#endif
            
            /* check if one of the sizes is too big to be multiplied */
            if ((sizeInt[0] | sizeInt[1] | sizeInt[2]) > 0xffffff) {
                bitSizeInt[0] = sizeofint(sizeInt[0]);
                bitSizeInt[1] = sizeofint(sizeInt[1]);
                bitSizeInt[2] = sizeofint(sizeInt[2]);
                bitSize = 0; /* flag the use of large sizes */
            } else {
                bitSize = sizeofints(3, sizeInt);
            }
            int smallIdx = FIRSTIDX;
            while (smallIdx < LASTIDX && magicInts[smallIdx] < minDiff) {
                smallIdx++;
            }
            *intOutputPtr++ = smallIdx;
            
            int maxIdx = std::min(LASTIDX, smallIdx + CHAR_BIT);
            int minIdx = maxIdx - CHAR_BIT; /* often this equal smallIdx */
            int smaller = magicInts[std::max(FIRSTIDX, smallIdx - 1)] / 2;
            int smallNum = magicInts[smallIdx] / 2;
            unsigned int sizeSmall[3];
            sizeSmall[0] = magicInts[smallIdx];
            sizeSmall[1] = magicInts[smallIdx];
            sizeSmall[2] = magicInts[smallIdx];
            int larger = magicInts[maxIdx] / 2;
            size_t i = 0;
            unsigned int *localUnsignedIntBufferPointer = reinterpret_cast<unsigned int *>(intBufferPoiner);
            int prevCoord[3] = {0, 0, 0};
            int prevRun = -1;
            while (i < numTriplets) {
                bool isSmall = false;
                int *thisCoord = reinterpret_cast<int *>(localUnsignedIntBufferPointer) + i * 3;
                int isSmaller;
                if (smallIdx < maxIdx && i >= 1 && std::abs(thisCoord[0] - prevCoord[0]) < larger
                    && std::abs(thisCoord[1] - prevCoord[1]) < larger
                    && std::abs(thisCoord[2] - prevCoord[2]) < larger) {
                    isSmaller = 1;
                } else if (smallIdx > minIdx) {
                    isSmaller = -1;
                } else {
                    isSmaller = 0;
                }
                if (i + 1 < numTriplets) {
                    if (std::abs(thisCoord[0] - thisCoord[3]) < smallNum
                        && std::abs(thisCoord[1] - thisCoord[4]) < smallNum
                        && std::abs(thisCoord[2] - thisCoord[5]) < smallNum) {
                        /* interchange first with second atom for better
                         * compression of water molecules
                         */
                        int tmp = thisCoord[0];
                        thisCoord[0] = thisCoord[3];
                        thisCoord[3] = tmp;
                        tmp = thisCoord[1];
                        thisCoord[1] = thisCoord[4];
                        thisCoord[4] = tmp;
                        tmp = thisCoord[2];
                        thisCoord[2] = thisCoord[5];
                        thisCoord[5] = tmp;
                        isSmall = true;
                    }
                }
                unsigned int tmpCoord[30];
                tmpCoord[0] = thisCoord[0] - minInt[0];
                tmpCoord[1] = thisCoord[1] - minInt[1];
                tmpCoord[2] = thisCoord[2] - minInt[2];
                if (bitSize == 0) {
                    sendbits(&buffer, bitSizeInt[0], tmpCoord[0]);
                    sendbits(&buffer, bitSizeInt[1], tmpCoord[1]);
                    sendbits(&buffer, bitSizeInt[2], tmpCoord[2]);
                } else {
                    sendints(&buffer, 3, bitSize, sizeInt, tmpCoord);
                }
                prevCoord[0] = thisCoord[0];
                prevCoord[1] = thisCoord[1];
                prevCoord[2] = thisCoord[2];
                thisCoord = thisCoord + 3;
                i++;
                
                int run = 0;
                if (!isSmall && isSmaller == -1) {
                    isSmaller = 0;
                }
                while (isSmall && run < CHAR_BIT * 3) {
                    if (isSmaller == -1
                        && (SQR(thisCoord[0] - prevCoord[0]) + SQR(thisCoord[1] - prevCoord[1])
                            + SQR(thisCoord[2] - prevCoord[2])
                            >= smaller * smaller)) {
                        isSmaller = 0;
                    }
                    
                    tmpCoord[run++] = thisCoord[0] - prevCoord[0] + smallNum;
                    tmpCoord[run++] = thisCoord[1] - prevCoord[1] + smallNum;
                    tmpCoord[run++] = thisCoord[2] - prevCoord[2] + smallNum;
                    
                    prevCoord[0] = thisCoord[0];
                    prevCoord[1] = thisCoord[1];
                    prevCoord[2] = thisCoord[2];
                    
                    i++;
                    thisCoord = thisCoord + 3;
                    isSmall = 0;
                    if (i < numTriplets && abs(thisCoord[0] - prevCoord[0]) < smallNum
                        && abs(thisCoord[1] - prevCoord[1]) < smallNum
                        && abs(thisCoord[2] - prevCoord[2]) < smallNum) {
                        isSmall = true;
                    }
                }
                if (run != prevRun || isSmaller != 0) {
                    prevRun = run;
                    sendbits(&buffer, 1, 1); /* flag the change in run-length */
                    sendbits(&buffer, 5, run + isSmaller + 1);
                } else {
                    sendbits(&buffer, 1, 0); /* flag the fact that runlength did not change */
                }
                for (int k = 0; k < run; k += 3) {
                    sendints(&buffer, 3, smallIdx, sizeSmall, &tmpCoord[k]);
                }
                if (isSmaller != 0) {
                    smallIdx += isSmaller;
                    if (isSmaller < 0) {
                        smallNum = smaller;
                        smaller = magicInts[smallIdx - 1] / 2;
                    } else {
                        smaller = smallNum;
                        smallNum = magicInts[smallIdx] / 2;
                    }
                    sizeSmall[0] = sizeSmall[1] = sizeSmall[2] = magicInts[smallIdx];
                }
            }
            if (buffer.lastbits != 0) {
                buffer.index++;
            }
            
            *(reinterpret_cast<uint64_t *>(intOutputPtr)) = buffer.index;
            intOutputPtr += sizeof(uint64_t) / sizeof(int);
            
            // Since this file is full of old code, and many signed-to-unsigned conversions, we
            // read data in batches if the smallest number that is a multiple of 4 that
            // fits in a signed integer to keep data access aligned if possible.
            size_t offset = 0;
            size_t remain = buffer.index;
            charOutputPtr = reinterpret_cast<unsigned char *>(intOutputPtr);
            do {
                // Max batch size is largest 4-tuple that fits in signed 32-bit int
                size_t batchSize = std::min(remain, static_cast<std::size_t>(2147483644));
                memcpy(charOutputPtr, buffer.data + offset, batchSize);
                charOutputPtr += batchSize;
                offset += batchSize;
                remain -= batchSize;
            } while (remain > 0);
            
            free(buffer.data);
            free(intBufferPoiner);
            
            size_t outputSize = charOutputPtr - bytes;

#ifdef DEBUG_OUTPUT
            printf("Finished encoding. Output length = %ld bytes, compression = %f.\n\n",
                   outputSize,
                   (float)size3 / outputSize);
#endif
            
            bytes = charOutputPtr;
            return outputSize;
        }
        
        void postprocess_encode() {}
        
        void preprocess_decode() {}
        
        /*! \brief Decompress 3d coordinates to memory.
         *
         * this routine decompresses a large number of compressed 3d coordinates.
         *
         */
        std::vector<T> decode(const unsigned char *&bytes, size_t targetLength) {
#ifdef DEBUG_OUTPUT
            printf("\nDecoding, targetLength: %ld\n", targetLength);
#endif
            std::vector<T> quantData(targetLength, 0);
            
            const unsigned char *inputBytesPointer = bytes;
            const int *inputIntPtr = reinterpret_cast<const int *>(inputBytesPointer);
            
            size_t bufferSize = targetLength * 1.2;
            struct DataBuffer buffer;
            buffer.data = reinterpret_cast<unsigned char *>(malloc(bufferSize * sizeof(int)));
            if (buffer.data == nullptr) {
                fprintf(stderr, "malloc failed\n");
            }
            buffer.index = 0;
            buffer.lastbits = 0;
            buffer.lastbyte = 0;
            
            int minInt[3];
            int maxInt[3];
            
            minInt[0] = *inputIntPtr++;
            minInt[1] = *inputIntPtr++;
            minInt[2] = *inputIntPtr++;
            maxInt[0] = *inputIntPtr++;
            maxInt[1] = *inputIntPtr++;
            maxInt[2] = *inputIntPtr++;

#ifdef DEBUG_OUTPUT
            printf("    minInt %d %d %d, maxInt %d %d %d\n",
                   minInt[0],
                   minInt[1],
                   minInt[2],
                   maxInt[0],
                   maxInt[1],
                   maxInt[2]);
#endif
            
            unsigned int sizeInt[3];
            sizeInt[0] = maxInt[0] - minInt[0] + 1;
            sizeInt[1] = maxInt[1] - minInt[1] + 1;
            sizeInt[2] = maxInt[2] - minInt[2] + 1;
            unsigned int bitSizeInt[3];
            int bitSize;
            /* check if one of the sizes is too big to be multiplied */
            if ((sizeInt[0] | sizeInt[1] | sizeInt[2]) > 0xffffff) {
                bitSizeInt[0] = sizeofint(sizeInt[0]);
                bitSizeInt[1] = sizeofint(sizeInt[1]);
                bitSizeInt[2] = sizeofint(sizeInt[2]);
                bitSize = 0; /* flag the use of large sizes */
            } else {
                bitSize = sizeofints(3, sizeInt);
            }
            
            int smallIdx = *inputIntPtr++;
            
            int smaller = magicInts[std::max(FIRSTIDX, smallIdx - 1)] / 2;
            int smallNum = magicInts[smallIdx] / 2;
            unsigned int sizeSmall[3];
            sizeSmall[0] = magicInts[smallIdx];
            sizeSmall[1] = magicInts[smallIdx];
            sizeSmall[2] = magicInts[smallIdx];
            
            size_t size3 = targetLength;
            bufferSize = size3 * 1.2;
            buffer.data = reinterpret_cast<unsigned char *>(malloc(bufferSize * sizeof(int)));
            buffer.index = *(reinterpret_cast<const uint64_t *>(inputIntPtr));
            inputIntPtr += sizeof(uint64_t) / sizeof(int);
            
            size_t offset = 0;
            size_t remain = buffer.index;
            inputBytesPointer = reinterpret_cast<const unsigned char *>(inputIntPtr);
            do {
                // Max batch size is largest 4-tuple that fits in signed 32-bit int
                size_t batchSize = std::min(remain, static_cast<std::size_t>(2147483644));
                memcpy(buffer.data + offset, inputBytesPointer, batchSize);
                inputBytesPointer += batchSize;
                offset += batchSize;
                remain -= batchSize;
            } while (remain > 0);
            
            buffer.index = 0;
            buffer.lastbits = 0;
            buffer.lastbyte = 0;
            
            int run = 0;
            size_t i = 0;
            int *intBufferPoiner = reinterpret_cast<int *>(malloc(size3 * sizeof(*intBufferPoiner)));
            int *localIntBufferPointer = intBufferPoiner;
            unsigned char *charOutputPtr = reinterpret_cast<unsigned char *>(quantData.data());
            int *intOutputPtr = reinterpret_cast<int *>(charOutputPtr);
            int prevCoord[3] = {0, 0, 0};
            uint64_t numTriplets = targetLength / 3;
            while (i < numTriplets) {
                int *thisCoord = reinterpret_cast<int *>(localIntBufferPointer) + i * 3;
                
                if (bitSize == 0) {
                    thisCoord[0] = receivebits(&buffer, bitSizeInt[0]);
                    thisCoord[1] = receivebits(&buffer, bitSizeInt[1]);
                    thisCoord[2] = receivebits(&buffer, bitSizeInt[2]);
                } else {
                    receiveints(&buffer, 3, bitSize, sizeInt, thisCoord);
                }
                
                i++;
                thisCoord[0] += minInt[0];
                thisCoord[1] += minInt[1];
                thisCoord[2] += minInt[2];
                
                prevCoord[0] = thisCoord[0];
                prevCoord[1] = thisCoord[1];
                prevCoord[2] = thisCoord[2];
                
                int flag = receivebits(&buffer, 1);
                int isSmaller = 0;
                if (flag == 1) {
                    run = receivebits(&buffer, 5);
                    isSmaller = run % 3;
                    run -= isSmaller;
                    isSmaller--;
                }
                if (run > 0) {
                    thisCoord += 3;
                    for (int k = 0; k < run; k += 3) {
                        receiveints(&buffer, 3, smallIdx, sizeSmall, thisCoord);
                        i++;
                        thisCoord[0] += prevCoord[0] - smallNum;
                        thisCoord[1] += prevCoord[1] - smallNum;
                        thisCoord[2] += prevCoord[2] - smallNum;
                        if (k == 0) {
                            /* interchange first with second atom for better
                             * compression of water molecules
                             */
                            int tmp = thisCoord[0];
                            thisCoord[0] = prevCoord[0];
                            prevCoord[0] = tmp;
                            tmp = thisCoord[1];
                            thisCoord[1] = prevCoord[1];
                            prevCoord[1] = tmp;
                            tmp = thisCoord[2];
                            thisCoord[2] = prevCoord[2];
                            prevCoord[2] = tmp;
                            *intOutputPtr++ = prevCoord[0];
                            *intOutputPtr++ = prevCoord[1];
                            *intOutputPtr++ = prevCoord[2];
                        } else {
                            prevCoord[0] = thisCoord[0];
                            prevCoord[1] = thisCoord[1];
                            prevCoord[2] = thisCoord[2];
                        }
                        *intOutputPtr++ = thisCoord[0];
                        *intOutputPtr++ = thisCoord[1];
                        *intOutputPtr++ = thisCoord[2];
                    }
                } else {
                    *intOutputPtr++ = thisCoord[0];
                    *intOutputPtr++ = thisCoord[1];
                    *intOutputPtr++ = thisCoord[2];
                }
                
                smallIdx += isSmaller;
                if (isSmaller < 0) {
                    smallNum = smaller;
                    if (smallIdx > FIRSTIDX) {
                        smaller = magicInts[smallIdx - 1] / 2;
                    } else {
                        smaller = 0;
                    }
                } else if (isSmaller > 0) {
                    smaller = smallNum;
                    smallNum = magicInts[smallIdx] / 2;
                }
                sizeSmall[0] = sizeSmall[1] = sizeSmall[2] = magicInts[smallIdx];
            }
            free(buffer.data);
            free(intBufferPoiner);

#ifdef DEBUG_OUTPUT
            printf("Decoded %llu triplets.\n", numTriplets);
            for (size_t i = 0; i < numTriplets; i++)
            {
                printf("Triplet %ld: %d %d %d\n", i, quantData[i * 3], quantData[i * 3 + 1], quantData[i * 3 + 2]);
            }
#endif
            
            return quantData;
        }
        
        void postprocess_decode() {}
        
        void save(uchar *&c) {}
        
        void load(const uchar *&c, size_t &remaining_length) {}
    };
    
} // namespace SZ3


#endif
