#ifndef _SZ_HUFFMAN_ENCODER_HPP
#define _SZ_HUFFMAN_ENCODER_HPP

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/utils/ska_hash/unordered_map.hpp"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>


namespace SZ {


    template<class T>
    class HuffmanEncoder : public concepts::EncoderInterface<T> {

    public:

        typedef struct node_t {
            struct node_t *left, *right;
            size_t freq;
            char t; //in_node:0; otherwise:1
            T c;
        } *node;

        typedef struct HuffmanTree {
            unsigned int stateNum;
            unsigned int allNodes;
            struct node_t *pool;
            node *qqq, *qq; //the root node of the HuffmanTree is qq[1]
            int n_nodes; //n_nodes is for compression
            int qend;
            uint64_t **code;
            unsigned char *cout;
            int n_inode; //n_inode is for decompression
            int maxBitCount;
        } HuffmanTree;


        HuffmanEncoder() {
            int x = 1;
            char *y = (char *) &x;
            if (*y == 1)
                sysEndianType = 0;
            else //=0
                sysEndianType = 1;
        }

        ~HuffmanEncoder() {
            SZ_FreeHuffman();
        }

        //build huffman tree
        HuffmanTree *createHuffmanTree(int stateNum) {
            HuffmanTree *huffmanTree = (HuffmanTree *) malloc(sizeof(HuffmanTree));
            memset(huffmanTree, 0, sizeof(HuffmanTree));
            huffmanTree->stateNum = stateNum;
            huffmanTree->allNodes = 2 * stateNum;

            huffmanTree->pool = (struct node_t *) malloc(huffmanTree->allNodes * 2 * sizeof(struct node_t));
            huffmanTree->qqq = (node *) malloc(huffmanTree->allNodes * 2 * sizeof(node));
            huffmanTree->code = (uint64_t **) malloc(huffmanTree->stateNum * sizeof(uint64_t *));
            huffmanTree->cout = (unsigned char *) malloc(huffmanTree->stateNum * sizeof(unsigned char));

            memset(huffmanTree->pool, 0, huffmanTree->allNodes * 2 * sizeof(struct node_t));
            memset(huffmanTree->qqq, 0, huffmanTree->allNodes * 2 * sizeof(node));
            memset(huffmanTree->code, 0, huffmanTree->stateNum * sizeof(uint64_t *));
            memset(huffmanTree->cout, 0, huffmanTree->stateNum * sizeof(unsigned char));
            huffmanTree->qq = huffmanTree->qqq - 1;
            huffmanTree->n_nodes = 0;
            huffmanTree->n_inode = 0;
            huffmanTree->qend = 1;

            return huffmanTree;
        }

        /**
         * build huffman tree using bins
         * @param bins
         * @param stateNum is no longer needed
         */
        void preprocess_encode(const std::vector<T> &bins, int stateNum) {
            preprocess_encode(bins.data(), bins.size(), stateNum);
        }

        /**
         * build huffman tree using bins
         * @param bins
         * @param num_bin
         * @param stateNum is no longer needed
         */
        void preprocess_encode(const T *bins, size_t num_bin, int stateNum) {
            nodeCount = 0;
            if (num_bin == 0) {
                printf("Huffman bins should not be empty\n");
                exit(0);
            }
            init(bins, num_bin);
            for (int i = 0; i < huffmanTree->stateNum; i++)
                if (huffmanTree->code[i]) nodeCount++;
            nodeCount = nodeCount * 2 - 1;
        }

        //save the huffman Tree in the compressed data
        uint save(uchar *&c) {
            auto cc = c;
            write(offset, c);
            int32ToBytes_bigEndian(c, nodeCount);
            c += sizeof(int);
            int32ToBytes_bigEndian(c, huffmanTree->stateNum / 2);
            c += sizeof(int);
            uint totalSize = 0;// = convert_HuffTree_to_bytes_anyStates(nodeCount, c);
            // std::cout << "nodeCount = " << nodeCount << std::endl;
            if (nodeCount <= 256)
                totalSize = convert_HuffTree_to_bytes_anyStates<unsigned char>(nodeCount, c);
            else if (nodeCount <= 65536)
                totalSize = convert_HuffTree_to_bytes_anyStates<unsigned short>(nodeCount, c);
            else
                totalSize = convert_HuffTree_to_bytes_anyStates<unsigned int>(nodeCount, c);
            c += totalSize;
            return c - cc;
        }

        size_t size_est() {
            size_t b = (nodeCount <= 256) ? sizeof(unsigned char) : ((nodeCount <= 65536) ? sizeof(unsigned short) : sizeof(unsigned int));
            return 1 + 2 * nodeCount * b + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T) + sizeof(int) + sizeof(int) + sizeof(T);
        }

        //perform encoding
        size_t encode(const std::vector<T> &bins, uchar *&bytes) {
            return encode(bins.data(), bins.size(), bytes);
        }

        //perform encoding
        size_t encode(const T *bins, size_t num_bin, uchar *&bytes) {
            size_t outSize = 0;
            size_t i = 0;
            unsigned char bitSize = 0, byteSize, byteSizep;
            int state;
            uchar *p = bytes + sizeof(size_t);
            int lackBits = 0;
            //int64_t totalBitSize = 0, maxBitSize = 0, bitSize21 = 0, bitSize32 = 0;
            for (i = 0; i < num_bin; i++) {
                state = bins[i] - offset;
                bitSize = huffmanTree->cout[state];

                if (lackBits == 0) {
                    byteSize = bitSize % 8 == 0 ? bitSize / 8 : bitSize / 8 +
                                                                1; //it's equal to the number of bytes involved (for *outSize)
                    byteSizep = bitSize / 8; //it's used to move the pointer p for next data
                    if (byteSize <= 8) {
                        int64ToBytes_bigEndian(p, (huffmanTree->code[state])[0]);
                        p += byteSizep;
                    } else //byteSize>8
                    {
                        int64ToBytes_bigEndian(p, (huffmanTree->code[state])[0]);
                        p += 8;
                        int64ToBytes_bigEndian(p, (huffmanTree->code[state])[1]);
                        p += (byteSizep - 8);
                    }
                    outSize += byteSize;
                    lackBits = bitSize % 8 == 0 ? 0 : 8 - bitSize % 8;
                } else {
                    *p = (*p) | (unsigned char) ((huffmanTree->code[state])[0] >> (64 - lackBits));
                    if (lackBits < bitSize) {
                        p++;

                        int64_t newCode = (huffmanTree->code[state])[0] << lackBits;
                        int64ToBytes_bigEndian(p, newCode);

                        if (bitSize <= 64) {
                            bitSize -= lackBits;
                            byteSize = bitSize % 8 == 0 ? bitSize / 8 : bitSize / 8 + 1;
                            byteSizep = bitSize / 8;
                            p += byteSizep;
                            outSize += byteSize;
                            lackBits = bitSize % 8 == 0 ? 0 : 8 - bitSize % 8;
                        } else //bitSize > 64
                        {
                            byteSizep = 7; //must be 7 bytes, because lackBits!=0
                            p += byteSizep;
                            outSize += byteSize;

                            bitSize -= 64;
                            if (lackBits < bitSize) {
                                *p = (*p) | (unsigned char) ((huffmanTree->code[state])[0] >> (64 - lackBits));
                                p++;
                                newCode = (huffmanTree->code[state])[1] << lackBits;
                                int64ToBytes_bigEndian(p, newCode);
                                bitSize -= lackBits;
                                byteSize = bitSize % 8 == 0 ? bitSize / 8 : bitSize / 8 + 1;
                                byteSizep = bitSize / 8;
                                p += byteSizep;
                                outSize += byteSize;
                                lackBits = bitSize % 8 == 0 ? 0 : 8 - bitSize % 8;
                            } else //lackBits >= bitSize
                            {
                                *p = (*p) | (unsigned char) ((huffmanTree->code[state])[0] >> (64 - bitSize));
                                lackBits -= bitSize;
                            }
                        }
                    } else //lackBits >= bitSize
                    {
                        lackBits -= bitSize;
                        if (lackBits == 0)
                            p++;
                    }
                }
            }
            *reinterpret_cast<size_t *>(bytes) = outSize;
            bytes += sizeof(size_t) + outSize;
            return outSize;
        }

        void postprocess_encode() {
            SZ_FreeHuffman();
        }

        void preprocess_decode() {};

        //perform decoding
        std::vector<T> decode(const uchar *&bytes, size_t targetLength) {
            node t = treeRoot;
            std::vector<T> out(targetLength);
            size_t i = 0, byteIndex = 0, count = 0;
            int r;
            node n = treeRoot;
            size_t encodedLength = *reinterpret_cast<const size_t *>(bytes);
            bytes += sizeof(size_t);
            if (n->t) //root->t==1 means that all state values are the same (constant)
            {
                for (count = 0; count < targetLength; count++)
                    out[count] = n->c + offset;
                return out;
            }

            for (i = 0; count < targetLength; i++) {
                byteIndex = i >> 3; //i/8
                r = i % 8;
                if (((bytes[byteIndex] >> (7 - r)) & 0x01) == 0)
                    n = n->left;
                else
                    n = n->right;

                if (n->t) {
                    out[count] = n->c + offset;
                    n = t;
                    count++;
                }
            }
            bytes += encodedLength;
            return out;
        }

        //empty function
        void postprocess_decode() {
            SZ_FreeHuffman();
        }

        //load Huffman tree
        void load(const uchar *&c, size_t &remaining_length) {
            read(offset, c, remaining_length);
            nodeCount = bytesToInt32_bigEndian(c);
            int stateNum = bytesToInt32_bigEndian(c + sizeof(int)) * 2;
            size_t encodeStartIndex;
            if (nodeCount <= 256)
                encodeStartIndex = 1 + 3 * nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);
            else if (nodeCount <= 65536)
                encodeStartIndex =
                        1 + 2 * nodeCount * sizeof(unsigned short) + nodeCount * sizeof(unsigned char) +
                        nodeCount * sizeof(T);
            else
                encodeStartIndex =
                        1 + 2 * nodeCount * sizeof(unsigned int) + nodeCount * sizeof(unsigned char) +
                        nodeCount * sizeof(T);

            huffmanTree = createHuffmanTree(stateNum);
            treeRoot = reconstruct_HuffTree_from_bytes_anyStates(c + sizeof(int) + sizeof(int), nodeCount);
            c += sizeof(int) + sizeof(int) + encodeStartIndex;
            loaded = true;
        }

        bool isLoaded() { return loaded; }

    private:
        HuffmanTree *huffmanTree = NULL;
        node treeRoot;
        unsigned int nodeCount = 0;
        uchar sysEndianType; //0: little endian, 1: big endian
        bool loaded = false;
        T offset;


        node reconstruct_HuffTree_from_bytes_anyStates(const unsigned char *bytes, uint nodeCount) {
            if (nodeCount <= 256) {
                unsigned char *L = (unsigned char *) malloc(nodeCount * sizeof(unsigned char));
                memset(L, 0, nodeCount * sizeof(unsigned char));
                unsigned char *R = (unsigned char *) malloc(nodeCount * sizeof(unsigned char));
                memset(R, 0, nodeCount * sizeof(unsigned char));
                T *C = (T *) malloc(nodeCount * sizeof(T));
                memset(C, 0, nodeCount * sizeof(T));
                unsigned char *t = (unsigned char *) malloc(nodeCount * sizeof(unsigned char));
                memset(t, 0, nodeCount * sizeof(unsigned char));
                // TODO: Endian type
                // unsigned char cmpSysEndianType = bytes[0];
                // if(cmpSysEndianType!=(unsigned char)sysEndianType)
                // {
                // 	unsigned char* p = (unsigned char*)(bytes+1+2*nodeCount*sizeof(unsigned char));
                // 	size_t i = 0, size = nodeCount*sizeof(unsigned int);
                // 	while(1)
                // 	{
                // 		symTransform_4bytes(p);
                // 		i+=sizeof(unsigned int);
                // 		if(i<size)
                // 			p+=sizeof(unsigned int);
                // 		else
                // 			break;
                // 	}
                // }
                memcpy(L, bytes + 1, nodeCount * sizeof(unsigned char));
                memcpy(R, bytes + 1 + nodeCount * sizeof(unsigned char), nodeCount * sizeof(unsigned char));
                memcpy(C, bytes + 1 + 2 * nodeCount * sizeof(unsigned char), nodeCount * sizeof(T));
                memcpy(t, bytes + 1 + 2 * nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T),
                       nodeCount * sizeof(unsigned char));
                node root = this->new_node2(C[0], t[0]);
                this->unpad_tree<uchar>(L, R, C, t, 0, root);
                free(L);
                free(R);
                free(C);
                free(t);
                return root;
            } else if (nodeCount <= 65536) {
                unsigned short *L = (unsigned short *) malloc(nodeCount * sizeof(unsigned short));
                memset(L, 0, nodeCount * sizeof(unsigned short));
                unsigned short *R = (unsigned short *) malloc(nodeCount * sizeof(unsigned short));
                memset(R, 0, nodeCount * sizeof(unsigned short));
                T *C = (T *) malloc(nodeCount * sizeof(T));
                memset(C, 0, nodeCount * sizeof(T));
                unsigned char *t = (unsigned char *) malloc(nodeCount * sizeof(unsigned char));
                memset(t, 0, nodeCount * sizeof(unsigned char));

                // TODO: Endian type
                // unsigned char cmpSysEndianType = bytes[0];
                // if(cmpSysEndianType!=(unsigned char)sysEndianType)
                // {
                // 	unsigned char* p = (unsigned char*)(bytes+1);
                // 	size_t i = 0, size = 3*nodeCount*sizeof(unsigned int);
                // 	while(1)
                // 	{
                // 		symTransform_4bytes(p);
                // 		i+=sizeof(unsigned int);
                // 		if(i<size)
                // 			p+=sizeof(unsigned int);
                // 		else
                // 			break;
                // 	}
                // }

                memcpy(L, bytes + 1, nodeCount * sizeof(unsigned short));
                memcpy(R, bytes + 1 + nodeCount * sizeof(unsigned short), nodeCount * sizeof(unsigned short));
                memcpy(C, bytes + 1 + 2 * nodeCount * sizeof(unsigned short), nodeCount * sizeof(T));

                memcpy(t, bytes + 1 + 2 * nodeCount * sizeof(unsigned short) + nodeCount * sizeof(T),
                       nodeCount * sizeof(unsigned char));

                node root = this->new_node2(0, 0);
                this->unpad_tree<unsigned short>(L, R, C, t, 0, root);
                free(L);
                free(R);
                free(C);
                free(t);
                return root;
            } else //nodeCount>65536
            {
                unsigned int *L = (unsigned int *) malloc(nodeCount * sizeof(unsigned int));
                memset(L, 0, nodeCount * sizeof(unsigned int));
                unsigned int *R = (unsigned int *) malloc(nodeCount * sizeof(unsigned int));
                memset(R, 0, nodeCount * sizeof(unsigned int));
                T *C = (T *) malloc(nodeCount * sizeof(T));
                memset(C, 0, nodeCount * sizeof(T));
                unsigned char *t = (unsigned char *) malloc(nodeCount * sizeof(unsigned char));
                memset(t, 0, nodeCount * sizeof(unsigned char));
                // TODO: Endian type
                // unsigned char cmpSysEndianType = bytes[0];
                // if(cmpSysEndianType!=(unsigned char)sysEndianType)
                // {
                // 	unsigned char* p = (unsigned char*)(bytes+1);
                // 	size_t i = 0, size = 3*nodeCount*sizeof(unsigned int);
                // 	while(1)
                // 	{
                // 		symTransform_4bytes(p);
                // 		i+=sizeof(unsigned int);
                // 		if(i<size)
                // 			p+=sizeof(unsigned int);
                // 		else
                // 			break;
                // 	}
                // }

                memcpy(L, bytes + 1, nodeCount * sizeof(unsigned int));
                memcpy(R, bytes + 1 + nodeCount * sizeof(unsigned int), nodeCount * sizeof(unsigned int));
                memcpy(C, bytes + 1 + 2 * nodeCount * sizeof(unsigned int), nodeCount * sizeof(T));

                memcpy(t, bytes + 1 + 2 * nodeCount * sizeof(unsigned int) + nodeCount * sizeof(T),
                       nodeCount * sizeof(unsigned char));

                node root = this->new_node2(0, 0);
                this->unpad_tree<unsigned int>(L, R, C, t, 0, root);
                free(L);
                free(R);
                free(C);
                free(t);
                return root;
            }
        }

        node new_node(size_t freq, T c, node a, node b) {
            node n = huffmanTree->pool + huffmanTree->n_nodes++;
            if (freq) {
                n->c = c;
                n->freq = freq;
                n->t = 1;
            } else {
                n->left = a;
                n->right = b;
                n->freq = a->freq + b->freq;
                n->t = 0;
                //n->c = 0;
            }
            return n;
        }

        node new_node2(T c, unsigned char t) {
            huffmanTree->pool[huffmanTree->n_nodes].c = c;
            huffmanTree->pool[huffmanTree->n_nodes].t = t;
            return huffmanTree->pool + huffmanTree->n_nodes++;
        }

        /* priority queue */
        void qinsert(node n) {
            int j, i = huffmanTree->qend++;
            while ((j = (i >> 1)))  //j=i/2
            {
                if (huffmanTree->qq[j]->freq <= n->freq) break;
                huffmanTree->qq[i] = huffmanTree->qq[j], i = j;
            }
            huffmanTree->qq[i] = n;
        }

        node qremove() {
            int i, l;
            node n = huffmanTree->qq[i = 1];
            node p;
            if (huffmanTree->qend < 2) return 0;
            huffmanTree->qend--;
            huffmanTree->qq[i] = huffmanTree->qq[huffmanTree->qend];

            while ((l = (i << 1)) < huffmanTree->qend) {  //l=(i*2)
                if (l + 1 < huffmanTree->qend && huffmanTree->qq[l + 1]->freq < huffmanTree->qq[l]->freq) l++;
                if (huffmanTree->qq[i]->freq > huffmanTree->qq[l]->freq) {
                    p = huffmanTree->qq[i];
                    huffmanTree->qq[i] = huffmanTree->qq[l];
                    huffmanTree->qq[l] = p;
                    i = l;
                } else {
                    break;
                }
            }
            return n;
        }

        /* walk the tree and put 0s and 1s */
        /**
         * @out1 should be set to 0.
         * @out2 should be 0 as well.
         * @index: the index of the byte
         * */
        void build_code(node n, int len, uint64_t out1, uint64_t out2) {
            if (n->t) {
                huffmanTree->code[n->c] = (uint64_t *) malloc(2 * sizeof(uint64_t));
                if (len <= 64) {
                    (huffmanTree->code[n->c])[0] = out1 << (64 - len);
                    (huffmanTree->code[n->c])[1] = out2;
                } else {
                    (huffmanTree->code[n->c])[0] = out1;
                    (huffmanTree->code[n->c])[1] = out2 << (128 - len);
                }
                huffmanTree->cout[n->c] = (unsigned char) len;
                return;
            }
            int index = len >> 6; //=len/64
            if (index == 0) {
                out1 = out1 << 1;
                out1 = out1 | 0;
                build_code(n->left, len + 1, out1, 0);
                out1 = out1 | 1;
                build_code(n->right, len + 1, out1, 0);
            } else {
                if (len % 64 != 0)
                    out2 = out2 << 1;
                out2 = out2 | 0;
                build_code(n->left, len + 1, out1, out2);
                out2 = out2 | 1;
                build_code(n->right, len + 1, out1, out2);
            }
        }

        /**
         * Compute the frequency of the data and build the Huffman tree
         * @param HuffmanTree* huffmanTree (output)
         * @param int *s (input)
         * @param size_t length (input)
         * */
        void init(const T *s, size_t length) {
            T max = s[0];
            offset = s[0]; //offset is min

            ska::unordered_map<T, size_t> frequency;
            for (size_t i = 0; i < length; i++) {
                frequency[s[i]]++;
            }

            for (const auto &kv: frequency) {
                auto k = kv.first;
                if (k > max) {
                    max = k;
                }
                if (k < offset) {
                    offset = k;
                }
            }

            int stateNum = max - offset + 2;
            huffmanTree = createHuffmanTree(stateNum);

            for (const auto &f: frequency) {
                qinsert(new_node(f.second, f.first - offset, 0, 0));
            }

            while (huffmanTree->qend > 2)
                qinsert(new_node(0, 0, qremove(), qremove()));

            build_code(huffmanTree->qq[1], 0, 0, 0);
            treeRoot = huffmanTree->qq[1];

        }

        template<class T1>
        void pad_tree(T1 *L, T1 *R, T *C, unsigned char *t, unsigned int i, node root) {
            C[i] = root->c;
            t[i] = root->t;
            node lroot = root->left;
            if (lroot != 0) {
                huffmanTree->n_inode++;
                L[i] = huffmanTree->n_inode;
                pad_tree(L, R, C, t, huffmanTree->n_inode, lroot);
            }
            node rroot = root->right;
            if (rroot != 0) {
                huffmanTree->n_inode++;
                R[i] = huffmanTree->n_inode;
                pad_tree(L, R, C, t, huffmanTree->n_inode, rroot);
            }
        }

        template<class T1>
        void unpad_tree(T1 *L, T1 *R, T *C, unsigned char *t, unsigned int i, node root) {
            //root->c = C[i];
            if (root->t == 0) {
                T1 l, r;
                l = L[i];
                if (l != 0) {
                    node lroot = new_node2(C[l], t[l]);
                    root->left = lroot;
                    unpad_tree(L, R, C, t, l, lroot);
                }
                r = R[i];
                if (r != 0) {
                    node rroot = new_node2(C[r], t[r]);
                    root->right = rroot;
                    unpad_tree(L, R, C, t, r, rroot);
                }
            }
        }

        template<class T1>
        unsigned int convert_HuffTree_to_bytes_anyStates(unsigned int nodeCount, unsigned char *out) {
            T1 *L = (T1 *) malloc(nodeCount * sizeof(T1));
            memset(L, 0, nodeCount * sizeof(T1));
            T1 *R = (T1 *) malloc(nodeCount * sizeof(T1));
            memset(R, 0, nodeCount * sizeof(T1));
            T *C = (T *) malloc(nodeCount * sizeof(T));
            memset(C, 0, nodeCount * sizeof(T));
            unsigned char *t = (unsigned char *) malloc(nodeCount * sizeof(unsigned char));
            memset(t, 0, nodeCount * sizeof(unsigned char));

            pad_tree(L, R, C, t, 0, huffmanTree->qq[1]);

            unsigned int totalSize =
                    1 + 2 * nodeCount * sizeof(T1) + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);
            //*out = (unsigned char*)malloc(totalSize);
            out[0] = (unsigned char) sysEndianType;
            memcpy(out + 1, L, nodeCount * sizeof(T1));
            memcpy(out + 1 + nodeCount * sizeof(T1), R, nodeCount * sizeof(T1));
            memcpy(out + 1 + 2 * nodeCount * sizeof(T1), C, nodeCount * sizeof(T));
            memcpy(out + 1 + 2 * nodeCount * sizeof(T1) + nodeCount * sizeof(T), t, nodeCount * sizeof(unsigned char));

            free(L);
            free(R);
            free(C);
            free(t);
            return totalSize;
        }

        void SZ_FreeHuffman() {
            if (huffmanTree != NULL) {
                size_t i;
                free(huffmanTree->pool);
                huffmanTree->pool = NULL;
                free(huffmanTree->qqq);
                huffmanTree->qqq = NULL;
                for (i = 0; i < huffmanTree->stateNum; i++) {
                    if (huffmanTree->code[i] != NULL)
                        free(huffmanTree->code[i]);
                }
                free(huffmanTree->code);
                huffmanTree->code = NULL;
                free(huffmanTree->cout);
                huffmanTree->cout = NULL;
                free(huffmanTree);
                huffmanTree = NULL;
            }
        }

    };
}

#endif
