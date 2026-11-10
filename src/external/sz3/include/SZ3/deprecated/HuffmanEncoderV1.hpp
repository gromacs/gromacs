#ifndef SZ3_HUFFMAN_ENCODER_V1_HPP
#define SZ3_HUFFMAN_ENCODER_V1_HPP

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/// Decodes the Huffman streams of data version 3.3.2 only.
template <class T>
class HuffmanEncoderV1 : public concepts::EncoderInterface<T> {
   public:
    typedef struct node_t {
        struct node_t *left, *right;
        char t;  // in_node:0; otherwise:1
        T c;
    } *node;

    typedef struct HuffmanTree {
        struct node_t *pool;
        int n_nodes;
    } HuffmanTree;

    HuffmanEncoderV1() = default;

    ~HuffmanEncoderV1() override { SZ_FreeHuffman(); }

    // Declared because the destructor above suppresses the implicit ones. Copying is shallow and
    // has to stay -- the SZAlgo* wirings pass encoders by value -- so do not copy one that has
    // built its tree: both copies would free huffmanTree. Today every copy happens before that.
    HuffmanEncoderV1(const HuffmanEncoderV1 &) = default;
    HuffmanEncoderV1 &operator=(const HuffmanEncoderV1 &) = default;
    HuffmanEncoderV1(HuffmanEncoderV1 &&) noexcept = default;
    HuffmanEncoderV1 &operator=(HuffmanEncoderV1 &&) noexcept = default;

    HuffmanTree *createHuffmanTree(int stateNum) {
        HuffmanTree *tree = static_cast<HuffmanTree *>(malloc(sizeof(HuffmanTree)));
        tree->pool = static_cast<struct node_t *>(malloc(4 * static_cast<size_t>(stateNum) * sizeof(struct node_t)));
        memset(tree->pool, 0, 4 * static_cast<size_t>(stateNum) * sizeof(struct node_t));
        tree->n_nodes = 0;
        return tree;
    }

    void preprocess_encode(const std::vector<T> &, int) override { throw std::logic_error(kDecodeOnly); }

    void save(uchar *&) override { throw std::logic_error(kDecodeOnly); }

    size_t size_est() override { throw std::logic_error(kDecodeOnly); }

    static size_t size_bound(size_t, size_t) { throw std::logic_error(kDecodeOnly); }

    size_t encode(const std::vector<T> &, uchar *&) override { throw std::logic_error(kDecodeOnly); }

    void postprocess_encode() override { throw std::logic_error(kDecodeOnly); }

    void preprocess_decode() override {}

    // perform decoding
    std::vector<T> decode(const uchar *&bytes, size_t targetLength, size_t &remaining_length) override {
        node t = treeRoot;
        std::vector<T> out(targetLength);
        size_t i = 0, byteIndex = 0, count = 0;
        int r;
        node n = treeRoot;
        uint64_t encodedLength = 0;
        read(encodedLength, bytes, remaining_length);
        if (n->t)  // root->t==1 means that all state values are the same (constant)
        {
            for (count = 0; count < targetLength; count++) out[count] = n->c + offset;
            return out;
        }

        if (encodedLength > remaining_length)
            throw std::out_of_range("SZ3 Huffman: encoded length exceeds compressed buffer");

        // Walk at most the bits the stream holds, and stop once targetLength symbols are out.
        const size_t maxBits = targetLength > 0 ? encodedLength * 8 : 0;
        for (i = 0; i < maxBits; i++) {
            byteIndex = i >> 3;  // i/8
            r = i % 8;
            if (((bytes[byteIndex] >> (7 - r)) & 0x01) == 0)
                n = n->left;
            else
                n = n->right;

            if (n->t) {
                out[count] = n->c + offset;
                n = t;
                if (++count == targetLength) break;
            }
        }
        if (count < targetLength) throw std::out_of_range("SZ3 Huffman: corrupted encoded stream");
        bytes += encodedLength;
        remaining_length -= encodedLength;
        return out;
    }

    void postprocess_decode() override { SZ_FreeHuffman(); }

    // load Huffman tree
    void load(const uchar *&c, size_t &remaining_length) override {
        read(offset, c, remaining_length);
        if (remaining_length < 2 * sizeof(int)) throw std::out_of_range("SZ3 Huffman: truncated tree header");
        nodeCount = bytesToInt32_bigEndian(c);
        // The stored state count is skipped: it sizes the encode-side code tables, which decoding never
        // touches. nodeCount is bounded before it sizes anything, or the encodeStartIndex arithmetic
        // below overflows and the tree is read past the buffer.
        if (nodeCount <= 0 || static_cast<size_t>(nodeCount) > remaining_length)
            throw std::out_of_range("SZ3 Huffman: invalid node count");
        size_t encodeStartIndex;
        if (nodeCount <= 256)
            encodeStartIndex = 1 + 3 * nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);
        else if (nodeCount <= 65536)
            encodeStartIndex =
                1 + 2 * nodeCount * sizeof(unsigned short) + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);
        else
            encodeStartIndex =
                1 + 2 * nodeCount * sizeof(unsigned int) + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);

        size_t tree_bytes = sizeof(int) + sizeof(int) + encodeStartIndex;
        if (tree_bytes > remaining_length) throw std::out_of_range("SZ3 Huffman: tree exceeds compressed buffer");

        // The pool is 4x what is asked for, and unpad_tree builds each of the nodeCount nodes once.
        huffmanTree = createHuffmanTree(nodeCount);
        treeRoot = reconstruct_HuffTree_from_bytes_anyStates(c + sizeof(int) + sizeof(int), nodeCount);
        c += tree_bytes;
        remaining_length -= tree_bytes;
    }

   private:
    static constexpr const char *kDecodeOnly = "SZ3 Huffman: HuffmanEncoderV1 only decodes data version 3.3.2";
    HuffmanTree *huffmanTree = nullptr;
    node treeRoot;
    unsigned int nodeCount = 0;
    T offset;

    node reconstruct_HuffTree_from_bytes_anyStates(const unsigned char *bytes, uint nodeCount_) {
        if (nodeCount_ <= 256) {
            std::vector<unsigned char> L(nodeCount_);
            std::vector<unsigned char> R(nodeCount_);
            std::vector<T> C(nodeCount_);
            std::vector<unsigned char> t(nodeCount_);
            memcpy(L.data(), bytes + 1, nodeCount_ * sizeof(unsigned char));
            memcpy(R.data(), bytes + 1 + nodeCount_ * sizeof(unsigned char), nodeCount_ * sizeof(unsigned char));
            memcpy(C.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned char), nodeCount_ * sizeof(T));
            memcpy(t.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned char) + nodeCount_ * sizeof(T),
                   nodeCount_ * sizeof(unsigned char));
            std::vector<bool> seen(nodeCount_, false);
            seen[0] = true;
            node root = this->new_node2(C[0], t[0]);
            this->unpad_tree<uchar>(L.data(), R.data(), C.data(), t.data(), 0, root, nodeCount_, seen);
            return root;
        } else if (nodeCount_ <= 65536) {
            std::vector<unsigned short> L(nodeCount_);
            std::vector<unsigned short> R(nodeCount_);
            std::vector<T> C(nodeCount_);
            std::vector<unsigned char> t(nodeCount_);
            memcpy(L.data(), bytes + 1, nodeCount_ * sizeof(unsigned short));
            memcpy(R.data(), bytes + 1 + nodeCount_ * sizeof(unsigned short), nodeCount_ * sizeof(unsigned short));
            memcpy(C.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned short), nodeCount_ * sizeof(T));
            memcpy(t.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned short) + nodeCount_ * sizeof(T),
                   nodeCount_ * sizeof(unsigned char));
            std::vector<bool> seen(nodeCount_, false);
            seen[0] = true;
            node root = this->new_node2(0, 0);
            this->unpad_tree<unsigned short>(L.data(), R.data(), C.data(), t.data(), 0, root, nodeCount_, seen);
            return root;
        } else  // nodeCount_>65536
        {
            std::vector<unsigned int> L(nodeCount_);
            std::vector<unsigned int> R(nodeCount_);
            std::vector<T> C(nodeCount_);
            std::vector<unsigned char> t(nodeCount_);
            memcpy(L.data(), bytes + 1, nodeCount_ * sizeof(unsigned int));
            memcpy(R.data(), bytes + 1 + nodeCount_ * sizeof(unsigned int), nodeCount_ * sizeof(unsigned int));
            memcpy(C.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned int), nodeCount_ * sizeof(T));
            memcpy(t.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned int) + nodeCount_ * sizeof(T),
                   nodeCount_ * sizeof(unsigned char));
            std::vector<bool> seen(nodeCount_, false);
            seen[0] = true;
            node root = this->new_node2(0, 0);
            this->unpad_tree<unsigned int>(L.data(), R.data(), C.data(), t.data(), 0, root, nodeCount_, seen);
            return root;
        }
    }

    node new_node2(T c, unsigned char t) {
        huffmanTree->pool[huffmanTree->n_nodes].c = c;
        huffmanTree->pool[huffmanTree->n_nodes].t = t;
        return huffmanTree->pool + huffmanTree->n_nodes++;
    }

    template <class T1>
    void unpad_tree(T1 *L, T1 *R, T *C, unsigned char *t, unsigned int i, node root, unsigned int nodeCount_,
                    std::vector<bool> &seen) {
        if (root->t == 0) {
            T1 l, r;
            l = L[i];
            if (l != 0) {
                // The encoder gave a child a higher index than its parent, so a valid index satisfies i < l <
                // nodeCount_. Enforcing it keeps L/R/C/t reads inside the pool and rules out a cycle.
                if (l <= i || l >= nodeCount_) throw std::out_of_range("SZ3 Huffman: invalid left child index in tree");
                // Increasing indices rule out a cycle but not two parents naming one child, which would
                // expand the tree exponentially instead of building nodeCount_ nodes.
                if (seen[l]) throw std::out_of_range("SZ3 Huffman: tree node reached twice");
                seen[l] = true;
                node lroot = new_node2(C[l], t[l]);
                root->left = lroot;
                unpad_tree(L, R, C, t, l, lroot, nodeCount_, seen);
            }
            r = R[i];
            if (r != 0) {
                if (r <= i || r >= nodeCount_)
                    throw std::out_of_range("SZ3 Huffman: invalid right child index in tree");
                if (seen[r]) throw std::out_of_range("SZ3 Huffman: tree node reached twice");
                seen[r] = true;
                node rroot = new_node2(C[r], t[r]);
                root->right = rroot;
                unpad_tree(L, R, C, t, r, rroot, nodeCount_, seen);
            }
            if (root->left == nullptr || root->right == nullptr) {
                throw std::out_of_range("SZ3 Huffman: internal tree node is missing a child");
            }
        }
    }

    void SZ_FreeHuffman() {
        if (huffmanTree != nullptr) {
            free(huffmanTree->pool);
            free(huffmanTree);
            huffmanTree = nullptr;
        }
    }
};
}  // namespace SZ3

#endif
