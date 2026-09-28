/**
 * @file lps_soa.h
 * @brief LCP deepening on per-field arrays (structure of arrays) for cores
 * that start as k-mers.
 *
 * `lps_soa` holds the cores of one level as separate arrays instead of an
 * array of `struct core`. Level-1 cores are k-mer cores: label = the k-mer
 * code (< 2^kmer_bits), bit_rep = 0x8000000000000000 | code and bit_size =
 * kmer_bits, so only the labels and spans are stored for level 1.
 *
 * lps_soa_deepen() gives exactly the cores (label, start, end, bit_rep,
 * bit_size) that lps_deepen() gives for the same level-1 cores, but runs the
 * compression, the core rules and the label hashes over contiguous arrays,
 * with SIMD where available (NEON, AVX2, SSE2/SSE4.1). After compression a
 * core is a 7-bit code (2 * first differing bit + that bit of the right
 * core), so the core rules become byte compares.
 *
 * Positions are stored as 32-bit start and end in one 64-bit span.
 *
 * @code
 *   struct lps_soa s;
 *   lps_soa_init(&s);
 *   lps_soa_reserve(&s, n);
 *   for (i = 0; i < n; i++) {
 *       s.label[i] = code[i];
 *       s.span[i] = start[i] | (uint64_t)end[i] << 32;
 *   }
 *   s.size = n;
 *   s.level = 1;
 *   lps_soa_deepen(&s, 2 * k, 3);
 *   // s.label[0, s.size), s.span[0, s.size)
 *   lps_soa_free(&s);
 * @endcode
 */

#ifndef LPS_SOA_H
#define LPS_SOA_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

struct lps_soa {
    int level;
    int64_t size;
    int64_t capacity;
    uint32_t *label;
    uint64_t *span;     /* start | end << 32 */
    uint64_t *bit_rep;  /* levels >= 2 */
    uint8_t *bit_size;  /* levels >= 2 */
    /* next level and scratch, owned by lps_soa */
    uint32_t *label_next;
    uint64_t *span_next;
    uint64_t *bit_rep_next;
    uint8_t *bit_size_next;
    uint8_t *code;      /* compressed cores, 64 bytes of padding */
    uint64_t *key;      /* label hash inputs: first | last << 32 */
    uint64_t *eq, *gt, *lt; /* pair masks, one bit per compressed core */
};

/**
 * @brief Initializes an empty lps_soa (no allocation).
 */
void lps_soa_init(struct lps_soa *s);

/**
 * @brief Makes room for n cores; keeps the current cores.
 *
 * @return 0 on success, -1 if out of memory (s is unchanged then).
 */
int lps_soa_reserve(struct lps_soa *s, int64_t n);

/**
 * @brief Frees all arrays and resets s to empty.
 */
void lps_soa_free(struct lps_soa *s);

/**
 * @brief Deepens to lcp_level, like lps_deepen().
 *
 * If s->level is 1, the cores must be k-mer cores of kmer_bits bits
 * (1 <= kmer_bits <= 32); kmer_bits is ignored otherwise.
 *
 * @return 1 if deepening was attempted, 0 if lcp_level <= s->level.
 */
int lps_soa_deepen(struct lps_soa *s, int kmer_bits, int lcp_level);

#ifdef __cplusplus
}
#endif

#endif
