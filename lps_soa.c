#include "lps_soa.h"

#include <stdlib.h>
#include <string.h>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#elif defined(__AVX2__)
#include <immintrin.h>
#elif defined(__SSE4_1__)
#include <smmintrin.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#endif

void lps_soa_init(struct lps_soa *s) {
    memset(s, 0, sizeof(*s));
    s->level = 1;
}

void lps_soa_free(struct lps_soa *s) {
    free(s->label);
    free(s->span);
    free(s->bit_rep);
    free(s->bit_size);
    free(s->label_next);
    free(s->span_next);
    free(s->bit_rep_next);
    free(s->bit_size_next);
    free(s->code);
    free(s->key);
    free(s->eq);
    free(s->gt);
    free(s->lt);
    lps_soa_init(s);
}

static int lps_soa_grow(void **p, size_t bytes) {
    void *q = realloc(*p, bytes ? bytes : 1);
    if (q == NULL)
        return -1;
    *p = q;
    return 0;
}

int lps_soa_reserve(struct lps_soa *s, int64_t n) {
    if (n <= s->capacity)
        return 0;
    const size_t c = (size_t)n;
    const size_t words = c / 64 + 2;
    if (lps_soa_grow((void **)&s->label, c * sizeof(uint32_t)) ||
        lps_soa_grow((void **)&s->span, c * sizeof(uint64_t)) ||
        lps_soa_grow((void **)&s->bit_rep, c * sizeof(uint64_t)) ||
        lps_soa_grow((void **)&s->bit_size, c) ||
        lps_soa_grow((void **)&s->label_next, c * sizeof(uint32_t)) ||
        lps_soa_grow((void **)&s->span_next, c * sizeof(uint64_t)) ||
        lps_soa_grow((void **)&s->bit_rep_next, c * sizeof(uint64_t)) ||
        lps_soa_grow((void **)&s->bit_size_next, c) ||
        lps_soa_grow((void **)&s->code, c + 64) ||
        lps_soa_grow((void **)&s->key, c * sizeof(uint64_t)) ||
        lps_soa_grow((void **)&s->eq, words * sizeof(uint64_t)) ||
        lps_soa_grow((void **)&s->gt, words * sizeof(uint64_t)) ||
        lps_soa_grow((void **)&s->lt, words * sizeof(uint64_t)))
        return -1;
    s->capacity = n;
    return 0;
}

/* MurmurHash3_32(key, 16, 42) of the 4 words (a, b, c, d), as core labels
 * use it (the body step rotates h1 by 15). */
static inline uint32_t lps_soa_rotl15(uint32_t x) { return (x << 15) | (x >> 17); }

static inline uint32_t lps_soa_murmur16(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    const uint32_t in[4] = {a, b, c, d};
    uint32_t h = 42;
    for (int i = 0; i < 4; i++) {
        uint32_t k = in[i] * 0xcc9e2d51u;
        k = lps_soa_rotl15(k) * 0x1b873593u;
        h = lps_soa_rotl15(h ^ k) * 5 + 0xe6546b64u;
    }
    h ^= 16u;
    h ^= h >> 16;
    h *= 0x85ebca6bu;
    h ^= h >> 13;
    h *= 0xc2b2ae35u;
    h ^= h >> 16;
    return h;
}

/* bit_size of a compressed core with code r (r < 128) */
static inline uint32_t lps_soa_code_bits(uint32_t r) {
    return r < 4 ? 2 : 32 - (uint32_t)__builtin_clz(r);
}

/* Compressed codes of k-mer cores: c[j] from cores j (left) and j + 1.
 * With x = label_j ^ label_j+1 != 0 the first differing bit is ctz(x) <
 * kmer_bits and the code is 2 * ctz(x) + that bit of the right core; equal
 * cores give 2 * kmer_bits (bit kmer_bits of 0x8000000000000000 | code is
 * 0). In the SIMD paths ctz(x) is the exponent of (float)(x & -x). */
static void lps_soa_compress_kmers(const uint32_t *label, int64_t m, uint32_t kmer_bits,
                                   uint8_t *c) {
    int64_t j = 0;
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    const uint32x4_t eqcode = vdupq_n_u32(2 * kmer_bits);
    for (; j + 8 <= m; j += 8) {
        uint16x4_t half[2];
        for (int h = 0; h < 2; h++) {
            const uint32x4_t l = vld1q_u32(label + j + 4 * h);
            const uint32x4_t r = vld1q_u32(label + j + 4 * h + 1);
            const uint32x4_t x = veorq_u32(l, r);
            const uint32x4_t low =
                vandq_u32(x, vreinterpretq_u32_s32(vnegq_s32(vreinterpretq_s32_u32(x))));
            const uint32x4_t f = vsubq_u32(
                vshrq_n_u32(vreinterpretq_u32_f32(vcvtq_f32_u32(low)), 23), vdupq_n_u32(127));
            const uint32x4_t bit = vshrq_n_u32(vtstq_u32(r, low), 31);
            const uint32x4_t v = vaddq_u32(vaddq_u32(f, f), bit);
            half[h] = vmovn_u32(vbslq_u32(vceqzq_u32(x), eqcode, v));
        }
        vst1_u8(c + j, vmovn_u16(vcombine_u16(half[0], half[1])));
    }
#elif defined(__AVX2__)
    const __m256i c127 = _mm256_set1_epi32(127), cff = _mm256_set1_epi32(0xFF);
    const __m256i eqcode = _mm256_set1_epi32((int)(2 * kmer_bits));
    const __m256i one = _mm256_set1_epi32(1), zero = _mm256_setzero_si256();
    for (; j + 8 <= m; j += 8) {
        const __m256i l = _mm256_loadu_si256((const __m256i *)(label + j));
        const __m256i r = _mm256_loadu_si256((const __m256i *)(label + j + 1));
        const __m256i x = _mm256_xor_si256(l, r);
        const __m256i low = _mm256_and_si256(x, _mm256_sub_epi32(zero, x));
        /* the float is negative for bit 31: drop the sign with the mask */
        const __m256i f = _mm256_sub_epi32(
            _mm256_and_si256(
                _mm256_srli_epi32(_mm256_castps_si256(_mm256_cvtepi32_ps(low)), 23), cff),
            c127);
        const __m256i bit = _mm256_andnot_si256(
            _mm256_cmpeq_epi32(_mm256_and_si256(r, low), zero), one);
        const __m256i v = _mm256_add_epi32(_mm256_add_epi32(f, f), bit);
        const __m256i w = _mm256_blendv_epi8(v, eqcode, _mm256_cmpeq_epi32(x, zero));
        const __m128i p16 =
            _mm_packs_epi32(_mm256_castsi256_si128(w), _mm256_extracti128_si256(w, 1));
        _mm_storel_epi64((__m128i *)(c + j), _mm_packus_epi16(p16, p16));
    }
#elif defined(__SSE2__)
    const __m128i c127 = _mm_set1_epi32(127), cff = _mm_set1_epi32(0xFF);
    const __m128i eqcode = _mm_set1_epi32((int)(2 * kmer_bits));
    const __m128i one = _mm_set1_epi32(1), zero = _mm_setzero_si128();
    for (; j + 8 <= m; j += 8) {
        __m128i w[2];
        for (int h = 0; h < 2; h++) {
            const __m128i l = _mm_loadu_si128((const __m128i *)(label + j + 4 * h));
            const __m128i r = _mm_loadu_si128((const __m128i *)(label + j + 4 * h + 1));
            const __m128i x = _mm_xor_si128(l, r);
            const __m128i low = _mm_and_si128(x, _mm_sub_epi32(zero, x));
            const __m128i f = _mm_sub_epi32(
                _mm_and_si128(_mm_srli_epi32(_mm_castps_si128(_mm_cvtepi32_ps(low)), 23), cff),
                c127);
            const __m128i bit =
                _mm_andnot_si128(_mm_cmpeq_epi32(_mm_and_si128(r, low), zero), one);
            const __m128i v = _mm_add_epi32(_mm_add_epi32(f, f), bit);
            const __m128i e = _mm_cmpeq_epi32(x, zero);
            w[h] = _mm_or_si128(_mm_and_si128(e, eqcode), _mm_andnot_si128(e, v));
        }
        const __m128i p16 = _mm_packs_epi32(w[0], w[1]);
        _mm_storel_epi64((__m128i *)(c + j), _mm_packus_epi16(p16, p16));
    }
#endif
    for (; j < m; j++) {
        const uint32_t x = label[j] ^ label[j + 1];
        c[j] = (uint8_t)(x ? 2 * (uint32_t)__builtin_ctz(x) + ((label[j + 1] >> __builtin_ctz(x)) & 1)
                           : 2 * kmer_bits);
    }
}

/* Compressed codes of general cores (levels >= 2), as core_compress(). */
static void lps_soa_compress(const uint64_t *br, const uint8_t *bs, int64_t m, uint8_t *c) {
    for (int64_t j = 0; j < m; j++) {
        const uint64_t l = br[j], r = br[j + 1];
        const uint32_t lb = bs[j], rb = bs[j + 1];
        const uint64_t x = l ^ r;
        uint32_t f = x ? (uint32_t)__builtin_ctzll(x) : rb;
        f = f < lb ? f : lb;
        f = f < rb ? f : rb;
        c[j] = (uint8_t)(2 * f + ((r >> f) & 1));
    }
}

/* Pair masks of c (readable up to c[64 * words]): bit i of eq/gt/lt is
 * c[i] ==/>/< c[i + 1]. Codes are < 128, so signed byte compares work. */
static void lps_soa_pair_masks(const uint8_t *c, int64_t words, uint64_t *eq, uint64_t *gt,
                               uint64_t *lt) {
    for (int64_t w = 0; w < words; w++) {
        const uint8_t *p = c + 64 * w;
        uint64_t e = 0, g = 0, l = 0;
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
        static const int8_t shift[16] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7};
        const int8x16_t sh = vld1q_s8(shift);
        uint8x16_t me[4], mg[4], ml[4];
        for (int q = 0; q < 4; q++) {
            const uint8x16_t a = vld1q_u8(p + 16 * q), b = vld1q_u8(p + 16 * q + 1);
            me[q] = vshlq_u8(vshrq_n_u8(vceqq_u8(a, b), 7), sh);
            mg[q] = vshlq_u8(vshrq_n_u8(vcgtq_u8(a, b), 7), sh);
            ml[q] = vshlq_u8(vshrq_n_u8(vcltq_u8(a, b), 7), sh);
        }
#define LPS_SOA_PACK64(v)                                                                  \
    vgetq_lane_u64(vreinterpretq_u64_u8(vpaddq_u8(                                         \
                       vpaddq_u8(vpaddq_u8(v[0], v[1]), vpaddq_u8(v[2], v[3])), vdupq_n_u8(0))), \
                   0)
        e = LPS_SOA_PACK64(me);
        g = LPS_SOA_PACK64(mg);
        l = LPS_SOA_PACK64(ml);
#undef LPS_SOA_PACK64
#elif defined(__AVX2__)
        for (int h = 0; h < 2; h++) {
            const __m256i a = _mm256_loadu_si256((const __m256i *)(p + 32 * h));
            const __m256i b = _mm256_loadu_si256((const __m256i *)(p + 32 * h + 1));
            e |= (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b)) << (32 * h);
            g |= (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpgt_epi8(a, b)) << (32 * h);
            l |= (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpgt_epi8(b, a)) << (32 * h);
        }
#elif defined(__SSE2__)
        for (int q = 0; q < 4; q++) {
            const __m128i a = _mm_loadu_si128((const __m128i *)(p + 16 * q));
            const __m128i b = _mm_loadu_si128((const __m128i *)(p + 16 * q + 1));
            e |= (uint64_t)(uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(a, b)) << (16 * q);
            g |= (uint64_t)(uint32_t)_mm_movemask_epi8(_mm_cmpgt_epi8(a, b)) << (16 * q);
            l |= (uint64_t)(uint32_t)_mm_movemask_epi8(_mm_cmpgt_epi8(b, a)) << (16 * q);
        }
#else
        for (int i = 0; i < 64; i++) {
            e |= (uint64_t)(p[i] == p[i + 1]) << i;
            g |= (uint64_t)(p[i] > p[i + 1]) << i;
            l |= (uint64_t)(p[i] < p[i + 1]) << i;
        }
#endif
        eq[w] = e;
        gt[w] = g;
        lt[w] = l;
    }
}

/* out[i] = label of the core over compressed [b, e), key[i] = b | e << 32:
 * murmur16(lab[b + 1], lab[e - 1], lab[e], e - b - 2). Lanes are
 * independent, so SIMD hides the multiply latency of one hash. */
static void lps_soa_labels(const uint64_t *key, const uint32_t *lab, uint32_t *out, int64_t n) {
    int64_t i = 0;
#if defined(__AVX2__) && !(defined(__ARM_NEON) || defined(__ARM_NEON__))
    const __m256i c1 = _mm256_set1_epi32((int)0xcc9e2d51u);
    const __m256i c2 = _mm256_set1_epi32(0x1b873593);
    const __m256i c5 = _mm256_set1_epi32(5);
    const __m256i ce = _mm256_set1_epi32((int)0xe6546b64u);
    const __m256i f1 = _mm256_set1_epi32((int)0x85ebca6bu);
    const __m256i f2 = _mm256_set1_epi32((int)0xc2b2ae35u);
    const __m256i one = _mm256_set1_epi32(1), two = _mm256_set1_epi32(2);
    const __m256i split = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);
    const int *L = (const int *)lab;
    for (; i + 8 <= n; i += 8) {
        const __m256i k0 =
            _mm256_permutevar8x32_epi32(_mm256_loadu_si256((const __m256i *)(key + i)), split);
        const __m256i k1 = _mm256_permutevar8x32_epi32(
            _mm256_loadu_si256((const __m256i *)(key + i + 4)), split);
        const __m256i B = _mm256_permute2x128_si256(k0, k1, 0x20);
        const __m256i E = _mm256_permute2x128_si256(k0, k1, 0x31);
        __m256i in[4];
        in[0] = _mm256_i32gather_epi32(L, _mm256_add_epi32(B, one), 4);
        in[1] = _mm256_i32gather_epi32(L, _mm256_sub_epi32(E, one), 4);
        in[2] = _mm256_i32gather_epi32(L, E, 4);
        in[3] = _mm256_sub_epi32(_mm256_sub_epi32(E, B), two);
        __m256i h = _mm256_set1_epi32(42);
        for (int q = 0; q < 4; q++) {
            __m256i k = _mm256_mullo_epi32(in[q], c1);
            k = _mm256_or_si256(_mm256_slli_epi32(k, 15), _mm256_srli_epi32(k, 17));
            k = _mm256_mullo_epi32(k, c2);
            h = _mm256_xor_si256(h, k);
            h = _mm256_or_si256(_mm256_slli_epi32(h, 15), _mm256_srli_epi32(h, 17));
            h = _mm256_add_epi32(_mm256_mullo_epi32(h, c5), ce);
        }
        h = _mm256_xor_si256(h, _mm256_set1_epi32(16));
        h = _mm256_xor_si256(h, _mm256_srli_epi32(h, 16));
        h = _mm256_mullo_epi32(h, f1);
        h = _mm256_xor_si256(h, _mm256_srli_epi32(h, 13));
        h = _mm256_mullo_epi32(h, f2);
        h = _mm256_xor_si256(h, _mm256_srli_epi32(h, 16));
        _mm256_storeu_si256((__m256i *)(out + i), h);
    }
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__SSE4_1__)
    /* 4 lanes; the inputs are loaded one by one (no gather). */
    for (; i + 4 <= n; i += 4) {
        uint32_t w[4][4];
        for (int t = 0; t < 4; t++) {
            const uint32_t b = (uint32_t)key[i + t], e = (uint32_t)(key[i + t] >> 32);
            w[0][t] = lab[b + 1];
            w[1][t] = lab[e - 1];
            w[2][t] = lab[e];
            w[3][t] = e - b - 2;
        }
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
        uint32x4_t h = vdupq_n_u32(42);
        for (int q = 0; q < 4; q++) {
            uint32x4_t k = vmulq_n_u32(vld1q_u32(w[q]), 0xcc9e2d51u);
            k = vsriq_n_u32(vshlq_n_u32(k, 15), k, 17);
            k = vmulq_n_u32(k, 0x1b873593u);
            h = veorq_u32(h, k);
            h = vsriq_n_u32(vshlq_n_u32(h, 15), h, 17);
            h = vaddq_u32(vmulq_n_u32(h, 5), vdupq_n_u32(0xe6546b64u));
        }
        h = veorq_u32(h, vdupq_n_u32(16));
        h = veorq_u32(h, vshrq_n_u32(h, 16));
        h = vmulq_n_u32(h, 0x85ebca6bu);
        h = veorq_u32(h, vshrq_n_u32(h, 13));
        h = vmulq_n_u32(h, 0xc2b2ae35u);
        h = veorq_u32(h, vshrq_n_u32(h, 16));
        vst1q_u32(out + i, h);
#else
        __m128i h = _mm_set1_epi32(42);
        for (int q = 0; q < 4; q++) {
            __m128i k = _mm_mullo_epi32(_mm_loadu_si128((const __m128i *)w[q]),
                                        _mm_set1_epi32((int)0xcc9e2d51u));
            k = _mm_or_si128(_mm_slli_epi32(k, 15), _mm_srli_epi32(k, 17));
            k = _mm_mullo_epi32(k, _mm_set1_epi32(0x1b873593));
            h = _mm_xor_si128(h, k);
            h = _mm_or_si128(_mm_slli_epi32(h, 15), _mm_srli_epi32(h, 17));
            h = _mm_add_epi32(_mm_mullo_epi32(h, _mm_set1_epi32(5)),
                              _mm_set1_epi32((int)0xe6546b64u));
        }
        h = _mm_xor_si128(h, _mm_set1_epi32(16));
        h = _mm_xor_si128(h, _mm_srli_epi32(h, 16));
        h = _mm_mullo_epi32(h, _mm_set1_epi32((int)0x85ebca6bu));
        h = _mm_xor_si128(h, _mm_srli_epi32(h, 13));
        h = _mm_mullo_epi32(h, _mm_set1_epi32((int)0xc2b2ae35u));
        h = _mm_xor_si128(h, _mm_srli_epi32(h, 16));
        _mm_storeu_si128((__m128i *)(out + i), h);
#endif
    }
#endif
    for (; i < n; i++) {
        const uint32_t b = (uint32_t)key[i], e = (uint32_t)(key[i] >> 32);
        out[i] = lps_soa_murmur16(lab[b + 1], lab[e - 1], lab[e], e - b - 2);
    }
}

/* Bits i + s of mask v at word w (s in {-1, 1, 2}); v[words] is 0. */
static inline uint64_t lps_soa_at(const uint64_t *v, int64_t w, int s) {
    if (s > 0)
        return (v[w] >> s) | (v[w + 1] << (64 - s));
    return (v[w] << 1) | (w ? v[w - 1] >> 63 : 0);
}

/* Bits of word w (positions 64w..64w+63) below position lim. */
static inline uint64_t lps_soa_below(int64_t w, int64_t lim) {
    const int64_t lo = 64 * w;
    if (lim <= lo)
        return 0;
    if (lim - lo >= 64)
        return ~(uint64_t)0;
    return ((uint64_t)1 << (lim - lo)) - 1;
}

/* One level: n >= 2 cores in the current arrays -> the next level's cores in
 * the *_next arrays; returns their count. Compressed core j is core j + 1
 * compressed against core j: its label and end are core j + 1's, its start
 * is core j's. The rules follow parse3() of lps.c. */
static int64_t lps_soa_level(struct lps_soa *s, int kmer_bits) {
    const int64_t m = s->size - 1;
    uint8_t *c = s->code;
    if (s->level == 1)
        lps_soa_compress_kmers(s->label, m, (uint32_t)kmer_bits, c);
    else
        lps_soa_compress(s->bit_rep, s->bit_size, m, c);
    memset(c + m, 0, 64);

    const int64_t words = (m + 63) / 64;
    uint64_t *E = s->eq, *G = s->gt, *Lt = s->lt;
    lps_soa_pair_masks(c, words, E, G, Lt);
    E[words] = G[words] = Lt[words] = 0;
    /* pairs (i, i + 1) exist for i <= m - 2 */
    const int64_t npairs = m > 1 ? m - 1 : 0;
    for (int64_t w = npairs / 64; w < words; w++) {
        const uint64_t keep = lps_soa_below(w, npairs);
        E[w] &= keep;
        G[w] &= keep;
        Lt[w] &= keep;
    }

    const uint64_t *span = s->span;
    uint64_t *obr = s->bit_rep_next, *ospan = s->span_next, *key = s->key;
    uint8_t *obs = s->bit_size_next;
    int64_t k = 0;
    /* new core k over compressed [b, e); bit_rep keeps whole codes from the
     * right while they fit in 64 bits (at most 21 bits for 3 codes) */
#define LPS_SOA_PUT(b_, e_)                                                          \
    do {                                                                             \
        const int64_t b = (b_), e = (e_);                                            \
        ospan[k] = (span[b] & 0xFFFFFFFFu) | (span[e] & ~(uint64_t)0xFFFFFFFFu);     \
        key[k] = (uint64_t)b | (uint64_t)e << 32;                                    \
        uint32_t size = 0;                                                           \
        uint64_t rep = 0;                                                            \
        if (e - b == 3) {                                                            \
            const uint32_t s0 = lps_soa_code_bits(c[b]), s1 = lps_soa_code_bits(c[b + 1]), \
                           s2 = lps_soa_code_bits(c[b + 2]);                        \
            size = s0 + s1 + s2;                                                     \
            rep = c[b + 2] | (uint64_t)c[b + 1] << s2 | (uint64_t)c[b] << (s2 + s1); \
        } else {                                                                     \
            uint32_t index = 0;                                                      \
            int full = 0;                                                            \
            for (int64_t t = e - 1; b <= t; t--) {                                   \
                const uint32_t sz = lps_soa_code_bits(c[t]);                         \
                size += sz;                                                          \
                if (!full && index + sz <= 64) {                                     \
                    rep |= (uint64_t)c[t] << index;                                  \
                    index += sz;                                                     \
                } else {                                                             \
                    full = 1;                                                        \
                }                                                                    \
            }                                                                        \
        }                                                                            \
        obr[k] = rep & 0x7FFFFFFFFFFFFFFFull;                                        \
        obs[k] = (uint8_t)(size < 63 ? size : 63);                                   \
        k++;                                                                         \
    } while (0)

    int64_t i2 = m; /* end of the last core; m: none yet */
    for (int64_t w = 0; w < words; w++) {
        const uint64_t e0 = E[w], e1 = lps_soa_at(E, w, 1);
        const uint64_t g0 = G[w], g1 = lps_soa_at(G, w, 1), gm = lps_soa_at(G, w, -1);
        const uint64_t l0 = Lt[w], l1 = lps_soa_at(Lt, w, 1), l2 = lps_soa_at(Lt, w, 2);
        const uint64_t rint = ~e0 & e1;
        const uint64_t lmin = g0 & l1;
        uint64_t lmax = l0 & g1 & ~gm & ~l2;
        if (w == 0)
            lmax &= ~(uint64_t)1;
        /* every core needs i + 2 < m, LMAX i + 3 < m */
        uint64_t starts =
            ((rint | lmin) & lps_soa_below(w, m - 2)) | (lmax & lps_soa_below(w, m - 3));
        while (starts) {
            const int64_t i = 64 * w + __builtin_ctzll(starts);
            starts &= starts - 1;
            /* The core ends one past the first unequal pair after i + 1: the
             * end of the run for RINT, pair i + 1 itself for LMIN/LMAX. */
            int64_t pw = (i + 1) / 64;
            uint64_t v = ~E[pw] & (~(uint64_t)0 << ((i + 1) % 64));
            while (v == 0)
                v = ~E[++pw];
            const int64_t t = 64 * pw + __builtin_ctzll(v) + 1;
            if (t >= m)
                continue; /* a run reaching the end is no core */
            if (i2 < i)
                LPS_SOA_PUT(i2 - 1, i + 1); /* the gap since the last core */
            i2 = t + 1;
            LPS_SOA_PUT(i, i2);
        }
    }
#undef LPS_SOA_PUT
    lps_soa_labels(key, s->label, s->label_next, k);
    return k;
}

int lps_soa_deepen(struct lps_soa *s, int kmer_bits, int lcp_level) {
    if (lcp_level <= s->level)
        return 0;
    while (s->level < lcp_level) {
        if (s->size < 2) { /* as lcp_dct(): too few cores, none left */
            s->size = 0;
            s->level++;
            break;
        }
        s->size = lps_soa_level(s, kmer_bits);
        uint32_t *l = s->label;
        s->label = s->label_next;
        s->label_next = l;
        uint64_t *sp = s->span;
        s->span = s->span_next;
        s->span_next = sp;
        uint64_t *br = s->bit_rep;
        s->bit_rep = s->bit_rep_next;
        s->bit_rep_next = br;
        uint8_t *bs = s->bit_size;
        s->bit_size = s->bit_size_next;
        s->bit_size_next = bs;
        s->level++;
    }
    return 1;
}
