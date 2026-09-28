#include "core.h"
#include "lps.h"
#include "lps_soa.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

void log(const std::string &message) {
	std::cout << message << std::endl;
}

// lps_soa_deepen() must give exactly the cores of lps_deepen() for the same
// level-1 k-mer cores, at every level. Codes come from small alphabets too,
// so equal neighbours (RINT cores) and gaps (SSEQ cores) are frequent.
static void check_sequence(const std::vector<uint32_t> &codes, int kmer_bits, int lcp_level) {
    const int n = (int)codes.size();

    struct lps ref;
    ref.level = 1;
    ref.size = n;
    ref.capacity = n > 0 ? n : 1;
    ref.cores = (struct core *)malloc(ref.capacity * sizeof(struct core));

    struct lps_soa soa;
    lps_soa_init(&soa);
    int rc = lps_soa_reserve(&soa, n);
    assert(rc == 0 && "lps_soa_reserve should succeed");
    (void)rc;

    uint64_t pos = 0;
    for (int i = 0; i < n; i++) {
        const uint64_t start = pos, end = pos + 3 + (codes[i] % 5);
        pos += 1 + (codes[i] % 3);
        init_core4(&ref.cores[i], kmer_bits, 0x8000000000000000 | codes[i], codes[i], start, end);
        soa.label[i] = codes[i];
        soa.span[i] = start | end << 32;
    }
    soa.size = n;

    lps_deepen(&ref, lcp_level);
    lps_soa_deepen(&soa, kmer_bits, lcp_level);

    assert(ref.level == soa.level && "levels should match");
    assert(ref.size == soa.size && "core counts should match");
    for (int i = 0; i < ref.size; i++) {
        const struct core *c = &ref.cores[i];
        assert(c->label == soa.label[i] && "labels should match");
        assert(c->start == (soa.span[i] & 0xFFFFFFFFu) && "starts should match");
        assert(c->end == (soa.span[i] >> 32) && "ends should match");
        if (lcp_level > 1) {
            assert(c->bit_rep == soa.bit_rep[i] && "bit_reps should match");
            assert(c->bit_size == soa.bit_size[i] && "bit_sizes should match");
        }
    }

    free(ref.cores);
    lps_soa_free(&soa);
}

void test_lps_soa_deepen() {

    std::mt19937_64 rng(12345);
    int checked = 0;
    for (int kmer_bits : {2, 8, 22, 31, 32}) {
        const uint64_t full = kmer_bits == 32 ? 0xFFFFFFFFull : (1ull << kmer_bits) - 1;
        for (uint64_t alphabet : {2ull, 3ull, 5ull, 64ull, 0ull}) {
            for (int n : {0, 1, 2, 3, 4, 5, 7, 64, 65, 127, 128, 129, 1000, 5000}) {
                std::vector<uint32_t> codes(n);
                for (auto &c : codes) {
                    // alphabet 0: any code; else few distinct codes
                    const uint64_t r = rng();
                    c = (uint32_t)(alphabet ? (r % alphabet) * (full / alphabet) : r) & full;
                }
                for (int lcp_level : {1, 2, 3, 4, 5, 7}) {
                    check_sequence(codes, kmer_bits, lcp_level);
                    checked++;
                }
            }
        }
    }

    // deepening in steps must equal deepening at once
    std::vector<uint32_t> codes(3000);
    for (auto &c : codes)
        c = (uint32_t)(rng() % 7) * 1234567u & 0x3FFFFF;
    struct lps_soa a, b;
    lps_soa_init(&a);
    lps_soa_init(&b);
    lps_soa_reserve(&a, 3000);
    lps_soa_reserve(&b, 3000);
    for (int i = 0; i < 3000; i++) {
        a.label[i] = b.label[i] = codes[i];
        a.span[i] = b.span[i] = (uint64_t)i | (uint64_t)(i + 11) << 32;
    }
    a.size = b.size = 3000;
    lps_soa_deepen(&a, 22, 5);
    lps_soa_deepen(&b, 22, 3);
    int again = lps_soa_deepen(&b, 22, 3);
    assert(again == 0 && "deepening to the same level should do nothing");
    (void)again;
    lps_soa_deepen(&b, 22, 5);
    assert(a.size == b.size && a.level == b.level && "sizes should match");
    for (int i = 0; i < a.size; i++)
        assert(a.label[i] == b.label[i] && a.span[i] == b.span[i] && "cores should match");
    lps_soa_free(&a);
    lps_soa_free(&b);

    log("...  test_lps_soa_deepen passed! (" + std::to_string(checked) + " sequences)");
}

int main() {

	log("Running test_lps_soa...");

	test_lps_soa_deepen();

	log("All tests in test_lps_soa completed successfully!");

	return 0;
}
