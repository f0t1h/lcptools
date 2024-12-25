#include "lps.h"

void reverse(const char *str, int len, char **rev) {
    *rev = (char*) malloc(len*sizeof(char));
    size_t left = 0;
    size_t right = len - 1;

    while (left < right) {
        (*rev)[left] = str[right];
        (*rev)[right] = str[left];

        left++;
        right--;
    }
    if (left == right) {
        (*rev)[left] = str[left];
    }
};

void init_lps(struct lps *nlps, const char *str, int len, int rev_comp) {   
    nlps->level = 1;
    nlps->size = 0;
    nlps->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));

    if (!rev_comp) {
        nlps->size = parse1(str, str+len, nlps->cores);
    } else {
        char *rev = NULL;
        reverse(str, len, &rev);
        nlps->size = parse2(rev, rev+len, nlps->cores);
        free(rev);
    }
};

void free_lps(struct lps *nlps) {
    for(int i=0; i<nlps->size; i++) {
        free(nlps->cores[i].bit_rep);
    }
    free(nlps->cores);
    nlps->size = 0;
};

int parse1(const char *begin, const char *end, struct core *cores) {

    const char *it1 = begin;
    const char *it2 = end;
    int core_index = 0;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (alphabet[(unsigned char)*it1] == alphabet[(unsigned char)*(it1+1)]) {
            continue;
        }

        // check for RINT core
        if (alphabet[(unsigned char)*(it1+1)] == alphabet[(unsigned char)*(it1+2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 + 2;
            while (temp < end && alphabet[(unsigned char)*(temp-1)] == alphabet[(unsigned char)*temp]) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1) {
                    init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1, it1-begin+1);
                    core_index++;
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                init_core1(&(cores[core_index]), it1, it2-it1, it1-begin, it2-begin);
                core_index++;

                continue;
            }
        }

        if (alphabet[(unsigned char)*it1] > alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] < alphabet[(unsigned char)*(it1+2)] ) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1, it1-begin+1);
                core_index++;
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, it2-it1, it1-begin, it2-begin);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            alphabet[(unsigned char)*it1] < alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] > alphabet[(unsigned char)*(it1+2)] &&
            alphabet[(unsigned char)*(it1-1)] <= alphabet[(unsigned char)*(it1)] &&
            alphabet[(unsigned char)*(it1+2)] >= alphabet[(unsigned char)*(it1+3)] ) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1, it1-begin+1);
                core_index++;
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, it2-it1, it1-begin, it2-begin);
            core_index++;

            continue;
        }
    }

    return core_index;
};

int parse2(const char *begin, const char *end, struct core *cores) {

    const char *it1 = begin;
    const char *it2 = end;
    int core_index = 0;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (rc_alphabet[(unsigned char)*it1] == rc_alphabet[(unsigned char)*(it1+1)]) {
            continue;
        }

        // check for RINT core
        if (rc_alphabet[(unsigned char)*(it1+1)] == rc_alphabet[(unsigned char)*(it1+2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 + 2;
            while (temp < end && rc_alphabet[(unsigned char)*(temp-1)] == rc_alphabet[(unsigned char)*temp]) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1) {
                    init_core2(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1, it1-begin+1);
                    core_index++;
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                init_core2(&(cores[core_index]), it1, it2-it1, it1-begin, it2-begin);
                core_index++;

                continue;
            }
        }

        if (rc_alphabet[(unsigned char)*it1] > rc_alphabet[(unsigned char)*(it1+1)] &&
            rc_alphabet[(unsigned char)*(it1+1)] < rc_alphabet[(unsigned char)*(it1+2)] ) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core2(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1, it1-begin+1);
                core_index++;
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core2(&(cores[core_index]), it1, it2-it1, it1-begin, it2-begin);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            rc_alphabet[(unsigned char)*it1] < rc_alphabet[(unsigned char)*(it1+1)] &&
            rc_alphabet[(unsigned char)*(it1+1)] > rc_alphabet[(unsigned char)*(it1+2)] &&
            rc_alphabet[(unsigned char)*(it1-1)] <= rc_alphabet[(unsigned char)*(it1)] &&
            rc_alphabet[(unsigned char)*(it1+2)] >= rc_alphabet[(unsigned char)*(it1+3)] ) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core2(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1, it1-begin+1);
                core_index++;
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core2(&(cores[core_index]), it1, it2-it1, it1-begin, it2-begin);
            core_index++;

            continue;
        }
    }

    return core_index;
};

int parse3(struct core *begin, struct core *end, struct core *cores) {

    struct core *it1 = begin + DCT_ITERATION_COUNT;
    struct core *it2 = end;
    int core_index = 0;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (core_eq(it1, it1+1)) {
            continue;
        }

        // check for RINT core
        if (core_eq(it1+1, it1+2)) {

            // count middle characters
            uint32_t middle_count = 1;
            struct core *temp = it1 + 2;
            while (temp < end && core_eq(temp-1, temp)) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1) {
                    init_core3(&(cores[core_index]), it2-1-DCT_ITERATION_COUNT, it1-it2+2+DCT_ITERATION_COUNT);
                    core_index++;
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                init_core3(&(cores[core_index]), it1-DCT_ITERATION_COUNT, it2-it1+DCT_ITERATION_COUNT);
                core_index++;

                continue;
            }
        }

        // check for LMIN
        if (core_gt(it1, it1+1) && core_lt(it1+1, it1+2)) {
            
            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core3(&(cores[core_index]), it2-1-DCT_ITERATION_COUNT, it1-it2+2+DCT_ITERATION_COUNT);
                core_index++;
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1-DCT_ITERATION_COUNT, it2-it1+DCT_ITERATION_COUNT);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            core_lt(it1, it1+1) &&
            core_gt(it1+1, it1+2) &&
            core_leq(it1-1, it1) &&
            core_geq(it1+2, it1+3)) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core3(&(cores[core_index]), it2-1-DCT_ITERATION_COUNT, it1-it2+2+DCT_ITERATION_COUNT);
                core_index++;
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1-DCT_ITERATION_COUNT, it2-it1+DCT_ITERATION_COUNT);
            core_index++;

            continue;
        }
    }
    return core_index;
};

int64_t lps_memsize(const struct lps *str) {
    uint64_t total = sizeof(struct lps);
    
    for(int i=0; i<str->size; i++) {
        total += core_memsize(&(str->cores[i]));
    }

    return total;
};

/**
 * @brief Performs Deterministic Coin Tossing (DCT) compression on binary sequences.
 *
 * This function is a central part of the LCP (Locally Consisted Parsing) algorithm. It identifies differences
 * between consecutive binary strings, compressing the information by focusing on the position and value of
 * the first divergent bit from the right-end of the strings. This difference is used to generate a compact
 * 'core' that encapsulates the unique elements of each sequence in a smaller binary form.
 *
 * This compression significantly reduces redundant information, making further analysis of the sequences
 * within the LCP framework more efficient and manageable.
 *
 * @return 0 if dct is performed, -1 if no enough cores are available for dct.
 */
int lcp_dct(struct lps *str) {

    // at least 2 cores are needed for compression
    if (str->size < DCT_ITERATION_COUNT + 2) {
        return -1;
    }

    for (size_t dct_index = 0; dct_index < DCT_ITERATION_COUNT; dct_index++) {
        struct core *it_left = str->cores + str->size - 2, *it_right = str->cores + str->size - 1;

        for (; str->cores + dct_index <= it_left; it_left--, it_right--) {
            core_compress(it_left, it_right);
        }
    }

    return 0;
};

int lps_deepen1(struct lps *str) {

    // compress cores
    if (lcp_dct(str) < 0) {
        return 0;
    }

    // find new cores
    int new_size = parse3(str->cores + DCT_ITERATION_COUNT, str->cores + str->size, str->cores);
    int temp = new_size;

    // remove old cores
    while(temp < str->size) {
        free(str->cores[temp].bit_rep);
        temp++;
    }
    str->size = new_size;

    str->level++;

    return 1;
};


int lps_deepen(struct lps *str, int lcp_level) {

    if (lcp_level <= str->level)
        return 0;

    while (str->level < lcp_level && lps_deepen1(str))
        ;

    return 1;
};

void print_lps(const struct lps *str) {
    printf("Level: %d \n", str->level);
    for(int i=0; i<str->size; i++) {
        print_core(&(str->cores[i]));
        printf(" ");
    }
};

int lps_eq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 1;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 1;
        }
    }

    return 0;
};

int lps_neq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 0;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 0;
        }
    }

    return 1;
};