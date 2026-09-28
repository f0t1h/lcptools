#ifndef _LCPTOOLS_HO_H_
#define _LCPTOOLS_HO_H_
/**
 * @file lps.h
 * @brief Defines the lps struct and its associated methods for handling locally consistent
 * parsing (LCP) of strings.
 *
 * The lps struct is responsible for performing LCP operations on input strings, constructing
 * cores, and supporting functionalities like parsing, compression (deepening), and memory
 * usage calculations. It includes methods for reading and writing the data to files as well
 * as deepening the LCP to higher levels of compression.
 *
 * The lps struct leverages various helper classes like core, encoding, and hash to manage the
 * string data and its compressed forms. Additionally, it supports both standard and
 * reverse-complement parsing for specialized string handling in bioinformatics and other fields.
 *
 * Key functionalities include:
 * - Parsing an input string or file to extract LCP cores.
 * - Performing multi-level compression of LCP cores (deepening).
 * - Saving and loading LCP cores from files.
 * - Calculating memory usage of the constructed LCP structure.
 *
 * Dependencies:
 * - Requires core.h, encoding.h, hash.h, and constant.h for auxiliary data structures and utilities.
 *
 * Example usage:
 * @code
 *   std::string sequence = "AGCTAGCTAG";
 *   lcp::lps parser(sequence);
 *   parser.deepen();
 *   parser.write("output.lps");
 * @endcode
 *
 * @see core.h
 * @see encoding.h
 * @see hash.h
 * @see constant.h
 *
 * @namespace lcp
 * @struct lps
 *
 * @note Destructor handles clean-up of allocated memory for cores.
 *
 * @author Akmuhammet Ashyralyyev
 * @version 1.0
 * @date 2024-09-14
 *
 */

#ifndef LPS_H
#define LPS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <math.h>

#define CONSTANT_FACTOR         1.5

struct lps {
    int level;
    int size;
    int capacity;
    struct core *cores;
};

/**
 * @brief Constructs an lps object from a string.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 */
void init_lps(struct lps *lps_ptr, const char *str, int len);

/**
 * @brief Constructs an lps object from a string.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 * @param offset The length of the offset in which each index will be shifted.
 */
void init_lps_offset(struct lps *lps_ptr, const char *str, int len, uint64_t offset);

/**
 * @brief Constructs an lps object from a string, with reverse complement
 * transformation.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 */
void init_lps2(struct lps *lps_ptr, const char *str, int len);
/**
 * @brief Initializes an lps object by reading its contents from a binary file.
 *
 * This function reads the level and size of the `lps` object from the provided file, 
 * allocates memory for the cores array if necessary, and then reads each core's data. 
 * If any error occurs during reading, it will print an error message and terminate 
 * the program or return control, depending on error-handling strategy.
 *
 * @param lps_ptr The `lps` object that will be initialized
 * @param in File pointer to the binary file containing the serialized lps data.
 * 
 * @note The caller must ensure that the `FILE *in` is a valid and open binary file 
 *       for reading. The function will allocate memory for `lps_ptr->cores` if 
 *       `lps_ptr->size > 0`. The caller is responsible for freeing this memory later.
 */
void init_lps3(struct lps *lps_ptr, FILE *in);

/**
 * @brief Constructs an lps object from a string, using split and merge paradigm.
 * The give string will be roughly split into the length of `chunk_size`
 * and processed independently. The final cores will be merged into a single array.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 * @param chunk_size The length of the chunks to be processed.
 */
void init_lps4(struct lps *lps_ptr, const char *str, int len, int lcp_level, int chunk_size);

/**
 * @brief Destructor for the lps object. Frees dynamically allocated memory for cores.
 */
void free_lps(struct lps *lps_ptr);

/**
 * @brief Serializes and writes an lps object to a binary file.
 *
 * This function writes the level, size, and all core objects of the `lps` object 
 * to the specified binary file. Each core's data, including its bit representation, 
 * is written sequentially to the file. The resulting file can later be read to 
 * reconstruct the `lps` object.
 *
 * @param lps_ptr The `lps` object that will be initialized
 * @param out File pointer to the binary file where the lps data will be written.
 * 
 * @note The caller must ensure that the `FILE *out` is a valid and open binary file 
 *       for writing. If the file cannot be written to, the function behavior is undefined.
 */
void write_lps(struct lps *lps_ptr, FILE *out);

/**
 * @brief Parses a sequence to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores.
 *
 * This function iterates over a sequence defined by iterators `begin` and `end` and identifies key
 * segments, called "cores," that represent the (LCP) regions. By analyzing
 * character relationships in the sequence (such as equality or relative order), it builds and stores
 * these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the sequence to parse.
 * @param end Iterator pointing to the end of the sequence to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @param offset The distance measure where the indecies of the core will be shifted by.
 * @return Size of the cores identified in the given string.
 */
int parse1(const char *begin, const char *end, struct core *cores, uint64_t offset);

/**
 * @brief Parses a sequence to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores using complement alphabet.
 *
 * This function iterates over a sequence defined by iterators `begin` and `end` and identifies key
 * segments, called "cores," that represent the (LCP) regions. By analyzing
 * character relationships in the sequence (such as equality or relative order based on complement), 
 * it builds and stores these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the sequence to parse.
 * @param end Iterator pointing to the end of the sequence to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @param offset The distance measure where the indecies of the core will be shifted by.
 * @return Size of the cores identified in the given string.
 */
int parse2(const char *begin, const char *end, struct core *cores, uint64_t offset);

/**
 * @brief Parses a array of cores to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores.
 *
 * This function iterates over a array of `core` structures defined by iterators `begin` and `end` and 
 * identifies key segments, called "cores," that represent the (LCP) regions. By analyzing
 * `core` structure relationships in the array (such as equality or relative order), it builds and stores
 * these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the `core` array to parse.
 * @param end Iterator pointing to the end of the `core` array to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @return Size of the cores identified in the given string.
 */
int parse3(struct core *begin, struct core *end, struct core *cores);

/**
 * @brief Calculates and returns the memory size used by the `lps` structure.
 *
 * @return The memory size (in bytes) used by the `lps` structure.
 */
int64_t lps_memsize(const struct lps *lps_ptr);

/**
 * @brief Deepens the compression level of the LCP structure. This method compresses the
 * existing cores and finds new cores.
 *
 * @param lps_ptr The `lps` object that will be parsed over.
 * @return 1 if successful in deepening the structure, 0 otherwise.
 */
int lps_deepen1(struct lps *lps_ptr);

/**
 * @brief Deepens the compression level of the LCP structure to a specific level.
 *
 * @param lps_ptr The `lps` object that will be parsed over.
 * @param lcp_level The target compression level to deepen to.
 * @return 1 if deepening was successful, 0 otherwise.
 */
int lps_deepen(struct lps *lps_ptr, int lcp_level);

/**
 * @brief Outputs the representation of a `lcp` pointer.
 *
 * This function iterates over the cores in the `cores` array and outputs them.
 * It prints each cores as bit representation, after initially printing the LCP level.
 *
 * @param lps_ptr A pointer to the `lps` object to be output.
 */
void print_lps(const struct lps *lps_ptr);

/**
 * @brief Equality operator for comparing two lcp::lps objects.
 *
 * This function compares the sizes of the core arrays of both objects, and then
 * compares each core element by dereferencing the core pointers. If all elements match,
 * the function returns true; otherwise, it returns false.
 *
 * @param lhs The left-hand side lps object to compare.
 * @param rhs The right-hand side lps object to compare.
 * @return 1 if both lps objects are equal, 0 otherwise.
 */
int lps_eq(const struct lps *lhs, const struct lps *rhs);

/**
 * @brief Inequality operator for comparing two lcp::lps objects.
 *
 * This function first checks if the sizes of the core arrays are different.
 * If they are, the objects are not equal. It then checks each core element.
 * If any element is different, the function returns true; otherwise, it returns false.
 *
 * @param lhs The left-hand side lps object to compare.
 * @param rhs The right-hand side lps object to compare.
 * @return 1 if the lps objects are not equal, 0 otherwise.
 */
int lps_neq(const struct lps *lhs, const struct lps *rhs);

/**
 * @brief Clears the lps object.
 *
 * @param lps_ptr The `lps` object to be cleared.
 */
void lps_clear(struct lps *lps_ptr);

#ifdef __cplusplus
}
#endif

#endif
/**
 * @file encoding.h
 * @brief Provides utilities for encoding and reverse complementing DNA
 * sequences using a character-to-bit mapping system.
 *
 * This file defines functions and variables for managing the encoding of DNA
 * sequences into bit-encoded values, as well as their reverse complement
 * encodings. It includes functionality for initializing encoding coefficients,
 * displaying the encoding summary, and loading custom encodings from
 * user-specified files or maps.
 *
 * Key functionalities include:
 *   - Initializes encoding coefficients for standard DNA bases (A, C, G, T) and
 * their reverse complements.
 *   - Supports custom encoding initialization via a map or file input.
 *   - Automatically calculates dictionary bit size based on the maximum
 * encoding values.
 *   - Provides reverse complement encodings for DNA sequences, a critical
 * feature in bioinformatics.
 *   - Displays an overview of the encoding scheme including coefficients and
 * bit size.
 *   - Loads encoding mappings from an external file, making it easy to extend
 * the encoding system for custom alphabets or symbols.
 *
 * Example usage:
 * @code
 *   // Initialize standard encoding and reverse complements
 *   lcp::init_coefficients(true);
 *
 *   // Initialize encoding from a file
 *   std::string encoding_file = "dna_encoding.txt";
 *   lcp::init_coefficients(encoding_file, true);
 * @endcode
 *
 * @see constant.h
 * @see core.h
 *
 * @namespace lcp
 *
 * @note Initialization functions throw std::invalid_argument exceptions for
 * invalid maps or file formats.
 *
 * @author Akmuhammet Ashyralyyev
 * @version 1.0
 * @date 2024-09-14
 *
 */

#ifndef ENCODING_H
#define ENCODING_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#define maximum(a, b) ((a) > (b) ? (a) : (b))

extern int alphabet[128];
extern int rc_alphabet[128];
extern char characters[128];

/**
 * @brief Displays the alphabet encoding summary including coefficients
 * and dictionary bit size.
 */
void LCP_SUMMARY(void);

/**
 * @brief Initializes the encoding coefficients for standard DNA bases
 * (A, C, G, T) and their reverse complements. Sets default values for
 * coefficients and dictionary bit size.
 */
void LCP_INIT(void);

/**
 * @brief Initializes the encoding coefficients for standard DNA bases
 * (A, C, G, T) and their reverse complements. Sets default values for
 * coefficients and dictionary bit size.
 * @param verbose If 0, prints the encoding summary after initialization.
 */
void LCP_INIT2(int verbose);

/**
 * @brief Initializes the encoding coefficients by reading them from a
 * file. The file must contain character, encoding, and reverse
 * complement encoding on each line.
 * @param filename Path to the file containing the character encodings.
 * @param verbose If true, prints the encoding summary after
 * initialization.
 * @return Always returns 0 upon successful initialization.
 * @throws std::invalid_argument if any invalid data is found in the
 * file.
 */
int LCP_INIT_FILE(const char *filename, int verbose);

#ifdef __cplusplus
}
#endif

#endif
/**
 * @file core.h
 * @brief Header file for the `core` struct, which represents bit-encoded
 * sequences and provides various utilities for encoding, compression, and file
 * I/O operations on bit sequences. (I/O operations will be supported later)
 *
 * This file contains the declaration of the `core` struct along with its member
 * functions, constructors, destructors, and operator overloads. It supports
 * operations for encoding strings into bit sequences, combining multiple
 * sequences, compressing data, and reading/writing from/to files.
 *
 * Key functionalities include:
 * - Encoding sequences of characters (e.g., ACGT) into a compact bit-encoded
 * form.
 * - Supporting reverse complement encoding of DNA sequences.
 * - Saving and loading core from files.
 * - Calculating memory usage of the constructed core structure.
 *
 * Dependencies:
 * - Requires constant.h and encoding.h for auxiliary data structures and
 * utilities.
 *
 * @see constant.h
 * @see encoding.h
 *
 * @namespace lcp
 * @struct core
 *
 *
 * @author Akmuhammet Ashyralyyev
 * @version 1.0
 * @date 2024-09-14
 *
 */

#ifndef CORE_H
#define CORE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> 

#define DCT_ITERATION_COUNT 1

#define minimum(a, b) ((a) < (b) ? (a) : (b))

typedef uint32_t ubit_size;
typedef uint32_t ulabel;

struct core {
    ubit_size bit_size;
    uint64_t bit_rep;
    ulabel label;
    uint64_t start;
    uint64_t end;
};


/**
 * @brief Computes the 32-bit MurmurHash3 hash for a given key.
 *
 * This function computes a 32-bit hash of the input data 'key' with the
 * specified length 'len' and an optional seed value. It processes the
 * input in blocks and handles any remaining bytes.
 *
 * @param key Pointer to the data to be hashed.
 * @param len The length of the data in bytes.
 * @param seed An initial seed value for the hash computation.
 * @return The resulting 32-bit hash value.
 */
uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed) ;

/**
 * @brief Initializes a core structure with the provided string data and index range.
 * 
 * This function processes a given substring starting at `begin` with a specified 
 * distance (length of the substring) and assigns start and end indices for tracking.
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param begin Pointer to the start of the string data.
 * @param distance Length of the substring to process.
 * @param start_index Start index of the substring within the data.
 * @param end_index End index of the substring within the data.
 */
void init_core1(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index);

/**
 * @brief Initializes a core structure with the provided string data and index range.
 * 
 * Similar to `init_core1`, but initalizes the core structure with reverse complement
 * alphabet encoding.
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param begin Pointer to the start of the string data.
 * @param distance Length of the substring to process.
 * @param start_index Start index of the substring within the data.
 * @param end_index End index of the substring within the data.
 */
void init_core2(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index);

/**
 * @brief Initializes a core structure by combining data from other core structures.
 * 
 * This function initializes a new core structure (`cr`) using a sequence of 
 * `core` objects starting from `begin` with the specified `distance` (number 
 * of `core` objects to process).
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param begin Pointer to the start of the sequence of core structures.
 * @param distance Number of core structures to process in the sequence.
 */
void init_core3(struct core *cr, struct core *begin, uint64_t distance);

/**
 * @brief Directly initializes a core structure with precomputed representation and metadata.
 * 
 * This function allows the initialization of a core structure when the bit 
 * representation, bit size, and label are already computed. Useful for 
 * deserializing or cloning a core structure.
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param bit_size Size of the bit representation in bits.
 * @param bit_rep Pointer to the precomputed bit representation array.
 * @param label Unique label assigned to the core structure.
 * @param start Start index of the substring or sequence represented by the core.
 * @param end End index of the substring or sequence represented by the core.
 */
void init_core4(struct core *cr, ubit_size bit_size, uint64_t bit_rep, ulabel label, uint64_t start, uint64_t end);

/**
 * @brief Compresses the right `core` object by comparing it with
 * left `core` object.
 *
 * This function compresses the core's bit sequence by identifying
 * common patterns between the right `core` object and left `core` object,
 * and reduces the size of the sequence accordingly.
 *
 * @param left_core The `core` object to compare against for compression.
 * @param right_core The `core` object that will be compressed.
 */
void core_compress(const struct core *left_core, struct core *right_core);

/**
 * @brief Output the bit representation of a `core` object.
 *
 * This function iterates over the bits in the `core` object and outputs
 * them to the given output stream. It prints each bit as either 0 or 1,
 * starting from the most significant bit to the least significant bit.
 *
 * @param cr The `core` object to be output.
 */
void print_core(const struct core *cr);

// core comparison operator implementation

/**
 * @brief Operator overload for equality comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the two objects are equal, 0 otherwise.
 */
int core_eq(const struct core *lhs, const struct core *rhs);

/**
 * @brief Operator overload for greater-than comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is greater, 0 otherwise.
 */
int core_neq(const struct core *lhs, const struct core *rhs);

/**
 * @brief Operator overload for smaller-than comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is smaller, 0 otherwise.
 */
int core_gt(const struct core *lhs, const struct core *rhs);

/**
 * @brief Operator overload for not-equal-to comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the two objects are not equal, 0 otherwise.
 */
int core_lt(const struct core *lhs, const struct core *rhs);

/**
 * @brief Operator overload for greater-than-or-equal-to comparison between
 * two `core` objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is greater than or equal, 0 otherwise.
 */
int core_geq(const struct core *lhs, const struct core *rhs);

/**
 * @brief Operator overload for smaller-than-or-equal-to comparison between
 * two `core` objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is smaller than or equal, 0 otherwise.
 */
int core_leq(const struct core *lhs, const struct core *rhs);

#ifdef __cplusplus
}
#endif

#endif
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
#ifdef LCPTOOLS_IMPL

void init_lps(struct lps *lps_ptr, const char *str, int len) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->capacity = (len/CONSTANT_FACTOR);
    lps_ptr->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1(str, str+len, lps_ptr->cores, 0);
}

void init_lps_offset(struct lps *lps_ptr, const char *str, int len, uint64_t offset) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->capacity = (len/CONSTANT_FACTOR);
    lps_ptr->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1(str, str+len, lps_ptr->cores, offset);
}

void init_lps2(struct lps *lps_ptr, const char *str, int len) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->capacity = (len/CONSTANT_FACTOR);
    lps_ptr->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse2(str, str+len, lps_ptr->cores, 0);
}

void init_lps3(struct lps *lps_ptr, FILE *in) {
    // read the level from the binary file
    if (fread(&(lps_ptr->level), sizeof(int), 1, in) != 1) {
        fprintf(stderr, "Error reading level from file\n");
        exit(EXIT_FAILURE);
    }

    // read the size (number of cores)
    if(fread(&(lps_ptr->size), sizeof(int), 1, in) != 1) {
        fprintf(stderr, "Error reading size from file\n");
        exit(EXIT_FAILURE);
    }

    lps_ptr->cores = NULL;

    if (lps_ptr->size) {
        // allocate memory for the cores array
        lps_ptr->capacity = (lps_ptr->size);
        lps_ptr->cores = (struct core *)malloc(lps_ptr->size * sizeof(struct core));
        if (fread(lps_ptr->cores, lps_ptr->size * sizeof(struct core), 1, in) != 1) {
            fprintf(stderr, "Error reading cores from file\n");
            exit(EXIT_FAILURE);
        }
    }
}

void init_lps4(struct lps *lps_ptr, const char *str, int len, int lcp_level, int chunk_size) {

    if (lcp_level < 1)
        return;

    lps_ptr->level = 1;
    lps_ptr->size = 0; 
    int estimated_size = (int)(len / pow((double)CONSTANT_FACTOR, lcp_level));
    lps_ptr->cores = (struct core *)malloc(estimated_size*sizeof(struct core));
    lps_ptr->capacity = (estimated_size);

    int str_index = 0, core_index = 0;

    {
        int str_len = minimum(chunk_size, len);
        struct lps temp_lps;
        init_lps_offset(&temp_lps, str, str_len, 0);
        lps_deepen(&temp_lps, lcp_level);

        if (temp_lps.size) {
            memcpy(lps_ptr->cores, temp_lps.cores, (temp_lps.size)*sizeof(struct core));
            core_index = (temp_lps.size);
            lps_ptr->size = (temp_lps.size);
            if (temp_lps.size>1)
                str_index = lps_ptr->cores[core_index-2].start;
            else 
                str_index = lps_ptr->cores[core_index-1].start;
        }
        free(temp_lps.cores);
    }

    while (str_index < len) {
        int str_len = minimum(chunk_size, len-str_index);
        struct lps temp_lps;
        init_lps_offset(&temp_lps, str+str_index, str_len, str_index);
        lps_deepen(&temp_lps, lcp_level);

        if (1<temp_lps.size) {
            int overlap = 2;
            while (0<overlap) {
                if (lps_ptr->cores[core_index-overlap].start == temp_lps.cores[0].start)
                    break;
                overlap--;
            }
            memcpy(lps_ptr->cores+core_index, temp_lps.cores+overlap, (temp_lps.size-overlap)*sizeof(struct core));
            core_index += (temp_lps.size-overlap);
            lps_ptr->size += (temp_lps.size-overlap);

            if ((uint64_t)str_index < lps_ptr->cores[core_index-2].start) {
                str_index = lps_ptr->cores[core_index-2].start;
                free(temp_lps.cores);
                continue;
            } 
        }
        
        // find next start point
        for(int i=str_index+str_len-1; str_index <= i; i--) {
            if (alphabet[(unsigned char)*(str+i)] == -1) {
                str_index = i+1;
                break;
            }
        }
        if (alphabet[(unsigned char)*(str+str_index)] != -1) { // all of the characters are valid, so not valid cores found
            str_index += str_len;
        }
        
        free(temp_lps.cores);
    }

    if (lps_ptr->size && lps_ptr->size >= lps_ptr->capacity){
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));
        lps_ptr->capacity = lps_ptr->size;
    }
}

void free_lps(struct lps *lps_ptr) {
    free(lps_ptr->cores);
    lps_ptr->size = 0;
}

void write_lps(struct lps *lps_ptr, FILE *out) {
    // write the level field
    fwrite(&(lps_ptr->level), sizeof(int), 1, out);

    // write the size (number of cores)
    fwrite(&(lps_ptr->size), sizeof(int), 1, out);

    // write each core object iteratively
    if (lps_ptr->size) {
        fwrite(lps_ptr->cores, lps_ptr->size*sizeof(struct core), 1, out);
    }
}

int parse1(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = begin;
    const char *it2 = end;
    int core_index = 0;
    int last_invalid_char_index = -1;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (alphabet[(unsigned char)*it1] == -1) {
            last_invalid_char_index = it1 - begin;
            continue;
        }

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
                if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                    init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                    core_index++;
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                init_core1(&(cores[core_index]), it1, 2+middle_count, it1-begin+offset, it2-begin+offset);
                core_index++;

                continue;
            }
        }

        if (alphabet[(unsigned char)*it1] > alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] < alphabet[(unsigned char)*(it1+2)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                core_index++;
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, 3, it1-begin+offset, it2-begin+offset);
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
            alphabet[(unsigned char)*(it1+2)] >= alphabet[(unsigned char)*(it1+3)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                core_index++;
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, 3, it1-begin+offset, it2-begin+offset);
            core_index++;

            continue;
        }
    }

    return core_index;
}

int parse2(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = end - 1;
    const char *it2 = begin - 1;
    int core_index = 0;

    // find lcp cores
    for (; begin <= it1 - 2; it1--) {

        // skip invalid character
        if (rc_alphabet[(unsigned char)*it1] == rc_alphabet[(unsigned char)*(it1-1)]) {
            continue;
        }

        // check for RINT core
        if (rc_alphabet[(unsigned char)*(it1-1)] == rc_alphabet[(unsigned char)*(it1-2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 - 2;
            while (begin <= temp && rc_alphabet[(unsigned char)*(temp+1)] == rc_alphabet[(unsigned char)*temp]) {
                temp--;
                middle_count++;
            }
            if (begin <= temp) {
                // check if there is any SSEQ cores left behind
                if (it1 < it2) {
                    init_core2(&(cores[core_index]), it2+1, it2-it1+2, end-it2-1+offset, end-it1-1+offset);
                    core_index++;
                }

                // create RINT core
                it2 = it1 - 2 - middle_count;
                init_core2(&(cores[core_index]), it1, 2+middle_count, end-it1-1+offset, end-it2-1+offset);
                core_index++;

                continue;
            }
        }

        if (rc_alphabet[(unsigned char)*it1] > rc_alphabet[(unsigned char)*(it1-1)] &&
            rc_alphabet[(unsigned char)*(it1-1)] < rc_alphabet[(unsigned char)*(it1-2)]) {

            // check if there is any SSEQ cores left behind
            if (it1 < it2) {
                init_core2(&(cores[core_index]), it2+1, it2-it1+2, end-it2-1+offset, end-it1-1+offset);
                core_index++;
            }

            // create LMIN core
            it2 = it1 - 3;
            init_core2(&(cores[core_index]), it1, 3, end-it1-1+offset, end-it2-1+offset);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (begin <= it1-3 &&
            rc_alphabet[(unsigned char)*it1] < rc_alphabet[(unsigned char)*(it1-1)] &&
            rc_alphabet[(unsigned char)*(it1-1)] > rc_alphabet[(unsigned char)*(it1-2)] &&
            rc_alphabet[(unsigned char)*(it1+1)] <= rc_alphabet[(unsigned char)*(it1)] &&
            rc_alphabet[(unsigned char)*(it1-2)] >= rc_alphabet[(unsigned char)*(it1-3)]) {

            // check if there is any SSEQ cores left behind
            if (it1 < it2) {
                init_core2(&(cores[core_index]), it2+1, it2-it1+2, end-it2-1+offset, end-it1-1+offset);
                core_index++;
            }

            // create LMAX core
            it2 = it1 - 3;
            init_core2(&(cores[core_index]), it1, 3, end-it1-1+offset, end-it2-1+offset);
            core_index++;

            continue;
        }
    }

    return core_index;
}

int parse3(struct core *begin, struct core *end, struct core *cores) {

    struct core *it1 = begin;
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
                    init_core3(&(cores[core_index]), it2-1, it1-it2+2);
                    core_index++;
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                init_core3(&(cores[core_index]), it1, it2-it1);
                core_index++;

                continue;
            }
        }

        // check for LMIN
        if (core_gt(it1, it1+1) && core_lt(it1+1, it1+2)) {
            
            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core3(&(cores[core_index]), it2-1, it1-it2+2);
                core_index++;
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1, it2-it1);
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
                init_core3(&(cores[core_index]), it2-1, it1-it2+2);
                core_index++;
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1, it2-it1);
            core_index++;

            continue;
        }
    }
    return core_index;
}

int64_t lps_memsize(const struct lps *lps_ptr) {
    return sizeof(struct lps) + lps_ptr->size * sizeof(struct core);
}

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
int lcp_dct(struct lps *lps_ptr) {

    // at least 2 cores are needed for compression
    if (lps_ptr->size < DCT_ITERATION_COUNT + 1) {
        return -1;
    }

    for (uint64_t dct_index = 0; dct_index < DCT_ITERATION_COUNT; dct_index++) {
        struct core *it_left = lps_ptr->cores + lps_ptr->size - 2, *it_right = lps_ptr->cores + lps_ptr->size - 1;

        for (; lps_ptr->cores + dct_index <= it_left; it_left--, it_right--) {
            core_compress(it_left, it_right);
        }
    }

    return 0;
}

int lps_deepen1(struct lps *lps_ptr) {

    // compress cores
    if (lcp_dct(lps_ptr) < 0) {
        lps_ptr->size = 0;
        lps_ptr->level++;
        return 0;
    }

    // find new cores
    int new_size = parse3(lps_ptr->cores + DCT_ITERATION_COUNT, lps_ptr->cores + lps_ptr->size, lps_ptr->cores);
    int temp = new_size;

    // remove old cores
    while(temp < lps_ptr->size) {
        temp++;
    }
    lps_ptr->size = new_size;

    lps_ptr->level++;

    if (lps_ptr->size && lps_ptr->size >= lps_ptr->capacity){
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));
        lps_ptr->capacity = lps_ptr->size;
    }

    return 1;
}

int lps_deepen(struct lps *lps_ptr, int lcp_level) {

    if (lcp_level <= lps_ptr->level)
        return 0;

    while (lps_ptr->level < lcp_level && lps_deepen1(lps_ptr))
        ;

    return 1;
}

void print_lps(const struct lps *lps_ptr) {
    printf("Level: %d \n", lps_ptr->level);
    for(int i=0; i<lps_ptr->size; i++) {
        print_core(&(lps_ptr->cores[i]));
        printf(" ");
    }
}

int lps_eq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 0;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 0;
        }
    }

    return 1;
}

int lps_neq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 1;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 1;
        }
    }

    return 0;
}

void lps_clear(struct lps *lps_ptr) {
    lps_ptr->size = 0;
}
/**
 * @file encoding.c
 * @brief Implementation of encoding functions.
 *
 * This file contains the implementation of encoding functions used to
 * initialize the alphabet with their corresponding coefficients. The encodings
 * support initialization with default coefficients, specific coefficients, or
 * by reading coefficients from a file.
 */


int alphabet[128];
int rc_alphabet[128];
char characters[128];

void LCP_SUMMARY(void) {
    printf("# Alphabet encoding summary\n");
    printf("# Coefficients: ");
    for (int i = 0; i < 128; i++) {
        if (alphabet[i] != -1) {
            printf("%c:%d ", i, alphabet[i]);
        }
    }
    printf("\n");
}

void LCP_INIT(void) {
    LCP_INIT2(0);
}

void LCP_INIT2(int verbose) {

    // init coefficients A/a=0, T/t=3, G/g=2, C/c=1
    for (int current_index = 0; current_index < 128; current_index++) {
        alphabet[current_index] = -1;
        characters[current_index] = 126;
    }
    alphabet['A'] = 0; alphabet['a'] = 0;
    alphabet['T'] = 3; alphabet['t'] = 3;
    alphabet['G'] = 2; alphabet['g'] = 2;
    alphabet['C'] = 1; alphabet['c'] = 1;

    rc_alphabet['A'] = 3; rc_alphabet['a'] = 3;
    rc_alphabet['T'] = 0; rc_alphabet['t'] = 0;
    rc_alphabet['G'] = 1; rc_alphabet['g'] = 1;
    rc_alphabet['C'] = 2; rc_alphabet['c'] = 2;

    characters[0] = 'A';
    characters[1] = 'C';
    characters[2] = 'G';
    characters[3] = 'T';

    if (verbose)
        LCP_SUMMARY();
}

int LCP_INIT_FILE(const char *encoding_file, int verbose) {
    
    FILE *encodings = fopen(encoding_file, "r");
    if (!encodings) {
        if (verbose) {
            fprintf(stderr, "Error: Could not open file %s\n", encoding_file);
        }
        return -1;
    }

    // clear arrays
    for (int current_index = 0; current_index < 128; current_index++) {
        alphabet[current_index] = -1;
        characters[current_index] = 126;
    }

    char character;
    int encoding, rev_encoding, mx = -1;
    while (fscanf(encodings, " %c %d %d", &character, &encoding, &rev_encoding) == 3) {
        alphabet[(unsigned char)character] = encoding;
        rc_alphabet[(unsigned char)character] = rev_encoding;

        mx = maximum(encoding, mx);
        mx = maximum(rev_encoding, mx);
    }

    fclose(encodings);

    int bit_count = 0;
    while (mx > 0) {
        bit_count++;
        mx = mx / 2;
    }

    if (bit_count != 2) {
        fprintf(stderr, "You alphabet has to have at most 2 binary digits in encoding. %d", bit_count);
        exit(EXIT_FAILURE);
    }

    return 0;
}
/**
 * @file core.c
 * @brief Implementation of the `core` struct and its associated functions.
 *
 * This file contains the implementation of the `core` struct, which is used to
 * represent a sequence of encoded bits for string data. The stuct supports
 * operations such as compression, comparison, and writing/reading to files.
 *
 * Key operations include:
 * - Encoding strings into bit arrays using coefficient-based encoding.
 * - Constructing `core` objects from strings or other `core` objects.
 * - Compressing bit representations to optimize memory usage.
 * - Writing and reading `core` objects to and from files.
 * - Comparing `core` objects with overloaded operators.
 * - Efficiently handling block-wise bit manipulations.
 *
 * @note The `STATS` macro is used to conditionally compile sections of the code
 * that track additional metadata such as `start` and `end` indices for
 * performance analysis.
 */

#include <inttypes.h>
/**
 * @brief Computes the 32-bit MurmurHash3 hash for a given key.
 *
 * This function computes a 32-bit hash of the input data 'key' with the
 * specified length 'len' and an optional seed value. It processes the
 * input in blocks and handles any remaining bytes.
 *
 * @param key Pointer to the data to be hashed.
 * @param len The length of the data in bytes.
 * @param seed An initial seed value for the hash computation.
 * @return The resulting 32-bit hash value.
 */
uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed) {
    const uint8_t *data = (const uint8_t *)key;
    const int nblocks = len / 4;

    uint32_t h1 = seed;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    // Body: Process blocks of 4 bytes at a time
    const uint32_t *blocks = (const uint32_t *)(data + nblocks * 4);

    for (int i = -nblocks; i; i++) {
        uint32_t k1 = blocks[i];

        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;

        h1 ^= k1;
        h1 = (h1 << 15) | (h1 >> (32 - 15));
        h1 = h1 * 5 + 0xe6546b64;
    }

    // Tail: Process remaining bytes
    const uint8_t *tail = (const uint8_t *)(data + nblocks * 4);

    uint32_t k1 = 0;

    switch (len & 3) {
    case 3:
        k1 ^= tail[2] << 16;
        break;
    case 2:
        k1 ^= tail[1] << 8;
        break;
    case 1:
        k1 ^= tail[0];
        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;
        h1 ^= k1;
    }

    // Finalization: Mix the hash to ensure the last few bits are fully mixed
    h1 ^= len;

    /* fmix32 */
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    return h1;
}

void init_core1(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index) {
    cr->start = start_index;
    cr->end = end_index;
    cr->label = 0;
    cr->label |= ((distance-2) << 6);
    cr->label |= (alphabet[(int)(*begin)] << 4);
    cr->label |= (alphabet[(int)(*(begin+distance-2))] << 2);
    cr->label |= (alphabet[(int)(*(begin+distance-1))]);
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * distance;
}

void init_core2(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index) {
    cr->start = start_index;
    cr->end = end_index;
    cr->label = 0;
    cr->label |= ((distance-2) << 6);
    cr->label |= (rc_alphabet[(int)(*(begin))] << 4);
    cr->label |= (rc_alphabet[(int)(*(begin-distance+2))] << 2);
    cr->label |= (rc_alphabet[(int)(*(begin-distance+1))]);
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * distance;
}

void init_core3(struct core *cr, struct core *begin, uint64_t distance) {
    cr->start = begin->start;
    cr->end = (begin+distance-1)->end;
    cr->bit_rep = 0;
    cr->bit_size = 0;

    for (struct core *it=begin; it<begin+distance; it++) {
        cr->bit_size += it->bit_size;
    }

    int index = 0;
    for (struct core *it = begin+distance-1; begin <= it && index + it->bit_size <= 64; it--) {
        cr->bit_rep |= (it->bit_rep << index);
        index += it->bit_size;
    }

    cr->bit_rep = 0x7FFFFFFFFFFFFFFF & cr->bit_rep;
    cr->bit_size = minimum(cr->bit_size, 63);

    ulabel data[4];
    data[0] = (begin)->label;
    data[1] = (begin+distance-2)->label;
    data[2] = (begin+distance-1)->label;
    data[3] = distance-2;
    cr->label = MurmurHash3_32((void*)data, 4 * sizeof(ulabel), 42);
}

void init_core4(struct core *cr, ubit_size bit_size, uint64_t bit_rep, ulabel label, uint64_t start, uint64_t end) {
    cr->bit_size = bit_size;
    cr->bit_rep = bit_rep;
    cr->label = label;
    cr->start = start;
    cr->end = end;
}

void core_compress(const struct core *left_core, struct core *right_core) {
#ifdef LCP_LEVEL1_BRANCH
    if (left_core->bit_rep & 0x8000000000000000) { // if compressing 1-level cores
        uint64_t left_core_3 = (left_core->bit_rep) & 3;
        uint64_t left_core_2 = (left_core->bit_rep >> 2) & 3;
        uint64_t left_core_middle_count = (left_core->bit_rep & 0x7FFFFFFFFFFFFFFF) >> 6;
        uint64_t left_core_1 = (left_core->bit_rep >> 4) & 3;

        uint64_t right_core_3 = (right_core->bit_rep) & 3;
        uint64_t right_core_2 = (right_core->bit_rep >> 2) & 3;
        uint64_t right_core_middle_count = (right_core->bit_rep & 0x7FFFFFFFFFFFFFFF) >> 6;
        uint64_t right_core_1 = (right_core->bit_rep >> 4) & 3;

        if (left_core_3 != right_core_3) { // if right characters mismatches
            if ((left_core_3 & 1) != (right_core_3 & 1)) {
                right_core->bit_rep = (right_core_3 & 1); // 0b00 + r % 2
            } else {
                right_core->bit_rep = 2 + ((right_core_3 >> 1) & 1); // 0b10 + r % 2
            }
            right_core->bit_size = 2;
        }
        else if (left_core_2 != right_core_2) { // if middle characters mismatches
            if ((left_core_2 & 1) != (right_core_2 & 1)) {
                right_core->bit_rep = 4 + (right_core_2 & 1); // 0b100 + r % 2
            } else {
                right_core->bit_rep = 6 + ((right_core_2 >> 1) & 1); // 0b110 + r % 2
            }
            right_core->bit_size = (64 - __builtin_clzll(right_core->bit_rep));
        }
        else if (left_core_middle_count != right_core_middle_count) { // middle character counts mismatches
            if (left_core_middle_count < right_core_middle_count) {
                // compare left_core_1 with right_core_2
                if ((left_core_1 & 1) != (right_core_2 & 1)) {
                    right_core->bit_rep = 4 * (left_core_middle_count + 1) + (right_core_2 & 1); // 2 * 2 * (mid + 1) + r % 2
                } else {
                    right_core->bit_rep = 2 * (2 * (left_core_middle_count + 1) + 1) + ((right_core_2 >> 1) & 1); // 2 * (2 * (mid + 1) + 1) + r % 2
                }
                right_core->bit_size = (64 - __builtin_clzll(right_core->bit_rep));
            } else {
                // compare left_core_2 with right_core_1
                if ((left_core_2 & 1) != (right_core_1 & 1)) {
                    right_core->bit_rep = 4 * (right_core_middle_count + 1) + (right_core_1 & 1);
                } else {
                    right_core->bit_rep = 2 * (2 * (right_core_middle_count + 1) + 1) + ((right_core_1 >> 1) & 1);
                }
                right_core->bit_size = (64 - __builtin_clzll(right_core->bit_rep));
            }
        }
        else if (left_core_1 != right_core_1) { // left characters mismatches
            if ((left_core_1 & 1) != (right_core_1 & 1)) {
                right_core->bit_rep = 4 * (left_core_middle_count + 1) + (right_core_1 & 1);
            } else {
                right_core->bit_rep = 2 * (2 * (left_core_middle_count + 1) + 1) + ((right_core_1 >> 1) & 1);
            }
            right_core->bit_size = (64 - __builtin_clzll(right_core->bit_rep));
        }
        else { // they are same
            right_core->bit_rep = 2 * right_core->bit_size;
            right_core->bit_size = (64 - __builtin_clzll(right_core->bit_rep));
        }
    } else
#endif /* LCP_LEVEL1_BRANCH */
    { // generic path: first differing bit from the LSB, bounded by bit_size
        ubit_size first_differing_index = 64;
        if (left_core->bit_rep != right_core->bit_rep) {
            first_differing_index = __builtin_ctzll(left_core->bit_rep ^ right_core->bit_rep); // trailing zero count (0-index)
        } else {
            first_differing_index = right_core->bit_size;
        }
        first_differing_index = minimum(first_differing_index, minimum(left_core->bit_size, right_core->bit_size));
        right_core->bit_rep = 2 * first_differing_index + ((right_core->bit_rep >> first_differing_index) & 1);
        right_core->bit_size = right_core->bit_rep == 0 ? 2 : (64 - __builtin_clzll(right_core->bit_rep));
        right_core->bit_size = right_core->bit_size < 2 ? 2 : right_core->bit_size;
    }

    // now, the right core is dependent on the left; hence, its coverage spans towards the left
    right_core->start = left_core->start;
}

void print_core(const struct core *cr) {
    if (cr->bit_rep & 0x8000000000000000) { // if printing 1-level cores
        uint64_t middle_count = (0x7FFFFFFFFFFFFFFF & cr->bit_rep) >> 6;
        uint64_t middle_val = (cr->bit_rep >> 2) & 3;
        printf("%" PRIu64, ((cr->bit_rep >> 5) & 1));
        printf("%" PRIu64, ((cr->bit_rep >> 4) & 1));
        for (uint64_t i=0; i<middle_count; i++) {
            printf("%" PRIu64, ((middle_val >> 1) & 1));
            printf("%" PRIu64, (middle_val & 1));
        }
        printf("%" PRIu64, ((cr->bit_rep >> 1) & 1));
        printf("%" PRIu64, (cr->bit_rep & 1));
    } else {
        for (ubit_size index = cr->bit_size - 1; 0 < index; index--) {
            printf("%" PRIu64, ((cr->bit_rep >> index) & 1));
        }
        printf("%" PRIu64, (cr->bit_rep & 1));
    }
}

// core comparison operator implementation

int core_eq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep == rhs->bit_rep;
}

int core_neq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep != rhs->bit_rep;
}

int core_gt(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep > rhs->bit_rep;
}

int core_lt(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep < rhs->bit_rep;
}

int core_geq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep >= rhs->bit_rep;
}

int core_leq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep <= rhs->bit_rep;
}

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
#endif
#endif
