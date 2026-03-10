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

#include <pthread.h>
#include <stdio.h>
#include <math.h>

#ifndef CONSTANT_FACTOR
#define CONSTANT_FACTOR     1.5
#endif

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
 * @brief Deepens the compression level of the LCP structure in parallel. This method
 * compresses the existing cores and finds new cores.
 *
 * @param lps_ptr The `lps` object that will be parsed over.
 * @param thread_number The number of theads to run the DCT.
 * @return 1 if successful in deepening the structure, 0 otherwise.
 */
int lps_deepen1_parallel(struct lps *lps_ptr, int thread_number);

/**
 * @brief Deepens the compression level of the LCP structure to a specific level in threads.
 *
 * @param lps_ptr The `lps` object that will be parsed over.
 * @param lcp_level The target compression level to deepen to.
 * @param thread_number The number of theads to run the DCT.
 * @return 1 if deepening was successful, 0 otherwise.
 */
int lps_deepen_parallel(struct lps *lps_ptr, int lcp_level, int thread_number);

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

typedef struct {
    uint64_t x;
    uint64_t y;
} uint128_t;
#if defined(__cplusplus)
    #define CL(type)      type
#else
    #define CL(type)      (type)
#endif


struct core {
    ubit_size bit_size;
    uint128_t bit_rep;
    ulabel label;
    uint64_t start;
    uint64_t end;
};

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
static inline int core_eq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep.x == rhs->bit_rep.x && lhs->bit_rep.y == rhs->bit_rep.y;
}

/**
 * @brief Operator overload for greater-than comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is greater, 0 otherwise.
 */
static inline int core_neq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep.x != rhs->bit_rep.x || lhs->bit_rep.y != rhs->bit_rep.y;
}

/**
 * @brief Operator overload for smaller-than comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is smaller, 0 otherwise.
 */
static inline int core_gt(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep.x > rhs->bit_rep.x || (lhs->bit_rep.x == rhs->bit_rep.x && lhs->bit_rep.y > rhs->bit_rep.y);
}

/**
 * @brief Operator overload for not-equal-to comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the two objects are not equal, 0 otherwise.
 */
static inline int core_lt(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep.x < rhs->bit_rep.x || (lhs->bit_rep.x == rhs->bit_rep.x && lhs->bit_rep.y < rhs->bit_rep.y);
}

/**
 * @brief Operator overload for greater-than-or-equal-to comparison between
 * two `core` objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is greater than or equal, 0 otherwise.
 */
static inline int core_geq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep.x > rhs->bit_rep.x || (lhs->bit_rep.x == rhs->bit_rep.x && lhs->bit_rep.y >= rhs->bit_rep.y);
}

/**
 * @brief Operator overload for smaller-than-or-equal-to comparison between
 * two `core` objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is smaller than or equal, 0 otherwise.
 */
static inline int core_leq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep.x < rhs->bit_rep.x || (lhs->bit_rep.x == rhs->bit_rep.x && lhs->bit_rep.y <= rhs->bit_rep.y);
}

/**
 * @brief Return the minimum of two ubit_size values.
 *
 * @param a First value.
 * @param b Second value.
 * @return The smaller of a and b.
 */
static inline ubit_size umin(ubit_size a, ubit_size b) {
    return a < b ? a : b;
}

/**
 * @brief Compute the bit-length of a 128-bit unsigned integer.
 *
 * Bit-length is defined as the position of the most significant set bit
 * (1-based). Returns 0 if x == 0.
 *
 * Uses __builtin_clzll for efficient leading-zero counting.
 *
 * @param x Input 128-bit unsigned integer.
 * @return Number of significant bits required to represent x.
 */
static inline ubit_size bitlen_u128(uint128_t x) {
    return x.x ?
    ((ubit_size)(64u - (ubit_size)__builtin_clzll(x.x)) + 64) :
    ( x.y ? (ubit_size)(64u - (ubit_size)__builtin_clzll(x.y)) : 0u);
}

/**
 * @brief Compute the bit-length of x with a minimum result of 2.
 *
 * Ensures the returned bit-length is at least 2, even if x has fewer
 * significant bits (or is zero).
 *
 * @param x Input 128-bit unsigned integer.
 * @return max(bitlen_u128(x), 2).
 */
static inline ubit_size bitlen_min2_128(uint128_t x) {
    ubit_size bl = bitlen_u128(x);
    return bl < 2u ? 2u : bl;
}

/**
 * @brief Extract a 2-bit symbol from the low 6 bits of a packed value.
 *
 * The low 6 bits are interpreted as three 2-bit slots:
 *   k = 0 → rightmost 2 bits
 *   k = 1 → middle 2 bits
 *   k = 2 → leftmost 2 bits
 *
 * @param rep Packed representation.
 * @param k   Slot index (0–2).
 * @return The 2-bit symbol at slot k.
 */
static inline uint128_t sym2_128(uint128_t rep, unsigned k) {
    return k ? CL(uint128_t){0, (rep.y >> k) & 3ull} : CL(uint128_t){0, rep.y & 3ull};
}

/**
 * @brief Extract the middle repetition count stored above the low 6 bits.
 *
 * Bits 0–5 are reserved for 2-bit symbols. Bits above 6 contain the
 * repetition count. The most significant bit is masked off before shifting.
 *
 * @param rep Packed representation.
 * @return Middle repetition count.
 */
static inline uint64_t mid_count_128(uint128_t rep) {
    return rep.y >> 6;
}

/**
 * @brief Emit an encoded index-bit value based on symbol comparison.
 *
 * Computes:
 *   result = i + selected_bit or 2 + i + selected_bit, depending on
 *  from where the bit is selected
 *
 * If the least significant bits of a2 and b2 differ, select (b2 & 1).
 * Otherwise, select ((b2 >> 1) & 1).
 *
 * @param a2 First 2-bit symbol.
 * @param b2 Second 2-bit symbol.
 * @param i  Base index multiplier.
 * @return Encoded value 2 + i plus chosen bit from b2.
 */
static inline uint128_t emit_idx_bit_128(uint128_t a2, uint128_t b2, uint64_t i) {
    if ((a2.y & 1) != (b2.y & 1)) {
        return CL(uint128_t){0, i + (b2.y & 1)};
    }
    return CL(uint128_t){0, 2 + i + ((b2.y >> 1) & 1)};
}

/**
 * @brief Perform level-1 compression of two adjacent cores.
 *
 * Implements the 3-symbol + middle-count encoding scheme.
 * The caller guarantees that both `left` and `right` are already
 * level-1 encoded.
 *
 * The function compares corresponding components of the packed
 * representations in the following priority order:
 *
 *   1. Rightmost 2-bit symbol (L3 vs R3)
 *   2. Leftmost 2-bit symbol (L2 vs R2)
 *   3. Middle repetition count (Lm vs Rm)
 *      - If counts differ, performs a boundary comparison:
 *          • If Lm < Rm: compare L1 with R2
 *          • Otherwise:  compare L2 with R1
 *        The emitted index depends on the smaller count.
 *   4. Final 2-bit symbol (L1 vs R1)
 *   5. If all components match, emit a “same” encoding.
 *
 * In each mismatch case, the emitted value is:
 *
 *     out = 2 * index + selected_bit
 *
 * where the selected bit is derived from the compared 2-bit symbols.
 *
 * Side effects:
 *   - Overwrites `right->bit_rep` with the compressed result.
 *   - Updates `right->bit_size` to bitlen_min2(out).
 *   - Propagates `left->start` into `right->start`.
 *
 * @param left   Pointer to left core (read-only, level-1 encoded).
 * @param right  Pointer to right core (level-1 encoded, updated in place).
 */
static inline void core_compress_level1(const struct core *left, struct core *right) {
    uint128_t L = left->bit_rep;
    uint128_t R = right->bit_rep;

    uint128_t L3 = sym2_128(L, 0), L2 = sym2_128(L, 2), L1 = sym2_128(L, 4);
    uint128_t R3 = sym2_128(R, 0), R2 = sym2_128(R, 2), R1 = sym2_128(R, 4);

    uint64_t Lm = mid_count_128(L);
    uint64_t Rm = mid_count_128(R);

    uint128_t out;

    if (L3.y != R3.y) {
        // base = 0 (2*i with i=0)
        out = emit_idx_bit_128(L3, R3, 0);
        right->bit_rep = out;
        right->bit_size = 2;
    }
    else if (L2.y != R2.y) {
        // Your original used base=4/6 -> i=2 (since 2*i = 4)
        out = emit_idx_bit_128(L2, R2, 4);
        right->bit_rep = out;
        right->bit_size = bitlen_min2_128(out);
    }
    else if (Lm != Rm) {
        // Preserve your “compare across boundary depending on which count is smaller”
        if (Lm < Rm) {
            // compare left L1 with right R2; index = Lm + 1
            out = emit_idx_bit_128(L1, R2, 4 * Lm + 4);
        } else {
            // compare left L2 with right R1; index = Rm + 1
            out = emit_idx_bit_128(L2, R1, 4 * Rm + 4);
        }
        right->bit_rep = out;
        right->bit_size = bitlen_min2_128(out);
    }
    else if (L1.y != R1.y) {
        // left mismatch: index = Lm + 1
        out = emit_idx_bit_128(L1, R1, 4 * Lm + 4);
        right->bit_rep = out;
        right->bit_size = bitlen_min2_128(out);
    }
    else {
        // same
        out = CL(uint128_t){0, 2ull * (uint64_t)right->bit_size};
        right->bit_rep = out;
        right->bit_size = bitlen_min2_128(out);
    }

    right->start = left->start;
}


/**
 * @brief Perform upper-level compression via bounded first-difference encoding.
 *
 * The caller guarantees that both `left` and `right` are already
 * upper-level encoded.
 *
 * The algorithm:
 *
 *   1. Compute a comparison bound:
 *        bound = min(left->bit_size, right->bit_size, 64).
 *
 *   2. If the bit representations are identical:
 *        idx = bound.
 *      Otherwise:
 *        - Compute x = left->bit_rep ^ right->bit_rep.
 *        - Find the least significant differing bit:
 *              idx = ctz(x).
 *        - Clamp idx to `bound`.
 *
 *   3. Emit:
 *        out = 2 * idx + bit_at_idx(right)
 *
 *      where bit_at_idx(right) is the bit of right->bit_rep at position idx.
 *
 * Side effects:
 *   - Overwrites `right->bit_rep` with the compressed result.
 *   - Updates `right->bit_size` to bitlen_min2(out).
 *   - Propagates `left->start` into `right->start`.
 *
 * This encoding effectively captures the first differing bit
 * (from the least significant side), bounded by the smaller size.
 *
 * @param left   Pointer to left core (read-only, upper-level encoded).
 * @param right  Pointer to right core (upper-level encoded, updated in place).
 */
static inline void core_compress_upper(const struct core *left, struct core *right) {
    ubit_size bound = umin(right->bit_size, umin(left->bit_size, 64u));

    ubit_size idx;
    if (left->bit_rep.x == right->bit_rep.x && left->bit_rep.y == right->bit_rep.y) {
        idx = bound;
    } else {
        if (left->bit_rep.y == right->bit_rep.y) {
            uint64_t x = left->bit_rep.x ^ right->bit_rep.x;   // nonzero here
            idx = (ubit_size)__builtin_ctzll(x) + 63;
            idx = umin(idx, bound);
        } else {
            uint64_t x = left->bit_rep.y ^ right->bit_rep.y;   // nonzero here
            idx = (ubit_size)__builtin_ctzll(x);
            idx = umin(idx, bound);
        }
    }

    uint128_t out = CL(uint128_t){0, 2ull * (uint64_t)idx + (( idx > 63 ? right->bit_rep.x >> (idx - 63) : right->bit_rep.y >> idx) & 1ull)};
    right->bit_rep = out;
    right->bit_size = bitlen_min2_128(out);
    right->start = left->start;
}

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

    if (lps_ptr->level == 1) {
        for (uint64_t dct_index = 0; dct_index < DCT_ITERATION_COUNT; dct_index++) {
            struct core *it_left = lps_ptr->cores + lps_ptr->size - 2, *it_right = lps_ptr->cores + lps_ptr->size - 1;

            for (; lps_ptr->cores + dct_index <= it_left; it_left--, it_right--) {
                core_compress_level1(it_left, it_right);
            }
        }
    } else {
        for (uint64_t dct_index = 0; dct_index < DCT_ITERATION_COUNT; dct_index++) {
            struct core *it_left = lps_ptr->cores + lps_ptr->size - 2, *it_right = lps_ptr->cores + lps_ptr->size - 1;

            for (; lps_ptr->cores + dct_index <= it_left; it_left--, it_right--) {
                core_compress_upper(it_left, it_right);
            }
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

typedef struct {
    struct core *cores;     // destination cores (real array)
    struct core dummy_left; // uncompressed core at begin-1
    int offset_begin;       // inclusive i
    int offset_end;         // exclusive i
    int flags;              // 1 bit (process head) 1 bit (level 1 or not)
} dct_worker_args_t;

/**
 * @brief Worker routine for parallel DCT compression.
 *
 * This function is executed by a single POSIX thread and performs
 * compression on a contiguous subrange of cores. The behavior depends
 * on the `flags` field:
 *
 * - flags == -1 : No work is performed (thread exits immediately).
 * - flags % 2   : Use core_compress_level1().
 * - otherwise   : Use core_compress_upper().
 *
 * If flags > 1, the worker also performs an additional boundary
 * compression using `dummy_left` and the first core in its range.
 * This is required when the thread's assigned range does not begin
 * at the global DCT starting index, ensuring correctness across
 * chunk boundaries.
 *
 * The function assumes that:
 * - Each thread operates on a disjoint range of cores.
 * - Any required boundary state is provided through `dummy_left`.
 *
 * @param argp Pointer to a dct_worker_args_t structure containing
 *             the thread parameters.
 *
 * @return Always returns NULL (required by pthread signature).
 */
static void *dct_worker(void *argp) {
    dct_worker_args_t *dct_params = (dct_worker_args_t *)argp;

    if (dct_params->flags == -1) return NULL;

    if (dct_params->flags % 2) {
        for (int i = dct_params->offset_end - 1; dct_params->offset_begin < i; i--) {
            core_compress_level1(dct_params->cores + (i - 1), dct_params->cores + i);
        }
        if (dct_params->flags > 1) {
            core_compress_level1(&(dct_params->dummy_left), dct_params->cores + dct_params->offset_begin);
        }
    } else {
        for (int i = dct_params->offset_end - 1; dct_params->offset_begin < i; i--) {
            core_compress_upper(dct_params->cores + (i - 1), dct_params->cores + i);
        }
        if (dct_params->flags > 1) {
            core_compress_upper(&(dct_params->dummy_left), dct_params->cores + dct_params->offset_begin);
        }
    }

    return NULL;
}

/**
 * @brief Executes one parallel DCT sweep over the cores array.
 *
 * This function divides the compression work starting at `dct_index`
 * into approximately equal chunks and distributes them among
 * `thread_number` threads. Each thread processes a contiguous
 * subrange of cores using the dct_worker() routine.
 *
 * Work partitioning:
 * - total = lps_ptr->size - dct_index
 * - chunk size is computed using ceiling division.
 * - Each thread receives a non-overlapping [begin, end) interval.
 *
 * Boundary handling:
 * - Threads whose range does not start at dct_index receive a
 *   copy of the left boundary core (`dummy_left`) to ensure
 *   correct cross-boundary compression.
 *
 * Memory management:
 * - Dynamically allocates thread handles and argument arrays.
 * - Joins all threads before returning.
 *
 * @param lps_ptr       Pointer to the LPS structure containing cores.
 * @param dct_index     Starting index for this DCT iteration.
 * @param thread_number Number of worker threads to spawn.
 *
 * @return 0 on success.
 * @return -1 on allocation failure.
 */
static int run_parallel_sweep(struct lps *lps_ptr, int dct_index, int thread_number) {

    const int total = (int)lps_ptr->size - dct_index;   // number of pairs (i, i+1)
    if (total <= 0) return 0;                               // nothing to do
    const int chunk = (total + thread_number - 1) / thread_number;

    pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * thread_number);
    dct_worker_args_t *args = (dct_worker_args_t *)malloc(sizeof(dct_worker_args_t) * thread_number);

    if (!threads || !args) {
        perror("LCP: couldn't allocate thread arguments\n");
        return -1;
    }

    for (int t = 0; t < thread_number; t++) {
        int begin = dct_index + t * chunk;
        int end = ((begin + chunk) <= lps_ptr->size ? (begin + chunk) : lps_ptr->size);

        args[t].cores = lps_ptr->cores;
        args[t].offset_begin = begin;
        args[t].offset_end = end;

        if ((int)lps_ptr->size <= begin || end <=begin) {
            args[t].flags = -1;
        } else if (begin != dct_index) {
            args[t].dummy_left = lps_ptr->cores[begin-1];
            args[t].flags = 2 + (lps_ptr->level == 1 ? 1 : 0);
        } else {
            args[t].flags = (lps_ptr->level == 1 ? 1 : 0);
        }

        pthread_create(&threads[t], NULL, dct_worker, &args[t]);
    }

    for (int t = 0; t < thread_number; t++) {
        pthread_join(threads[t], NULL);
    }

    free(threads);
    free(args);
    return 0;
}

/**
 * @brief Performs multi-iteration parallel Deterministic Coin Tossing (DCT).
 *
 * This function executes DCT_ITERATION_COUNT consecutive compression
 * sweeps over the cores stored in `lps_ptr`, using parallel execution
 * for each sweep.
 *
 * Each iteration reduces redundancy by compressing adjacent core pairs.
 * The process prepares the core sequence for subsequent parsing stages
 * in the LCP framework.
 *
 * Preconditions:
 * - At least DCT_ITERATION_COUNT + 1 cores must be available.
 *
 * @param lps_ptr       Pointer to the LPS structure.
 * @param thread_number Number of worker threads per sweep.
 *
 * @return 0 on success.
 * @return -1 if there are not enough cores for DCT.
 * @return Propagates non-zero errors from run_parallel_sweep().
 */
int lcp_dct_parallel(struct lps *lps_ptr, int thread_number) {
    // at least 2 cores are needed for compression
    if (lps_ptr->size < DCT_ITERATION_COUNT + 1) {
        return -1;
    }

    for (int dct_index = 0; dct_index < DCT_ITERATION_COUNT; dct_index++) {
        int rc = run_parallel_sweep(lps_ptr, dct_index, thread_number);
        if (rc != 0) {
            return rc;
        }
    }

    return 0;
}

/**
 * @brief Performs a single LPS deepening step using parallel DCT.
 *
 * This function advances the LPS structure by one level:
 *
 * 1. Executes parallel DCT compression.
 * 2. Parses the resulting cores using parse3() to identify new cores.
 * 3. Updates the core array and shrinks it via realloc().
 * 4. Increments the LPS level.
 *
 * If DCT cannot be performed (insufficient cores),
 * the structure is reset (size = 0) and the level is still incremented.
 *
 * @param lps_ptr       Pointer to the LPS structure to update.
 * @param thread_number Number of worker threads used during DCT.
 *
 * @return 1 if deepening produced a valid new level.
 * @return 0 if no further deepening was possible.
 */
int lps_deepen1_parallel(struct lps *lps_ptr, int thread_number) {

    // compress cores
    if (lcp_dct_parallel(lps_ptr, thread_number) < 0) {
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

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));

    return 1;
}

/**
 * @brief Deepens the LPS structure up to a target LCP level.
 *
 * Repeatedly invokes lps_deepen1_parallel() until either:
 * - The desired `lcp_level` is reached, or
 * - Further deepening is no longer possible.
 *
 * If the current level already satisfies the requested level,
 * no work is performed.
 *
 * @param lps_ptr       Pointer to the LPS structure.
 * @param lcp_level     Target LCP level to reach.
 * @param thread_number Number of worker threads used per deepening step.
 *
 * @return 1 if processing completed (even if stopped early).
 * @return 0 if no deepening was required.
 */
int lps_deepen_parallel(struct lps *lps_ptr, int lcp_level, int thread_number) {

    if (lcp_level <= lps_ptr->level)
        return 0;

    while (lps_ptr->level < lcp_level && lps_deepen1_parallel(lps_ptr, thread_number))
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
        h1 = (h1 << 15) | (h1 >> (32 - 15)); // it should be (h1 << 13) | (h1 >> (32 - 13)) but left at it is, let it be a legacy :)
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
    cr->bit_rep = CL(uint128_t){0x8000000000000000, cr->label};
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
    cr->bit_rep = CL(uint128_t){0x8000000000000000, cr->label};
    cr->bit_size = 2 * distance;
}

void init_core3(struct core *cr, struct core *begin, uint64_t distance) {
    cr->start = begin->start;
    cr->end = (begin+distance-1)->end;
    cr->bit_rep = CL(uint128_t){0, 0};
    cr->bit_size = 0;

    for (struct core *it=begin; it<begin+distance; it++) {
        cr->bit_size += it->bit_size;
    }

    int index = 0;
    for (struct core *it = begin+distance-1; begin <= it && index + it->bit_size <= 64; it--) {
        cr->bit_rep.y |= (it->bit_rep.y << index);
        index += it->bit_size;
    }

    cr->bit_rep.x = 0x7FFFFFFFFFFFFFFF & cr->bit_rep.x;
    cr->bit_size = minimum(cr->bit_size, 127);

    ulabel data[4];
    data[0] = (begin)->label;
    data[1] = (begin+distance-2)->label;
    data[2] = (begin+distance-1)->label;
    data[3] = distance-2;
    cr->label = MurmurHash3_32((void*)data, 4 * sizeof(ulabel), 42);
}

void init_core4(struct core *cr, ubit_size bit_size, uint64_t bit_rep, ulabel label, uint64_t start, uint64_t end) {
    cr->bit_size = bit_size;
    cr->bit_rep = CL(uint128_t){0, bit_rep};
    cr->label = label;
    cr->start = start;
    cr->end = end;
}

void print_core(const struct core *cr) {
    if (cr->bit_rep.x & 0x8000000000000000) { // if printing 1-level cores
        uint64_t middle_count = (cr->bit_rep.y) >> 6;
        uint64_t middle_val = (cr->bit_rep.y >> 2) & 3;
        printf("%" PRIu64, (uint64_t)((cr->bit_rep.y >> 5) & 1));
        printf("%" PRIu64, (uint64_t)((cr->bit_rep.y >> 4) & 1));
        for (uint64_t i=0; i<middle_count; i++) {
            printf("%" PRIu64, (uint64_t)((middle_val >> 1) & 1));
            printf("%" PRIu64, (uint64_t)(middle_val & 1));
        }
        printf("%" PRIu64, (uint64_t)((cr->bit_rep.y >> 1) & 1));
        printf("%" PRIu64, (uint64_t)(cr->bit_rep.y & 1));
    } else {
        for (ubit_size index = cr->bit_size - 1; 63 < index; index--) {
            printf("%" PRIu64, (uint64_t)((cr->bit_rep.x >> index) & 1));
        }
        for (ubit_size index = cr->bit_size - 1; 0 < index; index--) {
            printf("%" PRIu64, (uint64_t)((cr->bit_rep.y >> index) & 1));
        }
        printf("%" PRIu64, (uint64_t)(cr->bit_rep.y & 1));
    }
}
#endif
#endif
