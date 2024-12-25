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

#include "core.h"
#include "encoding.h"

#define CONSTANT_FACTOR         1.5

struct lps {
    int level;
    int size;
    struct core *cores;
};

void reverse(const char *str, int len, char **rev);

/**
 * @brief Constructs an lps object from a string, with an option to apply reverse complement
 * transformation.
 * 
 * @param nlps The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 * @param rev_comp Whether to apply reverse complement (default is false).
 */
void init_lps(struct lps *nlps, const char *str, int len, int rev_comp);

/**
 * @brief Destructor for the lps object. Frees dynamically allocated memory for cores.
 */
void free_lps(struct lps *nlps);

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
 * @return Size of the cores identified in the given string.
 */
int parse1(const char *begin, const char *end, struct core *cores);

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
 * @return Size of the cores identified in the given string.
 */
int parse2(const char *begin, const char *end, struct core *cores);

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
int64_t lps_memsize(const struct lps *str);

/**
 * @brief Deepens the compression level of the LCP structure. This method compresses the
 * existing cores and finds new cores.
 *
 * @param str The `lps` object that will be parsed over.
 * @return 1 if successful in deepening the structure, 0 otherwise.
 */
int lps_deepen1(struct lps *str);

/**
 * @brief Deepens the compression level of the LCP structure to a specific level.
 *
 * @param str The `lps` object that will be parsed over.
 * @param lcp_level The target compression level to deepen to.
 * @return 1 if deepening was successful, 0 otherwise.
 */
int lps_deepen(struct lps *str, int lcp_level);

/**
 * @brief Outputs the representation of a `lcp` pointer.
 *
 * This function iterates over the cores in the `cores` array and outputs them.
 * It prints each cores as bit representation, after initially printing the LCP level.
 *
 * @param element A pointer to the `lps` object to be output.
 */
void print_lps(const struct lps *str);

/**
 * @brief Equality operator for comparing two lcp::lps objects.
 *
 * This function compares the sizes of the core arrays of both objects, and then
 * compares each core element by dereferencing the core pointers. If all elements match,
 * the function returns true; otherwise, it returns false.
 *
 * @param lhs The left-hand side lps object to compare.
 * @param rhs The right-hand side lps object to compare.
 * @return 0 if both lps objects are equal, 1 otherwise.
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
 * @return 0 if the lps objects are not equal, 1 otherwise.
 */
int lps_neq(const struct lps *lhs, const struct lps *rhs);

#ifdef __cplusplus
}
#endif

#endif
