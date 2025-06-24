/**
 * @file encoding.c
 * @brief Implementation of encoding functions.
 *
 * This file contains the implementation of encoding functions used to
 * initialize the alphabet with their corresponding coefficients. The encodings
 * support initialization with default coefficients, specific coefficients, or
 * by reading coefficients from a file.
 */

#include "encoding.h"

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
