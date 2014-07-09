#ifndef __GOTOH_H__
#define __GOTOH_H__

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// dynamic programming arrays
typedef struct gth_Arr gth_Arr;
struct gth_Arr {
    size_t lenX;        // length of first sequence
    size_t lenY;        // length of second sequence
    int (**data)[3];    // arrays holding scores
    int (**path)[3];    // arrays holding backpointers
    int *gapX;          // gapX[i] = number of gaps between i-1 and i
    int *gapY;          // gapY[j] = number of gaps between j-1 and j
};

// sequence
typedef struct gth_Seq gth_Seq;
struct gth_Seq {
    size_t len;         // length
    char name[256];     // name
    char *res;          // residues
};

// substitution matrix
typedef struct gth_Sub gth_Sub;
struct gth_Sub {
    char alpha[27];     // alphabet
    int score[26][26];  // score matrix
};

// for the Gotoh algorithm
gth_Arr gth_init(size_t lenX, size_t lenY);
void    gth_set_sub(gth_Arr arr, const char *resX, const char *resY, int score[26][26]);
void    gth_set_gap(gth_Arr arr, int gap, int ext, int endgap, int endext);
int     gth_align(gth_Arr arr);
void    gth_free(gth_Arr arr);

// for input/output
gth_Seq gth_read_fasta(const char *filename);
gth_Sub gth_read_matrix(const char *filename);
void    gth_putseq(FILE *stream, const char *res, const int *gap);

#endif // __GOTOH_H__

