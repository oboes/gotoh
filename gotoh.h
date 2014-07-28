#ifndef __GOTOH_H__
#define __GOTOH_H__

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// data structure for the dynamic programming arrays
typedef struct gth_Arr gth_Arr;
struct gth_Arr {
    size_t lenX;        // length of first sequence
    size_t lenY;        // length of second sequence
    int (**data)[3];    // arrays holding scores
    int (**path)[3];    // arrays holding backpointers
    int *gapX;          // gapX[i] = number of gaps between i-1 and i
    int *gapY;          // gapY[j] = number of gaps between j-1 and j
};


/*  *******************************  EXPLANATION OF THE gth_Arr DATA STRUCTURE  *******************************
 *
 *  Sequence 1 is called seqX and its residues are numbered from 1 to lenX
 *  Sequence 2 is called seqY and its residues are numbered from 1 to lenY
 *  So indices start at 1 rather than 0. We use i for seqX residues and j for seqY residues.
 *
 *  A sequence alignment is something like this (here with lenX=9 and lenY=8):
 *      -123--456-789
 *      123-456--78--
 *  Where numerical indices are the residues and '-' symbols are the gaps
 *  A partial alignment between subsequences seqX[1..i] and seqY[1..j] is called a (i,j)-subalignment
 *
 *  The actual residues in the sequences (e.g. amino acids or nucleic acids)
 *  are not stored in the gth_Arr data structure: the algorithm doesn't need them
 *  Instead it uses gap penalties and substitution scores stored in the arr.data arrays
 *
 *  For all (0 <= i <= lenX,  0 <= j <= lenY) let us define:
 *      X[i][j]  =  arr.data[i][j][0]  (seqX insertion penalties = seqY gap penalties)
 *      Y[i][j]  =  arr.data[i][j][1]  (seqY insertion penalties = seqX gap penalties)
 *      Z[i][j]  =  arr.data[i][j][2]  (substitution scores)
 *  These three matrices must contain all the parameters used by the algorithm
 *
 *  Before calling gth_align you should set the following scores:
 *      X[0][j]  =  score for opening in seqY a gap between j and j+1  (0 <= j <= lenY)
 *      X[i][j]  =  score for matching i with a gap between j and j+1  (0 <= j <= lenY,  1 <= i <= lenX)
 *      Y[i][0]  =  score for opening in seqX a gap between i and i+1  (0 <= i <= lenX)
 *      Y[i][j]  =  score for matching j with a gap between i and i+1  (0 <= i <= lenX,  1 <= j <= lenY)
 *      Z[i][j]  =  score for matching i in seqX with j in seqY        (1 <= i <= lenX,  1 <= j <= lenY)
 *  You don't need to set the values of Z[0][j] and Z[i][0], they are not used by the algorithm
 *
 *  Once gth_align has been called, the matrices contain subalignments scores:
 *      X[i][j]  =  max score for (i,j)-subalignments ending with i matched to a gap
 *      Y[i][j]  =  max score for (i,j)-subalignments ending with j matched to a gap
 *      Z[i][j]  =  max score for (i,j)-subalignments ending with i matched to j
 *  Return value = optimal alignment score = max( X[lenX][lenY], Y[lenX][lenY], Z[lenX][lenY] )
 *
 *  Corresponding gap lengths for an optimal alignment are:
 *      arr.gapX[i]  =  number of gaps between i-1 and i in seqX  (0 <= i <= lenX)
 *      arr.gapY[j]  =  number of gaps between j-1 and j in seqY  (0 <= j <= lenY)
 *
 *  You don't need to pay attention to the arr.path arrays: they are internally used for
 *  storing indices that will be used for backtracking once the score arrays are computed.
 *  It should be possible to modify the code so that the path.arrays indices contain
 *  information on all possible backtracking directions (rather than one single arbitrary one),
 *  which could then be used to retrieve *all* optimal alignments having the same optimal score.
 *  But you will need to read and modify the actual code if you want to do such a thing.
 *
 */


// data structure for an ungapped sequence (read it from a FASTA file using gth_read_fasta)
typedef struct gth_Seq gth_Seq;
struct gth_Seq {
    size_t len;         // number of residues in the sequences (without the terminating NULL character)
    char name[256];     // name of the sequence
    char *res;          // sequence residues (should be A-Z letters); don't forget to free the memory once not needed anymore
};


// data structure for a substitution matrix (read it from a NCBI/EMBOSS file using gth_read_matrix)
typedef struct gth_Sub gth_Sub;
struct gth_Sub {
    char alpha[27];     // alphabet used in this matrix (letters are replaced by their position in the alphabet)
    int score[26][26];  // score[i][j] = substitution score for matching letter i with letter j
};


// functions implementing the Needleman-Wunsch-Gotoh algorithm
gth_Arr gth_init(size_t lenX, size_t lenY);
void    gth_set_sub(gth_Arr arr, const char *resX, const char *resY, int score[26][26]);
void    gth_set_gap(gth_Arr arr, int gap, int ext, int endgap, int endext);
int     gth_align(gth_Arr arr);
void    gth_free(gth_Arr arr);


// functions handling basic input/ouput
gth_Seq gth_read_fasta(const char *filename);
gth_Sub gth_read_matrix(const char *filename);
void    gth_putseq(FILE *stream, const char *res, const int *gap);

/*  See the GitHub repository for a working example on how to use the library!  */

#endif // __GOTOH_H__

