#include <stdio.h>
#include "gotoh.h"


const gth_Sub BLOSUM62;


int main(int argc, char **argv) {

    // arguments parsing
    
    int error = 0;
    int dump = 0;
    int quiet = 0;
    char *seqX_path = NULL;
    char *seqY_path = NULL;
    char *matrix_path = NULL;
    char *arr_paths[] = {NULL, NULL, NULL};
    double gapopen = 9.5;
    double gapextend = 0.5;
    double endopen = 0.0;
    double endextend = 0.0;

    for (size_t i=1 ; i<argc && !error ; i++) {
        if (!strcmp(argv[i], "-dump")) dump = 1;
        else if (!strcmp(argv[i], "-quiet")) quiet = 1;
        else if (argv[i][0] == '-') {
            if (i+1 >= argc) {
                fprintf(stderr, "ERROR: missing argument for option '%s'\n\n", argv[i]);
                error = 1;
            }
            else if (!strcmp(argv[i], "-matrix"))    matrix_path  = argv[i+1];
            else if (!strcmp(argv[i], "-gapopen"))   gapopen      = atof(argv[i+1]);
            else if (!strcmp(argv[i], "-gapextend")) gapextend    = atof(argv[i+1]);
            else if (!strcmp(argv[i], "-endopen"))   endopen      = atof(argv[i+1]);
            else if (!strcmp(argv[i], "-endextend")) endextend    = atof(argv[i+1]);
            else if (!strcmp(argv[i], "-arrxfile"))  arr_paths[0] = argv[i+1];
            else if (!strcmp(argv[i], "-arryfile"))  arr_paths[1] = argv[i+1];
            else if (!strcmp(argv[i], "-arrzfile"))  arr_paths[2] = argv[i+1];
            else {
                fprintf(stderr, "ERROR: unknown option '%s'\n\n", argv[i]);
                error = 1;
            }
            i++;
        }
        else if (seqX_path == NULL) seqX_path = argv[i];
        else if (seqY_path == NULL) seqY_path = argv[i];
        else {
            fprintf(stderr, "ERROR: too many arguments\n\n");
            error = 1;
        }
    }

    if (argc <= 1) error = 1;
    else if ((seqX_path == NULL || seqY_path == NULL) && !error) {
        fprintf(stderr, "ERROR: you must provide two sequence files\n\n");
        error = 1;
    }


    // help message

    if (error) {
        printf("Needleman-Wunsch global alignment of two sequences\n\n");
        printf("Usage: %s [OPTIONS] SEQUENCE1.fasta SEQUENCE2.fasta\n\n", argv[0]);
        printf("Options and their defaults:\n");
        printf("    -matrix     BLOSUM62    matrix file in NCBI/EMBOSS format\n");
        printf("    -gapopen    9.5         opening gap penalty\n");
        printf("    -gapextend  0.5         extending gap penalty\n");
        printf("    -endopen    0.0         opening end gap penalty\n");
        printf("    -endextend  0.0         extending end gap penalty\n");
        printf("    -quiet                  decrease verbosity\n\n");
        printf("Advanced:\n");
        printf("    -arrxfile   arrx.txt    load initial array X from this file\n");
        printf("    -arryfile   arry.txt    load initial array Y from this file\n");
        printf("    -arrzfile   arrz.txt    load initial array Z from this file\n");
        printf("    -dump                   output final arrays\n\n");
        printf("Default parameters are the same as those of the EMBOSS needle command.\n");
        return 0;
    }


    // load sequences

    gth_Seq seqX = gth_read_fasta(seqX_path);
    if (seqX.len == 0) {
        fprintf(stderr, "ERROR: problem reading file '%s'\n", seqX_path);
        return -1;
    }

    gth_Seq seqY = gth_read_fasta(seqY_path);
    if (seqY.len == 0) {
        fprintf(stderr, "ERROR: problem reading file '%s'\n", seqY_path);
        return -1;
    }


    // load matrix

    gth_Sub matrix = BLOSUM62;
    if (matrix_path != NULL) {
        matrix = gth_read_matrix(matrix_path);
        if (matrix.alpha[0] == '\0') {
            fprintf(stderr, "ERROR: problem reading file '%s'\n", matrix_path);
            return -1;
        }
    }
    for (size_t i=0 ; i<26 ; i++) {
        for (size_t j=0 ; j<26 ; j++) {
            matrix.score[i][j] *= 10;
        }
    }


    // create, fill, and backtrack the arrays

    gth_Arr array = gth_init(seqX.len, seqY.len);
    gth_set_sub(array, seqX.res, seqY.res, matrix.score);
    gth_set_gap(array, (int)(gapopen*10), (int)(gapextend*10), (int)(endopen*10), (int)(endextend*10));
    for (int k=0 ; k<3 ; k++) {
        if (arr_paths[k] == NULL) continue;
        FILE *file = fopen(arr_paths[k], "r");
        if (!file) {
            fprintf(stderr, "ERROR: problem reading file '%s'\n", arr_paths[k]);
            gth_free(array);
            free(seqX.res);
            free(seqY.res);
            return -1;
        }
        for (size_t i=0 ; i<=array.lenX ; i++) {
            for (size_t j=0 ; j<=array.lenY ; j++) {
                fscanf(file, "%d", &array.data[i][j][k]);
            }
            fscanf(file, "%*[^\n]");
        }
        fclose(file);
    }
    double score = (double)(gth_align(array)) / 10;


    // output

    if (!quiet) {
        printf("# Needleman-Wunsch global alignment of two sequences\n");
        printf("#\n");
        if (arr_paths[2] == NULL) {
            printf("# matrix file: %s\n", (matrix_path == NULL) ? "none specified, using BLOSUM62" : matrix_path);
        }
        else {
            printf("# substitution scores array provided by user\n");
        }
        if (arr_paths[0] == NULL && arr_paths[1] == NULL) {
            printf("# gap opening penalty: %.1f\n", gapopen);
            printf("# gap extending penalty: %.1f\n", gapextend);
            printf("# end gap opening penalty: %.1f\n", endopen);
            printf("# end gap extending penalty: %.1f\n", endextend);
        }
        else {
            printf("# gap scores array(s) provided by user\n");
        }
        printf("#\n");
        printf("# score: %.1f\n\n", score);
    }
    if (dump) {
        for (int k=0 ; k<3 ; k++) {
            printf("# scores in array %c:\n", 'X'+k);
            for (size_t i=0 ; i<=array.lenX ; i++) {
                for (size_t j=0 ; j<=array.lenY ; j++) {
                    printf("% 15d ", array.data[i][j][k]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }
    printf(">%s\n", seqX.name);
    gth_putseq(stdout, seqX.res, array.gapX);
    printf("\n\n>%s\n", seqY.name);
    gth_putseq(stdout, seqY.res, array.gapY);
    printf("\n");


    // cleanup

    gth_free(array);
    free(seqX.res);
    free(seqY.res);
    return 0;
}


const gth_Sub BLOSUM62 = {
    .alpha = "ARNDCQEGHILKMFPSTWYVBZX", .score = 
    {{  4, -2,  0, -2, -1, -2,  0, -2, -1, -4, -1, -1, -1, -2, -4, -1, -1, -1,  1,  0, -4,  0, -3,  0, -2, -1},
     { -2,  4, -3,  4,  1, -3, -1,  0, -3, -4,  0, -4, -3,  3, -4, -2,  0, -1,  0, -1, -4, -3, -4, -1, -3,  1},
     {  0, -3,  9, -3, -4, -2, -3, -3, -1, -4, -3, -1, -1, -3, -4, -3, -3, -3, -1, -1, -4, -1, -2, -2, -2, -3},
     { -2,  4, -3,  6,  2, -3, -1, -1, -3, -4, -1, -4, -3,  1, -4, -1,  0, -2,  0, -1, -4, -3, -4, -1, -3,  1},
     { -1,  1, -4,  2,  5, -3, -2,  0, -3, -4,  1, -3, -2,  0, -4, -1,  2,  0,  0, -1, -4, -2, -3, -1, -2,  4},
     { -2, -3, -2, -3, -3,  6, -3, -1,  0, -4, -3,  0,  0, -3, -4, -4, -3, -3, -2, -2, -4, -1,  1, -1,  3, -3},
     {  0, -1, -3, -1, -2, -3,  6, -2, -4, -4, -2, -4, -3,  0, -4, -2, -2, -2,  0, -2, -4, -3, -2, -1, -3, -2},
     { -2,  0, -3, -1,  0, -1, -2,  8, -3, -4, -1, -3, -2,  1, -4, -2,  0,  0, -1, -2, -4, -3, -2, -1,  2,  0},
     { -1, -3, -1, -3, -3,  0, -4, -3,  4, -4, -3,  2,  1, -3, -4, -3, -3, -3, -2, -1, -4,  3, -3, -1, -1, -3},
     { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4},
     { -1,  0, -3, -1,  1, -3, -2, -1, -3, -4,  5, -2, -1,  0, -4, -1,  1,  2,  0, -1, -4, -2, -3, -1, -2,  1},
     { -1, -4, -1, -4, -3,  0, -4, -3,  2, -4, -2,  4,  2, -3, -4, -3, -2, -2, -2, -1, -4,  1, -2, -1, -1, -3},
     { -1, -3, -1, -3, -2,  0, -3, -2,  1, -4, -1,  2,  5, -2, -4, -2,  0, -1, -1, -1, -4,  1, -1, -1, -1, -1},
     { -2,  3, -3,  1,  0, -3,  0,  1, -3, -4,  0, -3, -2,  6, -4, -2,  0,  0,  1,  0, -4, -3, -4, -1, -2,  0},
     { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4},
     { -1, -2, -3, -1, -1, -4, -2, -2, -3, -4, -1, -3, -2, -2, -4,  7, -1, -2, -1, -1, -4, -2, -4, -2, -3, -1},
     { -1,  0, -3,  0,  2, -3, -2,  0, -3, -4,  1, -2,  0,  0, -4, -1,  5,  1,  0, -1, -4, -2, -2, -1, -1,  3},
     { -1, -1, -3, -2,  0, -3, -2,  0, -3, -4,  2, -2, -1,  0, -4, -2,  1,  5, -1, -1, -4, -3, -3, -1, -2,  0},
     {  1,  0, -1,  0,  0, -2,  0, -1, -2, -4,  0, -2, -1,  1, -4, -1,  0, -1,  4,  1, -4, -2, -3,  0, -2,  0},
     {  0, -1, -1, -1, -1, -2, -2, -2, -1, -4, -1, -1, -1,  0, -4, -1, -1, -1,  1,  5, -4,  0, -2,  0, -2, -1},
     { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4},
     {  0, -3, -1, -3, -2, -1, -3, -3,  3, -4, -2,  1,  1, -3, -4, -2, -2, -3, -2,  0, -4,  4, -3, -1, -1, -2},
     { -3, -4, -2, -4, -3,  1, -2, -2, -3, -4, -3, -2, -1, -4, -4, -4, -2, -3, -3, -2, -4, -3, 11, -2,  2, -3},
     {  0, -1, -2, -1, -1, -1, -1, -1, -1, -4, -1, -1, -1, -1, -4, -2, -1, -1,  0,  0, -4, -1, -2, -1, -1, -1},
     { -2, -3, -2, -3, -2,  3, -3,  2, -1, -4, -2, -1, -1, -2, -4, -3, -1, -2, -2, -2, -4, -1,  2, -1,  7, -2},
     { -1,  1, -3,  1,  4, -3, -2,  0, -3, -4,  1, -3, -1,  0, -4, -1,  3,  0,  0, -1, -4, -2, -3, -1, -2,  4}}
};

