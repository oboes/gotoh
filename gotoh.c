#include "gotoh.h"
#include "inline.c"


// allocate memory for the internal arrays
gth_Arr gth_init(size_t lenX, size_t lenY) {
    gth_Arr arr;
    arr.lenX = lenX;
    arr.lenY = lenY;
    arr.data = alloc_arr(lenX+1, lenY+1);
    arr.path = alloc_arr(lenX+1, lenY+1);
    arr.gapX = (int*) malloc( sizeof(int) * (lenX+1) );
    arr.gapY = (int*) malloc( sizeof(int) * (lenY+1) );
    if (arr.data == NULL || arr.path == NULL || arr.gapX == NULL || arr.gapY == NULL) gth_free(arr);
    return arr;
}


// initialize the substitution array Z
void gth_set_sub(gth_Arr arr, const char *resX, const char *resY, int score[26][26]) {
    for (size_t i=1 ; i<=arr.lenX ; i++) {
        for (size_t j=1 ; j<=arr.lenY ; j++) {
            arr.data[i][j][Z] = score[resX[i-1] - 'A'][resY[j-1] - 'A'];
        }
    }
}


// initialize the gap arrays X and Y
void gth_set_gap(gth_Arr arr, int gap, int ext, int endgap, int endext) {

    fill_arr(-gap, arr.data, X, 0,          0, 1, arr.lenY-1);
    fill_arr(-gap, arr.data, Y, 1, arr.lenX-1, 0,          0);

    fill_arr(-ext, arr.data, X, 1, arr.lenX  , 1, arr.lenY-1);
    fill_arr(-ext, arr.data, Y, 1, arr.lenX-1, 1, arr.lenY  );

    arr.data[       0][       0][X] = -endgap;
    arr.data[       0][arr.lenY][X] = -endgap;
    arr.data[       0][       0][Y] = -endgap;
    arr.data[arr.lenX][       0][Y] = -endgap;

    fill_arr(-endext, arr.data, X,        1, arr.lenX,        0,        0);
    fill_arr(-endext, arr.data, X,        1, arr.lenX, arr.lenY, arr.lenY);
    fill_arr(-endext, arr.data, Y,        0,        0,        1, arr.lenY);
    fill_arr(-endext, arr.data, Y, arr.lenX, arr.lenX,        1, arr.lenY);
}


// fill the arrays, backtrack, and store the alignment in arr.gapX and arr.gapY
int gth_align(gth_Arr arr) {

    for (size_t i=1 ; i<=arr.lenX ; i++) {
        arr.gapX[i] = arr.data[i][0][Y];
        arr.data[i][0][X] += arr.data[i-1][0][X];
    }

    for (size_t j=1 ; j<=arr.lenY ; j++) {
        arr.gapY[j] = arr.data[0][j][X];
        arr.data[0][j][Y] += arr.data[0][j-1][Y];
    }

    fill_arr(X, arr.path, X, 1, arr.lenX, 0,        0);
    fill_arr(Y, arr.path, Y, 0,        0, 1, arr.lenY);

    fill_arr(-Inf, arr.data, Y, 1, arr.lenX, 0,        0);
    fill_arr(-Inf, arr.data, Z, 1, arr.lenX, 0,        0);
    fill_arr(-Inf, arr.data, X, 0,        0, 1, arr.lenY);
    fill_arr(-Inf, arr.data, Z, 0,        0, 1, arr.lenY);

    // recurrence
    for (size_t i=1 ; i<=arr.lenX ; i++) {
        for (size_t j=1 ; j<=arr.lenY ; j++) {

            arr.data[i][j][X] += maximum(
                arr.data[i-1][j  ][X],
                arr.data[i-1][j  ][Y] + arr.gapY[j],
                arr.data[i-1][j  ][Z] + arr.gapY[j],
            &arr.path[i][j][X]);

            arr.data[i][j][Y] += maximum(
                arr.data[i  ][j-1][X] + arr.gapX[i],
                arr.data[i  ][j-1][Y],
                arr.data[i  ][j-1][Z] + arr.gapX[i],
            &arr.path[i][j][Y]);

            arr.data[i][j][Z] += maximum(
                arr.data[i-1][j-1][X],
                arr.data[i-1][j-1][Y],
                arr.data[i-1][j-1][Z],
            &arr.path[i][j][Z]);
        }
    }

    for (size_t i=0 ; i<=arr.lenX ; i++) arr.gapX[i] = 0;
    for (size_t j=0 ; j<=arr.lenY ; j++) arr.gapY[j] = 0;

    // backtracking
    int K;
    int max = maximum(
        arr.data[arr.lenX][arr.lenY][X],
        arr.data[arr.lenX][arr.lenY][Y],
        arr.data[arr.lenX][arr.lenY][Z],
    &K);
    size_t i = arr.lenX, j = arr.lenY;
    while (i > 0 || j > 0) {

        if (K == X) {
            arr.gapY[j]++;
            K = arr.path[i--][j  ][X];
        }
        else if (K == Y) {
            arr.gapX[i]++;
            K = arr.path[i  ][j--][Y];
        }
        else if (K == Z) {
            K = arr.path[i--][j--][Z];
        }
    }

    return max;
}


// free memory
void gth_free(gth_Arr arr) {
    free_arr(arr.data, arr.lenX+1);
    free_arr(arr.path, arr.lenX+1);
    free(arr.gapX);
    free(arr.gapY);
}


// read a sequence from a FASTA file
gth_Seq gth_read_fasta(const char *filename) {

    gth_Seq seq;
    seq.len = 0;
    seq.name[0] = '\0';
    seq.res = NULL;
    FILE *file = fopen(filename, "r");
    if (!file) return seq;

    // default sequence name is filename
    strncpy(seq.name, filename, sizeof(seq.name)-1);
    seq.name[sizeof(seq.name)-1] = '\0';

    // read sequence name
    char format[16];
    sprintf(format, " >%%%lus%%*[^\n]", sizeof(seq.name)-1);
    fscanf(file, format, seq.name);

    // get sequence size and allocate memory
    long offset = ftell(file);
    while (!feof(file)) {
        int c = fgetc(file);
        if (isalpha(c)) seq.len++;
    }
    if (seq.len == 0) return seq;
    seq.res = (char*) malloc(sizeof(char) * (seq.len+1));

    // read sequence
    if (seq.res != NULL) {
        fseek(file, offset, SEEK_SET);

        size_t i = 0;
        while (!feof(file) && i<seq.len) {
            int c = fgetc(file);
            if (isalpha(c)) seq.res[i++] = toupper(c);
        }
        seq.res[i] = '\0';
    }
    else seq.len = 0;

    fclose(file);
    return seq;
}


// read a substitution matrix from a NCBI/BAST file
gth_Sub gth_read_matrix(const char *filename) {

    gth_Sub matrix;
    matrix.alpha[0] = '\0';
    for (int i=0 ; i<sizeof(matrix.alpha-1) ; i++) {
        for (int j=0 ; j<sizeof(matrix.alpha-1) ; j++) {
            matrix.score[i][j] = 0;
        }
    }
    FILE *file = fopen(filename, "r");
    if (!file) return matrix;

    // skip lines beginning with a '#'
    char line[1024];
    while (fgets(line, sizeof(line), file) != NULL && line[0] == '#');
    
    // read the line with the alphabet
    size_t size = 0;
    for (char *p = strtok(line, " \t") ; p && isalpha(*p) && size<sizeof(matrix.alpha)-1 ; size++) {
        matrix.alpha[size] = toupper(*p);
        p = strtok(NULL, " \t");
    }

    // read the score matrix
    for (int i=0 ; i<size ; i++) {
        fscanf(file, " %*c");
        for (int j=0 ; j<size ; j++) {
            fscanf(file, "%d", &matrix.score[matrix.alpha[i] - 'A'][matrix.alpha[j] - 'A']);
        }
        fscanf(file, "%*[^\n]");
    }

    return matrix;
}


// print sequence res/gap in stream
void gth_putseq(FILE *stream, const char *res, const int *gap) {
    if (gap == NULL) {
        fputs(res, stream);
        return;
    }
    for (size_t k=0 ; k<gap[0] ; k++) fputc('-', stream);
    for (size_t i=0, len=strlen(res) ; i<len ; i++) {
        fputc(res[i], stream);
        for (size_t k=0 ; k<gap[i+1] ; k++) fputc('-', stream);
    }
}

