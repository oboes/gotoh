// this file contains some useful inline functions, mostly for managing memory

typedef int Vec[3];

static const size_t X = 0;
static const size_t Y = 1;
static const size_t Z = 2;
static const int Inf = 1024*1024*1024;

static inline Vec** alloc_arr(size_t x, size_t y) {
    Vec **arr = (Vec**) malloc( x * sizeof(Vec*) );
    if (arr == NULL) return NULL;
    for (size_t i=0 ; i<x ; i++) {
        arr[i] = (Vec*) malloc( y * sizeof(Vec) );
    }
    return arr;
}

static inline void free_arr(Vec **arr, size_t size) {
    if (arr != NULL) {
        for (size_t i=0 ; i<size ; i++) free(arr[i]);
        free(arr);
    }
}

static inline void fill_arr(int value, Vec **arr, size_t pos, size_t imin, size_t imax, size_t jmin, size_t jmax) {
    for (size_t i=imin ; i<=imax ; i++) {
        for (size_t j=jmin ; j<=jmax ; j++) {
            arr[i][j][pos] = value;
            arr[i][j][pos] = value;
            arr[i][j][pos] = value;
        }
    }
}

static inline int maximum(int arg0, int arg1, int arg2, int *idx) {

    if (arg0 >= arg1 && arg0 >= arg2) {
        *idx = X;
        return arg0;
    }
    if (arg1 >= arg0 && arg1 >= arg2) {
        *idx = Y;
        return arg1;
    }

    *idx = Z;
    return arg2;
}

