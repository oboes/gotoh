minimalist implementation of the Gotoh algorithm
================================================

A very short (around 500 lines) implementation of the Gotoh algorithm, also known as the
[Needleman-Wunsch](http://en.wikipedia.org/wiki/Needleman-Wunsch) algorithm with affine gap penalties.
It allows arbitrary initial values in the three dynamic programming arrays, and so it is possible
to use the algorithm for e.g. :
* gap opening penalties depending on the residues surrounding the gap
* substitution scores depending on the matched residues AND their positions in the sequence
* doing something else than biological sequence alignment, such as dynamic time warping
With its default parameters, the program will give the exact same result as BLAST needle


### Compiling ###
```
gcc -std=c99 gotoh.c main.c -o gotoh
```



