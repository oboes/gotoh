Gotoh
=====


Short (around 500 lines) implementation of the Gotoh algorithm, also known as the
[Needleman-Wunsch](http://en.wikipedia.org/wiki/Needleman-Wunsch) algorithm with affine gap penalties.
The main difference with classical Needleman-Wunsch is that Gotoh requires three different
dynamic programming arrays.


### Features ###

* C99-compliant, only uses `stdio.h`, `stdlib.h`, `string.h`, and `ctype.h`
* basic command-line interface
* read FASTA files
* read NCBI/EMBOSS scoring matrix files
* exact same results as [EMBOSS Needle](http://www.ebi.ac.uk/Tools/psa/emboss_needle/)
* initial values in the dynamic programming arrays can be loaded from files

Due to the last point, it should be possible to experiment with different scoring schemes, for example:
* gap opening penalties depending on the residues surrounding the gap
* substitution scores depending on the matched residues *and* their positions in the sequence
* doing something else than biological sequence alignment, such as dynamic time warping


### Usage ###

Compile with this command:
```
gcc -std=c99 gotoh.c main.c -o align
```
Then type `./align` and read the help message.


### Reference ###
Durbin, Eddy, Krogh, Mitchison. [Biological Sequence Analysis](http://books.google.com/books?id=R5P2GlJvigQC),
*Cambridge University Press*, 1998.

