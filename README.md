Short C implementation (around 500 lines of codes + comments) of the Gotoh algorithm, also known as the
[Needleman-Wunsch](http://en.wikipedia.org/wiki/Needleman-Wunsch) algorithm with affine gap penalties.

This software was developed for my 2014 bioinformatics master thesis, [available here](./doc/thesis.pdf).
See also the [defense slides](./doc/slides.pdf) for a high-level summary.

Also included in this repo is this [small interactive explanation](http://oboes.github.io/gotoh/doc/needleman-wunsch.html) of the Needleman-Wunsch algorithm.


### Features ###

* C99-compliant, only uses `stdio.h`, `stdlib.h`, `string.h`, and `ctype.h`
* basic command-line interface
* read FASTA files
* read NCBI/EMBOSS scoring matrix files
* exact same results as [EMBOSS Needle](http://www.ebi.ac.uk/Tools/psa/emboss_needle/)
* initial values in the dynamic programming arrays can be loaded from files


Due to the last point, the program can be used for experimenting with different scoring schemes, for example:
* gap opening penalties depending on the residues surrounding the gap
* substitution scores depending on the matched residues *and* their positions in the sequence
* doing something else than biological sequence alignment, such as dynamic time warping


### Usage ###

Compile with this command:
```
gcc -std=c99 gotoh.c main.c -o align
```
Then type `./align` and read the help message.

Basic example showing how to use the code:
```C
gth_Seq sequenceX = gth_read_fasta("first_sequence.fasta");
gth_Seq sequenceY = gth_read_fasta("second_sequence.fasta");
gth_Sub matrix = gth_read_matrix("BLOSUM62");
int gapopen=10, gapextend=1, endopen=0, endextend=0;

gth_Arr array = gth_init(sequenceX.len, sequenceY.len);
gth_set_sub(array, sequenceX.res, sequenceY.res, matrix.score);
gth_set_gap(array, gapopen, gapextend, endopen, endextend);
int score = gth_align(array);

printf("Score: %d\n", score);
gth_putseq(stdout, sequenceX.res, array.gapX);
printf("\n");
gth_putseq(stdout, sequenceY.res, array.gapY);
printf("\n");

gth_free(array);
free(sequenceX.res);
free(sequenceY.res);
```


### Citing this work ###

If you use this software, you can cite my bioinformatics master thesis:

```
@mastersthesis{oboes-dynamine,
    author    = "Olivier Boes",
    title     = "Improving the Needleman-Wunsch algorithm with the DynaMine predictor",
    school    = "Université Libre de Bruxelles",
    year      = "2014"
}
```
