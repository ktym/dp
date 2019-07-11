# Pairwise sequence alignment algorithms

Currently implements Needleman-Wunsch, Smith-Waterman
and Needleman-Wunsch-Gotoh algorithms based on
the Dynamic Programming (DP).

## Coordinates

To treat (i,j) as 1-based coordinates (instead of Ruby's 0-based),
use the following accessors:

* DP::Matrix#get(i,j), DP::Matrix#set(i,j)
* DP::Sequence#[i], DP::Sequence#[j]

```
M(i,j) | Target A      C      C      A      G      T
-------+--------------------------------------------------
 Query | (0,0)  i=1    i=2    i=3    i=4    i=5    i=6
     A | j=1    (1,1)  (2,1)  (3,1)  (4,1)  (5,1)  (6,1)
     C | j=2    (1,2)  (2,2)  (3,2)  (4,2)  (5,2)  (6,2)
     A | j=3    (1,3)  (2,3)  (3,3)  (4,3)  (5,3)  (6,3)
     G | j=4    (1,4)  (2,4)  (3,4)  (4,4)  (5,4)  (6,4)
     C | j=5    (1,5)  (2,5)  (3,5)  (4,5)  (5,5)  (6,5)
```

## Accessors

Each sub-class has the following accessors to modify parameters
before the calculation of DP:

```
aligner = DP.new(target, query)
aligner.gap_penalty = 2            # gap open penalty
aligner.ext_penalty = 1            # gap extention penalty
aligner.match = 1
aligner.mis_match = -1
```

## Calculation

In the DP#dp method, init_matrix, calc_matrix and trace_back functions
are called which are required to be implemented in each sub-class.

```
aligner.dp
# => call init_matrix(), calc_matrix() and trace_back() then print alignment
```

## References

* Needleman SB, Wunsch CD: A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol 1970, 48:443-53.
* Smith TF, Waterman MS: Identification of common molecular subsequences. J Mol Biol 1981, 147:195-7.
* Gotoh O: An improved algorithm for matching biological sequences. J Mol Biol 1982, 162:705-8.

## Notes

* [Lecture notes by Terai](http://asailab.cb.k.u-tokyo.ac.jp/wp-content/uploads/2014/06/20181002_terai_pairwisealignment.mod2_.pdf)
* [Lecture notes by Asai](http://asailab.cb.k.u-tokyo.ac.jp/asai/lecture/H30Genome3_Alignment.pdf)
* [Lecture notes by Yokoyama](https://paper.dropbox.com/doc/--Af886MpyNXIKPsDWulpka5t1Ag-dqPm9SbT1DoLXnSyzu6lo)
* [DP basics in Ruby](https://www.jabba.cloud/20161020172918/)

## Copyright

MIT license (C) 2019 Toshiaki Katayama <ktym@dbcls.jp> and Toshiyuki Yokoyama
